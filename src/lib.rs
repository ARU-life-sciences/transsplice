pub mod assign;
pub mod candidate;
pub mod template;

use anyhow::{Context, Result};
use assign::{best_assignment, best_assignments_beam, classify_junction, Assignment};
use candidate::{build_candidates, Candidate, RawFragment};
use orfedit::align::{align, ColumnCall};
use orfedit::codon::AA_ORDER;
use orfedit::fasta::{codons_for_frame, FastaIndex};
use orfedit::{scan_gene, ScanParams};
use template::GeneTemplate;

pub struct ExonResult {
    pub slot: usize, // 1-indexed exon number
    pub provenance: String,
    pub contig: String,
    pub start: usize,
    pub end: usize,
    pub strand: char,
    pub score: f64,
    /// Real, genomically-located predicted edit sites for this exon
    /// (straight from `orfedit::scan_gene`'s own `ScanOutcome.edits` -
    /// not re-derived, not just a count).
    pub edits: Vec<orfedit::EditCall>,
    pub via_denovo: bool,
    /// GFF3 CDS phase: bases at the exon's 5' end (transcript orientation)
    /// that complete a codon begun in the previous exon.
    pub phase: usize,
}

/// Bonus added to the whole-gene score for each intron end matching the
/// group II consensus: 5' `GNGCG` (15/18 Arabidopsis mito cis introns),
/// 3' `AY`. Calibrated on Arabidopsis (nad1/2/4/5/7, ccmFc, cox2, rpl2,
/// rps3 vs RefSeq): 2.0 got 9/32 exons base-exact - a split codon's two
/// readings differ by one residue, and the profile often prefers the wrong
/// one by more than a tiebreak; 20.0 got 21/32 and 40/100 changed nothing,
/// with protein identity identical throughout (a chance motif within the
/// few-bp refinement window is rare).
pub const SPLICE_MOTIF_BONUS: f64 = 20.0;

/// How far (in codons) each side of a junction may move from the
/// codon-aligned exon match during junction refinement.
const JUNCTION_CODON_SHIFT: i64 = 1;

/// One assigned exon before junction refinement: the codon-aligned match
/// from `scan_gene` (its "core"), which refinement extends by `ext5`/`ext3`
/// bases at the transcript 5'/3' ends (negative = trimmed back into the core).
struct Placed {
    slot: usize,
    contig: String,
    strand: char,
    core_start: i64,
    core_end: i64,
    score: f64,
    edits: Vec<orfedit::EditCall>,
    provenance: String,
    via_denovo: bool,
}

impl Placed {
    fn span(&self, ext5: i64, ext3: i64) -> (i64, i64) {
        if self.strand == '-' {
            (self.core_start - ext3, self.core_end + ext5)
        } else {
            (self.core_start - ext5, self.core_end + ext3)
        }
    }
}

fn sense_bases(fasta_idx: &FastaIndex, contig: &str, start: i64, end: i64, strand: char) -> Vec<u8> {
    fasta_idx
        .extract_window(contig, start, end, strand)
        .map(|(_, _, s)| s.to_ascii_uppercase())
        .unwrap_or_default()
}

/// Group II intron consensus at a junction between exon `a` (whose 3' end
/// the intron follows) and exon `b` (whose 5' end it precedes). Holds for
/// trans-splicing too: each half of a split intron keeps its own end.
fn splice_motif_bonus(fasta_idx: &FastaIndex, a: &Placed, a_span: (i64, i64), b: &Placed, b_span: (i64, i64)) -> f64 {
    let five = if a.strand == '-' {
        sense_bases(fasta_idx, &a.contig, a_span.0 - 5, a_span.0, '-')
    } else {
        sense_bases(fasta_idx, &a.contig, a_span.1, a_span.1 + 5, '+')
    };
    let three = if b.strand == '-' {
        sense_bases(fasta_idx, &b.contig, b_span.1, b_span.1 + 2, '-')
    } else {
        sense_bases(fasta_idx, &b.contig, b_span.0 - 2, b_span.0, '+')
    };
    let mut bonus = 0.0;
    if five.len() == 5 && five[0] == b'G' && &five[2..5] == b"GCG" {
        bonus += SPLICE_MOTIF_BONUS;
    }
    if three.len() == 2 && three[0] == b'A' && (three[1] == b'C' || three[1] == b'T') {
        bonus += SPLICE_MOTIF_BONUS;
    }
    bonus
}

pub struct JunctionResult {
    pub from_slot: usize,
    pub to_slot: usize,
    pub kind: &'static str,
    pub gap_bp: i64,
}

pub struct ReconstructionOutcome {
    pub n_slots: usize,
    pub n_filled: usize,
    pub complete: bool,
    pub exons: Vec<ExonResult>,
    pub junctions: Vec<JunctionResult>,
    pub whole_gene_score: f64,
    pub whole_gene_n_edits: usize,
    pub protein: String,
}

/// Reconstruct one (species, gene) call from its raw `.ctg.bed` fragments.
/// See transsplice's README for the algorithm; this ties together
/// candidate generation, per-slot scoring (via `orfedit::scan_gene`, reused
/// as-is - each candidate-vs-slot score IS a boundary-corrected alignment,
/// not just a raw score), assignment, a de novo fallback scan for slots no
/// existing fragment fills, junction classification, and a final whole-gene
/// alignment pass for an aggregate score/protein.
pub fn reconstruct_gene(
    fasta_idx: &FastaIndex,
    fragments: &[RawFragment],
    template: &GeneTemplate,
    params: &ScanParams,
    merge_distance: i64,
    junction_distance_threshold: i64,
    exon_flank: i64,
) -> Result<ReconstructionOutcome> {
    let k = template.slots.len();
    let mut candidates = build_candidates(fragments, merge_distance);

    let exon_params = ScanParams { flank: exon_flank, ..*params };
    let mut scores: Vec<Vec<f64>> = candidates
        .iter()
        .map(|cand| {
            template
                .slots
                .iter()
                .map(|slot| {
                    scan_gene(fasta_idx, &cand.contig, cand.start, cand.end, cand.strand, &slot.profile, &exon_params)
                        .map(|o| o.score)
                        .unwrap_or(f64::NEG_INFINITY)
                })
                .collect()
        })
        .collect();

    let floors: Vec<f64> = template.slots.iter().map(|s| s.min_self_score).collect();
    let mut assignment = best_assignment(&scores, &floors);

    // De novo fallback: any slot no candidate could fill gets a whole-contig
    // scan (flank=0 => the whole requested region is the search window),
    // restricted to contigs already known to carry this gene in this
    // species - not a genome-wide scan.
    let mut contigs: Vec<String> = fragments.iter().map(|f| f.contig.clone()).collect();
    contigs.sort();
    contigs.dedup();
    let denovo_params = ScanParams { flank: 0, ..*params };
    let mut denovo_candidate_indices: std::collections::HashSet<usize> = std::collections::HashSet::new();

    for (si, slot) in template.slots.iter().enumerate() {
        if assignment.slot_to_candidate[si].is_some() {
            continue;
        }
        let mut best: Option<(f64, i64, i64, String, char)> = None;
        for contig in &contigs {
            let Some(len) = fasta_idx.sequences.get(contig).map(|s| s.len() as i64) else { continue };
            for strand in ['+', '-'] {
                if let Ok(o) = scan_gene(fasta_idx, contig, 0, len, strand, &slot.profile, &denovo_params) {
                    if o.score >= slot.min_self_score && best.as_ref().is_none_or(|b| o.score > b.0) {
                        best = Some((o.score, o.corrected_start as i64, o.corrected_end as i64, contig.clone(), strand));
                    }
                }
            }
        }
        if let Some((score, start, end, contig, strand)) = best {
            let idx = candidates.len();
            candidates.push(Candidate { contig, start, end, strand, provenance: "denovo".to_string() });
            for row in scores.iter_mut() {
                row.push(f64::NEG_INFINITY);
            }
            let mut row = vec![f64::NEG_INFINITY; k];
            row[si] = score;
            scores.push(row);
            assignment.slot_to_candidate[si] = Some(idx);
            assignment.total_score += score;
            denovo_candidate_indices.insert(idx);
        }
    }

    // Per-slot-independent scoring can rank a wrong permutation of several
    // real, decently-scoring candidates above the right one (confirmed
    // against ground truth - see README) - so evaluate every plausible
    // assignment in a small beam by its effect on the WHOLE protein, not
    // just accept the single locally-greedy one.
    // Beam members are compared with template-default junctions; only the
    // winner gets junction refinement (the costly part - many whole-gene
    // alignments), which fine-tunes boundaries rather than changing which
    // candidates are exons.
    let beam = best_assignments_beam(&scores, &floors, 2);
    let mut best: Option<(usize, (Vec<ExonResult>, Vec<JunctionResult>, f64, usize, String))> = None;
    for (bi, candidate_assignment) in beam.iter().enumerate() {
        let Ok(result) = evaluate_assignment(
            fasta_idx, &candidates, template, candidate_assignment, &exon_params, params,
            &denovo_candidate_indices, junction_distance_threshold, false,
        ) else {
            continue;
        };
        if best.as_ref().is_none_or(|(_, prev)| result.2 > prev.2) {
            best = Some((bi, result));
        }
    }
    let (bi, unrefined) = best.context("no candidate assignment (not even the empty one) produced a result")?;
    let (exons, junctions, whole_gene_score, whole_gene_n_edits, protein) = evaluate_assignment(
        fasta_idx, &candidates, template, &beam[bi], &exon_params, params,
        &denovo_candidate_indices, junction_distance_threshold, true,
    )
    .unwrap_or(unrefined);

    Ok(ReconstructionOutcome {
        n_slots: k,
        n_filled: exons.len(),
        complete: exons.len() == k,
        exons,
        junctions,
        whole_gene_score,
        whole_gene_n_edits,
        protein,
    })
}

/// Build the exon list, junction classifications, and whole-gene alignment
/// for one candidate assignment - the thing `reconstruct_gene`'s beam
/// search evaluates once per assignment it's choosing between.
#[allow(clippy::too_many_arguments)]
fn evaluate_assignment(
    fasta_idx: &FastaIndex,
    candidates: &[Candidate],
    template: &GeneTemplate,
    assignment: &Assignment,
    exon_params: &ScanParams,
    params: &ScanParams,
    denovo_candidate_indices: &std::collections::HashSet<usize>,
    junction_distance_threshold: i64,
    refine_junctions: bool,
) -> Result<(Vec<ExonResult>, Vec<JunctionResult>, f64, usize, String)> {
    let k = template.slots.len();
    let mut placed: Vec<Placed> = Vec::new();
    for (si, slot) in template.slots.iter().enumerate() {
        let Some(ci) = assignment.slot_to_candidate[si] else { continue };
        let cand = &candidates[ci];
        let outcome = scan_gene(fasta_idx, &cand.contig, cand.start, cand.end, cand.strand, &slot.profile, exon_params)?;
        let (mut core_start, mut core_end) = (outcome.corrected_start as i64, outcome.corrected_end as i64);
        let mut edits = outcome.edits;
        // scan_gene folds in a following stop codon (possibly edit-created).
        // Only the gene's last exon ends in one; after any other exon that
        // "codon" is split-codon/intron sequence, and CGA/CAA/CAG read as
        // edit-created stops often enough to extend exons wrongly.
        if outcome.includes_stop && si + 1 != k {
            if cand.strand == '-' {
                core_start += 3;
            } else {
                core_end -= 3;
            }
            edits.retain(|e| e.resulting_aa != '*');
        }
        placed.push(Placed {
            slot: si + 1,
            contig: cand.contig.clone(),
            strand: cand.strand,
            core_start,
            core_end,
            score: outcome.score,
            edits,
            provenance: cand.provenance.clone(),
            via_denovo: denovo_candidate_indices.contains(&ci),
        });
    }

    // Two slots can't be the same stretch of genome. Assignment only stops a
    // candidate being reused, not two overlapping candidates (a fragment and
    // a merged run containing it) filling different slots - seen with a
    // short slot profile matching inside its neighbour's exon (Arabidopsis
    // cox2 exon 2, 27 aa). Keep the better-scoring exon; the slot empties.
    let mut dropped = vec![false; placed.len()];
    for i in 0..placed.len() {
        for j in (i + 1)..placed.len() {
            let (a, b) = (&placed[i], &placed[j]);
            if a.contig == b.contig && a.core_start < b.core_end && b.core_start < a.core_end {
                if a.score >= b.score { dropped[j] = true } else { dropped[i] = true }
            }
        }
    }
    let placed: Vec<Placed> = placed.into_iter().zip(dropped).filter(|(_, d)| !d).map(|(p, _)| p).collect();

    // Junction refinement. Slot profiles cover whole codons only, so each
    // core is short of its true exon by the split-codon bases (template
    // lead/trail); adding them back keeps the joined CDS in frame. The
    // split point itself varies a little between species, so at every
    // junction between adjacent slots try each split-codon distribution
    // (3' bases t on the upstream exon, (3-t)%3 on the downstream one) and
    // a one-codon shift either side, keeping the placement the whole-gene
    // alignment scores best (intron-motif bonus as tiebreak).
    let mut ext5: Vec<i64> = placed.iter().map(|p| template.slots[p.slot - 1].lead as i64).collect();
    let mut ext3: Vec<i64> = placed.iter().map(|p| template.slots[p.slot - 1].trail as i64).collect();
    for i in 0..placed.len().saturating_sub(1) {
        if placed[i + 1].slot == placed[i].slot + 1 {
            ext5[i + 1] = (3 - ext3[i]).rem_euclid(3);
        }
    }

    let concat_of = |ext5: &[i64], ext3: &[i64]| -> Vec<u8> {
        let mut concat = Vec::new();
        for (i, p) in placed.iter().enumerate() {
            let (s, e) = p.span(ext5[i], ext3[i]);
            if e > s {
                concat.extend(sense_bases(fasta_idx, &p.contig, s, e, p.strand));
            }
        }
        concat
    };
    let score_of = |ext5: &[i64], ext3: &[i64]| -> f64 {
        let codons = codons_for_frame(&concat_of(ext5, ext3), 0);
        align(&template.whole_gene_profile, &codons, params.gap_open, params.gap_extend, params.edit_penalty).score
    };
    let junction_bonus = |i: usize, ext5: &[i64], ext3: &[i64]| -> f64 {
        let (a, b) = (&placed[i], &placed[i + 1]);
        splice_motif_bonus(fasta_idx, a, a.span(ext5[i], ext3[i]), b, b.span(ext5[i + 1], ext3[i + 1]))
    };

    // Coordinate descent per junction: split-codon distribution first, then
    // a codon shift on the upstream exon, then on the downstream one (9
    // whole-gene alignments per junction instead of all 27 combinations).
    for i in 0..placed.len().saturating_sub(1) {
        if !refine_junctions || placed[i + 1].slot != placed[i].slot + 1 {
            continue; // a missing slot between them: no shared codon to place
        }
        let floor_a = (-(placed[i].core_end - placed[i].core_start) + 3).max(-3); // keep >= one codon
        let floor_b = (-(placed[i + 1].core_end - placed[i + 1].core_start) + 3).max(-3);
        let eval = |t: i64, da: i64, db: i64| -> Option<f64> {
            let (e3, e5) = (t + 3 * da, (3 - t).rem_euclid(3) + 3 * db);
            if e3 < floor_a || e5 < floor_b {
                return None;
            }
            let (mut t5, mut t3) = (ext5.clone(), ext3.clone());
            t3[i] = e3;
            t5[i + 1] = e5;
            Some(score_of(&t5, &t3) + junction_bonus(i, &t5, &t3))
        };
        let argmax = |opts: Vec<(i64, Option<f64>)>, default: i64| -> i64 {
            opts.into_iter()
                .filter_map(|(v, s)| s.map(|s| (v, s)))
                .fold(None::<(i64, f64)>, |b, (v, s)| if b.is_none_or(|(_, bs)| s > bs) { Some((v, s)) } else { b })
                .map_or(default, |(v, _)| v)
        };
        let t0 = ext3[i].rem_euclid(3);
        let t = argmax((0..3).map(|t| (t, eval(t, 0, 0))).collect(), t0);
        let da = argmax((-JUNCTION_CODON_SHIFT..=JUNCTION_CODON_SHIFT).map(|d| (d, eval(t, d, 0))).collect(), 0);
        let db = argmax((-JUNCTION_CODON_SHIFT..=JUNCTION_CODON_SHIFT).map(|d| (d, eval(t, da, d))).collect(), 0);
        if eval(t, da, db).is_some() {
            ext3[i] = t + 3 * da;
            ext5[i + 1] = (3 - t).rem_euclid(3) + 3 * db;
        }
    }

    let exons: Vec<ExonResult> = placed
        .into_iter()
        .enumerate()
        .map(|(i, p)| {
            let (s, e) = p.span(ext5[i], ext3[i]);
            ExonResult {
                slot: p.slot,
                provenance: p.provenance,
                contig: p.contig,
                start: s.max(0) as usize,
                end: e.max(0) as usize,
                strand: p.strand,
                score: p.score,
                edits: p.edits,
                via_denovo: p.via_denovo,
                phase: ext5[i].rem_euclid(3) as usize,
            }
        })
        .collect();

    let mut junctions = Vec::new();
    for w in exons.windows(2) {
        let (a, b) = (&w[0], &w[1]);
        let (kind, gap) = classify_junction(
            &a.contig, a.start as i64, a.end as i64, a.strand,
            &b.contig, b.start as i64, b.end as i64, b.strand,
            junction_distance_threshold,
        );
        junctions.push(JunctionResult { from_slot: a.slot, to_slot: b.slot, kind: kind.as_str(), gap_bp: gap });
    }

    // Whole-gene pass: concatenate the corrected exon sequences (in slot
    // order, each already strand-corrected/frame-trimmed by its own
    // scan_gene call above) and align once against the whole-gene profile -
    // this is the arbiter the beam search picks an assignment by.
    let mut concat: Vec<u8> = Vec::new();
    for e in &exons {
        let (_, _, seq) = fasta_idx.extract_window(&e.contig, e.start as i64, e.end as i64, e.strand)?;
        concat.extend(seq);
    }
    let codons = codons_for_frame(&concat, 0);
    let whole = align(&template.whole_gene_profile, &codons, params.gap_open, params.gap_extend, params.edit_penalty);
    let mut protein = String::with_capacity(whole.columns.len());
    let mut whole_gene_n_edits = 0usize;
    for call in &whole.columns {
        match call {
            ColumnCall::Aligned { aa, n_edits, .. } => {
                protein.push(*aa as char);
                if *n_edits > 0 {
                    whole_gene_n_edits += 1;
                }
            }
            ColumnCall::Deleted => protein.push('-'),
        }
    }
    let _ = AA_ORDER; // re-exported for downstream convenience if needed

    Ok((exons, junctions, whole.score, whole_gene_n_edits, protein))
}

#[cfg(test)]
mod tests {
    use super::*;
    use orfedit::fasta::FastaIndex;
    use orfedit::profile::Profile;
    use std::collections::HashMap;
    use template::SlotTemplate;

    fn mkv_profile(gene: &str) -> Profile {
        use orfedit::codon::AA_ORDER;
        let columns = [b'M', b'K', b'V']
            .iter()
            .map(|&target| {
                let mut col = [-5.0f64; 20];
                col[AA_ORDER.iter().position(|&x| x == target).unwrap()] = 5.0;
                col
            })
            .collect();
        Profile { gene: gene.to_string(), columns }
    }

    #[test]
    fn exon_edit_sites_carry_through_to_the_result() {
        // "TTT" junk, "ACG" (Thr, unedited - needs 1 edit to become Met),
        // "AAA" (Lys), "GTT" (Val), "CCC" junk.
        let seq = b"TTTACGAAAGTTCCC".to_vec();
        let mut sequences = HashMap::new();
        sequences.insert("ctg1".to_string(), seq);
        let fasta_idx = FastaIndex { sequences };

        let template = GeneTemplate {
            gene: "test".to_string(),
            slots: vec![SlotTemplate { profile: mkv_profile("test_exon1"), min_self_score: -100.0, lead: 0, trail: 0 }],
            whole_gene_profile: mkv_profile("test_whole"),
        };
        let fragments = vec![RawFragment { contig: "ctg1".to_string(), start: 3, end: 12, strand: '+' }];
        let params = ScanParams::default();

        let outcome = reconstruct_gene(&fasta_idx, &fragments, &template, &params, 0, 1000, 0).unwrap();

        assert_eq!(outcome.n_filled, 1);
        let exon = &outcome.exons[0];
        assert_eq!(exon.edits.len(), 1, "expected exactly one edit (the ACG->ATG rescue)");
        assert_eq!(exon.edits[0].resulting_aa, 'M');
        assert_eq!(exon.edits[0].genomic_pos, 4); // the C in "ACG" at offset 4
    }

    #[test]
    fn split_codon_junction_stays_in_frame() {
        // Gene MKVF = ATG|A  AA|GTT TTT: a phase-1 junction splits Lys's
        // codon (AAA) 1+2 across a group II intron (GTGCG...AT). Slot
        // profiles cover full codons only (M; V F), so the codon-aligned
        // matches leave out the split codon's bases - joining them as-is
        // (the pre-0.2 behaviour) gave ATG GTT TTT = MVF, dropping K and,
        // in real genes, frame-shifting everything downstream.
        let mut seq = b"CCC".to_vec();
        seq.extend(b"ATGA"); // exon1: 3..7
        seq.extend(b"GTGCGTTTTTTTTTTTAT"); // intron: 7..25
        seq.extend(b"AAGTTTTT"); // exon2: 25..33
        seq.extend(b"CCC");
        let mut sequences = HashMap::new();
        sequences.insert("ctg1".to_string(), seq);
        let fasta_idx = FastaIndex { sequences };

        let template = GeneTemplate {
            gene: "test".to_string(),
            slots: vec![
                SlotTemplate {
                    profile: custom_profile("e1", &[&[(b'M', 10.0)]]),
                    min_self_score: -100.0,
                    lead: 0,
                    trail: 1,
                },
                SlotTemplate {
                    profile: custom_profile("e2", &[&[(b'V', 10.0)], &[(b'F', 10.0)]]),
                    min_self_score: -100.0,
                    lead: 2,
                    trail: 0,
                },
            ],
            whole_gene_profile: custom_profile(
                "whole",
                &[&[(b'M', 10.0)], &[(b'K', 10.0)], &[(b'V', 10.0)], &[(b'F', 10.0)]],
            ),
        };
        let fragments = vec![
            RawFragment { contig: "ctg1".to_string(), start: 3, end: 6, strand: '+' },
            RawFragment { contig: "ctg1".to_string(), start: 27, end: 33, strand: '+' },
        ];
        let outcome = reconstruct_gene(&fasta_idx, &fragments, &template, &ScanParams::default(), 0, 1000, 0).unwrap();

        assert_eq!(outcome.protein, "MKVF");
        assert_eq!((outcome.exons[0].start, outcome.exons[0].end), (3, 7));
        assert_eq!((outcome.exons[1].start, outcome.exons[1].end), (25, 33));
        assert_eq!((outcome.exons[0].phase, outcome.exons[1].phase), (0, 2));
    }

    #[test]
    fn internal_exon_does_not_absorb_a_following_edit_created_stop() {
        // exon1 = ATG AAA; the next "codon" (CGA, an edit-created TGA stop)
        // is intron sequence. scan_gene folds it in as a stop; for a
        // non-final exon that must be undone.
        let mut seq = b"CCC".to_vec();
        seq.extend(b"ATGAAA"); // exon1: 3..9
        seq.extend(b"CGATTTTTTTTTTTTTTT"); // intron: 9..27
        seq.extend(b"GTTTTT"); // exon2: 27..33
        seq.extend(b"CCC");
        let mut sequences = HashMap::new();
        sequences.insert("ctg1".to_string(), seq);
        let fasta_idx = FastaIndex { sequences };
        let template = GeneTemplate {
            gene: "test".to_string(),
            slots: vec![
                SlotTemplate {
                    profile: custom_profile("e1", &[&[(b'M', 10.0)], &[(b'K', 10.0)]]),
                    min_self_score: -100.0,
                    lead: 0,
                    trail: 0,
                },
                SlotTemplate {
                    profile: custom_profile("e2", &[&[(b'V', 10.0)], &[(b'F', 10.0)]]),
                    min_self_score: -100.0,
                    lead: 0,
                    trail: 0,
                },
            ],
            whole_gene_profile: custom_profile(
                "whole",
                &[&[(b'M', 10.0)], &[(b'K', 10.0)], &[(b'V', 10.0)], &[(b'F', 10.0)]],
            ),
        };
        let fragments = vec![
            RawFragment { contig: "ctg1".to_string(), start: 3, end: 9, strand: '+' },
            RawFragment { contig: "ctg1".to_string(), start: 27, end: 33, strand: '+' },
        ];
        let outcome = reconstruct_gene(&fasta_idx, &fragments, &template, &ScanParams::default(), 0, 1000, 0).unwrap();

        assert_eq!(outcome.protein, "MKVF");
        assert_eq!((outcome.exons[0].start, outcome.exons[0].end), (3, 9));
        assert!(outcome.exons[0].edits.iter().all(|e| e.resulting_aa != '*'));
    }

    #[test]
    fn overlapping_candidates_never_fill_two_slots() {
        // One exon ATG AAA GTT TTT (MKVF). Slot 2's short profile (VF) also
        // matches *inside* it, via a second fragment covering its 3' half.
        // Both slots must not claim the same bases.
        let mut seq = b"CCC".to_vec();
        seq.extend(b"ATGAAAGTTTTT"); // 3..15
        seq.extend(b"CCCCCCCCCCCC");
        let mut sequences = HashMap::new();
        sequences.insert("ctg1".to_string(), seq);
        let fasta_idx = FastaIndex { sequences };
        let template = GeneTemplate {
            gene: "test".to_string(),
            slots: vec![
                SlotTemplate {
                    profile: custom_profile("e1", &[&[(b'M', 10.0)], &[(b'K', 10.0)], &[(b'V', 10.0)], &[(b'F', 10.0)]]),
                    min_self_score: -100.0, lead: 0, trail: 0,
                },
                SlotTemplate {
                    profile: custom_profile("e2", &[&[(b'V', 10.0)], &[(b'F', 10.0)]]),
                    min_self_score: -100.0, lead: 0, trail: 0,
                },
            ],
            whole_gene_profile: custom_profile("whole", &[&[(b'M', 10.0)], &[(b'K', 10.0)], &[(b'V', 10.0)], &[(b'F', 10.0)]]),
        };
        let fragments = vec![
            RawFragment { contig: "ctg1".to_string(), start: 3, end: 15, strand: '+' },
            RawFragment { contig: "ctg1".to_string(), start: 9, end: 15, strand: '+' },
        ];
        let outcome = reconstruct_gene(&fasta_idx, &fragments, &template, &ScanParams::default(), 0, 1000, 0).unwrap();
        for (i, a) in outcome.exons.iter().enumerate() {
            for b in &outcome.exons[i + 1..] {
                assert!(!(a.contig == b.contig && a.start < b.end && b.start < a.end), "exons overlap: {:?}-{:?}",
                        (a.start, a.end), (b.start, b.end));
            }
        }
        assert_eq!(outcome.protein, "MKVF");
    }

    fn custom_profile(gene: &str, columns: &[&[(u8, f64)]]) -> Profile {
        use orfedit::codon::AA_ORDER;
        let cols = columns
            .iter()
            .map(|specs| {
                let mut col = [-5.0f64; 20];
                for &(aa, score) in specs.iter() {
                    col[AA_ORDER.iter().position(|&x| x == aa).unwrap()] = score;
                }
                col
            })
            .collect();
        Profile { gene: gene.to_string(), columns: cols }
    }

    #[test]
    fn whole_gene_score_resolves_a_real_permutation_ambiguity() {
        // Mirrors the real nad2 misassignment found via ground-truth
        // validation (Arabidopsis_thaliana): 3 candidates (M/K/V codons,
        // well separated on one contig) and 3 slots whose profiles overlap
        // enough that the naive per-slot-greedy assignment (M->0, K->1,
        // V->2, local score 10+10+10=30) beats the biologically-correct
        // one (M->0, V->1, K->2, local score 10+9+9=28) on local scores
        // alone - but the correct one is what the whole-gene profile
        // (which wants M-V-K, not M-K-V) actually recognises. The fix
        // must pick the correct one despite its lower local-only score.
        let seq = b"TTTATGCCCNNNNNNNNNNNTTTAAACCCNNNNNNNNNNNTTTGTTCCC".to_vec();
        assert_eq!(seq.len(), 49);
        let mut sequences = HashMap::new();
        sequences.insert("ctg1".to_string(), seq);
        let fasta_idx = FastaIndex { sequences };

        let slot0 = custom_profile("slot0", &[&[(b'M', 10.0), (b'K', 9.0), (b'V', 1.0)]]);
        let slot1 = custom_profile("slot1", &[&[(b'K', 10.0), (b'V', 9.0), (b'M', 1.0)]]);
        let slot2 = custom_profile("slot2", &[&[(b'V', 10.0), (b'K', 9.0), (b'M', 1.0)]]);
        let whole = custom_profile("whole", &[&[(b'M', 10.0)], &[(b'V', 10.0)], &[(b'K', 10.0)]]);

        let template = GeneTemplate {
            gene: "test".to_string(),
            slots: vec![
                SlotTemplate { profile: slot0, min_self_score: -100.0, lead: 0, trail: 0 },
                SlotTemplate { profile: slot1, min_self_score: -100.0, lead: 0, trail: 0 },
                SlotTemplate { profile: slot2, min_self_score: -100.0, lead: 0, trail: 0 },
            ],
            whole_gene_profile: whole,
        };
        let fragments = vec![
            RawFragment { contig: "ctg1".to_string(), start: 0, end: 9, strand: '+' },   // M
            RawFragment { contig: "ctg1".to_string(), start: 20, end: 29, strand: '+' }, // K
            RawFragment { contig: "ctg1".to_string(), start: 40, end: 49, strand: '+' }, // V
        ];
        let params = ScanParams::default();

        let outcome = reconstruct_gene(&fasta_idx, &fragments, &template, &params, 0, 1000, 0).unwrap();

        assert_eq!(outcome.n_filled, 3);
        assert_eq!(outcome.protein, "MVK", "must pick the whole-gene-coherent permutation, not the naive per-slot-greedy one (which would give MKV)");
    }
}
