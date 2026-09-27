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
    let beam = best_assignments_beam(&scores, &floors, 2);
    let mut best_outcome: Option<(Vec<ExonResult>, Vec<JunctionResult>, f64, usize, String)> = None;
    for candidate_assignment in &beam {
        let Ok(result) = evaluate_assignment(
            fasta_idx, &candidates, template, candidate_assignment, &exon_params, params,
            &denovo_candidate_indices, junction_distance_threshold,
        ) else {
            continue;
        };
        let better = best_outcome.as_ref().is_none_or(|(_, _, prev, _, _)| result.2 > *prev);
        if better {
            best_outcome = Some(result);
        }
    }
    let (exons, junctions, whole_gene_score, whole_gene_n_edits, protein) =
        best_outcome.context("no candidate assignment (not even the empty one) produced a result")?;

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
) -> Result<(Vec<ExonResult>, Vec<JunctionResult>, f64, usize, String)> {
    let mut exons: Vec<ExonResult> = Vec::new();
    for (si, slot) in template.slots.iter().enumerate() {
        let Some(ci) = assignment.slot_to_candidate[si] else { continue };
        let cand = &candidates[ci];
        let outcome = scan_gene(fasta_idx, &cand.contig, cand.start, cand.end, cand.strand, &slot.profile, exon_params)?;
        exons.push(ExonResult {
            slot: si + 1,
            provenance: cand.provenance.clone(),
            contig: cand.contig.clone(),
            start: outcome.corrected_start,
            end: outcome.corrected_end,
            strand: cand.strand,
            score: outcome.score,
            edits: outcome.edits,
            via_denovo: denovo_candidate_indices.contains(&ci),
        });
    }

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
            slots: vec![SlotTemplate { profile: mkv_profile("test_exon1"), min_self_score: -100.0 }],
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
                SlotTemplate { profile: slot0, min_self_score: -100.0 },
                SlotTemplate { profile: slot1, min_self_score: -100.0 },
                SlotTemplate { profile: slot2, min_self_score: -100.0 },
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
