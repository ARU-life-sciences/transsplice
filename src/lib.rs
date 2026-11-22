pub mod assign;
pub mod candidate;
pub mod template;

use anyhow::{Context, Result};
use assign::{assignments_from_shortlists, best_assignment, classify_junction, Assignment};
use candidate::{build_candidates, Candidate, RawFragment};
use orfedit::align::{NtColumnCall, align_nt_frame0, align_nt_frame0_score};
use orfedit::fasta::FastaIndex;
use orfedit::{frameshifts_supported, scan_gene, ScanParams};
use std::cell::RefCell;
use std::collections::HashMap;
use template::GeneTemplate;

/// What `evaluate_assignment` needs from one `scan_gene` call, kept so the
/// same (candidate, slot) pair is scanned once per gene however many beam
/// assignments use it.
#[derive(Clone)]
struct ScanSummary {
    score: f64,
    corrected_start: usize,
    corrected_end: usize,
    includes_stop: bool,
    edits: Vec<(usize, usize, char, bool)>, // genomic_pos, profile_column, resulting_aa, reverse
}

impl ScanSummary {
    fn from_outcome(o: orfedit::ScanOutcome) -> Self {
        ScanSummary {
            score: o.score,
            corrected_start: o.corrected_start,
            corrected_end: o.corrected_end,
            includes_stop: o.includes_stop,
            edits: o.edits.iter().map(|e| (e.genomic_pos, e.profile_column, e.resulting_aa, e.reverse)).collect(),
        }
    }
}

/// Per-gene memo tables: scans by (candidate index, slot index), always
/// with the exon-scan parameters; whole-gene alignment scores by the
/// joined exon sequence, always with the gene's frameshift penalty.
#[derive(Default)]
struct GeneCache {
    scans: RefCell<HashMap<(usize, usize), Option<ScanSummary>>>,
    whole_scores: RefCell<HashMap<Vec<u8>, f64>>,
}

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
/// How far (bases) junction refinement slides a junction either way; 0,
/// the default, for no slide. A setting of the whole run (`--junction-slide`).
/// Off by default because it is not a clear gain everywhere: at 6, 135 of
/// 144 junctions in 11 ferns agree with intron-spanning RNA reads against
/// 126, but held out on 24 angiosperm RefSeq genomes 425 of 477 exons are
/// placed right against 427 (two plastid exons moved off).
static JUNCTION_SLIDE: std::sync::atomic::AtomicI64 = std::sync::atomic::AtomicI64::new(0);

pub fn set_junction_slide(bases: i64) {
    JUNCTION_SLIDE.store(bases.max(0), std::sync::atomic::Ordering::Relaxed);
}
/// Whole-gene score a slid junction must gain to be taken.
const JUNCTION_SLIDE_MARGIN: f64 = 1.0;
pub const SPLICE_MOTIF_BONUS: f64 = 20.0;

/// An exon slot whose profile has at most this many columns (90 bp) is
/// "short": an HMM search rarely hits it and a chance match can clear its
/// floor, so whatever candidate fills it is always checked against a
/// whole-contig search and a search next to the neighbouring exons, the
/// whole-gene alignment (plus intron motifs) deciding between them.
pub const SHORT_SLOT_COLUMNS: usize = 30;

/// How far (bp) from a placed neighbouring exon a short exon is looked for
/// on the same strand (a cis intron; group II introns are 0.5-4 kb).
const NEARBY_WINDOW: i64 = 5000;

/// How many of a short slot's best whole-contig hits are tried.
const DENOVO_ALTERNATIVES: usize = 3;

/// How many different places for a short exon the assignment search
/// considers (it considers two candidates for any other slot).
const SHORT_SLOT_PLACES: usize = 3;

/// Gene ends. Profiles are aligned end to end, so where a species' protein
/// starts later or stops sooner than the references' the first exon began
/// upstream of the start codon and the last ran on past the stop codon (a
/// stop is just a skipped codon to the aligner). If the first aligned codon
/// isn't Met, the first exon is cut back to the first ATG within this many
/// profile columns that sits on a column where the references have Met...
const START_SEARCH_COLUMNS: usize = 15;
/// ...and the last exon ends at the first in-frame stop codon: one inside
/// it (provided at least half the profile is aligned before it), else one
/// within this many codons of its last aligned codon.
const STOP_SEARCH_CODONS: usize = 15;

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
fn splice_motifs(
    fasta_idx: &FastaIndex,
    (a_contig, a_strand, a_span): (&str, char, (i64, i64)),
    (b_contig, b_strand, b_span): (&str, char, (i64, i64)),
) -> (bool, bool) {
    let five = if a_strand == '-' {
        sense_bases(fasta_idx, a_contig, a_span.0 - 5, a_span.0, '-')
    } else {
        sense_bases(fasta_idx, a_contig, a_span.1, a_span.1 + 5, '+')
    };
    let three = if b_strand == '-' {
        sense_bases(fasta_idx, b_contig, b_span.1, b_span.1 + 2, '-')
    } else {
        sense_bases(fasta_idx, b_contig, b_span.0 - 2, b_span.0, '+')
    };
    (
        five.len() == 5 && five[0] == b'G' && &five[2..5] == b"GCG",
        three.len() == 2 && three[0] == b'A' && (three[1] == b'C' || three[1] == b'T'),
    )
}

/// The junction-refinement bonus: SPLICE_MOTIF_BONUS per intron end found.
fn splice_motif_bonus(fasta_idx: &FastaIndex, a: (&str, char, (i64, i64)), b: (&str, char, (i64, i64))) -> f64 {
    let (five, three) = splice_motifs(fasta_idx, a, b);
    (five as u8 + three as u8) as f64 * SPLICE_MOTIF_BONUS
}

/// What a junction's intron ends say about its two exons being real, as
/// log2 odds (the profile scores' unit) - for choosing between different
/// places for an exon, where SPLICE_MOTIF_BONUS (sized to settle a
/// junction's exact position) would let a chance `AY`, one site in eight,
/// outweigh the protein. 5' GNGCG follows 75% of mito and 61% of plastid
/// exons in the reference genomes against 0.4% of random sites; 3' AY is
/// taken as 85% against 12.5%.
fn splice_motif_evidence(fasta_idx: &FastaIndex, a: (&str, char, (i64, i64)), b: (&str, char, (i64, i64))) -> f64 {
    let (five, three) = splice_motifs(fasta_idx, a, b);
    (if five { 7.0 } else { -1.5 }) + (if three { 2.5 } else { -2.5 })
}

/// A change of reading frame inside the reconstructed CDS (see orfedit's
/// frameshifts): the 1 or 2 bases skipped, on exon `slot`, genomic 0-based
/// half-open [start, end), before whole-gene profile column `profile_column`.
#[derive(Debug, Clone, PartialEq)]
pub struct GeneFrameshift {
    pub slot: usize,
    pub contig: String,
    pub start: usize,
    pub end: usize,
    pub n_bases: usize,
    pub profile_column: usize,
    /// Genomic [start, end) the frame change could equally well sit in
    /// (orfedit's frameshift window), clipped to the exon.
    pub window_start: usize,
    pub window_end: usize,
}

/// A frame change this close (bases) to an exon end is taken as the
/// junction being placed a base or two off, not a frameshift in the exon.
const FRAMESHIFT_JUNCTION_MARGIN: usize = 3;

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
    pub frameshifts: Vec<GeneFrameshift>,
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
    let (prepared, tasks) = prepare_gene(fasta_idx, fragments, template, params, merge_distance, exon_flank);
    let scans: Vec<DenovoScan> = tasks.iter().map(|t| denovo_scan(fasta_idx, template, params, t)).collect();
    finish_gene(fasta_idx, prepared, &tasks, &scans, template, params, junction_distance_threshold, exon_flank)
}

/// One whole-contig de novo scan a gene needs: an empty slot's profile
/// against one strand of one contig.
pub struct DenovoTask {
    pub slot: usize, // 0-indexed
    pub contig: String,
    pub len: i64,
    pub strand: char,
}

impl DenovoTask {
    /// DP size (profile columns x bases): what the scan's cost scales with.
    pub fn cells(&self, template: &GeneTemplate) -> u64 {
        template.slots[self.slot].profile.len() as u64 * self.len as u64
    }
}

/// A de novo scan's score and corrected [start, end), or None if it failed.
pub type DenovoScan = Option<(f64, i64, i64)>;

/// A gene after candidate scoring and first assignment, waiting for its
/// de novo scans (`finish_gene`).
pub struct Prepared {
    candidates: Vec<Candidate>,
    scores: Vec<Vec<f64>>,
    assignment: Assignment,
    cache: GeneCache,
}

/// Stage 1 of `reconstruct_gene`: build and score the candidates, assign
/// them, and list the de novo scans the empty slots need (in slot, contig,
/// strand order). The three stages are separate so a batch can run every
/// gene's de novo scans - nearly all of the run time, and very uneven
/// between genes - as one flat pool of parallel tasks.
pub fn prepare_gene(
    fasta_idx: &FastaIndex,
    fragments: &[RawFragment],
    template: &GeneTemplate,
    params: &ScanParams,
    merge_distance: i64,
    exon_flank: i64,
) -> (Prepared, Vec<DenovoTask>) {
    let k = template.slots.len();
    let candidates = build_candidates(fragments, merge_distance);

    let exon_params = ScanParams { flank: exon_flank, ..*params };
    let cache = GeneCache::default();
    let scores: Vec<Vec<f64>> = candidates
        .iter()
        .enumerate()
        .map(|(ci, cand)| {
            template
                .slots
                .iter()
                .enumerate()
                .map(|(si, slot)| {
                    let summary =
                        scan_gene(fasta_idx, &cand.contig, cand.start, cand.end, cand.strand, &slot.profile, &exon_params)
                            .ok()
                            .map(ScanSummary::from_outcome);
                    let score = summary.as_ref().map_or(f64::NEG_INFINITY, |o| o.score);
                    cache.scans.borrow_mut().insert((ci, si), summary);
                    score
                })
                .collect()
        })
        .collect();

    let floors: Vec<f64> = template.slots.iter().map(|s| s.min_self_score).collect();
    let assignment = best_assignment(&scores, &floors);

    // De novo: any slot no candidate could fill, and every short slot
    // (SHORT_SLOT_COLUMNS) whether filled or not, gets a whole-contig scan
    // (flank=0 => the whole requested region is the search window),
    // restricted to contigs already known to carry this gene in this
    // species - not a genome-wide scan.
    let mut contigs: Vec<String> = fragments.iter().map(|f| f.contig.clone()).collect();
    contigs.sort();
    contigs.dedup();
    let mut tasks: Vec<DenovoTask> = Vec::new();
    for si in (0..k).filter(|&si| assignment.slot_to_candidate[si].is_none() || is_short(&template.slots[si])) {
        for contig in &contigs {
            let Some(len) = fasta_idx.sequences.get(contig).map(|s| s.len() as i64) else { continue };
            for strand in ['+', '-'] {
                tasks.push(DenovoTask { slot: si, contig: contig.clone(), len, strand });
            }
        }
    }
    (Prepared { candidates, scores, assignment, cache }, tasks)
}

fn is_short(slot: &template::SlotTemplate) -> bool {
    slot.profile.len() <= SHORT_SLOT_COLUMNS
}

/// Stage 2: one de novo scan.
pub fn denovo_scan(fasta_idx: &FastaIndex, template: &GeneTemplate, params: &ScanParams, task: &DenovoTask) -> DenovoScan {
    let denovo_params = ScanParams { flank: 0, ..*params };
    scan_gene(fasta_idx, &task.contig, 0, task.len, task.strand, &template.slots[task.slot].profile, &denovo_params)
        .ok()
        .map(|o| (o.score, o.corrected_start as i64, o.corrected_end as i64))
}

/// Stage 3: take the de novo scans (`scans[i]` for `tasks[i]`, as
/// `prepare_gene` listed them, so ties resolve as in a serial loop), then
/// choose among assignments by whole-gene alignment and refine junctions.
#[allow(clippy::too_many_arguments)]
pub fn finish_gene(
    fasta_idx: &FastaIndex,
    prepared: Prepared,
    tasks: &[DenovoTask],
    scans: &[DenovoScan],
    template: &GeneTemplate,
    params: &ScanParams,
    junction_distance_threshold: i64,
    exon_flank: i64,
) -> Result<ReconstructionOutcome> {
    let k = template.slots.len();
    let Prepared { mut candidates, mut scores, mut assignment, cache } = prepared;
    let exon_params = ScanParams { flank: exon_flank, ..*params };
    let floors: Vec<f64> = template.slots.iter().map(|s| s.min_self_score).collect();
    let mut denovo_candidate_indices: std::collections::HashSet<usize> = std::collections::HashSet::new();

    // A slot's whole-contig hits at or above its floor, best first (ties
    // in task order, as a serial loop keeping the first best would).
    let hits_for = |si: usize| -> Vec<(f64, i64, i64, &DenovoTask)> {
        let mut hits: Vec<(f64, i64, i64, &DenovoTask)> = tasks
            .iter()
            .zip(scans)
            .filter(|(task, _)| task.slot == si)
            .filter_map(|(task, scan)| scan.map(|(score, start, end)| (score, start, end, task)))
            .filter(|h| h.0 >= floors[si])
            .collect();
        hits.sort_by(|x, y| y.0.partial_cmp(&x.0).unwrap_or(std::cmp::Ordering::Equal));
        hits
    };
    // Add a found exon as a candidate for slot `si` (or find it, if that
    // very span is one already); its index. A fragment whose own scan
    // lands on the same exon is not reused: it may be needed for another
    // slot (one HMM hit often spans two cis exons).
    #[allow(clippy::too_many_arguments)]
    fn add_candidate(
        candidates: &mut Vec<Candidate>,
        scores: &mut Vec<Vec<f64>>,
        found: &mut std::collections::HashSet<usize>,
        new: Candidate,
        si: usize,
        score: f64,
        k: usize,
    ) -> usize {
        let same = |c: &Candidate| c.contig == new.contig && c.strand == new.strand && c.start == new.start && c.end == new.end;
        if let Some(idx) = candidates.iter().position(same) {
            return idx;
        }
        let idx = candidates.len();
        candidates.push(new);
        for row in scores.iter_mut() {
            row.push(f64::NEG_INFINITY);
        }
        let mut row = vec![f64::NEG_INFINITY; k];
        row[si] = score;
        scores.push(row);
        found.insert(idx);
        idx
    }

    let empty_slots: Vec<usize> = (0..k).filter(|&si| assignment.slot_to_candidate[si].is_none()).collect();
    for &si in &empty_slots {
        if let Some(&(score, start, end, task)) = hits_for(si).first() {
            let new = Candidate { contig: task.contig.clone(), start, end, strand: task.strand, provenance: "denovo".to_string() };
            let idx = add_candidate(&mut candidates, &mut scores, &mut denovo_candidate_indices, new, si, score, k);
            assignment.slot_to_candidate[si] = Some(idx);
            assignment.total_score += score;
        }
    }
    // A short slot's best whole-contig hits join the candidates whether or
    // not something fills it, so the beam weighs them together with the
    // other slots: with a wrong exon in a short slot (or none) its
    // neighbours can be out of frame, and are then left out too.
    let short_slots: Vec<usize> = (0..k).filter(|&si| is_short(&template.slots[si])).collect();
    for &si in &short_slots {
        for (score, start, end, task) in hits_for(si).into_iter().take(DENOVO_ALTERNATIVES) {
            let new = Candidate { contig: task.contig.clone(), start, end, strand: task.strand, provenance: "denovo".to_string() };
            add_candidate(&mut candidates, &mut scores, &mut denovo_candidate_indices, new, si, score, k);
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
    // Per slot, its best two candidates at or above the floor; for a short
    // slot, its candidates at the best three places - a whole-contig hit
    // is often the very exon a fragment's own scan lands on, and the two
    // must not crowd a different place out, nor can one stand in for the
    // other (the fragment may be needed for the next slot).
    let shortlists: Vec<Vec<(usize, f64)>> = (0..k)
        .map(|si| {
            let mut ranked: Vec<(usize, f64)> =
                (0..candidates.len()).map(|ci| (ci, scores[ci][si])).filter(|c| c.1 >= floors[si]).collect();
            ranked.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
            if !short_slots.contains(&si) {
                ranked.truncate(2);
                return ranked;
            }
            let mut places: Vec<Option<(&str, char, usize, usize)>> = Vec::new();
            ranked.retain(|&(ci, _)| {
                let c = &candidates[ci];
                let place = cached_scan(fasta_idx, &cache, c, ci, si, &template.slots[si].profile, &exon_params)
                    .map(|o| (c.contig.as_str(), c.strand, o.corrected_start, o.corrected_end));
                if places.contains(&place) {
                    return true;
                }
                if places.len() == SHORT_SLOT_PLACES {
                    return false;
                }
                places.push(place);
                true
            });
            ranked
        })
        .collect();
    let beam = assignments_from_shortlists(&shortlists, candidates.len());
    let mut best: Option<(usize, Evaluated)> = None;
    for (bi, candidate_assignment) in beam.iter().enumerate() {
        let Ok(result) = evaluate_assignment(
            fasta_idx, &candidates, template, candidate_assignment, &exon_params, params,
            &denovo_candidate_indices, junction_distance_threshold, false, &cache,
        ) else {
            continue;
        };
        if best.as_ref().is_none_or(|(_, prev)| result.whole_gene_score > prev.whole_gene_score) {
            best = Some((bi, result));
        }
    }
    let (bi, unrefined) = best.context("no candidate assignment (not even the empty one) produced a result")?;
    let mut chosen = beam[bi].clone();
    let mut refined = evaluate_assignment(
        fasta_idx, &candidates, template, &chosen, &exon_params, params,
        &denovo_candidate_indices, junction_distance_threshold, true, &cache,
    )
    .ok();

    // Short slots: whatever fills one (or nothing), try its best
    // whole-contig hits and the best match next to each neighbouring exon
    // in its place, and keep a replacement only if the refined whole-gene
    // alignment plus intron motifs scores higher with it. A short profile's
    // own score can't tell the exon from a chance match; the protein it
    // completes and the intron ends beside it can.
    let denovo_params = ScanParams { flank: 0, ..*params };
    for &si in &short_slots {
        let mut alternatives: Vec<usize> = Vec::new();
        for (score, start, end, task) in hits_for(si).into_iter().take(DENOVO_ALTERNATIVES) {
            let new = Candidate { contig: task.contig.clone(), start, end, strand: task.strand, provenance: "denovo".to_string() };
            alternatives.push(add_candidate(&mut candidates, &mut scores, &mut denovo_candidate_indices, new, si, score, k));
        }
        if let Some(current) = &refined {
            for (contig, strand, start, end) in nearby_windows(&current.exons, si + 1, 3 * template.slots[si].profile.len() as i64) {
                let Ok(o) = scan_gene(fasta_idx, &contig, start, end, strand, &template.slots[si].profile, &denovo_params) else {
                    continue;
                };
                if o.score >= floors[si] {
                    let new = Candidate {
                        contig,
                        start: o.corrected_start as i64,
                        end: o.corrected_end as i64,
                        strand,
                        provenance: "nearby".to_string(),
                    };
                    alternatives.push(add_candidate(
                        &mut candidates, &mut scores, &mut denovo_candidate_indices, new, si, o.score, k,
                    ));
                }
            }
        }
        let mut tried: Vec<(String, char, usize, usize)> = Vec::new();
        let locus = |ci: usize| -> Option<(String, char, usize, usize)> {
            let c = &candidates[ci];
            cached_scan(fasta_idx, &cache, c, ci, si, &template.slots[si].profile, &exon_params)
                .map(|o| (c.contig.clone(), c.strand, o.corrected_start, o.corrected_end))
        };
        if let Some(l) = chosen.slot_to_candidate[si].and_then(&locus) {
            tried.push(l);
        }
        for ci in alternatives {
            if chosen.slot_to_candidate.contains(&Some(ci)) {
                continue;
            }
            // the same stretch of genome as one already tried: same result
            let Some(l) = locus(ci) else { continue };
            if tried.contains(&l) {
                continue;
            }
            tried.push(l);
            let mut alternative = chosen.clone();
            alternative.slot_to_candidate[si] = Some(ci);
            let Ok(result) = evaluate_assignment(
                fasta_idx, &candidates, template, &alternative, &exon_params, params,
                &denovo_candidate_indices, junction_distance_threshold, true, &cache,
            ) else {
                continue;
            };
            if refined.as_ref().is_none_or(|r| result.total() > r.total()) {
                chosen = alternative;
                refined = Some(result);
            }
        }
    }
    let Evaluated { exons, junctions, whole_gene_score, whole_gene_n_edits, protein, frameshifts, .. } =
        refined.unwrap_or(unrefined);

    Ok(ReconstructionOutcome {
        n_slots: k,
        n_filled: exons.len(),
        complete: exons.len() == k,
        exons,
        junctions,
        whole_gene_score,
        whole_gene_n_edits,
        protein,
        frameshifts,
    })
}

/// Where to look for exon `slot` (1-indexed) beside its placed neighbours:
/// up to NEARBY_WINDOW after the previous exon and before the next one, on
/// their strand, not reaching past the other neighbour when both lie on
/// one strand of one contig. (contig, strand, start, end); windows shorter
/// than `min_len` are left out.
fn nearby_windows(exons: &[ExonResult], slot: usize, min_len: i64) -> Vec<(String, char, i64, i64)> {
    let prev = exons.iter().find(|e| e.slot + 1 == slot);
    let next = exons.iter().find(|e| e.slot == slot + 1);
    let mut windows: Vec<(String, char, i64, i64)> = Vec::new();
    if let Some(p) = prev {
        let (ps, pe) = (p.start as i64, p.end as i64);
        let (mut s, mut e) = if p.strand == '-' { (ps - NEARBY_WINDOW, ps) } else { (pe, pe + NEARBY_WINDOW) };
        if let Some(n) = next.filter(|n| n.contig == p.contig && n.strand == p.strand) {
            if p.strand == '-' && (n.end as i64) <= ps {
                s = s.max(n.end as i64);
            } else if p.strand != '-' && (n.start as i64) >= pe {
                e = e.min(n.start as i64);
            }
        }
        windows.push((p.contig.clone(), p.strand, s.max(0), e));
    }
    if let Some(n) = next {
        let (ns, ne) = (n.start as i64, n.end as i64);
        let (mut s, mut e) = if n.strand == '-' { (ne, ne + NEARBY_WINDOW) } else { (ns - NEARBY_WINDOW, ns) };
        if let Some(p) = prev.filter(|p| p.contig == n.contig && p.strand == n.strand) {
            if n.strand == '-' && (p.start as i64) >= ne {
                e = e.min(p.start as i64);
            } else if n.strand != '-' && (p.end as i64) <= ns {
                s = s.max(p.end as i64);
            }
        }
        windows.push((n.contig.clone(), n.strand, s.max(0), e));
    }
    windows.retain(|w| w.3 - w.2 >= min_len);
    windows.dedup();
    windows
}

/// `scan_gene` of candidate `ci` against slot `si`'s profile with the
/// exon-scan parameters, done once per gene. Same contig/coordinates/
/// strand/profile/parameters => same result, so the scan done for the
/// score matrix (or for an earlier assignment, for candidates added after
/// it) is reused.
fn cached_scan(
    fasta_idx: &FastaIndex,
    cache: &GeneCache,
    cand: &Candidate,
    ci: usize,
    si: usize,
    profile: &orfedit::profile::Profile,
    exon_params: &ScanParams,
) -> Option<ScanSummary> {
    if let Some(o) = cache.scans.borrow().get(&(ci, si)) {
        return o.clone();
    }
    let o = scan_gene(fasta_idx, &cand.contig, cand.start, cand.end, cand.strand, profile, exon_params)
        .ok()
        .map(ScanSummary::from_outcome);
    cache.scans.borrow_mut().insert((ci, si), o.clone());
    o
}

/// One assignment's exons, junctions and whole-gene alignment.
struct Evaluated {
    exons: Vec<ExonResult>,
    junctions: Vec<JunctionResult>,
    whole_gene_score: f64,
    whole_gene_n_edits: usize,
    protein: String,
    frameshifts: Vec<GeneFrameshift>,
    /// Intron-motif evidence (splice_motif_evidence) summed over the
    /// junctions between adjacent slots.
    motif_bonus: f64,
}

impl Evaluated {
    /// What alternatives for a short slot are compared by.
    fn total(&self) -> f64 {
        self.whole_gene_score + self.motif_bonus
    }
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
    cache: &GeneCache,
) -> Result<Evaluated> {
    let k = template.slots.len();
    let mut placed: Vec<Placed> = Vec::new();
    for (si, slot) in template.slots.iter().enumerate() {
        let Some(ci) = assignment.slot_to_candidate[si] else { continue };
        let cand = &candidates[ci];
        let outcome =
            cached_scan(fasta_idx, cache, cand, ci, si, &slot.profile, exon_params).context("exon scan failed")?;
        let (mut core_start, mut core_end) = (outcome.corrected_start as i64, outcome.corrected_end as i64);
        let mut edits: Vec<orfedit::EditCall> = outcome
            .edits
            .iter()
            .map(|&(genomic_pos, profile_column, resulting_aa, reverse)| orfedit::EditCall { genomic_pos, profile_column, resulting_aa, reverse })
            .collect();
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
    // Best first, so an exon that is itself dropped doesn't knock out a
    // third one it overlaps.
    let mut order: Vec<usize> = (0..placed.len()).collect();
    order.sort_by(|&i, &j| placed[j].score.partial_cmp(&placed[i].score).unwrap_or(std::cmp::Ordering::Equal).then(i.cmp(&j)));
    let mut kept = vec![false; placed.len()];
    for &i in &order {
        let a = &placed[i];
        let clash = (0..placed.len()).any(|j| {
            let b = &placed[j];
            kept[j] && a.contig == b.contig && a.core_start < b.core_end && b.core_start < a.core_end
        });
        kept[i] = !clash;
    }
    let placed: Vec<Placed> = placed.into_iter().zip(kept).filter(|(_, keep)| *keep).map(|(p, _)| p).collect();

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
        // Coordinate descent revisits placements (the chosen split is
        // re-scored as the zero shift, and so on): score each distinct
        // joined sequence once.
        let concat = concat_of(ext5, ext3);
        if let Some(&s) = cache.whole_scores.borrow().get(&concat) {
            return s;
        }
        let s = align_nt_frame0_score(
            &template.whole_gene_profile,
            &concat,
            params.gap_open,
            params.gap_extend,
            params.edit_costs(),
            params.frameshift_penalty,
        );
        cache.whole_scores.borrow_mut().insert(concat, s);
        s
    };
    let junction_bonus = |i: usize, ext5: &[i64], ext3: &[i64]| -> f64 {
        let (a, b) = (&placed[i], &placed[i + 1]);
        splice_motif_bonus(
            fasta_idx,
            (&a.contig, a.strand, a.span(ext5[i], ext3[i])),
            (&b.contig, b.strand, b.span(ext5[i + 1], ext3[i + 1])),
        )
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
        let eval_ext = |e3: i64, e5: i64| -> Option<f64> {
            if e3 < floor_a || e5 < floor_b {
                return None;
            }
            let (mut t5, mut t3) = (ext5.clone(), ext3.clone());
            t3[i] = e3;
            t5[i + 1] = e5;
            Some(score_of(&t5, &t3) + junction_bonus(i, &t5, &t3))
        };
        let eval = |t: i64, da: i64, db: i64| -> Option<f64> { eval_ext(t + 3 * da, (3 - t).rem_euclid(3) + 3 * db) };
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
        if let Some(found) = eval(t, da, db) {
            let (e3, e5) = (t + 3 * da, (3 - t).rem_euclid(3) + 3 * db);
            // Then slide the junction along: k bases off one exon and onto
            // the other, the joined length (and the frame) unchanged. The
            // descent above cannot get there - a base moved across is a
            // change of split and of both exons' lengths at once, and each
            // of those alone breaks the frame. Seen against intron-spanning
            // RNA reads in 11 ferns: nad5's exons 3-4 (a 22 bp exon) sat 1-2
            // bases along from where the reads put them in all ten.
            // Judged on the protein alone, and only for a clear gain: the
            // intron-end bonus, sized to settle a split codon, moved 13
            // junctions the reads confirm onto a chance motif nearby. (A
            // gain of 1 to 4 gives the same result; at 8 nothing moves.)
            let _ = found;
            let protein_of = |e3: i64, e5: i64| -> Option<f64> {
                if e3 < floor_a || e5 < floor_b {
                    return None;
                }
                let (mut t5, mut t3) = (ext5.clone(), ext3.clone());
                t3[i] = e3;
                t5[i + 1] = e5;
                Some(score_of(&t5, &t3))
            };
            let base = protein_of(e3, e5).unwrap_or(f64::NEG_INFINITY);
            let reach = JUNCTION_SLIDE.load(std::sync::atomic::Ordering::Relaxed);
            let (mut best, mut slide) = (base + JUNCTION_SLIDE_MARGIN, 0);
            for k in (-reach..=reach).filter(|&k| k != 0) {
                if let Some(s) = protein_of(e3 + k, e5 - k) {
                    if s > best {
                        (best, slide) = (s, k);
                    }
                }
            }
            ext3[i] = e3 + slide;
            ext5[i + 1] = e5 - slide;
        }
    }

    let mut exons: Vec<ExonResult> = placed
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

    // Whole-gene pass: concatenate the corrected exon sequences (in slot
    // order, each already strand-corrected/frame-trimmed by its own
    // scan_gene call above) and align once against the whole-gene profile -
    // this is the arbiter the beam search picks an assignment by. Frame
    // changes are allowed (params.frameshift_penalty) and kept only with
    // sound alignment either side (orfedit::frameshifts_supported), as in
    // orfedit's own single-exon scan.
    let align_exons = |exons: &[ExonResult]| -> Result<(Vec<u8>, Vec<(usize, usize)>, orfedit::align::NtAlignResult)> {
        let mut concat: Vec<u8> = Vec::new();
        let mut offsets = Vec::with_capacity(exons.len()); // (concat start, exon index)
        for (i, e) in exons.iter().enumerate() {
            let (_, _, seq) = fasta_idx.extract_window(&e.contig, e.start as i64, e.end as i64, e.strand)?;
            offsets.push((concat.len(), i));
            concat.extend(seq);
        }
        let align_whole = |penalty: f64| {
            align_nt_frame0(&template.whole_gene_profile, &concat, params.gap_open, params.gap_extend, params.edit_costs(), penalty)
        };
        let mut whole = align_whole(params.frameshift_penalty);
        cache.whole_scores.borrow_mut().entry(concat.clone()).or_insert(whole.score);
        if !whole.frameshifts.is_empty() && !frameshifts_supported(&template.whole_gene_profile, &whole, params) {
            whole = align_whole(f64::INFINITY);
        }
        Ok((concat, offsets, whole))
    };
    let (concat, mut offsets, mut whole) = align_exons(&exons)?;
    // Gene ends (the chosen assignment only): start the first exon at the
    // start codon and end the last at the stop codon, then align again.
    if refine_junctions && !exons.is_empty() {
        let (start_trim, end_change) = terminal_adjustment(fasta_idx, template, &exons, &concat, &offsets, &whole, params.edit_costs().reverse());
        let mut changed = false;
        if start_trim > 0 {
            let e = &mut exons[0];
            if e.end - e.start > start_trim {
                if e.strand == '-' { e.end -= start_trim } else { e.start += start_trim }
                changed = true;
            }
        }
        if end_change != 0 {
            let e = exons.last_mut().unwrap();
            let (start, end) = if e.strand == '-' {
                (e.start as i64 - end_change, e.end as i64)
            } else {
                (e.start as i64, e.end as i64 + end_change)
            };
            let contig_len = fasta_idx.sequences.get(&e.contig).map_or(0, |s| s.len() as i64);
            if start >= 0 && end <= contig_len && end - start >= 3 {
                (e.start, e.end) = (start as usize, end as usize);
                // the exon now ends at a literal stop, not one made by an edit
                e.edits.retain(|x| x.resulting_aa != '*');
                changed = true;
            }
        }
        if changed {
            for e in exons.iter_mut() {
                let (start, end) = (e.start, e.end);
                e.edits.retain(|x| x.genomic_pos >= start && x.genomic_pos < end);
            }
            (_, offsets, whole) = align_exons(&exons)?;
        }
    }

    let mut junctions = Vec::new();
    let mut motif_bonus = 0.0;
    for w in exons.windows(2) {
        let (a, b) = (&w[0], &w[1]);
        let (kind, gap) = classify_junction(
            &a.contig, a.start as i64, a.end as i64, a.strand,
            &b.contig, b.start as i64, b.end as i64, b.strand,
            junction_distance_threshold,
        );
        junctions.push(JunctionResult { from_slot: a.slot, to_slot: b.slot, kind: kind.as_str(), gap_bp: gap });
        if b.slot == a.slot + 1 {
            motif_bonus += splice_motif_evidence(
                fasta_idx,
                (&a.contig, a.strand, (a.start as i64, a.end as i64)),
                (&b.contig, b.strand, (b.start as i64, b.end as i64)),
            );
        }
    }

    let mut protein = String::with_capacity(whole.columns.len());
    let mut whole_gene_n_edits = 0usize;
    for call in &whole.columns {
        match call {
            NtColumnCall::Aligned { aa, n_edits, .. } => {
                protein.push(*aa as char);
                if *n_edits > 0 {
                    whole_gene_n_edits += 1;
                }
            }
            NtColumnCall::Deleted => protein.push('-'),
        }
    }
    let mut frameshifts = Vec::new();
    for f in &whole.frameshifts {
        // the exon holding the skipped bases, and their offset within it
        let Some(&(off, i)) = offsets.iter().rev().find(|(o, _)| *o <= f.sense_pos) else { continue };
        let e = &exons[i];
        let local = f.sense_pos - off;
        let len = e.end - e.start;
        if local < FRAMESHIFT_JUNCTION_MARGIN || local + f.n_bases + FRAMESHIFT_JUNCTION_MARGIN > len {
            continue;
        }
        // concat offsets [a, b) within this exon -> genomic [start, end)
        let to_genomic = |a: usize, b: usize| -> (usize, usize) {
            if e.strand == '-' { (e.end - b, e.end - a) } else { (e.start + a, e.start + b) }
        };
        let (start, end) = to_genomic(local, local + f.n_bases);
        let (wa, wb) = (f.window.0.max(off) - off, (f.window.1 - off).min(len));
        let (window_start, window_end) = to_genomic(wa.min(local), wb.max(local + f.n_bases));
        frameshifts.push(GeneFrameshift {
            slot: e.slot,
            contig: e.contig.clone(),
            start,
            end,
            n_bases: f.n_bases,
            profile_column: f.after_column,
            window_start,
            window_end,
        });
    }

    Ok(Evaluated { exons, junctions, whole_gene_score: whole.score, whole_gene_n_edits, protein, frameshifts, motif_bonus })
}

/// How far to move the gene's two ends, from the whole-gene alignment
/// `whole` of the joined exons `concat` (see START_SEARCH_COLUMNS):
/// bases to cut from the first exon's 5' end, and bases to add to (or,
/// negative, cut from) the last exon's 3' end. (0, 0) when the first/last
/// slot isn't placed or nothing needs moving.
fn terminal_adjustment(
    fasta_idx: &FastaIndex,
    template: &GeneTemplate,
    exons: &[ExonResult],
    concat: &[u8],
    offsets: &[(usize, usize)],
    whole: &orfedit::align::NtAlignResult,
    reverse_editing: bool,
) -> (usize, i64) {
    let profile = &template.whole_gene_profile;
    let n_columns = profile.len();
    // (profile column, concat position of its codon, amino acid called)
    let aligned: Vec<(usize, usize, u8)> = whole
        .columns
        .iter()
        .enumerate()
        .filter_map(|(c, call)| match call {
            NtColumnCall::Aligned { sense_start, aa, .. } => Some((c, *sense_start, *aa)),
            NtColumnCall::Deleted => None,
        })
        .collect();
    let (Some(&(c0, s0, aa0)), Some(&(c_last, s_last, _))) = (aligned.first(), aligned.last()) else {
        return (0, 0);
    };
    let codon = |q: usize| -> &[u8] { &concat[q..q + 3] };
    let is_stop = |c: &[u8]| matches!(c, b"TAA" | b"TAG" | b"TGA");

    let mut start_trim = 0usize;
    if exons[0].slot == 1 {
        let first_len = offsets.get(1).map_or(concat.len(), |o| o.0);
        let met = orfedit::codon::aa_index(b'M');
        let at = if aa0 == b'M' {
            s0
        } else {
            aligned
                .iter()
                .take_while(|a| a.0 <= c0 + START_SEARCH_COLUMNS)
                .find(|a| codon(a.1) == b"ATG" && met.is_some_and(|m| profile.columns[a.0][m] > 0.0))
                .map_or(0, |a| a.1)
        };
        if at > 0 && at + 3 <= first_len && !whole.frameshifts.iter().any(|f| f.sense_pos < at) {
            start_trim = at;
        }
    }

    let mut end_change = 0i64;
    let last = exons.last().unwrap();
    let last_off = offsets.last().map_or(0, |o| o.0);
    if last.slot == template.slots.len() && s_last >= last_off {
        // a stop codon inside the last exon, among the codons the alignment
        // skipped between two aligned ones
        let mut interior: Option<(usize, usize)> = None; // (concat position, column aligned before it)
        'search: for w in aligned.windows(2) {
            let ((ca, sa, _), (_, sb, _)) = (w[0], w[1]);
            if sb < sa + 6 || (sb - sa) % 3 != 0 {
                continue;
            }
            for q in (sa + 3..sb).step_by(3) {
                if q >= last_off && is_stop(codon(q)) {
                    interior = Some((q, ca));
                    break 'search;
                }
            }
        }
        // With U->C editing a stop codon inside the gene may be a residue
        // the alignment chose not to read (the edit cost more than the
        // column was worth): it says nothing about where the gene ends.
        if reverse_editing {
            interior = None;
        }
        if let Some((q, column)) = interior {
            // an early stop is a pseudogene's (or needs reverse editing):
            // leave the exon as the alignment has it
            if (column + 1) * 2 >= n_columns {
                end_change = (q + 3) as i64 - concat.len() as i64;
            }
        } else if n_columns - 1 - c_last <= STOP_SEARCH_CODONS {
            // no stop inside: the first one within reach of the last
            // aligned codon, in the exon's unaligned tail or beyond it
            let tail = &concat[s_last + 3..];
            let n = 3 * STOP_SEARCH_CODONS as i64;
            let beyond = if last.strand == '-' {
                sense_bases(fasta_idx, &last.contig, last.start as i64 - n, last.start as i64, '-')
            } else {
                sense_bases(fasta_idx, &last.contig, last.end as i64, last.end as i64 + n, '+')
            };
            let mut seq = tail.to_vec();
            seq.extend(beyond);
            // A tail of one CAA/CAG/CGA is a stop made by editing (orfedit
            // folded it in) - unless a literal stop follows it directly:
            // then it is the last residue and that the stop.
            let edit_stop = tail.len() == 3 && matches!(tail, b"CAA" | b"CAG" | b"CGA");
            let reach = if edit_stop { 2 } else { STOP_SEARCH_CODONS };
            if let Some(i) = (0..reach).find(|i| 3 * i + 3 <= seq.len() && is_stop(&seq[3 * i..3 * i + 3])) {
                end_change = (3 * i + 3) as i64 - tail.len() as i64;
            }
        }
    }
    (start_trim, end_change)
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

    fn one_contig(seq: &str) -> FastaIndex {
        let mut sequences = HashMap::new();
        sequences.insert("ctg1".to_string(), seq.as_bytes().to_vec());
        FastaIndex { sequences }
    }

    fn single_exon_template(residues: &[u8]) -> GeneTemplate {
        let columns: Vec<Vec<(u8, f64)>> = residues.iter().map(|&aa| vec![(aa, 10.0)]).collect();
        let columns: Vec<&[(u8, f64)]> = columns.iter().map(|c| c.as_slice()).collect();
        GeneTemplate {
            gene: "test".to_string(),
            slots: vec![SlotTemplate { profile: custom_profile("e1", &columns), min_self_score: -100.0, lead: 0, trail: 0 }],
            whole_gene_profile: custom_profile("whole", &columns),
        }
    }

    #[test]
    fn last_exon_ends_at_a_stop_codon_inside_it() {
        // The references' protein is MKVF; this gene stops after MKV, and
        // the aligner, to which a stop is just a skipped codon, carries on
        // to the F beyond it.
        let fasta_idx = one_contig("CCCATGAAAGTTTAATTTCCC");
        let template = single_exon_template(b"MKVF");
        let fragments = vec![RawFragment { contig: "ctg1".to_string(), start: 3, end: 18, strand: '+' }];
        let outcome = reconstruct_gene(&fasta_idx, &fragments, &template, &ScanParams::default(), 0, 1000, 0).unwrap();
        let exon = &outcome.exons[0];
        assert_eq!((exon.start, exon.end), (3, 15), "ATG AAA GTT TAA");
        assert_eq!(outcome.protein, "MKV-");
    }

    #[test]
    fn last_exon_is_extended_to_a_stop_codon_a_few_codons_on() {
        // profile MKV, gene MKVGP*: the last two residues aren't in the
        // profile, so the match ends two codons short of the stop
        for (seq, start, end, strand, expect) in [
            ("CCCATGAAAGTTGGGCCCTAAGGG", 3, 12, '+', (3, 21)),
            ("CCCTTAGGGCCCAACTTTCATGGG", 12, 21, '-', (3, 21)),
        ] {
            let fasta_idx = one_contig(seq);
            let template = single_exon_template(b"MKV");
            let fragments = vec![RawFragment { contig: "ctg1".to_string(), start, end, strand }];
            let outcome = reconstruct_gene(&fasta_idx, &fragments, &template, &ScanParams::default(), 0, 1000, 0).unwrap();
            let exon = &outcome.exons[0];
            assert_eq!((exon.start, exon.end), expect, "strand {strand}");
            assert_eq!(outcome.protein, "MKV");
        }
    }

    #[test]
    fn first_exon_starts_at_the_start_codon() {
        // The references start M K M V F; this gene starts at the second
        // Met, and the profile's first two columns were aligned to the
        // codons upstream of it (CTT AAA).
        let fasta_idx = one_contig("CCCCTTAAAATGGTTTTTTAACCC");
        let template = single_exon_template(b"MKMVF");
        let fragments = vec![RawFragment { contig: "ctg1".to_string(), start: 3, end: 18, strand: '+' }];
        let outcome = reconstruct_gene(&fasta_idx, &fragments, &template, &ScanParams::default(), 0, 1000, 0).unwrap();
        let exon = &outcome.exons[0];
        assert_eq!((exon.start, exon.end), (9, 21), "ATG GTT TTT TAA");
    }

    #[test]
    fn short_exon_is_taken_from_beside_its_neighbour_not_from_a_chance_match() {
        // Gene MKVF|WW. The only candidate handed in for the two-codon
        // second exon is a chance TGGTGG far from exon 1; the real one
        // follows exon 1's intron (GTGCG...AT) and is found by the search
        // every short slot gets. Both score the same against the profile;
        // the intron's 3' end decides.
        let mut seq = String::from("CCC");
        seq.push_str("ATGAAAGTTTTT"); // exon 1: 3..15
        seq.push_str("GTGCGTTTTTTTAT"); // intron: 15..29
        seq.push_str("TGGTGG"); // exon 2: 29..35
        seq.push_str("CCCCCCCCCCCC"); // 35..47
        seq.push_str("TGGTGG"); // chance match: 47..53
        seq.push_str("CCC");
        let fasta_idx = one_contig(&seq);
        let cols = |s: &[u8]| -> Vec<Vec<(u8, f64)>> { s.iter().map(|&aa| vec![(aa, 10.0)]).collect() };
        let profile = |name: &str, s: &[u8]| {
            let c = cols(s);
            custom_profile(name, &c.iter().map(|c| c.as_slice()).collect::<Vec<_>>())
        };
        let template = GeneTemplate {
            gene: "test".to_string(),
            slots: vec![
                SlotTemplate { profile: profile("e1", b"MKVF"), min_self_score: -100.0, lead: 0, trail: 0 },
                SlotTemplate { profile: profile("e2", b"WW"), min_self_score: 10.0, lead: 0, trail: 0 },
            ],
            whole_gene_profile: profile("whole", b"MKVFWW"),
        };
        let fragments = vec![
            RawFragment { contig: "ctg1".to_string(), start: 3, end: 15, strand: '+' },
            RawFragment { contig: "ctg1".to_string(), start: 47, end: 53, strand: '+' },
        ];
        let outcome = reconstruct_gene(&fasta_idx, &fragments, &template, &ScanParams::default(), 0, 1000, 0).unwrap();
        assert_eq!(outcome.n_filled, 2);
        let exon = &outcome.exons[1];
        assert_eq!((exon.start, exon.end), (29, 35));
        assert!(exon.via_denovo);
        assert_eq!(outcome.protein, "MKVFWW");
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

    #[test]
    fn frameshift_in_an_exon_is_reported_on_its_slot() {
        // exon1 = (MKVWF)x3; group II intron; exon2 = (MKVWF)x6 + stop with
        // one extra A after its 15th codon
        let half = "ATGAAAGTTTGGTTT".repeat(3);
        let mut seq = b"CCC".to_vec();
        seq.extend(half.as_bytes()); // exon1: 3..48
        seq.extend(b"GTGCGTTTTTTTTTTTTTTTTTTTTAT"); // intron: 48..75
        seq.extend(half.as_bytes()); // exon2: 75..
        seq.push(b'A'); // 120
        seq.extend(half.as_bytes());
        seq.extend(b"TAACCC");
        let mut sequences = HashMap::new();
        sequences.insert("ctg1".to_string(), seq);
        let fasta_idx = FastaIndex { sequences };
        let cols = |n: usize| -> Vec<&'static [(u8, f64)]> {
            let one: [&'static [(u8, f64)]; 5] =
                [&[(b'M', 5.0)], &[(b'K', 5.0)], &[(b'V', 5.0)], &[(b'W', 5.0)], &[(b'F', 5.0)]];
            (0..n).map(|i| one[i % 5]).collect()
        };
        let template = GeneTemplate {
            gene: "test".to_string(),
            slots: vec![
                SlotTemplate { profile: custom_profile("e1", &cols(15)), min_self_score: -100.0, lead: 0, trail: 0 },
                SlotTemplate { profile: custom_profile("e2", &cols(30)), min_self_score: -100.0, lead: 0, trail: 0 },
            ],
            whole_gene_profile: custom_profile("whole", &cols(45)),
        };
        let fragments = vec![
            RawFragment { contig: "ctg1".to_string(), start: 3, end: 48, strand: '+' },
            RawFragment { contig: "ctg1".to_string(), start: 75, end: 169, strand: '+' },
        ];
        let params = ScanParams { frameshift_penalty: 10.0, ..ScanParams::default() };
        let outcome = reconstruct_gene(&fasta_idx, &fragments, &template, &params, 0, 1000, 0).unwrap();
        assert_eq!(outcome.protein, "MKVWF".repeat(9));
        assert_eq!(outcome.frameshifts.len(), 1);
        let f = &outcome.frameshifts[0];
        assert_eq!((f.slot, f.start, f.end, f.n_bases), (2, 120, 121, 1));

        // single-frame: no frameshift, and exon 2 is cut short
        let outcome = reconstruct_gene(&fasta_idx, &fragments, &template, &ScanParams::default(), 0, 1000, 0).unwrap();
        assert!(outcome.frameshifts.is_empty());
    }
}
