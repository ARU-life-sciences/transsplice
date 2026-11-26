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
use std::collections::{HashMap, HashSet};
use std::path::Path;
use template::GeneTemplate;

/// The circular contigs of a FASTA and their lengths (`wrap_circular`). A
/// gene's exons either side of such a contig's origin lie together (see
/// `classify_junction`), an exon is looked for next to its neighbour
/// across the origin, and an exon may itself cross it.
pub type Circular = HashMap<String, i64>;

/// Names of the circular contigs in a FASTA: oatk's headers say
/// `circular=true` (or false).
pub fn circular_contigs(path: &Path) -> Result<HashSet<String>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    Ok(text
        .lines()
        .filter_map(|l| l.strip_prefix('>'))
        .filter(|h| h.split_whitespace().skip(1).any(|w| w == "circular=true"))
        .filter_map(|h| h.split_whitespace().next().map(str::to_string))
        .collect())
}

/// How far (bp) past its end a circular contig's sequence is continued
/// from its start: room for an exon across the origin, the window it is
/// looked for in beside its neighbour (NEARBY_WINDOW), and two fragments
/// merged across it (`--merge-distance`).
pub const WRAP_BP: usize = 20_000;

/// Continue each named contig's sequence past its end with its first
/// WRAP_BP bases (all of it, if shorter), so that a window can run round
/// the origin: a position p >= length is position p - length. Returns the
/// contigs' own lengths, which everything that needs a contig's length
/// takes from here, not from the longer sequence.
pub fn wrap_circular(fasta_idx: &mut FastaIndex, names: &HashSet<String>) -> Circular {
    let mut out = Circular::new();
    for name in names {
        if let Some(seq) = fasta_idx.sequences.get_mut(name) {
            let n = seq.len();
            seq.extend_from_within(..WRAP_BP.min(n));
            out.insert(name.clone(), n as i64);
        }
    }
    out
}

/// A circular contig's own length.
fn circular_len(circular: &Circular, contig: &str) -> Option<i64> {
    circular.get(contig).copied()
}

/// Whether any exon runs past its circular contig's end, across the origin.
fn crosses_origin(exons: &[ExonResult], circular: &Circular) -> bool {
    exons.iter().any(|e| circular_len(circular, &e.contig).is_some_and(|n| e.end as i64 > n))
}

/// A contig's own length: a circular one's, else its sequence's.
fn contig_len(fasta_idx: &FastaIndex, circular: &Circular, contig: &str) -> Option<i64> {
    circular_len(circular, contig).or_else(|| fasta_idx.sequences.get(contig).map(|s| s.len() as i64))
}

/// [start, end) found past a circular contig's end (in its continuation,
/// `wrap_circular`) moved back by its length; one across the origin keeps
/// its end past the length.
fn canonical(circular: &Circular, contig: &str, start: i64, end: i64) -> (i64, i64) {
    match circular_len(circular, contig) {
        Some(n) if start >= n => (start - n, end - n),
        _ => (start, end),
    }
}

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

/// What a junction that is not cis - exons far apart or on opposite
/// strands of one contig - costs when an exon on it is weighed against one
/// next to its neighbour, or against none (the pass for exons of two
/// copies, after the beam), times the share of the template's references
/// in which that junction is cis (`junctions.tsv`; a template without it
/// costs nothing).
/// A setting of the whole run (`--trans-junction-penalty`); 0, the default,
/// for none. Picking each exon's best match anywhere joins a gene from two
/// or three diverged copies of it (plastid clpP in *Aquilegia*, without
/// introns, in two copies 16 kb apart), or takes a chance match of a short
/// exon tens of kb away (plastid rps12's 8-codon exon 3): either way a
/// junction the references have cis comes out trans. The cost lets a copy
/// next to its neighbour win when it scores within it, and leaves a slot
/// empty where a far match adds less; a junction trans in the references
/// (nad1, nad2, nad5, rps12 exon 1) costs nearly nothing. A junction to
/// another contig costs the same, so that moving an exon onto its
/// neighbour's contig is not a loss, but no exon is left out for one alone:
/// a fragmented assembly splits genes.
static TRANS_JUNCTION_PENALTY: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

pub fn set_trans_junction_penalty(cost: f64) {
    let cost = if cost.is_finite() { cost.max(0.0) } else { 0.0 };
    TRANS_JUNCTION_PENALTY.store(cost.to_bits(), std::sync::atomic::Ordering::Relaxed);
}

fn trans_junction_penalty() -> f64 {
    f64::from_bits(TRANS_JUNCTION_PENALTY.load(std::sync::atomic::Ordering::Relaxed))
}

/// What `junctions` cost under `--trans-junction-penalty` (see there).
fn junction_cost(template: &GeneTemplate, junctions: &[JunctionResult]) -> f64 {
    let cost = trans_junction_penalty();
    if cost == 0.0 {
        return 0.0;
    }
    junctions
        .iter()
        .filter(|j| j.kind != "cis")
        .filter_map(|j| template.cis_share(j.from_slot, j.to_slot))
        .map(|share| cost * share)
        .sum()
}
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
/// it - read through at a price (`--internal-stop-penalty`), or, where
/// U-to-C editing lets an alignment step over a stop, skipped - provided
/// at least half the profile is aligned before it, else one within this
/// many codons of its last aligned codon.
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
    /// Within `--end-margin` columns of the gene's start or end: reported,
    /// as orfedit reports its own, but the gene is read in one frame
    /// there (it ends at its own stop).
    pub near_end: bool,
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
    let circular = Circular::new();
    let (prepared, tasks) = prepare_gene(fasta_idx, fragments, template, params, merge_distance, exon_flank, &circular);
    let scans: Vec<DenovoScan> = tasks.iter().map(|t| denovo_scan(fasta_idx, template, params, t)).collect();
    finish_gene(fasta_idx, prepared, &tasks, &scans, template, params, junction_distance_threshold, exon_flank, &circular)
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
    circular: &Circular,
) -> (Prepared, Vec<DenovoTask>) {
    let k = template.slots.len();
    let candidates = build_candidates(fragments, merge_distance, circular);

    let exon_params = ScanParams { flank: exon_flank, end_margin: 0, ..*params }; // (the margin is the whole gene's)
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
            // (the contig itself: not round its origin again)
            let Some(len) = contig_len(fasta_idx, circular, contig) else { continue };
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
    let denovo_params = ScanParams { flank: 0, end_margin: 0, ..*params };
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
    circular: &Circular,
) -> Result<ReconstructionOutcome> {
    let k = template.slots.len();
    let Prepared { mut candidates, mut scores, mut assignment, cache } = prepared;
    let exon_params = ScanParams { flank: exon_flank, end_margin: 0, ..*params }; // (the margin is the whole gene's)
    // The end margin reads the gene chosen, at the end; the exons are chosen
    // without it (see there).
    let margin_params = params;
    let search_params = ScanParams { end_margin: 0, ..*params };
    let params = &search_params;
    let floors: Vec<f64> = template.slots.iter().map(|s| s.min_self_score).collect();
    let mut denovo_candidate_indices: std::collections::HashSet<usize> = std::collections::HashSet::new();

    // A slot's whole-contig hits at or above its floor, best first (ties
    // in task order, as a serial loop keeping the first best would).
    let hits_for = |si: usize| -> Vec<(f64, i64, i64, &DenovoTask)> {
        let mut hits: Vec<(f64, i64, i64, &DenovoTask)> = tasks
            .iter()
            .zip(scans)
            .filter(|(task, _)| task.slot == si)
            .filter_map(|(task, scan)| {
                scan.map(|(score, start, end)| {
                    let (start, end) = canonical(circular, &task.contig, start, end);
                    (score, start, end, task)
                })
            })
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
    // the whole-contig hit taken for each slot the fragments left empty
    let mut found_for_slot: Vec<Option<usize>> = vec![None; k];
    for &si in &empty_slots {
        if let Some(&(score, start, end, task)) = hits_for(si).first() {
            let new = Candidate { contig: task.contig.clone(), start, end, strand: task.strand, provenance: "denovo".to_string() };
            let idx = add_candidate(&mut candidates, &mut scores, &mut denovo_candidate_indices, new, si, score, k);
            found_for_slot[si] = Some(idx);
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
                // The whole-contig hit of a slot the fragments left empty
                // stays in. Where a species has lost introns, one fragment
                // holds several exons end to end; it and the merged run
                // score every one of those slots as well as the hit does,
                // and at two candidates a slot they would crowd it out -
                // yet each can fill one slot only (liverwort cox1 without
                // four of its nine introns: exons 6 to 9 are one genomic
                // exon, and 6 and 7 stayed empty).
                let found = found_for_slot[si].filter(|&d| ranked.iter().skip(2).any(|c| c.0 == d));
                ranked.truncate(2);
                if let Some(d) = found {
                    ranked.push((d, scores[d][si]));
                }
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
            &denovo_candidate_indices, junction_distance_threshold, false, &cache, circular,
        ) else {
            continue;
        };
        // On a tie, a gene within its contig over one across the origin:
        // where one copy of a plastid's inverted repeat runs across the
        // origin (an exon cut by it, read whole since 0.3.3), the two copies
        // make the same protein, and the one across the origin is the
        // harder to use (as 0.3.2 chose: rpl2 in 152 plastomes).
        let better = |prev: &Evaluated| {
            let eps = 1e-9 * prev.whole_gene_score.abs().max(1.0);
            result.whole_gene_score > prev.whole_gene_score + eps
                || (result.whole_gene_score >= prev.whole_gene_score - eps
                    && crosses_origin(&prev.exons, circular)
                    && !crosses_origin(&result.exons, circular))
        };
        if best.as_ref().is_none_or(|(_, prev)| better(prev)) {
            best = Some((bi, result));
        }
    }
    let (bi, unrefined) = best.context("no candidate assignment (not even the empty one) produced a result")?;
    let mut chosen = beam[bi].clone();
    let mut refined = evaluate_assignment(
        fasta_idx, &candidates, template, &chosen, &exon_params, params,
        &denovo_candidate_indices, junction_distance_threshold, true, &cache, circular,
    )
    .ok();

    // Short slots: whatever fills one (or nothing), try its best
    // whole-contig hits and the best match next to each neighbouring exon
    // in its place, and keep a replacement only if the refined whole-gene
    // alignment plus intron motifs scores higher with it. A short profile's
    // own score can't tell the exon from a chance match; the protein it
    // completes and the intron ends beside it can.
    let denovo_params = ScanParams { flank: 0, end_margin: 0, ..*params };
    for &si in &short_slots {
        let mut alternatives: Vec<usize> = Vec::new();
        for (score, start, end, task) in hits_for(si).into_iter().take(DENOVO_ALTERNATIVES) {
            let new = Candidate { contig: task.contig.clone(), start, end, strand: task.strand, provenance: "denovo".to_string() };
            alternatives.push(add_candidate(&mut candidates, &mut scores, &mut denovo_candidate_indices, new, si, score, k));
        }
        if let Some(current) = &refined {
            for (contig, strand, start, end) in nearby_windows(fasta_idx, circular, &current.exons, si + 1, 3 * template.slots[si].profile.len() as i64) {
                let Ok(o) = scan_gene(fasta_idx, &contig, start, end, strand, &template.slots[si].profile, &denovo_params) else {
                    continue;
                };
                if o.score >= floors[si] {
                    let (start, end) = canonical(circular, &contig, o.corrected_start as i64, o.corrected_end as i64);
                    let new = Candidate { contig, start, end, strand, provenance: "nearby".to_string() };
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
                &denovo_candidate_indices, junction_distance_threshold, true, &cache, circular,
            ) else {
                continue;
            };
            if refined.as_ref().is_none_or(|r| result.total() > r.total()) {
                chosen = alternative;
                refined = Some(result);
            }
        }
    }
    // Exons of two copies. An exon in a repeat - a plastid's inverted
    // repeat, a repeat of a mitogenome - has a twin that scores as it does,
    // the beam weighs two candidates a slot, and one fragment holding two
    // exons fills only one slot: so a gene could be joined from exon 1 of
    // one copy of the inverted repeat and exon 2 of the other, tens of kb
    // apart on opposite strands (plastid ndhB and rpl2 in hundreds of
    // assemblies). For an exon on a junction that is not cis, look next to
    // its neighbours (as for a short slot, and round the origin of a
    // circular contig) and take what is found there if the gene scores no
    // lower and more of its junctions are cis. "No lower" is net of what
    // its junctions cost (`--trans-junction-penalty`): a diverged copy of
    // an exon next to its neighbour scores a little less than the best
    // copy elsewhere, and gives way only to that. Where nothing next to
    // the neighbours will do, the exon is left out if the gene scores
    // higher, net of that cost, without it: a chance match tens of kb away
    // adds less than the junction costs. (Only this pass weighs that cost;
    // the beam and the short-slot search before it weigh whole-gene scores
    // alone, as before. Where an intron is lost, one fragment holds both
    // exons and fills one slot only, so the copy that is whole is not among
    // the beam's assignments: charged there, a chimera of two copies gave
    // way to a worse copy that this pass never saw. Charged in the
    // short-slot search, a short exon 3 of mito cox2 was paired with a far
    // exon 2, which a pass moving one exon at a time cannot undo: Vaccinium
    // cox2, intronless and whole at 949, came out trans at 574.)
    if let Some(mut current) = refined.take() {
        let apart = |e: &Evaluated| e.junctions.iter().filter(|j| j.kind != "cis").count();
        // Each round weighs every exon on a junction that is not cis and
        // takes the best change (the first that would do, slot by slot,
        // took exon 1 of one copy of the inverted repeat across the origin
        // where exon 2 of the other copy, beside exon 1, scored higher):
        // an exon found next to its neighbours, else, failing any, an exon
        // left out (only one far from a neighbour on its own contig).
        loop {
            if apart(&current) == 0 {
                break;
            }
            let eps = 1e-6 * current.whole_gene_score.abs().max(1.0);
            // the better of two changes: net of the junction cost, and on a
            // tie the one within its contig (as in the beam: two copies of
            // the inverted repeat beside each other's exons, one across the
            // origin)
            let preferable = |r: &Evaluated, b: &Evaluated| {
                r.net() > b.net() + eps
                    || (r.net() >= b.net() - eps && crosses_origin(&b.exons, circular) && !crosses_origin(&r.exons, circular))
            };
            let on_apart = |slot: usize| current.junctions.iter().any(|j| j.kind != "cis" && (j.from_slot == slot || j.to_slot == slot));
            // apart on one contig: far, or on the other strand
            let on_far = |slot: usize| {
                current.junctions.iter().any(|j| j.kind != "cis" && j.gap_bp >= 0 && (j.from_slot == slot || j.to_slot == slot))
            };
            // (afresh each round: what was no better beside a neighbour may
            // be once another exon has moved)
            let mut tried: Vec<(usize, String, char, usize, usize)> = Vec::new();
            let mut best_swap: Option<(Assignment, Evaluated)> = None;
            for si in (0..k).filter(|&si| on_apart(si + 1)) {
                let slot = si + 1;
                let profile = &template.slots[si].profile;
                for (contig, strand, start, end) in nearby_windows(fasta_idx, circular, &current.exons, slot, 3 * profile.len() as i64) {
                    let Ok(o) = scan_gene(fasta_idx, &contig, start, end, strand, profile, &denovo_params) else {
                        continue;
                    };
                    if o.score < floors[si] {
                        continue;
                    }
                    let (start, end) = canonical(circular, &contig, o.corrected_start as i64, o.corrected_end as i64);
                    let place = (si, contig.clone(), strand, start as usize, end as usize);
                    if tried.contains(&place) || current.exons.iter().any(|e| e.slot == slot && e.contig == contig && e.strand == strand && e.start == start as usize && e.end == end as usize) {
                        continue;
                    }
                    tried.push(place);
                    let new = Candidate { contig, start, end, strand, provenance: "nearby".to_string() };
                    let ci = add_candidate(&mut candidates, &mut scores, &mut denovo_candidate_indices, new, si, o.score, k);
                    if chosen.slot_to_candidate.contains(&Some(ci)) {
                        continue;
                    }
                    let mut alternative = chosen.clone();
                    alternative.slot_to_candidate[si] = Some(ci);
                    let Ok(result) = evaluate_assignment(
                        fasta_idx, &candidates, template, &alternative, &exon_params, params,
                        &denovo_candidate_indices, junction_distance_threshold, true, &cache, circular,
                    ) else {
                        continue;
                    };
                    if result.exons.len() >= current.exons.len()
                        && result.net() >= current.net() - eps
                        && apart(&result) < apart(&current)
                        && best_swap.as_ref().is_none_or(|(_, b)| preferable(&result, b))
                    {
                        best_swap = Some((alternative, result));
                    }
                }
            }
            let mut best = best_swap;
            // A run of exons moved together. Where no exon beside a
            // neighbour will do - one copy of a repeat lacks it - the other
            // copy may hold the whole gene, and moving one exon at a time
            // only moves the break: Dryopteris nad2 and nad7 took exons 1-3
            // from one copy of a repeat and exon 4 from the other, which
            // holds all four and, given alone, makes the same protein all
            // cis. For each junction apart, the exons on one side of it are
            // found again one by one beside the other side's, outward from
            // the junction, each beside the one found before; the run is
            // taken by the same rule as one exon.
            if best.is_none() {
                let bare = |e: &ExonResult| ExonResult {
                    slot: e.slot,
                    provenance: e.provenance.clone(),
                    contig: e.contig.clone(),
                    start: e.start,
                    end: e.end,
                    strand: e.strand,
                    score: e.score,
                    edits: Vec::new(),
                    via_denovo: e.via_denovo,
                    phase: e.phase,
                };
                let filled: Vec<usize> = current.exons.iter().map(|e| e.slot).collect();
                for j in current.junctions.iter().filter(|j| j.kind != "cis") {
                    for keep_after in [true, false] {
                        // the side kept, and the slots found again, outward
                        let (mut placed, order): (Vec<ExonResult>, Vec<usize>) = if keep_after {
                            (
                                current.exons.iter().filter(|e| e.slot >= j.to_slot).map(bare).collect(),
                                filled.iter().copied().filter(|&s| s <= j.from_slot).rev().collect(),
                            )
                        } else {
                            (
                                current.exons.iter().filter(|e| e.slot <= j.from_slot).map(bare).collect(),
                                filled.iter().copied().filter(|&s| s >= j.to_slot).collect(),
                            )
                        };
                        if order.len() < 2 {
                            continue; // one exon: the swaps above
                        }
                        let mut alternative = chosen.clone();
                        let mut whole = true;
                        for &slot in &order {
                            let si = slot - 1;
                            let profile = &template.slots[si].profile;
                            let mut found: Option<(f64, String, char, i64, i64)> = None;
                            for (contig, strand, start, end) in nearby_windows(fasta_idx, circular, &placed, slot, 3 * profile.len() as i64) {
                                let Ok(o) = scan_gene(fasta_idx, &contig, start, end, strand, profile, &denovo_params) else {
                                    continue;
                                };
                                if o.score < floors[si] || found.as_ref().is_some_and(|f| f.0 >= o.score) {
                                    continue;
                                }
                                let (s, e) = canonical(circular, &contig, o.corrected_start as i64, o.corrected_end as i64);
                                found = Some((o.score, contig, strand, s, e));
                            }
                            let Some((score, contig, strand, s, e)) = found else {
                                whole = false;
                                break;
                            };
                            placed.push(ExonResult {
                                slot,
                                provenance: "nearby".to_string(),
                                contig: contig.clone(),
                                start: s as usize,
                                end: e as usize,
                                strand,
                                score,
                                edits: Vec::new(),
                                via_denovo: true,
                                phase: 0,
                            });
                            let new = Candidate { contig, start: s, end: e, strand, provenance: "nearby".to_string() };
                            alternative.slot_to_candidate[si] =
                                Some(add_candidate(&mut candidates, &mut scores, &mut denovo_candidate_indices, new, si, score, k));
                        }
                        let used: Vec<usize> = alternative.slot_to_candidate.iter().flatten().copied().collect();
                        let injective = used.iter().enumerate().all(|(i, c)| !used[..i].contains(c));
                        if !whole || !injective || alternative.slot_to_candidate == chosen.slot_to_candidate {
                            continue;
                        }
                        let Ok(result) = evaluate_assignment(
                            fasta_idx, &candidates, template, &alternative, &exon_params, params,
                            &denovo_candidate_indices, junction_distance_threshold, true, &cache, circular,
                        ) else {
                            continue;
                        };
                        if result.exons.len() >= current.exons.len()
                            && result.net() >= current.net() - eps
                            && apart(&result) < apart(&current)
                            && best.as_ref().is_none_or(|(_, b)| preferable(&result, b))
                        {
                            best = Some((alternative, result));
                        }
                    }
                }
            }
            if best.is_none() && trans_junction_penalty() > 0.0 {
                for si in (0..k).filter(|&si| chosen.slot_to_candidate[si].is_some() && on_far(si + 1)) {
                    let mut alternative = chosen.clone();
                    alternative.slot_to_candidate[si] = None;
                    let Ok(result) = evaluate_assignment(
                        fasta_idx, &candidates, template, &alternative, &exon_params, params,
                        &denovo_candidate_indices, junction_distance_threshold, true, &cache, circular,
                    ) else {
                        continue;
                    };
                    if result.net() > current.net() + eps
                        && apart(&result) < apart(&current)
                        && best.as_ref().is_none_or(|(_, b)| result.net() > b.net())
                    {
                        best = Some((alternative, result));
                    }
                }
            }
            let Some((alternative, result)) = best else { break };
            chosen = alternative;
            current = result;
        }
        refined = Some(current);
    }

    // The gene chosen, read in one frame where its only frame changes lie
    // within `end_margin` columns of either end (those still reported, as
    // near_end). Applied to every alignment the search weighed, the margin
    // changed which exons were chosen: of 894 1-bp indels planted in spliced
    // plastid genes 612 were found, against 707 without it, and calls
    // elsewhere in those genes went from 14 to 69.
    if margin_params.end_margin > 0 {
        if let Ok(r) = evaluate_assignment(
            fasta_idx, &candidates, template, &chosen, &exon_params, margin_params,
            &denovo_candidate_indices, junction_distance_threshold, true, &cache, circular,
        ) {
            refined = Some(r);
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
fn nearby_windows(fasta_idx: &FastaIndex, circular: &Circular, exons: &[ExonResult], slot: usize, min_len: i64) -> Vec<(String, char, i64, i64)> {
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
        windows.push((p.contig.clone(), p.strand, s, e));
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
        windows.push((n.contig.clone(), n.strand, s, e));
    }
    // On a circular contig a window past either end runs on round the
    // origin, as one window into the contig's continuation (`wrap_circular`;
    // what is found there goes back by `canonical`). Cut in two at the
    // origin, as before 0.3.3, it could not hold an exon across the origin
    // itself: fern ndhB, whose exon 2 crosses it beside exon 1, came from
    // the other copy of the inverted repeat.
    let mut out: Vec<(String, char, i64, i64)> = Vec::new();
    for (contig, strand, s, e) in windows {
        let seq_len = fasta_idx.sequences.get(&contig).map_or(i64::MAX, |q| q.len() as i64);
        let (s, e) = match circular_len(circular, &contig) {
            Some(len) if s < 0 => (s + len, e + len),
            _ => (s, e),
        };
        out.push((contig, strand, s.max(0), e.min(seq_len)));
    }
    out.retain(|w| w.3 - w.2 >= min_len);
    out.dedup();
    out
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
    /// What its junctions cost (`--trans-junction-penalty`).
    junction_cost: f64,
}

impl Evaluated {
    /// What alternatives for a short slot are compared by.
    fn total(&self) -> f64 {
        self.whole_gene_score + self.motif_bonus
    }

    /// What exons of two copies are compared by.
    fn net(&self) -> f64 {
        self.whole_gene_score - self.junction_cost
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
    circular: &Circular,
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
            // an exon wholly past a circular contig's end, in its
            // continuation, is the one at its start
            let (s2, e2) = canonical(circular, &p.contig, s.max(0), e.max(0));
            let shift = (s.max(0) - s2) as usize;
            let mut edits = p.edits;
            for x in edits.iter_mut() {
                x.genomic_pos -= shift;
            }
            ExonResult {
                slot: p.slot,
                provenance: p.provenance,
                contig: p.contig,
                start: s2 as usize,
                end: e2 as usize,
                strand: p.strand,
                score: p.score,
                edits,
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
    let align_exons = |exons: &[ExonResult]| -> Result<(Vec<u8>, Vec<(usize, usize)>, orfedit::align::NtAlignResult, Vec<orfedit::align::Frameshift>)> {
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
        // Frame changes within `end_margin` profile columns of the gene's
        // start or end buy only its first or last few columns, and cannot
        // be told from a gene that starts or ends differently there (see
        // orfedit): if those are all it has, the gene is read in one frame.
        // They are still reported (`near_end`), as orfedit's are: dropped
        // altogether, 7 and 13 points of planted-indel recall in spliced mito
        // and plastid genes went with them.
        let n = template.whole_gene_profile.len();
        let margin = params.end_margin.min(n / 4);
        let mut near_end = Vec::new();
        if !whole.frameshifts.is_empty() && whole.frameshifts.iter().all(|f| f.after_column < margin || f.after_column + margin > n) {
            near_end = std::mem::take(&mut whole.frameshifts);
            whole = align_whole(f64::INFINITY);
        }
        Ok((concat, offsets, whole, near_end))
    };
    let (concat, mut offsets, mut whole, mut near_end) = align_exons(&exons)?;
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
            (_, offsets, whole, near_end) = align_exons(&exons)?;
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
            circular_len(circular, &a.contig),
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
    for (f, near_end) in whole.frameshifts.iter().map(|f| (f, false)).chain(near_end.iter().map(|f| (f, true))) {
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
            near_end,
        });
    }

    let junction_cost = junction_cost(template, &junctions);
    Ok(Evaluated { exons, junctions, whole_gene_score: whole.score, whole_gene_n_edits, protein, frameshifts, motif_bonus, junction_cost })
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
        // ...or one the alignment read through as an internal stop
        // (orfedit's EditCosts::internal_stop), whichever comes first
        let read_through = aligned.iter().find(|a| a.2 == orfedit::codon::STOP && a.1 >= last_off).map(|a| (a.1, a.0.saturating_sub(1)));
        if let Some(r) = read_through {
            if interior.is_none_or(|i| r.0 < i.0) {
                interior = Some(r);
            }
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
        Profile { gene: gene.to_string(), columns, edited: None }
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
            junction_cis: Vec::new(),
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
            junction_cis: Vec::new(),
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
            junction_cis: Vec::new(),
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
            junction_cis: Vec::new(),
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
            junction_cis: Vec::new(),
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
            junction_cis: Vec::new(),
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
        Profile { gene: gene.to_string(), columns: cols, edited: None }
    }

    /// A profile scoring each residue of `protein` 10 in its column.
    fn protein_profile(gene: &str, protein: &[u8]) -> Profile {
        let cols: Vec<[(u8, f64); 1]> = protein.iter().map(|&aa| [(aa, 10.0)]).collect();
        let refs: Vec<&[(u8, f64)]> = cols.iter().map(|c| &c[..]).collect();
        custom_profile(gene, &refs)
    }

    #[test]
    fn exons_of_two_copies_give_way_to_the_copy_next_to_the_neighbour() {
        // Exon 2 has two copies, as in a plastid's inverted repeat. The
        // fragments name the one on the other strand, mid-contig; the other
        // lies just after exon 1, round the origin of the circular contig.
        let e1 = [b"ATG".to_vec(), b"AAATGGCAT".repeat(10)].concat(); // M (K W H)x10
        let e2 = b"TATTTTTGGTGT".repeat(8); // (Y F W C)x8
        let mut seq = vec![b'C'; 600];
        seq[9..9 + e2.len()].copy_from_slice(&e2);
        let mut rc = e2.clone();
        orfedit::fasta::revcomp_in_place(&mut rc);
        seq[250..250 + rc.len()].copy_from_slice(&rc);
        seq[480..480 + e1.len()].copy_from_slice(&e1);
        let mut sequences = HashMap::new();
        sequences.insert("ctg1".to_string(), seq);
        let linear_idx = FastaIndex { sequences };
        let mut circular_idx = FastaIndex { sequences: linear_idx.sequences.clone() };
        let circular = wrap_circular(&mut circular_idx, &["ctg1".to_string()].into_iter().collect());

        let p1 = [b"M".to_vec(), b"KWH".repeat(10)].concat();
        let p2 = b"YFWC".repeat(8);
        let template = GeneTemplate {
            gene: "test".to_string(),
            slots: vec![
                SlotTemplate { profile: protein_profile("e1", &p1), min_self_score: 100.0, lead: 0, trail: 0 },
                SlotTemplate { profile: protein_profile("e2", &p2), min_self_score: 100.0, lead: 0, trail: 0 },
            ],
            whole_gene_profile: protein_profile("whole", &[p1.clone(), p2.clone()].concat()),
            junction_cis: Vec::new(),
        };
        let fragments = vec![
            RawFragment { contig: "ctg1".to_string(), start: 480, end: 480 + e1.len() as i64, strand: '+' },
            RawFragment { contig: "ctg1".to_string(), start: 250, end: 250 + e2.len() as i64, strand: '-' },
        ];
        let params = ScanParams::default();
        let run = |fasta_idx: &FastaIndex, circular: &Circular| {
            let (prepared, tasks) = prepare_gene(fasta_idx, &fragments, &template, &params, 0, 0, circular);
            let scans: Vec<DenovoScan> = tasks.iter().map(|t| denovo_scan(fasta_idx, &template, &params, t)).collect();
            finish_gene(fasta_idx, prepared, &tasks, &scans, &template, &params, 100, 0, circular).unwrap()
        };

        // a linear contig: nothing after exon 1, so the far copy stays
        let linear = run(&linear_idx, &Circular::new());
        assert_eq!(linear.n_filled, 2);
        assert_eq!((linear.exons[1].start, linear.exons[1].strand), (250, '-'));
        assert_eq!(linear.junctions[0].kind, "trans");

        // circular: exon 2 found round the origin, the same protein, cis
        let circ = run(&circular_idx, &circular);
        assert_eq!(circ.protein, linear.protein);
        let (a, b) = (&circ.exons[1], &linear.exons[1]);
        assert_eq!(a.strand, '+');
        assert!(a.start >= 9 && a.end <= 9 + e2.len(), "exon 2 at {}..{}", a.start, a.end);
        assert_eq!(a.end - a.start, b.end - b.start);
        assert_eq!(circ.junctions[0].kind, "cis");
        assert_eq!(circ.junctions[0].gap_bp, 600 - circ.exons[0].end as i64 + a.start as i64);
    }

    #[test]
    fn an_exon_across_the_origin_is_found_beside_its_neighbour() {
        // As fern ndhB: exon 2's copy beside exon 1 runs across the origin
        // of the circular contig (bases 590..600 then 0..86); the fragments
        // name only the other copy, on the other strand mid-contig.
        let e1 = [b"ATG".to_vec(), b"AAATGGCAT".repeat(10)].concat(); // M (K W H)x10
        let e2 = b"TATTTTTGGTGT".repeat(8); // (Y F W C)x8, 96 bp
        let mut seq = vec![b'C'; 600];
        seq[590..600].copy_from_slice(&e2[..10]);
        seq[..86].copy_from_slice(&e2[10..]);
        let mut rc = e2.clone();
        orfedit::fasta::revcomp_in_place(&mut rc);
        seq[250..250 + rc.len()].copy_from_slice(&rc);
        seq[480..480 + e1.len()].copy_from_slice(&e1);
        let mut sequences = HashMap::new();
        sequences.insert("ctg1".to_string(), seq);
        let mut fasta_idx = FastaIndex { sequences };
        let circular = wrap_circular(&mut fasta_idx, &["ctg1".to_string()].into_iter().collect());
        assert_eq!(circular["ctg1"], 600);

        let p1 = [b"M".to_vec(), b"KWH".repeat(10)].concat();
        let p2 = b"YFWC".repeat(8);
        let template = GeneTemplate {
            gene: "test".to_string(),
            slots: vec![
                SlotTemplate { profile: protein_profile("e1", &p1), min_self_score: 100.0, lead: 0, trail: 0 },
                SlotTemplate { profile: protein_profile("e2", &p2), min_self_score: 100.0, lead: 0, trail: 0 },
            ],
            whole_gene_profile: protein_profile("whole", &[p1.clone(), p2.clone()].concat()),
            junction_cis: Vec::new(),
        };
        let fragments = vec![
            RawFragment { contig: "ctg1".to_string(), start: 480, end: 480 + e1.len() as i64, strand: '+' },
            RawFragment { contig: "ctg1".to_string(), start: 250, end: 250 + e2.len() as i64, strand: '-' },
        ];
        let params = ScanParams::default();
        let run = |fasta_idx: &FastaIndex, circular: &Circular| {
            let (prepared, tasks) = prepare_gene(fasta_idx, &fragments, &template, &params, 0, 0, circular);
            let scans: Vec<DenovoScan> = tasks.iter().map(|t| denovo_scan(fasta_idx, &template, &params, t)).collect();
            finish_gene(fasta_idx, prepared, &tasks, &scans, &template, &params, 100, 0, circular).unwrap()
        };
        let out = run(&fasta_idx, &circular);

        assert_eq!(out.n_filled, 2);
        let e = &out.exons[1];
        assert_eq!(e.strand, '+');
        // across the origin: it ends past the contig's length
        assert!(e.start >= 590 && e.start < 600 && e.end > 600 && e.end <= 686, "exon 2 at {}..{}", e.start, e.end);
        assert_eq!(out.junctions[0].kind, "cis");

        // the same as with exon 2 whole beside exon 1, on a linear contig
        let mut seq = fasta_idx.sequences["ctg1"][..600].to_vec();
        seq[..86].fill(b'C');
        seq.extend(vec![b'C'; 100]);
        seq[590..686].copy_from_slice(&e2);
        let mut sequences = HashMap::new();
        sequences.insert("ctg1".to_string(), seq);
        let control = run(&FastaIndex { sequences }, &Circular::new());
        assert_eq!(control.junctions[0].kind, "cis");
        assert_eq!(out.protein, control.protein);
        assert_eq!((e.start, e.end), (control.exons[1].start, control.exons[1].end));
    }

    #[test]
    fn a_run_of_exons_moves_to_the_copy_that_holds_the_whole_gene() {
        // Two copies of a repeat, further apart than an exon is looked for
        // beside its neighbour: A (at 100) holds exons 1 and 2, B (at 8000)
        // all three. The fragments name exons 1-2 of A and exon 3 of B, so
        // one exon at a time only moves the break.
        let e1 = [b"ATG".to_vec(), b"AAATGGCAT".repeat(10)].concat(); // M (K W H)x10
        let e2 = b"TATTTTTGGTGT".repeat(8); // (Y F W C)x8
        let e3 = b"GATGAACCTGGT".repeat(8); // (D E P G)x8
        let intron = b"GTAAGTTTTTTTTTTTTTTTTTTTTTTTAG".to_vec();
        let copy_a = [e1.clone(), intron.clone(), e2.clone()].concat();
        let copy_b = [e1.clone(), intron.clone(), e2.clone(), intron.clone(), e3.clone()].concat();
        let mut seq = vec![b'C'; 9000];
        seq[100..100 + copy_a.len()].copy_from_slice(&copy_a);
        seq[8000..8000 + copy_b.len()].copy_from_slice(&copy_b);
        let mut sequences = HashMap::new();
        sequences.insert("ctg1".to_string(), seq);
        let fasta_idx = FastaIndex { sequences };
        let (p1, p2, p3) = ([b"M".to_vec(), b"KWH".repeat(10)].concat(), b"YFWC".repeat(8), b"DEPG".repeat(8));
        let template = GeneTemplate {
            gene: "test".to_string(),
            slots: vec![
                SlotTemplate { profile: protein_profile("e1", &p1), min_self_score: 100.0, lead: 0, trail: 0 },
                SlotTemplate { profile: protein_profile("e2", &p2), min_self_score: 100.0, lead: 0, trail: 0 },
                SlotTemplate { profile: protein_profile("e3", &p3), min_self_score: 100.0, lead: 0, trail: 0 },
            ],
            whole_gene_profile: protein_profile("whole", &[p1.clone(), p2.clone(), p3.clone()].concat()),
            junction_cis: Vec::new(),
        };
        let b3 = 8000 + (e1.len() + e2.len() + 2 * intron.len()) as i64;
        let fragments = vec![
            RawFragment { contig: "ctg1".to_string(), start: 100, end: 100 + e1.len() as i64, strand: '+' },
            RawFragment { contig: "ctg1".to_string(), start: 100 + (e1.len() + intron.len()) as i64, end: 100 + copy_a.len() as i64, strand: '+' },
            RawFragment { contig: "ctg1".to_string(), start: b3, end: b3 + e3.len() as i64, strand: '+' },
        ];
        let params = ScanParams::default();
        let circular = Circular::new();
        let (prepared, tasks) = prepare_gene(&fasta_idx, &fragments, &template, &params, 0, 0, &circular);
        let scans: Vec<DenovoScan> = tasks.iter().map(|t| denovo_scan(&fasta_idx, &template, &params, t)).collect();
        let out = finish_gene(&fasta_idx, prepared, &tasks, &scans, &template, &params, 100, 0, &circular).unwrap();
        assert_eq!(out.n_filled, 3);
        assert!(out.junctions.iter().all(|j| j.kind == "cis"), "{:?}", out.junctions.iter().map(|j| j.kind).collect::<Vec<_>>());
        assert!(out.exons.iter().all(|e| e.start >= 8000), "{:?}", out.exons.iter().map(|e| e.start).collect::<Vec<_>>());
    }

    #[test]
    fn fragments_of_an_exon_across_the_origin_are_merged() {
        // the same exon 2 across the origin, named by its two pieces
        let e2 = b"TATTTTTGGTGT".repeat(8);
        let mut seq = vec![b'C'; 600];
        seq[590..600].copy_from_slice(&e2[..10]);
        seq[..86].copy_from_slice(&e2[10..]);
        let mut sequences = HashMap::new();
        sequences.insert("ctg1".to_string(), seq);
        let mut fasta_idx = FastaIndex { sequences };
        let circular = wrap_circular(&mut fasta_idx, &["ctg1".to_string()].into_iter().collect());
        let p2 = b"YFWC".repeat(8);
        let template = GeneTemplate {
            gene: "test".to_string(),
            slots: vec![SlotTemplate { profile: protein_profile("e2", &p2), min_self_score: 100.0, lead: 0, trail: 0 }],
            whole_gene_profile: protein_profile("whole", &p2),
            junction_cis: Vec::new(),
        };
        let fragments = vec![
            RawFragment { contig: "ctg1".to_string(), start: 0, end: 86, strand: '+' },
            RawFragment { contig: "ctg1".to_string(), start: 590, end: 600, strand: '+' },
        ];
        let params = ScanParams::default();
        let (prepared, tasks) = prepare_gene(&fasta_idx, &fragments, &template, &params, 100, 0, &circular);
        let scans: Vec<DenovoScan> = tasks.iter().map(|t| denovo_scan(&fasta_idx, &template, &params, t)).collect();
        let out = finish_gene(&fasta_idx, prepared, &tasks, &scans, &template, &params, 100, 0, &circular).unwrap();
        assert_eq!(out.n_filled, 1);
        assert_eq!((out.exons[0].start, out.exons[0].end), (590, 686));
        assert_eq!(out.protein, String::from_utf8(p2).unwrap());
    }

    /// Runs the three stages on one gene (as a batch would), junctions
    /// cis within 100 bp.
    fn run_gene(fasta_idx: &FastaIndex, fragments: &[RawFragment], template: &GeneTemplate) -> ReconstructionOutcome {
        let params = ScanParams::default();
        let (prepared, tasks) = prepare_gene(fasta_idx, fragments, template, &params, 0, 0, &Circular::new());
        let scans: Vec<DenovoScan> = tasks.iter().map(|t| denovo_scan(fasta_idx, template, &params, t)).collect();
        finish_gene(fasta_idx, prepared, &tasks, &scans, template, &params, 100, 0, &Circular::new()).unwrap()
    }

    // Every test that sets the run-wide cost sets it to this, as tests run
    // side by side (a template without junction shares costs nothing).
    const TEST_TRANS_JUNCTION_PENALTY: f64 = 50.0;

    #[test]
    fn a_diverged_copy_next_to_the_neighbour_beats_a_better_one_far_away() {
        // Exon 2 has two copies: a diverged one just after exon 1 (two of
        // its 32 residues differ: 30 less) and an exact one on the other
        // strand 1 kb away, which the fragments name - plastid clpP in
        // Aquilegia, without introns, in two copies 16 kb apart.
        set_trans_junction_penalty(TEST_TRANS_JUNCTION_PENALTY);
        let e1 = [b"ATG".to_vec(), b"AAATGGCAT".repeat(10)].concat(); // M (K W H)x10
        let e2 = b"TATTTTTGGTGT".repeat(8); // (Y F W C)x8
        let mut near = e2.clone();
        near[9..12].copy_from_slice(b"AAA"); // C -> K
        near[45..48].copy_from_slice(b"AAA");
        let mut seq = vec![b'C'; 1500];
        seq[100..100 + e1.len()].copy_from_slice(&e1);
        let near_at = 100 + e1.len() + 60;
        seq[near_at..near_at + near.len()].copy_from_slice(&near);
        let mut rc = e2.clone();
        orfedit::fasta::revcomp_in_place(&mut rc);
        seq[1200..1200 + rc.len()].copy_from_slice(&rc);
        let mut sequences = HashMap::new();
        sequences.insert("ctg1".to_string(), seq);
        let fasta_idx = FastaIndex { sequences };

        let p1 = [b"M".to_vec(), b"KWH".repeat(10)].concat();
        let p2 = b"YFWC".repeat(8);
        let template = |junction_cis: Vec<Option<f64>>| GeneTemplate {
            gene: "test".to_string(),
            slots: vec![
                SlotTemplate { profile: protein_profile("e1", &p1), min_self_score: 100.0, lead: 0, trail: 0 },
                SlotTemplate { profile: protein_profile("e2", &p2), min_self_score: 100.0, lead: 0, trail: 0 },
            ],
            whole_gene_profile: protein_profile("whole", &[p1.clone(), p2.clone()].concat()),
            junction_cis,
        };
        let fragments = vec![
            RawFragment { contig: "ctg1".to_string(), start: 100, end: 100 + e1.len() as i64, strand: '+' },
            RawFragment { contig: "ctg1".to_string(), start: 1200, end: 1200 + e2.len() as i64, strand: '-' },
        ];

        // nothing known of the junction: the better copy, far away
        let unknown = run_gene(&fasta_idx, &fragments, &template(Vec::new()));
        assert_eq!((unknown.exons[1].start, unknown.exons[1].strand), (1200, '-'));
        assert_eq!(unknown.junctions[0].kind, "trans");

        // cis in the references: the copy next to exon 1
        let cis = run_gene(&fasta_idx, &fragments, &template(vec![Some(1.0)]));
        assert_eq!(cis.n_filled, 2);
        let e = &cis.exons[1];
        assert_eq!(e.strand, '+');
        assert!(e.start >= near_at && e.end <= near_at + near.len(), "exon 2 at {}..{}", e.start, e.end);
        assert_eq!(cis.junctions[0].kind, "cis");
        assert!(cis.whole_gene_score < unknown.whole_gene_score);

        // trans in the references: the better copy stays
        let trans = run_gene(&fasta_idx, &fragments, &template(vec![Some(0.0)]));
        assert_eq!((trans.exons[1].start, trans.exons[1].strand), (1200, '-'));
    }

    #[test]
    fn an_exon_on_another_contig_is_not_left_out() {
        // As a_short_exon_far_away_..., but the only match for exon 2 is on
        // another contig: a fragmented assembly splits genes, so it stays.
        set_trans_junction_penalty(TEST_TRANS_JUNCTION_PENALTY);
        let e1 = [b"ATG".to_vec(), b"AAATGGCAT".repeat(10)].concat();
        let mut seq1 = vec![b'C'; 400];
        seq1[100..100 + e1.len()].copy_from_slice(&e1);
        let mut seq2 = vec![b'C'; 400];
        seq2[200..206].copy_from_slice(b"TGGTGG");
        let mut sequences = HashMap::new();
        sequences.insert("ctg1".to_string(), seq1);
        sequences.insert("ctg2".to_string(), seq2);
        let fasta_idx = FastaIndex { sequences };
        let p1 = [b"M".to_vec(), b"KWH".repeat(10)].concat();
        let template = GeneTemplate {
            gene: "test".to_string(),
            slots: vec![
                SlotTemplate { profile: protein_profile("e1", &p1), min_self_score: 100.0, lead: 0, trail: 0 },
                SlotTemplate { profile: protein_profile("e2", b"WW"), min_self_score: 10.0, lead: 0, trail: 0 },
            ],
            whole_gene_profile: protein_profile("whole", &[p1.clone(), b"WW".to_vec()].concat()),
            junction_cis: vec![Some(1.0)],
        };
        let fragments = vec![
            RawFragment { contig: "ctg1".to_string(), start: 100, end: 100 + e1.len() as i64, strand: '+' },
            RawFragment { contig: "ctg2".to_string(), start: 150, end: 250, strand: '+' },
        ];
        let outcome = run_gene(&fasta_idx, &fragments, &template);
        assert_eq!(outcome.n_filled, 2);
        assert_eq!(outcome.exons[1].contig, "ctg2");
        assert_eq!((outcome.junctions[0].kind, outcome.junctions[0].gap_bp), ("trans", -1));
    }

    #[test]
    fn a_copy_without_the_intron_beats_a_diverged_copy_with_it() {
        // Copy A has lost the intron: one fragment holds both exons, and
        // fills one slot only, so A's two exons together are not among the
        // beam's assignments. Copy B, 700 bp on, keeps the intron and
        // differs at two residues of each exon (plastid rpl2 in Cerastium).
        // The beam takes exon 1 of A and exon 2 of B, and the pass for two
        // copies then finds A's exon 2 beside its exon 1. Charging the
        // junction in the beam instead chose B whole.
        set_trans_junction_penalty(TEST_TRANS_JUNCTION_PENALTY);
        let e1 = [b"ATG".to_vec(), b"AAATGGCAT".repeat(10)].concat(); // M (K W H)x10
        let e2 = b"TATTTTTGGTGT".repeat(8); // (Y F W C)x8
        let (mut b1, mut b2) = (e1.clone(), e2.clone());
        b1[3..6].copy_from_slice(b"GAA"); // K -> E
        b1[30..33].copy_from_slice(b"GAA");
        b2[9..12].copy_from_slice(b"AAA"); // C -> K
        b2[45..48].copy_from_slice(b"AAA");
        let a = [e1.clone(), e2.clone()].concat();
        let b = [b1.clone(), vec![b'C'; 60], b2.clone()].concat();
        let mut seq = vec![b'C'; 1600];
        seq[100..100 + a.len()].copy_from_slice(&a);
        seq[800..800 + b.len()].copy_from_slice(&b);
        let mut sequences = HashMap::new();
        sequences.insert("ctg1".to_string(), seq);
        let fasta_idx = FastaIndex { sequences };

        let p1 = [b"M".to_vec(), b"KWH".repeat(10)].concat();
        let p2 = b"YFWC".repeat(8);
        let template = GeneTemplate {
            gene: "test".to_string(),
            slots: vec![
                SlotTemplate { profile: protein_profile("e1", &p1), min_self_score: 100.0, lead: 0, trail: 0 },
                SlotTemplate { profile: protein_profile("e2", &p2), min_self_score: 100.0, lead: 0, trail: 0 },
            ],
            whole_gene_profile: protein_profile("whole", &[p1.clone(), p2.clone()].concat()),
            junction_cis: vec![Some(1.0)],
        };
        let b2_at = 800 + b1.len() as i64 + 60;
        let fragments = vec![
            RawFragment { contig: "ctg1".to_string(), start: 100, end: 100 + a.len() as i64, strand: '+' },
            RawFragment { contig: "ctg1".to_string(), start: 800, end: 800 + b1.len() as i64, strand: '+' },
            RawFragment { contig: "ctg1".to_string(), start: b2_at, end: b2_at + b2.len() as i64, strand: '+' },
        ];
        let outcome = run_gene(&fasta_idx, &fragments, &template);
        assert_eq!(outcome.n_filled, 2);
        for e in &outcome.exons {
            assert!(e.start >= 100 && e.end <= 100 + a.len(), "exon {} at {}..{}", e.slot, e.start, e.end);
        }
        assert_eq!(outcome.junctions[0].kind, "cis");
        assert_eq!(outcome.protein, String::from_utf8([p1, p2].concat()).unwrap());
    }

    #[test]
    fn a_short_exon_far_away_is_left_out_where_the_references_join_it_cis() {
        // Exon 2 is two codons; nothing matches it near exon 1, and a
        // chance TGGTGG 1 kb on fills it - plastid rps12's exon 3 in some
        // liverworts and horsetails. Its 20 points do not pay for a
        // junction the references have cis.
        set_trans_junction_penalty(TEST_TRANS_JUNCTION_PENALTY);
        let e1 = [b"ATG".to_vec(), b"AAATGGCAT".repeat(10)].concat();
        let mut seq = vec![b'C'; 1500];
        seq[100..100 + e1.len()].copy_from_slice(&e1);
        seq[1200..1206].copy_from_slice(b"TGGTGG");
        let mut sequences = HashMap::new();
        sequences.insert("ctg1".to_string(), seq);
        let fasta_idx = FastaIndex { sequences };

        let p1 = [b"M".to_vec(), b"KWH".repeat(10)].concat();
        let template = |junction_cis: Vec<Option<f64>>| GeneTemplate {
            gene: "test".to_string(),
            slots: vec![
                SlotTemplate { profile: protein_profile("e1", &p1), min_self_score: 100.0, lead: 0, trail: 0 },
                SlotTemplate { profile: protein_profile("e2", b"WW"), min_self_score: 10.0, lead: 0, trail: 0 },
            ],
            whole_gene_profile: protein_profile("whole", &[p1.clone(), b"WW".to_vec()].concat()),
            junction_cis,
        };
        let fragments = vec![RawFragment { contig: "ctg1".to_string(), start: 100, end: 100 + e1.len() as i64, strand: '+' }];

        let unknown = run_gene(&fasta_idx, &fragments, &template(Vec::new()));
        assert_eq!(unknown.n_filled, 2);
        assert_eq!(unknown.exons[1].start, 1200);
        assert_eq!(unknown.junctions[0].kind, "trans");

        let cis = run_gene(&fasta_idx, &fragments, &template(vec![Some(1.0)]));
        assert_eq!(cis.n_filled, 1);
        assert!(!cis.complete);
        assert!(cis.junctions.is_empty());
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
            junction_cis: Vec::new(),
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
            junction_cis: Vec::new(),
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

    #[test]
    fn a_frameshift_near_the_gene_end_is_reported_and_not_read_through() {
        // as above, the extra A after exon 2's 27th codon: 3 columns from the end
        let half = "ATGAAAGTTTGGTTT".repeat(3);
        let mut seq = b"CCC".to_vec();
        seq.extend(half.as_bytes()); // exon1: 3..48
        seq.extend(b"GTGCGTTTTTTTTTTTTTTTTTTTTAT"); // intron: 48..75
        seq.extend("ATGAAAGTTTGGTTT".repeat(5).as_bytes()); // exon2 codons 1-25: 75..150
        seq.extend(b"ATGAA"); // 150..155, then the extra base
        seq.push(b'A'); // 155
        seq.extend(b"AGTTTGGTTTTAACCC");
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
            junction_cis: Vec::new(),
        };
        let fragments = vec![
            RawFragment { contig: "ctg1".to_string(), start: 3, end: 48, strand: '+' },
            RawFragment { contig: "ctg1".to_string(), start: 75, end: 172, strand: '+' },
        ];
        // no margin: the frame change is followed and not marked
        let params = ScanParams { frameshift_penalty: 10.0, ..ScanParams::default() };
        let followed = reconstruct_gene(&fasta_idx, &fragments, &template, &params, 0, 1000, 0).unwrap();
        assert_eq!(followed.frameshifts.len(), 1);
        assert!(!followed.frameshifts[0].near_end);
        // a margin: reported as near the end, and the gene read in one frame
        let params = ScanParams { frameshift_penalty: 10.0, end_margin: 10, ..ScanParams::default() };
        let near = reconstruct_gene(&fasta_idx, &fragments, &template, &params, 0, 1000, 0).unwrap();
        assert_eq!(near.frameshifts.len(), 1);
        let f = &near.frameshifts[0];
        assert!(f.near_end);
        assert_eq!((f.slot, f.start, f.end, f.n_bases), (followed.frameshifts[0].slot, followed.frameshifts[0].start, followed.frameshifts[0].end, 1));
        assert_ne!(near.protein, followed.protein);
    }
}
