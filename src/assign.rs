//! Assigning candidate exon fragments to exon-number "slots" (1..K), and
//! classifying each junction between assigned slots as cis or trans for the
//! species at hand - never assumed from a reference topology (see README:
//! only exon *identity/order* is fixed across lineages, not which specific
//! junction is cis vs trans).

#[derive(Clone, Debug, Default)]
pub struct Assignment {
    /// slot_to_candidate[slot] = Some(candidate_index) or None if unfilled.
    pub slot_to_candidate: Vec<Option<usize>>,
    pub total_score: f64,
}

/// Best injective partial matching of candidates to slots maximizing total
/// score, subject to a per-slot floor (candidates scoring below a slot's
/// floor are never eligible for it). Brute-force recursive search - fine
/// given K (exon count) is at most 5 here and the candidate pool is capped
/// by the caller before this is called.
pub fn best_assignment(scores: &[Vec<f64>], floors: &[f64]) -> Assignment {
    let k = floors.len();
    let m = scores.len();
    let mut used = vec![false; m];
    let mut current: Vec<Option<usize>> = vec![None; k];
    let mut best = Assignment { slot_to_candidate: vec![None; k], total_score: 0.0 };

    fn recurse(
        slot: usize,
        k: usize,
        scores: &[Vec<f64>],
        floors: &[f64],
        used: &mut [bool],
        current: &mut Vec<Option<usize>>,
        current_score: f64,
        best: &mut Assignment,
    ) {
        if slot == k {
            if current_score >= best.total_score {
                best.total_score = current_score;
                best.slot_to_candidate = current.clone();
            }
            return;
        }
        // Leave this slot unfilled.
        recurse(slot + 1, k, scores, floors, used, current, current_score, best);
        // Or try every eligible, unused candidate.
        for c in 0..scores.len() {
            if used[c] {
                continue;
            }
            let s = scores[c][slot];
            if s < floors[slot] {
                continue;
            }
            used[c] = true;
            current[slot] = Some(c);
            recurse(slot + 1, k, scores, floors, used, current, current_score + s, best);
            current[slot] = None;
            used[c] = false;
        }
    }

    if m > 0 {
        recurse(0, k, scores, floors, &mut used, &mut current, 0.0, &mut best);
    }
    best
}

/// Every plausible injective assignment constructible from each slot's
/// top-`beam_width` candidates (plus always the option to leave a slot
/// unfilled) - NOT just the single greedy-optimal one. Per-slot-independent
/// scoring can't distinguish "these 3 real, decently-scoring fragments go
/// to slots 3,4,5" from "...go to slots 5,4,3" when the profiles for those
/// slots aren't individually discriminating enough - confirmed directly
/// against real ground truth (a genuine `nad2` misassignment in
/// `Arabidopsis_thaliana` - see `transsplice`'s README). The caller is
/// expected to score each returned assignment some other way (this crate's
/// whole-gene alignment pass) and pick among them - that's the only
/// signal that can actually tell a permutation error from a correct
/// assignment when local per-slot scores are ambiguous between them.
///
/// Bounded and cheap: `(beam_width + 1)^k` combinations at most (the `+1`
/// is "leave unfilled"), trivial for the gene sizes this project handles
/// (k <= 5).
pub fn best_assignments_beam(scores: &[Vec<f64>], floors: &[f64], beam_width: usize) -> Vec<Assignment> {
    let k = floors.len();
    let m = scores.len();

    let mut shortlists: Vec<Vec<(usize, f64)>> = Vec::with_capacity(k);
    for slot in 0..k {
        let mut cands: Vec<(usize, f64)> = (0..m)
            .filter_map(|c| {
                let s = scores[c][slot];
                if s >= floors[slot] { Some((c, s)) } else { None }
            })
            .collect();
        cands.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        cands.truncate(beam_width);
        shortlists.push(cands);
    }

    let mut results = Vec::new();
    let mut used = vec![false; m];
    let mut current: Vec<Option<usize>> = vec![None; k];

    fn recurse(
        slot: usize,
        k: usize,
        shortlists: &[Vec<(usize, f64)>],
        used: &mut [bool],
        current: &mut Vec<Option<usize>>,
        current_score: f64,
        results: &mut Vec<Assignment>,
    ) {
        if slot == k {
            results.push(Assignment { slot_to_candidate: current.clone(), total_score: current_score });
            return;
        }
        // Leave this slot unfilled.
        recurse(slot + 1, k, shortlists, used, current, current_score, results);
        // Or try each shortlisted candidate.
        for &(c, s) in &shortlists[slot] {
            if used[c] {
                continue;
            }
            used[c] = true;
            current[slot] = Some(c);
            recurse(slot + 1, k, shortlists, used, current, current_score + s, results);
            current[slot] = None;
            used[c] = false;
        }
    }

    recurse(0, k, &shortlists, &mut used, &mut current, 0.0, &mut results);
    results
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Junction {
    Cis,
    Trans,
}

impl Junction {
    pub fn as_str(&self) -> &'static str {
        match self {
            Junction::Cis => "cis",
            Junction::Trans => "trans",
        }
    }
}

/// Genomic gap between two non-overlapping regions on the same contig (0 if
/// they overlap - shouldn't happen for correctly-assigned distinct exons,
/// but never negative/panicking if it does).
fn gap_bp(a_start: i64, a_end: i64, b_start: i64, b_end: i64) -> i64 {
    if a_end <= b_start {
        b_start - a_end
    } else if b_end <= a_start {
        a_start - b_end
    } else {
        0
    }
}

/// Classify the junction between two assigned exons *for this species*:
/// different contigs or a strand switch is a hard trans call (physically
/// impossible within one contiguous transcript); same contig + same strand
/// falls back to a distance threshold.
pub fn classify_junction(
    a_contig: &str,
    a_start: i64,
    a_end: i64,
    a_strand: char,
    b_contig: &str,
    b_start: i64,
    b_end: i64,
    b_strand: char,
    distance_threshold: i64,
) -> (Junction, i64) {
    if a_contig != b_contig {
        return (Junction::Trans, -1); // -1 = not a meaningful same-molecule gap
    }
    if a_strand != b_strand {
        let gap = gap_bp(a_start, a_end, b_start, b_end);
        return (Junction::Trans, gap);
    }
    let gap = gap_bp(a_start, a_end, b_start, b_end);
    if gap > distance_threshold {
        (Junction::Trans, gap)
    } else {
        (Junction::Cis, gap)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn picks_best_scoring_injective_matching() {
        // 3 candidates, 2 slots. Candidate 0 is great for slot 0 only;
        // candidate 1 is decent for both; candidate 2 is great for slot 1.
        // Optimal: 0->slot0, 2->slot1 (10+10=20), not 1->both (impossible)
        // or 1 used for either slot alone (worse than 20).
        let scores = vec![
            vec![10.0, 1.0],
            vec![5.0, 5.0],
            vec![1.0, 10.0],
        ];
        let floors = vec![0.0, 0.0];
        let a = best_assignment(&scores, &floors);
        assert_eq!(a.slot_to_candidate, vec![Some(0), Some(2)]);
        assert_eq!(a.total_score, 20.0);
    }

    #[test]
    fn beam_offers_alternatives_not_just_the_greedy_optimum() {
        // 3 candidates, 3 slots, mirroring the real nad2 misassignment
        // this was built to catch: per-slot-independent scoring alone
        // strictly prefers assigning 0->slot0, 1->slot1, 2->slot2 (sum
        // 10+10+10=30), but the "reversed" assignment 0->slot0, 2->slot1,
        // 1->slot2 (sum 10+9+9=28) is also a real, injective possibility
        // that a beam_width=2 search must surface for a caller (here,
        // whole-gene alignment) to be able to choose instead.
        let scores = vec![
            vec![10.0, 1.0, 1.0],
            vec![1.0, 10.0, 9.0],
            vec![1.0, 9.0, 10.0],
        ];
        let floors = vec![0.0, 0.0, 0.0];
        let beam = best_assignments_beam(&scores, &floors, 2);
        let greedy = vec![Some(0), Some(1), Some(2)];
        let reversed = vec![Some(0), Some(2), Some(1)];
        assert!(beam.iter().any(|a| a.slot_to_candidate == greedy), "greedy optimum must be in the beam");
        assert!(beam.iter().any(|a| a.slot_to_candidate == reversed), "the real alternative must be in the beam too");
    }

    #[test]
    fn slot_left_unfilled_if_nothing_clears_the_floor() {
        let scores = vec![vec![1.0, 1.0]];
        let floors = vec![5.0, 5.0];
        let a = best_assignment(&scores, &floors);
        assert_eq!(a.slot_to_candidate, vec![None, None]);
        assert_eq!(a.total_score, 0.0);
    }

    #[test]
    fn strand_switch_is_hard_trans_regardless_of_distance() {
        let (j, _) = classify_junction("ctg1", 100, 200, '+', "ctg1", 250, 300, '-', 10_000);
        assert_eq!(j, Junction::Trans);
    }

    #[test]
    fn nearby_same_strand_is_cis() {
        let (j, gap) = classify_junction("ctg1", 100, 200, '+', "ctg1", 300, 400, '+', 10_000);
        assert_eq!(j, Junction::Cis);
        assert_eq!(gap, 100);
    }

    #[test]
    fn far_same_strand_is_trans() {
        let (j, gap) = classify_junction("ctg1", 100, 200, '+', "ctg1", 500_000, 500_100, '+', 10_000);
        assert_eq!(j, Junction::Trans);
        assert_eq!(gap, 499_800);
    }

    #[test]
    fn different_contig_is_trans() {
        let (j, gap) = classify_junction("ctg1", 100, 200, '+', "ctg2", 100, 200, '+', 10_000);
        assert_eq!(j, Junction::Trans);
        assert_eq!(gap, -1);
    }
}
