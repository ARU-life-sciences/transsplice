//! Candidate exon fragments: the raw `.ctg.bed` hits for one (species, gene)
//! pair, plus merged-run candidates for the case where oatk's HMM splits one
//! real exon into two-or-more adjacent hits (confirmed in real data - see
//! transsplice's README).

use std::collections::HashMap;

/// How near (bp) a circular contig's two ends two fragments must come to be
/// taken as one exon cut by the origin (an HMM hit stops at the end of the
/// sequence it is given).
const ORIGIN_TOUCH: i64 = 30;

#[derive(Clone, Debug)]
pub struct RawFragment {
    pub contig: String,
    pub start: i64,
    pub end: i64,
    pub strand: char,
}

#[derive(Clone, Debug)]
pub struct Candidate {
    pub contig: String,
    pub start: i64,
    pub end: i64,
    pub strand: char,
    /// Human-readable provenance, e.g. "90240-90313" or a merged run
    /// "merged:670208-670631+671812-671964".
    pub provenance: String,
}

/// Every original fragment, plus one extra candidate per maximal run of
/// same-contig, same-strand fragments whose consecutive gaps are all
/// <= `merge_distance` (bp). Both the pieces and the merged span are kept -
/// scoring against each exon-slot profile decides which represents the
/// real exon (see README: a wrongly-merged span reads as a worse match to
/// any single slot than the correct pieces do to their own distinct slots).
/// On a circular contig (`circular`: its own length), a fragment that runs
/// to the contig's end and one that starts at its beginning (within
/// ORIGIN_TOUCH) are one exon cut by the origin (fern ndhB exon 2): they
/// are merged across it, the span ending past the contig's length. Only
/// those two: merging the whole last run with the whole first one made
/// one candidate of every exon of a liverwort cox1 lying either side of
/// the origin, which crowded an exon's own candidate out of its slot, and
/// moved plastid rpl2 onto the copy of the inverted repeat across the
/// origin, whose intron there was already measured round it.
pub fn build_candidates(fragments: &[RawFragment], merge_distance: i64, circular: &HashMap<String, i64>) -> Vec<Candidate> {
    let mut out: Vec<Candidate> = fragments
        .iter()
        .map(|f| Candidate {
            contig: f.contig.clone(),
            start: f.start,
            end: f.end,
            strand: f.strand,
            provenance: format!("{}-{}", f.start, f.end),
        })
        .collect();

    let mut groups: HashMap<(String, char), Vec<&RawFragment>> = HashMap::new();
    for f in fragments {
        groups.entry((f.contig.clone(), f.strand)).or_default().push(f);
    }

    // HashMap order is random per process: fix the order of merged-run
    // candidates so equal-scoring ties resolve the same way every run.
    let mut groups: Vec<((String, char), Vec<&RawFragment>)> = groups.into_iter().collect();
    groups.sort_by(|a, b| a.0.cmp(&b.0));
    for ((contig, strand), mut frags) in groups {
        frags.sort_by_key(|f| f.start);
        let mut i = 0;
        while i < frags.len() {
            let mut j = i;
            while j + 1 < frags.len() && frags[j + 1].start - frags[j].end <= merge_distance {
                j += 1;
            }
            if j > i {
                let span_start = frags[i].start;
                let span_end = frags[i..=j].iter().map(|f| f.end).max().unwrap();
                let provenance = frags[i..=j]
                    .iter()
                    .map(|f| format!("{}-{}", f.start, f.end))
                    .collect::<Vec<_>>()
                    .join("+");
                out.push(Candidate {
                    contig: contig.clone(),
                    start: span_start,
                    end: span_end,
                    strand,
                    provenance: format!("merged:{}", provenance),
                });
            }
            i = j + 1;
        }
        // one exon cut by a circular contig's origin
        if let Some(&n) = circular.get(&contig) {
            let last = frags.iter().max_by_key(|f| f.end).unwrap();
            let first = frags[0];
            if !std::ptr::eq(*last, first) && n - last.end <= ORIGIN_TOUCH && first.start <= ORIGIN_TOUCH {
                out.push(Candidate {
                    contig: contig.clone(),
                    start: last.start,
                    end: n + first.end,
                    strand,
                    provenance: format!("merged:{}-{}+{}-{}", last.start, last.end, first.start, first.end),
                });
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frag(contig: &str, start: i64, end: i64, strand: char) -> RawFragment {
        RawFragment { contig: contig.to_string(), start, end, strand }
    }

    #[test]
    fn adjacent_same_strand_fragments_get_a_merged_candidate_too() {
        let frags = vec![frag("ctg1", 670208, 670631, '-'), frag("ctg1", 671812, 671964, '-')];
        let cands = build_candidates(&frags, 5000, &HashMap::new());
        assert_eq!(cands.len(), 3); // 2 originals + 1 merged
        let merged = cands.iter().find(|c| c.provenance.starts_with("merged:")).unwrap();
        assert_eq!(merged.start, 670208);
        assert_eq!(merged.end, 671964);
    }

    #[test]
    fn distant_fragments_are_not_merged() {
        let frags = vec![frag("ctg1", 1000, 1200, '+'), frag("ctg1", 500000, 500200, '+')];
        let cands = build_candidates(&frags, 5000, &HashMap::new());
        assert_eq!(cands.len(), 2); // no merge candidate added
    }

    #[test]
    fn fragments_either_side_of_a_circular_origin_get_a_merged_candidate() {
        // an exon across the origin of a 148,588 bp circular contig
        let frags = vec![frag("ctg1", 0, 110, '-'), frag("ctg1", 945, 1662, '-'), frag("ctg1", 147917, 148588, '-')];
        let circular: HashMap<String, i64> = [("ctg1".to_string(), 148588)].into_iter().collect();
        let cands = build_candidates(&frags, 500, &circular);
        let merged: Vec<&Candidate> = cands.iter().filter(|c| c.provenance.starts_with("merged:")).collect();
        assert_eq!(merged.len(), 1);
        assert_eq!((merged[0].start, merged[0].end), (147917, 148588 + 110));
        // linear: nothing across the ends
        assert!(build_candidates(&frags, 500, &HashMap::new()).iter().all(|c| !c.provenance.starts_with("merged:")));
        // two exons either side of the origin, an intron between: not merged
        let frags = vec![frag("ctg1", 775, 1209, '+'), frag("ctg1", 146946, 147337, '+')];
        assert!(build_candidates(&frags, 5000, &circular).iter().all(|c| !c.provenance.starts_with("merged:")));
    }

    #[test]
    fn different_strand_fragments_are_never_merged() {
        let frags = vec![frag("ctg1", 1000, 1200, '+'), frag("ctg1", 1300, 1500, '-')];
        let cands = build_candidates(&frags, 5000, &HashMap::new());
        assert_eq!(cands.len(), 2);
    }
}
