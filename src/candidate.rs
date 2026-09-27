//! Candidate exon fragments: the raw `.ctg.bed` hits for one (species, gene)
//! pair, plus merged-run candidates for the case where oatk's HMM splits one
//! real exon into two-or-more adjacent hits (confirmed in real data - see
//! transsplice's README).

use std::collections::HashMap;

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
pub fn build_candidates(fragments: &[RawFragment], merge_distance: i64) -> Vec<Candidate> {
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
        let cands = build_candidates(&frags, 5000);
        assert_eq!(cands.len(), 3); // 2 originals + 1 merged
        let merged = cands.iter().find(|c| c.provenance.starts_with("merged:")).unwrap();
        assert_eq!(merged.start, 670208);
        assert_eq!(merged.end, 671964);
    }

    #[test]
    fn distant_fragments_are_not_merged() {
        let frags = vec![frag("ctg1", 1000, 1200, '+'), frag("ctg1", 500000, 500200, '+')];
        let cands = build_candidates(&frags, 5000);
        assert_eq!(cands.len(), 2); // no merge candidate added
    }

    #[test]
    fn different_strand_fragments_are_never_merged() {
        let frags = vec![frag("ctg1", 1000, 1200, '+'), frag("ctg1", 1300, 1500, '-')];
        let cands = build_candidates(&frags, 5000);
        assert_eq!(cands.len(), 2);
    }
}
