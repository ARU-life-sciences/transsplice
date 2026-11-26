//! Loading a gene's reconstruction template: one profile+floor per exon
//! slot, plus a whole-gene profile for the final aggregate pass. Built by
//! `analysis/trans_splicing/src/00_build_exon_profiles.py` from verified
//! GenBank reference exons (see that script and `reference/raw/`).
//!
//! Directory layout expected (one per gene):
//!   exon1.pssm, exon1.threshold, [exon1.frame], exon2.pssm, ...
//!   whole_gene.pssm
//!
//! `exonN.frame` (`lead<TAB>trail`) says how many bases at the exon's 5'
//! end complete a codon begun in the previous exon, and how many at its 3'
//! end begin one the next exon completes - the junction phase. Slot
//! profiles cover only the exon's full codons, so these split-codon bases
//! are exactly what a codon-aligned exon match leaves out. Absent => 0/0.
//!
//! `junctions.tsv` (optional; `exon<TAB>n_cis<TAB>n_trans`, a header line
//! allowed) says, for the junction from exon N to exon N+1, in how many of
//! the template's references it is cis and in how many trans. A junction
//! that is cis in nearly all of them is unlikely to be trans in the genome
//! at hand (see `set_trans_junction_penalty`). Absent => nothing known.

use anyhow::{Context, Result};
use orfedit::profile::Profile;
use std::fs;
use std::path::Path;

pub struct SlotTemplate {
    pub profile: Profile,
    pub min_self_score: f64,
    pub lead: usize,
    pub trail: usize,
}

pub struct GeneTemplate {
    pub gene: String,
    pub slots: Vec<SlotTemplate>,
    pub whole_gene_profile: Profile,
    /// `junction_cis[i]`: the share of the references in which the junction
    /// from slot i+1 to slot i+2 is cis, if `junctions.tsv` gives it.
    pub junction_cis: Vec<Option<f64>>,
}

impl GeneTemplate {
    pub fn load(dir: &Path, gene: &str) -> Result<Self> {
        let mut slot_indices: Vec<usize> = Vec::new();
        for entry in fs::read_dir(dir).with_context(|| format!("reading template dir {}", dir.display()))? {
            let entry = entry?;
            let name = entry.file_name();
            let name = name.to_string_lossy();
            if let Some(rest) = name.strip_prefix("exon").and_then(|s| s.strip_suffix(".pssm")) {
                if let Ok(n) = rest.parse::<usize>() {
                    slot_indices.push(n);
                }
            }
        }
        slot_indices.sort_unstable();
        if slot_indices.is_empty() {
            anyhow::bail!("{}: no exonN.pssm files found", dir.display());
        }

        let mut slots = Vec::with_capacity(slot_indices.len());
        for n in &slot_indices {
            let pssm_path = dir.join(format!("exon{n}.pssm"));
            let threshold_path = dir.join(format!("exon{n}.threshold"));
            let profile = Profile::load(&pssm_path, &format!("{gene}_exon{n}"))?;
            let min_self_score: f64 = fs::read_to_string(&threshold_path)
                .with_context(|| format!("reading {}", threshold_path.display()))?
                .trim()
                .parse()
                .with_context(|| format!("bad threshold in {}", threshold_path.display()))?;
            let frame_path = dir.join(format!("exon{n}.frame"));
            let (lead, trail) = if frame_path.exists() {
                let text = fs::read_to_string(&frame_path)
                    .with_context(|| format!("reading {}", frame_path.display()))?;
                let v: Vec<usize> = text
                    .split_whitespace()
                    .map(|x| x.parse())
                    .collect::<Result<_, _>>()
                    .with_context(|| format!("bad lead/trail in {}", frame_path.display()))?;
                match v.as_slice() {
                    [l, t] if *l < 3 && *t < 3 => (*l, *t),
                    _ => anyhow::bail!("{}: expected `lead<TAB>trail`, each 0-2", frame_path.display()),
                }
            } else {
                (0, 0)
            };
            slots.push(SlotTemplate { profile, min_self_score, lead, trail });
        }

        let whole_gene_profile = Profile::load(&dir.join("whole_gene.pssm"), gene)?;

        let junctions_path = dir.join("junctions.tsv");
        let junction_cis = if junctions_path.exists() {
            let text = fs::read_to_string(&junctions_path).with_context(|| format!("reading {}", junctions_path.display()))?;
            junction_shares(&text, &slot_indices).with_context(|| format!("in {}", junctions_path.display()))?
        } else {
            vec![None; slot_indices.len().saturating_sub(1)]
        };

        Ok(GeneTemplate { gene: gene.to_string(), slots, whole_gene_profile, junction_cis })
    }

    /// The share of the references in which exons in slots `from` and `to`
    /// (1-indexed, from < to) are joined cis: with slots between them left
    /// empty, the least of the junctions spanned (one often trans among
    /// them explains the gap). None if any of those is not known.
    pub fn cis_share(&self, from: usize, to: usize) -> Option<f64> {
        if from == 0 || to <= from || to - 1 > self.junction_cis.len() {
            return None;
        }
        self.junction_cis[from - 1..to - 1].iter().try_fold(1.0f64, |m, x| x.map(|f| m.min(f)))
    }
}

/// `junctions.tsv`'s text -> the cis share of each junction between the
/// slots (exon numbers `slot_indices`, sorted): entry i for slot i+1 to i+2.
pub fn junction_shares(text: &str, slot_indices: &[usize]) -> Result<Vec<Option<f64>>> {
    let mut by_exon: std::collections::HashMap<usize, f64> = std::collections::HashMap::new();
    for line in text.lines().map(str::trim).filter(|l| !l.is_empty() && !l.starts_with('#')) {
        let fields: Vec<&str> = line.split('\t').collect();
        let Ok(exon) = fields[0].parse::<usize>() else { continue }; // the header
        let counts: Vec<f64> =
            fields[1..].iter().map(|x| x.parse()).collect::<Result<_, _>>().with_context(|| format!("bad counts: {line}"))?;
        let [n_cis, n_trans] = counts[..] else {
            anyhow::bail!("expected `exon<TAB>n_cis<TAB>n_trans`, got {line}");
        };
        if n_cis + n_trans > 0.0 {
            by_exon.insert(exon, n_cis / (n_cis + n_trans));
        }
    }
    // between two slots whose exon numbers are not consecutive (a slot
    // with too few references was left out of the template), the least
    // of the junctions spanned
    Ok(slot_indices
        .windows(2)
        .map(|w| (w[0]..w[1]).try_fold(1.0f64, |m, n| by_exon.get(&n).map(|&f| m.min(f))))
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn junction_shares_are_read_per_junction() {
        let text = "exon\tn_cis\tn_trans\n1\t2\t19\n2\t21\t0\n4\t0\t0\n";
        let shares = junction_shares(text, &[1, 2, 3, 4, 5]).unwrap();
        assert_eq!(shares.len(), 4);
        assert!((shares[0].unwrap() - 2.0 / 21.0).abs() < 1e-12);
        assert_eq!(shares[1], Some(1.0));
        assert_eq!(shares[2], None); // not given
        assert_eq!(shares[3], None); // no references
        assert!(junction_shares("1\t2\n", &[1, 2]).is_err());
        // exon 2's slot missing from the template: 1 -> 3 spans two junctions
        assert_eq!(junction_shares("1\t3\t1\n2\t4\t0\n", &[1, 3]).unwrap(), vec![Some(0.75)]);
    }

    #[test]
    fn a_gap_over_an_empty_slot_takes_the_least_cis_junction() {
        let profile = || Profile { gene: "x".to_string(), columns: vec![[0.0; 20]], edited: None };
        let slot = || SlotTemplate { profile: profile(), min_self_score: 0.0, lead: 0, trail: 0 };
        let template = GeneTemplate {
            gene: "x".to_string(),
            slots: vec![slot(), slot(), slot()],
            whole_gene_profile: profile(),
            junction_cis: vec![Some(1.0), Some(0.25)],
        };
        assert_eq!(template.cis_share(1, 2), Some(1.0));
        assert_eq!(template.cis_share(2, 3), Some(0.25));
        assert_eq!(template.cis_share(1, 3), Some(0.25));
        assert_eq!(template.cis_share(2, 4), None);
        let unknown = GeneTemplate { junction_cis: vec![Some(1.0), None], ..template };
        assert_eq!(unknown.cis_share(1, 3), None);
    }
}
