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

        Ok(GeneTemplate { gene: gene.to_string(), slots, whole_gene_profile })
    }
}
