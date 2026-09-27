//! Loading a gene's reconstruction template: one profile+floor per exon
//! slot, plus a whole-gene profile for the final aggregate pass. Built by
//! `analysis/trans_splicing/src/00_build_exon_profiles.py` from verified
//! GenBank reference exons (see that script and `reference/raw/`).
//!
//! Directory layout expected (one per gene):
//!   exon1.pssm, exon1.threshold, exon2.pssm, exon2.threshold, ...
//!   whole_gene.pssm

use anyhow::{Context, Result};
use orfedit::profile::Profile;
use std::fs;
use std::path::Path;

pub struct SlotTemplate {
    pub profile: Profile,
    pub min_self_score: f64,
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
            slots.push(SlotTemplate { profile, min_self_score });
        }

        let whole_gene_profile = Profile::load(&dir.join("whole_gene.pssm"), gene)?;

        Ok(GeneTemplate { gene: gene.to_string(), slots, whole_gene_profile })
    }
}
