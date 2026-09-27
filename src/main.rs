//! `transsplice` - reconstructs trans-spliced plant mitochondrial genes
//! (nad1/nad2/nad5/rps3) from oatk's scattered per-exon `.ctg.bed` hits.
//! See README.md and `transsplice`'s design notes for the algorithm; reuses
//! `orfedit`'s edit-tolerant alignment engine as a library throughout.

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use orfedit::fasta::FastaIndex;
use orfedit::ScanParams;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;
use transsplice::candidate::RawFragment;
use transsplice::template::GeneTemplate;
use transsplice::reconstruct_gene;

#[derive(Parser)]
#[command(name = "transsplice", version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Reconstruct every (species, gene) group in a manifest TSV in parallel.
    ScanBatch {
        /// TSV, header: species organelle gene fasta template_dir contig start end strand
        /// (one row per raw `.ctg.bed` fragment - multiple rows share the
        /// same species+organelle+gene+fasta+template_dir key).
        #[arg(long)]
        manifest: PathBuf,
        #[arg(long)]
        genes_out: PathBuf,
        #[arg(long)]
        exons_out: PathBuf,
        #[arg(long)]
        junctions_out: PathBuf,
        #[arg(long)]
        edits_out: PathBuf,
        #[arg(long, default_value_t = 4)]
        threads: usize,
        #[arg(long, default_value_t = 8.0)]
        gap_open: f64,
        #[arg(long, default_value_t = 1.0)]
        gap_extend: f64,
        #[arg(long, default_value_t = 2.5)]
        edit_penalty: f64,
        #[arg(long, default_value_t = 5000)]
        merge_distance: i64,
        #[arg(long, default_value_t = 10000)]
        junction_distance_threshold: i64,
        #[arg(long, default_value_t = 60)]
        exon_flank: i64,
    },
}

struct ManifestRow {
    species: String,
    organelle: String,
    gene: String,
    fasta: String,
    template_dir: String,
    contig: String,
    start: i64,
    end: i64,
    strand: char,
}

fn parse_manifest(path: &PathBuf) -> Result<Vec<ManifestRow>> {
    let text = fs::read_to_string(path).with_context(|| format!("reading manifest {}", path.display()))?;
    let mut lines = text.lines();
    let header = lines.next().context("empty manifest")?;
    let cols: Vec<&str> = header.split('\t').collect();
    let expected = ["species", "organelle", "gene", "fasta", "template_dir", "contig", "start", "end", "strand"];
    if cols != expected {
        anyhow::bail!("manifest header must be {:?}, got {:?}", expected, cols);
    }
    let mut rows = Vec::new();
    for (lineno, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != 9 {
            anyhow::bail!("{}:{}: expected 9 columns, got {}", path.display(), lineno + 2, f.len());
        }
        rows.push(ManifestRow {
            species: f[0].to_string(),
            organelle: f[1].to_string(),
            gene: f[2].to_string(),
            fasta: f[3].to_string(),
            template_dir: f[4].to_string(),
            contig: f[5].to_string(),
            start: f[6].parse().with_context(|| format!("{}:{}: bad start", path.display(), lineno + 2))?,
            end: f[7].parse().with_context(|| format!("{}:{}: bad end", path.display(), lineno + 2))?,
            strand: f[8].chars().next().context("empty strand")?,
        });
    }
    Ok(rows)
}

#[derive(Clone, Hash, Eq, PartialEq)]
struct GroupKey {
    species: String,
    organelle: String,
    gene: String,
    fasta: String,
    template_dir: String,
}

struct GeneRow {
    species: String,
    organelle: String,
    gene: String,
    n_slots: usize,
    n_filled: usize,
    complete: bool,
    whole_gene_score: f64,
    whole_gene_n_edits: usize,
    protein: String,
}

struct ExonRow {
    species: String,
    organelle: String,
    gene: String,
    slot: usize,
    contig: String,
    start: usize,
    end: usize,
    strand: char,
    score: f64,
    n_edits: usize,
    via_denovo: bool,
    provenance: String,
}

struct EditRow {
    species: String,
    organelle: String,
    gene: String,
    slot: usize,
    genomic_pos: usize,
    profile_column: usize,
    resulting_aa: char,
}

struct JunctionRow {
    species: String,
    organelle: String,
    gene: String,
    from_slot: usize,
    to_slot: usize,
    kind: &'static str,
    gap_bp: i64,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::ScanBatch {
            manifest,
            genes_out,
            exons_out,
            junctions_out,
            edits_out,
            threads,
            gap_open,
            gap_extend,
            edit_penalty,
            merge_distance,
            junction_distance_threshold,
            exon_flank,
        } => {
            rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().ok();

            let rows = parse_manifest(&manifest)?;
            eprintln!("[info] transsplice: {} manifest rows", rows.len());

            let mut groups: HashMap<GroupKey, Vec<&ManifestRow>> = HashMap::new();
            for r in &rows {
                groups
                    .entry(GroupKey {
                        species: r.species.clone(),
                        organelle: r.organelle.clone(),
                        gene: r.gene.clone(),
                        fasta: r.fasta.clone(),
                        template_dir: r.template_dir.clone(),
                    })
                    .or_default()
                    .push(r);
            }

            // Cache templates (small: 4 genes total) and fastas (one load
            // per distinct species+organelle fasta path, shared across
            // whichever genes reference it) up front so the parallel
            // section below only reads, never races on first-load.
            let mut template_paths: Vec<&str> = groups.keys().map(|k| k.template_dir.as_str()).collect();
            template_paths.sort_unstable();
            template_paths.dedup();
            let mut templates: HashMap<String, GeneTemplate> = HashMap::new();
            for tdir in template_paths {
                let gene = groups.keys().find(|k| k.template_dir == tdir).map(|k| k.gene.clone()).unwrap_or_default();
                match GeneTemplate::load(&PathBuf::from(tdir), &gene) {
                    Ok(t) => {
                        templates.insert(tdir.to_string(), t);
                    }
                    Err(e) => eprintln!("[warn] transsplice: failed to load template {}: {:#}", tdir, e),
                }
            }

            let mut fasta_paths: Vec<&str> = groups.keys().map(|k| k.fasta.as_str()).collect();
            fasta_paths.sort_unstable();
            fasta_paths.dedup();
            let mut fastas: HashMap<String, FastaIndex> = HashMap::new();
            for fpath in fasta_paths {
                match FastaIndex::load(&PathBuf::from(fpath)) {
                    Ok(idx) => {
                        fastas.insert(fpath.to_string(), idx);
                    }
                    Err(e) => eprintln!("[warn] transsplice: failed to load fasta {}: {:#}", fpath, e),
                }
            }

            let params = ScanParams { gap_open, gap_extend, edit_penalty, flank: 0 };
            let group_vec: Vec<(&GroupKey, &Vec<&ManifestRow>)> = groups.iter().collect();

            let results: Vec<(GeneRow, Vec<ExonRow>, Vec<JunctionRow>, Vec<EditRow>)> = group_vec
                .par_iter()
                .filter_map(|(key, frags)| {
                    let template = templates.get(&key.template_dir)?;
                    let fasta_idx = fastas.get(&key.fasta)?;
                    let fragments: Vec<RawFragment> = frags
                        .iter()
                        .map(|f| RawFragment { contig: f.contig.clone(), start: f.start, end: f.end, strand: f.strand })
                        .collect();
                    match reconstruct_gene(
                        fasta_idx, &fragments, template, &params,
                        merge_distance, junction_distance_threshold, exon_flank,
                    ) {
                        Ok(outcome) => {
                            let gene_row = GeneRow {
                                species: key.species.clone(), organelle: key.organelle.clone(), gene: key.gene.clone(),
                                n_slots: outcome.n_slots, n_filled: outcome.n_filled, complete: outcome.complete,
                                whole_gene_score: outcome.whole_gene_score, whole_gene_n_edits: outcome.whole_gene_n_edits,
                                protein: outcome.protein,
                            };
                            let exon_rows = outcome.exons.iter().map(|e| ExonRow {
                                species: key.species.clone(), organelle: key.organelle.clone(), gene: key.gene.clone(),
                                slot: e.slot, contig: e.contig.clone(), start: e.start, end: e.end, strand: e.strand,
                                score: e.score, n_edits: e.edits.len(), via_denovo: e.via_denovo, provenance: e.provenance.clone(),
                            }).collect();
                            let junction_rows = outcome.junctions.iter().map(|j| JunctionRow {
                                species: key.species.clone(), organelle: key.organelle.clone(), gene: key.gene.clone(),
                                from_slot: j.from_slot, to_slot: j.to_slot, kind: j.kind, gap_bp: j.gap_bp,
                            }).collect();
                            let edit_rows = outcome.exons.iter().flat_map(|e| {
                                e.edits.iter().map(move |ed| EditRow {
                                    species: key.species.clone(), organelle: key.organelle.clone(), gene: key.gene.clone(),
                                    slot: e.slot, genomic_pos: ed.genomic_pos, profile_column: ed.profile_column,
                                    resulting_aa: ed.resulting_aa,
                                })
                            }).collect();
                            Some((gene_row, exon_rows, junction_rows, edit_rows))
                        }
                        Err(e) => {
                            eprintln!("[warn] transsplice: {} {} {}: {:#}", key.species, key.organelle, key.gene, e);
                            None
                        }
                    }
                })
                .collect();

            write_outputs(&genes_out, &exons_out, &junctions_out, &edits_out, &results)?;
            let n_complete = results.iter().filter(|(g, _, _, _)| g.complete).count();
            eprintln!(
                "[info] transsplice: {} genes reconstructed, {} complete -> {}",
                results.len(), n_complete, genes_out.display()
            );
            Ok(())
        }
    }
}

fn write_outputs(
    genes_out: &PathBuf,
    exons_out: &PathBuf,
    junctions_out: &PathBuf,
    edits_out: &PathBuf,
    results: &[(GeneRow, Vec<ExonRow>, Vec<JunctionRow>, Vec<EditRow>)],
) -> Result<()> {
    for p in [genes_out, exons_out, junctions_out, edits_out] {
        if let Some(parent) = p.parent() {
            fs::create_dir_all(parent)?;
        }
    }

    let mut g = String::from("species\torganelle\tgene\tn_slots\tn_filled\tcomplete\twhole_gene_score\twhole_gene_n_edits\tprotein\n");
    for (row, _, _, _) in results {
        g.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\n",
            row.species, row.organelle, row.gene, row.n_slots, row.n_filled, row.complete,
            row.whole_gene_score, row.whole_gene_n_edits, row.protein
        ));
    }
    let tmp = genes_out.with_extension("tsv.tmp");
    fs::write(&tmp, g)?;
    fs::rename(&tmp, genes_out)?;

    let mut e = String::from("species\torganelle\tgene\tslot\tcontig\tstart\tend\tstrand\tscore\tn_edits\tvia_denovo\tprovenance\n");
    for (_, exon_rows, _, _) in results {
        for row in exon_rows {
            e.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\n",
                row.species, row.organelle, row.gene, row.slot, row.contig, row.start, row.end, row.strand,
                row.score, row.n_edits, row.via_denovo, row.provenance
            ));
        }
    }
    let tmp = exons_out.with_extension("tsv.tmp");
    fs::write(&tmp, e)?;
    fs::rename(&tmp, exons_out)?;

    let mut j = String::from("species\torganelle\tgene\tfrom_slot\tto_slot\tkind\tgap_bp\n");
    for (_, _, junction_rows, _) in results {
        for row in junction_rows {
            j.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                row.species, row.organelle, row.gene, row.from_slot, row.to_slot, row.kind, row.gap_bp
            ));
        }
    }
    let tmp = junctions_out.with_extension("tsv.tmp");
    fs::write(&tmp, j)?;
    fs::rename(&tmp, junctions_out)?;

    let mut ed = String::from("species\torganelle\tgene\tslot\tgenomic_pos\tprofile_column\tresulting_aa\n");
    for (_, _, _, edit_rows) in results {
        for row in edit_rows {
            ed.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                row.species, row.organelle, row.gene, row.slot, row.genomic_pos, row.profile_column, row.resulting_aa
            ));
        }
    }
    let tmp = edits_out.with_extension("tsv.tmp");
    fs::write(&tmp, ed)?;
    fs::rename(&tmp, edits_out)?;

    Ok(())
}
