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
use transsplice::{denovo_scan, finish_gene, prepare_gene, DenovoScan, DenovoTask, Prepared};

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
        /// Also write one row per frameshift found in a reconstructed gene
        #[arg(long)]
        frameshifts_out: Option<PathBuf>,
        #[arg(long, default_value_t = 4)]
        threads: usize,
        #[arg(long, default_value_t = 8.0)]
        gap_open: f64,
        #[arg(long, default_value_t = 1.0)]
        gap_extend: f64,
        #[arg(long, default_value_t = 2.5)]
        edit_penalty: f64,
        /// Cost of a U-to-C ("reverse") RNA edit; "inf", the default, for
        /// C-to-U editing only. Give a finite cost for ferns, lycophytes
        /// and hornworts (see orfedit)
        #[arg(long, default_value_t = f64::INFINITY)]
        reverse_edit_penalty: f64,
        /// Cost of reading a genomic in-frame stop codon as Gln or Arg by a
        /// U-to-C edit; "inf", the default, for none (see orfedit)
        #[arg(long, default_value_t = f64::INFINITY)]
        stop_edit_penalty: f64,
        #[arg(long, default_value_t = 5000)]
        merge_distance: i64,
        #[arg(long, default_value_t = 10000)]
        junction_distance_threshold: i64,
        #[arg(long, default_value_t = 60)]
        exon_flank: i64,
        /// Cost of a frame change (skipping 1 or 2 bases) inside a gene;
        /// "inf" keeps every exon and the joined CDS in one reading frame
        #[arg(long, default_value_t = 20.0)]
        frameshift_penalty: f64,
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

#[derive(Clone, Hash, Eq, PartialEq, PartialOrd, Ord)]
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
    phase: usize,
}

struct EditRow {
    species: String,
    organelle: String,
    gene: String,
    slot: usize,
    genomic_pos: usize,
    profile_column: usize,
    resulting_aa: char,
    reverse: bool,
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
            reverse_edit_penalty,
            stop_edit_penalty,
            merge_distance,
            junction_distance_threshold,
            exon_flank,
            frameshifts_out,
            frameshift_penalty,
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

            let params = ScanParams { gap_open, gap_extend, edit_penalty, reverse_edit_penalty, stop_edit_penalty, flank: 0, frameshift_penalty, ..ScanParams::default() };
            let mut group_vec: Vec<(&GroupKey, &Vec<&ManifestRow>)> = groups.iter().collect();
            group_vec.sort_by(|a, b| a.0.cmp(b.0)); // stable output row order

            // Three stages, each a flat pool of parallel tasks (see
            // transsplice::prepare_gene): score and assign every gene's
            // candidates; run every de novo scan, largest first; finish
            // every gene. One task per item throughout: costs differ by
            // orders of magnitude, and a multi-item chunk can't be stolen
            // once a thread has started it.
            let prepared: Vec<Option<(Prepared, Vec<DenovoTask>)>> = group_vec
                .par_iter()
                .with_max_len(1)
                .map(|(key, frags)| {
                    let template = templates.get(&key.template_dir)?;
                    let fasta_idx = fastas.get(&key.fasta)?;
                    let fragments: Vec<RawFragment> = frags
                        .iter()
                        .map(|f| RawFragment { contig: f.contig.clone(), start: f.start, end: f.end, strand: f.strand })
                        .collect();
                    Some(prepare_gene(fasta_idx, &fragments, template, &params, merge_distance, exon_flank))
                })
                .collect();
            let (prepared, tasks): (Vec<Option<Prepared>>, Vec<Vec<DenovoTask>>) = prepared
                .into_iter()
                .map(|p| match p {
                    Some((prepared, tasks)) => (Some(prepared), tasks),
                    None => (None, Vec::new()),
                })
                .unzip();
            let mut jobs: Vec<(u64, usize, usize)> = Vec::new(); // (cells, gene, task)
            for (gi, gene_tasks) in tasks.iter().enumerate() {
                let template = &templates[&group_vec[gi].0.template_dir];
                for (ti, t) in gene_tasks.iter().enumerate() {
                    jobs.push((t.cells(template), gi, ti));
                }
            }
            jobs.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
            let scanned: Vec<(usize, usize, DenovoScan)> = jobs
                .par_iter()
                .with_max_len(1)
                .map(|&(_, gi, ti)| {
                    let key = group_vec[gi].0;
                    (gi, ti, denovo_scan(&fastas[&key.fasta], &templates[&key.template_dir], &params, &tasks[gi][ti]))
                })
                .collect();
            let mut scans: Vec<Vec<DenovoScan>> = tasks.iter().map(|t| vec![None; t.len()]).collect();
            for (gi, ti, scan) in scanned {
                scans[gi][ti] = scan;
            }
            let work: Vec<(&GroupKey, Prepared, Vec<DenovoTask>, Vec<DenovoScan>)> = group_vec
                .iter()
                .zip(prepared)
                .zip(tasks)
                .zip(scans)
                .filter_map(|((((key, _), prepared), tasks), scans)| Some((*key, prepared?, tasks, scans)))
                .collect();

            let results: Vec<(GeneRow, Vec<ExonRow>, Vec<JunctionRow>, Vec<EditRow>, Vec<String>)> = work
                .into_par_iter()
                .with_max_len(1)
                .filter_map(|(key, prepared, tasks, scans)| {
                    let template = templates.get(&key.template_dir)?;
                    let fasta_idx = fastas.get(&key.fasta)?;
                    match finish_gene(
                        fasta_idx, prepared, &tasks, &scans, template, &params,
                        junction_distance_threshold, exon_flank,
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
                                phase: e.phase,
                            }).collect();
                            let junction_rows = outcome.junctions.iter().map(|j| JunctionRow {
                                species: key.species.clone(), organelle: key.organelle.clone(), gene: key.gene.clone(),
                                from_slot: j.from_slot, to_slot: j.to_slot, kind: j.kind, gap_bp: j.gap_bp,
                            }).collect();
                            let edit_rows = outcome.exons.iter().flat_map(|e| {
                                e.edits.iter().map(move |ed| EditRow {
                                    species: key.species.clone(), organelle: key.organelle.clone(), gene: key.gene.clone(),
                                    slot: e.slot, genomic_pos: ed.genomic_pos, profile_column: ed.profile_column,
                                    resulting_aa: ed.resulting_aa, reverse: ed.reverse,
                                })
                            }).collect();
                            let shift_rows = outcome.frameshifts.iter().map(|f| format!(
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                                key.species, key.organelle, key.gene, f.slot, f.contig, f.start, f.end, f.n_bases,
                                f.profile_column, f.window_start, f.window_end
                            )).collect();
                            Some((gene_row, exon_rows, junction_rows, edit_rows, shift_rows))
                        }
                        Err(e) => {
                            eprintln!("[warn] transsplice: {} {} {}: {:#}", key.species, key.organelle, key.gene, e);
                            None
                        }
                    }
                })
                .collect();

            write_outputs(&genes_out, &exons_out, &junctions_out, &edits_out, &results)?;
            if let Some(path) = &frameshifts_out {
                let mut w = String::from(
                    "species\torganelle\tgene\tslot\tcontig\tstart\tend\tn_bases\tprofile_column\twindow_start\twindow_end\n",
                );
                for (_, _, _, _, rows) in &results {
                    rows.iter().for_each(|r| w.push_str(r));
                }
                if let Some(parent) = path.parent() {
                    fs::create_dir_all(parent)?;
                }
                let tmp = path.with_extension("tsv.tmp");
                fs::write(&tmp, w)?;
                fs::rename(&tmp, path)?;
            }
            let n_complete = results.iter().filter(|(g, ..)| g.complete).count();
            eprintln!(
                "[info] transsplice: {} genes reconstructed, {} complete, {} with frameshifts -> {}",
                results.len(), n_complete, results.iter().filter(|r| !r.4.is_empty()).count(), genes_out.display()
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
    results: &[(GeneRow, Vec<ExonRow>, Vec<JunctionRow>, Vec<EditRow>, Vec<String>)],
) -> Result<()> {
    for p in [genes_out, exons_out, junctions_out, edits_out] {
        if let Some(parent) = p.parent() {
            fs::create_dir_all(parent)?;
        }
    }

    let mut g = String::from("species\torganelle\tgene\tn_slots\tn_filled\tcomplete\twhole_gene_score\twhole_gene_n_edits\tprotein\n");
    for (row, ..) in results {
        g.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\n",
            row.species, row.organelle, row.gene, row.n_slots, row.n_filled, row.complete,
            row.whole_gene_score, row.whole_gene_n_edits, row.protein
        ));
    }
    let tmp = genes_out.with_extension("tsv.tmp");
    fs::write(&tmp, g)?;
    fs::rename(&tmp, genes_out)?;

    let mut e = String::from("species\torganelle\tgene\tslot\tcontig\tstart\tend\tstrand\tscore\tn_edits\tvia_denovo\tprovenance\tphase\n");
    for (_, exon_rows, ..) in results {
        for row in exon_rows {
            e.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\n",
                row.species, row.organelle, row.gene, row.slot, row.contig, row.start, row.end, row.strand,
                row.score, row.n_edits, row.via_denovo, row.provenance, row.phase
            ));
        }
    }
    let tmp = exons_out.with_extension("tsv.tmp");
    fs::write(&tmp, e)?;
    fs::rename(&tmp, exons_out)?;

    let mut j = String::from("species\torganelle\tgene\tfrom_slot\tto_slot\tkind\tgap_bp\n");
    for (_, _, junction_rows, ..) in results {
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

    let mut ed = String::from("species\torganelle\tgene\tslot\tgenomic_pos\tprofile_column\tresulting_aa\tedit\n");
    for (_, _, _, edit_rows, _) in results {
        for row in edit_rows {
            ed.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                row.species, row.organelle, row.gene, row.slot, row.genomic_pos, row.profile_column, row.resulting_aa,
                if row.reverse { "U>C" } else { "C>U" }
            ));
        }
    }
    let tmp = edits_out.with_extension("tsv.tmp");
    fs::write(&tmp, ed)?;
    fs::rename(&tmp, edits_out)?;

    Ok(())
}
