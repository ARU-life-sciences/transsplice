# transsplice

Reconstructs trans-spliced genes in plant mitochondrial genomes from
scattered per-exon gene-model hits, reusing
[`orfedit`](https://github.com/tolkit/orfedit)'s edit-tolerant alignment
engine.

## The problem

`nad1`, `nad2`, `nad5`, and (in some lineages) `rps3` are trans-spliced in
angiosperm mitochondria: their exons are transcribed as separate RNA
molecules from genomic loci that can be hundreds of kb apart, sometimes on
opposite strands, then spliced together at the RNA level. A targeted gene
caller (oatk, or any HMM/profile search) sees this as several independent
per-exon hits under the same gene name - nothing joins them into one gene.

## How it works

Given a gene's raw exon fragments (from any upstream caller) and a
reference "template" (one small profile per exon slot 1..K, plus a
whole-gene profile - see `Profile format` below and the per-slot
`.threshold` files), for one (species, gene):

1. **Merge** same-strand fragments within a configurable distance into an
   additional candidate spanning the whole run (a real gene model can
   split one true exon into adjacent hits) - both the pieces and the
   merged span are kept, scoring decides which is real.
2. **Score** every candidate against every exon-slot profile via
   `orfedit::scan_gene` reused as-is, so each score is already a
   boundary-corrected, edit-tolerant alignment, not a raw heuristic.
3. **De novo fallback**: if a slot has no candidate scoring above its own
   floor, scan the contig(s) the other exons were found on (bounded, not
   genome-wide) for that slot's profile directly.
4. **Assign** candidates to slots via the best-scoring injective partial
   matching (brute-force search - the exon count is always small).
5. **Concatenate** assigned exons in slot order, each independently
   strand-corrected, classify each junction cis/trans for this species
   specifically (never inherited from a reference topology - see below),
   and run one final whole-gene alignment pass for an aggregate score.
6. **Refine junctions** (0.2.0) on the winning assignment. Slot profiles
   cover each exon's full codons only, so a codon-aligned match is short
   of the true exon by the bases of any codon split across a junction
   (`exonN.frame`'s lead/trail). Those are added back, then each junction
   between adjacent slots is fine-tuned - split-codon distribution, then
   a one-codon shift either side - by whole-gene alignment score plus a
   bonus for each group II intron end (5' `GNGCG`, 3' `AY`). A stop codon
   `orfedit` folds onto any exon but the last is removed (after an
   internal exon that "codon" is intron sequence, and edit-created stops -
   CGA/CAA/CAG - are common there).

0.2.1: the intron-motif bonus is 20 (was 2). A split codon's two readings
differ by one residue and the profile often prefers the wrong one by more
than a tiebreak; calibrated on Arabidopsis (nad1/2/4/5/7, ccmFc, cox2,
rpl2, rps3 vs RefSeq) exons placed base-exact went 9/32 -> 21/32 (40 and
100 changed nothing; protein identity unchanged). Also: two slots can no
longer be filled by overlapping candidates (a short slot profile had
matched inside its neighbour's exon - Arabidopsis cox2); the lower-scoring
one's slot is left empty. Works for cis-spliced genes unchanged - the
junction classifier never assumed trans.

Before 0.2.0 split codons were dropped at every junction, frame-shifting
the joined CDS downstream of the first phase-1/2 junction: on
Arabidopsis thaliana (a training genome) nad2 came out 50% identical to
its RefSeq protein despite all five exons being placed within 2 bp.
0.2.0 gives nad1 90%, nad2 95%, nad5 95%, rps3 99% (the rest is mostly
RNA-editing differences and the tiny-exon limitation below).

**Known limitation - tiny exons.** A 7-19 codon exon (nad5 exon 3, 22 bp;
nad1 exon 4, 59 bp) is too short for an upstream HMM search to hit and
its profile too short to place it confidently de novo; its slot can end
up filled from an unrelated fragment. A neighbourhood search around
cis-adjacent exons was tried and did not fix it (a spurious hit nearby
outscored the real exon).

Critically: the algorithm never assumes a fixed cis/trans *pattern* per
junction, only fixed exon *identity/order* - real reference genomes show
the same gene's exon count is conserved across lineages while *which*
junction is spliced in cis vs. trans is not (independently gained/lost
multiple times - see the calling repo's notes). Cis/trans classification
is a per-species output label (genomic distance + strand check), computed
after assignment, never an input to it.

## Profile format

Same PSSM format as `orfedit` (plain TSV, 20 amino acid columns, one row
per alignment column) - a template directory holds `exon1.pssm`,
optionally `exon1.frame` (`lead<TAB>trail`: bases of a split codon at the
exon's 5'/3' end, 0-2 each; absent = 0/0 - profiles must then be built
from the exon's in-frame full codons only),
`exon1.threshold` (a `min_self_score` floor - not an invented cutoff, see
usage below), `exon2.pssm`, ..., and `whole_gene.pssm`.

## Usage

```
transsplice scan-batch \
  --manifest manifest.tsv \
  --genes-out reconstructed_genes.tsv \
  --exons-out reconstructed_exons.tsv \
  --junctions-out reconstructed_junctions.tsv \
  --threads 16
```

`manifest.tsv` header: `species organelle gene fasta template_dir contig
start end strand` - one row per raw exon fragment (multiple rows share the
same species+gene; all of them are candidates, not just the best hit).

Tunable: `--merge-distance` (default 5000bp), `--junction-distance-threshold`
(default 10000bp), `--exon-flank` (default 60bp), plus the same
`--gap-open`/`--gap-extend`/`--edit-penalty` scoring knobs as `orfedit`.

`min_self_score` thresholds should be computed by actually scoring a
slot's own training sequences against the profile built from them via a
real `orfedit scan-batch` call (same code path this tool uses at runtime),
not guessed - see the calling repo's `00_build_exon_profiles.py` for a
worked example.

## Output

`<genes-out>`: one row per species x gene - `n_filled`/`n_slots`,
`complete`, aggregate `whole_gene_score`/`whole_gene_n_edits`, joined
protein (`-` = unfilled/deleted position). `<exons-out>`: one row per
assigned exon (slot, corrected coordinates, strand, its own score,
`via_denovo` flag, raw-fragment provenance). `<junctions-out>`: one row
per junction between consecutive assigned exons (`cis`/`trans`, gap in bp).

## Testing

`cargo test` covers the assignment solver (including a merge-worthy
adjacent-fragment case) and the cis/trans junction classifier
(strand-switch and distance-threshold logic) in isolation. `orfedit`
itself has separate, more extensive tests for the alignment core this
tool depends on.
