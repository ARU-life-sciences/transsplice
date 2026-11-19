# transsplice

Reconstructs trans-spliced genes in plant mitochondrial genomes from
scattered per-exon gene-model hits, reusing
[`orfedit`](https://github.com/ARU-life-sciences/orfedit)'s edit-tolerant alignment
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
3. **De novo search**: if a slot has no candidate scoring above its own
   floor, scan the contig(s) the other exons were found on (bounded, not
   genome-wide) for that slot's profile directly. A short slot (profile
   of 30 columns or fewer) is always searched for this way, filled or
   not - see step 7.
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

0.3.0: frameshifts. With orfedit 0.2.0's frame-change alignment, each
exon's scan and the whole-gene pass over the joined exons may change
reading frame (`--frameshift-penalty`, default 20; `inf` for none), so
an exon with an indel is still placed whole and the joined CDS aligned
through it. Each frame change in the whole-gene alignment is mapped back
to its exon and genomic position and written to `--frameshifts-out`
(`slot`, `contig`, `start`/`end` of the 1 or 2 skipped bases, `n_bases`,
whole-gene `profile_column`, and `window_start`/`window_end`, the
stretch it could equally well sit in, clipped to the exon). A frame
change within 3 bp of an exon end is taken as the junction being a base
or two off, not reported. The joined exons are now aligned in frame 0
only (0.2.1 let that alignment start in any of the three frames), which
changes results even with `inf`: on the 8,993 DToL mito genes, 528
differ from 0.2.1, all nad5 bar four; 174 nad5 become complete (none
stop being complete), and the median number of profile columns left
unaligned in the changed proteins falls from 50 to 1. Works for any
organelle's templates: `poa` builds plastid ones (atpF, clpP, ndhA,
ndhB, rpl2, rpoC1, rps12, rps16, ycf3) from 13 RefSeq plastomes.

Speed and reproducibility (0.3.0). A batch runs in three stages, each
a pool of parallel tasks: score and assign every gene's candidates, run
every whole-contig search for an empty slot (largest first), then choose
among assignments and refine junctions. Each (candidate, slot) scan and
each distinct joined sequence is scored once per gene. With orfedit's
faster alignment, the 8,993 DToL mito genes take 13 minutes and 0.7 GB
at 16 threads (2 h 54 min and 9.8 GB before), with the same genes table
and, exon by exon, the same scores and edits; 22 of 28,191 exons sit on
a different, equal-scoring copy of a repeat. Candidates and output rows
are now in a fixed order, so two runs give byte-identical tables (before,
equal-scoring repeat copies won at random).

7. **Short exons** (0.3.0). A short exon (nad5 exon 3, 22 bp; nad1 exon
   4, 59 bp; plastid rps12 exon 3, 26 bp; rps16 exon 1, 40 bp) is rarely
   hit by an upstream HMM search, and a chance match elsewhere can clear
   its floor - so its slot used to be filled from an unrelated fragment,
   or left empty when that fragment overlapped a neighbouring exon. Now
   every slot with a profile of 30 columns or fewer gets, whatever fills
   it, its three best whole-contig matches as extra candidates, which the
   assignment search weighs with the other slots (the candidates at a
   short slot's best three places; two candidates for any other slot) - a wrong short exon, or none, can put
   its neighbours out of frame, so they have to be chosen together. Then
   those matches and the best match within 5 kb after the previous exon
   and before the next one are each tried in the short slot's place: the
   assignment is re-evaluated (junctions refined), and the alternative
   kept only if the whole-gene alignment score plus the intron-end
   evidence (5' `GNGCG` +7 / absent -1.5, 3' `AY` +2.5 / absent -2.5,
   as log2 odds) is higher. Exons found this way have `via_denovo` true
   and provenance `denovo` or `nearby`.
8. **Gene ends** (0.3.0). Profiles are aligned end to end, so where a
   species' protein starts later or stops sooner than the references'
   the first exon began upstream of the start codon and the last ran on
   past the stop (to the aligner a stop is just a skipped codon). After
   refinement the first exon is cut back to the start codon - the first
   aligned codon if it is Met, else the first ATG within 15 profile
   columns that sits on a column where the references have Met - and the
   last exon ends at the first in-frame stop codon: one inside it, if at
   least half the profile is aligned before it, else one within 15
   codons of its last aligned codon (a `CAA`/`CAG`/`CGA` right at the end
   still counts as an edit-made stop, unless a literal stop follows it
   directly).

Overlapping exons are now dropped best-first (before, an exon that was
itself dropped for overlapping a better one could still knock out a
third).

Measured with `poa`'s templates (`build/exon_impl` there). On the 24
RefSeq genomes the templates come from, each with its own genome left
out of the templates: 89% of 477 exons within 3 bp of GenBank at both
ends (74% with the same templates before steps 7-8; 68% base-exact
against 53%), exons of 30 bp or less 15 of 16 (8), last exons ending
exactly at the stop codon 142 of 169 (92), first exons starting exactly
at the start codon 166 of 178 (151), and no exon that was right made
wrong. On the DToL assemblies: nad5 exon 3 is within 20% edit distance
of a reference exon in 861 of 927 mito assemblies (none of 775 placed
before), nad1 exon 4 in 928 of 928 (79 of 834), plastid rps12 exon 3 in
970 of 1,139 (0 of 837) and rps16 exon 1 in 816 of 818 (2 of 777).
Complete mito genes go from 5,957 to 6,925 of 8,566 and plastid ones
from 9,373 to 9,792 of 9,954; genes with a frameshift call fall from
1,521 to 1,004 (mito) and 362 to 98 (plastid; rps12 136 to 0). These
are with `poa`'s slot floors at a quarter of the lowest reference
self-score; at half, as they were, divergent lineages on the angiosperm
templates (ferns, lycophytes) lost most of their exons. A full mito run
takes 17 minutes at 16 threads (13 before the searches for short exons).

Still wrong: a gene whose exon count differs from the template's (cox2
with one intron or two; an intron lost so two template exons are one
stretch of genome, which one candidate can't fill twice - nad7 in
Rosaceae, nad4 in liverworts); an intron-side exon end more than a codon
from where the slot profile ends; and a last exon whose stop is made by
editing somewhere the profile doesn't end (ccmFc).

**The whole-gene profile decides step 7, so it has to be right for the
short exon's own columns.** A whole-gene profile refined from earlier
reconstructions in which the short exon was misplaced (poa's
`02_bootstrap_refine.py`) votes for the misplaced sequence.

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
