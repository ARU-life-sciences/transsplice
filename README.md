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
whole-gene `profile_column`, `window_start`/`window_end`, the
stretch it could equally well sit in, clipped to the exon, and, since
0.3.3, `near_end`: true for one `--end-margin` keeps out of the gene). A frame
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
templates (ferns, lycophytes) lost most of their exons, and at 0.075,
where poa now has them, mito is 7,423 complete and 762 with a call. A
full mito run takes 17 minutes at 16 threads (13 before the searches for
short exons).

Reverse editing. `--stop-edit-penalty` and `--reverse-edit-penalty` (both
`inf`, off, by default) are orfedit's: the cost of reading a genomic
in-frame stop codon as Gln or Arg by a U-to-C edit, and of any other
U-to-C edit; see orfedit's README for why they are separate and how they
were calibrated. They apply to every exon scan and to the whole-gene
alignment, and with either on, a stop codon inside the last exon no
longer ends the gene there. For ferns, quillworts and hornworts: on 249
spliced mito genes from 30 DToL fern assemblies, `--stop-edit-penalty
2.5` takes genes aligned over 99% of their profile from 0 to 80, and
plastid ones from 74 to 117 of 233. With both off the output is
unchanged. The edits table has an `edit` column, `C>U` or `U>C`.

Editing sites. `--known-edit-penalty` and `--novel-edit-penalty` are
orfedit's too: for templates whose profiles say where their reference
proteins are edited (a 21st `edit` column; see orfedit's README), the
cost of a C-to-U edit at a column edited in them and at any other, each
`--edit-penalty` unless given. They apply to the exon scans, where a
gene's edits are called, and to the whole-gene alignment. On 107 DToL
angiosperms with RNA-seq, held out by genus, 1.5 at the known columns
and 10 away from them takes the edits called in spliced plastid genes
from 59% to 94% in the RNA, and those in spliced mitochondrial genes
from 96% to 98%, with more of them right in both (822 against 721, 7,140
against 6,963).
Templates without the column are read as before.

Stops inside a gene. An alignment cannot step over a stop codon (see
orfedit, "Where a gene ends"), which the whole-gene alignment used to do
for a stop inside an exon, and which now, with nothing else done, comes
out as two frame changes round the stop. `--internal-stop-penalty`
(orfedit's; `inf`, off, by default) aligns such a stop to a profile
column at that cost: `*` in the protein, counted in the genes table's
`n_internal_stops`. One in the last exon, with at least half the profile
aligned before it, is the gene's stop, as a skipped one was. Of 1-bp
indels planted in real spliced genes, as many are found with
`--internal-stop-penalty 20` as before any of this (mito 693 of 768,
plastid 691 of 894), with as few calls elsewhere in those genes (32 and
14; without the option, 102 and 39, the pairs round stops).
`--end-margin N` (0 by default, and poa leaves it there): a gene whose
only frame changes lie within N profile columns of its start or end is
read in one frame. At 25 it found fewer of the planted indels (mito
83%, plastid 64%) and called more elsewhere, so it is not the remedy
here that it is for single-exon genes in orfedit. Since 0.3.3 the margin
reads only the gene chosen, after the search (applied to every alignment
the search weighed, it changed which exons were chosen), and the frame
changes it keeps out are still reported, `near_end` true, as orfedit's
are. Plastid planted indels found at 25: 691 of 894, against 612 before
and 707 without a margin, with as few calls elsewhere as without (14):
still short of no margin, since a frame change can fall outside the exon
once the gene ends at its own stop.

Junction slide. `--junction-slide N` (0, off, by default) lets junction
refinement move each exon junction up to N bases along - bases off one
exon and onto its neighbour, the joined length unchanged - where the
whole-gene alignment gains by more than 1. The refinement before it tries
each split of the shared codon and a codon either side one at a time, and
cannot reach a slide, which changes all three at once. It is judged on
the protein alone: with the intron-end bonus it moved 13 junctions that
RNA reads confirm onto a chance motif nearby. Checked against
intron-spanning RNA-seq reads in 11 DToL ferns: at 6, 135 of 144
junctions give the transcript the reads show against 126 (nad5's 22 bp
exon 3 was a base or two off in all ten). It is off by default because
held out on 24 angiosperm RefSeq genomes it places 425 of 477 exons right
against 427. With it off the output is unchanged.

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
after assignment, and an input to it only between copies of an exon that
score the same, and, given what the references say of each junction,
between an exon far from its neighbour and one beside it (below).

**Exons of two copies (0.3.1).** An exon in a repeat has a twin that
scores as it does: both copies of a plastid's inverted repeat, a repeat of
a mitogenome. The beam weighs two candidates a slot, and one fragment
holding two exons fills only one slot, so a gene could be joined from
exon 1 of one copy of the inverted repeat and exon 2 of the other, tens of
kb apart on opposite strands. Now, for an exon on a junction that is not
cis, the exon is looked for next to its neighbours (as for a short slot)
and taken from there if the gene scores no lower and more of its
junctions are cis. On a circular contig (oatk's `circular=true` header)
that search, and the junction distance, go round the origin: an intron
across the origin is cis. Over one dataset of 22,160 reconstructed genes
(about 1,300 plastomes and 1,100 mitogenomes) against 0.3.0: plastid ndhB
with a trans junction in 426 species -> 8, plastid rpl2 329 -> 0, cox2
142 -> 20, nad7 63 -> 21, liverwort atp1 30 -> 5; plastid rps12's cis
junction (exon 2 to 3) 386 -> 1,238 of 1,326. No gene scores lower; 125
proteins changed, every one to a higher score, most of them (cox2, atp1)
from an exon that had been only partly aligned (median 50 unaligned
columns -> 2). Planted 1-bp indels found as before (mito 693 -> 692 of
768, plastid 691 of 894). The 8 ndhB that stay are ferns whose exon 2
copy next to exon 1 runs across the origin itself, which a search window
could not (it can since 0.3.3, below).

**Junctions the references have cis (0.3.2).** Two causes left genes
joined across loci where the references join them cis. A gene in two or
three diverged copies (plastid clpP without introns in Cyperaceae and
*Aquilegia*): the best exon per slot, from any copy, outscores every real
copy, so the pass above never takes the copy beside the neighbour. And a
short exon's chance match far away (plastid rps12's 8-codon exon 3, in
lineages that have lost its intron): its floor is near zero or below, and
it cannot be told from a real short exon by score (nad5's 7-codon exon 3
is trans on both sides). What tells them apart is the references. A
template may hold `junctions.tsv` (`exon<TAB>n_cis<TAB>n_trans`, for the
junction from exon N to N+1), and `--trans-junction-penalty P` charges a
junction that is not cis - far apart, on the other strand, or on another
contig - P times the share of references that have it cis. The charge
applies in the pass for exons apart only, after the beam and the
short-slot search, which weigh whole-gene scores alone as before (charged
there, a chimera gave way to a worse copy, and a short exon was paired
with the wrong neighbour). In that pass an exon beside its neighbour is
taken if the gene scores no lower net of the charge; failing any, an exon
far from a neighbour on its own contig is left out if the gene scores
higher without it; each round weighs every exon apart and takes the best
change; an exon is never left out for a junction to another contig alone.
A junction trans in the references (nad1, nad2, nad5, plastid rps12 exon 1)
costs nearly nothing. Without a `junctions.tsv`, or with P = 0, nothing
changes. Over the same dataset, at P = 60, against 0.3.1: junctions the
references have cis that came out trans within one contig 332 -> 162
(plastid rps12 exon 2 to 3: 87 -> 0, clpP 52 -> 8; mito rps3 24 -> 6);
junctions the references have trans 5,034 -> 5,034; no gene has more
junctions apart. 195 proteins changed: 194 to a lower whole-gene score
(median 15.6 lower, at most 61), one higher. Exon placement held out on
24 RefSeq genomes is unchanged (523 of 578 exons); planted 1-bp indels in
spliced plastid genes are found more often (691 -> 707 of 894), and in
mito ones 692 -> 697 of 768.

**Round the origin, and copies holding the whole gene (0.3.3).** A
circular contig's sequence is now continued past its end with its first
20 kb (`WRAP_BP`), so a window beside a neighbouring exon runs on round
the origin as one window, and an exon can cross the origin itself. Cut
in two at the origin, as before, a window could not hold such an exon:
fern ndhB took exon 2 from the other copy of the inverted repeat, 65 kb
away on the other strand, because its copy beside exon 1 crosses the
origin. Two fragments that reach the contig's two ends (within 30 bp) are
merged into one candidate across the origin; only those two, since
merging whole runs of fragments either side made one candidate of every
exon of a liverwort cox1 and crowded an exon's own out of its slot. An
exon found wholly in the continuation is moved back to the contig's start;
one across the origin keeps its end past the contig's length, as poa
writes such genes (GFF3's convention for circular sequences). On a tie,
a gene within its contig is preferred over one across the origin: the
two copies of a plastid's inverted repeat make the same protein, and the
one across the origin is the harder to use.

The pass for exons apart also moves a run of exons together. Where no
exon beside a neighbour will do - one copy of a repeat lacks it - the
other copy may hold the whole gene, and moving one exon at a time only
moves the break (Dryopteris nad2 and nad7: exons 1-3 from one copy of a
mitochondrial repeat and the last from the other, which holds all of
them). For each junction apart, the exons on one side are found again
one by one beside the other side's, outward from the junction, and the
run is taken by the same rule as one exon. Of 139 genes in the dataset
below with a junction apart that the references have cis, 98 have no
whole cis copy in their assembly at all (tblastn chains of >= 70%
identity covering 90% of the protein), so the assembly holds them in
pieces; given only the cis copy's fragments, 15 of the other 41 made an
all-cis gene scoring the same, higher, or within the junction cost.

Over the same dataset against 0.3.2: junctions the references have cis
that came out trans within one contig 162 -> 133 (fern plastid ndhB 6 ->
0, mito nad7 25 -> 16, nad2 25 -> 22, nad1 14 -> 11, cox2 12 -> 9,
plastid clpP 8 -> 6); junctions the references have trans 5,034 ->
5,031 (three Geranium nad1, intronless from exon 2 to 5, now read from
that one locus); complete genes 21,616 as before, and none with fewer
exons. 59 proteins changed, 54 to a higher whole-gene score (median 149
higher: exons cut at the origin, now whole - mito ccmFc, nad1, nad2 and
nad4 of mosses and liverworts, Sphagnum inundatum nad2 exon 2 925 ->
1629), five lower, each now from one copy: Carex rariflora and C.
vaginata clpP from their intronless copy (55 and 11 lower), Vandenboschia
cox2 (39 lower), and two horsetail rps12, exon 3 beside exon 2, its
intron lost. Planted 1-bp indels found
as before (plastid 707 of 894, mito 697 of 768); exon placement held out
on 24 RefSeq genomes as before (525 of 578 with the options poa passes).

## Profile format

Same PSSM format as `orfedit` (plain TSV, 20 amino acid columns, one row
per alignment column) - a template directory holds `exon1.pssm`,
optionally `exon1.frame` (`lead<TAB>trail`: bases of a split codon at the
exon's 5'/3' end, 0-2 each; absent = 0/0 - profiles must then be built
from the exon's in-frame full codons only),
`exon1.threshold` (a `min_self_score` floor - not an invented cutoff, see
usage below), `exon2.pssm`, ..., and `whole_gene.pssm`; optionally
`junctions.tsv` (`exon<TAB>n_cis<TAB>n_trans`, a header line allowed: in
how many references the junction from exon N to N+1 is cis and in how
many trans; read by `--trans-junction-penalty`).

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
(default 10000bp), `--exon-flank` (default 60bp), `--trans-junction-penalty`
(default 0, none; poa uses 60), plus the same
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
