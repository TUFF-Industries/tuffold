# Cross-checks

Dev scripts that check Tuffold's validation against established tools and
yardsticks it isn't trained on. All read projects without writing them;
each script's docstring has its usage.

| script | checks | needs |
|---|---|---|
| `synthetic.py` | mokapot on the seeded synthetic data of `tuffold-validate`'s tests | mokapot |
| `project.py` | mokapot on a real project's rank-1 PSMs with Tuffold's features; with R, its entrapment estimate | mokapot 0.10, pandas |
| `entrapment_fasta.py` | builds a database with a shuffled twin of every (or every chosen) protein, prints r | stdlib |
| `entrapment.py` | entrapment estimates for Tuffold projects (its q-values, or TDC) and TPP pepXML | stdlib |
| `embed_pepxml.py` | embedded entrapment: renames twins to their real protein in a pepXML, so pooled evidence (siblings, NSP) can lift them | stdlib |
| `embed_tuff.py` | embedded entrapment inside a Tuffold project (no new search): re-points twins to their real protein (`same`) or to another real protein (`cross`) | stdlib |
| `mokapot_siblings.py` | Tuffold's peptide-level pooling inputs in mokapot, with decoys' sibling counts none, unpaired, paired, or imputed; scored by `entrapment_seq.py`'s `tsv` mode | mokapot 0.10, pandas |
| `entrapment_seq.py` | entrapment estimates by sequence (twins sharing names), for TPP pepXML (decoy cut or iProphet's own model) and Tuffold; with a decoy cut also twins per decoy, split into near-copy/far and shadow/free twins | stdlib |
| `proteometools.py` | ground-truth FDP on ProteomeTools synthetic peptide pools, for TPP pepXML and Tuffold | stdlib |
| `headroom.py` | correct matches a validation leaves below its 1% cut (targets minus decoys), i.e. whether rescoring is saturated | stdlib |
| `rt_consistency.py` | whether the extra spectra of one validation over another sit at their peptide's retention time | stdlib |
| `sweep.py` | validates copies of projects several ways (deep features off and on by default) and tabulates them | a `tuffold` binary |

## An entrapment check end to end

```
python3 -I tools/crosscheck/entrapment_fasta.py yeast_ups.fasta entrap.fasta --only OS=Saccharomyces
tuffold search runs/*.mzML --fasta entrap.fasta -o entrap.tuff
python3 -I tools/crosscheck/sweep.py --tuffold target/release/tuffold --out scratch entrap.tuff
python3 -I tools/crosscheck/entrapment.py tuff 1.0452 scratch/entrap-*.tuff --q 0.005 0.01 0.02
python3 -I tools/crosscheck/rt_consistency.py scratch/entrap-plain.tuff scratch/entrap-deep.tuff
python tools/crosscheck/project.py scratch/entrap-plain.tuff 1.0452
```

(1.0452 is the r `entrapment_fasta.py` printed.) An external engine's
results work the same way: search the entrapment FASTA with decoys added,
import, and sweep.

A peptide counts as entrapment when every non-decoy protein it maps to is a
twin, so a twin peptide that is also some decoy's sequence counts. Before
these scripts, one scratch version counted such peptides as real, which put
the PXD006447 Comet baseline at 1.02% instead of 1.11% (47 spectra); the
earliest deep-feature numbers used that version.

Estimates are only as good as their assumptions: entrapment can't see a
wrong spectrum on a present peptide (that's what `rt_consistency.py` is
for), and twins that differ from their real peptide in a few positions share
most fragments, which near-copies in any real database do too.

## Embedded entrapment

Standard twins are proteins of their own, absent from the sample, so
nothing that pools evidence over a protein ever lifts them. A wrong peptide
of a *present* protein is what sibling evidence lifts, and decoys (whose
proteins are never present) aren't lifted unless paired. To measure that,
give each twin its real protein's name: search the twin database as usual,
then rewrite the pepXML (`embed_pepxml.py`) before iProphet or a Tuffold
import, and score by sequence (`entrapment_seq.py`). The PSMs are the same
as the standard design's, so both can be compared on one search. Report the
lower bound beside the combined estimate here: a wrong spectrum on a present
real peptide is no wrong peptide, but on a twin it always is, so with
pooling the combined estimate reads high (a little: weighting each twin by
its protein's coverage moved it from 1.39% to 1.28% on PXD006447).

Every real protein has a twin, so an embedded twin lands in a protein with
identified peptides as often as a wrong real match does (PXD006447 Comet:
11.5% of low-scoring twins have a sibling, 10.3% of real wrong targets). So
this is the natural mix, not a worst case. Standard twins (0.9% there)
understate what pooling adds.

Then read the twins-per-decoy columns. If decoys count what twins stand for,
twins pass at r/(1 + r) of the decoys (1.00x). Shuffled twins keep their
real peptide's composition and mass, which creates two kinds no wrong real
peptide has:
- near-copies, which differ from the real peptide in a few positions;
- shadow twins, whose real peptide passes, so they can take its spectra.

Sibling evidence lifts both kinds, and reversed decoys can't stand for
them. Wrong real peptides are rarely like that inside their own protein:
same-composition near-copies of passing yeast peptides sit in another
protein 755 times in 756. So far and free twins are the ones decoys must
count; the other columns show what the twin design adds.

To take the shadow twins out by construction, embed each twin in *another*
real protein: `embed_tuff.py cross` renames the twin of Yᵢ (and its decoys)
to the next real protein by length, inside a validated standard-design
project, so a twin still lands in a present protein at the natural rate but
the protein whose spectra it takes is not the one whose siblings it gets.
`embed_tuff.py same` is the control and reproduces `embed_pepxml.py`'s
results exactly. On PXD006447 (Comet, four twin seeds) twins per decoy go
from 1.41–1.72× (same) to 0.91–1.10× (cross) with paired siblings and from
1.67–2.06× to 0.94–1.17× imputed.

## Ground truth: ProteomeTools pools

ProteomeTools (PRIDE PXD004732) measured pools of ~1,000 synthetic human
peptides, so what is in a run is known. Searched against a whole proteome,
every other identification is wrong, including wrong peptides of the
pool peptides' own proteins, the case pooled evidence lifts. Per pool,
download the run's `.raw` and its MaxQuant `-tryptic.zip` (its
`peptides.txt` is the pool's peptides as found by a search of the pool
alone). Then:

```
python3 -I tools/crosscheck/proteometools.py truth sprot.fasta truth.tsv 61=pool61/peptides.txt 62=pool62/peptides.txt ...
# search the runs (decoys added, reverse_ prefix), run TPP and import into Tuffold as for entrapment
python3 -I tools/crosscheck/proteometools.py score truth.tsv sprot.fasta tuff project.tuff Tuffold
python3 -I tools/crosscheck/proteometools.py score truth.tsv sprot.fasta tpp ipro.pep.xml interprophet iProphet
python3 -I tools/crosscheck/proteometools.py bins truth.tsv interact.pep.xml
```

The first-pool DDA runs fragment every precursor twice (Orbitrap HCD and ion
trap CID); `activation_method = HCD` in Comet keeps one. Comet needs indexed
mzML (ThermoRawFileParser `-f=2`). The truth list misses some real analytes
(carryover from earlier pools, MaxQuant's missed ~2%, synthesis by-products
such as `…K` read as `…KK`), so every method shows a floor above its nominal
FDR; `bins` shows where it comes from (decoys and unlisted targets agree
where errors dominate, and only unlisted targets sit at the top). Compare
methods on one search rather than reading the floor as their error.
