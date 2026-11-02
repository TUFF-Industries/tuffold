# Tuffold user guide

Tuffold validates MS/MS proteomics results and lets you browse and quantify
them. You bring search results (or spectra and a FASTA, and Tuffold
searches them); Tuffold estimates how many of the identifications are wrong
with target-decoy competition, groups proteins, and shows everything linked
together, down to annotated spectra. It works offline: nothing leaves your
computer.

To install it, see [install.md](install.md). Every step below also has a
command (`tuffold …`), listed at the end of each section.

## Contents

1. [Projects](#1-projects)
2. [Bringing data in](#2-bringing-data-in)
3. [Validation and thresholds](#3-validation-and-thresholds)
4. [Proteins](#4-proteins)
5. [Samples](#5-samples)
6. [Quant: comparing sample groups](#6-quant-comparing-sample-groups)
7. [Exports](#7-exports)
8. [GO terms](#8-go-terms)
9. [Undo](#9-undo)
10. [Coming from Scaffold](#10-coming-from-scaffold)
11. [AI assistants](#11-ai-assistants)
12. [When something is left out](#12-when-something-is-left-out)

## 1. Projects

Everything lives in one project file, `name.tuff`: the identifications,
scores, spectra, validation results, sample groups, and your marks. It is a
SQLite database that Tuffold versions; newer Tuffolds open older projects,
and an older Tuffold refuses a newer one rather than misread it. Opening a
project never changes it.

- **New project from data…** (⌘N) picks files and creates a project from them.
- **Open project…** (⌘O) opens a `.tuff`.
- **Open the demo project** builds one from synthetic data, marked
  *SYNTHETIC DATA* in the title bar, to try every view without lab data.

With a project open, one bar at the top holds its name, the views (Overview,
Proteins, Samples, Quant; ⌘1–4), and its actions: + adds data, the magnifier
searches spectra, *Export* opens the export menu, and the folder opens another
project. Hover a button for its name and shortcut.

`tuffold demo -o demo.tuff`, `tuffold info study.tuff`

## 2. Bringing data in

**Add data** (+ in the bar, ⌘I) adds files to the open project. Tuffold reads:

| Source | What to pick |
|---|---|
| FragPipe | the results folder (MSFragger's own pepXML files are used when present) |
| MaxQuant | the `txt` folder (`msms.txt`, or `evidence.txt`) |
| Mascot | `.dat` |
| Proteome Discoverer 1.4–3.x | `.pdResult` or `.msf` |
| Scaffold 4.0–5.3 | `.sf3` |
| Scaffold 5, Scaffold DDA 6.5 | `.sfdb` |
| Any engine | pepXML (`.pep.xml`, `.pepXML`) or mzIdentML (`.mzid`), gzipped or not |
| Spectra | mzML, MGF, or a Bruker timsTOF `.d` folder (DDA-PASEF) |

**Thermo .raw files** aren't read directly: convert them to mzML first,
then add or search the `.mzML`. Either converter works:

- [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser):
  `ThermoRawFileParser -i run01.raw -o out/ -f 2` (`-f 2` is indexed mzML).
- ProteoWizard `msconvert run01.raw --mzML` (Windows), with peak picking
  on if the MS/MS spectra are profile (Tuffold centroids profile spectra
  itself, but the files are much larger).

Adding a `.raw` file says so instead of importing it.

Each search engine in a file becomes its own import. Spectra files attach
their peaks to the identified spectra (matched by run and scan number) and,
when they hold MS1 scans, let Tuffold measure precursor intensities. A
timsTOF `.d` gives one MS/MS spectrum per precursor (its PASEF frames
summed), and its MS1 frames give intensities, each precursor measured only
near its own ion mobility (±0.015 1/K0), so other ions at the same m/z
don't add to it. Match-between-runs looks for a transfer at the mobility
the other runs saw it at. An MS1 scan that is damaged (a timsTOF frame
whose data won't decode, or an mzML scan holding peaks no instrument could
have measured, such as an m/z below a proton's or an intensity above 10^15)
is left out of intensities and purity, and the import says how many. The
MS/MS spectra that name it as their precursor scan are still read.

After an import, the Overview lists what came in and, under **Findings**,
everything that didn't, with the reason: a file it couldn't read, PSMs
whose mass doesn't fit their spectrum, unknown modifications. Nothing is
dropped silently.

**Give Tuffold unfiltered results.** Validation needs the decoy matches.
Many programs write only what passed their own FDR filter, and decoys are
then missing or rare: set the PSM FDR to 100% in MaxQuant, export
Proteome Discoverer results before its validator, and keep FragPipe's
MSFragger pepXML files. Tuffold says when a file looks filtered.

### Searching spectra yourself

**Search…** (⌘⇧F) searches mzML or MGF files, or Bruker timsTOF `.d` folders,
against a FASTA with the
search core of [Sage](https://github.com/lazear/sage): pick the FASTA, then
the spectra. The app uses common settings (trypsin, 2 missed cleavages,
20 ppm precursor and fragment, carbamidomethyl C fixed, oxidized M
variable). Decoys are generated (reversed proteins); any decoys in the
FASTA are dropped. Before the search starts, Tuffold estimates the memory
its index needs and refuses a search that wouldn't fit, rather than let the
computer swap.

For other settings use the command, e.g. a semi-tryptic phospho search:

```sh
tuffold search run1.mzML run2.mzML --fasta human.fasta -o study.tuff \
    --semi --variable S:Phospho --variable T:Phospho --fragment-ppm 10
```

Enzymes: trypsin, trypsin/p, lys-c, lys-n, arg-c, glu-c, asp-n,
chymotrypsin. Modifications are `SITE:MOD`, by Unimod name or mass
(`K:TMT6plex`, `C:57.021464`; `^`/`$` are peptide termini, `[`/`]` protein
termini); `--label tmt` or `tmtpro` fixes an isobaric label.
`--max-memory 12` allows the index 12 GB. For low-resolution (ion-trap)
MS/MS give the fragment tolerance in Da, `--fragment-da 0.5`, instead of
`--fragment-ppm` (one or the other).

`tuffold import … -o study.tuff`, `tuffold search …`

## 3. Validation and thresholds

Click **Validate now** on the Overview (or run `tuffold validate`); after
new data is added the Overview says the results are out of date. For each engine Tuffold:

1. keeps each spectrum's top-ranked match;
2. **rescores** it (Percolator style): a linear model trained on the
   engine's own scores, cross-validated in three folds so no match is
   scored by a model that saw it. When every match has a retention time,
   one feature is how far it eluted from when its sequence says it should
   (a per-run model from amino-acid composition, trained on confident
   matches the same cross-validated way). With the deep-learning models
   installed (below), each match is also compared with the spectrum
   AlphaPeptDeep predicts for its peptide, and a deep retention-time model,
   fine-tuned on each run's confident matches, replaces the composition
   model. A second pass adds how many other spectra confidently identify
   the same peptide (a decoy counts those of the target it was made from,
   so decoys keep measuring the errors this adds); when the project's decoys
   can't be matched to their targets, it's left out and the Overview says so. The model is kept only if it
   beats the engine's best single score; the Overview says what each
   engine was ranked by;
3. lets engines **vote** per spectrum when several searched the same data
   (the strongest calibrated vote wins; votes are never added up);
4. scores each **peptide** from all its spectra (iProphet style): its best
   and second-best spectrum, how many other spectra and charges match it,
   and how many other peptides of its protein are identified, by the same
   cross-validated rescoring, kept only if it beats the best spectrum
   alone. A decoy peptide counts the identified peptides of the target
   protein it was made from, so decoys still count the wrong peptides this
   lifts; when decoys can't be matched to their targets (a database whose
   decoys are numbered rather than named after their targets), protein
   evidence is left out and the Overview says so. Modified forms are not
   pooled, and spectra never are: a wrong spectrum on a peptide that is
   really there would be boosted by the right ones, and decoys cannot count
   such errors;
5. estimates **q-values** by target-decoy competition at the spectrum,
   peptide, and protein-group level, and a **probability** (1 − PEP) for each.

Proteins are grouped by parsimony: proteins with the same peptides form one
group, peptides shared between groups go to the group with the most
peptides (razor peptides; ties go to the better-scoring group), and groups are thresholded with picked
target-decoy competition.

**Thresholds** change without validating again. By **FDR** (q-value; 0.1%,
1%, or 5% per level) or by **Probability** (e.g. 95% peptides, 99%
protein groups, as in Scaffold), plus a minimum number of peptides per
group. In probability mode the decoy-estimated FDR of what passes is shown
too: a 99% cut says each item is likely right, not how many in the list
are wrong.

Other values: `tuffold thresholds study.tuff --protein-fdr 2% --min-peptides 2`
or `--by probability --peptide-prob 90%`.

The Overview's engine table shows how many spectra each engine finds alone
at 1% FDR, and notes say when something limits the result: an engine
without decoys (not validated), too few decoys to train a model, decoys so
rare the file was probably filtered before Tuffold saw it.

### PTM site localization

When a passing peptide carries a variable modification that could sit on
more than one residue (phosphorylation on S/T/Y, oxidation on M,
acetylation on K, and so on), Tuffold scores every candidate placement from
the spectrum's site-determining b and y ions (Beausoleil et al.'s Ascore,
weighted across peak depths 1–10), turns the score differences into
placement probabilities, and sums them into a probability per site. Ions
are matched at 20 ppm, or at 0.5 Da when the spectra are low-resolution
(Tuffold decides from the data, and a validation note, also shown by
`tuffold localize`, says which it used and the matches that decided it). To check how often a confident site call
is still wrong, Tuffold also lets the modification onto Alanine as a decoy
site and reports the decoy-estimated false localization rate (FLR) on the
Overview and in `tuffold localize`. The estimate is checked at the default
95% site probability and above; at lower cuts it can run low (two or three
times, on synthetic spectra at 90%) where candidate sites sit next to each
other, since one stray ion then decides between them. Spectra without attached peaks, with
more than 256 candidate isoforms, or whose isoforms can't be fragmented are
counted and reported, never dropped silently. A peptide's sites are those
of its best-localized passing spectrum. Attaching spectra files later needs
validating again before their spectra are localized.

```
tuffold localize study.tuff --site-prob 95% --export localization.tsv
```

### Deep-learning models

Tuffold can rescore with AlphaPeptDeep's pretrained models (Mann lab,
Apache-2.0). They aren't bundled: download `pretrained_models_v3.zip`
from AlphaPeptDeep's GitHub releases (about 29 MB) and run

```
tuffold models --install pretrained_models_v3.zip
```

`tuffold models` says whether they're installed. From then on validation
uses them whenever it rescores a project whose spectra are attached (it
takes a few minutes for a large project); `tuffold validate --no-deep`
skips them. The Overview's notes say which collision energy, instrument,
and fragment tolerance (20 ppm, or 0.5 Da for ion-trap spectra) fit your
spectra, or that the model doesn't describe them, in which case only the
retention-time model is used. Without the models installed, the notes say
that validation ran without them.

## 4. Proteins

The **Proteins** view (⌘2) lists the protein groups passing the thresholds.
Select one to see:

- its header: accession, description, peptides, spectra, coverage, its GO
  slim terms, and **☆ Star** to mark it (stars are kept in the project);
- **Similarity**: other passing groups sharing its peptides, as a peptide ×
  group matrix, and the size of its protein cluster;
- its **peptides** (unique, razor, or shared), with localized PTM site
  probabilities in the **Sites** column, each linked to its **spectra**;
- the selected spectrum, **annotated** with b and y ions at 20 ppm or
  0.5 Da, the fragment ladder over the sequence (with localized PTM sites
  badged and their site probabilities and Ascore shown under the spectrum
  header), and every engine's scores.

Tables sort by any column. A peptide with very many spectra lists the best
500 and says how many there are.

**Filters** above the list narrow what is shown (never what is validated,
quantified, or exported): text in the accession or description, a minimum
number of spectra, organism (UniProt's `OS=`), a modification, a GO slim
term, and starred groups only. Type in the organism and GO menus to search
them. When the group on show is filtered out, the first listed group is
shown instead.

## 5. Samples

The **Samples** view (⌘3) shows the same groups with one column per
biosample. **Show** picks the measure:

- **Spectra**: spectral counts;
- **Exclusive spectra** and **exclusive peptides**: from peptides unique to
  the group;
- **NSAF**: counts normalized by protein length and the sample's total;
- **emPAI**: 10^(observed / observable peptides) − 1, with observable
  peptides fully tryptic, 6–30 residues, no missed cleavage.

Double-click a group to open it in Proteins. The filters apply here too.
Top 3 intensities (MS1) are in the exports.

## 6. Quant: comparing sample groups

On the **Quant** view (⌘4), first assign biosamples to **categories** under
*Sample groups* (e.g. control, treated) and **Save groups**; a biosample
without a category is left out. Then:

- **Compare**: *Two categories* (a baseline and the one compared to it) or
  *All categories* (one-way ANOVA, with a clustered heat map of what differs),
  or *Sample QC* (below).
- **Measure**: *Spectral counts*, *MS1 intensity* when spectra files with
  MS1 scans are attached, or *Reporter ions* for isobaric-labelled pools.
  Intensities trace each identified precursor
  through its elution peak; proteins are summarized by median polish, which compares
  samples through the peptides they share. Missing values stay missing (no
  imputation), except where match-between-runs measured them (below).
- **Test** (two categories, counts): **Welch's t-test** across replicates
  (the default when both categories have two or more samples), **Fisher's
  exact test** on pooled counts (for samples without replicates; it
  overstates significance when replicates exist), or **Mann–Whitney U**
  (rank-based; with few replicates it cannot reach significance, and the
  view says so). Intensities and reporter ions use Welch's t-test.

Results show as a volcano plot (log2 fold against −log10 q) and a table,
with Benjamini–Hochberg q-values. **Export TSV…** writes the table.

`tuffold samples study.tuff --set ctrl-1=control --set trt-1=treated`,
`tuffold quant study.tuff --compare control treated --export quant.tsv`,
`--by intensity`, `--by reporters`, `--anova`

### Match between runs

A peptide identified in some runs is often missing from others only
because no MS/MS spectrum happened to pick it there, though its MS1 peak is
present. After validating, press **Match between runs…** on the Quant tab
(with *Measure: MS1 intensity* chosen) and pick the runs' spectra files again,
with their MS1 scans (mzML files, or timsTOF `.d` folders; each named like its
run), or run `tuffold mbr study.tuff run1.mzML run2.mzML …`. It looks for each
passing precursor of the other runs in every given run: where it should
elute (each other run's retention times are mapped onto this one's through
the precursors both identified) and at its m/z. Every peak found is checked
against a decoy, the same m/z looked for at the wrong time, and transfers
are kept at 1% FDR by target-decoy competition, scored by a model trained on
the run's own identified precursors (retention-time, mass and isotope
errors, intensity, peak width, and how far the intensity is from what the
other runs measured). A wrong identification in another run whose ion is
also in this one is transferred like a right one: no decoy can tell them
apart, so the identification FDR is what limits those. Accepted transfers fill MS1 intensity values
only where the run identified nothing of that precursor and the peak is as
strong (intensity, scans) as all but the weakest 5% of the run's own
identified precursors, and the Quant view says how many values came from
them (the line under the toolbar: "N values matched between runs", which
stays when the project is reopened). When matching finishes the tab shows,
for each run, how many of the precursors looked for were transferred, and says
if damaged MS1 scans were left out or too few known peaks trained no model
(then retention-time error alone ranks the peaks). A file that is no run of
the project, or a project that isn't validated, is reported as an error. The
button is greyed out until the project has MS1 data to put the transfers
beside. Transfers are quantified but never used to normalize samples.
Matching cannot be undone from the app (it takes no undo point); importing new
data clears the transfers, and you match again after validating.

### Sample QC

**Compare: Sample QC** asks whether the samples look like what they're labelled
as, on whichever **Measure** is chosen, before you compare anything:

- **PCA**, over the groups measured in every sample (nothing is imputed).
  Replicates should sit together; a sample among the wrong group may be
  swapped, and groups split by run date point at a batch effect.
- **Correlation**: Pearson's r of every pair of samples, over the groups both
  measured. On spectral counts abundance raises every r, so look at which are
  lowest rather than how high they are.
- **Samples**: groups measured, median, and mean r with the others, where a
  failed run or a poorly loaded channel shows.

Samples without a category are included here. **Export TSV…** writes it all,
or `tuffold quant study.tuff --qc --export qc.tsv` (with `--by intensity` or
`--by reporters`).

### TMT and iTRAQ (isobaric labels)

An isobaric run holds every sample at once, one per reporter channel. Search
with the label (`tuffold search … --label tmt`, `tmtpro`, `itraq4`, or `itraq8`,
or import an engine's results and the run's spectra) and validate. On the
**Quant** view, choose **Measure: Reporter ions**, pick the reagent set under
*Plex · channel samples*, name the biosample in each channel, and click **Save plex**,
then give those channel biosamples categories under *Sample groups* and click
**Save groups**. Or from the command line:

`tuffold plex study.tuff pool-1 --label tmt10 126=ctrl-1 127N=ctrl-2 …`

Reagent sets: `tmt6`, `tmt6-original`, `tmt10`, `tmt11`, `tmtpro16`,
`tmtpro18`, `itraq4` (114–117), `itraq8` (113–119, 121). iTRAQ reporters are
summed within 60 ppm; TMT's take the top peak within 20. TMT6 kits sold until about 2012 carry their 127 and 129
reporters 6 mDa heavier than today's (`tmt6-original`); when another reagent
set completes more spectra, the Quant view and `tuffold quant` say so.
Without
`--run`, the plex takes every run not already in another plex (the
fractions of one pool). Channel samples are ordinary biosamples, so give
them categories with `tuffold samples` and compare them (two categories or
`--anova` across all categories):

`tuffold quant study.tuff --by reporters --compare control treated --export tmt.tsv`

Each passing spectrum's reporter ions are read from its peaks (the
strongest peak within 20 ppm of each channel), or, for SPS-MS3 methods
(Orbitrap Fusion, Lumos, Eclipse, Ascend), from the MS3 scan acquired from
it: attach the run's mzML (ThermoRawFileParser writes which MS/MS scan each
MS3 came from) and Tuffold stores each MS3 with its spectrum, and the Quant
view says how many reporters came from MS3 scans. MGF files don't carry that
link. When the run's mzML has its MS1 scans, Tuffold also measures each
spectrum's precursor purity, the share of its isolation window that was
the identified peptide, and leaves out spectra under 50% (co-isolated
peptides add their own reporters and pull ratios toward 1:1);
`--min-purity 0.75` is stricter, `--min-purity none` keeps every spectrum.
The cut applies to SPS-MS3 reporters too (it measures the MS/MS window the
MS3 fragments came from), as Proteome Discoverer's co-isolation threshold
does.
A spectrum missing any
channel is left out and counted. Proteins sum their unique and razor
spectra per channel, channels are normalized on the median protein (equal
loading), and categories are compared with Welch's t-test or one-way ANOVA
(with the same clustered heat map). Ratios are
relative within the plex.

Each lot comes with a certificate of analysis listing, per channel, the
percentage of its reporter ions that carry one or two heavy atoms more or
less than they should, and so turn up in a neighbouring channel. Copy it
into a text file, one line per channel of the kit (every channel, zeros for
a pure one): the channel's name, then its four percentages for −2, −1, +1,
+2 Da (TMT6/10/11 and iTRAQ certificates), or eight for TMTpro's −2C, −N−C, −C, −N,
+N, +C, +N+C, +2C columns; `NA` means none, `#` starts a comment.

```
# TMT10 lot …
126   0.0  0.0   5.09  0.0
127N  0.0  0.25  5.27  0.0
…
```

Load it with *Certificate…* in the plex panel, or
`tuffold plex study.tuff pool-1 --label tmt10 --impurities lot.txt 126=…`
(`--impurities none` removes it). Each spectrum's reporters are then unmixed
before anything else, and the Quant view says so. A plex saved again
without a certificate keeps the one it has.

**Several plexes.** An experiment with more samples than channels runs as
several plexes, one per labelled pool. Put the same pooled reference sample
(a mix of all samples) in a channel of each, define each plex with its runs,
and mark the reference channel. In the app, click **New plex**, name the
plex, pick its runs under *Runs* (runs already in another plex are disabled;
leaving every run unselected takes all runs for the first plex only), mark
the reference channel with **ref**, and click **Save plex**; use the *Plex*
buttons to switch between plexes, or **Remove plex** (twice to confirm) to
drop one. Or from the command line:

`tuffold plex study.tuff pool-1 --label tmt10 --run A1 --run A2 --reference 126 126=ref-1 127N=ctrl-1 …`
`tuffold plex study.tuff pool-2 --label tmt10 --run B1 --run B2 --reference 126 126=ref-2 127N=ctrl-4 …`

`tuffold quant --by reporters` (and the Quant view) then puts the plexes on
one scale: for each protein, every plex is shifted so its reference reads the
plexes' average (internal reference scaling), and the reference channels are
left out of the comparison. A plex whose reference has no value for a protein
can't be bridged, so its values for that protein are left out and counted.
`--plex pool-1` still quantifies one plex alone, and `tuffold plex study.tuff
pool-1 --remove` removes a plex.

Interference MS1 can't see still compresses ratios somewhat; SPS-MS3
methods remove most of it.

## 7. Exports

| Export | What | Menu / command |
|---|---|---|
| Excel workbook | Summary, Methods, Proteins, Peptides, Spectra, and ANOVA sheets at the current thresholds (targets only) | *Export* → *Excel workbook…*, `tuffold export study.tuff --to study.xlsx` |
| Methods text | a paragraph for a paper, written from what the project did, with references | Methods sheet, `tuffold methods study.tuff` |
| Spectral library | NIST `.msp`: the best passing spectrum of each peptide ion, b/y-annotated, for SpectraST, Skyline, or NIST MS Search | *Export* → *Spectral library (.msp)…*, `--to library.msp` |
| mzIdentML 1.2 | Tuffold's results (q-values, PEPs, protein groups) for PRIDE or other tools | *Export* → *mzIdentML…*, `--to results.mzid` |
| TSV tables | quant comparisons, ANOVA, and PTM localization | Quant view, `tuffold quant … --export`, `tuffold localize … --export` |

`tuffold export` never overwrites a file; in the app the save dialog asks
first. Excel holds 1,048,575 rows per
sheet; a longer sheet is cut and the Summary says so (the TSVs carry
everything).

## 8. GO terms

Tuffold doesn't download anything, so GO annotation uses files you fetch:

- the ontology, `go-basic.obo` (geneontology.org);
- a GO slim, e.g. `goslim_generic.obo` (geneontology.org);
- annotations for your organisms as GAF, e.g. `goa_human.gaf.gz`
  (EBI GOA).

**Add GO Annotations…** takes all of them at once (a file named `*slim*.obo`
is the slim). Each protein with a UniProt accession gets the slim terms its
annotations reach; they appear in the group header and in the GO filter.
Adding again replaces the earlier annotations.

`tuffold go study.tuff --ontology go-basic.obo --slim goslim_generic.obo --annotations goa_human.gaf.gz`

## 9. Undo

**Undo Last Change** (⌘Z) undoes the last import, search, plex, or sample-group
edit; undoing again redoes it. The previous version is kept beside the
project as `name.tuff.undo`. Threshold changes and validation aren't undone
(they are recomputed from the data, so just change them back).

`tuffold undo study.tuff`

## 10. Coming from Scaffold

Add the `.sf3` or `.sfdb` file. Biosamples, categories, runs, and every engine's
matches come in, and Tuffold validates the engines' own scores; Scaffold's
probabilities are not imported. Expect different numbers from Scaffold's:
the method differs, and Tuffold's FDR is estimated from the decoys.

Some things to know:

- Scaffold may store only the matches it accepted. Tuffold then sees few
  decoys, says so, and its FDR checks little; import the original search
  results if you have them.
- X! Tandem score names were worked out from the values and give their
  position in the file instead.
- From an `.sfdb`, Percolator's PEP and Scaffold's probabilities are left
  out like Scaffold's own; a Scaffold 5 `.sfdb` converted from an `.sf3`
  may hold no decoys at all, and the import says so.
- Peaks come in with the matches. Scaffold keeps at most 300 per spectrum,
  so attaching the runs' mzML, MGF, or `.d` files replaces them with
  complete spectra (and adds MS1 intensities).

Scaffold users will recognise: probability thresholds, exclusive counts,
NSAF and emPAI, Fisher and t-tests, TMT and iTRAQ reporter-ion quant, sample
PCA and correlations, PTM site
probabilities and Ascores, the similarity view, GO filters, stars, and the
Excel export. Not yet: DIA.

## 11. AI assistants

`tuffold-mcp` lets an assistant (Claude Desktop, Claude Code, a local model
through Ollama) open projects, list and compare proteins, and run imports,
searches, and validation. Privacy mode shows it aliases instead of sample,
file, and protein names, and a local model keeps everything on your
computer. What each tool sends and how to connect: [mcp.md](mcp.md).

## 12. When something is left out

Tuffold never invents values and never drops rows without saying so.
Unknown q-values are shown as "—"; proteins without a sequence have no
coverage; an engine without decoys is listed but not validated. Every
import, validation, and export reports what it left out and why, on the
Overview (Findings and notes) and in the command's output.

Problems or questions: open an issue on the repository.
