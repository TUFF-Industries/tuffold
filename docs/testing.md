# Trying to break Tuffold

Tuffold is young, and it has been tested mostly by the person who wrote it, so
real data from other labs is what it needs most. Wrong numbers matter more than
crashes.

## Before relying on a result

Checked on real data so far: FDR and rescoring on several public datasets, and
MS1 protein intensities against FlashLFQ on iPRG 2016. **Not yet checked:**

- iTRAQ, which has only run on synthetic spectra
- TMT SPS-MS3 against MS2 reporters on a real dataset
- MS1 intensities against MaxQuant
- PTM site probabilities against a ground-truth set
- any head-to-head with Scaffold on the same ground truth

Known gaps: no SILAC, no DIA, no Windows build, and no reading Thermo `.raw`
directly (convert with ThermoRawFileParser). The macOS package isn't notarized
([install notes](install.md)).

## Scaffold feature count

The README's badge is this list: 15 of 16, scaled so 95% is the top, because there's
always more to do (15/16 × 95 = 89%). The list is our own, from Scaffold's public
descriptions, and nobody has compared the two side by side, since we don't have a
license. Having a feature isn't the same as its results matching Scaffold's.

| Feature | Tuffold |
|---|---|
| Importing results from many search engines | Yes |
| Spectrum, peptide, and protein thresholds, by FDR or probability | Yes |
| Combining engines | Yes |
| Protein grouping | Yes |
| Annotated spectra and sequence coverage | Yes |
| Spectral counts, NSAF, emPAI, Top 3 | Yes |
| Fisher, t-test, Mann–Whitney, ANOVA | Yes |
| Heat maps | Yes |
| Similarity view | Yes |
| GO terms | Yes |
| PTM site localization | Yes |
| TMT and iTRAQ | Yes |
| MS1 label-free intensity | Yes |
| PCA and sample QC | Yes |
| Excel, spectral library, and mzIdentML exports | Yes |
| SILAC | **No** |

## Comparing with another tool

The fair comparison is one search, validated two ways:

1. Search your spectra with one engine and one FASTA containing decoys.
2. Validate that search in the other tool the way you normally would, and note
   its accepted spectra, peptides, and protein groups at 1% FDR.
3. Import the same search results and the spectra into Tuffold, and validate:

   ```
   tuffold import results.pep.xml run1.mzML -o study.tuff
   tuffold validate study.tuff
   tuffold info study.tuff
   ```

4. Compare the same level at the same FDR.

More accepted isn't automatically better: a tool can accept more by
under-controlling its errors. To test that, search against a FASTA with added
entrapment proteins (from a species that can't be in the sample, or shuffled)
and count how many accepted identifications hit them, in each tool.

Other useful checks: known mixtures or spike-ins for fold changes, and the same
data run twice, which should give the same answer.

## Reporting

Use the [issue form](https://github.com/TUFF-Industries/tuffold/issues/new/choose).
**Save Diagnostics…** in the Tuffold menu (or `tuffold diagnostics study.tuff
--to report.txt`) writes a text file with the version, system, models, project
counts, and any panic logs. It has no sample, run, or protein names, but read it
before sending it, since a crash message can name a file. Nothing is ever sent
automatically; crash logs stay in `~/Library/Logs/Tuffold` (macOS) or
`~/.local/state/tuffold/logs` (Linux).

Useful things to throw at it: a file from an engine or version it hasn't seen,
huge projects, odd modifications, empty or truncated files, samples with no
identifications, and anything it refuses that you think it should read.
