# Proteome Discoverer result files: what Tuffold reads

Tuffold imports Thermo Proteome Discoverer **2.x and 3.x** results: the
consensus `.pdResult` and the per-search `.msf` files those versions write.
Both are SQLite databases with the same PSM tables. PD **1.x** `.msf` files
have their own layout (a `Peptides` table); see *PD 1.x* below.

Everything below was checked on two public PD 3.1 files: an `.msf`
(PRIDE PXD082486, Sequest HT, 28,604 PSMs) and a consensus `.pdResult`
(PRIDE PXD083435, Sequest HT, 71,769 PSMs of which 23,956 decoys; it
validates to 712 protein groups at 1% FDR). Table and column names
were also cross-read against open-source readers (CompOmics psm_utils,
pyproteome); no PD code was used.

## Tables

| table | what Tuffold takes |
|---|---|
| `SchemaInfo` | `SoftwareVersion` of the `Result` row → the software version |
| `WorkflowInputFiles` | `FileID` → `FileName` (a raw file path) → run name. A `.pdResult` also lists its input `.msf` files here, so only files that PSMs name become runs |
| `TargetPsms`, `DecoyPsms` | one row per PSM: `WorkflowID` + `PeptideID` key; `IdentifyingNodeName`/`IdentifyingNodeTypeName` (the engine); `Sequence`; `Charge`; `MassOverCharge` (experimental m/z); `TheoreticalMass` (theoretical **MH+**); `RetentionTime` (**minutes**); `FirstScan`, `LastScan`; `SpectrumFileId`; `SearchEngineRank` (the engine's rank; `Rank` is PD's); score columns |
| `ProcessingNodeScores` | `ScoreName` (= the PSM table's column), `DisplayName`, `IsMainScore` |
| `FoundModifications` | `ModificationID`, `Name`, `DeltaMonoisotopicMass`, `PositionType`, `UnimodAccession` |
| `TargetPsmsFoundModifications`, `DecoyPsmsFoundModifications` | PSM key → modification id and `Position` |
| `TargetProteins`, `DecoyProteins` | `UniqueSequenceID`, `Accession`, `Description`, `Sequence`, `FastaTitleLines` |
| `TargetProteinsTargetPsms`, `DecoyProteinsDecoyPsms` | protein ↔ PSM |

## Conventions (checked)

- Modification `Position` is **1-based** on residues (`S2(Phospho)` →
  position 2). `PositionType` 0 = residue; 1 and 3 = (protein) N-terminus,
  with `Position` 1; 2 and 4 = C-terminus.
- Sequence + modification deltas + proton equals `TheoreticalMass` for
  every PSM of both samples (100,373, targets and decoys, to 0.01 Da). The importer applies this
  check to every PSM and refuses the file on any failure.
- Decoy PSMs live in `DecoyPsms` (PD's target-decoy search); their proteins
  in `DecoyProteins` **under the target's accession**, so Tuffold names them
  `DECOY_<accession>` with an explicit decoy flag.
- Accessions need not be unique: the sample's FASTA headers (`>protein 1
  sequence`, `>protein 2 sequence`, ...) all give `protein`. Proteins that
  share an accession are named by their FASTA title line.
- Score columns are engine-specific columns of the PSM tables (`XCorr` for
  Sequest HT). Tuffold takes those named in `ProcessingNodeScores` plus
  `DeltaCn` and `DeltaScore`, when both PSM tables have them; PD's
  `MatchConfidence` and validator outputs are not imported.
- A workflow may run several nodes of one engine type (the sample: Sequest
  HT nodes 6 and 8 over the same 9,969 spectra, decoys only from node 8);
  each node is its own engine, named by node.

## Target and decoy competition

PD searches the target and decoy databases separately and stores each
search's own ranks (`SearchEngineRank`), so a spectrum usually has a rank-1
target and a rank-1 decoy. Tuffold imports both as rank 1, and validation
competes them per spectrum on the rescored score (target-decoy competition
on the final score, as Percolator does on separate searches with
`--post-processing-tdc`). They are not ranked against each other at import.

That was tried and measured (2026-10-27, PRIDE PXD071927,
`TcNMT_TMT6plex_6-16-2025.msf`, PD 2.5, Sequest HT, 24 runs, 106,726
spectra; 103,440 `TargetPsms` and 88,132 `DecoyPsms` rows with
`SearchEngineRank` 1). Ranking by XCorr at import, as a concatenated
search would, passes the same 1,838 spectra without rescoring. With
rescoring it passes 3,064 against 3,534 (peptides 509 against 687,
protein groups 131 against 143), because a spectrum whose decoy outscored
its target is gone before rescoring can weigh them (inferred: the file has
no truth). On seeded synthetic separate searches both ways keep the true
FDP under 1% (`crates/tuffold-validate/tests/separate_search.rs`).

In that file 7,757 spectra have a target and a decoy with exactly the same
XCorr: 949 of them the same peptide in both searches, the others different
sequences of that precursor, many differing only by I/L or residue order.
Validation gives such a tie to the decoy. The `.sfdb` Scaffold made from
the same search, as Tuffold imports it, ranks targets and decoys together
and breaks equal scores by PSM id, which favours the target there. That is
why it passes more (2,452 spectra on XCorr alone, 3,947 rescored without
deep features, which the `.msf` lacks peaks for): ranking the `.msf`
jointly with ties to the target gives the same rank 1s and 2,452.

PD's `DeltaCn` is measured from the same search's rank 1, so it is 0 for
every rank-1 PSM (all 191,572 in that file).

## PD 1.x

Checked on PRIDE PXD001296 (`AF_secretome_B21-22_1.msf`, PD 1.4.0.288,
Mascot and Sequest HT nodes plus Percolator; 21,702 PSMs).
`crates/tuffold-io/src/pd/v1.rs`.

| table | what Tuffold takes |
|---|---|
| `Peptides`, `Peptides_decoy` | `ProcessingNodeNumber` + `PeptideID` key, `SpectrumID`, `SearchEngineRank`, `Sequence`. Target and decoy searches are separate, so a spectrum has a rank 1 in each; they compete in validation (see *Target and decoy competition*) |
| `ProcessingNodes` | `FriendlyName` of each node with peptides → the engine |
| `PeptideScores(_decoy)`, `ProcessingNodeScores` | score id → name (`XCorr`, `IonScore`): each engine stores only its main score |
| `SpectrumHeaders` + `MassPeaks` | scans, `RetentionTime` (minutes), `Charge`, `Mass` = precursor **MH+**; `MassPeaks.Mass` = precursor m/z, `FileID` |
| `FileInfos` | `FileID` → `FileName` (raw path) → run |
| `PeptidesAminoAcidModifications(_decoy)` | modification id and **0-based** `Position` |
| `PeptidesTerminalModifications(_decoy)` | terminal modification id; its `PositionType` 1, 3 = N-term, 2, 4 = C-term |
| `AminoAcidModifications` | name, `DeltaMass`, `PositionType`, `UnimodAccession` |
| `PeptidesProteins(_decoy)`, `Proteins`, `ProteinAnnotations` | protein links, sequences, FASTA title lines (accession = first word) |

Checked: every one of 24,013 residue modifications (targets and decoys)
sits on a residue its type allows when positions are 0-based (7,869 don't
when read 1-based). Fixed modifications are listed per PSM. 1.x stores no
theoretical mass, so the importer checks each PSM against its precursor
MH+ (isotope peaks −1..+3, 0.05 Da): all 21,702 fit within 20 ppm.

**Protein links are incomplete.** 1.x stores only proteins some target
matched (5,460 of the FASTA's 13,501) and links a decoy only when its
target is stored: 4,564 of 10,798 decoy peptides (and 709 targets, mostly
Mascot's) have no link. A PSM without proteins counts as a target in
Tuffold, so each is linked by sequence to stored proteins containing its
peptide (decoys: reversed, whole or but the last residue; 954 found), and a
decoy still without one is filed under a placeholder decoy protein
`DECOY_unstored:<peptide>` (4,062). Before this, 1,700 Sequest decoys
counted as targets.

Percolator's q-values and PEPs (`CustomDataPeptides`) are not imported.
The sample is a weak search (target beats decoy on 52% of spectra; XCorr
means 0.92 vs 0.88): PD's Percolator passes 95 Sequest PSMs at 1% and 153
spectra at 5%; Tuffold passes none at 1% (with the +1 correction 1% needs
100 targets before the first decoy; it reaches 1.8%) and 114 spectra at 5%.

## Still unknown

- PD 3.x `.pdResultDetails` spectra.
- Peak data (`MassSpectrumItems`/blobs) — not read; attach mzML instead.
