# sf3 analysis tools

Python 3 (standard library only) for working out Scaffold's `.sf3` layout
from sample files. Findings: `docs/sf3-format.md`. Never commit `.sf3`
files or what these tools extract.

    python3 -I dump.py FILE.sf3 OUTDIR          # unpack the gzip blobs
    python3 -I explore.py tools/sf3 OUTDIR/Experiment
    python3 -I desc.py tools/sf3 OUTDIR         # record-file descriptors
    python3 -I profile.py OUTDIR/PeptideIdentification0 160 4 35
    python3 -I strings.py OUTDIR/SpectrumName0
    python3 -I tree.py tools/sf3 OUTDIR/ModificationTypes 6   # any Java stream

The checks below run from inside OUTDIR (`-I` keeps it off the import
path; the tools directory is passed explicitly):

    python3 -I $TOOLS/mods.py $TOOLS    # modifications sit on allowed residues
    python3 -I $TOOLS/mass.py $TOOLS    # peptide mass vs [M+H]+ for every PSM
    python3 -I $TOOLS/prot.py $TOOLS    # PSM peptide inside its protein's sequence
    python3 -I $TOOLS/lookup.py $TOOLS  # protein -> peptides map vs sequences
    python3 -I $TOOLS/scores.py $TOOLS  # score field ranges per engine
    python3 -I $TOOLS/rankauc.py $TOOLS # target/decoy AUC ranking by each field
    python3 -I $TOOLS/ident.py $TOOLS   # exact identities between score fields
    python3 -I $TOOLS/rank1.py $TOOLS   # what the offset-147 flag picks
    python3 -I $TOOLS/peaks.py $TOOLS   # stored peaks: chains, overlap, counts

`jser.py` reads Java serialization streams; `bbser.py` the block-based
object dialect and `bbload.py` its record files; `modlib.py` reads
fixed-block record files across parts and the modification table.
