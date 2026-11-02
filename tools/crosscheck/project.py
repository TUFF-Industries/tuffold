"""Run mokapot on a validated Tuffold project's rank-1 PSMs, with the
feature set Tuffold's rescoring uses (the engine scores reported on every
PSM, -log10 for e-values; charge indicators, |mass error| to the nearest
isotope peak, length, missed cleavages, modifications), and compare the
1% FDR lists at spectrum and peptide level. Read-only on the project.

Usage (needs mokapot 0.10 and pandas; NumPy 2 is patched over below):
    python tools/crosscheck/project.py PROJECT.tuff [R]

With R (from `entrapment_fasta.py`), mokapot's lists also get the entrapment
estimate `entrapment.py` gives Tuffold's, at 0.5%, 1% and 2%.

Validate the project with `tuffold validate` first. Several engines in one
project are pooled here, so use it on single-engine projects.
"""

import sqlite3
import sys

import numpy as np, pandas as pd
np.float_ = np.float64  # mokapot 0.10 predates NumPy 2
import mokapot

path = sys.argv[1]
db = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
psms = pd.read_sql_query("""
    select p.id psm, p.spectrum_id spectrum, p.charge, p.calculated_mz calc, p.experimental_mz exp,
           pe.proforma, pe.sequence, d.is_decoy decoy
    from psms p join peptides pe on pe.id = p.peptide_id join psm_decoys d on d.psm_id = p.id
    where p.rank = 1""", db)
scores = pd.read_sql_query("select psm_id psm, name, value from psm_scores", db)
scores = scores[scores.psm.isin(psms.psm)].pivot(index="psm", columns="name", values="value")
scores = scores.dropna(axis=1)  # only scores on every PSM, as Tuffold
for name in list(scores.columns):
    if any(m in name.lower() for m in ("expect", "evalue", "e-value", "pvalue", "p-value")):
        scores[f"-log10({name})"] = -np.log10(scores.pop(name).clip(lower=1e-300))
df = psms.merge(scores, left_on="psm", right_index=True)
for z in sorted(df.charge.clip(upper=5).unique()):
    df[f"charge{z}"] = (df.charge.clip(upper=5) == z).astype(float)
iso = 1.003354835
df["masserr"] = np.min([np.abs((df.exp - (df.calc + k * iso / df.charge)) / df.calc * 1e6)
                        for k in range(-1, 4)], axis=0)
df["length"] = df.sequence.str.len()
df["missed"] = df.sequence.apply(lambda s: sum(1 for a, b in zip(s, s[1:]) if a in "KR" and b != "P"))
df["mods"] = df.proforma.str.count(r"\[")
features = [c for c in df.columns if c not in ("psm", "spectrum", "charge", "calc", "exp", "proforma", "sequence", "decoy")]
features = [c for c in features if df[c].nunique() > 1]
df["target"] = df.decoy == 0
print(f"{path.rsplit('/',1)[-1]}: {len(df)} rank-1 PSMs ({int(df.decoy.sum())} decoy); features: {', '.join(features)}")
data = mokapot.LinearPsmDataset(psms=df, target_column="target", spectrum_columns=["spectrum"],
                                peptide_column="proforma", feature_columns=features, copy_data=False)
results, _ = mokapot.brew(data, rng=1)
mk = results.confidence_estimates["psms"]
mk_spec = set(mk.loc[mk["mokapot q-value"] <= 0.01, "spectrum"])
mkp = results.confidence_estimates["peptides"]
mk_pep = set(mkp.loc[mkp["mokapot q-value"] <= 0.01, "proforma"])
tu = pd.read_sql_query("select spectrum_id spectrum, q_value q, is_decoy from spectrum_results", db)
tu_spec = set(tu.loc[(tu.q <= 0.01) & (tu.is_decoy == 0), "spectrum"])
tp = pd.read_sql_query("select pe.proforma, pr.q_value q, pr.is_decoy from peptide_results pr join peptides pe on pe.id = pr.peptide_id", db)
tu_pep = set(tp.loc[(tp.q <= 0.01) & (tp.is_decoy == 0), "proforma"])
print(f"spectra at 1%: mokapot {len(mk_spec)}, Tuffold {len(tu_spec)}, both {len(mk_spec & tu_spec)}, "
      f"mokapot only {len(mk_spec - tu_spec)}, Tuffold only {len(tu_spec - mk_spec)}")
print(f"peptides at 1%: mokapot {len(mk_pep)}, Tuffold {len(tu_pep)}, both {len(mk_pep & tu_pep)}")
if len(sys.argv) > 2:
    r = float(sys.argv[2])
    accessions = {}
    for peptide, accession in db.execute("""SELECT DISTINCT pe.proforma, pr.accession FROM psms ps
            JOIN peptides pe ON pe.id = ps.peptide_id JOIN psm_proteins x ON x.psm_id = ps.id
            JOIN proteins pr ON pr.id = x.protein_id WHERE ps.rank = 1"""):
        accessions.setdefault(peptide, []).append(accession)
    entrap = {pf for pf, accs in accessions.items()
              if (real := [a for a in accs if not a.startswith(("reverse_", "rev_"))])
              and all(a.startswith("ENTRAP_") for a in real)}
    for q in (0.005, 0.01, 0.02):
        accepted = mk.loc[mk["mokapot q-value"] <= q].drop_duplicates("spectrum")
        ne = int(accepted.proforma.isin(entrap).sum())
        print(f"mokapot q≤{q:g}: {len(accepted)} spectra, entrapment estimate {ne * (1 + 1 / r) / max(len(accepted), 1):.2%}")
