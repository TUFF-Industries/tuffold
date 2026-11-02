# Agent guide

Tuffold is a free, open-source (Apache-2.0) desktop app for validating and
viewing MS/MS proteomics results, in the spirit of Scaffold. It's all Rust: a
headless engine with three thin front ends, a GPUI desktop app, a CLI, and an
MCP server. It implements its own FDR, rescoring, protein inference, and
quantitation, checked against established tools.

Start with [README.md](README.md) for what it does and [docs/guide.md](docs/guide.md)
for how it's used.

## Layout

| Path | What's there |
|---|---|
| `crates/tuffold-core` | Data model and chemistry; no I/O |
| `crates/tuffold-io` | Importers and spectra readers; each produces an `ImportBatch` |
| `crates/tuffold-validate` | Target-decoy FDR, rescoring, engine merging, protein inference, PTM localization, as pure functions |
| `crates/tuffold-quant` | Spectral counts, MS1 intensity, reporter ions, statistics, sample QC, as pure functions |
| `crates/tuffold-search` | The built-in search, over Sage's `sage-core` |
| `crates/tuffold-predict` | AlphaPeptDeep's models for deep-learning rescoring features |
| `crates/tuffold-store` | The `.tuff` project file (SQLite): schema, migrations, queries |
| `crates/tuffold-engine` | Ties the crates into import, validate, quantify, and export |
| `crates/tuffold-cli`, `crates/tuffold-mcp` | The `tuffold` command and the MCP server |
| `app` | The desktop app: view models in `project.rs`, `browse.rs`, `quant.rs`; rendering in `views/`; design tokens in `theme.rs` |
| `docs` | User guide, install notes, how to test, MCP, and notes on the file formats read |
| `tools` | Cross-checks against other tools, packaging, format analysis |

Front ends call the engine crates directly, with no IPC layer, and long work
stays off the UI thread. Keep the engine headless.

## Commands

```
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo test -p tuffold-app --test snapshots
```

Plain `cargo test` covers only the headless crates; the app needs GPUI's
platform libraries. The snapshot test renders the main screens in light and dark
to `target/tmp/snapshots/`, on macOS. Look at them after a UI change.

## Rules

- **Never fabricate data.** Unknown q-values stay `None`. Missing sequences, PSMs,
  or groups are reported, never invented, and rows are never silently dropped.
- **Validation is held to its error rate.** Changes must keep
  `crates/tuffold-validate/tests/calibration.rs` passing and stay deterministic,
  with no unseeded randomness. An engine without decoys or a usable score is
  reported unscored, not ranked by guesswork. Opening a project never writes to it.
- **Spectrum identity is `(run, SpectrumKey)`.** Every importer derives keys the
  same way, or engines and spectra files won't meet.
- **Diagnostics are specific.** `IssueCode` values are a stable vocabulary. A message
  names the file, line, field, and what it means for the loaded data.
- **Fixtures are checked by hand.** Expectations are derived from the fixture files
  and written inline, never by running the code under test. Files in
  `fixtures/third_party/` stay byte-identical to upstream.
- **Keep data out of the repository:** no `.tuff` files, datasets, generated exports,
  or machine-specific paths.
- **Scaffold:** reading its files for import is fine. Never write its formats or copy
  its branding or assets. Tuffold validates engines itself rather than importing
  Scaffold's probabilities.
- **The UI is a desktop instrument:** blue only for real actions, mono for data and
  sans for words, no filler text, and exclusive choices on a track like the tabs.
- **No TypeScript or JavaScript front ends.** Solve GPUI problems in Rust.

Keep `docs/guide.md` and the README current with what you change.
