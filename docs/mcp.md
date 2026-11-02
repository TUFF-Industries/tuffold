# The MCP server and your data

Tuffold works fully offline. The desktop app and the `tuffold` command make
no network connections, and nothing in a project ever leaves the computer
unless you connect an AI assistant to `tuffold-mcp`. This page says exactly
what the assistant then sees, how to show it less, and how to keep even the
assistant on your own machine.

## What leaves, and where it goes

`tuffold-mcp` itself also makes no network connections. It answers the
assistant's tool calls over standard input and output, and the assistant
app (Claude Desktop, Claude Code, an Ollama bridge, …) sends those answers to
whichever model it uses. If that model runs in the cloud, everything below
goes to its provider; if it runs locally, nothing leaves.

Only replies to tool calls are sent, and only for the tools the assistant
calls. What each reply holds:

| Tool | Reply |
|---|---|
| `project_summary` | counts; run, sample, category, and plex names; import file paths, software and versions; validation settings, PTM localization summary (scored and localized PSM/site counts, decoy-AA FLR), and notes |
| `list_protein_groups` | accessions, descriptions, peptide and spectrum counts (all and exclusive) per sample, q-values, probabilities, coverage, GO slim terms |
| `protein_group` | member accessions and descriptions, peptides (sequence, modifications, localized peptide, site probabilities, and Ascore), sequence length, coverage |
| `peptide_spectra` | run and sample names, scan numbers, charge, precursor m/z, retention time, q-values, localized peptide, site probabilities, and Ascore |
| `compare_categories`, `compare_all_categories` | accessions, category and sample names, means, fold changes, p- and q-values (counts, MS1 intensity, or TMT reporters) |
| `set_categories` | sample and category names |
| `set_plex`, `remove_plex`, `quantify_reporters` | plex name, reagent set, run and channel sample names, whether reporters are corrected for the lot's isotope impurities (`set_plex` takes the certificate's local path; its contents are not sent back), reference channels and whether several plexes were bridged through them, the precursor purity cut and spectra left out by it, normalization shifts, PSM counts, better-fitting reagent set if any, per-channel reporter intensities, and optional two-category comparison |
| `match_between_runs` | the run names and counts of precursors looked for, transferred, and won by decoys, and of damaged MS1 scans left out |
| `import_data`, `search_spectra` | the paths the assistant gave, the search settings (enzyme, tolerances in ppm or Da, modifications), counts (including MS3 scans stored for SPS-MS3 reporters and damaged MS1 scans left out), warnings (which can quote file names) |
| `validate`, `set_thresholds` | settings, counts, engine names, how peptides were scored, PTM localization summary, notes |
| `export_workbook`, `export_library`, `export_mzid` | the file's path and what went in (the file itself stays on disk) |
| `methods_text` | the methods paragraph: software, settings, counts, references |
| `undo_last_change`, `create_demo` | what was undone; the demo's validation |

Never sent by any tool: spectra and peaks, protein sequences (only their
length and coverage), raw or mzML files, FASTA files, or anything outside the
project the assistant names. Error messages can quote names and paths.

The assistant also sees what you type to it and the paths it is asked to
open: a path like `~/patient-cohort/run1.mzML` tells it that much
whatever mode the server is in.

## Privacy mode

Start the server with `--private` and the assistant sees aliases instead of
names; the server maps them back when the assistant uses them, so every tool
works as before:

| Level | Hidden |
|---|---|
| `--private` (or `--private=names`) | sample, category, run, plex, spectra-file, and import names, as `sample-1`, `group-2`, `run-3`, `plex-1`, `file-1`, `import-2` |
| `--private=proteins` | the above, plus protein accessions (`protein-7`) and descriptions (left out) |

`TUFFOLD_MCP_PRIVATE=names` (or `proteins`) does the same where a client
can't pass arguments. The level is set by whoever starts the server; no tool
can change it.

Every string in every reply is rewritten, error messages included: a value
that is a name becomes its alias, and names of three characters or more are
replaced inside text too, at word boundaries. Names of one or two
characters (a category called `A`) are hidden only as whole values. Aliases
are numbered per project the first time the server sees each name and stay
the same while it runs, so `sample-2` means one sample throughout a
conversation (a restarted server may number them differently).

Still visible in any mode: numbers, software and settings, and **peptide
sequences**, which can identify a protein. `proteins` hides which protein a
group is called, not what it is; when that matters, use a local model.

## A local model

With a model running on your own computer nothing leaves at all.
`tools/local-model/ollama_mcp.py` connects an [Ollama](https://ollama.com)
model to `tuffold-mcp` (Python standard library only; it refuses an Ollama
host other than this computer):

```sh
ollama pull qwen2.5                      # any model with tool calling
cargo build --release -p tuffold-mcp
python3 tools/local-model/ollama_mcp.py --model qwen2.5 \
    --server target/release/tuffold-mcp --private names \
    "Create the demo project at /tmp/demo.tuff and summarize it."
```

It prints every tool call and reply, so you can see what the model asked for
and got. Checked with Ollama 0.30.9 and `qwen2.5` (7B) on the demo project:
the model created it, summarized it, and listed protein groups, seeing only
aliases. Small models use the tools fine but summarize loosely; larger ones
(`qwen2.5:14b` and up) answer better.

Any MCP client that can use a local model works the same way: point it at
`tuffold-mcp` with the arguments above.

## Cloud assistants

With a cloud model, the replies above go to its provider under your
account's terms. Before connecting real data, check what those terms say
about training on and retaining your inputs: business, API, and enterprise
agreements usually exclude training and can offer zero data retention;
personal accounts may have a setting for it. Privacy mode reduces what is
identifiable; it does not make a cloud model local.

Claude Code, with privacy mode:

```sh
claude mcp add tuffold -- "$PWD/target/release/tuffold-mcp" --private
```

Claude Desktop (`claude_desktop_config.json`):

```json
{
  "mcpServers": {
    "tuffold": {
      "command": "/absolute/path/to/tuffold-mcp",
      "args": ["--private=proteins"]
    }
  }
}
```
