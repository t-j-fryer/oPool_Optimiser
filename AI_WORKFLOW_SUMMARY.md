# AI Workflow Summary

Last updated: 2026-03-06

This file is a handoff guide for AI coding tools and new contributors working on `oPool_Optimiser`.

## Project Purpose

`oPool_Optimiser` packages a notebook workflow that:
1. Reverse-translates amino acid sequences and codon-optimizes DNA.
2. Assigns Golden Gate overhangs and fragments for oPool assembly.
3. Adds primer pairs and Type IIS sites to produce order-ready oligos.
4. Exports ordering and reference outputs.

## Canonical Repo Layout

- `notebooks/oPool_Cloning_Notebook.ipynb`: main end-to-end workflow.
- `data/AAseq_dTF001_dTF016.csv`: example amino acid input (renamed from original merged CSV).
- `data/overhangs.csv`: example overhang inventory.
- `data/orthogonal_oligos.csv`: example primer inventory.
- `outputs/`: generated outputs (git-ignored except `.gitkeep`).
- `requirements.txt`: pip environment lock.
- `environment.yml`: conda environment option.

## Machine-Readable Context

```yaml
project:
  name: oPool_Optimiser
  notebook: notebooks/oPool_Cloning_Notebook.ipynb
runtime:
  python: "3.11"
  kernel_display_name: "Python 3 (opool-cloning)"
  kernel_name: "python3"
env_options:
  pip_requirements: requirements.txt
  conda_env: environment.yml
example_inputs:
  aa_sequences_csv: data/AAseq_dTF001_dTF016.csv
  overhangs_csv: data/overhangs.csv
  orthogonal_oligos_csv: data/orthogonal_oligos.csv
default_outputs_dir: outputs/
path_policy:
  absolute_paths_allowed: false
  use_repo_relative_paths: true
git:
  default_branch: main
  backup_branch_pre_replace: backup/pre-replace-2026-03-06
```

## Setup Commands

### Option A: venv (recommended)

```bash
python3.11 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
python -m ipykernel install --user --name opool-cloning --display-name "Python 3 (opool-cloning)"
jupyter lab
```

### Option B: conda

```bash
conda env create -f environment.yml
conda activate opool-cloning
python -m ipykernel install --user --name opool-cloning --display-name "Python 3 (opool-cloning)"
jupyter lab
```

## Expected Input Contracts

### `data/AAseq_*.csv`

- Intended format in this notebook: 2 columns, no header.
- Column 1: sequence name (unique ID).
- Column 2: amino-acid sequence (1-letter code).
- The code loads this using `pd.read_csv(..., names=["name", "aa_seq"])`.
- If you include a header row, it will be treated as data unless code is changed.

### `data/overhangs.csv`

- Parsed as a single cell containing comma-separated 4-nt overhangs.
- Must not include forbidden vector overhangs for the selected mode.

### `data/orthogonal_oligos.csv`

- Parsed with `header=None`.
- First column is primer name/ID.
- Last column is primer sequence (5'->3').
- Rows are used as `[F, R, F, R, ...]` in `unique_pairs` mode.

## Notebook Workflow and Key Settings

The notebook has multiple configurable blocks. These are the key variables AI tools should preserve and reason about.

### 1) Gene Optimization Block

Defaults:
- `INPUT_PDB = "dTF001_dTF016"`
- `input_path = "data/AAseq_dTF001_dTF016.csv"`
- `output_path = "outputs/dTF001_dTF016_Optimised.csv"`
- `AVOID_SEQS = ["GGTCTC", "AAAAA", "GGGGG", "CCCCC", "TTTTT"]`
- `GC_MIN = 0.30`
- `GC_MAX = 0.70`
- `GC_WINDOW = 50`
- `CODON_SPECIES = "e_coli"`
- `CODON_METHOD = "match_codon_usage"`

Purpose:
- Reverse-translate amino acids to DNA.
- Remove motifs and optimize codons with DnaChisel.

### 2) Pool Assignment Block

Defaults:
- `OPOOL_LENGTH = 300`
- `VECTOR_OH1_MODE = "fixed"` (`"fixed"` or `"ser_windows"`)
- `VECTOR_OVERHANG_1_FIXED = "GGCA"`
- `VECTOR_OVERHANG_2 = "TCGG"`
- `VEC1_WINDOW0_FRAG1_PENALTY_NT = 2`
- `SHORT_POOL_MAX_SIZE = None`
- `STRIP_NTERM_MET = True`
- `WRITE_STRIP_LOG = True`
- `STRIP_LOG_CSV = "outputs/dTF001_dTF016_stripped_ATG_log.csv"`
- `OVERHANGS_FILE = "data/overhangs.csv"`
- `INPUT_FILE = output_path`
- `OUTPUT_FILE = "outputs/dTF001_dTF016_Assigned.csv"`
- `UNASSIGNED_LOG = "outputs/dTF001_dTF016_unassigned.csv"`

Purpose:
- Split constructs into fragment designs based on pool length and overhang compatibility.

### 3) Primer + Type IIS Addition Block

Defaults:
- `OPOOL_LENGTH = 300`
- `ADD_STUFFER = True`
- `STUFFER_GC_MIN = 0.40`
- `STUFFER_GC_MAX = 0.60`
- `STUFFER_MAX_HOMOPOLYMER = 4`
- `STUFFER_MAX_TRIES = 10000`
- `STUFFER_SEED = 123`
- `TYPEIIS_CUT_SITE = "GGTCTCA"` (BsaI-compatible pattern including chosen `N`)
- `PRIMER_ASSIGN_MODE = "unique_pairs"` (`"unique_pairs"` or `"combinatorial"`)
- `PRIMER_START_AT = None`
- `BLOCK_BASE = 1`
- `PRIMER_CSV = "data/orthogonal_oligos.csv"`
- `OUTPUT_UNUSED_PRIMERS = "outputs/orthogonal_oligos_unused.csv"`
- `WRITE_UNUSED_PRIMERS = True`
- `GENERATE_PAIR_TABLE_ONLY = False`
- `ALL_PAIRS_CSV = "outputs/orthogonal_oligos_pairs_ALL.csv"`
- `PAIRS_CSV = "outputs/orthogonal_oligos_pairs_unused.csv"`
- `OUTPUT_UNUSED_PAIRS = "outputs/orthogonal_oligos_pairs_unused.csv"`
- `INPUT_ASSEMBLY_CSV = OUTPUT_FILE`
- `OUTPUT_WITH_PRIMERS = "outputs/dTF001_dTF016_FULL_INFO.csv"`
- `OUTPUT_FASTA = "outputs/dTF001_dTF016_references.fasta"`
- `OUTPUT_FRAGMENTS_CSV = "outputs/dTF001_dTF016_oPool_Order_Fragments.csv"`
- `MAKE_OUTPUT_DIRS = True`
- `FALLBACK_VECTOR_OH1 = "GGCA"`
- `FALLBACK_VECTOR_OH2 = "TCGG"`
- `VEC1_NEEDS_2NT_SOURCE_LABELS = {"W0"}`

Purpose:
- Add vector overhang context, primer handles, Type IIS sites, and optional stuffer.
- Produce files needed for synthesis orders and references.

### 4) Optional Merge Block

Defaults:
- `ROOT_DIR = Path("outputs")`
- `OUT_CSV = ROOT_DIR / "opTF003_oPool_Order_Fragments.csv"`

Purpose:
- Merge multiple `*oPool_Order_Fragments.csv` files into one table.

## Core Outputs to Expect

- `outputs/dTF001_dTF016_Optimised.csv`
- `outputs/dTF001_dTF016_Assigned.csv`
- `outputs/dTF001_dTF016_FULL_INFO.csv`
- `outputs/dTF001_dTF016_references.fasta`
- `outputs/dTF001_dTF016_oPool_Order_Fragments.csv`
- `outputs/dTF001_dTF016_unassigned.csv` (only when assignment fails for some entries)
- `outputs/dTF001_dTF016_stripped_ATG_log.csv` (if enabled)

## Recommended Operator Workflow

1. Create environment and install kernel.
2. Open notebook in Jupyter Lab and select `Python 3 (opool-cloning)`.
3. Confirm top-level variables point to intended `data/` and `outputs/` files.
4. Run cells in order from top to bottom.
5. Inspect generated CSV/FASTA files in `outputs/`.
6. Commit notebook changes and any intentional input updates.
7. Do not commit generated files in `outputs/` unless explicitly required.

## AI Tool Guardrails

When an AI tool edits this repo, it should:
1. Keep all file paths repo-relative (`data/...`, `outputs/...`).
2. Avoid reintroducing machine-specific absolute paths.
3. Preserve input contracts unless intentionally migrating formats.
4. Keep kernel/environment instructions synchronized across README and config files.
5. Avoid changing biology-critical defaults silently; surface diffs explicitly.
6. Clear large notebook outputs before commit unless needed for review.
7. Prefer additive docs/scripts over destructive rewrites.

## FAQ

Q: Which settings most often need changing for a new experiment?
A: `input_path`, `output_path`, `OPOOL_LENGTH`, vector overhang settings, and primer mode/settings.

Q: How do I skip gene optimization if I already have optimized DNA?
A: Point `INPUT_FILE` (pool assignment block) to your pre-optimized CSV and ensure expected columns exist.

Q: My AA input has a header row. What should I do?
A: Either remove the header row or change the load call to use `header=0` and explicit column mapping.

Q: What causes unassigned sequences?
A: Overhang compatibility constraints, vector overhang conflicts, fragment-length constraints, or restricted available overhang sets.

Q: When should I use `ser_windows` mode?
A: Use it when per-gene unique `VectorOH1` assignment from Ser-Ser window-derived sets is required.

Q: Why is `TYPEIIS_CUT_SITE` set to `GGTCTCA` instead of `GGTCTCN`?
A: The script requires concrete DNA bases, so a specific base is chosen for `N`.

Q: Where should generated outputs go?
A: `outputs/` only. The repository is configured to ignore this directory in git by default.

Q: How can I recover the old GitHub repo content from before replacement?
A: Use branch `backup/pre-replace-2026-03-06` as the restore source.

## Troubleshooting Checklist

1. Kernel mismatch: verify Python 3.11 and correct Jupyter kernel selected.
2. Missing package errors: reinstall with `pip install -r requirements.txt`.
3. File not found: check that paths are relative to repo root and files exist in `data/`.
4. Invalid primer format: ensure first column is name and last column is sequence.
5. Unexpected AA parsing: confirm no header row unless code is updated accordingly.
6. Oligo length assertion failures: adjust `OPOOL_LENGTH`, stuffer options, or fragment settings.
7. Pair assignment failures: verify primer inventory size is sufficient for number of blocks.

## Maintainer Notes

- This handoff file is intended to be the first document an AI tool reads before editing code or notebook settings.
- If defaults or file schema change, update this file in the same commit.
