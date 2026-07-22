# oPool Optimiser

Fast, user-friendly workflows for optimizing genes and preparing pooled Golden Gate cloning oligos.

For AI assistant handoff, settings, workflow conventions, and FAQ, see `AI_WORKFLOW_SUMMARY.md`.

## Repository Layout

- `notebooks/oPool_Cloning_Notebook_Fast_Pool_Assignment.ipynb`: modular notebook with the faster long-gene pool search
- `notebooks/oPool_Cloning_Notebook_Simple.ipynb`: edit one input cell, then choose **Run All**
- `scripts/opool_cli.py`: terminal command-line interface
- `scripts/opool_workflow.py`: shared implementation used by the CLI and simple notebook
- `data/orthogonal_oligos.csv`: example primer inventory
- `data/overhangs.csv`: example overhang list
- `data/AAseq_dTF001_dTF016.csv`: example AA input table (copied from `dTF001_dTF016_merged.csv`)
- `outputs/`: generated outputs from notebook runs

## Python Environment + Kernel

### Option A: `venv` (recommended)

```bash
python3.11 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
python -m ipykernel install --user --name opool-cloning --display-name "Python 3 (opool-cloning)"
```

Then in Jupyter, choose kernel: `Python 3 (opool-cloning)`.

### Option B: Conda/Mamba

```bash
conda env create -f environment.yml
conda activate opool-cloning
python -m ipykernel install --user --name opool-cloning --display-name "Python 3 (opool-cloning)"
```

## Simplest notebook workflow

```bash
jupyter lab
```

Open `notebooks/oPool_Cloning_Notebook_Simple.ipynb`, edit its single **User inputs** cell, and choose **Run All**.

## Terminal workflow

The input can be either a two-column amino-acid CSV or an existing `*_Optimised.csv`; the CLI detects which one it received.

```bash
source .venv/bin/activate
python scripts/opool_cli.py \
  --input "/path/to/opTF010_Optimised.csv" \
  --overhangs "/path/to/overhangs_bgal.csv" \
  --opool-length 350 \
  --vector-oh1 TATG \
  --vector-oh2 GGAT \
  --genes-per-subpool 1
```

Only `--input` is universally required. Repository overhang and primer inventories, 250-nt oligos, `GCTT`/`AGTG` vector overhangs, automatic pool packing, combinatorial primers, and BsaI are the defaults. Run this for all options:

```bash
python scripts/opool_cli.py --help
```

Existing outputs are protected by default. Choose a different `--run-name`/`--output-dir`, or explicitly pass `--force` to replace them.

### Codon-optimization species

For amino-acid inputs, set `CODON_SPECIES` in the simple or fast notebook, or pass `--codon-species` to the CLI. The default is `e_coli`. The built-in species keywords supplied by `python_codon_tables` are:

| Species | Short keyword | Full table name |
| --- | --- | --- |
| *Bacillus subtilis* | `b_subtilis` | `b_subtilis_1423` |
| *Caenorhabditis elegans* | `c_elegans` | `c_elegans_6239` |
| *Drosophila melanogaster* | `d_melanogaster` | `d_melanogaster_7227` |
| *Escherichia coli* | `e_coli` | `e_coli_316407` |
| *Gallus gallus* | `g_gallus` | `g_gallus_9031` |
| *Homo sapiens* | `h_sapiens` | `h_sapiens_9606` |
| *Mus musculus* | `m_musculus` | `m_musculus_10090` |
| *Mus musculus domesticus* | `m_musculus_domesticus` | `m_musculus_domesticus_10092` |
| *Saccharomyces cerevisiae* | `s_cerevisiae` | `s_cerevisiae_4932` |

Either the short keyword or full table name is accepted. Alternatively, use a numeric NCBI taxonomy ID; DnaChisel will then need internet access to retrieve its codon table. This setting is ignored when the input already contains optimized DNA.

## Notes

- Notebook defaults are set to the bundled example input files in `data/` and write outputs into `outputs/`.
- If you switch to your own datasets, update the single user-input cell in the simple notebook or the top configuration cell in the fast notebook.
