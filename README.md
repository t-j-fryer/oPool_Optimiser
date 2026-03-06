# oPool Cloning Notebook

Portable GitHub-ready version of `oPool_Cloning_Notebook.ipynb` with example inputs and reproducible Python/Jupyter kernel setup.

For AI assistant handoff, settings, workflow conventions, and FAQ, see `AI_WORKFLOW_SUMMARY.md`.

## Repository Layout

- `notebooks/oPool_Cloning_Notebook.ipynb`: main notebook (paths patched to be repo-relative)
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

## Running

```bash
jupyter lab
```

Open `notebooks/oPool_Cloning_Notebook.ipynb` and run cells in order.

## Notes

- Notebook defaults are set to the bundled example input files in `data/` and write outputs into `outputs/`.
- If you switch to your own datasets, update the top config variables in the notebook.
