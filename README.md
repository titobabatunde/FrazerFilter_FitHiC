# FrazerFilter_FitHiC - Neighbour-Support Filtering of Fit-Hi-C Loops

This package filters Fit-Hi-C loop calls by FDR and by neighbour significance: a
loop is kept only if both of its anchors are corroborated by significant
interactions with the bins flanking the opposing anchor. The effect is to drop
isolated single-pixel calls and keep loops whose signal is supported by their
local neighbourhood.

Each replicate is filtered independently - there are no cross-replicate
comparisons and no union files.

## Origin and Citation

The filtering criterion is taken from the methods of:

> W. W. Greenwald, N. Li, P. Benaglio, D. Jakubosky, H. Matsui, A. Schmitt, S. Selvaraj, M. D'Antonio, A. D'Antonio-Chronowska, E. N. Smith, K. A. Frazer, Subtle changes in chromatin loop contact propensity are associated with differential gene regulation and expression. *Nat. Commun.* **10**, 1054 (2019). doi: [10.1038/s41467-019-08940-5](https://doi.org/10.1038/s41467-019-08940-5)

That study quantifies chromatin loop *contact propensity*; this repo implements
only the loop-filtering step from it, applied to **fithic** output. It does not
compute contact propensity.

Please cite the paper above if you use this.

## Overview

1. Filter each replicate by FDR threshold
2. Keep only interactions whose anchors have enough significant neighbours
   (the Frazer filter, below)
3. Write the survivors in both fithic and BEDPE format, the latter for APAs

## Files

- **`0.0_create_frazerTB_env.sh`** - Builds the `frazerTB` mamba environment and verifies its imports
- **`frazerTB_env.yml`** - Conda/mamba environment spec
- **`1.1_filter_fithic_replicate_frazer.py`** - The filter, one replicate and chromosome per invocation
- **`1.1_filter_fithic_replicates_frazer.sh`** - Generates SLURM jobs per replicate per chromosome

## Setup

Requires `mamba` (or `conda`). Everything installs from `conda-forge`; the
filter is pure pandas/numpy, so there are no version-sensitive pins.

```bash
bash 0.0_create_frazerTB_env.sh
mamba activate frazerTB
```

The script verifies that `pandas` and `numpy` import afterwards - a successful
solve alone is not proof the environment works.

Options:

```bash
ENV_NAME=frazerTB2 bash 0.0_create_frazerTB_env.sh              # build side-by-side
ENV_PREFIX=/path/to/envs/frazerTB bash 0.0_create_frazerTB_env.sh
```

## The Frazer Filter

For an interaction A ↔ B to pass, with defaults `min_neighbors=3` of
`total_neighbors=5`:

1. **Anchor A test** - A must have ≥3 significant interactions among the bins
   flanking B: A ↔ B-1, A ↔ B-2, … A ↔ B+1, A ↔ B+2, …
2. **Anchor B test** - B must have ≥3 significant interactions among the bins
   flanking A: A-1 ↔ B, A-2 ↔ B, … A+1 ↔ B, A+2 ↔ B, …
3. **Both must pass.** Either anchor failing drops the interaction.

## Workflow

### Generate SLURM Scripts

```bash
bash 1.1_filter_fithic_replicates_frazer.sh
```

This will:
- Find every replicate directory matching each cell type under `$INPUT_DIR`
- Locate each replicate's fithic significances file (`.txt` or `.txt.gz`)
- Generate SLURM scripts in `qshs/{date}_filter_fithic_frazer_replicate_fdr{fdr}/`
- Create one script per replicate per chromosome

Then submit:

```bash
sbatch qshs/*/filter_frazer_*.sh
```

#### Arguments

Set as environment variables. The defaults point at the original project, so the
required ones must be changed for your own data.

**Required:**
- `CELL_TYPES`: Space-separated list of cell types / conditions to process
- `INPUT_DIR`: Directory of per-replicate fithic output, as `<dir>/{replicate}/fithic/{resolution}/{replicate}.{FITHIC_TEMPLATE}`

**Optional:**
- `CHROMS`: Space-separated chromosomes to write jobs for (default: `chr1` through `chr19`)
- `SKIP_REPLICATES`: Replicate names to exclude (default: `Th2-2`; set to empty to skip none)
- `FITHIC_TEMPLATE`: Fithic filename after the leading `{replicate}.` (default: `L20000.U3000000.p2.b200.spline_pass2.res{resolution}.significances.txt`)
- `RESULTS_DIR`: Output root (default: `$PWD/results`)
- `ENV_NAME`: Environment the generated jobs activate (default: `frazerTB`)
- `RESOLUTION`: Resolution in base pairs (default: 10000)
- `FDR_THRESHOLD`: FDR threshold for significance (default: 0.0005)
- `MIN_NEIGHBORS` / `TOTAL_NEIGHBORS`: Filter stringency (default: 3 of 5)
- `WORKING_DIR`: Directory containing the `.py` file (default: `$PWD`)

### Run the Python Script Directly

```bash
python3 1.1_filter_fithic_replicate_frazer.py \
    --input_file rep1.significances.txt \
    --replicate_name rep1 \
    --chromosome chr1 \
    --fdr_threshold 0.0005 \
    --resolution 10000 \
    --output_dir /path/to/output \
    --verbose
```

#### Arguments

**Required:**
- `--input_file`: Input fithic significance file (`.txt` or `.txt.gz`) for a single replicate
- `--replicate_name`: Replicate name (e.g., `Th2-1`)
- `--resolution`: Resolution in base pairs (e.g., 10000)
- `--output_dir`: Output directory for results

**Optional:**
- `--chromosome`: Chromosome to process (default: `chr1`)
- `--fdr_threshold`: FDR threshold for significance (default: 0.0001)
- `--min_neighbors`: Minimum number of significant neighbours required (default: 3)
- `--total_neighbors`: Total number of neighbours to check (default: 5)
- `--verbose`: Enable verbose output

Note the two FDR defaults differ: the Python defaults to `0.0001`, the generator
passes `0.0005`. The generator's value wins for SLURM runs.

## Running the scripts

Run from inside the repo - `workingDir` defaults to `$PWD`.

```bash
sbatch 1.1_filter_fithic_replicates_frazer.sh     # recommended
bash   1.1_filter_fithic_replicates_frazer.sh     # fine: only writes job scripts
```

`bash` is fine for the generator, which takes seconds. The jobs it writes are
the heavy part and must go through `sbatch`.

Override any input variable without editing the file:

```bash
CELL_TYPES="condA condB" INPUT_DIR=/path/to/matrix \
  bash 1.1_filter_fithic_replicates_frazer.sh

sbatch --export=ALL,CELL_TYPES="condA condB",INPUT_DIR=/path/to/matrix \
  1.1_filter_fithic_replicates_frazer.sh
```

## Output

For each replicate and chromosome, under `{output_dir}/frazer-chrs/`:

- `{replicate_name}.{chromosome}.frazer.fdr{fdr}.txt` - survivors, fithic format
- `{replicate_name}.{chromosome}.frazer.fdr{fdr}.bedpe` - same, BEDPE format for APAs

With the generator, `{output_dir}` is
`{RESULTS_DIR}/fithic_frazer_replicate_fdr{fdr}/{cellType}`, so the files land in
`…/{cellType}/frazer-chrs/`.

## Example Workflow

`{condition}` below is one entry of `CELL_TYPES`, `{replicate}` one of its
replicate directories, and `{date}` the `YYMMDD` stamp the generator applies.

```bash
# 1. Set up environment
bash 0.0_create_frazerTB_env.sh
mamba activate frazerTB

# 2. Generate SLURM scripts
bash 1.1_filter_fithic_replicates_frazer.sh

# 3. Submit jobs
cd qshs/{date}_filter_fithic_frazer_replicate_fdr0.0005/
sbatch filter_frazer_{condition}_{replicate}_chr1.sh
sbatch filter_frazer_{condition}_{replicate}_chr2.sh
# ... etc, one per replicate x chromosome
```
