# FrazerFilter_FitHiC - Neighbour-Support Filtering of Fit-Hi-C Loops

This package filters Fit-Hi-C loop calls by FDR and by neighbour significance: a loop is kept only if both of its anchors are corroborated by significant interactions with the bins flanking the opposing anchor. The effect is to drop isolated single-pixel calls and keep loops whose signal is supported by their local neighbourhood.

Each replicate is filtered independently 

## Origin and Citation

The filtering criterion is taken from the methods of:

> W. W. Greenwald, N. Li, P. Benaglio, D. Jakubosky, H. Matsui, A. Schmitt, S. Selvaraj, M. D'Antonio, A. D'Antonio-Chronowska, E. N. Smith, K. A. Frazer, Subtle changes in chromatin loop contact propensity are associated with differential gene regulation and expression. *Nat. Commun.* **10**, 1054 (2019). doi: [10.1038/s41467-019-08940-5](https://doi.org/10.1038/s41467-019-08940-5)

That study quantifies chromatin loop *contact propensity*; this repo implements only the loop-filtering step from it, applied to **fithic** output. It does not compute contact propensity.

Please cite the paper above if you use this.

## Overview

1. Filter each replicate by FDR threshold
2. Keep only interactions whose anchors have enough significant neighbours
   (the Frazer filter, below)
3. Write the survivors in both fithic and BEDPE format, the latter for APAs
4. Concatenate the per-chromosome output into one genome-wide file per replicate

## Files

- **`0.0_create_frazerTB_env.sh`** - Builds the `frazerTB` mamba environment and verifies its imports
- **`frazerTB_env.yml`** - Conda/mamba environment spec
- **`1.1_filter_fithic_replicate_frazer.py`** - The filter, one replicate and chromosome per invocation
- **`1.1_filter_fithic_replicates_frazer.sh`** - Generates SLURM jobs per replicate per chromosome
- **`1.2_concat_frazer_chrs.sh`** - Concatenates the per-chromosome output into one file per replicate

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

"Three of the five immediately upstream or downstream bins" is the paper's
wording: three on *one* side is enough, and each side is checked separately, so
`total_neighbors` is the count per direction.

Two filters are applied before this one, inside
`load_and_sort_fithic_data`:

- `q-value < fdr_threshold`
- **anchors at least 32,000 bp apart** - hardcoded, not currently a parameter.
  It removes very short-range interactions; be aware it is stricter than the
  2 kb self-ligation cutoff used in the paper.

## Workflow

### Step 1: Filter Each Replicate

#### Generate SLURM Scripts

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
- `CHROMS`: Space-separated chromosomes to write jobs for (default: `chr1` through `chr19`). Assumes mouse genome.
- `SKIP_REPLICATES`: Replicate names to exclude (default: none)
- `FITHIC_TEMPLATE`: Fithic filename after the leading `{replicate}.` (default: `L20000.U3000000.p2.b200.spline_pass2.res{resolution}.significances.txt`)
- `RESULTS_DIR`: Output root (default: `$PWD/results`)
- `ENV_NAME`: Environment the generated jobs activate (default: `frazerTB`)
- `RESOLUTION`: Resolution in base pairs (default: 10000)
- `FDR_THRESHOLD`: q-value threshold for significance (default: 0.01, the paper's value)
- `MIN_NEIGHBORS` / `TOTAL_NEIGHBORS`: Filter stringency (default: 3 of 5)
- `WORKING_DIR`: Directory containing the `.py` file (default: `$PWD`)

#### Run the Python Script Directly

```bash
python3 1.1_filter_fithic_replicate_frazer.py \
    --input_file rep1.significances.txt \
    --replicate_name rep1 \
    --chromosome chr1 \
    --fdr_threshold 0.01 \
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
- `--fdr_threshold`: q-value threshold for significance (default: 0.01, the paper's value)
- `--min_neighbors`: Minimum number of significant neighbours required (default: 3)
- `--total_neighbors`: Total number of neighbours to check (default: 5)
- `--verbose`: Enable verbose output

### Step 2: Concatenate Across Chromosomes

Step 1 writes one file per replicate per chromosome. Step 2 joins them into a
single genome-wide file per replicate, in both formats.

```bash
bash 1.2_concat_frazer_chrs.sh
```

Unlike step 1 this does the work itself rather than generating jobs.

| | |
|---|---|
| Input | `<results>/fithic_frazer_replicate_fdr<fdr>/{condition}/frazer-chrs/{replicate}.{chrom}.frazer.fdr<fdr>.{txt,bedpe}` |
| Output | `<results>/fithic_frazer_replicate_fdr<fdr>/{condition}/{replicate}.frazer.fdr<fdr>.{txt,bedpe}` |

Replicates are discovered from what step 1 actually produced rather than by
rescanning the fithic input, so a partial run concatenates what exists and
reports what is missing. Check the `N chromosomes, M missing` line before using
the result.

An existing output file is left alone; pass `FORCE=1` to rebuild it.

**BEDPE loop ids are renumbered.** Step 1 names loops `loop_0`…`loop_N` within
each chromosome, so a naive concatenation would repeat ids. Step 2 renumbers the
`name` column genome-wide.

#### Arguments

**Required:** none beyond matching what step 1 was run on.

**Optional:**
- `CELL_TYPES`: Conditions to concatenate (default: the same list as step 1)
- `CHROMS`: Chromosomes, in output order (default: `chr1` through `chr19`)
- `RESULTS_DIR`: Output root; must match what step 1 used (default: `$PWD/results`)
- `FDR_THRESHOLD`: Must match step 1, since it is part of the filenames (default: 0.01)
- `FORCE`: Overwrite existing concatenated files (default: 0)
- `WORKING_DIR`: This repo (default: `$PWD`)

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

After step 2, the genome-wide files sit one level up, alongside `frazer-chrs/`:

- `…/{cellType}/{replicate_name}.frazer.fdr{fdr}.txt`
- `…/{cellType}/{replicate_name}.frazer.fdr{fdr}.bedpe`

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
cd qshs/{date}_filter_fithic_frazer_replicate_fdr0.01/
sbatch filter_frazer_{condition}_{replicate}_chr1.sh
sbatch filter_frazer_{condition}_{replicate}_chr2.sh
# ... etc, one per replicate x chromosome
# Wait for all of them to finish...

# 4. Concatenate the per-chromosome output
cd -
bash 1.2_concat_frazer_chrs.sh
```

Step 2 must not start until every chromosome job has finished. It will happily
concatenate a partial set and only warn about what is missing.
