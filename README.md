# FrazerTB - Fit-Hi-C Individual Replicate Filtering with Frazer Approach

This package implements the Frazer filtering approach for filtering individual Fit-Hi-C loop replicate files. It filters interactions by FDR threshold and neighbor significance criteria for each replicate independently.

## Overview

The Frazer filtering approach:
1. Filters each replicate file by FDR threshold
2. Filters each replicate by requiring significant interactions with neighbors around opposing anchors (Frazer filtering)
3. Outputs filtered files for each replicate independently in both fithic and BEDPE formats for APAs

Each replicate is processed independently - no cross-replicate comparisons or union files are created.

## Files

- **`1.1_filter_fithic_replicate_frazer.py`** - Main Python script for filtering individual fithic replicates
- **`1.1_filter_fithic_replicates_frazer.sh`** - Bash script to generate SLURM jobs per replicate and chromosome
- **`create_frazerTB_env.sh`** - Script to create the `frazerTB` mamba environment
- **`frazerTB_env.yml`** - Conda/mamba environment file with all dependencies

## Setup

### 1. Create the Environment

First, create the `frazerTB` mamba environment:

```bash
bash create_frazerTB_env.sh
```

This will:
- Create a new mamba environment called `frazerTB` with required packages (pandas, numpy)


### 2. Activate the Environment

```bash
mamba activate frazerTB
```

## Usage

### Option 1: Generate SLURM Scripts (Recommended)

The bash script automatically generates SLURM scripts for each replicate and chromosome combination:

```bash
bash 1.1_filter_fithic_replicates_frazer.sh
```

This will:
- Find all replicates for each cell type
- Generate SLURM scripts in `qshs/{date}_filter_fithic_frazer_replicate_fdr{fdr}/`
- Create scripts for each replicate × chromosome combination (e.g., 3 replicates × 19 chromosomes = 57 scripts per cell type)

Then submit the jobs:

```bash
sbatch qshs/*/filter_frazer_*.sh
```

### Option 2: Run Python Script Directly

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

## Python Script Arguments

### Required Arguments

- `--input_file`: Input fithic significance file (.txt or .txt.gz) for a single replicate
- `--replicate_name`: Replicate name (e.g., Th2-1)
- `--resolution`: Resolution in base pairs (e.g., 10000)
- `--output_dir`: Output directory for results

### Optional Arguments

- `--chromosome`: Chromosome to process (default: chr1)
- `--fdr_threshold`: FDR threshold for significance (default: 0.0001)
- `--min_neighbors`: Minimum number of significant neighbors required (default: 3)
- `--total_neighbors`: Total number of neighbors to check (default: 5)
- `--verbose`: Enable verbose output

## Frazer Filtering Algorithm

For an interaction A ↔ B to pass the Frazer filter:

1. **Anchor A test**: Anchor A must have ≥`min_neighbors` significant interactions with upstream/downstream bins around anchor B
   - Checks: A ↔ B-1, A ↔ B-2, ..., A ↔ B+1, A ↔ B+2, ...
   
2. **Anchor B test**: Anchor B must have ≥`min_neighbors` significant interactions with upstream/downstream bins around anchor A
   - Checks: A-1 ↔ B, A-2 ↔ B, ..., A+1 ↔ B, A+2 ↔ B, ...

3. **Both anchors must pass**: The interaction is kept only if BOTH anchors pass their respective tests.

## Output Files

For each replicate, the script generates:

1. **Fithic format file**:
   - `{replicate_name}.{chromosome}.frazer.fdr{fdr}.txt` - Filtered significant interactions that passed FDR and neighbor significance tests

2. **BEDPE format file**:
   - `{replicate_name}.{chromosome}.frazer.fdr{fdr}.bedpe` - Same interactions in BEDPE format for visualization/APAs

All files are saved in: `{output_dir}/frazer-chrs/`

## Example Workflow

```bash
# 1. Set up environment
bash create_frazerTB_env.sh
mamba activate frazerTB

# 2. Generate SLURM scripts
bash 1.1_filter_fithic_replicates_frazer.sh

# 3. Submit jobs
cd qshs/250101_filter_fithic_frazer_replicate_fdr0.0005/
sbatch filter_frazer_npTh17_Th2-1_chr1.sh
sbatch filter_frazer_npTh17_Th2-1_chr2.sh
sbatch filter_frazer_npTh17_Th2-3_chr1.sh
# ... etc (one script per replicate × chromosome combination)
```
