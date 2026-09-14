#!/bin/bash
#SBATCH --job-name=filter_fithic_frazer_replicate
#SBATCH --output=1.1_filter_fithic_frazer_replicate_%j.out
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=150g

# Writes one SLURM job script per replicate per chromosome for the Frazer
# neighbour filter. This script generates jobs; it does not filter anything.
#
# ---------------------------------------------------------------------------
# HOW TO RUN   (workingDir defaults to $PWD)
#
#   sbatch 1.1_filter_fithic_replicates_frazer.sh
#   bash   1.1_filter_fithic_replicates_frazer.sh
#
# bash is fine for THIS script: it only writes job scripts and takes seconds.
# It is NOT fine for the jobs it generates -- those are per-chromosome filter
# runs and must be submitted with sbatch.
#
# Override any INPUT VARIABLE below on the command line:
#
#   CELL_TYPES="condA condB" INPUT_DIR=/path/to/matrix \
#     bash 1.1_filter_fithic_replicates_frazer.sh
#
#   sbatch --export=ALL,CELL_TYPES="condA condB",INPUT_DIR=/path/to/matrix \
#     1.1_filter_fithic_replicates_frazer.sh
# ---------------------------------------------------------------------------

source ~/.bashrc

# ===========================================================================
# INPUT VARIABLES
# ===========================================================================

# Cell types / conditions to process.
read -r -a cellTypes <<< "${CELL_TYPES:-npTh17 pTh17-1 Th0 Th1 Th2 Treg}"

# Replicates to skip, by exact name. 
read -r -a skipReplicates <<< "${SKIP_REPLICATES:-}"

# Chromosomes to write jobs for. Mouse autosomes by default; use chr1..chr22
# for human, and add chrX / chrY if you want them.
read -r -a chroms <<< "${CHROMS:-chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19}"

# Directory holding per-replicate fithic output.
inputDir="${INPUT_DIR:-/mnt/BioAdHoc/Groups/vd-ay/bbabatunde/projects/25-06-Kuchroo-Ay/yard/250818_HiCPro/results/hicpro/hic_results/matrix}"

# Filename of the fithic call table, after the leading "<replicate>.".
fithicTemplate="${FITHIC_TEMPLATE:-L20000.U3000000.p2.b200.spline_pass2.res${RESOLUTION:-10000}.significances.txt}"

# Where filtered output goes.
resultsRoot="${RESULTS_DIR:-$(pwd)/results}"

# This repo, holding the .py file. Defaults to the directory you run from.
workingDir="${WORKING_DIR:-$(pwd)}"

# Conda/mamba environment the generated jobs activate. Built by
# 0.0_create_frazerTB_env.sh.
envName="${ENV_NAME:-frazerTB}"

resolution="${RESOLUTION:-10000}"
fdr="${FDR_THRESHOLD:-0.0001}"

# Frazer filter parameters: an interaction is kept only if BOTH anchors have at
# least minNeighbors significant partners among the totalNeighbors bins
# flanking the opposing anchor.
minNeighbors="${MIN_NEIGHBORS:-3}"
totalNeighbors="${TOTAL_NEIGHBORS:-5}"

# ===========================================================================
# Derived - no need to edit below here
# ===========================================================================
curr_date=$(date +"%y%m%d")
pythonFile="${workingDir}/1.1_filter_fithic_replicate_frazer.py"
scriptsDir="${workingDir}/qshs/${curr_date}_filter_fithic_frazer_replicate_fdr${fdr}"
outputDir="${resultsRoot}/fithic_frazer_replicate_fdr${fdr}"

mkdir -p ${scriptsDir}
mkdir -p ${outputDir}

# Check if Python script exists
if [ ! -f "${pythonFile}" ]; then
    echo "Error: Python script not found: ${pythonFile}"
    echo "  Run this from inside the repo, or set WORKING_DIR=."
    exit 1
fi

if [ ! -d "${inputDir}" ]; then
    echo "Error: fithic input directory not found: ${inputDir}"
    echo "  Set INPUT_DIR= to where your per-replicate fithic output lives."
    exit 1
fi

echo "=================================================="
echo "Generating Frazer replicate filtering scripts"
echo "=================================================="
echo "Input directory:   ${inputDir}"
echo "Scripts directory: ${scriptsDir}"
echo "Output directory:  ${outputDir}"
echo "Filter:            >=${minNeighbors} of ${totalNeighbors} neighbours, both anchors, FDR<${fdr}"
echo ""

for cellType in ${cellTypes[@]}; do
    echo "Processing: ${cellType}"

    fithicFiles=()
    replicateNames=()

    while IFS= read -r folder; do
        replicate=$(basename ${folder})

        # Skip any replicate named in SKIP_REPLICATES
        skip=0
        for s in ${skipReplicates[@]}; do
            [[ "${replicate}" == "${s}" ]] && { skip=1; break; }
        done
        [ "${skip}" -eq 1 ] && { echo "  Skipping ${replicate} (in SKIP_REPLICATES)"; continue; }

        fithicFile="${folder}/fithic/${resolution}/${replicate}.${fithicTemplate}"

        if [ ! -f "${fithicFile}" ]; then
            if [ -f "${fithicFile}.gz" ]; then
                fithicFile="${fithicFile}.gz"
            else
                echo "  Warning: Fithic file not found: ${fithicFile}[.gz]"
                continue
            fi
        fi

        fithicFiles+=("${fithicFile}")
        replicateNames+=("${replicate}")
    done < <(find ${inputDir} -maxdepth 1 -type d -name "${cellType}*" | sort)

    if [ ${#fithicFiles[@]} -eq 0 ]; then
        echo "  Warning: No fithic files found for ${cellType}, skipping..."
        continue
    fi

    echo "  Found ${#fithicFiles[@]} replicates"

    # The python appends "/frazer-chrs" itself, so pass only the cell type here.
    # Upstream passed "${cellType}-frazer-chrs", producing a doubled
    # <cellType>-frazer-chrs/frazer-chrs/ path.
    cellTypeOutputDir="${outputDir}/${cellType}"

    for i in "${!fithicFiles[@]}"; do
        fithicFile="${fithicFiles[$i]}"
        replicateName="${replicateNames[$i]}"

        for chrom in ${chroms[@]}; do
            scriptFile="${scriptsDir}/filter_frazer_${cellType}_${replicateName}_${chrom}.sh"

            cat <<EOF > ${scriptFile}
#!/bin/bash
#SBATCH --job-name=frazer-${cellType}-${replicateName}-${chrom}-fdr${fdr}
#SBATCH --output=${scriptsDir}/frazer-${cellType}-${replicateName}-${chrom}-fdr${fdr}_%j.out
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --mem=250g

source ~/.bashrc
mamba activate ${envName}
cd ${workingDir}

python3 ${pythonFile} \\
    --input_file ${fithicFile} \\
    --replicate_name ${replicateName} \\
    --chromosome ${chrom} \\
    --fdr_threshold ${fdr} \\
    --resolution ${resolution} \\
    --min_neighbors ${minNeighbors} \\
    --total_neighbors ${totalNeighbors} \\
    --output_dir ${cellTypeOutputDir} \\
    --verbose

EOF
            chmod +x ${scriptFile}
        done
    done

    echo "  ✓ Created ${#fithicFiles[@]} replicates x ${#chroms[@]} chromosomes = $(( ${#fithicFiles[@]} * ${#chroms[@]} )) scripts"
    echo ""
done

echo "=================================================="
echo "Generated scripts in: ${scriptsDir}"
echo "Submit jobs using:"
echo "  sbatch ${scriptsDir}/filter_frazer_*.sh"
echo "=================================================="
