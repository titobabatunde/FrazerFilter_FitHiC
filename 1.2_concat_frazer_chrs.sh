#!/bin/bash

#SBATCH --job-name=concat_frazer_chrs
#SBATCH --output=1.2_concat_frazer_chrs_%j.out
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16g

# Concatenates the per-chromosome output of 1.1 into one genome-wide file per
# replicate, in both fithic and BEDPE format. Unlike 1.1 this does the work
# itself rather than generating jobs.
#
# ---------------------------------------------------------------------------
# HOW TO RUN   (run it from inside the repo -- workingDir defaults to $PWD)
#
#   sbatch 1.2_concat_frazer_chrs.sh
#   bash   1.2_concat_frazer_chrs.sh
#
# Prefer sbatch. Usually quick, but it reads every per-chromosome table for
# every replicate, which on a large run is minutes of I/O on a login node.
#
# Override any INPUT VARIABLE below on the command line:
#
#   CELL_TYPES="condA condB" FORCE=1 bash 1.2_concat_frazer_chrs.sh
#
#   sbatch --export=ALL,CELL_TYPES="condA condB",FORCE=1 1.2_concat_frazer_chrs.sh
# ---------------------------------------------------------------------------
#
# Input  (from 1.1): <results>/fithic_frazer_replicate_fdr<fdr>/<cellType>/
#                        frazer-chrs/<replicate>.<chrom>.frazer.fdr<fdr>.{txt,bedpe}
#
# Output:            <results>/fithic_frazer_replicate_fdr<fdr>/<cellType>/
#                        <replicate>.frazer.fdr<fdr>.{txt,bedpe}
#
# Replicates are discovered from what 1.1 actually produced, not by rescanning
# the fithic input, so a partial run concatenates what exists and says what is
# missing.

source ~/.bashrc

# ===========================================================================
# INPUT VARIABLES
# ===========================================================================
# THESE MUST MATCH WHAT 1.1 WAS RUN ON.

# Cell types / conditions to concatenate.
read -r -a cellTypes <<< "${CELL_TYPES:-npTh17 pTh17-1 Th0 Th1 Th2 Treg}"

# Chromosomes, in output order. Mouse autosomes by default.
read -r -a chroms <<< "${CHROMS:-chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19}"

# Output root, must match what 1.1 used.
resultsRoot="${RESULTS_DIR:-$(pwd)/results}"

# This repo. Defaults to the directory you run from.
workingDir="${WORKING_DIR:-$(pwd)}"

fdr="${FDR_THRESHOLD:-0.01}"

# Set FORCE=1 to overwrite existing concatenated files instead of skipping.
FORCE="${FORCE:-0}"

# ===========================================================================
# Derived - no need to edit below here
# ===========================================================================
outputDir="${resultsRoot}/fithic_frazer_replicate_fdr${fdr}"

if [ ! -d "${outputDir}" ]; then
    echo "Error: 1.1 output directory not found: ${outputDir}"
    echo "  Run 1.1_filter_fithic_replicates_frazer.sh first, or set RESULTS_DIR="
    exit 1
fi

echo "=================================================="
echo "Concatenating Frazer-filtered loops across chromosomes"
echo "=================================================="
echo "Directory:   ${outputDir}"
echo "Chromosomes: ${#chroms[@]} (${chroms[0]}..${chroms[$((${#chroms[@]}-1))]})"
echo ""

for cellType in ${cellTypes[@]}; do
    chrDir="${outputDir}/${cellType}/frazer-chrs"

    if [ ! -d "${chrDir}" ]; then
        echo "Warning: No 1.1 output for ${cellType}, skipping..."
        continue
    fi

    # Discover replicates from the filenames 1.1 wrote:
    #   <replicate>.<chrom>.frazer.fdr<fdr>.txt
    # Strip ".<chrom>.frazer.<everything>" to leave the replicate name.
    replicates=$(ls "${chrDir}"/*.frazer.fdr${fdr}.txt 2>/dev/null \
                 | xargs -n1 basename 2>/dev/null \
                 | sed 's/\.[^.]*\.frazer\..*$//' | sort -u)

    if [ -z "${replicates}" ]; then
        echo "Warning: No filtered files found for ${cellType}, skipping..."
        continue
    fi

    echo "Processing: ${cellType}"

    for replicate in ${replicates}; do
        outTxt="${outputDir}/${cellType}/${replicate}.frazer.fdr${fdr}.txt"
        outBedpe="${outputDir}/${cellType}/${replicate}.frazer.fdr${fdr}.bedpe"

        if [ -f "${outTxt}" ] && [ "${FORCE}" != "1" ]; then
            echo "  ${replicate}: already concatenated (FORCE=1 to overwrite)"
            continue
        fi

        # Write through temp files so an interrupted run cannot leave a
        # truncated table that the FORCE check above would then skip.
        tmpTxt="${outTxt}.tmp.$$"
        tmpBedpe="${outBedpe}.tmp.$$"
        : > "${tmpTxt}"; : > "${tmpBedpe}"

        nFound=0; nMissing=0; headerTxt=0; headerBedpe=0

        for chrom in ${chroms[@]}; do
            fTxt="${chrDir}/${replicate}.${chrom}.frazer.fdr${fdr}.txt"
            fBedpe="${chrDir}/${replicate}.${chrom}.frazer.fdr${fdr}.bedpe"

            if [ ! -f "${fTxt}" ]; then
                nMissing=$((nMissing + 1))
                continue
            fi

            # fithic format: header once, then rows
            if [ "${headerTxt}" -eq 0 ]; then
                head -1 "${fTxt}" > "${tmpTxt}"; headerTxt=1
            fi
            tail -n +2 "${fTxt}" >> "${tmpTxt}"

            # BEDPE: same, but the name column is loop_0..loop_N numbered per
            # chromosome, so it is renumbered genome-wide below.
            if [ -f "${fBedpe}" ]; then
                if [ "${headerBedpe}" -eq 0 ]; then
                    head -1 "${fBedpe}" > "${tmpBedpe}"; headerBedpe=1
                fi
                tail -n +2 "${fBedpe}" >> "${tmpBedpe}"
            fi

            nFound=$((nFound + 1))
        done

        if [ "${nFound}" -eq 0 ]; then
            echo "  ${replicate}: no chromosome files found, nothing written"
            rm -f "${tmpTxt}" "${tmpBedpe}"
            continue
        fi

        mv "${tmpTxt}" "${outTxt}"

        # Renumber the BEDPE name column so loop ids are unique genome-wide.
        if [ "${headerBedpe}" -eq 1 ]; then
            awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} {$7="loop_" (n++); print}' \
                "${tmpBedpe}" > "${outBedpe}"
            rm -f "${tmpBedpe}"
        else
            rm -f "${tmpBedpe}"
        fi

        nLoops=$(( $(wc -l < "${outTxt}") - 1 ))
        echo "  ${replicate}: ${nFound} chromosomes, ${nMissing} missing, ${nLoops} loops"
        echo "    ${outTxt}"
        [ -f "${outBedpe}" ] && echo "    ${outBedpe}"
    done
    echo ""
done

echo "=================================================="
echo "Concatenated files in: ${outputDir}/<cellType>/"
echo "=================================================="
