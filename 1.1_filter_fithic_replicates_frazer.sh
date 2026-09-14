#!/bin/bash

#SBATCH --job-name=filter_fithic_frazer_replicate
#SBATCH --output=1.1_filter_fithic_frazer_replicate_%j.out
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=150g

# Script to generate SLURM scripts for filtering fithic replicates using Frazer replicate approach
# Creates example scripts in qshs directory

source ~/.bashrc

# Parameters
curr_date=$(date +"%y%m%d")
resolution=10000
fdr=0.0005

# Directories
baseDir="/mnt/BioAdHoc/Groups/vd-ay/bbabatunde/projects/25-06-Kuchroo-Ay"
workingDir="/home/bbabatunde/packages/25-11-frazer/Frazer_TB"
inputDir="${baseDir}/yard/250818_HiCPro/results/hicpro/hic_results/matrix"
pythonFile="${workingDir}/1.1_filter_fithic_replicate_frazer.py"
scriptsDir="${workingDir}/qshs/${curr_date}_filter_fithic_frazer_replicate_fdr${fdr}"
outputDir="${baseDir}/results/fithic_frazer_replicate_fdr${fdr}"

# Cell types to process
cellTypes=("npTh17" "pTh17-1" "Th0" "Th1" "Th2" "Treg")
# cellTypes=("Th2")

mkdir -p ${scriptsDir}
mkdir -p ${outputDir}

# Check if Python script exists
if [ ! -f "${pythonFile}" ]; then
    echo "Error: Python script not found: ${pythonFile}"
    exit 1
fi

echo "=================================================="
echo "Generating Frazer replicate filtering scripts"
echo "=================================================="
echo "Scripts directory: ${scriptsDir}"
echo "Output directory: ${outputDir}"
echo ""

for cellType in ${cellTypes[@]}; do
    echo "Processing: ${cellType}"
    
    # Collect replicate names and fithic files
    fithicFiles=()
    replicateNames=()
    
    while IFS= read -r folder; do
        replicate=$(basename ${folder})
        
        # Skip Th2-2 if present
        if [[ "${replicate}" == "Th2-2" ]]; then
            continue
        fi
        
        # Construct fithic file path
        fithicFolder="${folder}/fithic/${resolution}"
        fithicFile="${fithicFolder}/${replicate}.L20000.U3000000.p2.b200.spline_pass2.res${resolution}.significances.txt"
        
        # Check if file exists
        if [ ! -f "${fithicFile}" ]; then
            echo "  Warning: Fithic file not found: ${fithicFile}"
            continue
        fi
        
        fithicFiles+=("${fithicFile}")
        replicateNames+=("${replicate}")
        
        if [ -n "${verbose}" ]; then
            echo "  Found replicate: ${replicate}"
        fi
    done < <(find ${inputDir} -maxdepth 1 -type d -name "${cellType}*" | sort)
    
    if [ ${#fithicFiles[@]} -eq 0 ]; then
        echo "  Warning: No fithic files found for ${cellType}, skipping..."
        continue
    fi
    
    echo "  Found ${#fithicFiles[@]} replicates"
    
    # Create output directory for this cell type
    cellTypeOutputDir="${outputDir}/${cellType}-frazer-chrs"
    
    # Create scripts for each replicate and chromosome combination
    for i in "${!fithicFiles[@]}"; do
        fithicFile="${fithicFiles[$i]}"
        replicateName="${replicateNames[$i]}"
        
        for chrom in {1..19}; do
            scriptFile="${scriptsDir}/filter_frazer_${cellType}_${replicateName}_chr${chrom}.sh"
            
            cat <<EOF > ${scriptFile}
#!/bin/bash
#SBATCH --job-name=frazer-${cellType}-${replicateName}-chr${chrom}-fdr${fdr}
#SBATCH --output=${scriptsDir}/frazer-${cellType}-${replicateName}-chr${chrom}-fdr${fdr}_%j.out
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --mem=250g

source ~/.bashrc
mamba activate your_environment
cd ${workingDir}

python3 ${pythonFile} \\
    --input_file ${fithicFile} \\
    --replicate_name ${replicateName} \\
    --chromosome chr${chrom} \\
    --fdr_threshold ${fdr} \\
    --resolution ${resolution} \\
    --output_dir ${cellTypeOutputDir} \\
    --verbose

EOF
            chmod +x ${scriptFile}
        done
    done
    
    echo "  ✓ Created scripts for ${#fithicFiles[@]} replicates × 19 chromosomes = $(( ${#fithicFiles[@]} * 19 )) scripts"
    echo ""
done

echo "=================================================="
echo "Generated scripts in: ${scriptsDir}"
echo "Submit jobs using:"
echo "  sbatch ${scriptsDir}/filter_frazer_*.sh"
echo "=================================================="

