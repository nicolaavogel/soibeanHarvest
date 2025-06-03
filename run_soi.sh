#!/bin/bash
#SBATCH --job-name=soibean_runs
#SBATCH --partition=compregular
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
#SBATCH --array=0-47%10
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=10G

# Define directories
TREE_DIR="/projects/wintherpedersen/people/bfj994/marineGraphs/tree_dir"
BELUGA_DIR="/projects/wintherpedersen/people/bfj994/marineGraphs/Beluga"
SEAL_DIR="/projects/wintherpedersen/people/bfj994/marineGraphs/harpSeal"
NARWHAL_DIR="/projects/wintherpedersen/people/bfj994/marineGraphs/Narwhal"
BOW_DIR="/projects/wintherpedersen/people/bfj994/marineGraphs/bowHeadWhale"

# Get FASTQ from filelist based on SLURM array index (1-based)
FASTQ=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" fastq_filelist.txt)

# Safety check
if [[ -z "$FASTQ" ]]; then
    echo "❌ Error: No FASTQ file found for index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

BASENAME=$(basename "$FASTQ")
OUTPUT_NAME=$(echo "$BASENAME" | sed 's/_fin_.*//')

# Determine prefix and input directory
if [[ "$BASENAME" == harpSeal* ]]; then
    DBPREFIX="harpSeal"
    SOIBEAN_DIR="$SEAL_DIR"
elif [[ "$BASENAME" == Beluga* ]]; then
    DBPREFIX="Beluga"
    SOIBEAN_DIR="$BELUGA_DIR"
elif [[ "$BASENAME" == Narwhal* ]]; then
    DBPREFIX="Narwhal"
    SOIBEAN_DIR="$NARWHAL_DIR"
elif [[ "$BASENAME" == bowHeadWhale* ]]; then
    DBPREFIX="bowHeadWhale"
    SOIBEAN_DIR="$BOW_DIR"
else
    echo "❌ Error: Unknown species prefix in filename $BASENAME"
    exit 1
fi

# Run vgan
/projects/wintherpedersen/apps/vgan/bin/vgan soibean \
    -fq1 "$FASTQ" \
    --dbprefix "$DBPREFIX" \
    --soibean_dir "$SOIBEAN_DIR" \
    --tree_dir "$TREE_DIR" \
    -o "${SOIBEAN_DIR}/${OUTPUT_NAME}" \
    -t 20 \
    -M "${DBPREFIX}.rec" \
	-k 3
