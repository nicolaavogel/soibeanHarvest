#!/bin/bash

# Paths
TREE_DIR="/projects/wintherpedersen/people/bfj994/marineGraphs/tree_dir"
BELUGA_DIR="/projects/wintherpedersen/people/bfj994/marineGraphs/Beluga"
SEAL_DIR="/projects/wintherpedersen/people/bfj994/marineGraphs/harpSeal/"
NARWHAL_DIR="/projects/wintherpedersen/people/bfj994/marineGraphs/Narwhal"
BOW_DIR="/projects/wintherpedersen/people/bfj994/marineGraphs/bowHeadWhale"

# Read each FASTQ path from the file list
while IFS= read -r FASTQ; do
    BASENAME=$(basename "$FASTQ")
    OUTPUT_NAME=$(echo "$BASENAME" | sed 's/_fin_.*//')

    # Determine dbprefix and soibean_dir
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
        echo "⚠️  Skipping unknown species file: $BASENAME"
        continue
    fi

    # Echo the constructed command
    echo "/projects/wintherpedersen/apps/vgan/bin/vgan soibean \\"
    echo "  -fq1 \"$FASTQ\" \\"
    echo "  --dbprefix \"$DBPREFIX\" \\"
    echo "  --soibean_dir \"$SOIBEAN_DIR\" \\"
    echo "  --tree_dir \"$TREE_DIR\" \\"
    echo "  -o \"$OUTPUT_NAME\ -M \"$DBPREFIX\".rec"
    echo ""
done < fastq_filelist.txt

