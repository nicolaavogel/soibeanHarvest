#!/bin/bash

# Check if both species and fqlist arguments are provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <species> <fqlist>"
    exit 1
fi

species=$1        # Set the species from the first command-line argument
fqlist=$2         # Set the fqlist file from the second command-line argument

# Loop through each entry in the specified fqlist file
for name in $(cat "$fqlist"); do
    # Extract the name between the first '/' and '.fq'
    base_name=$(echo "$name" | sed 's|.*/||; s|\.fq.*||')
    echo "Processing file: $base_name for species: $species"

    # Run the vgan command with the extracted name and specified species
    nice -19 /home/projects/animalia_mito/vgan/bin/vgan soibean -fq1 ${name} -o ${base_name} --soibean_dir /home/projects2/animalia_mito_external/costumeGraphs/${species}/ --tree_dir /home/projects2/animalia_mito_external/costumeGraphs/tree_dir --dbprefix ${species} --iter 1000000 --burnin 150000 -k 4 -t 5 2> ${base_name}.log

done
