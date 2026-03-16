for name in $(cat fqlist); do

    species=$(echo "$name" | cut -d'/' -f1)
    echo "$species"
    nice -19 /home/projects/animalia_mito/vgan/bin/vgan soibean -fq1 ${name}.fq.gz -o ${name} --soibean_dir /home/projects2/animalia_mito_external/costumeGraphs/${species}/ --tree_dir /home/projects2/animalia_mito_external/costumeGraphs/tree_dir/ --dbprefix ${species} -t 20 -k 3 2> ${name}.log

done
