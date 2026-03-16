/home/ctools/vg_1.44.0/bin/vg paths -Q $1 -F -x $2 > ${3}.fa

samtools faidx ${3}.fa

/home/ctools/gargammel/src/fragSim -n $4 -f size_freq.tsv  --circ $1 ${3}.fa > reads_${3}.fa
/home/ctools/gargammel/src/deamSim -mapdamage misincorporation.txt double reads_${3}.fa > reads_${3}High.fa
/home/ctools/gargammel/src/adptSim -l 140 -artp reads_${3}HighAdapt.fa reads_${3}High.fa
/home/ctools/gargammel/art_src_MountRainier_Linux/art_illumina -ss HS25 -amp -na -p -l 140 -c 1 -i reads_${3}HighAdapt.fa -o ${3}High
/home/ctools/leehom-1.2.17 --ancientdna -fq1 ${3}High1.fq -fq2 ${3}High2.fq -fqo ${3}High_fin_${4}

