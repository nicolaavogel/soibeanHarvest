samtools faidx ${1}.fa

/projects/wintherpedersen/apps/gargammel/src/fragSim -n $2 --loc 3.7344 --scale 0.35 ${1}.fa > reads_${1}.fa
/projects/wintherpedersen/apps/gargammel/src/deamSim -matfile dhigh reads_${1}.fa > reads_${1}High.fa
/projects/wintherpedersen/apps/gargammel/src/adptSim -l 140 -artp reads_${1}HighAdapt.fa reads_${1}High.fa
/projects/wintherpedersen/apps/gargammel/art_src_MountRainier_Linux/art_illumina -ss HS25 -amp -na -p -l 140 -c 1 -i reads_${1}HighAdapt.fa -o ${1}High
/projects/wintherpedersen/apps/gargammel/leeHom --ancientdna -fq1 ${1}High1.fq -fq2 ${1}High2.fq -fqo ${1}High_fin_${2}


