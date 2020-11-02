#!/bin/bash

for SAMPLE in 1P_S20 2_II_S16 2_S6 3P_S21 4P_S17 5_I_S15 5P_S19 6_S10 8_I_S9 10_S8 11_II_S12 12_I_S7 12_S13 15_I_S14 16_S2 22_S1 25_S11 26_S4 C_S3 E_S5 2P_S18 

do

echo "performing STAR on" ${SAMPLE}


STAR --runThreadN 4 --genomeDir /mnt/c/Users/Jgriffin/Desktop/star_genome_idx/ --chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimOutJunctionFormat 1 --twopassMode Basic --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --alignSJDBoverhangMin 3 --chimMultimapNmax 20 --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 5 --alignSplicedMateMapLminOverLmate 0 --sjdbOverhang 75 --chimScoreDropMax 37 --quantMode GeneCounts --sjdbGTFfile /mnt/c/Users/Jgriffin/Desktop/gencode.vH34.annotation.gtf --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --readFilesIn /mnt/e/fastq/${SAMPLE}_R1_001.fastq /mnt/e/fastq/${SAMPLE}_R2_001.fastq --outFileNamePrefix /mnt/e/human_aligned_star_3/${SAMPLE}_human

done
