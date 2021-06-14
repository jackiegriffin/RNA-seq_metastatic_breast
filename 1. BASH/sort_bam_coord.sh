#!/bin/bash

for SAMPLE in 1P_S20 2_II_S16 2_S6 3P_S21 4P_S17 5_I_S15 5P_S19 6_S10 8_I_S9 10_S8 11_II_S12 12_I_S7 12_S13 15_I_S14 16_S2 22_S1 25_S11 26_S4 C_S3 E_S5 2P_S18 

do

echo "performing samtools on" ${SAMPLE}


samtools sort /mnt/e/human_aligned_star_3/${SAMPLE}_humanAligned.sortedByCoord.out.bam > /mnt/e/human_aligned_star_3/samtools_sorted_human_3/${SAMPLE}_human.sorted.bam

done


