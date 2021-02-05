#!/bin/bash

for SAMPLE in 2_S6 12_I_S7 12_S13 E_S5 5_I_S15

do

echo "performing bwa mem on" ${SAMPLE}

bwa mem hg38.fa -T 19 Griffin_RNAseq_1-15-20/${SAMPLE}_R1_001.fastq Griffin_RNAseq_1-15-20/${SAMPLE}_R2_001.fastq  > ${SAMPLE}-pe.sam

done