# RNA_seq_metastases
Objective: Process and analyze RNA-seq (Ribo-) data from matched primary (breast) and metastatic (bone, adrenal, & brain) solid tumors. Identify therapeutic opportunities to prolong and/or treat metastatic disease in ER+ breast cancer patients. 

Sequencing depth range of 19 - 28 million. Paired-end, short (76 bp) reads were assessed for quality using FASTQC and mapped to both human (GRCh38 + gencode v33 annotation) and mouse (GRCm38.p6 release 24 + gencode vM24 annotation) genomes using STAR aligner with relaxed parameters; --outFilterMatchNminOverLread 0.3 and –outFilterScoreMinOverLread 0.3. Computational deconvolution of mouse and human reads was done using XenofilteR. Human only filtered mapped reads were annotated using featureCounts with built-in SAF hg38 annotation file. Differential expression analysis was done using Bioconductor package edgeR. 


XenofilteR requires coordinate sorted bam file with NM-tag input file (o	--outSAMattributes NM)
https://github.com/PeeperLab/XenoFilteR/releases
XenofilteR limitation: Potential loss of conserved genes with reads of an equal edit distance to mouse and human genome are not assigned. Potential loss of conserved genes.

To analyze differential expression (edgeR) using non biological replicate data, a coefficient of variation was extrapolated from analysis of biologically comparable data (3 reps of primary and brain) and applied to individual case studies, generating dispersion estimates and fold change values to estimate differential expression between primary tumor and matched metastasis.





featureCounts -d 20 -M -p -B -J -F SAF -a /mnt/c/Users/Jgriffin/Desktop/subread-2.0.0-source/annotation/hg38_RefSeq_exon.txt -O -o human_filtered_featCounts.txt _human_sorted_Filtered.bam  
