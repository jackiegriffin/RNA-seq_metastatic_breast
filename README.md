## RNA-seq data from matched primary (ER+ breast) and metastatic, patient-derived, solid tumors. 
#### Objective: Using edgeR, perform differential expression analysis and identify uniquely upregulated transcripts driving metastases to the bone, brain, and adrenals.  



> FASTQC -> 1. STAR -> 2. XenofilteR -> 3. featurecounts -> **4. edgeR** 


1. Paired-end, short (76 bp) reads were mapped to:
   + Human (GRCh38 + gencode v33 annotation)
   + Mouse (GRCm38.p6 release 24 + gencode vM24 annotation) 
2. Computational deconvolution of mouse and human reads was done using XenofilteR.
   + limitation: Potential loss of conserved genes; equal edit distance genes are not assigned.
3. Human only filtered reads were annotated using featureCounts (SAF hg38).
4. Counts were normalized using the weighted trimmed mean of log expression ratios (TMM) and significance threshold of foldchange 2 was set.
5. Functional enrichment of DEGs was done using Bioconductor packages enrichplot, clusterProfiler, and pathview and GO, KEGG databases:
   + breast cancer = hsa05224
   + Transcriptional misregulation in cancer = hsa05202
   + MicroRNAs in cancer = hsa05206
   + MSigDB oncogenic C6
   
   
   
   
   
   

