## RNA-seq data from matched primary (ER+ breast) and metastatic patient-derived, solid tumors. 
#### Objective: Perform differential expression analysis to identify transcripts uniquely upregulated in bone, brain, and adrenal metastases. Integrate public datasets to stratify patient outcome. Select clinically relevant, druggable targets for pre-clinical investigation in prolonging/preventing metastatic recurrence in ER+ breast cancer.

#
>
> FASTQC -> 1. STAR -> 2. XenofilteR -> 3. featurecounts -> **4. edgeR -> 5. Functional enrichment -> 6. cBioportal/GEO data integration** -> 7. Druggable target identification

#
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
6. Cox proportional hazard regression with mRNA expression of DEGs as predictors to stratify patient outcome.
7. Target identification was prioritized to repurposing FDA approved drugs. Dgidb.org was used to search drug-gene interactions + FDA status and PubChem and clinicaltrials.gov for specificity, toxicity, clinical history.
