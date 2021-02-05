# RNA_seq_metastases
**Objective**: Process and analyze RNA-seq (Ribo-) data from matched primary (breast) and metastatic (bone, adrenal, & brain) solid tumors. Identify differentially expressed transcripits in metastatic tumors, determine whether their expression can stratify patient outcome, and probe for druggable targets.

**Goal**: Identify therapeutic opportunities for prolonging and/or treating metastatic disease in ER+ breast cancer patients. 

Sequencing depth range of 19 - 28 million. Paired-end, short (76 bp) reads were assessed for quality using FASTQC and mapped to both human (GRCh38 + gencode v33 annotation) and mouse (GRCm38.p6 release 24 + gencode vM24 annotation) genomes using STAR aligner with relaxed parameters; --outFilterMatchNminOverLread 0.3 and –outFilterScoreMinOverLread 0.3. Computational deconvolution of mouse and human reads was done using XenofilteR. Human only filtered mapped reads were annotated using featureCounts with built-in SAF hg38 annotation file. Differential expression analysis was done using Bioconductor package edgeR. Counts were TMM normalized and classified significant at a threshold of foldchange>2. Functional enrichment of DEGs was done using Bioconductor packages enrichplot, clusterProfiler, and pathview and GO, KEGG, and MSigDB signature databases.


KEGG Pathway Maps:
       breast cancer = hsa05224
       PI3K-Akt signaling pathway = hsa04151
       Transcriptional misregulation in cancer = hsa05202
       MicroRNAs in cancer = hsa05206



Notes:
XenofilteR requires coordinate sorted bam file with NM-tag input file (o	--outSAMattributes NM)
https://github.com/PeeperLab/XenoFilteR/releases
XenofilteR limitation: Potential loss of conserved genes with reads of an equal edit distance to mouse and human genome are not assigned. Potential loss of conserved genes.

To analyze non-replicate data, coefficient of variation was extrapolated from biologically comparable data (3 reps of primary and brain) and applied to individual case analysis, generating dispersion estimates and fold change values to estimate differential expression between primary tumor and matched metastasis.



