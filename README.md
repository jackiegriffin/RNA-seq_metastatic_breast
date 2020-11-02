# RNA_seq_metastases
Ribodepleted RNA-seq was performed on total RNA extracted from primary and matched metastases in mice. Sequencing depth ranged from 1.9 to 2.8e7 and stranded (opposing strand) paired-end reads were 76 bp long. Read quality was assessed using FASTQC.
14 metastases and their matched primary breast tumors were analyzed for differential expression analysis to determine transcriptional changes driving metastatic progression:
1.	Paired-end reads were mapped to both the human genome GRCh38 with gencode v33 annotation and mouse genome GRCm38.p6 release 24 with gencode vM24 annotation, using the two-pass STAR1 alignment protocol1 v2.7.3a.
o	--outSAMattributes NM argument needed or XenoFilteR input
2.	Computational deconvolution of mouse and human reads was done using XenoFilteR
o	Download package: https://github.com/PeeperLab/XenoFilteR/releases
o	INPUT = coordinate sorted .bam files with ‘NM’-tag
3.	Number of mapped reads to annotated genes were counted using featurecounts.
o	featureCounts built-in SAF hg38 annotation file outputs entrez ids as meta features
4.	Differentially expression genes (DEGs) were quantified using the EdgeR package in R.
---------------------------------------------------------------------------------------------------------------------------
EDGER
Analyze non biological replicate data; primary and matched metastatic tumors. 
  - coefficient of variation extrapolated from similar data to get dispersion estimates and subsequent FC values to estimate DE expression between primary tumor and matched metastasis with single replicate data. 
  - Three biological replicates are available for estrogen deprived ZR75-1 primary and matched brain metastasis treatment group – use coefficient of variation from this analysis and apply for matched primary and metastatic tumors for other treatment groups where only 1 biological replicate is available (I.e., Estrogen deprived, MCF-7 primary and matched bone metastasis).

---------------------------------------------------------------------------------------------------------------------------

STAR
call HUMAN GENOME (hg38.fa): -- sjdbGTFfile human.gtf
call MOUSE GENOME: --sjdbGTFfile /mnt/c/Users/Jgriffin/Desktop/gencode.vM24.annotation.gtf
------------------------------------------------------------------------------------------------------------------------------------------
XenoFilteR
- Sequence reads that only map to a single reference genome are classified to that specific organism.
- For reads that map to both the human and mouse reference genome the edit distance (the number of base pairs different between the sequence read and the reference genome) is calculated by summing soft clips, insertions (CIGAR string) and the number of mismatches (‘NM’-tag) (F and R averaged for paired-end). Reads having a lower edit distance for the reference genome of a species are classified as originating from that species

Limitation: Potential loss of conserved genes with reads of an equal edit distance to mouse as well as human are not assigned. Poteltioal loss of conserved genes.


---------------------------------------------------------------------------------------------------------------------------
FEATURECOUNTS 

http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf :

featureCounts -d 20 -M -p -B -J -F SAF -a /mnt/c/Users/Jgriffin/Desktop/subread-2.0.0-source/annotation/hg38_RefSeq_exon.txt -O -o human_filtered_featCounts.txt _human_sorted_Filtered.bam  

- -d <int>; Minimum fragment/template length, 50 by default. This option must be used together with -p and -P
- -M; If specified, multi-mapping reads/fragments will be counted. The program uses the ‘NH’ tag to find multi-mapping reads. Each alignment reported for a multi-mapping read will be counted individually. Each alignment will carry 1 count or a fractional count (--fraction).
- -p; (isPairedEnd) If specified, fragments (or templates) will be counted instead of reads. This option is only applicable for paired-end reads.
- -B; If specified, only fragments that have both ends successfully aligned will be considered for summarization
- -J; Count the number of reads supporting each exon-exon junction. Junctions are identified from those exon-spanning reads (containing ‘N’ in CIGAR string) in input data (note that options ‘–splitOnly’ and ‘–nonSplitOnly’ are not considered by this parameter). The output result includes names of primary and secondary genes that overlap at least one of the two splice sites of a junction. Only one primary gene is reported, but there might be more than one secondary gene reported. Secondary genes do not overlap more splice sites than the primary gene. When the primary and secondary genes overlap same number of splice sites, the gene with the smallest leftmost base position is selected as the primary gene. Also included in the output result are the position information for the left splice site (‘Site1’) and the right splice site (‘Site2’) of a junction. These include chromosome name, coordinate and strand of the splice site. In the last columns of the output, number of supporting reads is provided for each junction for each library.
- -F (isGTFAnnotationFile); Specify the format of the annotation file. Acceptable formats include ‘GTF’ and ‘SAF’ (see Section 6.2.2 for details). By default, C version of featureCounts program accepts a GTF format annotation and R version accepts a SAF format annotation. In-built annotations in SAF format are provided.
- -O; (allowMultiOverlap) If specified, reads (or fragments if -p is specified) will be allowed to be assigned to more than one matched meta-feature (or feature if -f is specified). Reads/fragments overlapping with more than one meta-feature/feature will be counted more than once. Note that when performing meta-feature level summarization, a read (or fragment) will still be counted once if it overlaps with multiple features within the same meta-feature (as long as it does not overlap with other metafeatures). Also note that this parameter is applied to each individual alignment when there are more than one alignment reported for a read (ie. multi-mapping read).


---------------------------------------------------------------------------------------------------------------------------



1STAR: ultrafast universal RNA-seq aligner, dobin 2013

---------------------------------------------------------------------------------------------------------------------------

