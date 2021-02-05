# Libraries ----
library("pathview")

#   KEGG Pathway Maps ----
#       breast cancer = hsa05224
#       PI3K-Akt signaling pathway = hsa04151
#       Transcriptional misregulation in cancer = hsa05202
#       MicroRNAs in cancer = hsa05206
#      

##    COMPLETE --------------------------------------------------------
  complete <- read.csv("Pathway_Enrichment_Analysis/Data_files_INPUT/gene.list_sig_complete_DEGs.csv")
  geneListcomplete <- complete[,2]
  names(geneListcomplete) <- as.character(complete[,1])
  geneListcomplete <- sort(geneListcomplete, decreasing = TRUE)
  gene <- names(geneListcomplete)[abs(geneListcomplete) > 2] 
  hsa03100 <- pathview(gene.data  = geneListcomplete,
                       pathway.id = "hsa03100",
                       species    = "hsa",
                       limit      = list(gene=max(abs(geneListcomplete)), cpd=1),
                       out.suffix = "complete")

##    BONE ONLY --------------------------------------------------------
  bone <- read.csv("Pathway_Enrichment_Analysis/Data_files_INPUT/bone_gene_list.csv")
  geneListbone <- bone[,2]
  names(geneListbone) <- as.character(bone[,1])
  geneListbone <- sort(geneListbone, decreasing = TRUE)
  gene <- names(geneListbone)[abs(geneListbone) > 2] 
  hsa03100 <- pathview(gene.data  = geneListbone,
                              pathway.id = "hsa03100",
                              species    = "hsa",
                              limit      = list(gene=max(abs(geneListbone)), cpd=1),
                              out.suffix = "bone")
  
  ##    BRAIN ONLY --------------------------------------------------------
  brain <- read.csv("Pathway_Enrichment_Analysis/Data_files_INPUT/brain_gene_list.csv")
  geneListbrain <- brain[,2]
  names(geneListbrain) <- as.character(brain[,1])
  geneListbrain <- sort(geneListbrain, decreasing = TRUE)
  gene <- names(geneListbrain)[abs(geneListbrain) > 2] 
  hsa04512 <- pathview(gene.data  = geneListbrain,
                            pathway.id = "hsa04512",
                            species    = "hsa",
                            limit      = list(gene=max(abs(geneListbrain)), cpd=1),
                            out.suffix = "brain")
  
  ##    ADRENAL ONLY --------------------------------------------------------
  adrenal <- read.csv("Pathway_Enrichment_Analysis/Data_files_INPUT/adrenal_gene_list.csv")
  geneListadrenal <- adrenal[,2]
  names(geneListadrenal) <- as.character(adrenal[,1])
  geneListadrenal <- sort(geneListadrenal, decreasing = TRUE)
  gene <- names(geneListadrenal)[abs(geneListadrenal) > 2] 
  hsa05206 <- pathview(gene.data  = geneListadrenal,
                       pathway.id = "hsa05206",
                       species    = "hsa",
                       limit      = list(gene=max(abs(geneListadrenal)), cpd=1),
                       out.suffix = "adrenal",
                       kegg.native = FALSE)
  
  
  
  
  