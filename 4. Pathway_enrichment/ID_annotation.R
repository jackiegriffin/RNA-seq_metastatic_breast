library(clusterProfiler)
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)


d <- read.csv("4. Pathway_enrichment/DEGs_organ_stratified/bone.2.11.2021.csv")

## Format data ----
geneList <- d[,2]
names(geneList) <- as.character(d[,1])## feature 2: named vector
geneList <- sort(geneList, decreasing = TRUE)## feature 3: decreasing order
gene <- names(geneList)[abs(geneList) > 2] # set genelist to include those with FC > 2
