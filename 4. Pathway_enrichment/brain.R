# libraries ----
library(clusterProfiler)
library(magrittr)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)

# brain data ----
d <- read.csv("4. Pathway_enrichment/DEGs_organ_stratified/brain.2.11.2021.csv")

# Format data ----
geneList <- d[,2]
names(geneList) <- as.character(d[,1])## feature 2: named vector
geneList <- sort(geneList, decreasing = TRUE)## feature 3: decreasing order

gene <- names(geneList)[abs(geneList) > 1]
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)

ego3_brain <- gseGO(geneList     = geneList,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "CC",
                   nPerm        = 1000,
                   minGSSize    = 100,
                   maxGSSize    = 500,
                   # pvalueCutoff = 0.05,
                   verbose      = FALSE)
head(ego3_brain)


ego3_brain_dot <- dotplot(ego3_brain, showCategory=20) + ggtitle("gseGO (Brain)")
plot(ego3_brain_dot)

###-------



enrich_KEGG <- enrichKEGG(gene = gene, organism='hsa',pvalueCutoff = 0.05)
KEGG_brain <- dotplot(enrich_KEGG, showCategory=20) + ggtitle("KEGG ORA for brain metastases")
plot(KEGG_brain)


plot_grid(KEGG_bone, p2, ncol=2)
# GO ORA ---
ego_brain <- enrichGO(gene          = gene,
                     universe      = names(geneList),
                     OrgDb         = org.Hs.eg.db,
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
head(ego_brain)

GO_ORA_bar_brain <- barplot(ego_brain, showCategory=40)
plot(GO_ORA_bar_brain)
GO_ORA_cnet_brain <- cnetplot(ego_brain, foldChange=geneList)
plot(GO_ORA_cnet_brain)



# GO Gene set enrighment ----
ggo_brain <- groupGO(gene     = gene,
                    OrgDb    = org.Hs.eg.db,
                    ont      = "CC",
                    level    = 3,
                    readable = TRUE)

head(ggo_brain)

go_group_brain <-barplot(ggo_brain, showCategory=40)
plot(go_group_brain)




