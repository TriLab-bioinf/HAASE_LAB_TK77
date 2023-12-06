# load library
library(org.Mm.eg.db)
library(clusterProfiler)
library(dplyr)

# import data
data <- read.delim("E16.5_vs_oxistress_E16.5.bed",sep="\t",header=T)
data2 <- read.delim("E16.5_vs_oxistress_oxistress.bed",sep="\t",header=T)
data3 <- read.delim("E16.5_vs_oxistress_oxistress.bed_E16.5.bed",sep="\t",header=T)

res <- rbind(data,data2,data3)

res$Parent.files <- gsub(".bed","",res$Parent.files)

res$gene <- case_when(
    is.na(res[,10]) ~ res[,9],
    is.na(res[,9]) ~ res[,10],
    .default = paste0(res[,9],",",res[,10])
         )

comparisons <- unique(res$Parent.files)

comparisons

enrich_gene_list_for_test <- list()

# convert symbols to ENTREZID
for (key in comparisons){ 
    symbols <- res[which(res$Parent.files==key),]$gene
    symbols <- gsub("--[0-9]*",'',symbols)
    unique_symbols <- unique(unlist(strsplit(symbols, split = ","), recursive = FALSE))
    df <- bitr(unique_symbols, fromType = "SYMBOL",
        toType = c("ENSEMBL", "ENTREZID"),
        OrgDb = org.Mm.eg.db)
    enrich_gene_list_for_test[[ key ]] <- df$ENTREZID
    #enrich_gene_list_for_test[[ key ]] <- unique_symbols
}

enrich_GO <- compareCluster(geneCluster = enrich_gene_list_for_test, 
                                   ont           = "ALL",
                                   OrgDb = org.Mm.eg.db,
                                   keyType = 'ENTREZID',
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.05,
                                   fun =  enrichGO)

write.table(enrich_GO,"oxistress_vs_E16.5_GO_enrichment.txt",sep="\t",quote=F)

options(repr.plot.width = 8, repr.plot.height = 8, repr.plot.res = 300)
p1 <- dotplot(enrich_GO,label_format = 100 ,showCategory=20, font.size=8)
pdf("oxistress_vs_E16.5_GO_enrichment.pdf",height=8,width=8)
p1
dev.off()
p1

enrich_kegg <- compareCluster(geneCluster = enrich_gene_list_for_test, 
                                   organism = "mmu",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.05,
                                   fun =  enrichKEGG)

write.table(enrich_kegg,"oxistress_vs_E16.5_KEGG_enrichment.txt",sep="\t",quote=F)

options(repr.plot.width = 9, repr.plot.height = 8, repr.plot.res = 300)
p2 <- dotplot(enrich_kegg,label_format = 100 ,showCategory=20, font.size=8)
pdf("oxistress_vs_E16.5_KEGG_enrichment.pdf",height=8,width=9)
p2
dev.off()
p2


