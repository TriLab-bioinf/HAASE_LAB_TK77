library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(vioplot)
library(ggdist)
library(ggforce)
library(dplyr)
library(tidyr)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(ggrepel)

files1 <- list.files("/gpfs/gsfs12/users/wangy80/TK77/ARTDeco_stress_output/dogs",pattern="*out.dogs.bed")
files2 <- list.files("/gpfs/gsfs12/users/wangy80/TK77/ARTDeco_Output_v2/dogs",pattern="*out.dogs.bed")

dir1 <- "/gpfs/gsfs12/users/wangy80/TK77/ARTDeco_stress_output/dogs/"
dir2 <- "/gpfs/gsfs12/users/wangy80/TK77/ARTDeco_Output_v2/dogs/"

files<- c(paste0(dir1,files1),paste0(dir2,files2))

load(paste0(getwd(), "/MILI.prepach.clusters_forThenia.Rdata")) #MILI.prepach.clusters
MILI.prepach.clusters <- MILI.prepach.clusters$regions

write.table(as.data.frame(MILI.prepach.clusters),"MILI.prepach.clusters.bed",sep="\t",quote=F,row.names=F)

print(paste("Number of MILI prepachytene clusters:", length(MILI.prepach.clusters)))

#sort and rank prepachytene clusters
MILI.prepach.clusters <- MILI.prepach.clusters[order(-MILI.prepach.clusters$all_reads_primary_alignments_FPM)]
MILI.prepach.clusters$Rank <- seq(1: length(MILI.prepach.clusters))

top1000 <- MILI.prepach.clusters[1:1000,]

sta <- data.frame()

for(i in 1:length(files)){
    samplename <- gsub("_vSplGTF","",gsub("^.*totalRNA_mouse_","",gsub("^.*mmu_NIH3TR_","",gsub("_Aligned.sortedByCoord.out.dogs.bed","",files[i]))))
    data <- import(files[i], format="bed")
    data$length <- width(data)
    data$sample <- samplename
    print(paste(samplename, length(data)))
    sta[i,1] <- samplename
    sta[i,2] <- length(data)
    ## overlap with PiRNA clusters
    df <- as.data.frame(data)
    df$ID <- paste0(df$seqnames,":",df$start,"-",df$end)
    overlap <- as.data.frame(subsetByOverlaps(data, top1000,type=c("any"), ignore.strand=FALSE, invert=FALSE))
    overlap$ID <- paste0(overlap$seqnames,":",overlap$start,"-",overlap$end)
    df$PiCluster <- ifelse(df$ID %in% overlap$ID, "yes", "no")
    assign(samplename,df)
}

### number of DoGS for each sample

sta$V2 <- as.numeric(sta$V2)
options(repr.plot.width = 5, repr.plot.height = 3, repr.plot.res = 300)

p1<-ggplot(data=sta, aes(x=V1, y=V2)) +
  geom_bar(stat="identity",fill="steelblue") + 
  coord_flip() +
  geom_text(aes(label=V2), hjust=1.2, size=3,color="white")+
  xlab("sample") +
  ylab("count") +
  ggtitle("number of DoGs") +
  theme_minimal()

pdf("number_of_DoGs_barplot.pdf",height=3,width=5)
p1
dev.off()
p1

## linux
# for i in `ls *bam`; do file=`basename $i _Aligned.sortedByCoord.out.bam`; n=`samtools view $i|wc -l`; echo -e "$file $n" >>sta.txt; done

reads <- read.delim("/data/wangy80/TK77/FranziToYuejunReadthrough/sta.txt",sep=" ",header=F)
names(reads) <- c("filename","readcount")

reads$sample <- gsub("_vSplGTF","",gsub("totalRNA_mouse_","",gsub("mmu_NIH3TR_","",reads$filename)))

options(repr.plot.width = 5, repr.plot.height = 3, repr.plot.res = 300)
p2<-ggplot(data=reads, aes(x=sample, y=readcount)) +
  geom_bar(stat="identity",fill="steelblue") + 
  coord_flip() +
  geom_text(aes(label=readcount), hjust=0.6, size=1.5,color="black")+
  xlab("sample") +
  ylab("readcount") +
  ggtitle("number of reads") +
  theme_minimal()

pdf("number_of_reads_barplot.pdf",height=3,width=5)
p2
dev.off()
p2

#load gene.exp.fpkm.txt file
file_geneExprFPKM <- "/gpfs/gsfs12/users/wangy80/TK77/ARTDeco_stress_output/quantification/gene.exp.fpkm.txt"
trExprFPKM <- read.delim(file_geneExprFPKM)
nrow(trExprFPKM[order(trExprFPKM$ID, decreasing=TRUE), ])

#load gtf file and extract gene_name
gtf.file = "/gpfs/gsfs12/users/wangy80/TK77/Mm10_refSeq3_copies_annotated3.sorted.gtf"
gtf.gr = rtracklayer::import(gtf.file) # creates a GRanges object

gtf <- as.data.frame(gtf.gr)

gtf$ID <- gtf$transcript_id
anno <- unique(gtf[,c("gene_id","ID")])

trExprFPKM_2 <- left_join(trExprFPKM,anno,by="ID")

#rename "ID" column to transcript_id
colnames(trExprFPKM_2)[colnames(trExprFPKM_2) == "ID"] <- "transcript_id"
trExprFPKM_2 <- trExprFPKM_2[order(trExprFPKM_2$transcript_id, decreasing=TRUE), ]
nrow(trExprFPKM_2)

trExprFPKM.m <- melt(trExprFPKM_2,id.vars=c("transcript_id","Length","gene_id"))

GeneExprFPKM.m <- trExprFPKM.m %>% 
    group_by(gene_id,variable) %>%
    arrange(-value) %>%
    slice(1)

dim(GeneExprFPKM.m)

dim(trExprFPKM.m)

#load gene.exp.fpkm.txt file
file_geneExprFPKM2 <- "/gpfs/gsfs12/users/wangy80/TK77/ARTDeco_Output_v2/quantification/gene.exp.fpkm.txt"
trExprFPKM2 <- read.delim(file_geneExprFPKM2)
nrow(trExprFPKM2[order(trExprFPKM2$ID, decreasing=TRUE), ])

trExprFPKM2_2 <- left_join(trExprFPKM2,anno,by="ID")

#rename "ID" column to transcript_id
colnames(trExprFPKM2_2)[colnames(trExprFPKM2_2) == "ID"] <- "transcript_id"
trExprFPKM2_2 <- trExprFPKM2_2[order(trExprFPKM2_2$transcript_id, decreasing=TRUE), ]
nrow(trExprFPKM2_2)

trExprFPKM2.m <- melt(trExprFPKM2_2,id.vars=c("transcript_id","Length","gene_id"))

GeneExprFPKM2.m <- trExprFPKM2.m %>% 
    group_by(gene_id,variable) %>%
    arrange(-value) %>%
    slice(1)

dim(trExprFPKM2.m)

allGeneExprFPKM.m <- rbind(GeneExprFPKM.m,GeneExprFPKM2.m)

allGeneExprFPKM.m %>% 
group_by(variable) %>%
filter(value>0.1) %>%
count()

allGeneExprFPKM.m$sample <- gsub("_vSplGTF","",gsub("totalRNA_mouse_","",gsub("mmu_NIH3TR_","",gsub("_Aligned.sortedByCoord.out","",allGeneExprFPKM.m$variable))))

all_dogs <- rbind(
    heatshock_SRR5558715,
    heatshock_SRR5558716,
    osmstress_SRR5558719,
    osmstress_SRR5558720,
    oxistress_SRR5558717,
    oxistress_SRR5558718,
    untreated_SRR5558713,
    untreated_SRR5558714,
    control_P42_SRR765631,
    pach_P14_SRR7760359,
    prepach_E16.5_SRR11916388
)

all_dogs$gene_id <- all_dogs$name

res <- left_join(allGeneExprFPKM.m,all_dogs,by=c("sample","gene_id"))

res$DoGs <- ifelse(!is.na(res$name),"Yes","No")
res$expr <- ifelse(res$value>0.1,"expressed","non-expressed")

res$group <- case_when(res$expr == "expressed" & res$DoGs=="Yes" ~ "expressed DoGs",
                       res$expr == "expressed" & res$DoGs=="No" ~ "expressed genes",
                       res$expr == "non-expressed" & res$DoGs=="Yes" ~ "non-expressed DoGs",
                       res$expr == "non-expressed" & res$DoGs=="No" ~ "non-expressed genes"
)

res2 <- res %>% 
group_by(sample,group)  %>%
summarize(count=n())

res3 <- res2 %>% pivot_wider(names_from = group, values_from = count)
res3

# enrichment analysis of DoGs genes
rownames(res3) <- res3$sample
apply(res3[,-1], 1, function(x) fisher.test(matrix(as.numeric(x[1:4]), ncol=2, byrow=T)))

options(repr.plot.width = 7, repr.plot.height = 5, repr.plot.res = 300)
pie <- ggplot(res2, aes(x=sample, y=count, fill=group)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) +
    scale_fill_manual(values=c("#56B4E9", "#E69F00", "#999999","#009966"), name="Category") +
    theme_minimal() +
    geom_text(aes(label = count),
        position = position_stack(vjust = 0.5),
        show.legend = FALSE, 
        size = 1) 

pdf("pie.pdf",height=5,width=7)
pie
dev.off()
pie

options(repr.plot.width = 7, repr.plot.height = 3, repr.plot.res = 300)
bar <- ggplot(res2, aes(x=count, y=sample,fill=group))+
  geom_col() +
  geom_text(aes(label = count,x=count),
           position = position_stack(vjust = 0.5),
           size=2)+
  scale_fill_brewer(palette="Set3") +
  theme_minimal()

pdf("bar.pdf",height=3,width=7)
bar
dev.off()
bar

options(repr.plot.width = 5, repr.plot.height = 5, repr.plot.res = 300)
library(scales) # to access break formatting functions

graph1 <- ggplot(all_dogs, aes(x=sample, y=log10(length), fill=sample)) + 
            geom_violin(width=1, scale="count") +
            geom_boxplot(width=0.2, color="grey36", alpha=0.6, outlier.shape = NA) +
            annotation_logticks(
                base = 10,
                sides = "l", 
                colour = "black") +
            theme_classic() + 
            coord_cartesian(ylim = c(2.7, 5.2)) +
            scale_fill_brewer(palette="Set3") +
            labs(title = "Log10 Transformed Violin Plot of nucleotide length of DoGs", y = "Log10(length of DoG)") + 
            theme(text=element_text(size=10,  family="Helvetica"), legend.position = "none",
                  axis.text.x = element_text(angle=90,hjust=1,vjust=0.5,size=6)
                 )+
            stat_compare_means(ref.group="prepach_E16.5_SRR11916388",aes(label = ..p.signif..))

pdf("DoGs_lengths.pdf",height=5,width=5)
graph1
dev.off()
graph1

res4 <- res %>% filter(DoGs=="Yes")

options(repr.plot.width = 8, repr.plot.height = 4, repr.plot.res = 300)
graph2 <- ggplot(res4, aes(x=sample, y=log10(value), fill=PiCluster)) + 
            #geom_violin() +
            geom_boxplot() +
            annotation_logticks(
                base = 10,
                sides = "l", 
                colour = "black") +
            theme_classic() + 
            #coord_cartesian(ylim = c(2.7, 5.2)) +
            scale_fill_brewer(palette="Set3") +
            labs(title = "Log10 Transformed Expression", y = "Log10(Expression)") +
            theme(text=element_text(size=10,  family="Helvetica"), 
                  axis.text.x = element_text(angle=90,hjust=1,vjust=0.5,size=6)
                 ) +
            stat_compare_means(aes(group = PiCluster,label = ..p.signif..))

graph2

options(repr.plot.width = 8, repr.plot.height = 5, repr.plot.res = 300)
graph3 <- ggplot(res4, aes(x=sample, y=log10(length), fill=PiCluster)) + 
            #geom_violin() +
            geom_boxplot() +
            annotation_logticks(
                base = 10,
                sides = "l", 
                colour = "black") +
            theme_classic() + 
            #coord_cartesian(ylim = c(2.7, 5.2)) +
            scale_fill_brewer(palette="Set3") +
            labs(title = "Log10 Transformed Length", y = "Log10(Length)") +
            theme(text=element_text(size=10,  family="Helvetica"), 
                  axis.text.x = element_text(angle=90,hjust=1,vjust=0.5,size=6)
                 ) +
            stat_compare_means(aes(group = PiCluster,label = ..p.signif..))

pdf("DoGs_overlapped_with_top1000_Piclusters.pdf",height=5,width=8)
graph3
dev.off()
graph3

res4$group2 <- case_when(res4$expr == "expressed" & res4$PiCluster=="yes" ~ "expressed DoGs in PiClusters",
                       res4$expr == "expressed" & res4$PiCluster=="no" ~ "expressed DoGs not in PiClusters",
                       res4$expr == "non-expressed" & res4$PiCluster=="yes" ~ "non-expressed DoGs in PiClusters",
                       res4$expr == "non-expressed" & res4$PiCluster=="no" ~ "non-expressed DoGs not in PiClusters"
)

res5 <- res4 %>% 
group_by(sample,group2)  %>%
summarize(count=n())

res6 <- res5 %>% pivot_wider(names_from = group2, values_from = count)
res6

# enrichment analysis of DoGs genes
rownames(res6) <- res6$sample
apply(res6[,-1], 1, function(x) fisher.test(matrix(as.numeric(x[1:4]), ncol=2, byrow=T)))

