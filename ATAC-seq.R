library(DiffBind)
library(dplyr)
library(RUVSeq)


samples <- read.csv("SH_ATAC_sample_all_2.csv")
result1 <- dba(sampleSheet = "SH_ATAC_sample_all_2.csv")
result1
plot(result1)
result2 <- dba.count(result1, summits = 250)
result2
plot(result2)
#batch correction-------------
sample1 <- result2$peaks[[1]]
for( i in c(2:20) ){
  assign(paste0("sample",i),result2$peaks[[i]])
}

Read_all <- cbind(sample1$Reads,sample2$Reads,sample3$Reads,sample4$Reads,sample5$Reads,
                  sample6$Reads,sample7$Reads,sample8$Reads,sample9$Reads,sample10$Reads,
                  sample11$Reads,sample12$Reads,sample13$Reads,sample14$Reads,sample15$Reads,
                  sample16$Reads,sample17$Reads,sample18$Reads,sample19$Reads,sample20$Reads)

colnames(Read_all) <- as.character(samples$SampleID)

x1 <- as.factor(c(rep(1,2),rep(2,4),rep(1,6),rep(2,8)))
set <- newSeqExpressionSet(as.matrix(Read_all),phenoData
                           = data.frame(x1,row.names=colnames(Read_all)))

library(RColorBrewer)
library(ggrepel)

colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x1])
plotPCA(set,col=colors[x1], cex=1)

set2 <- betweenLaneNormalization(set, which="full")
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x1])
plotPCA(set2, col=colors[x1], cex=1,labels=TRUE)
pData(set2)

dat2 <- normCounts(set2)


library(DESeq2)
sample_info <- read.table("meta_ATAC_2.txt", header = TRUE, as.is = FALSE, row.names = 1)
dds <- DESeqDataSetFromMatrix(countData = dat2, #counts(set2),
                              colData = sample_info,
                              design = ~type)

vsd <- varianceStabilizingTransformation(dds)
plotPCA(vsd,intgroup="type")
pcaData <- plotPCA(vsd,intgroup="type", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(y = PC2 + 0.5, label = type), size = 3, vjust = 0)+
  theme_bw()

