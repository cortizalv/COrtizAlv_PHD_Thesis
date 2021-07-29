setwd("Path/to/datasets")
library(DESeq2)
library("dplyr")
library("ggplot2")
library("RColorBrewer")
library(ppcor)
library(qgraph)

gene.counts = as.matrix(read.csv(file = "gene_counts_PR_bees_NFO_only.csv", row.names = "GeneID"))

sample.names = as.factor(c("C1-12-1",
                           "C1-12-2",
                           "C1-12-3",
                           "C1-12-4",
                           "C1-21-1",
                           "C1-21-2",
                           "C1-21-3",
                           "C1-21-4",
                           "C3-12-1",
                           "C3-12-2",
                           "C3-12-3",
                           "C3-12-4",
                           "C3-21-1",
                           "C3-21-2",
                           "C3-21-3",
                           "C3-21-4",
                           "C4-12-1",
                           "C4-12-2",
                           "C4-21-1",
                           "C4-21-2",
                           "T1-70-1",
                           "T1-70-2",
                           "T1-70-3",
                           "T1-70-4",
                           "T3-70-3",
                           "T3-70-4",
                           "T4-70-1",
                           "T4-70-2"))

trt.comb = as.factor(c("nurse",
                       "nurse",
                       "nurse",
                       "nurse",
                       "forager",
                       "forager",
                       "forager",
                       "forager",
                       "nurse",
                       "nurse",
                       "nurse",
                       "nurse",
                       "forager",
                       "forager",
                       "forager",
                       "forager",
                       "nurse",
                       "nurse",
                       "forager",
                       "forager",
                       "old",
                       "old",
                       "old",
                       "old",
                       "old",
                       "old",
                       "old",
                       "old"))

sample.info.comb = data.frame(sample.names, 
                              trt.comb,
                              row.names = sample.names)

ddsMat.comb = DESeqDataSetFromMatrix(countData = gene.counts, 
                                     colData = sample.info.comb, 
                                     design = ~ trt.comb)

nrow(ddsMat.comb) #number of rows pre removal of small counts
keep.comb <- rowSums(counts(ddsMat.comb)) > 1 #keep all rows which count total sums to more than 1
ddsMat.comb <- ddsMat.comb[keep.comb,] #new ddsMat object with removed rows
nrow(ddsMat.comb) #number of rows after removal of small counts

vsd.comb <- vst(ddsMat.comb, blind = FALSE) #blind=FALSE as so that the variables in the design do not contribute to the expected variance-mean trend of the experiment.
head(assay(vsd.comb), 3)
colData(vsd.comb)

######                                         ######
######                RESULTS                  ######
######                                         ######

pcaData <-  plotPCA(vsd.comb, intgroup = c("trt.comb"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = trt.comb)) +
  geom_point(size =2) +
  expand_limits(x=c(-35,35)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  stat_ellipse() +
  scale_color_manual(breaks = c("nurse", "forager", "old"), 
                    values=c("#FF934F", "#058ED9", "#2D3142")) +
  theme(legend.position = c(1.05, 0.5), panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "white"), legend.key = element_rect(fill = NA, color = NA)) +
  labs(colour = "Worker type", title = "PCA of overall gene counts \n across worker type") 

ddsMatres.comb = DESeq(ddsMat.comb)
res.comb <- results(ddsMatres.comb, alpha = 0.05, contrast=c("trt.comb","nurse","forager")) #call as is will save the results from the last variable of the design formula. In this case we only have one: trt.grps
res.comb

clock.genes = c("GeneID_410197", "GeneID_550811", "GeneID_410757", "GeneID_725614", "GeneID_410253", "GeneID_408449", "GeneID_406112", "GeneID_726204")
vsdclock =  subset(vsd.comb, rownames(vsd.comb) %in% clock.genes)

pcaData.clock <-  plotPCA(vsdclock, intgroup = c("trt.comb"), returnData = TRUE)
percentVar.clock <- round(100 * attr(pcaData.clock, "percentVar"))

ggplot(pcaData.clock, aes(x = PC1, y = PC2, color = trt.comb)) +
  geom_point(size =2) +
  xlab(paste0("PC1: ", percentVar.clock[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar.clock[2], "% variance")) +
  coord_fixed() +
  stat_ellipse() +
  scale_color_manual(breaks = c("nurse", "forager", "old"), 
                     values=c("#FF934F", "#058ED9", "#2D3142")) +
  theme(legend.position = c(1.10, 0.5), panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "white"), legend.key = element_rect(fill = NA, color = NA)) +
  labs(colour = "Worker type", title = "PCA of clock gene counts \n across worker type")
  
###Partial correlation Netwirk Analysis
data.clock.nurse2 = read.csv(file = "count_data_clock_genes_nurse_inverse.csv", row.names = "SampleID")
data.clock.forager2 = read.csv(file = "count_data_clock_genes_foragers_inverse2.csv", row.names = "SampleID")
data.clock.old2 = read.csv(file = "count_data_clock_genes_old_inverse2.csv", row.names = "SampleID")

x.nurse2 = pcor(data.clock.nurse2, method = c("spearman"))
y.forager2 = pcor(data.clock.forager2, method = c("spearman"))
z.old2 = pcor(data.clock.old2, method = c("spearman"))

qgraph(x.nurse2$estimate, threshold = 0.80, minimum = 0.80, theme = "colorblind")
qgraph(y.forager2$estimate, threshold = 0.80, minimum = 0.80, theme = "colorblind")
qgraph(z.old2$estimate, threshold = 0.80, minimum = 0.80, theme = "colorblind")
