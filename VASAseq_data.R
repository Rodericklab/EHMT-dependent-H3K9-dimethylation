# Install required packages if not already installed

if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  install.packages("SingleCellExperiment")
}
if (!requireNamespace("scater", quietly = TRUE)) {
  install.packages("scater")
}
if (!requireNamespace("platetools", quietly = TRUE)) {
  install.packages("platetools")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("viridis", quietly = TRUE)) {
  install.packages("viridis")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
if (!requireNamespace("matrixStats", quietly = TRUE)) {
  install.packages("matrixStats")
}
if (!requireNamespace("edgeR", quietly = TRUE)) {
  install.packages("edgeR")
}
if (!requireNamespace("limma", quietly = TRUE)) {
  install.packages("limma")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}


library(Seurat)
library(SingleCellExperiment)
library(scater)
library(platetools)
library(ggplot2)
library(viridis)
library(dplyr)
library(patchwork)
library(matrixStats)
library(edgeR)
library(limma)
library(pheatmap)

#load VASAseq MI data from the GEOXXX

data_dir <- "~"
data_file <- file.path(data_dir, "WT_KO_MI.csv")

# Read the data
MI <- read.csv(data_file, header = TRUE)
row.names(MI) <- MI$X
MI <- MI[, -c(1)]

Treatment <- factor(c("WT","WT","WT","KO","KO","KO"))

test <- data.frame(Sample=colnames(MI),Treatment)
test$Treatment <- relevel(test$Treatment, ref="WT")
rownames(test) <- test$Sample

edgeR.DGElist1 <- DGEList(counts = MI)
head(edgeR.DGElist1$counts)
edgeR.DGElist1$samples
hist(log2(rowSums(edgeR::cpm(edgeR.DGElist1))))

edgeR.DGElist1 <- calcNormFactors(edgeR.DGElist1 , method = "TMM")
edgeR.DGElist1$samples
plotMDS(edgeR.DGElist1, col=c(rep("Red",3), rep("Green",3)))

Ca <- round(cor(cpm(edgeR.DGElist1, log = TRUE)), 6)
library(gplots)
heatmap.2(Ca, symm=TRUE, margin=c(12,12), trace="none", col=colorRampPalette(c("white", "grey", "black"))(150))


design <- model.matrix(~test$Treatment)
rownames(design) <- colnames(MI)

edgeR.DGElist1 <- estimateGLMCommonDisp(edgeR.DGElist1, design)
edgeR.DGElist1 <- estimateGLMTrendedDisp(edgeR.DGElist1, design)
edgeR.DGElist1 <- estimateGLMTagwiseDisp(edgeR.DGElist1, design)


fit <- glmQLFit(edgeR.DGElist1, design, robust=TRUE)
plotQLDisp(fit)
colnames(fit)
qlf.2vs1 <- glmQLFTest(fit, coef=2)
topTags(qlf.2vs1, n=90)

results_MI <- as.data.frame(topTags(qlf.2vs1,n = 27121))

# Adjust FDR 0 for plotting
small_value <- 1.1e-190
results_MI$FDR <- ifelse(results_MI$FDR == 0, results_MI$FDR + small_value, results_MI$FDR)

results_MI <- subset(results_MI, FDR < 0.05)
up <- subset(results_MI, results_MI$logFC > 0.25)
down <- subset(results_MI, results_MI$logFC < -0.25)
results_MI <- rbind(up, down)
results_MI$gene <- rownames(results_MI)



#volcano plot

genes_to_label <- c('Anln', 'Aurka', 'Cdc25b', 'Cdk1', 'Kif11', 'Kif23', 'Mki67', 'Nuf2', 'Smc4', 'Top2a', 'Myh7', 'Col1a1', 'Col1a2','Nppa', 'Nppb', 'Kcnj3', 'Ndufb9', 'Ndufa11')
label_data <- results_MI[results_MI$gene %in% genes_to_label, ]

# Add a column to classify genes as up or downregulated
results_MI$regulation <- ifelse(results_MI$logFC > 0.25, "Upregulated",
                                 ifelse(results_MI$logFC < -0.25, "Downregulated", "Neutral"))


library(ggplot2)
library(ggrepel)


# Create the plot
volcano_plot <- ggplot(results_MI, aes(x = logFC, y = -log10(FDR), color = regulation)) +
  geom_point(size = 2) +
  geom_text_repel(
    data = label_data,
    aes(label = gene),
    box.padding   = 0.35, 
    point.padding = 0.2,
    segment.color = 'black', # Set segment color
    color = "black",          # Set text color to black
    size = 5,                 # Adjust text size
    max.overlaps = Inf        # Allow infinite overlaps
  ) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  labs(title = "",
       x = "log2(fold change)",
       y = "log10(pvalue)") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, colour = 'black', face = 'bold'), # Adjust axis text size
    legend.title = element_blank(),
    legend.position = "none" # Adjust legend position if necessary
  ) +
  xlim(-8, 8) + # Set x-axis limits
  ylim(0, 350) +
  guides(color = guide_legend(title = "Regulation")) # Adjust legend title

# Print the plot
print(volcano_plot)


#heatmap of some interest cell cycle genes 

int_cc <- c('Aurka', 'Mki67', 'Anln', 'Cdk1', 'Top2a', 'Gtse1', 'Ckap2', 'Nuf2', 'Kif11', 'Cdc25b', 'Cdc25c', 'Kif23', 'Smc4')

logcpm <- cpm(edgeR.DGElist1, log = TRUE)
x  <- as.matrix(logcpm)
x <- x[(rownames(x) %in% int_cc),]
pheatmap(x, cluster_cols = T,cluster_rows = F, clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", clustering_method = "complete", scale = 'row', color = colorRampPalette(c("navy", "white", "firebrick3"))(50))





#load VASAseq Sham data from the GEOXXX

data_dir <- "~"
data_file <- file.path(data_dir, "WT_KO_Sham.csv")

# Read the data
Sham <- read.csv(data_file, header = TRUE)
row.names(Sham) <- Sham$X
Sham <- Sham[, -c(1)]

Treatment <- factor(c("WT","WT","KO","KO","KO"))

test <- data.frame(Sample=colnames(Sham),Treatment)
test$Treatment <- relevel(test$Treatment, ref="WT")
rownames(test) <- test$Sample

edgeR.DGElist1 <- DGEList(counts = Sham)
head(edgeR.DGElist1$counts)
edgeR.DGElist1$samples
hist(log2(rowSums(edgeR::cpm(edgeR.DGElist1))))

edgeR.DGElist1 <- calcNormFactors(edgeR.DGElist1 , method = "TMM")
edgeR.DGElist1$samples
plotMDS(edgeR.DGElist1, col=c(rep("Red",2), rep("Green",3)))

Ca <- round(cor(cpm(edgeR.DGElist1, log = TRUE)), 5)
library(gplots)
heatmap.2(Ca, symm=TRUE, margin=c(12,12), trace="none", col=colorRampPalette(c("white", "grey", "black"))(150))


design <- model.matrix(~test$Treatment)
rownames(design) <- colnames(Sham)

edgeR.DGElist1 <- estimateGLMCommonDisp(edgeR.DGElist1, design)
edgeR.DGElist1 <- estimateGLMTrendedDisp(edgeR.DGElist1, design)
edgeR.DGElist1 <- estimateGLMTagwiseDisp(edgeR.DGElist1, design)

edgeR.DGElist1$common.dispersion
plotBCV(edgeR.DGElist1)

fit <- glmQLFit(edgeR.DGElist1, design, robust=TRUE)
plotQLDisp(fit)
colnames(fit)
qlf.2vs1 <- glmQLFTest(fit, coef=2)
topTags(qlf.2vs1, n=90)

results_Sham <- as.data.frame(topTags(qlf.2vs1,n = 27121))

small_value <- 10e-305
results_Sham$FDR <- ifelse(results_Sham$FDR == 0, results_Sham$FDR + small_value, results_Sham$FDR)

results_Sham <- subset(results_Sham, FDR < 0.05)
up <- subset(results_Sham, results_Sham$logFC > 0.25)
down <- subset(results_Sham, results_Sham$logFC < -0.25)
results_Sham <- rbind(up, down)
results_Sham$gene <- rownames(results_Sham)


#volcano plot

genes_to_label <- c('Myh7', 'Aurkb', 'Bub1', 'Ccnb2', 'Cdc25b', 'Cdc25c', 'Dlgap5', 'G2e3', 'Hells', 'Kif23', 'Ndc80', 'Prim1', 'Top2a', 'Col3a1', 'Col1a1', 'Col1a2', 'Runx1', 'Nppa', 'Nppb', 'Nppc')
label_data <- results_Sham[results_Sham$gene %in% genes_to_label, ]

# Add a column to classify genes as up or downregulated
results_Sham$regulation <- ifelse(results_Sham$logFC > 0.25, "Upregulated",
                                   ifelse(results_Sham$logFC < -0.25, "Downregulated", "Neutral"))


library(ggplot2)
library(ggrepel)

# label_data should be a subset of results_Sham that contains only the points you want to label.

# Create the plot
volcano_plot <- ggplot(results_Sham, aes(x = logFC, y = -log10(FDR), color = regulation)) +
  geom_point(size = 2) +
  geom_text_repel(
    data = label_data,
    aes(label = gene),
    box.padding   = 0.35, 
    point.padding = 0.2,
    segment.color = 'black', # Set segment color
    color = "black",          # Set text color to black
    size = 5,                 # Adjust text size
    max.overlaps = Inf        # Allow infinite overlaps
  ) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  labs(title = "",
       x = "log2(fold change)",
       y = "log10(pvalue)") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, colour = 'black', face = 'bold'), # Adjust axis text size
    legend.title = element_blank(),
    legend.position = "none" # Adjust legend position if necessary
  ) +
  xlim(-8, 8) + # Set x-axis limits
  ylim(0, 350) +
  guides(color = guide_legend(title = "Regulation")) # Adjust legend title

# Print the plot
print(volcano_plot)


#heatmap of some interest cell cycle genes 

int_cc <- c('Aurkb', 'Ccnb2', 'G2e3', 'Top2a', 'Cdc25b', 'Cdc25c', 'Prim1', 'Gmnn', 'Hells', 'Kif23', 'Bub1', 'Ube2c', 'Tmpo', 'Dlgap5', 'Tlk', 'Ndc80')


logcpm <- cpm(edgeR.DGElist1, log = TRUE)
x  <- as.matrix(logcpm)
x <- x[(rownames(x) %in% int_cc),]
pheatmap(x, cluster_cols = T,cluster_rows = F, clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", clustering_method = "complete", scale = 'row', color = colorRampPalette(c("navy", "white", "firebrick3"))(50))



