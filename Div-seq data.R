# Install required packages if not already installed
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
if (!requireNamespace("edgeR", quietly = TRUE)) {
  install.packages("edgeR")
}
if (!requireNamespace("enrichR", quietly = TRUE)) {
  install.packages("enrichR")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}
if (!requireNamespace("forcats", quietly = TRUE)) {
  install.packages("forcats")
}
if (!requireNamespace("gplots", quietly = TRUE)) {
  install.packages("gplots")
}



# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)
library(edgeR)
library(enrichR)
library(RColorBrewer)
library(forcats) # for fct_reorder()
library(gplots)


#load Div seq mouse CM nuclei (Edu +ve/-ve A366) data from the GSE261924

data_dir <- "~"
data_file <- file.path(data_dir, "Mouse_merged_counts_Edu+ve_vs_Edu-ve_A366.csv")

# Read the data
A366 <- read.csv(data_file, header = TRUE)
rownames(A366) <- A366$X
A366 <- A366[, -c(1)]

Animal <- factor(c("F","N","O","F","N","O"))
Treatment <- factor(c("Edu_pos","Edu_pos","Edu_pos","Edu_neg","Edu_neg","Edu_neg"))

test <- data.frame(Sample=colnames(A366),Animal,Treatment)
test$Treatment <- relevel(test$Treatment, ref="Edu_neg")
rownames(test) <- test$Sample

edgeR.DGElist1 <- DGEList(counts = A366)
head(edgeR.DGElist1$counts)
edgeR.DGElist1$samples
hist(log2(rowSums(cpm(edgeR.DGElist1))))

keep <- filterByExpr(edgeR.DGElist1)
table(keep)
edgeR.DGElist1 <- edgeR.DGElist1[keep, , keep.lib.sizes=FALSE]
hist(log2(rowSums(cpm(edgeR.DGElist1))))


edgeR.DGElist1 <- calcNormFactors(edgeR.DGElist1 , method = "TMM")
edgeR.DGElist1$samples
plotMDS(edgeR.DGElist1, col=c(rep("Red",3), rep("Green",3)))


Ca <- round(cor(cpm(edgeR.DGElist1, log = TRUE)), 6)

heatmap.2(Ca, symm=TRUE, margin=c(12,12), trace="none", col=colorRampPalette(c("white", "grey", "black"))(150))


design <- model.matrix(~test$Treatment)
rownames(design) <- colnames(A366)

edgeR.DGElist1 <- estimateGLMCommonDisp(edgeR.DGElist1, design)
edgeR.DGElist1 <- estimateGLMTrendedDisp(edgeR.DGElist1, design)
edgeR.DGElist1 <- estimateGLMTagwiseDisp(edgeR.DGElist1, design)

fit <- glmQLFit(edgeR.DGElist1, design, robust=TRUE)
plotQLDisp(fit)
colnames(fit)
qlf.2vs1 <- glmQLFTest(fit, coef=2)
topTags(qlf.2vs1, n=20)

results_A366 <- as.data.frame(topTags(qlf.2vs1,n = 10864))
results_A366$gene <- row.names(results_A366)
results_A366 <- subset(results_A366, results_A366$FDR < 0.05)


#volcano plot
genes_to_label <- c('Myh6', 'Tnnt2', 'Atp1a2', 'Obscn', 'Etfdh', 'Kcnh2', 'Cenpa', 'Cenpe', 'Mcm8', 'Cdc7', 'Kif11', 'Ttk', 'Kif14', 'Cdca8', 'Pmf1', 'Top2a', 'Dmd', 'Ryr2', 'Slc8a1', 'kcnq1', 'Slc4a3', 'Myh7')
label_data <- results_A366[results_A366$gene %in% genes_to_label, ]

# Add a column to classify genes as up or downregulated
results_A366$regulation <- ifelse(results_A366$logFC > 0.25, "Upregulated",
                            ifelse(results_A366$logFC < -0.25, "Downregulated", "Neutral"))

results_A366 <- subset(results_A366, regulation %in% c('Upregulated', 'Downregulated'))


# Create the Volcano plot
# Adjusting the very near to 0 FDR values

small_value <- 8.919876e-9

results_A366$FDR <- ifelse(results_A366$FDR == 0, results_A366$FDR + small_value, results_A366$FDR)
results_A366$FDR <- ifelse(results_A366$FDR == 9.5687e-280, results_A366$FDR + small_value, results_A366$FDR)
results_A366$FDR <- ifelse(results_A366$FDR == 8.919876e-80, results_A366$FDR + small_value, results_A366$FDR)


volcano_plot <- ggplot(results_A366, aes(x = logFC, y = -log10(FDR), color = regulation)) +
  geom_point(size = 4) +
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
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, colour = 'black', face = 'bold'), # Adjust axis text size
    legend.title = element_blank(),
    legend.position = "none" # Adjust legend position if necessary
  ) +
  xlim(-1.5, 1.5) + # Set x-axis limits
  ylim(0, 80) # set y-axis limits

# Print the plot
print(volcano_plot)

#GO analysis

dbs <- c("GO_Cellular_Component_2021","GO_Biological_Process_2021" , "GO_Molecular_Function_2021","KEGG_2018")
up <- subset(results_A366, logFC > 0.25)
down <- subset(results_A366, logFC < -0.25)

up <- up$gene 
down <- down$gene

up <- enrichr(up, dbs)
down <- enrichr(down, dbs)

up <- up$GO_Biological_Process_2021
down <- down$GO_Biological_Process_2021


up$type <- "Up-regulated"
down$type <- "Down-regulated"

up <- up[order(up$Adjusted.P.value), ]  # sort
up <- subset(up, up$Term %in% c('microtubule cytoskeleton organization involved in mitosis (GO:1902850)','mitotic spindle organization (GO:0007052)', 'positive regulation of cell cycle process (GO:0090068)', 'G2/M transition of mitotic cell cycle (GO:0000086)', 'chromosome condensation (GO:0030261)'))

down <- down[order(down$Adjusted.P.value), ]  # sort
down <- subset(down, down$Term %in% c('response to muscle stretch (GO:0035994)', 'cardiac muscle contraction (GO:0060048)', 'actomyosin structure organization (GO:0031032)', 'cardiac muscle cell action potential (GO:0086001)', 'dephosphorylation (GO:0016311)'))

gos <- rbind(up,down)
gos$Term <- factor(gos$Term, levels=gos$Term)


# Filter the data for only "Up-regulated" genes
upregulated_gos <- gos[gos$type == "Up-regulated", ]

# Remove '(GO:' and any subsequent characters from the Term labels
upregulated_gos$Term <- gsub("\\(GO:.*$", "", upregulated_gos$Term)

# Order the terms based on the -log10(qvalue) from high to low
upregulated_gos <- upregulated_gos %>%
  mutate(Term = fct_reorder(Term, -log10(Adjusted.P.value)))

# Convert Overlap to numeric after removing the part after '/'
upregulated_gos$Overlap <- as.numeric(sub("/.*", "", upregulated_gos$Overlap))


plot <- ggplot(upregulated_gos, aes(y= Term, x = -log10(Adjusted.P.value), fill = Overlap)) + 
  geom_bar(stat = 'identity', width = .7, color = "black", linewidth = 1.05) + # Add black borders to bars
  scale_fill_gradientn(
    colors = c('orange', 'Brown'),  # Custom color gradient
    name = "Counts", 
    limits = range(upregulated_gos$Overlap, na.rm = TRUE),
    breaks = c(13, 22, 47),  # Specify breaks at 5, 6, and 8
    labels = c("13", "22", "47")  # Corresponding labels for these breaks
  ) +
  labs(x = "-log10(qvalue)", y = "") + 
  theme_bw() + # Start with theme_bw
  theme(axis.text = element_text(color = "black", size = 22, face = 'bold'), # Set axis text color to black
        axis.title = element_text(color = "black", size = 24, face = 'bold'), # Set axis title color to black
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        legend.position = c(0.83, 0.3),# Place legend inside the plot at the bottom right
        panel.border = element_rect(colour = "black", fill=NA, size=3.5), # Plot border
        legend.text = element_text(size= 22, color = 'black', face = 'bold'), # Legend text
        legend.title = element_text(size= 22, color = 'black', face = 'bold'),
        legend.key.size = unit(0.65, "cm"),
        legend.spacing.x = unit(0.5, "cm"),  # Horizontal spacing
        legend.spacing.y = unit(0.5, "cm"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt")))

# Print the plot
print(plot)


# Down_reg

downregulated_gos <- gos[gos$type == "Down-regulated", ]

# Remove '(GO:' and any subsequent characters from the Term labels
downregulated_gos$Term <- gsub("\\(GO:.*$", "", downregulated_gos$Term)

# Order the terms based on the -log10(qvalue) from high to low
downregulated_gos <- downregulated_gos %>%
  mutate(Term = fct_reorder(Term, -log10(Adjusted.P.value)))

# Convert Overlap to numeric after removing the part after '/'
downregulated_gos$Overlap <- as.numeric(sub("/.*", "", downregulated_gos$Overlap))


plot <- ggplot(downregulated_gos, aes(y= Term, x = -log10(Adjusted.P.value), fill = Overlap)) + 
  geom_bar(stat = 'identity', width = .7, color = "black", linewidth = 1.05) + # Add black borders to bars
  scale_fill_gradientn(
    colors = c('#B3E5FC', '#311B92'),  # Custom color gradient
    name = "Counts", 
    limits = range(downregulated_gos$Overlap, na.rm = TRUE),
    breaks = c(4, 5, 7),  # Specify breaks at 5, 6, and 8
    labels = c("4", "5", "7")  # Corresponding labels for these breaks
  ) +
  labs(x = "-log10(qvalue)", y = "") + 
  theme_bw() + # Start with theme_bw
  theme(axis.text = element_text(color = "black", size = 22, face = 'bold'), # Set axis text color to black
        axis.title = element_text(color = "black", size = 24, face = 'bold'), # Set axis title color to black
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        legend.position = c(0.83, 0.3),# Place legend inside the plot at the bottom right
        panel.border = element_rect(colour = "black", fill=NA, size=3.5), # Plot border
        legend.text = element_text(size= 22, color = 'black', face = 'bold'), # Legend text
        legend.title = element_text(size= 22, color = 'black', face = 'bold'),
        legend.key.size = unit(0.65, "cm"),
        legend.spacing.x = unit(0.5, "cm"),  # Horizontal spacing
        legend.spacing.y = unit(0.5, "cm"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt")))

# Print the plot
print(plot)



#NO_A366 data from the GEOXXX

data_dir <- "~"
data_file <- file.path(data_dir, "Mouse_merged_counts_Edu+ve_vs_Edu-ve_NO_A366.csv")

# Read the data
NO_A366 <- read.csv(data_file, header = TRUE)
rownames(NO_A366) <- NO_A366$X
NO_A366 <- NO_A366[, -c(1)]


Animal <- factor(c("H","Q","U","H","Q","U"))
Treatment <- factor(c("Edu_pos","Edu_pos","Edu_pos","Edu_neg","Edu_neg","Edu_neg"))

test <- data.frame(Sample=colnames(NO_A366),Animal,Treatment)
test$Treatment <- relevel(test$Treatment, ref="Edu_neg")
rownames(test) <- test$Sample

edgeR.DGElist2 <- DGEList(counts = NO_A366)
head(edgeR.DGElist2$counts)
edgeR.DGElist2$samples
hist(log2(rowSums(cpm(edgeR.DGElist2))))


edgeR.DGElist2 <- calcNormFactors(edgeR.DGElist2 , method = "TMM")
edgeR.DGElist2$samples
plotMDS(edgeR.DGElist2, col=c(rep("Red",3), rep("Green",3)))


Ca <- round(cor(cpm(edgeR.DGElist2, log = TRUE)), 6)

heatmap.2(Ca, symm=TRUE, margin=c(12,12), trace="none", col=colorRampPalette(c("white", "grey", "black"))(150))


design <- model.matrix(~test$Treatment)
rownames(design) <- colnames(NO_A366)

edgeR.DGElist2 <- estimateGLMCommonDisp(edgeR.DGElist2, design)
edgeR.DGElist2 <- estimateGLMTrendedDisp(edgeR.DGElist2, design)
edgeR.DGElist2 <- estimateGLMTagwiseDisp(edgeR.DGElist2, design)

fit <- glmQLFit(edgeR.DGElist2, design, robust=TRUE)
plotQLDisp(fit)
colnames(fit)
qlf.2vs1 <- glmQLFTest(fit, coef=2)
topTags(qlf.2vs1, n=20)

results_NO_A366 <- as.data.frame(topTags(qlf.2vs1,n = 10864))
results_NO_A366$gene <- row.names(results_NO_A366)
results_NO_A366 <- subset(results_NO_A366, results_NO_A366$FDR < 0.05)

#volcano plot
genes_to_label <- unique(c('Mcm10', 'Vrk1', 'Tlr5', 'Pelp1', 'Fancd2', 'Blm', 'Cacna1g', 'Rpl23a', 'Fgf14', 'Cenpa', 'Cenpe', 'Mcm8', 'Cdc7', 'Kif11', 'Ttk', 'Kif14', 'Cdca8', 'Pmf1', 'Top2a'))

label_data <- results_NO_A366[results_NO_A366$gene %in% genes_to_label, ]

# Add a column to classify genes as up or downregulated
results_NO_A366$regulation <- ifelse(results_NO_A366$logFC > 0.25, "Upregulated",
                                  ifelse(results_NO_A366$logFC < -0.25, "Downregulated", "Neutral"))

results_NO_A366 <- subset(results_NO_A366, regulation %in% c('Upregulated', 'Downregulated'))


# Create the Volcano plot
# Adjusting the very near to 0 FDR values

small_value <- 8.919876e-9

results_NO_A366$FDR <- ifelse(results_NO_A366$FDR == 0, results_NO_A366$FDR + small_value, results_NO_A366$FDR)


volcano_plot <- ggplot(results_NO_A366, aes(x = logFC, y = -log10(FDR), color = regulation)) +
  geom_point(size = 4) +
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
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, colour = 'black', face = 'bold'), # Adjust axis text size
    legend.title = element_blank(),
    legend.position = "none" # Adjust legend position if necessary
  ) +
  xlim(-1.5, 1.5) + # Set x-axis limits
  ylim(0, 30) # set y-axis limits

# Print the plot
print(volcano_plot)

#GO analysis

dbs <- c("GO_Cellular_Component_2021","GO_Biological_Process_2021" , "GO_Molecular_Function_2021","KEGG_2018")
up <- subset(results_NO_A366, logFC > 0.25)
down <- subset(results_NO_A366, logFC < -0.25)

up <- up$gene 
down <- down$gene

up <- enrichr(up, dbs)
down <- enrichr(down, dbs)

up <- up$GO_Biological_Process_2021
down <- down$GO_Biological_Process_2021


up$type <- "Up-regulated"
down$type <- "Down-regulated"

up <- up[order(up$Adjusted.P.value), ]  # sort
up <- subset(up, up$Term %in% c('mitotic sister chromatid segregation (GO:0000070)', 'mitotic spindle organization (GO:0007052)', 'mitotic cytokinesis (GO:0000281)', 'positive regulation of cell cycle process (GO:0090068)', 'centromere complex assembly (GO:0034508)'))

down <- down[order(down$Adjusted.P.value), ]  # sort
down <- subset(down, down$Term %in% c('regulation of intrinsic apoptotic signaling pathway by p53 class mediator (GO:1902253)', 'catecholamine metabolic process (GO:0006584)', 'mitochondrion organization (GO:0007005)', 'regulation of lipid transport (GO:0032368)', 'negative regulation of reactive oxygen species metabolic process (GO:2000378)'))

gos <- rbind(up,down)
gos$Term <- factor(gos$Term, levels=gos$Term)

# Filter the data for only "Up-regulated" genes
upregulated_gos <- gos[gos$type == "Up-regulated", ]

# Remove '(GO:' and any subsequent characters from the Term labels
upregulated_gos$Term <- gsub("\\(GO:.*$", "", upregulated_gos$Term)

# Order the terms based on the -log10(qvalue) from high to low
upregulated_gos <- upregulated_gos %>%
  mutate(Term = fct_reorder(Term, -log10(Adjusted.P.value)))

# Convert Overlap to numeric after removing the part after '/'
upregulated_gos$Overlap <- as.numeric(sub("/.*", "", upregulated_gos$Overlap))


plot <- ggplot(upregulated_gos, aes(y= Term, x = -log10(Adjusted.P.value), fill = Overlap)) + 
  geom_bar(stat = 'identity', width = .7, color = "black", linewidth = 1.05) + # Add black borders to bars
  scale_fill_gradientn(
    colors = c('orange', 'Brown'),  # Custom color gradient
    name = "Counts", 
    limits = range(upregulated_gos$Overlap, na.rm = TRUE),
    breaks = c(10, 15, 23),  # Specify breaks at 5, 6, and 8
    labels = c("10", "15", "23")  # Corresponding labels for these breaks
  ) +
  labs(x = "-log10(qvalue)", y = "") + 
  theme_bw() + # Start with theme_bw
  theme(axis.text = element_text(color = "black", size = 22, face = 'bold'), # Set axis text color to black
        axis.title = element_text(color = "black", size = 24, face = 'bold'), # Set axis title color to black
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        legend.position = c(0.83, 0.3),# Place legend inside the plot at the bottom right
        panel.border = element_rect(colour = "black", fill=NA, size=3.5), # Plot border
        legend.text = element_text(size= 22, color = 'black', face = 'bold'), # Legend text
        legend.title = element_text(size= 22, color = 'black', face = 'bold'),
        legend.key.size = unit(0.65, "cm"),
        legend.spacing.x = unit(0.5, "cm"),  # Horizontal spacing
        legend.spacing.y = unit(0.5, "cm"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt")))

# Print the plot
print(plot)


# Down_reg

downregulated_gos <- gos[gos$type == "Down-regulated", ]

# Remove '(GO:' and any subsequent characters from the Term labels
downregulated_gos$Term <- gsub("\\(GO:.*$", "", downregulated_gos$Term)

# Order the terms based on the -log10(qvalue) from high to low
downregulated_gos <- downregulated_gos %>%
  mutate(Term = fct_reorder(Term, -log10(P.value)))

# Convert Overlap to numeric after removing the part after '/'
downregulated_gos$Overlap <- as.numeric(sub("/.*", "", downregulated_gos$Overlap))


plot <- ggplot(downregulated_gos, aes(y= Term, x = -log10(P.value), fill = Overlap)) + 
  geom_bar(stat = 'identity', width = .7, color = "black", linewidth = 1.05) + # Add black borders to bars
  scale_fill_gradientn(
    colors = c('#B3E5FC', '#311B92'),  # Custom color gradient
    name = "Counts", 
    limits = range(downregulated_gos$Overlap, na.rm = TRUE),
    breaks = c(2, 3, 4),  # Specify breaks at 5, 6, and 8
    labels = c("2", "3", '4')  # Corresponding labels for these breaks
  ) +
  labs(x = "-log10(qvalue)", y = "") + 
  theme_bw() + # Start with theme_bw
  theme(axis.text = element_text(color = "black", size = 22, face = 'bold'), # Set axis text color to black
        axis.title = element_text(color = "black", size = 24, face = 'bold'), # Set axis title color to black
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        legend.position = c(0.83, 0.3),# Place legend inside the plot at the bottom right
        panel.border = element_rect(colour = "black", fill=NA, size=3.5), # Plot border
        legend.text = element_text(size= 22, color = 'black', face = 'bold'), # Legend text
        legend.title = element_text(size= 22, color = 'black', face = 'bold'),
        legend.key.size = unit(0.65, "cm"),
        legend.spacing.x = unit(0.5, "cm"),  # Horizontal spacing
        legend.spacing.y = unit(0.5, "cm"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt"))) + xlim(0, 5)

# Print the plot
print(plot)



# Assuming downregulated_gos is already defined and contains "Down-regulated" genes

plot <- ggplot(downregulated_gos, aes(y= Term, x = -log10(P.value), fill = Overlap)) + 
  geom_bar(stat = 'identity', width = .7, color = "black", linewidth = 1.05) + # Add black borders to bars
  scale_fill_gradientn(
    colors = c('#B3E5FC', '#311B92'),  # Custom color gradient
    name = "Counts", 
    limits = range(downregulated_gos$Overlap, na.rm = TRUE),
    breaks = c(2, 3, 4),  # Specify breaks at 5, 6, and 8
    labels = c("2", "3", '4')  # Corresponding labels for these breaks
  ) +
  labs(x = "-log10(qvalue)", y = "") + 
  theme_bw() + # Start with theme_bw
  theme(axis.text = element_text(color = "black", size = 22, face = 'bold'), # Set axis text color to black
        axis.title = element_text(color = "black", size = 24, face = 'bold'), # Set axis title color to black
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        legend.position = c(0.83, 0.3),# Place legend inside the plot at the bottom right
        panel.border = element_rect(colour = "black", fill=NA, size=3.5), # Plot border
        legend.text = element_text(size= 22, color = 'black', face = 'bold'), # Legend text
        legend.title = element_text(size= 22, color = 'black', face = 'bold'),
        legend.key.size = unit(0.65, "cm"),
        legend.spacing.x = unit(0.5, "cm"),  # Horizontal spacing
        legend.spacing.y = unit(0.5, "cm"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt"))) + xlim(0, 5)

# Print the plot
print(plot)



# Shared DEGs between A366 and NO_A366

shared <- results_A366[results_A366$gene %in% results_NO_A366$gene, ]

#volcano plot
genes_to_label <- c('Kif15', 'Kif14', 'Cenpe', 'Cep70', 'Mki67', 'Cenpa', 'Top2a','Cdca2', 'Ppara', 'Ehd4', 'Ddc')
label_data <- shared[shared$gene %in% genes_to_label, ]

# Add a column to classify genes as up or downregulated
shared$regulation <- ifelse(shared$logFC > 0.25, "Upregulated",
                            ifelse(shared$logFC < -0.25, "Downregulated", "Neutral"))


# Create the plot
volcano_plot <- ggplot(shared, aes(x = logFC, y = -log10(FDR), color = regulation)) +
  geom_point(size = 4) +
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
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, colour = 'black', face = 'bold'), # Adjust axis text size
    legend.title = element_blank(),
    legend.position = "none" # Adjust legend position if necessary
  ) +
  xlim(-1.5, 1.5) + # Set x-axis limits
  ylim(0, 80) +
  guides(color = guide_legend(title = "Regulation")) # Adjust legend title

# Print the plot
print(volcano_plot)

# GO

dbs <- c("GO_Cellular_Component_2021","GO_Biological_Process_2021" , "GO_Molecular_Function_2021","KEGG_2018")
up <- subset(shared, logFC > 0.25)
down <- subset(shared, logFC < -0.25)

up <- up$gene 
down <- down$gene

up <- enrichr(up, dbs)
down <- enrichr(down, dbs)

up <- up$GO_Biological_Process_2021
down <- down$GO_Biological_Process_2021


up$type <- "Up-regulated"
down$type <- "Down-regulated"

up <- up[order(up$Adjusted.P.value), ]  # sort
up <- subset(up, up$Term %in% c('mitotic sister chromatid segregation (GO:0000070)', 'mitotic spindle organization (GO:0007052)', 'positive regulation of cell cycle process (GO:0090068)', 'mitotic cytokinesis (GO:0000281)', 'regulation of cell cycle G2/M phase transition (GO:1902749)'))

down <- down[order(down$Adjusted.P.value), ]  # sort
down <- subset(down, down$Term %in% c('organelle assembly (GO:0070925)', 'positive regulation of fatty acid beta-oxidation (GO:0032000)', 'negative regulation of appetite (GO:0032099)', 'indolalkylamine biosynthetic process (GO:0046219)', 'catechol-containing compound biosynthetic process (GO:0009713)'))

gos <- rbind(up,down)
gos$Term <- factor(gos$Term, levels=gos$Term)


# Filter the data for only "Up-regulated" genes
upregulated_gos <- gos[gos$type == "Up-regulated", ]

# Remove '(GO:' and any subsequent characters from the Term labels
upregulated_gos$Term <- gsub("\\(GO:.*$", "", upregulated_gos$Term)

# Order the terms based on the -log10(qvalue) from high to low
upregulated_gos <- upregulated_gos %>%
  mutate(Term = fct_reorder(Term, -log10(Adjusted.P.value)))


# Convert Overlap to numeric after removing the part after '/'
upregulated_gos$Overlap <- as.numeric(sub("/.*", "", upregulated_gos$Overlap))


plot <- ggplot(upregulated_gos, aes(y= Term, x = -log10(Adjusted.P.value), fill = Overlap)) + 
  geom_bar(stat = 'identity', width = .7, color = "black", linewidth = 1.05) + # Add black borders to bars
  scale_fill_gradientn(
    colors = c('orange', 'Brown'),  # Custom color gradient
    name = "Counts", 
    limits = range(upregulated_gos$Overlap, na.rm = TRUE),
    breaks = c(11, 14, 23),  # Specify breaks
    labels = c("11", "14", "23")  # Corresponding labels for these breaks
  ) +
  labs(x = "-log10(qvalue)", y = "") + 
  theme_bw() + # Start with theme_bw
  theme(axis.text = element_text(color = "black", size = 22, face = 'bold'), # Set axis text color to black
        axis.title = element_text(color = "black", size = 24, face = 'bold'), # Set axis title color to black
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        legend.position = c(0.83, 0.3),# Place legend inside the plot at the bottom right
        panel.border = element_rect(colour = "black", fill=NA, size=3.5), # Plot border
        legend.text = element_text(size= 22, color = 'black', face = 'bold'), # Legend text
        legend.title = element_text(size= 22, color = 'black', face = 'bold'),
        legend.key.size = unit(0.65, "cm"),
        legend.spacing.x = unit(0.5, "cm"),  # Horizontal spacing
        legend.spacing.y = unit(0.5, "cm"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt")))

# Print the plot
print(plot)


# Down_reg

downregulated_gos <- gos[gos$type == "Down-regulated", ]

# Remove '(GO:' and any subsequent characters from the Term labels
downregulated_gos$Term <- gsub("\\(GO:.*$", "", downregulated_gos$Term)

# Order the terms based on the -log10(qvalue) from high to low
downregulated_gos <- downregulated_gos %>%
  mutate(Term = fct_reorder(Term, -log10(P.value)))

# Convert Overlap to numeric after removing the part after '/'
downregulated_gos$Overlap <- as.numeric(sub("/.*", "", downregulated_gos$Overlap))

plot <- ggplot(downregulated_gos, aes(y= Term, x = -log10(P.value), fill = Overlap)) + 
  geom_bar(stat = 'identity', width = .7, color = "black", linewidth = 1.05) + # Add black borders to bars
  scale_fill_gradientn(
    colors = c('#B3E5FC', '#311B92'),  # Custom color gradient
    name = "Counts", 
    limits = range(downregulated_gos$Overlap, na.rm = TRUE),
    breaks = c(1, 2),  # Specify breaks at 5, 6, and 8
    labels = c("1", "2")  # Corresponding labels for these breaks
  ) +
  labs(x = "-log10(qvalue)", y = "") + 
  theme_bw() + # Start with theme_bw
  theme(axis.text = element_text(color = "black", size = 22, face = 'bold'), # Set axis text color to black
        axis.title = element_text(color = "black", size = 24, face = 'bold'), # Set axis title color to black
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        legend.position = c(0.83, 0.3),# Place legend inside the plot at the bottom right
        panel.border = element_rect(colour = "black", fill=NA, size=3.5), # Plot border
        legend.text = element_text(size= 22, color = 'black', face = 'bold'), # Legend text
        legend.title = element_text(size= 22, color = 'black', face = 'bold'),
        legend.key.size = unit(0.65, "cm"),
        legend.spacing.x = unit(0.5, "cm"),  # Horizontal spacing
        legend.spacing.y = unit(0.5, "cm"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt"))) + xlim(0, 4)

# Print the plot
print(plot)









