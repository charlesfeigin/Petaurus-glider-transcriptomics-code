# Load DESeq2
library(magrittr)
#library(apeglm)
library(DESeq2)

# Read in feature counts table
#read.counts <- read.table("featureCounts_results_prelim1_all.counts", header = TRUE)
read.counts <- read.table("bat_read_counts.txt", header = TRUE)

# Set the row names to gene IDs.
row.names(read.counts) <- read.counts$Geneid

# Remove extraneous rows before counts
read.counts <- read.counts[, -c(1:6)]

# Extract each tissue type from the sample name (Pat = patagium, Dors = dorsal skin)
tissueType <- gsub("[0-9].*", "", colnames(read.counts))

# Make a data frame of sample info
sample_info <- data.frame(sample = colnames(read.counts), tissue = tissueType)

row.names(sample_info) = colnames(read.counts)

# Generate DESeq2 datasets
dds_pat_vs_dors <- DESeqDataSetFromMatrix(countData = read.counts, colData = sample_info, design = ~ tissue)

## Remove genes that have effectively zero read counts
dds_pat_vs_dors <- dds_pat_vs_dors[rowMeans(counts(dds_pat_vs_dors)) >= 10, ]

# Relevel deseq datasets to set comparison tissue as reference level
dds_pat_vs_dors$tissue <- relevel(dds_pat_vs_dors$tissue, ref = "Dors")

## Run DESeq2 analysis
dds_pat_vs_dors <- DESeq(dds_pat_vs_dors, minReplicatesForReplace = 3)

# Collect the results of the DESeq2 analysis
dds_pat_vs_dors.results <- results(dds_pat_vs_dors, independentFiltering = TRUE, alpha = 0.01)

# Write them to a table
write.table(as.data.frame(dds_pat_vs_dors.results), file = "dds_pat_vs_dors.results.txt", col.names = TRUE, quote = FALSE, sep = "\t")

# Collect differentially expressed genes
dds_pat_vs_dors.results.degs <- subset(dds_pat_vs_dors.results, dds_pat_vs_dors.results$padj < 0.01 & abs(dds_pat_vs_dors.results$log2FoldChange) >= 0.584962501)

# Write them to a table
write.table(as.data.frame(dds_pat_vs_dors.results.degs), file = "dds_pat_vs_dors_1.5fold.degs.txt", col.names = TRUE, quote = FALSE, sep = "\t")



### Volcano plots ###

library("ggplot2")
library("ggrepel")
library(grid)
library(gridExtra)
#library("ggfortify")
#library("ggbiplot")
library("extrafont")
library("dplyr")
library("tibble")
library("ggthemes")
font_import()
loadfonts()



### Candidate lists are the set of genes differentially expressed in the given pairwise DGE comparison.
### and belonging to enriched limb-morphogenesis GO terms


candidates <- c("HOXD11",
"WNT5A",
"FGF10",
"GRHL2",
"EXT1",
"TBX5",
"TP63",
"HAND2",
"HOXD9",
"GREM1",
"OSR1",
"HOXD13",
"MAP3K20",
"TBX3")


# dds_stg2_pat_vs_dls LIMB GENES HIGHLIGHTED
volcanoPlotData <- as.data.frame(dds_pat_vs_dors.results)
volcanoPlotData <- volcanoPlotData[volcanoPlotData$baseMean < 25000, ]
volcanoPlotData <- volcanoPlotData %>% rownames_to_column('geneName')
volcanoPlotData <- volcanoPlotData %>% mutate(candidate = (geneName %in% candidates))
#write.table(as.data.frame(volcanoPlotData), file = "volcanoPlotData.txt", col.names = TRUE, quote = FALSE, sep = "\t")
candidateSubset <- volcanoPlotData %>% filter (candidate == TRUE)
candidateSubset <- candidateSubset %>% filter (padj < .01)
upreg <- candidateSubset %>% filter (log2FoldChange >= 0.584962501)
downreg <- candidateSubset %>% filter (log2FoldChange <= -0.584962501)
upreg <- candidateSubset
# Sanity check
#volcanoPlotData %>% filter (candidate == TRUE)

myplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean)) +
    geom_point(color = "grey48", alpha = 0.75, size = 1.5) +
	geom_text_repel(data = downreg, color = "black", ylim = c(50,Inf), xlim = c(-4.75,-4.75), size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", downreg$geneName))) +
	geom_text_repel(data = upreg, color = "black", ylim = c(20,Inf), xlim = c(4,7.75), force = 128, size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", upreg$geneName))) +
	geom_point(data = upreg, color = "darkgreen", alpha = 1, size = 2.75) + 
	geom_point(data = downreg, color = "purple", alpha = 1, size = 2.75) +
	xlim(-5,8) +
	xlab(paste0("Log 2-fold change")) +
	ylab(paste0("Mean expression")) +
	geom_rangeframe() +
    theme(text=element_text(family="Arial", size=7), legend.position = "none") +
	theme_classic()

ggsave("dds_bat_pat_vs_dors_LIMB.png",
	width = 8.76,
	height = 6.62,
	myplot,
	dpi = 600)
	
	
ggsave("dds_bat_pat_vs_dors_LIMB.svg",
	width = 8.76,
	height = 6.62,
	myplot,
	dpi = 600)






































# dds_stg1_pat_vs_shd
volcanoPlotData <- as.data.frame(dds_pat_vs_dors.results)
volcanoPlotData <- volcanoPlotData[volcanoPlotData$baseMean < 2500, ]
volcanoPlotData <- volcanoPlotData %>% rownames_to_column('geneName')
volcanoPlotData <- volcanoPlotData %>% mutate(candidate = (geneName %in% candidates))
candidateSubset <- volcanoPlotData %>% filter (candidate == TRUE)
candidateSubset <- candidateSubset %>% filter (padj < .01)
upreg <- candidateSubset %>% filter (log2FoldChange > 0)
downreg <- candidateSubset %>% filter (log2FoldChange < 0)

# Sanity check
#volcanoPlotData %>% filter (candidate == TRUE)

myplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean)) +
    geom_point(color = "grey48", alpha = 0.75, size = 1.5) +
	geom_text_repel(data = downreg, color = "black", ylim = c(50,Inf), xlim = c(-4.75,-4.75), force = 1, size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", downreg$geneName))) +
	geom_text_repel(data = upreg, color = "black", ylim = c(20,Inf), xlim = c(4,7.75), force = 16, size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", upreg$geneName))) +
	geom_point(data = upreg, color = "darkgreen", alpha = 1, size = 2.75) + 
	geom_point(data = downreg, color = "purple", alpha = 1, size = 2.75) +
	xlab(paste0("Log 2-fold change")) +
	ylab(paste0("Mean expression")) +
	geom_rangeframe() +
    theme(text=element_text(family="Arial", size=7), legend.position = "none") +
	theme_classic()

ggsave("dds_stg1_pat_vs_shd.png",
	width = 8.76,
	height = 6.62,
	myplot,
	dpi = 600)
	
ggsave("dds_stg1_pat_vs_shd.svg",
	width = 8.76,
	height = 6.62,
	myplot,
	dpi = 600)



# dds_stg2_pat_vs_shd
volcanoPlotData <- as.data.frame(dds_stg2_pat_vs_shd.results)
volcanoPlotData <- volcanoPlotData[volcanoPlotData$baseMean < 2500, ]
volcanoPlotData <- volcanoPlotData %>% rownames_to_column('geneName')
volcanoPlotData <- volcanoPlotData %>% mutate(candidate = (geneName %in% candidates))
candidateSubset <- volcanoPlotData %>% filter (candidate == TRUE)
candidateSubset <- candidateSubset %>% filter (padj < .01)
upreg <- candidateSubset %>% filter (log2FoldChange > 0)
downreg <- candidateSubset %>% filter (log2FoldChange < 0)

# Sanity check
#volcanoPlotData %>% filter (candidate == TRUE)

myplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean)) +
    geom_point(color = "grey48", alpha = 0.75, size = 1.5) +
	geom_text_repel(data = downreg, color = "black", ylim = c(50,Inf), xlim = c(-4,-4), size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", downreg$geneName))) +
	geom_text_repel(data = upreg, color = "black", ylim = c(20,Inf), xlim = c(4,7.75), force = 16, size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", upreg$geneName))) +
	geom_point(data = upreg, color = "darkgreen", alpha = 1, size = 2.75) + 
	geom_point(data = downreg, color = "purple", alpha = 1, size = 2.75) +
	xlab(paste0("Log 2-fold change")) +
	ylab(paste0("Mean expression")) +
	geom_rangeframe() +
    theme(text=element_text(family="Arial", size=7), legend.position = "none") +
	theme_classic()

ggsave("dds_stg2_pat_vs_shd.png",
	width = 8.76,
	height = 6.62,
	myplot,
	dpi = 600)
	


ggsave("dds_stg2_pat_vs_shd.svg",
	width = 8.76,
	height = 6.62,
	myplot,
	dpi = 600)

# dds_stg3_pat_vs_shd
#volcanoPlotData <- as.data.frame(dds_stg3_pat_vs_shd.results)
#volcanoPlotData <- volcanoPlotData[volcanoPlotData$baseMean < 2500, ]
#volcanoPlotData <- volcanoPlotData %>% rownames_to_column('geneName')
#volcanoPlotData <- volcanoPlotData %>% mutate(candidate = (geneName %in% candidates))
#candidateSubset <- volcanoPlotData %>% filter (candidate == TRUE)
#candidateSubset <- candidateSubset %>% filter (padj < .01)
#upreg <- candidateSubset %>% filter (log2FoldChange > 0)
#downreg <- candidateSubset %>% filter (log2FoldChange < 0)

# Sanity check
#volcanoPlotData %>% filter (candidate == TRUE)

#myplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean)) +
#    geom_point(color = "grey48", alpha = 0.75, size = 1.5) +
#	geom_text_repel(data = downreg, color = "black", ylim = c(50,Inf), xlim = c(-4.75,-4.75), size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", downreg$geneName))) +
#	geom_text_repel(data = upreg, color = "black", ylim = c(20,Inf), xlim = c(4,7), force = 24, size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", upreg$geneName))) +
#	geom_point(data = upreg, color = "darkgreen", alpha = 1, size = 2.75) + 
#	geom_point(data = downreg, color = "purple", alpha = 1, size = 2.75) +
#	xlab(paste0("Log 2-fold change")) +
#	ylab(paste0("Mean expression")) +
#	geom_rangeframe() +
#   theme(text=element_text(family="Arial", size=7), legend.position = "none") +
#	theme_classic()


#ggsave("dds_stg3_pat_vs_shd.svg",
#	width = 8.76,
#	height = 6.62,
#	myplot,
#	dpi = 600)



# dds_stg1_pat_vs_dls
volcanoPlotData <- as.data.frame(dds_stg1_pat_vs_dls.results)
volcanoPlotData <- volcanoPlotData[volcanoPlotData$baseMean < 2500, ]
volcanoPlotData <- volcanoPlotData %>% rownames_to_column('geneName')
volcanoPlotData <- volcanoPlotData %>% mutate(candidate = (geneName %in% candidates))
candidateSubset <- volcanoPlotData %>% filter (candidate == TRUE)
candidateSubset <- candidateSubset %>% filter (padj < .01)
upreg <- candidateSubset %>% filter (log2FoldChange > 0)
downreg <- candidateSubset %>% filter (log2FoldChange < 0)

# Sanity check
#volcanoPlotData %>% filter (candidate == TRUE)

myplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean)) +
    geom_point(color = "grey48", alpha = 0.75, size = 1.5) +
	geom_text_repel(data = downreg, color = "black", ylim = c(50,Inf), xlim = c(-4.75,-4.75), size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", downreg$geneName))) +
	geom_text_repel(data = upreg, color = "black", ylim = c(20,Inf), xlim = c(4,7.75), force = 16, size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", upreg$geneName))) +
	geom_point(data = upreg, color = "darkgreen", alpha = 1, size = 2.75) + 
	geom_point(data = downreg, color = "purple", alpha = 1, size = 2.75) +
	xlab(paste0("Log 2-fold change")) +
	ylab(paste0("Mean expression")) +
	geom_rangeframe() +
    theme(text=element_text(family="Arial", size=7), legend.position = "none") +
	theme_classic()

ggsave("dds_stg1_pat_vs_dls.png",
	width = 8.76,
	height = 6.62,
	myplot,
	dpi = 600)
	
	
ggsave("dds_stg1_pat_vs_dls.svg",
	width = 8.76,
	height = 6.62,
	myplot,
	dpi = 600)




# dds_stg2_pat_vs_dls
volcanoPlotData <- as.data.frame(dds_stg2_pat_vs_dls.results)
volcanoPlotData <- volcanoPlotData[volcanoPlotData$baseMean < 2500, ]
volcanoPlotData <- volcanoPlotData %>% rownames_to_column('geneName')
volcanoPlotData <- volcanoPlotData %>% mutate(candidate = (geneName %in% candidates))
candidateSubset <- volcanoPlotData %>% filter (candidate == TRUE)
candidateSubset <- candidateSubset %>% filter (padj < .01)
upreg <- candidateSubset %>% filter (log2FoldChange > 0)
downreg <- candidateSubset %>% filter (log2FoldChange < 0)

# Sanity check
#volcanoPlotData %>% filter (candidate == TRUE)

myplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean)) +
    geom_point(color = "grey48", alpha = 0.75, size = 1.5) +
	geom_text_repel(data = downreg, color = "black", ylim = c(50,Inf), xlim = c(-4.75,-4.75), size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", downreg$geneName))) +
	geom_text_repel(data = upreg, color = "black", ylim = c(20,Inf), xlim = c(4,7.75), force = 16, size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", upreg$geneName))) +
	geom_point(data = upreg, color = "darkgreen", alpha = 1, size = 2.75) + 
	geom_point(data = downreg, color = "purple", alpha = 1, size = 2.75) +
	xlab(paste0("Log 2-fold change")) +
	ylab(paste0("Mean expression")) +
	geom_rangeframe() +
    theme(text=element_text(family="Arial", size=7), legend.position = "none") +
	theme_classic()

ggsave("dds_stg2_pat_vs_dls.png",
	width = 8.76,
	height = 6.62,
	myplot,
	dpi = 600)
	
	
ggsave("dds_stg2_pat_vs_dls.svg",
	width = 8.76,
	height = 6.62,
	myplot,
	dpi = 600)











candidates <- c("EBF3|XM_020978347.1|MultiHit|fragGeneModel-490-gene9-frag1",
"OSR1|XM_020967455.1|SingleHit|maker-559-snap-gene-8.20",
"MSC|XM_020971030.1|SingleHit|maker-38-augustus-gene-278.2",
"TBX3|XM_021007970.1|SingleHit|maker-99-snap-gene-60.8",
"EMX2|XM_020988183.1|SingleHit|maker-37-augustus-gene-315.2",
"GLIS3|XM_021006269.1|MultiHit|fragGeneModel-187-gene20-frag1",
"IRX6|XM_020984923.1|SingleHit|maker-171-snap-gene-6.16",
"ZFHX3|XM_020963589.1|MultiHit|fragGeneModel-176-gene2-frag1",
"PAX1|XM_020986943.1|MultiHit|fragGeneModel-90919-gene24-frag1",
"BARX1|XM_020980965.1|SingleHit|maker-174-augustus-gene-42.12")






























##########################




myplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean)) +
    geom_point(color = "grey48", alpha = 0.75, size = 1.5) +
	geom_text_repel(data = downreg, color = "black", ylim = c(50,Inf), xlim = c(-4.75,-4.75), size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", downreg$geneName))) +
	geom_text_repel(data = upreg, color = "black", ylim = c(20,Inf), xlim = c(6,7.75), force = 16, size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", upreg$geneName))) +
	geom_point(data = upreg, color = "darkgreen", alpha = 0.8, size = 2.75) + 
	geom_point(data = downreg, color = "purple", alpha = 0.8, size = 2.75) +
	xlab(paste0("Log 2-fold change")) +
	ylab(paste0("Mean expression")) +
	geom_rangeframe() +
    theme(text=element_text(family="Arial", size=7), legend.position = "none") +
	theme_classic()



ggsave("test.png",
	width = 8.76,
	height = 6.62,
	myplot,
	dpi = 2400)



##########################
volcanoPlotData <- as.data.frame(dds_stg1_pat_vs_dls.results.lfcshrink)
volcanoPlotData$candidate <- ifelse(rownames(volcanoPlotData) %in% candidates & rownames(volcanoPlotData) %in% rownames(dds_stg1_pat_vs_dls.results.degs), "red", "black")
volcanoPlotData$degs <- "nondeg"

for (row in 1:nrow(volcanoPlotData)) {
	if (volcanoPlotData[row,]$padj < .01 & volcanoPlotData[row,]$log2FoldChange > 0) {
	volcanoPlotData[row,]$degs <- "upreg"
	} else if (volcanoPlotData[row,]$padj < .01 & volcanoPlotData[row,]$log2FoldChange < 0) {
	volcanoPlotData[row,]$degs <- "downreg"
	}
}

volcanoPlotData <- volcanoPlotData[volcanoPlotData$baseMean < 10000, ]
#volcanoPlotData <- volcanoPlotData[volcanoPlotData$padj != "NA",] 
candidateLabels <- subset(volcanoPlotData, rownames(volcanoPlotData) %in% candidates & rownames(volcanoPlotData) %in% rownames(dds_stg1_pat_vs_dls.results.degs))
#volcano.plot <- 
testplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean, color = degs)) +
			geom_point(size = 1.5, alpha = 0.5) + labs(color = "Significance") + xlab(paste0("Log 2-fold change")) + ylab(paste0("Mean expression")) +
			geom_point(data = candidateLabels, aes(color = "Red"), size = 2.0) +
 			geom_text_repel(data = candidateLabels, ylim = c(1000,Inf), color = "red", size=3, force = 32, segment.size = 0.5, segment.alpha = 1, segment.color = "Red", aes(label=gsub("\\|.*", "", rownames(candidateLabels))))
testplot <- testplot + scale_colour_manual(values = c("Blue", "Black", "Red", "darkgoldenrod3")) + theme(panel.background = element_blank())
ggsave("stg1_pat_vs_dls_volcano_labeled.png",
	width =8.76,
	height =6.62,
	testplot,
	dpi = 2400)
		
volcanoPlotData <- volcanoPlotData[volcanoPlotData$baseMean < 10000, ]
#volcanoPlotData <- volcanoPlotData[volcanoPlotData$padj != "NA",] 
candidateLabels <- subset(volcanoPlotData, rownames(volcanoPlotData) %in% candidates & rownames(volcanoPlotData) %in% rownames(dds_stg1_pat_vs_dls.results.degs))
#volcano.plot <- 
testplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean, color = degs)) +
			geom_point(size = 1.5, alpha = 0.5) + labs(color = "Significance") + xlab(paste0("Log 2-fold change")) + ylab(paste0("Mean expression")) +
			geom_point(data = candidateLabels, aes(color = "Red"), size = 2.0)
testplot <- testplot + scale_colour_manual(values = c("Blue", "Black", "Red", "darkgoldenrod3")) + theme(panel.background = element_blank())
ggsave("stg1_pat_vs_dls_volcano_unlabeled.png",
	width =8.76,
	height =6.62,
	testplot,
	dpi = 2400)


##########################
volcanoPlotData <- as.data.frame(dds_stg2_pat_vs_dls.results.lfcshrink)
volcanoPlotData$candidate <- ifelse(rownames(volcanoPlotData) %in% candidates & rownames(volcanoPlotData) %in% rownames(dds_stg2_pat_vs_dls.results.degs), "red", "black")
volcanoPlotData$degs <- "nondeg"

for (row in 1:nrow(volcanoPlotData)) {
	if (volcanoPlotData[row,]$padj < .01 & volcanoPlotData[row,]$log2FoldChange > 0) {
	volcanoPlotData[row,]$degs <- "upreg"
	} else if (volcanoPlotData[row,]$padj < .01 & volcanoPlotData[row,]$log2FoldChange < 0) {
	volcanoPlotData[row,]$degs <- "downreg"
	}
}



volcanoPlotData <- volcanoPlotData[volcanoPlotData$baseMean < 10000, ]
#volcanoPlotData <- volcanoPlotData[volcanoPlotData$padj != "NA",] 
candidateLabels <- subset(volcanoPlotData, rownames(volcanoPlotData) %in% candidates & rownames(volcanoPlotData) %in% rownames(dds_stg2_pat_vs_dls.results.degs))
#volcano.plot <- 
testplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean, color = degs)) +
			geom_point(size = 1.5, alpha = 0.5) + labs(color = "Significance") + xlab(paste0("Log 2-fold change")) + ylab(paste0("Mean expression")) +
			geom_point(data = candidateLabels, aes(color = "Red"), size = 2.0) +
 			geom_text_repel(data = candidateLabels, ylim = c(1000,Inf), color = "red", size=3, force = 32, segment.size = 0.5, segment.alpha = 1, segment.color = "Red", aes(label=gsub("\\|.*", "", rownames(candidateLabels))))
testplot <- testplot + scale_colour_manual(values = c("Blue", "Black", "Red", "darkgoldenrod3")) + theme(panel.background = element_blank())
ggsave("stg2_pat_vs_dls_volcano_labeled.png",
	width =8.76,
	height =6.62,
	testplot,
	dpi = 2400)
		
volcanoPlotData <- volcanoPlotData[volcanoPlotData$baseMean < 10000, ]
#volcanoPlotData <- volcanoPlotData[volcanoPlotData$padj != "NA",] 
candidateLabels <- subset(volcanoPlotData, rownames(volcanoPlotData) %in% candidates & rownames(volcanoPlotData) %in% rownames(dds_stg2_pat_vs_dls.results.degs))
#volcano.plot <- 
testplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean, color = degs)) +
			geom_point(size = 1.5, alpha = 0.5) + labs(color = "Significance") + xlab(paste0("Log 2-fold change")) + ylab(paste0("Mean expression")) +
			geom_point(data = candidateLabels, aes(color = "Red"), size = 2.0)
testplot <- testplot + scale_colour_manual(values = c("Blue", "Black", "Red", "darkgoldenrod3")) + theme(panel.background = element_blank())
ggsave("stg2_pat_vs_dls_volcano_unlabeled.png",
	width =8.76,
	height =6.62,
	testplot,
	dpi = 2400)
		
		

##########################
volcanoPlotData <- as.data.frame(dds_stg1_pat_vs_shd.results.lfcshrink)
volcanoPlotData$candidate <- ifelse(rownames(volcanoPlotData) %in% candidates & rownames(volcanoPlotData) %in% rownames(dds_stg1_pat_vs_shd.results.degs), "red", "black")
volcanoPlotData$degs <- "nondeg"

for (row in 1:nrow(volcanoPlotData)) {
	if (volcanoPlotData[row,]$padj < .01 & volcanoPlotData[row,]$log2FoldChange > 0) {
	volcanoPlotData[row,]$degs <- "upreg"
	} else if (volcanoPlotData[row,]$padj < .01 & volcanoPlotData[row,]$log2FoldChange < 0) {
	volcanoPlotData[row,]$degs <- "downreg"
	}
}

volcanoPlotData <- volcanoPlotData[volcanoPlotData$baseMean < 10000, ]
#volcanoPlotData <- volcanoPlotData[volcanoPlotData$padj != "NA",] 
candidateLabels <- subset(volcanoPlotData, rownames(volcanoPlotData) %in% candidates & rownames(volcanoPlotData) %in% rownames(dds_stg1_pat_vs_shd.results.degs))
#volcano.plot <- 
testplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean, color = degs)) +
			geom_point(size = 1.5, alpha = 0.5) + labs(color = "Significance") + xlab(paste0("Log 2-fold change")) + ylab(paste0("Mean expression")) +
			geom_point(data = candidateLabels, aes(color = "Red"), size = 2.0) +
 			geom_text_repel(data = candidateLabels, ylim = c(1000,Inf), color = "red", size=3, force = 32, segment.size = 0.5, segment.alpha = 1, segment.color = "Red", aes(label=gsub("\\|.*", "", rownames(candidateLabels))))
testplot <- testplot + scale_colour_manual(values = c("Blue", "Black", "Red", "darkgoldenrod3")) + theme(panel.background = element_blank())
ggsave("stg1_pat_vs_shd_volcano_labeled.png",
	width =8.76,
	height =6.62,
	testplot,
	dpi = 2400)
		
volcanoPlotData <- volcanoPlotData[volcanoPlotData$baseMean < 10000, ]
#volcanoPlotData <- volcanoPlotData[volcanoPlotData$padj != "NA",] 
candidateLabels <- subset(volcanoPlotData, rownames(volcanoPlotData) %in% candidates & rownames(volcanoPlotData) %in% rownames(dds_stg1_pat_vs_shd.results.degs))
#volcano.plot <- 
testplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean, color = degs)) +
			geom_point(size = 1.5, alpha = 0.5) + labs(color = "Significance") + xlab(paste0("Log 2-fold change")) + ylab(paste0("Mean expression")) +
			geom_point(data = candidateLabels, aes(color = "Red"), size = 2.0)
testplot <- testplot + scale_colour_manual(values = c("Blue", "Black", "Red", "darkgoldenrod3")) + theme(panel.background = element_blank())
ggsave("stg1_pat_vs_shd_volcano_unlabeled.png",
	width =8.76,
	height =6.62,
	testplot,
	dpi = 2400)
		
		
##########################
volcanoPlotData <- as.data.frame(dds_stg2_pat_vs_shd.results.lfcshrink)
volcanoPlotData$candidate <- ifelse(rownames(volcanoPlotData) %in% candidates & rownames(volcanoPlotData) %in% rownames(dds_stg2_pat_vs_shd.results.degs), "red", "black")
volcanoPlotData$degs <- "nondeg"

for (row in 1:nrow(volcanoPlotData)) {
	if (volcanoPlotData[row,]$padj < .01 & volcanoPlotData[row,]$log2FoldChange > 0) {
	volcanoPlotData[row,]$degs <- "upreg"
	} else if (volcanoPlotData[row,]$padj < .01 & volcanoPlotData[row,]$log2FoldChange < 0) {
	volcanoPlotData[row,]$degs <- "downreg"
	}
}

volcanoPlotData <- volcanoPlotData[volcanoPlotData$baseMean < 10000, ]
#volcanoPlotData <- volcanoPlotData[volcanoPlotData$padj != "NA",] 
candidateLabels <- subset(volcanoPlotData, rownames(volcanoPlotData) %in% candidates & rownames(volcanoPlotData) %in% rownames(dds_stg2_pat_vs_shd.results.degs))
#volcano.plot <- 
testplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean, color = degs)) +
			geom_point(size = 1.5, alpha = 0.5) + labs(color = "Significance") + xlab(paste0("Log 2-fold change")) + ylab(paste0("Mean expression")) +
			geom_point(data = candidateLabels, aes(color = "Red"), size = 2.0) +
 			geom_text_repel(data = candidateLabels, ylim = c(1000,Inf), color = "red", size=3, force = 32, segment.size = 0.5, segment.alpha = 1, segment.color = "Red", aes(label=gsub("\\|.*", "", rownames(candidateLabels))))
testplot <- testplot + scale_colour_manual(values = c("Blue", "Black", "Red", "darkgoldenrod3")) + theme(panel.background = element_blank())
ggsave("stg2_pat_vs_shd_volcano_labeled.png",
	width =8.76,
	height =6.62,
	testplot,
	dpi = 2400)
		
volcanoPlotData <- volcanoPlotData[volcanoPlotData$baseMean < 10000, ]
#volcanoPlotData <- volcanoPlotData[volcanoPlotData$padj != "NA",] 
candidateLabels <- subset(volcanoPlotData, rownames(volcanoPlotData) %in% candidates & rownames(volcanoPlotData) %in% rownames(dds_stg2_pat_vs_shd.results.degs))
#volcano.plot <- 
testplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean, color = degs)) +
			geom_point(size = 1.5, alpha = 0.5) + labs(color = "Significance") + xlab(paste0("Log 2-fold change")) + ylab(paste0("Mean expression")) +
			geom_point(data = candidateLabels, aes(color = "Red"), size = 2.0)
testplot <- testplot + scale_colour_manual(values = c("Blue", "Black", "Red", "darkgoldenrod3")) + theme(panel.background = element_blank())
ggsave("stg2_pat_vs_shd_volcano_unlabeled.png",
	width =8.76,
	height =6.62,
	testplot,
	dpi = 2400)
		
		
##########################
volcanoPlotData <- as.data.frame(dds_stg3_pat_vs_shd.results.lfcshrink)
volcanoPlotData$candidate <- ifelse(rownames(volcanoPlotData) %in% candidates & rownames(volcanoPlotData) %in% rownames(dds_stg3_pat_vs_shd.results.degs), "red", "black")
volcanoPlotData$degs <- "nondeg"

for (row in 1:nrow(volcanoPlotData)) {
	if (volcanoPlotData[row,]$padj < .01 & volcanoPlotData[row,]$log2FoldChange > 0) {
	volcanoPlotData[row,]$degs <- "upreg"
	} else if (volcanoPlotData[row,]$padj < .01 & volcanoPlotData[row,]$log2FoldChange < 0) {
	volcanoPlotData[row,]$degs <- "downreg"
	}
}

volcanoPlotData <- volcanoPlotData[volcanoPlotData$baseMean < 10000, ]
#volcanoPlotData <- volcanoPlotData[volcanoPlotData$padj != "NA",] 
candidateLabels <- subset(volcanoPlotData, rownames(volcanoPlotData) %in% candidates & rownames(volcanoPlotData) %in% rownames(dds_stg3_pat_vs_shd.results.degs))
#volcano.plot <- 
testplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean, color = degs)) +
			geom_point(size = 1.5, alpha = 0.5) + labs(color = "Significance") + xlab(paste0("Log 2-fold change")) + ylab(paste0("Mean expression")) +
			geom_point(data = candidateLabels, aes(color = "Red"), size = 2.0) +
 			geom_text_repel(data = candidateLabels, ylim = c(1000,Inf), color = "red", size=3, force = 32, segment.size = 0.5, segment.alpha = 1, segment.color = "Red", aes(label=gsub("\\|.*", "", rownames(candidateLabels))))
testplot <- testplot + scale_colour_manual(values = c("Blue", "Black", "Red", "darkgoldenrod3")) + theme(panel.background = element_blank())
ggsave("stg3_pat_vs_shd_volcano_labeled.png",
	width =8.76,
	height =6.62,
	testplot,
	dpi = 2400)
		
volcanoPlotData <- volcanoPlotData[volcanoPlotData$baseMean < 10000, ]
#volcanoPlotData <- volcanoPlotData[volcanoPlotData$padj != "NA",] 
candidateLabels <- subset(volcanoPlotData, rownames(volcanoPlotData) %in% candidates & rownames(volcanoPlotData) %in% rownames(dds_stg3_pat_vs_shd.results.degs))
#volcano.plot <- 
testplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean, color = degs)) +
			geom_point(size = 1.5, alpha = 0.5) + labs(color = "Significance") + xlab(paste0("Log 2-fold change")) + ylab(paste0("Mean expression")) +
			geom_point(data = candidateLabels, aes(color = "Red"), size = 2.0)
testplot <- testplot + scale_colour_manual(values = c("Blue", "Black", "Red", "darkgoldenrod3")) + theme(panel.background = element_blank())
ggsave("stg3_pat_vs_shd_volcano_unlabeled.png",
	width =8.76,
	height =6.62,
	testplot,
	dpi = 2400)		






















# PCA


vst <- dds_pat_vs_dors.vst.bl
intgroup <- "tissue"
ntop <- 200

	rv <- rowVars(assay(dds_pat_vs_dors.vst.bl))

	select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
	
	df2 <- as.data.frame(t(assay(dds_pat_vs_dors.vst.bl)[select, ]))
	
	pca_vst <- prcomp(df2, center=TRUE,scale=TRUE)
	
	pca_df <- as.data.frame(pca_vst$x)
	
	combined <- pca_df
	combined$sample <- sample_info$sample
	combined$tissue <- sample_info$tissue
	combined$stage <- sample_info$stage
	combined$weight <- sample_info$weight
	
	PC1v2 <- subset(combined, select=c(PC1, PC2, sample, tissue, stage, weight))
	
	png("pat_vs_dors_vst_Blind_200_PC1_vs_PC2.png")
	p <- ggplot(PC1v2, aes(x=PC1, y=PC2, color = "tissue"))
	p <- p + geom_point( size = 3, aes(shape=factor(sample_info$stage), colour=factor(sample_info$tissue)))
	p
	dev.off()


rld <- dds_pat_vs_dors.rld.bl
intgroup <- "tissue"
ntop <- 200

	rv <- rowVars(assay(dds_pat_vs_dors.rld.bl))

	select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
	
	df2 <- as.data.frame(t(assay(dds_pat_vs_dors.rld.bl)[select, ]))
	
	pca_rld <- prcomp(df2, center=TRUE,scale=TRUE)
	
	pca_df <- as.data.frame(pca_rld$x)
	
	combined <- pca_df
	combined$sample <- sample_info$sample
	combined$tissue <- sample_info$tissue
	combined$stage <- sample_info$stage
	combined$weight <- sample_info$weight
	
	PC1v2 <- subset(combined, select=c(PC1, PC2, sample, tissue, stage, weight))
	
	png("pat_vs_dors_rld_Blind_200_PC1_vs_PC2.png")
	p <- ggplot(PC1v2, aes(x=PC1, y=PC2, color = "tissue"))
	p <- p + geom_point( size = 3, aes(shape=factor(sample_info$stage), colour=factor(sample_info$tissue)))
	p
	dev.off()




