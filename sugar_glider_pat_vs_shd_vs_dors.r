library(magrittr)
library(DESeq2)

# Read in feature counts table
read.counts <- read.table("pbrev_read_counts.txt", header = TRUE)

# Set the row names to gene IDs.
row.names(read.counts) <- read.counts$Geneid

# Remove extraneous rows before counts
read.counts <- read.counts[, -c(1:6)]

# Get rid of "X" character added to names
betterNames <- gsub( "X", "", names(read.counts))

# Update names
names(read.counts) <- betterNames

# Extract each tissue type from the sample name (pat = patagium, dors = dorsal skin, shd = shoulder skin)
tissueType <- gsub("^[0-9]+[A-Z]_", "", betterNames)
tissueType <- gsub("_stg[1-3]+_[0-9]+$", "", tissueType)

# Extract each developmental stage from the sample name (stage 1 is prior to outgrowth, stage 2 is during early outgrowth)
devStage <- gsub("^[0-9]+[A-Z]_[a-z]+_", "", betterNames)
devStage <- gsub("_[0-9]+$", "", devStage)

# Load in animals weights separately
sampleWeight <- read.table("pbrev_sample_weights.txt", header=TRUE)

# Make a data frame of sample info
sample_info <- data.frame(sample = betterNames, tissue = tissueType, stage = devStage, weight = sampleWeight$weight)

# Add matching names to sample_info
row.names(sample_info) = betterNames


# Subset data for different tests.
#pat vs dors
stg1_pat_vs_dors_info <- sample_info[c(c(grep("pat_stg1", row.names(sample_info))), c(grep("dors_stg1", row.names(sample_info)))),]
stg1_pat_vs_dors_counts <- read.counts[c(c(grep("pat_stg1", colnames(read.counts))), c(grep("dors_stg1", colnames(read.counts))))]
stg2_pat_vs_dors_info <- sample_info[c(c(grep("pat_stg2", row.names(sample_info))), c(grep("dors_stg2", row.names(sample_info)))),]
stg2_pat_vs_dors_counts <- read.counts[c(c(grep("pat_stg2", colnames(read.counts))), c(grep("dors_stg2", colnames(read.counts))))]
#pat vs shd
stg1_pat_vs_shd_info <- sample_info[c(c(grep("pat_stg1", row.names(sample_info))), c(grep("shd_stg1", row.names(sample_info)))),]
stg1_pat_vs_shd_counts <- read.counts[c(c(grep("pat_stg1", colnames(read.counts))), c(grep("shd_stg1", colnames(read.counts))))]

## Generate DESeq2 datasets
#pat vs dors
dds_stg1_pat_vs_dors <- DESeqDataSetFromMatrix(countData = stg1_pat_vs_dors_counts, colData = stg1_pat_vs_dors_info, design = ~ tissue)
dds_stg2_pat_vs_dors <- DESeqDataSetFromMatrix(countData = stg2_pat_vs_dors_counts, colData = stg2_pat_vs_dors_info, design = ~ tissue)
#pat vs shd
dds_stg1_pat_vs_shd <- DESeqDataSetFromMatrix(countData = stg1_pat_vs_shd_counts, colData = stg1_pat_vs_shd_info, design = ~ tissue)

## Remove genes that have effectively zero read counts
#pat vs dors
dds_stg1_pat_vs_dors <- dds_stg1_pat_vs_dors[rowMeans(counts(dds_stg1_pat_vs_dors)) >= 10, ]
dds_stg2_pat_vs_dors <- dds_stg2_pat_vs_dors[rowMeans(counts(dds_stg2_pat_vs_dors)) >= 10, ]
#pat vs shd
dds_stg1_pat_vs_shd <- dds_stg1_pat_vs_shd[rowMeans(counts(dds_stg1_pat_vs_shd)) >= 10, ]


## Relevel deseq datasets to set comparison tissue as reference level
#pat vs dors
dds_stg1_pat_vs_dors$tissue <- relevel(dds_stg1_pat_vs_dors$tissue, ref = "dors")
dds_stg2_pat_vs_dors$tissue <- relevel(dds_stg2_pat_vs_dors$tissue, ref = "dors")
#pat vs shd
dds_stg1_pat_vs_shd$tissue <- relevel(dds_stg1_pat_vs_shd$tissue, ref = "shd")


## Run DESeq2 analysis
#pat vs dors
dds_stg1_pat_vs_dors <- DESeq(dds_stg1_pat_vs_dors, minReplicatesForReplace = 3)
dds_stg2_pat_vs_dors <- DESeq(dds_stg2_pat_vs_dors, minReplicatesForReplace = 3)
#pat vs shd
dds_stg1_pat_vs_shd <- DESeq(dds_stg1_pat_vs_shd, minReplicatesForReplace = 3)

# Collect the results of the DESeq2 analysis
# pat vs dors
dds_stg1_pat_vs_dors.results <- results(dds_stg1_pat_vs_dors, independentFiltering = TRUE, alpha = 0.01)
dds_stg2_pat_vs_dors.results <- results(dds_stg2_pat_vs_dors, independentFiltering = TRUE, alpha = 0.01)
# pat vs shd
dds_stg1_pat_vs_shd.results <- results(dds_stg1_pat_vs_shd, independentFiltering = TRUE, alpha = 0.01)

# Write them to a table
write.table(as.data.frame(dds_stg1_pat_vs_dors.results), file = "dds_stg1_pat_vs_dors.results_5-22.txt", col.names = TRUE, quote = FALSE, sep = "\t")
write.table(as.data.frame(dds_stg2_pat_vs_dors.results), file = "dds_stg2_pat_vs_dors.results_5-22.txt", col.names = TRUE, quote = FALSE, sep = "\t")
write.table(as.data.frame(dds_stg1_pat_vs_shd.results), file = "dds_stg1_pat_vs_shd.results_5-22.txt", col.names = TRUE, quote = FALSE, sep = "\t")

# Collect differentially expressed genes (log2FoldChange of 0.584962501 corresponds to fold change of ~1.5).
dds_stg1_pat_vs_dors.results.degs <- subset(dds_stg1_pat_vs_dors.results, dds_stg1_pat_vs_dors.results$padj < 0.01 & abs(dds_stg1_pat_vs_dors.results$log2FoldChange) >= 0.584962501)
dds_stg2_pat_vs_dors.results.degs <- subset(dds_stg2_pat_vs_dors.results, dds_stg2_pat_vs_dors.results$padj < 0.01 & abs(dds_stg2_pat_vs_dors.results$log2FoldChange) >= 0.584962501)
dds_stg1_pat_vs_shd.results.degs <- subset(dds_stg1_pat_vs_shd.results, dds_stg1_pat_vs_shd.results$padj < 0.01 & abs(dds_stg1_pat_vs_shd.results$log2FoldChange) >= 0.584962501)



# Write them to a table
write.table(as.data.frame(dds_stg1_pat_vs_shd.results.degs), file = "dds_stg1_pat_vs_shd.degs_1.5fold_5-22.txt", col.names = TRUE, quote = FALSE, sep = "\t")
write.table(as.data.frame(dds_stg2_pat_vs_dors.results.degs), file = "dds_stg2_pat_vs_dors.degs_1.5fold_5-22.txt", col.names = TRUE, quote = FALSE, sep = "\t")
write.table(as.data.frame(dds_stg1_pat_vs_dors.results.degs), file = "dds_stg1_pat_vs_dors.degs_1.5fold_5-22.txt", col.names = TRUE, quote = FALSE, sep = "\t")






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


### Candidate lists are the set of PM8 genes belonging to relevant Wnt-signaling GO Term
### and also differentially expressed in the given pairwise DGE comparison.

candidates <- c("SULF2|XM_020984167.1|maker-90847-snap-gene-5.13",
"WNT5A|XM_020970602.1|maker-174-augustus-gene-221.10",
"FGF2|XM_021008271.1|maker-689-augustus-gene-11.25",
"BARX1|XM_020980965.1|maker-174-augustus-gene-42.12",
"SULF1|XM_021002089.1|maker-38-augustus-gene-270.5",
"G3BP1|XM_020964239.1|fragGeneModel-551-gene5-frag1",
"CCNY|XM_020990387.1|fragGeneModel-178-gene27-frag1",
"WNT11|XM_020989533.1|maker-90823-augustus-gene-20.15",
"APCDD1|XM_020985551.1|maker-38-augustus-gene-73.9",
"MLLT3|XM_021004579.1|maker-184-snap-gene-6.3",
"CDH2|XM_020992563.1|fragGeneModel-38-gene52-frag1",
"DEPDC1B|XM_021007291.1|maker-229-augustus-gene-76.9",
"ADGRA2|XM_020985273.1|maker-90919-augustus-gene-138.18",
"SHISA2|XM_020997241.1|maker-91133-augustus-gene-196.6",
"ISL1|XM_020993916.1|maker-226-augustus-gene-13.8",
"FGF9|XM_020968777.1|maker-91133-augustus-gene-184.3",
"MDFI|XM_020999646.1|augustus_masked-90352-processed-gene-3.0",
"SFRP1|XM_020969593.1|maker-90919-augustus-gene-126.4",
"GREM1|XM_020995779.1|augustus_masked-37-processed-gene-15.3",
"DKK2|XM_021001334.1|maker-126-augustus-gene-114.11",
"ILK|XM_020975852.1|fragGeneModel-90294-gene2-frag1")

# dds_stg1_pat_vs_shd
volcanoPlotData <- as.data.frame(dds_stg1_pat_vs_shd.results)
volcanoPlotData <- volcanoPlotData[volcanoPlotData$baseMean < 2500, ]
volcanoPlotData <- volcanoPlotData %>% rownames_to_column('geneName')
volcanoPlotData <- volcanoPlotData %>% mutate(candidate = (geneName %in% candidates))
candidateSubset <- volcanoPlotData %>% filter (candidate == TRUE)
candidateSubset <- candidateSubset %>% filter (padj < .01)
upreg <- candidateSubset %>% filter (log2FoldChange >= 0.584962501)
downreg <- candidateSubset %>% filter (log2FoldChange <= -0.584962501)

# Sanity check
#volcanoPlotData %>% filter (candidate == TRUE)

myplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean)) +
    geom_point(color = "grey48", alpha = 0.75, size = 1.5) +
	geom_text_repel(data = downreg, color = "black", ylim = c(50,Inf), xlim = c(-3.0,-4.75), force = 1, size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", downreg$geneName))) +
	geom_text_repel(data = upreg, color = "black", ylim = c(20,Inf), xlim = c(4,7.75), force = 16, size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", upreg$geneName))) +
	geom_point(data = upreg, color = "darkgreen", alpha = 1, size = 2.75) + 
	geom_point(data = downreg, color = "purple", alpha = 1, size = 2.75) +
	xlim(-4.5,7.5) +
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



# dds_stg1_pat_vs_dors
volcanoPlotData <- as.data.frame(dds_stg1_pat_vs_dors.results)
volcanoPlotData <- volcanoPlotData[volcanoPlotData$baseMean < 2500, ]
volcanoPlotData <- volcanoPlotData %>% rownames_to_column('geneName')
volcanoPlotData <- volcanoPlotData %>% mutate(candidate = (geneName %in% candidates))
candidateSubset <- volcanoPlotData %>% filter (candidate == TRUE)
candidateSubset <- candidateSubset %>% filter (padj < .01)
upreg <- candidateSubset %>% filter (log2FoldChange >= 0.584962501)
downreg <- candidateSubset %>% filter (log2FoldChange <= -0.584962501)

# Sanity check
#volcanoPlotData %>% filter (candidate == TRUE)

myplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean)) +
    geom_point(color = "grey48", alpha = 0.75, size = 1.5) +
	geom_text_repel(data = downreg, color = "black", ylim = c(50,Inf), xlim = c(-4.75,-4.75), size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", downreg$geneName))) +
	geom_text_repel(data = upreg, color = "black", ylim = c(20,Inf), xlim = c(3,7.75), force = 256, size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", upreg$geneName))) +
	geom_point(data = upreg, color = "darkgreen", alpha = 1, size = 2.75) + 
	geom_point(data = downreg, color = "purple", alpha = 1, size = 2.75) +
	xlab(paste0("Log 2-fold change")) +
	ylab(paste0("Mean expression")) +
	geom_rangeframe() +
    theme(text=element_text(family="Arial", size=7), legend.position = "none") +
	theme_classic()

ggsave("dds_stg1_pat_vs_dors.png",
	width = 8.76,
	height = 6.62,
	myplot,
	dpi = 600)
	
	
ggsave("dds_stg1_pat_vs_dors.svg",
	width = 8.76,
	height = 6.62,
	myplot,
	dpi = 600)




# dds_stg2_pat_vs_dors
volcanoPlotData <- as.data.frame(dds_stg2_pat_vs_dors.results)
volcanoPlotData <- volcanoPlotData[volcanoPlotData$baseMean < 2500, ]
volcanoPlotData <- volcanoPlotData %>% rownames_to_column('geneName')
volcanoPlotData <- volcanoPlotData %>% mutate(candidate = (geneName %in% candidates))
candidateSubset <- volcanoPlotData %>% filter (candidate == TRUE)
candidateSubset <- candidateSubset %>% filter (padj < .01)
upreg <- candidateSubset %>% filter (log2FoldChange >= 0.584962501)
downreg <- candidateSubset %>% filter (log2FoldChange <= -0.584962501)

# Sanity check
#volcanoPlotData %>% filter (candidate == TRUE)

myplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean)) +
    geom_point(color = "grey48", alpha = 0.75, size = 1.5) +
	geom_text_repel(data = downreg, color = "black", ylim = c(50,Inf), xlim = c(-4.75,-4.75), size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", downreg$geneName))) +
	geom_text_repel(data = upreg, color = "black", ylim = c(20,Inf), xlim = c(4,7.75), force = 16, size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", upreg$geneName))) +
	geom_point(data = upreg, color = "darkgreen", alpha = 1, size = 2.75) + 
	geom_point(data = downreg, color = "purple", alpha = 1, size = 2.75) +
	xlim(-5,8) +
	xlab(paste0("Log 2-fold change")) +
	ylab(paste0("Mean expression")) +
	geom_rangeframe() +
    theme(text=element_text(family="Arial", size=7), legend.position = "none") +
	theme_classic()

ggsave("dds_stg2_pat_vs_dors.png",
	width = 8.76,
	height = 6.62,
	myplot,
	dpi = 600)
	
	
ggsave("dds_stg2_pat_vs_dors.svg",
	width = 8.76,
	height = 6.62,
	myplot,
	dpi = 600)



### Candidate lists are the set of genes differentially expressed in the given pairwise DGE comparison.
### and belonging to enriched limb-morphogenesis GO terms


candidates <- c("SHOX2|XM_020988878.1|maker-173-snap-gene-61.3",
"RARG|XM_020966147.1|maker-454-augustus-gene-5.35",
"PRRX1|XM_020970510.1|maker-220-augustus-gene-80.11",
"WNT5A|XM_020970602.1|maker-174-augustus-gene-221.10",
"OSR2|XM_020969132.1|maker-98-augustus-gene-213.15",
"FREM2|XM_020969996.1|fragGeneModel-90639-gene12-frag1",
"TBX5|XM_021008018.1|maker-99-augustus-gene-61.9",
"HAND2|XM_020997884.1|maker-126-augustus-gene-232.3",
"GREM1|XM_020995779.1|augustus_masked-37-processed-gene-15.3",
"OSR1|XM_020967455.1|maker-559-snap-gene-8.20",
"TBX3|XM_021007970.1|maker-99-snap-gene-60.8",
"COL2A1|XM_020977216.1|fragGeneModel-191-gene13-frag1")


# dds_stg2_pat_vs_dors LIMB GENES HIGHLIGHTED
volcanoPlotData <- as.data.frame(dds_stg2_pat_vs_dors.results)
volcanoPlotData <- volcanoPlotData[volcanoPlotData$baseMean < 2500, ]
volcanoPlotData <- volcanoPlotData %>% rownames_to_column('geneName')
volcanoPlotData <- volcanoPlotData %>% mutate(candidate = (geneName %in% candidates))
candidateSubset <- volcanoPlotData %>% filter (candidate == TRUE)
candidateSubset <- candidateSubset %>% filter (padj < .01)
upreg <- candidateSubset %>% filter (log2FoldChange >= 0.584962501)
downreg <- candidateSubset %>% filter (log2FoldChange <= -0.584962501)

# Sanity check
#volcanoPlotData %>% filter (candidate == TRUE)

myplot <- ggplot(data = volcanoPlotData, aes(x = log2FoldChange, y = baseMean)) +
    geom_point(color = "grey48", alpha = 0.75, size = 1.5) +
	geom_text_repel(data = downreg, color = "black", ylim = c(50,Inf), xlim = c(-4.75,-4.75), size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", downreg$geneName))) +
	geom_text_repel(data = upreg, color = "black", ylim = c(20,Inf), xlim = c(4,7.75), force = 256, size=4, segment.size = 0.5, segment.alpha = 1, segment.color = "black", aes(label=gsub("\\|.*", "", upreg$geneName))) +
	geom_point(data = upreg, color = "darkgreen", alpha = 1, size = 2.75) + 
	geom_point(data = downreg, color = "purple", alpha = 1, size = 2.75) +
	xlim(-5,8) +
	xlab(paste0("Log 2-fold change")) +
	ylab(paste0("Mean expression")) +
	geom_rangeframe() +
    theme(text=element_text(family="Arial", size=7), legend.position = "none") +
	theme_classic()

ggsave("dds_stg2_pat_vs_dors_LIMB.png",
	width = 8.76,
	height = 6.62,
	myplot,
	dpi = 600)
	
	
ggsave("dds_stg2_pat_vs_dors_LIMB.svg",
	width = 8.76,
	height = 6.62,
	myplot,
	dpi = 600)











































