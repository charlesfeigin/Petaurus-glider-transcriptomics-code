library(WGCNA)
library(gplots)

options(stringsAsFactors = FALSE);

################# Prepare input data

# Load in sample info
sample_info <- read.csv("sample_info_pat.txt", header = TRUE)
row.names(sample_info) <- sample_info[,1]
sample_info <- sample_info[,3:4]
sample_info <- sample_info[c(grep("pat", rownames(sample_info))),]

# Load in normalized read count data
# Make sure rows/columns are arranged the right way
# Pre-remove non-count columns
allData <- read.table("rld_matrix_pat.txt", header = TRUE)
colnames(allData) <- gsub(x = colnames(allData), pattern = 'X', replacement = "")

# transpose rows/columns
allData <- t(allData)

# Load WGCNA
library(WGCNA)

# Said as mandatory in WGCNA tutorial without explanation...ok!
options(stringsAsFactors = FALSE);


# If above returns FALSE (it will), remove the bad genes
if (!gsg$allOK)
{
 # Optionally, print the gene and sample names that were removed:
 if (sum(!gsg$goodGenes)>0) 
   #printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
 if (sum(!gsg$goodSamples)>0) 
   #printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
 # Remove the offending genes and samples from the data:
 allData <- allData[gsg$goodSamples, gsg$goodGenes]
}

collectGarbage()

#enableWGCNAThreads()
allowWGCNAThreads(nThreads=6)

powers <- c(1:30)

# Find soft threshold
sft = pickSoftThreshold(allData, powerVector = powers, 
			networkType = "signed", corFnc = "cor",
			corOptions = list(use = 'p'), blockSize = 20000,
			verbose = 5)

# plot to visualize
cex1 = 0.9;
pdf("pat_scalefree_cor_bigblock_pat.pdf")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
	main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	labels=powers,cex=cex1,col="red");
abline(h=0.8,col="red")
dev.off()

#pdf("pat_scalefree_bicor_bigblock_pat.pdf")


net1 = blockwiseModules(allData, power = 11,
	TOMType = "unsigned", networkType = "signed",
	deepslit = 2, minKMEtoStay = 0, maxBlockSize = 20000,
	minModuleSize = 50, mergeCutHeight = 0.30,
	#reassignThreshold = 0, corType = "bicor",
	#corOptions = list(use = 'p', maxPOutliers = 0.1),
	reassignThreshold = 0, corType = "pearson",
	corOptions = list(use = 'p'),
	numericLabels = TRUE, pamRespectsDendro = FALSE,
#	saveTOMs = TRUE, saveTOMFileBase = "pat-NetworkConstruction-auto_pearson_sft14",
	verbose = 3)


table(net1$colors)


## Write network genes to files
# Load stringr for string replacement
library(stringr)

# Make a range representing the number of modules
x <- c(0:13)
# Write each module gene list to a file
for (i in x) {
	write(names(as.data.frame(allData))[net1$colors == i], str_replace("pat_network_module_Z.txt", "Z", toString(i)))
}


# Parse out relevant info for subsequent analysis
moduleLabels = net1$colors
moduleColors = labels2colors(net1$colors)
MEs = net1$MEs

# Calculate module-trait correlations
moduleTraitCor <- cor(MEs, sample_info, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Reorganize module eigengene data
MEs0b <- moduleEigengenes(allData, moduleLabels)$eigengenes
row.names(MEs0b) <- row.names(allData)
MEsb <- orderMEs(MEs0b)
MEsc <- t(MEsb)

# Set order of modules to be displayed
x <- c("ME8", "ME4", "ME6", "ME2", "ME13", "ME9", "ME10", "ME5", "ME12", "ME1", "ME11", "ME7", "ME3")
MEsd <- MEsc[x,,drop=TRUE]

# Plot module eigengene expression against age-sorted samples
pdf(file="pat_network_ME_expression_vs_samples.pdf")
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = MEsd,
	xLabels = colnames(MEsd),
	yLabels = row.names(MEsd),
	ySymbols = row.names(MEsd),
	colorLabels = FALSE,
	colors = colorRampPalette(c("purple", "black", "yellow"))(10000),
	setStdMargins = FALSE,
	cex.lab=.5,
	cex.text = 0.75,
	main = paste("Module-Sample relationships"))
dev.off()
