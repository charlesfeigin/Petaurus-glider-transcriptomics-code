library(magrittr)
#library(apeglm)

##################################################################################
#########################################
# Prep input filtered read counts
#########################################

# Read in feature counts table (only pat reads)
read.counts <- read.table("pbrev_read_counts.txt", header = TRUE)

# Set the row names 
row.names(read.counts) <- read.counts$Geneid

# Remove irrelevant data
read.counts <- read.counts[, -c(1:6)]

# Get rid of "X" character added to sample names
betterNames <- gsub( "X", "", names(read.counts))

# Update names
names(read.counts) <- betterNames

# Get tissue type from name
tissueType <- gsub("^[0-9]+[A-Z]_", "", betterNames)
tissueType <- gsub("_stg[1-3]+_[0-9]+$", "", tissueType)

# Get dev stage from name
devStage <- gsub("^[0-9]+[A-Z]_[a-z]+_", "", betterNames)
devStage <- gsub("_[0-9]+$", "", devStage)

# Load in animals weights separately
sampleWeight <- read.table("pbrev_pat_sample_weights.txt", header=TRUE)


# Recode tissue as a number (pat = 1, shd = 2, dls = 3, latter two will just be tossed)
count <- 1
for (val in tissueType) {
	if (val == "pat") {
	tissueType[count] <- 3
	} else if (val == "shd") {
	tissueType[count] <- 1
	} else if (val == "dls") {
	tissueType[count] <- 2
	}
	count <- count + 1
}
rm(count)

# Recode stage as a number: stg1 = 1, stg2 = 2, stg3 = 3 (stg3 = samples after patagium has outgrown)
count <- 1
for (val in devStage) {
	if (val == "stg1") {
	devStage[count] <- 1
	} else if (val == "stg2") {
	devStage[count] <- 2
	} else if (val == "stg3") {
	devStage[count] <- 3
	}
	count <- count + 1
}

# Make a data frame of sample info
sample_info <- data.frame(sample = betterNames, weight = sampleWeight$weight)

row.names(sample_info) = betterNames

# Load DESeq2
library(DESeq2)

# Convert read counts dataframe into a matrix
sample_matrix <- data.matrix(read.counts)

# Remove genes that are effectively off
sample_matrix <- sample_matrix[rowMeans(sample_matrix) >= 10, ]

# Remove non-patagium samples for this analysis
pat_sample_matrix <- sample_matrix[,c(grep("pat", colnames(sample_matrix)))]
pat_sample_info <- sample_info[c(grep("pat", rownames(sample_info))),]

# Remove genes with mean counts < 10. Do above, need all to have same # genes
#pat_sample_matrix <- pat_sample_matrix[rowMeans(pat_sample_matrix) >= 10, ]

# Perform rld on counts matrix
pat_rld_matrix <- rlog(pat_sample_matrix, blind=TRUE)
rownames(pat_rld_matrix) <- rownames(pat_sample_matrix)

# Write filtered rld matrix and pat sample info to files
write.table(pat_rld_matrix, file="rld_matrix_pat.txt", sep = "\t")
write.table(pat_sample_info, file="sample_info_pat.txt", sep = ",", row.names = FALSE, quote = FALSE)

