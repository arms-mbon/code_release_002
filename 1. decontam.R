# This code is used to perform Blank Correction using the decontam package
# This is STEP 1 after running PEMA
# For each batch, we need: extended final tables + taxonomy assignments + a file sorting samples into "True Sample" or "Control Sample"

library(phyloseq)
BiocManager::install("decontam")
library(decontam)
library(ggplot2)
install.packages("gtable")
library(tidyr)
library(dplyr)
library(Biostrings)

#----------------COI HCMR data-------------------------

setwd("~/Datapaper/1.blankCorrection/COI")

# Do the following for each batch:
# _predecontam suffixes means that empty samples have been removed and recorded
# Reformatting of PEMA final table
ASVtable_COI <- read.table("extenedFinalTable_Sep20_COI_200924_predecontam.tsv", header = TRUE, sep = "\t")
ASVtable_COI_format <- ASVtable_COI %>%
  separate(col = "ASV_number.amplicon", into = c("ASV_number", "ID"), sep = ":")
OTUtable <- select(ASVtable_COI_format, -c("ASV_number", "Classification", "TAXON.NCBI_TAX_ID"))

# add IDs as row names
rownames(OTUtable) <- OTUtable[[1]]  # This sets the first column as row names
OTUtable <- OTUtable[,-1]   

# Read file containing sample info (True Sample or Control Sample)
Sample <- read.csv("sample_data_Sep20_COI.csv", row.names = 1)

# OUTSIDE OF R: Reformat tax_assignments from PEMA's outputs
# add columns names, remove empty columns + confidence levels, remove tails from ASV number
# Then, read it here
Tax <- read.table("tax_assignments_Sep20_COI_200924.tsv", header = TRUE, sep = "\t") 
Tax <- Tax %>%
  separate(col = "ASV", into = c("ASV", "Number"), sep = "_")
Tax <- select(Tax, -c("Number"))
rownames(Tax) <- Tax[, 1] 
Tax <- select(Tax, -c("ASV"))

# BUild phyloseq object
OTU <- otu_table(OTUtable, taxa_are_rows = TRUE)
sample_data <- sample_data(Sample)
Tax <- as.matrix(Tax) # Ensure Tax is a matrix
taxonomy <- tax_table(Tax)

ps <- phyloseq(OTU, sample_data, taxonomy)

# Remove samples with zero reads
ps <- prune_samples(sample_sums(ps) > 0, ps)


df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

# Optional: Filter by threshold (0.5 example shown)
contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples (prevalence method is based on presence/absence of each sequence, and not on DNA concentration)
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

write.table(df.pa, "contaminantResults_Sep20_COI_200924.csv")

# Check which ASV/OTU are classified as contaminants. 
# Remove them from extended final table outside of R. Save into "extenedFinalTable_..._postdecontam.tsv". Proceed to next step "2. rename_samples". 

#----------------18S HCMR data-------------------------


setwd("~/Datapaper/1.blankCorrection/18S")

# Do the following for each batch:
# _predecontam means that empty samples have been removed and recorded
# Reformatting of PEMA final table
OTUtable <- read.table("extenedFinalTable_Aug23_18S_pr2_0624_predecontam.tsv", header = TRUE, sep = "\t")
rownames(OTUtable) <- OTUtable[[1]]  # This sets the first column as row names
OTUtable <- OTUtable[,-1] 
OTUtable <- select(OTUtable, -c("Classification", "TAXON.NCBI_TAX_ID"))

Sample <- read.csv("sample_data_Aug23_18S.csv", row.names = 1)
Tax <- read.table("tax_Aug23_18S_pr2.txt",sep = "\t", header=T, row.names=1)

# BUild phyloseq object
OTU <- otu_table(OTUtable, taxa_are_rows = TRUE)
sample_data <- sample_data(Sample)
Tax <- as.matrix(Tax) # Ensure Tax is a matrix
taxonomy <- tax_table(Tax)

ps <- phyloseq(OTU, sample_data, taxonomy)

df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)


contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)


# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

write.table(df.pa, "contaminantResults_Aug23_18S.csv")

# Check which ASV/OTU are classified as contaminants. 
# Remove them from extended final table outside of R. Save into "extenedFinalTable_..._postdecontam.tsv". Proceed to next step "2. rename_samples". 


#----------------ITS HCMR data-------------------------

setwd("~/Datapaper/1.blankCorrection/ITS")


# Do the following for each batch:
# _predecontam means that empty samples have been removed and recorded
# Reformatting of PEMA final table
OTUtable <- read.table("extenedFinalTable_Sep20_ITS_210924_predecontam.tsv", header = TRUE, sep = "\t")
rownames(OTUtable) <- OTUtable[[1]]  # This sets the first column as row names
OTUtable <- select(OTUtable, -c("OTU", "Classification", "TAXON.NCBI_TAX_ID"))

# Load sample info (classification into true sample or negative control)
Sample <- read.csv("sample_data_Sep20_ITS.csv", row.names = 1)

# Load tax classification obtained from extended final table, manually reformatted to have the different levels
Tax <- read.table("tax_Sep20_ITS_210924.tsv",sep = "\t", header=T, row.names=1)

# BUild phyloseq object
OTU <- otu_table(OTUtable, taxa_are_rows = TRUE)
sample_data <- sample_data(Sample)
Tax <- as.matrix(Tax) # Ensure Tax is a matrix
taxonomy <- tax_table(Tax)

ps <- phyloseq(OTU, sample_data, taxonomy)

df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)


contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)


# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

write.table(df.pa, "contaminantResults_Sep20_ITS_210924.csv")

# Check which ASV/OTU are classified as contaminants. 
# Remove them from extended final table outside of R. Save into "extenedFinalTable_..._postdecontam.tsv". Proceed to next step "2. rename_samples". 


