# This code is used to merge all COI extended final tables after renaming of samples
# This is STEP 3

library(tidyr)
library(tibble)
library(Biostrings)
library(dplyr)                                                 
library(plyr)                                                  
library(readxl)
library(openxlsx)

setwd("~/Datapaper/3.mergeGenes/COI")


# read Extended_final_tables of each sequencing run
# separate ASV_xy prefixes from long ID code in first column

ASVcounts_April2021 <- read.csv("extenedFinalTable_Apr21_COI_200924_renamed.csv",header=T)
ASVcounts_April2021 <- separate(data = ASVcounts_April2021, col = "OTU", into = c("ASV_number", "ID"), sep = "\\:")

ASVcounts_Sep2020 <- read.csv("extenedFinalTable_Sep20_COI_200924_renamed.csv",header=T)
ASVcounts_Sep2020 <- separate(data = ASVcounts_Sep2020, col = "OTU", into = c("ASV_number", "ID"), sep = "\\:")

ASVcounts_August2023 <- read.csv("extenedFinalTable_Aug23_COI_170924_renamed.csv",header=T)
ASVcounts_August2023 <- separate(data = ASVcounts_August2023, col = "OTU", into = c("ASV_number", "ID"), sep = "\\:")

# read tax_assignment files of each sequencing run 
# /!\ PEMA version (with confidence levels), with contaminating ASVs removed

Taxfile_April2021 <- read.table("tax_assignments_Apr21_COI_200924_postdecontam.tsv",sep="\t")
Taxfile_Sep2020 <- read.table("tax_assignments_Sep20_COI_200924_postdecontam.tsv",sep="\t")
Taxfile_Aug2023 <- read.table("tax_assignments_Aug23_COI_170924_postdecontam.tsv",sep="\t")

# Make lists of Taxfiles and ASVcounts

all_Taxfile<-tibble::lst(Taxfile_April2021,Taxfile_Sep2020, Taxfile_Aug2023) # lst from tibble will preserve object names
all_ASVcounts<-tibble::lst(ASVcounts_April2021,ASVcounts_Sep2020, ASVcounts_August2023) # lst from tibble will preserve object names

# set taxonomy assignments with confidence below 0.8 to NA and subset to ASVs in ASVcounts

for(i in 1:length(all_Taxfile)) { 
  colnames(all_Taxfile[[i]])[1]<-"ID"
  all_Taxfile[[i]] <- separate(data = all_Taxfile[[i]], col = ID, into = c("ID", "Read_number"), sep = "\\_")
  all_Taxfile[[i]] <-all_Taxfile[[i]][,-c(2,4,7,10,13,16,19,22)]
  all_Taxfile[[i]][,2]<-ifelse(all_Taxfile[[i]][,3]<=0.8,NA,all_Taxfile[[i]][,2])
  all_Taxfile[[i]][,4]<-ifelse(all_Taxfile[[i]][,5]<=0.8,NA,all_Taxfile[[i]][,4])
  all_Taxfile[[i]][,6]<-ifelse(all_Taxfile[[i]][,7]<=0.8,NA,all_Taxfile[[i]][,6])
  all_Taxfile[[i]][,8]<-ifelse(all_Taxfile[[i]][,9]<=0.8,NA,all_Taxfile[[i]][,8])
  all_Taxfile[[i]][,10]<-ifelse(all_Taxfile[[i]][,11]<=0.8,NA,all_Taxfile[[i]][,10])
  all_Taxfile[[i]][,12]<-ifelse(all_Taxfile[[i]][,13]<=0.8,NA,all_Taxfile[[i]][,12])
  all_Taxfile[[i]][,14]<-ifelse(all_Taxfile[[i]][,15]<=0.8,NA,all_Taxfile[[i]][,14])
  all_Taxfile[[i]]<-all_Taxfile[[i]][,-c(3,5,7,9,11,13,15)]
  colnames(all_Taxfile[[i]])<-c("ID","Kingdom","Phylum","Class","Order","Family","Genus","Species")
  all_Taxfile[[i]] <- all_Taxfile[[i]] %>% filter(ID %in% all_ASVcounts[[i]]$ID) # Just precautionary measure in case taxonomy files still contain ASVs that were removed during blank-correction
}  

# read fasta files to get sequences 
# all_sequences_grouped_XX_COI.fa files contain all clustered ASV sequences, including ASVs present in blank samples and ASVs not classified taxonomically
# the following fasta files therefore contain more sequences than ASVs present in the Extended_final_table_XX_COI_noBlank.xlsx files

fasta_April2021 <- readDNAStringSet("all_sequences_grouped_Apr21_COI_200924.fa")

fasta_Sep2020 <- readDNAStringSet("all_sequences_grouped_Sep20_COI_200924.fa")

fasta_Aug2023 <- readDNAStringSet("all_sequences_grouped_Aug23_COI_170924.fa")

# Make list of fastas

all_fasta<-tibble::lst(fasta_April2021,fasta_Sep2020,fasta_Aug2023)

# Make fasta files into dataframes and remove read counts from sequence identifiers

IDs<-list()
sequences<-list()
all_seqs<-list()
for(i in 1:length(all_fasta)) {
  IDs[[i]]<-names(all_fasta[[i]])
  sequences[[i]]<-paste(all_fasta[[i]])
  all_seqs[[i]]<-data.frame(IDs[[i]],sequences[[i]])
  colnames(all_seqs[[i]])<-c("ID","seqs")
  all_seqs[[i]]$ID<-gsub("\\_.*","",all_seqs[[i]]$ID)
}  

# Get number of unique ASVs prior to any curation (NOTE THIS NUMBER!)

length(unique(rbind(all_seqs[[1]],all_seqs[[2]],all_seqs[[3]])$seqs))

# Subset to ASVs found in count tables (these are the ASVs with a taxonomic classification remaining after blank-curation)

for(i in 1:length(all_seqs)) {
  all_seqs[[i]] <- all_seqs[[i]] %>% filter(ID %in% all_ASVcounts[[i]]$ID)
}  

# Add sequences to Taxfiles

for(i in 1:length(all_Taxfile)) {
  for(j in 1:nrow(all_Taxfile[[i]])) {
    all_Taxfile[[i]]$seqs[all_Taxfile[[i]]$ID == all_seqs[[i]]$ID[j]] <- all_seqs[[i]][,2][j]
  }
}

# Add sequences to ASVcounts

for(i in 1:length(all_ASVcounts)) {
  for(j in 1:nrow(all_ASVcounts[[i]])) {
    all_ASVcounts[[i]]<-data.frame(all_ASVcounts[[i]])
    all_ASVcounts[[i]]$seqs[all_ASVcounts[[i]]$ID == all_seqs[[i]]$ID[j]] <- all_seqs[[i]][,2][j]
  }
}

# Bring sequence columns to front, remove ID and other columns

all_Taxfile<-lapply(all_Taxfile,function(x) x %>% relocate(seqs))
all_Taxfile<-lapply(all_Taxfile,function(x) x %>% select(-ID))

all_ASVcounts<-lapply(all_ASVcounts,function(x) x %>% relocate(seqs))
all_ASVcounts<-lapply(all_ASVcounts,function(x) x %>% select(-c(ASV_number,ID,Classification,TAXON.NCBI_TAX_ID)))

# Combine Taxfiles and ASVcounts among each other, set NAs to zero

ASVcounts<-rbind.fill(all_ASVcounts[[1]],all_ASVcounts[[2]],all_ASVcounts[[3]])
ASVcounts[is.na(ASVcounts)]<-0

Taxfile<-rbind.fill(all_Taxfile[[1]],all_Taxfile[[2]],all_Taxfile[[3]])

# Aggregate ASV counts based on unique sequences

ASVcounts<-aggregate(.~ seqs,data = ASVcounts,FUN=sum)

# Keep only distinct sequence-taxonomy combinations

Taxfile<-Taxfile %>% distinct()

# Sort Taxfile based on order in ASVcounts

Taxfile<-Taxfile[order(match(Taxfile[,1],ASVcounts[,1])),]

# In the Species column

Taxfile$Species <- ifelse(grepl("*_sp$", Taxfile$Species), NA, Taxfile$Species) # sets all Species entries that end with just _sp as NA

# Make fasta file of unique sequences with short ASV names and save to file

fasta<-data.frame(ASV=paste0(">ASV",1:nrow(ASVcounts)),seqs=ASVcounts$seqs) |>
  pivot_longer(everything()) |> 
  subset(select=-name)

write.table(fasta,"COI_ASVs.fa",row.names = F,col.names = F,quote=F)

# Write ASV table to file

ASVcounts$ASV<-paste0("ASV",1:nrow(ASVcounts))
ASVcounts<-ASVcounts %>% relocate(ASV) %>% select(-seqs)
write.table(ASVcounts,"COI_ASV_table.txt",sep="\t",row.names = F)

# Write taxonomy table to file

Taxfile$ASV<-paste0("ASV",1:nrow(Taxfile))
Taxfile<-Taxfile %>% relocate(ASV) %>% select(-seqs)
write.table(Taxfile,"COI_tax_table.txt",sep="\t",row.names = F)

# get number of unique, classified, blank-corrected ASVs (NOTE THIS NUMBER!)

length(ASVcounts$ASV)

# Get number of reads (NOTE THIS NUMBER!)

sum(colSums(ASVcounts[,-1]))

