# This is the main script for analyses, post-cleaning of the data
# This is STEP 4

library(dplyr)
library(ggplot2)
library(egg)
library(grafify)
library(phyloseq)
library(vegan)
library(grid) # had to be called, even though it should come with ggplot??
library(ggpubr)
library(scales)
library(tidyr)
library(xlsx)
library(readxl)
library(writexl)
library(UpSetR)
install.packages("openxlsx")
library(openxlsx)
library(data.table)


setwd("~/Datapaper")

# Read sample metadata and format a bit

samples<-read.csv("sample_data.csv",header=T,check.names=F,sep=",") 
rownames(samples) <- samples$MaterialSampleID
colnames(samples)[1]<-"MaterialSampleID"
samples$MaterialSampleID<-gsub("_r1","",samples$MaterialSampleID) # Remove the "_r1" strings in MaterialSampleID column
samples$MaterialSampleID<-gsub("_r2","",samples$MaterialSampleID) # Remove the "_r2" strings in MaterialSampleID column

### Create phyloseq objects, do first quick data assessments and get info on recovered phyla species-level assignments ###

#----------------COI-------------------------

# Read ASV counts

ASVcountsCOI<-read.table("COI_ASV_table.txt",header=T,check.names=F, sep="\t",row.names = 1)

# Read ASV taxonomy (needs to be read as matrix, may cause problems otherwise when phyloseq object will be created)

ASVtaxaCOI<-as.matrix(read.table("COI_tax_table.txt",header=T,check.names=F, sep="\t",row.names = 1))

# Sort count table based on order in tax table (precautionary measure)

ASVcountsCOI<-ASVcountsCOI[order(match(rownames(ASVcountsCOI),rownames(ASVtaxaCOI))),]

# Create phyloseq object

psCOI <- phyloseq(otu_table(ASVcountsCOI,taxa_are_rows = TRUE), sample_data(samples), tax_table(ASVtaxaCOI))

# Get number of samples that produced ASVs through PEMA processing

psCOI

# Remove certain erroneous ASVs

to_remove_taxa_COI<-c("Homo","Bos","Canis","Felis", "Sus", "Gynaikothrips","Dorypteryx","Fannia","Bactrocera", "Aleochara", "Larus","Entomobrya", "Hypogastrura","Psylla","Tanytarsus","Ptinus","Paratanytarsus")

psCOI<-subset_taxa(psCOI,!Genus %in% to_remove_taxa_COI) # Make phyloseq object without these sequences

# Remove samples with are left with a read number of zero

psCOI <- prune_samples(sample_sums(psCOI) > 0, psCOI)


# Make data.table for plot
# Violin plot based on sequencing events
read_sums_COI<- data.table(as(sample_data(psCOI), "data.frame"),
                       TotalReads = sample_sums(psCOI), keep.rownames = TRUE)
setnames(read_sums_COI,"rn","SampleID")
reads_plot_COI <- ggplot(read_sums_COI, aes(y=TotalReads,x=Sequenced,color=Sequenced)) + 
  geom_violin(na.rm = TRUE) + 
  geom_jitter(shape=16, position=position_jitter(0.2),size=1,na.rm = TRUE) +
  scale_y_continuous(labels = scales::comma)+
  #scale_x_discrete(limits=c("20-Sep", "21-Apr", "23-Aug"))+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("Sequencing Depth")

reads_plot_COI
ggsave("reads_plot_COI.png", reads_plot_COI, width = 6, height = 4, dpi = 300)



# Check number of ASVs after further curation and filtering

psCOI

# Get total number of reads

sum(sample_sums(psCOI))

# Get number of ASVs assigned to species level

subset_taxa(psCOI,!is.na(Species))

# Get number of unique species identified with Linnean name 

length(unique(tax_table(subset_taxa(psCOI,!is.na(Species)))[,ncol(tax_table(psCOI))]))

# Get number of species observations for occurrences with a minimum of 2 reads (i.e. all presence-absence occurrences of all ASVs classified to species level across all samples)

psCOI_species<-subset_taxa(psCOI,!is.na(Species))

psCOI_species_obs<-psCOI_species

otu_table(psCOI_species_obs)[otu_table(psCOI_species_obs)<2]<-0

otu_table(psCOI_species_obs)[otu_table(psCOI_species_obs)>1]<-1

sum(sample_sums(psCOI_species_obs))

# Get data set agglomerated at phylum level and with relative abundances #

psCOI.phylum <- tax_glom(psCOI, taxrank = "Phylum",NArm=F) 

phylum_COI<-as.data.frame(cbind(as.data.frame(tax_table(psCOI.phylum))[,2],taxa_sums(psCOI.phylum)))

colnames(phylum_COI)<-c("Phylum","reads")

phylum_COI$reads<-as.numeric(phylum_COI$reads)

phylum_COI[is.na(phylum_COI)]<-"NA"

phylum_COI<-aggregate(.~ Phylum,data = phylum_COI,FUN=sum)

phylum_COI<-phylum_COI[order(phylum_COI$reads, decreasing = TRUE),]  

phylum_COI$reads<-phylum_COI$reads/sum(phylum_COI$reads)

# Change classification in the phylum level column to correct phylum name
# This was done based on web-based search using WoRMS or other scientific publications

phylum_COI$Phylum # check all unique entries in the Phylum column of the COI data set

phylum_COI$Phylum<-ifelse(grepl("Ectocarpales|Fucales|Eustigmatales|Laminariales|Cutleriales|Dictyotales|Chattonellales|Chromulinales|Parmales|Sphacelariales|Florenciellales|Goniochloridales|Dictyochales", phylum_COI$Phylum),"Ochrophyta",phylum_COI$Phylum)
phylum_COI$Phylum<-ifelse(grepl("Discosea|Evosea|Tubulinea", phylum_COI$Phylum),"Amoebozoa",phylum_COI$Phylum)
phylum_COI$Phylum<-ifelse(grepl("Peridiniales|Syndiniales|Suessiales|Gymnodiniales|Gonyaulacales", phylum_COI$Phylum),"Myzozoa",phylum_COI$Phylum)
phylum_COI$Phylum<-ifelse(grepl("Myzocytiopsidales|Saprolegniales|Anisolpidiales|Peronosporales|Pythiales|Olpidiopsidales", phylum_COI$Phylum),"Oomycota",phylum_COI$Phylum)
phylum_COI$Phylum<-ifelse(grepl("Imbricatea", phylum_COI$Phylum),"Cercozoa",phylum_COI$Phylum)
phylum_COI$Phylum<-ifelse(grepl("Apusomonadidae", phylum_COI$Phylum),"Apusozoa",phylum_COI$Phylum)
phylum_COI$Phylum<-ifelse(grepl("Ministeria", phylum_COI$Phylum),"Choanozoa",phylum_COI$Phylum)
phylum_COI$Phylum<-ifelse(grepl("Thraustochytrida|Bicosoecida", phylum_COI$Phylum),"Bigyra",phylum_COI$Phylum)
phylum_COI$Phylum<-ifelse(grepl("Mucoromycota", phylum_COI$Phylum),"Zygomycota",phylum_COI$Phylum)

phylum_COI<-aggregate(.~ Phylum,data = phylum_COI,FUN=sum)

phylum_COI<-phylum_COI[order(phylum_COI$reads, decreasing = TRUE),]  

# get number and names of identified phyla
phyla_identified_COI<-phylum_COI
length(which(phyla_identified_COI$Phylum!="NA")) # Number of recovered phyla
phyla_identified_COI[phyla_identified_COI$Phylum=="NA",] # Percentage of reads unclassified at phylum level
write.table(phyla_identified_COI,"phyla_identified_COI.txt",sep="\t",row.names = F,quote = F)

phylum_COI$Phylum[12:nrow(phylum_COI)]<-"Other" # keep only top 10 most abundant phyla and NA and set rest to "Other"

phylum_COI<-aggregate(.~ Phylum,data = phylum_COI,FUN=sum)

# For ASVs with species level classification, get data set agglomerated at phylum  level 

# Extract taxonomic table from phyloseq object
taxa_COI <- as.data.frame(tax_table(psCOI_species))

# Ensure distinct species level classification
specs_COI <- distinct(taxa_COI, Species, .keep_all = TRUE)

# Create abundance column and set to 1 for each species
specs_COI$abundance <- rep(1, nrow(specs_COI))

# Select relevant columns
specs_COI <- specs_COI %>% select(c(Phylum, abundance))

# Aggregate data at Phylum level by summing abundance
specs_COI <- aggregate(. ~ Phylum, data = specs_COI, FUN = sum)

# Order by abundance in decreasing order
specs_COI <- specs_COI[order(specs_COI$abundance, decreasing = TRUE),]

# Change classification in the phylum level column to correct phylum name
# This was done based on web-based search using WoRMS or other scientific publications
# We just use the code from above here, because the species-assigned data set is just a subset of the data set above

specs_COI$Phylum<-ifelse(grepl("Ectocarpales|Fucales|Eustigmatales|Laminariales|Cutleriales|Dictyotales|Chattonellales|Chromulinales|Parmales|Sphacelariales|Florenciellales|Goniochloridales|Dictyochales", specs_COI$Phylum),"Ochrophyta",specs_COI$Phylum)
specs_COI$Phylum<-ifelse(grepl("Discosea|Evosea|Tubulinea", specs_COI$Phylum),"Amoebozoa",specs_COI$Phylum)
specs_COI$Phylum<-ifelse(grepl("Peridiniales|Syndiniales|Suessiales|Gymnodiniales|Gonyaulacales", specs_COI$Phylum),"Myzozoa",specs_COI$Phylum)
specs_COI$Phylum<-ifelse(grepl("Myzocytiopsidales|Saprolegniales|Anisolpidiales|Peronosporales|Pythiales|Olpidiopsidales", specs_COI$Phylum),"Oomycota",specs_COI$Phylum)
specs_COI$Phylum<-ifelse(grepl("Imbricatea", specs_COI$Phylum),"Cercozoa",specs_COI$Phylum)
specs_COI$Phylum<-ifelse(grepl("Apusomonadidae", specs_COI$Phylum),"Apusozoa",specs_COI$Phylum)
specs_COI$Phylum<-ifelse(grepl("Ministeria", specs_COI$Phylum),"Choanozoa",specs_COI$Phylum)
specs_COI$Phylum<-ifelse(grepl("Thraustochytrida|Bicosoecida", specs_COI$Phylum),"Bigyra",specs_COI$Phylum)
specs_COI$Phylum<-ifelse(grepl("Mucoromycota", specs_COI$Phylum),"Zygomycota",specs_COI$Phylum)

specs_COI<-aggregate(.~ Phylum,data = specs_COI,FUN=sum)

specs_COI<-specs_COI[order(specs_COI$abundance, decreasing = TRUE),]  

# get number and names of identified phyla
phyla_specs_identified_COI<-specs_COI
length(phyla_specs_identified_COI$Phylum) # Number of recovered phyla
write.table(phyla_specs_identified_COI,"phyla_specs_identified_COI.txt",sep="\t",row.names = F,quote = F)

specs_COI$Phylum<-ifelse(specs_COI$abundance>=10,specs_COI$Phylum,"Other") # keep only phyla with at least 10 species and set rest to "Other"

specs_COI<-aggregate(.~ Phylum,data = specs_COI,FUN=sum)

#----------------18S-------------------------

# Read OTU counts

OTUcounts18S<-read.table("18S_OTU_table.txt",header=T,check.names=F, sep="\t",row.names = 1)

# Read OTU taxonomy (needs to be read as matrix, may cause problems otherwise when phyloseq object will be created)

OTUtaxa18S<-as.matrix(read.table("18S_tax_table.txt",header=T,check.names=F, sep="\t",row.names = 1))

# Sort count table based on order in tax table (precautionary measure)

OTUcounts18S<-OTUcounts18S[order(match(rownames(OTUcounts18S),rownames(OTUtaxa18S))),]

# Create phyloseq object

ps18S <- phyloseq(otu_table(OTUcounts18S,taxa_are_rows = TRUE), sample_data(samples), tax_table(OTUtaxa18S))

# Get number of samples that produced OTUs through PEMA processing

ps18S

# Remove certain erroneous OTUs

to_remove_taxa_18S<-c("Zea_mays","Stegobium_paniceum", "Umbelopsis_isabellina", "Backusella_ctenidia", "Davidiella_tassiana")

ps18S<-subset_taxa(ps18S,!Species %in% to_remove_taxa_18S) # Make phyloseq object without these sequences
ps18S<-subset_taxa(ps18S,!Level_XXX %in% "Drosophila") # Make phyloseq object without these sequences

# Remove samples with are left with a read number of zero

ps18S <- prune_samples(sample_sums(ps18S) > 0, ps18S)
ps18S

# Make data.table for plot
# Violin plot based on sequencing events
read_sums_18S<- data.table(as(sample_data(ps18S), "data.frame"),
                       TotalReads = sample_sums(ps18S), keep.rownames = TRUE)
setnames(read_sums_18S,"rn","SampleID")
reads_plot_18S <- ggplot(read_sums_18S, aes(y=TotalReads,x=Sequenced,color=Sequenced)) + 
  geom_violin() + 
  geom_jitter(shape=16, position=position_jitter(0.2),size=1) +
  scale_y_continuous(labels = scales::comma)+
  #scale_x_discrete(limits=c("20-Sep", "21-Apr","23-Aug"))+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("Sequencing Depth")

reads_plot_18S
ggsave("reads_plot_18S.png", reads_plot_18S, width = 6, height = 4, dpi = 300)


# Remove OTUs with zeros reads

ps18S<-prune_taxa(rowSums(otu_table(ps18S))>0,ps18S)

# Check number of OTUs

ps18S

# Get total number of reads

sum(sample_sums(ps18S))

# Get number of OTUs assigned to species level

subset_taxa(ps18S,!is.na(Species))

# Get number of unique species identified with Linnean name 

length(unique(tax_table(subset_taxa(ps18S,!is.na(Species)))[,ncol(tax_table(ps18S))]))

# Get number of species observations for occurrences with a minimum of 2 reads (i.e. all presence-absence occurrences of all OTUs classified to species level across all samples)

ps18S_species<-subset_taxa(ps18S,!is.na(Species))

ps18S_species_obs<-ps18S_species

otu_table(ps18S_species_obs)[otu_table(ps18S_species_obs)<2]<-0

otu_table(ps18S_species_obs)[otu_table(ps18S_species_obs)>1]<-1

sum(sample_sums(ps18S_species_obs))

# get data set agglomerated at phylum level and with relative abundances #

ps18S.phylum <- tax_glom(ps18S, taxrank = "Phylum.Class",NArm=F) 

phylum_18S<-as.data.frame(cbind(as.data.frame(tax_table(ps18S.phylum))[,4],taxa_sums(ps18S.phylum)))

colnames(phylum_18S)<-c("Phylum","reads")

phylum_18S$reads<-as.numeric(phylum_18S$reads)

phylum_18S[is.na(phylum_18S)]<-"NA"

phylum_18S<-aggregate(.~ Phylum,data = phylum_18S,FUN=sum)

phylum_18S<-phylum_18S[order(phylum_18S$reads, decreasing = TRUE),]  

phylum_18S$reads<-phylum_18S$reads/sum(phylum_18S$reads)

# Change classification in the phylum/class level column to correct phylum name
# This was done based on web-based search using WoRMS or other scientific publications

phylum_18S$Phylum # check all unique entries in the Phylum/Class column of the 18S data set

phylum_18S$Phylum<-ifelse(grepl("Urochordata|Craniata|Cephalochordata", phylum_18S$Phylum),"Chordata",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Gregarinomorphea|Dinophyceae|Coccidiomorphea|Syndiniales|Colpodellidea|Perkinsida|Ellobiophyceae|Gregarinomorphea", phylum_18S$Phylum),"Myzozoa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Phaeophyceae|Dictyochophyceae|Chrysophyceae|Synurophyceae|Chrysomerophyceae|Raphidophyceae|Pelagophyceae|Eustigmatophyceae|MOCH-5|Pinguiophyceae|Xanthophyceae", phylum_18S$Phylum),"Ochrophyta",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Florideophyceae|Bangiophyceae|Compsopogonophyceae|Rhodellophyceae|Porphyridiophyceae", phylum_18S$Phylum),"Rhodophyta",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Ulvophyceae|Trebouxiophyceae|Pyramimonadophyceae|Embryophyceae|Chlorophyceae|Mamiellophyceae|Chloropicophyceae|Chlorodendrophyceae|Pedinophyceae|Nephroselmidophyceae|Scotinosphaera", phylum_18S$Phylum),"Chlorophyta",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Spirotrichea|Phyllopharyngea|Oligohymenophorea|Heterotrichea|Litostomatea|CONTH|CONThreeP|Cariacotrichea|Colpodea|Protocruziidae|Prostomatea|Karyorelictea|Nassophorea|Plagiopylea", phylum_18S$Phylum),"Ciliophora",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Labyrinthulomycetes|Bicoecea|Placidideae", phylum_18S$Phylum),"Bigyra",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Developea", phylum_18S$Phylum),"Gyrista",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Thecofilosea|Endomyxa-Ascetosporea|Imbricatea|Endomyxa|Phytomyxea|Filosa-Sarcomonadea|Filosa-Granofilosea|Chlorarachniophyceae|Filosa|Novel-clade-10-12|Metromonadea|Filosa-Thecofilosea|Filosa-Imbricatea|Filosa-Metromonadea", phylum_18S$Phylum),"Cercozoa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Ichthyosporea|Choanoflagellatea", phylum_18S$Phylum),"Choanozoa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Monothalamids|Globothalamea|Allogromida|Tubothalamea", phylum_18S$Phylum),"Foraminifera",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Telonemia-Group-2|Telonemia-Group-1|Katablepharidaceae|Cryptophyceae", phylum_18S$Phylum),"Cryptophyta",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Group-1", phylum_18S$Phylum),"Apusozoa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Myxozoa", phylum_18S$Phylum),"Cnidaria",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Loxomitra", phylum_18S$Phylum),"Entoprocta",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Crasiella", phylum_18S$Phylum),"Gastrotricha",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Prymnesiophyceae|HAP5|Pavlovophyceae", phylum_18S$Phylum),"Haptophyta",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Planomonadida|Subulatomonas-lineage|Mantamonadida|YS16Ec34-lineage", phylum_18S$Phylum),"Sulcozoa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Tubulinea|Discosea-Flabellinia|Stygamoebida|LKM74-lineage|Variosea|Centramoebida|Breviata-lineage|Lobosa-G1", phylum_18S$Phylum),"Amoebozoa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Diplonemea|Kinetoplastea|Euglenida|Symbiontida", phylum_18S$Phylum),"Euglenozoa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Prasinodermophyceae", phylum_18S$Phylum),"Prasinodermophyta",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Pterocystida", phylum_18S$Phylum),"Heliozoa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Zygnematophyceae|Zygnemophyceae|Coleochaetophyceae|Klebsormidiophyceae", phylum_18S$Phylum),"Charophyta",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Pirsonia", phylum_18S$Phylum),"Hyphochytridiomycota",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("RAD-B|Polycystinea", phylum_18S$Phylum),"Radiozoa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Baseodiscus", phylum_18S$Phylum),"Nemertea",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Heterolobosea", phylum_18S$Phylum),"Percolozoa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Preaxostyla", phylum_18S$Phylum),"Metamonada",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Gnosonesima", phylum_18S$Phylum),"Platyhelminthes",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("twista", phylum_18S$Phylum),"Alveidia",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Hyphochytridiomycota", phylum_18S$Phylum),"Hyphochytriomyceta",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Endomyxa-Phytomyxea", phylum_18S$Phylum),"Endomyxa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Group-1", phylum_18S$Phylum),"Sulcozoa",phylum_18S$Phylum)

phylum_18S<-aggregate(.~ Phylum,data = phylum_18S,FUN=sum)
phylum_18S<-phylum_18S[order(phylum_18S$reads, decreasing = TRUE),]  

# get number and names of identified phyla
# do not count the MAST clades as unique phyla, they either belong to Bigyra or Gyrista which are already present in the table
phyla_identified_18S<-phylum_18S
length(phyla_identified_18S[!grepl("MAST-",phyla_identified_18S$Phylum),]$Phylum)-1 # Number of recovered phyla; -1 to remove "NA" as phylum
phyla_identified_18S[phyla_identified_18S$Phylum=="NA",] # Percentage of reads unclassified at phylum level
write.table(phyla_identified_18S,"phyla_identified_18S.txt",sep="\t",row.names = F,quote = F)

phylum_18S$Phylum[12:nrow(phylum_18S)]<-"Other" # keep only top 10 most abundant phyla and NA and set rest to "Other"

phylum_18S<-aggregate(.~ Phylum,data = phylum_18S,FUN=sum)

# For OTUs with species level classification, get data set agglomerated at phylum/class level 

taxa_18S<-as.data.frame(tax_table(ps18S_species))

specs_18S<-distinct(taxa_18S,Species,.keep_all=T)

specs_18S$abundance<-rep(1,nrow(specs_18S))

specs_18S<-specs_18S %>% select(c(Phylum.Class,abundance))

specs_18S<-aggregate(.~ Phylum.Class,data = specs_18S,FUN=sum)

# Change classification in the phylum/class level column to correct phylum name
# This was done based on web-based search using WoRMS or other scientific publications
# We just use the code from above here, because the species-assigned data set is just a subset of the data set above

specs_18S$Phylum.Class<-ifelse(grepl("Urochordata|Craniata|Cephalochordata", specs_18S$Phylum.Class),"Chordata",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Gregarnimorphea|Dinophyceae|Coccidiomorphea|Syndiniales|Colpodellidea|Perkinsida|Ellobiophyceae|Gregarinomorphea", specs_18S$Phylum.Class),"Myzozoa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Phaeophyceae|Dictyochophyceae|Chrysophyceae|Synurophyceae|Chrysomerophyceae|Raphidophyceae|Pelagophyceae|Eustigmatophyceae|MOCH-5|Pinguiophyceae|Xanthophyceae", specs_18S$Phylum.Class),"Ochrophyta",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Florideophyceae|Bangiophyceae|Compsopogonophyceae|Rhodellophyceae|Porphyridiophyceae", specs_18S$Phylum.Class),"Rhodophyta",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Ulvophyceae|Trebouxiophyceae|Pyramimonadophyceae|Embryophyceae|Chlorophyceae|Mamiellophyceae|Chloropicophyceae|Chlorodendrophyceae|Pedinophyceae|Nephroselmidophyceae|Scotinosphaera", specs_18S$Phylum.Class),"Chlorophyta",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Spirotrichea|Phyllopharyngea|Oligohymenophorea|Heterotrichea|Litostomatea|CONTH|CONThreeP|Cariacotrichea|Colpodea|Protocruziidae|Prostomatea|Karyorelictea|Nassophorea|Plagiopylea", specs_18S$Phylum.Class),"Ciliophora",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Labyrinthulomycetes|Bicoecea|Placidideae", specs_18S$Phylum.Class),"Bigyra",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Developea", specs_18S$Phylum.Class),"Gyrista",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Thecofilosea|Endomyxa-Ascetosporea|Imbricatea|Endomyxa|Phytomyxea|Filosa-Sarcomonadea|Filosa-Granofilosea|Chlorarachniophyceae|Filosa|Novel-clade-10-12|Metromonadea|Filosa-Thecofilosea|Filosa-Imbricatea|Filosa-Metromonadea", specs_18S$Phylum.Class),"Cercozoa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Ichthyosporea|Choanoflagellatea", specs_18S$Phylum.Class),"Choanozoa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Monothalamids|Globothalamea|Allogromida|Tubothalamea", specs_18S$Phylum.Class),"Foraminifera",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Telonemia-Group-2|Telonemia-Group-1|Katablepharidaceae|Cryptophyceae", specs_18S$Phylum.Class),"Cryptophyta",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Group-1", specs_18S$Phylum.Class),"Apusozoa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Myxozoa", specs_18S$Phylum.Class),"Cnidaria",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Loxomitra", specs_18S$Phylum.Class),"Entoprocta",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Crasiella", specs_18S$Phylum.Class),"Gastrotricha",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Prymnesiophyceae|HAP5|Pavlovophyceae", specs_18S$Phylum.Class),"Haptophyta",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Planomonadida|Subulatomonas-lineage|Mantamonadida|YS16Ec34-lineage", specs_18S$Phylum.Class),"Sulcozoa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Tubulinea|Discosea-Flabellinia|Stygamoebida|LKM74-lineage|Variosea|Centramoebida|Breviata-lineage|Lobosa-G1", specs_18S$Phylum.Class),"Amoebozoa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Diplonemea|Kinetoplastea|Euglenida|Symbiontida", specs_18S$Phylum.Class),"Euglenozoa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Prasinodermophyceae", specs_18S$Phylum.Class),"Prasinodermophyta",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Pterocystida", specs_18S$Phylum.Class),"Heliozoa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Zygnematophyceae|Zygnemophyceae|Klebsormidiophyceae", specs_18S$Phylum.Class),"Charophyta",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Pirsonia", specs_18S$Phylum.Class),"Hyphochytridiomycota",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("RAD-B|Polycystinea", specs_18S$Phylum.Class),"Radiozoa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Baseodiscus", specs_18S$Phylum.Class),"Nemertea",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Heterolobosea", specs_18S$Phylum.Class),"Percolozoa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Preaxostyla", specs_18S$Phylum.Class),"Metamonada",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Gnosonesima", specs_18S$Phylum.Class),"Platyhelminthes",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("twista", specs_18S$Phylum.Class),"Ancoracystidae",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Hyphochytridiomycota", specs_18S$Phylum.Class),"Hyphochytriomyceta",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Endomyxa-Phytomyxea", specs_18S$Phylum.Class),"Endomyxa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Group-1", specs_18S$Phylum.Class),"Sulcozoa",specs_18S$Phylum.Class)

specs_18S<-aggregate(.~ Phylum.Class,data = specs_18S,FUN=sum)

specs_18S<-specs_18S[order(specs_18S$abundance, decreasing = TRUE),]

# get number and names of identified phyla
phyla_specs_identified_18S<-specs_18S
length(phyla_specs_identified_18S$Phylum) # Number of recovered phyla
write.table(phyla_specs_identified_18S,"phyla_specs_identified_18S.txt",sep="\t",row.names = F,quote = F)

specs_18S$Phylum.Class<-ifelse(specs_18S$abundance>2,specs_18S$Phylum.Class,"Other") # Keep only phyla with at least 3 species and set rest to "Other"

specs_18S<-aggregate(.~ Phylum.Class,data = specs_18S,FUN=sum)



#----------------ITS-------------------------

# Read ASV counts

ASVcountsITS<-read.table("ITS_ASV_table.txt",header=T,check.names=F, sep="\t",row.names = 1)

# Read ASV taxonomy (needs to be read as matrix, may cause problems otherwise when phyloseq object will be created)

ASVtaxaITS<-as.matrix(read.table("ITS_tax_table.txt",header=T,check.names=F, sep="\t",row.names = 1))

# Sort count table based on order in tax table (precautionary measure)

ASVcountsITS<-ASVcountsITS[order(match(rownames(ASVcountsITS),rownames(ASVtaxaITS))),]

# Create phyloseq object

psITS <- phyloseq(otu_table(ASVcountsITS,taxa_are_rows = TRUE), sample_data(samples), tax_table(ASVtaxaITS))

# Get number of samples that produced ASVs through PEMA processing

psITS


# Remove ASVs which have a total abundance of zero after removing samples during the previous step
# This also removes ASVs that only were found in the blank

psITS<-prune_taxa(rowSums(otu_table(psITS))>0,psITS)

# Check number of ASVs

psITS

# Violin plot based on sequencing events
read_sums_ITS<- data.table(as(sample_data(psITS), "data.frame"),
                       TotalReads = sample_sums(psITS), keep.rownames = TRUE)
setnames(read_sums_ITS,"rn","SampleID")
reads_plot_ITS <- ggplot(read_sums_ITS, aes(y=TotalReads,x=Sequenced,color=Sequenced)) + 
  geom_violin(na.rm = TRUE) + 
  geom_jitter(shape=16, position=position_jitter(0.2),size=1,na.rm = TRUE) +
  scale_y_continuous(labels = scales::comma)+
  #scale_x_discrete(limits=c("20-Sep", "21-Apr", "23-Aug"))+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("Sequencing Depth")

reads_plot_ITS
ggsave("reads_plot_ITS.png", reads_plot_ITS, width = 6, height = 4, dpi = 300)


# Get total number of reads

sum(sample_sums(psITS))

# Get number of ASVs assigned to species level

subset_taxa(psITS,!is.na(Species))

# Get number of unique species identified with Linnean name 

length(unique(tax_table(subset_taxa(psITS,!is.na(Species)))[,ncol(tax_table(psITS))]))

# Get number of species observations for occurrences with a minimum of 2 reads (i.e. all presence-absence occurrences of all ASVs classified to species level across all samples)

psITS_species<-subset_taxa(psITS,!is.na(Species))

psITS_species_obs<-psITS_species

otu_table(psITS_species_obs)[otu_table(psITS_species_obs)<2]<-0

otu_table(psITS_species_obs)[otu_table(psITS_species_obs)>1]<-1

sum(sample_sums(psITS_species_obs))

# get data set agglomerated at phylum level and with relative abundances

psITS.phylum <- tax_glom(psITS, taxrank = "Phylum",NArm=F) 

phylum_ITS<-as.data.frame(cbind(as.data.frame(tax_table(psITS.phylum))[,4],taxa_sums(psITS.phylum)))

colnames(phylum_ITS)<-c("Phylum","reads")

phylum_ITS$reads<-as.numeric(phylum_ITS$reads)

phylum_ITS[is.na(phylum_ITS)]<-"NA"

phylum_ITS<-aggregate(.~ Phylum,data = phylum_ITS,FUN=sum)

phylum_ITS<-phylum_ITS[order(phylum_ITS$reads, decreasing = TRUE),]  

phylum_ITS$reads<-phylum_ITS$reads/sum(phylum_ITS$reads)

# All entries in the Phylum level column represent correct phylum names. No correction needed.
phylum_ITS$Phylum 

# get number and names of identified phyla
phyla_identified_ITS<-phylum_ITS
length(which(phyla_identified_ITS$Phylum!="NA")) # Get number of phyla 
phyla_identified_ITS[phyla_identified_ITS$Phylum=="NA",] # Percentage of reads unclassified at phylum level
write.table(phyla_identified_ITS,"phyla_identified_ITS.txt",sep="\t",row.names = F,quote = F)

phylum_ITS$Phylum[7:nrow(phylum_ITS)]<-"Other" # Keep only top 5 most abundant phyla and set rest to "Other"

phylum_ITS<-aggregate(.~ Phylum,data = phylum_ITS,FUN=sum)


# For ASVs with species level classification, get data set agglomerated at class  level.
# we chose class level here for a better resolution in the plot.

taxa_ITS<-as.data.frame(tax_table(psITS_species))

specs_ITS<-distinct(taxa_ITS,Species,.keep_all=T)

specs_ITS$abundance<-rep(1,nrow(specs_ITS))

specs_ITS<-specs_ITS %>% select(c(Class,abundance))

specs_ITS<-aggregate(.~ Class,data = specs_ITS,FUN=sum)

specs_ITS<-specs_ITS[order(specs_ITS$abundance, decreasing = TRUE),]

# All entries in the class level column represent correct class names. No correction needed.

# get number and names of identified classes and manually add phylum classification for each class
phyla_specs_identified_ITS<-specs_ITS
length(phyla_specs_identified_ITS$Class) # Number of recovered classes
phyla_specs_identified_ITS$Phylum<-"unknown" # Prepare to add respective phylum classification for the classes
phyla_specs_identified_ITS$Phylum<-ifelse(grepl("Eurotiomycetes|Dothideomycetes|Sordariomycetes|Saccharomycetes|Leotiomycetes|Lecanoromycetes|Pezizomycetes|Taphrinomycetes|Arthoniomycetes", phyla_specs_identified_ITS$Class),"Ascomycota",phyla_specs_identified_ITS$Phylum)
phyla_specs_identified_ITS$Phylum<-ifelse(grepl("Agaricomycetes|Tremellomycetes|Agaricostilbomycetes|Microbotryomycetes|Wallemiomycetes|Malasseziomycetes|Cystobasidiomycetes", phyla_specs_identified_ITS$Class),"Basidiomycota",phyla_specs_identified_ITS$Phylum)
phyla_specs_identified_ITS$Phylum<-ifelse(grepl("Mortierellomycetes", phyla_specs_identified_ITS$Class),"Mortierellomycota",phyla_specs_identified_ITS$Phylum)
phyla_specs_identified_ITS$Phylum<-ifelse(grepl("Mucoromycetes", phyla_specs_identified_ITS$Class),"Mucoromycota",phyla_specs_identified_ITS$Phylum)
write.table(phyla_specs_identified_ITS,"phyla_class_specs_identified_ITS.txt",sep="\t",row.names = F,quote = F)

specs_ITS$Class<-ifelse(specs_ITS$abundance>2,specs_ITS$Class,"Other") # Keep only phyla with at least 3 species and set rest to "Other"

specs_ITS<-aggregate(.~ Class,data = specs_ITS,FUN=sum)

#--------merge reads plots----

library(cowplot)

combined_reads_plot <- plot_grid(reads_plot_COI, reads_plot_18S, reads_plot_ITS, ncol = 1, labels = c("A", "B", "C"))

# Display combined plot
combined_reads_plot

# Save the combined figure
ggsave("combined_reads_plot.png", combined_reads_plot, width = 12, height = 12, dpi = 300)


#--------Plot overall info on phyla and species recovered ###------------------------

phyla<-unique(c(phyla_identified_COI$Phylum,phyla_identified_18S$Phylum,phyla_identified_ITS$Phylum)) # vector with all unique phyla of all genes
length(which(!grepl("MAST-",phyla)))-1 # Number of recovered phyla excluding the MASTs (they are already included in either Bigyra or Gyrista phyla); -1 to remove "NA" as phylum

species<-unique(c(phyla_specs_identified_COI$Phylum,phyla_specs_identified_18S$Phylum,phyla_specs_identified_ITS$Phylum)) # vector with all unique species of all genes
length(species)

# Plots of phyla recovered and phylum/class classification of species recovered

# Create table with colors for each Phylum/Class 

phylum_palette<-as.data.frame(c(phylum_COI[,1],phylum_18S[,1],phylum_ITS[,1],specs_ITS[,1],specs_18S[,1],specs_COI[,1]))

phylum_palette <- as.data.frame(
  c(as.character(phylum_COI[,1]),
    as.character(phylum_18S[,1]),
    as.character(phylum_ITS[,1]),
    as.character(specs_ITS[,1]),
    as.character(specs_18S[,1]),
    as.character(specs_COI[,1]))
)


phylum_palette<-unique(phylum_palette)

write.table(phylum_palette,"palette_to_be_formatted.txt",sep = "\t",col.names = F,row.names = F,quote=F)

# Check kelly palette of grafify package for colorblind-friendly codes

graf_palettes$kelly

# outside of R, add color codes as a second column to palette_to_be_formatted.txt in Excel. Give the fungal classes colors according to the phylum they belong to
# the kelly palette contains less colors than needed, choose a few more suitable hexacodes and add them manually
# Save file as palette.txt once ready
# try out plotting with the colors and switch hexacodes around as needed for a nice order of colors in the plots

# read prepared palette table back into R

phylum_palette<-read.table("palette.txt",sep="\t",stringsAsFactors = FALSE,na.strings = "",comment.char = "")

colnames(phylum_palette)<-c("Phylum","color")

# Plot with gene names on the left side

# Reorder levels to place "NA" at the end
sorted_phylum_COI <- phylum_COI %>% 
  arrange(desc(reads)) %>% 
  filter(!(Phylum %in% c("Other", "NA")))
sorted_phylum_COI <- rbind(sorted_phylum_COI, phylum_COI %>% filter(Phylum == "Other"), phylum_COI %>% filter(Phylum == "NA"))

sorted_phylum_18S <- phylum_18S %>% 
  arrange(desc(reads)) %>% 
  filter(!(Phylum %in% c("Other", "NA")))
sorted_phylum_18S <- rbind(sorted_phylum_18S, phylum_18S %>% filter(Phylum == "Other"), phylum_18S %>% filter(Phylum == "NA"))

sorted_phylum_ITS <- phylum_ITS %>% 
  arrange(desc(reads)) %>% 
  filter(!(Phylum %in% c("Other", "NA")))
sorted_phylum_ITS <- rbind(sorted_phylum_ITS, phylum_ITS %>% filter(Phylum == "Other"), phylum_ITS %>% filter(Phylum == "NA"))

sorted_specs_COI <- specs_COI %>% 
  arrange(desc(abundance)) %>% 
  filter(Phylum != "Other")
sorted_specs_COI <- rbind(sorted_specs_COI, specs_COI %>% filter(Phylum == "Other"))

sorted_specs_18S <- specs_18S %>% 
  arrange(desc(abundance)) %>% 
  filter(Phylum.Class != "Other")
sorted_specs_18S <- rbind(sorted_specs_18S, specs_18S %>% filter(Phylum.Class == "Other"))

sorted_specs_ITS <- specs_ITS %>% 
  arrange(desc(abundance)) %>% 
  filter(Class != "Other")
sorted_specs_ITS <- rbind(sorted_specs_ITS, specs_ITS %>% filter(Class == "Other"))


# Update factor levels for Phylum
phylum_COI$Phylum <- factor(phylum_COI$Phylum, levels = sorted_phylum_COI$Phylum)
phylum_18S$Phylum <- factor(phylum_18S$Phylum, levels = sorted_phylum_18S$Phylum)
phylum_ITS$Phylum <- factor(phylum_ITS$Phylum, levels = sorted_phylum_ITS$Phylum)
specs_COI$Phylum <- factor(specs_COI$Phylum, levels = sorted_specs_COI$Phylum)
specs_18S$Phylum.Class <- factor(specs_18S$Phylum.Class, levels = sorted_specs_18S$Phylum.Class)
specs_ITS$Class <- factor(specs_ITS$Class, levels = sorted_specs_ITS$Class)

cols1<-phylum_palette[phylum_palette$Phylum %in% phylum_COI$Phylum,]
cols1<-cols1[order(match(cols1[,1],phylum_COI$Phylum)),]
cols1<-cols1$color
phylum_COI_plot<-ggplot(phylum_COI, aes(x=Phylum, y=reads, fill=Phylum))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=levels(phylum_COI$Phylum))+
  scale_y_continuous(limits=c(0,0.65),breaks = c(0,.1,.2,.3,.4,.5,.6))+
  scale_fill_manual(breaks=phylum_COI$Phylum,values=cols1)+
  ggtitle("A")+
  annotate("text", y = 0.25, x = -4, label = "COI", size = 7, fontface =2)+
  coord_cartesian(xlim = c(1,12),clip="off")+
  theme_bw()+
  theme(plot.margin=unit(c(.2,.5,0,2.2), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",axis.text.y=element_text(color="black",size=10),axis.text.x=element_text(angle=45,hjust=1,color="black",size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank())

cols2<-phylum_palette[phylum_palette$Phylum %in% phylum_18S$Phylum,]
cols2<-cols2[order(match(cols2[,1],phylum_18S$Phylum)),]
cols2<-cols2$color
phylum_18S_plot<-ggplot(phylum_18S, aes(x=Phylum, y=reads, fill=Phylum))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=levels(phylum_18S$Phylum))+
  scale_y_continuous(limits=c(0,0.25),breaks = c(0,.1,.2))+
  scale_fill_manual(breaks=phylum_18S$Phylum,values=cols2)+
  ggtitle("C")+
  annotate("text", y = 0.15, x = -4, label = "18S", size = 7, fontface =2)+
  coord_cartesian(xlim = c(1,12),clip="off")+
  labs(y="Relative abundance")+
  theme_bw()+
  theme(plot.margin=unit(c(-0.3,.5,-0.3,2.2), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.x=element_blank(),legend.position="none",axis.text.y=element_text(color="black",size=10),axis.text.x=element_text(angle=45,hjust=1,color="black",size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank(),axis.title.y=element_text(size=18))

cols3<-phylum_palette[phylum_palette$Phylum %in% phylum_ITS$Phylum,]
cols3<-cols3[order(match(cols3[,1],phylum_ITS$Phylum)),]
cols3<-cols3$color
phylum_ITS_plot<-ggplot(phylum_ITS, aes(x=Phylum, y=reads, fill=Phylum))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=levels(phylum_ITS$Phylum))+
  scale_y_continuous(limits=c(0,0.6),breaks = c(0,0.2,0.4,0.6))+
  scale_fill_manual(breaks=phylum_ITS$Phylum,values=cols3)+
  ggtitle("E")+
  annotate("text", y = 0.3, x = -2.2, label = "ITS", size = 7, fontface =2)+
  coord_cartesian(xlim = c(1,7),clip="off")+
  theme_bw()+
  theme(plot.margin=unit(c(0,.5,.2,2.2), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position="none",axis.text.y=element_text(color="black",size=10),axis.text.x=element_text(angle=45,hjust=1,color="black",size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank())

cols4<-phylum_palette[phylum_palette$Phylum %in% specs_COI$Phylum,]
cols4<-cols4[order(match(cols4[,1],specs_COI$Phylum)),]
cols4<-cols4$color
specs_COI_plot<-ggplot(specs_COI, aes(x=Phylum, y=abundance, fill=Phylum))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=levels(specs_COI$Phylum))+
  scale_y_continuous(limits=c(0,130),breaks = c(0,30,60,90,120))+
  scale_fill_manual(breaks=specs_COI$Phylum,values=cols4)+
  ggtitle("B")+
  theme_bw()+
  theme(plot.margin=unit(c(.2,.2,0,0), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",axis.text.y=element_text(color="black",size=10),axis.text.x=element_text(angle=45,hjust=1,color="black",size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank())

cols5<-phylum_palette[phylum_palette$Phylum %in% specs_18S$Phylum,]
cols5<-cols5[order(match(cols5[,1],specs_18S$Phylum)),]
cols5<-cols5$color
specs_18S_plot<-ggplot(specs_18S, aes(x=Phylum.Class, y=abundance, fill=Phylum.Class))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=levels(specs_18S$Phylum.Class))+
  scale_y_continuous(limits=c(0,40),breaks = c(0,10,20,30,40))+
  labs(y="Number of species")+
  ggtitle("D")+
  scale_fill_manual(breaks=specs_18S$Phylum,values=cols5)+
  theme_bw()+
  theme(plot.margin=unit(c(-0.3,.2,-0.3,0), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.x=element_blank(),legend.position="none",axis.text.y=element_text(color="black",size=10),axis.text.x=element_text(angle=45,hjust=1,color="black",size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank(),axis.title.y=element_text(size=18))

cols6<-phylum_palette[phylum_palette$Phylum %in% specs_ITS$Class,]
cols6<-cols6[order(match(cols6[,1],specs_ITS$Class)),]
cols6<-cols6$color
specs_ITS_plot<-ggplot(specs_ITS, aes(x=Class, y=abundance, fill=Class))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=levels(specs_ITS$Class))+
  scale_fill_manual(breaks=specs_ITS$Class,values=cols6)+
  ggtitle("F")+
  theme_bw()+
  theme(plot.margin=unit(c(0,.2,.2,0), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",axis.text.y=element_text(color="black",size=10),axis.text.x=element_text(angle=45,hjust=1,color="black",size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank())

phylum_plot_genes<-egg::ggarrange(phylum_COI_plot,specs_COI_plot,phylum_18S_plot,specs_18S_plot,phylum_ITS_plot,specs_ITS_plot,nrow=3)
phylum_plot_genes

ggsave(phylum_plot_genes,file="Figures/Figure_2.pdf",height=8,width=8,device = cairo_pdf,dpi=600)

# Plot without gene names

cols1<-phylum_palette[phylum_palette$Phylum %in% phylum_COI$Phylum,]
cols1<-cols1[order(match(cols1[,1],phylum_COI$Phylum)),]
cols1<-cols1$color
phylum_COI_plot<-ggplot(phylum_COI, aes(x=Phylum, y=reads, fill=Phylum))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=phylum_COI$Phylum)+
  scale_y_continuous(limits=c(0,0.5),breaks = c(0,.1,.2,.3,.4,.5))+
  scale_fill_manual(breaks=phylum_COI$Phylum,values=cols1)+
  ggtitle("A")+
  theme_bw()+
  theme(plot.margin=unit(c(0,0,-0.3,0), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",axis.text.y=element_text(color="black"),axis.text.x=element_text(angle=45,hjust=1,color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank())

cols2<-phylum_palette[phylum_palette$Phylum %in% phylum_18S$Phylum,]
cols2<-cols2[order(match(cols2[,1],phylum_18S$Phylum)),]
cols2<-cols2$color
phylum_18S_plot<-ggplot(phylum_18S, aes(x=Phylum, y=reads, fill=Phylum))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=phylum_18S$Phylum)+
  scale_y_continuous(limits=c(0,0.3),breaks = c(0,.1,.2,.3))+
  scale_fill_manual(breaks=phylum_18S$Phylum,values=cols2)+
  ggtitle("C")+
  labs(y="Relative abundance")+
  theme_bw()+
  theme(plot.margin=unit(c(0,0,-0.3,0), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.x=element_blank(),legend.position="none",axis.text.y=element_text(color="black"),axis.text.x=element_text(angle=45,hjust=1,color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank(),axis.title.y=element_text(size=15))

cols3<-phylum_palette[phylum_palette$Phylum %in% phylum_ITS$Phylum,]
cols3<-cols3[order(match(cols3[,1],phylum_ITS$Phylum)),]
cols3<-cols3$color
phylum_ITS_plot<-ggplot(phylum_ITS, aes(x=Phylum, y=reads, fill=Phylum))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=phylum_ITS$Phylum)+
  scale_y_continuous(limits=c(0,0.45),breaks = c(0,0.2,0.4))+
  scale_fill_manual(breaks=phylum_ITS$Phylum,values=cols3)+
  ggtitle("E")+
  theme_bw()+
  theme(plot.margin=unit(c(0,0,0,0), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position="none",axis.text.y=element_text(color="black"),axis.text.x=element_text(angle=45,hjust=1,color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank())

cols4<-phylum_palette[phylum_palette$Phylum %in% specs_COI$Phylum,]
cols4<-cols4[order(match(cols4[,1],specs_COI$Phylum)),]
cols4<-cols4$color
specs_COI_plot<-ggplot(specs_COI, aes(x=Phylum, y=abundance, fill=Phylum))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=specs_COI$Phylum)+
  scale_y_continuous(limits=c(0,130),breaks = c(0,30,60,90,120))+
  scale_fill_manual(breaks=specs_COI$Phylum,values=cols4)+
  ggtitle("B")+
  theme_bw()+
  theme(plot.margin=unit(c(0,0,-0.3,.5), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",axis.text.y=element_text(color="black"),axis.text.x=element_text(angle=45,hjust=1,color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank())

cols5<-phylum_palette[phylum_palette$Phylum %in% specs_18S$Phylum,]
cols5<-cols5[order(match(cols5[,1],specs_18S$Phylum)),]
cols5<-cols5$color
specs_18S_plot<-ggplot(specs_18S, aes(x=Phylum.Class, y=abundance, fill=Phylum.Class))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=specs_18S$Phylum.Class)+
  scale_y_continuous(limits=c(0,40),breaks = c(0,10,20,30,40))+
  labs(y="Number of species")+
  ggtitle("D")+
  scale_fill_manual(breaks=specs_18S$Phylum,values=cols5)+
  theme_bw()+
  theme(plot.margin=unit(c(0,0,-0.3,.5), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.x=element_blank(),legend.position="none",axis.text.y=element_text(color="black"),axis.text.x=element_text(angle=45,hjust=1,color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank(),axis.title.y=element_text(size=15))

cols6<-phylum_palette[phylum_palette$Phylum %in% specs_ITS$Class,]
cols6<-cols6[order(match(cols6[,1],specs_ITS$Class)),]
cols6<-cols6$color
specs_ITS_plot<-ggplot(specs_ITS, aes(x=Class, y=abundance, fill=Class))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=specs_ITS$Class)+
  scale_fill_manual(breaks=specs_ITS$Class,values=cols6)+
  ggtitle("F")+
  theme_bw()+
  theme(plot.margin=unit(c(0,0,0,.5), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",axis.text.y=element_text(color="black"),axis.text.x=element_text(angle=45,hjust=1,color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank())

phylum_plot<-egg::ggarrange(phylum_COI_plot,specs_COI_plot,phylum_18S_plot,specs_18S_plot,phylum_ITS_plot,specs_ITS_plot,nrow=3)

ggsave(phylum_plot,file="Figures/phylum_plot.png",height=8,width=8,type="cairo")

# ---- UpSet plot for number of species identified which are shared between the marker gene data sets ###----

# Get unique species occurrences for each gene

list_coi<-as.data.frame(unique(as.data.frame(tax_table(psCOI_species))$Species))
list_coi$COI<-1
colnames(list_coi)[1]<-"Species"

list_18s<-as.data.frame(unique(as.data.frame(tax_table(ps18S_species))$Species))
list_18s$'18S'<-1
colnames(list_18s)[1]<-"Species"

list_its<-as.data.frame(unique(as.data.frame(tax_table(psITS_species))$Species))
list_its$ITS<-1
colnames(list_its)[1]<-"Species"

# merge tables by species name, set NAs to zero and set species names as rownames

list_all<-merge(list_coi,list_18s,by="Species",all=T)
list_all<-merge(list_all,list_its,by="Species",all=T)

list_all[is.na(list_all)] <- 0
rownames(list_all) <- list_all$Species
list_all$Species <- NULL

# Create UpSet plot and save it to file

png(file="upset.png",height = 1500, width = 1200,res=300,type="cairo")
upset(list_all, main.bar.color = "#009E73", sets.bar.color = "lightgray", matrix.color = "black", order.by = "freq", decreasing = T, set_size.show = T,set_size.scale_max=890,
      mainbar.y.label = "No. of species", sets.x.label = "No. of species",text.scale = c(1.3, 1.3, 1, 1, 1.3, 1.5))
dev.off()

## at the Genus level
# Get unique genus occurrences for each gene

list_coi_genus <- as.data.frame(unique(as.data.frame(tax_table(psCOI_species))$Genus))
list_coi_genus$COI <- 1
colnames(list_coi_genus)[1] <- "Genus"

list_18s_genus <- as.data.frame(tax_table(ps18S_species))
list_18s_genus$Genus <- sapply(strsplit(as.character(list_18s_genus$Species), "_"), `[`, 1)
list_18s_genus <- as.data.frame(unique(list_18s_genus$Genus))
list_18s_genus$'18S' <- 1
colnames(list_18s_genus)[1] <- "Genus"

list_its_genus <- as.data.frame(unique(as.data.frame(tax_table(psITS_species))$Genus))
list_its_genus$ITS <- 1
colnames(list_its_genus)[1] <- "Genus"

# Merge tables by genus name, set NAs to zero, and set genus names as rownames

list_all_genus <- merge(list_coi_genus, list_18s_genus, by = "Genus", all = TRUE)
list_all_genus <- merge(list_all_genus, list_its_genus, by = "Genus", all = TRUE)

list_all_genus[is.na(list_all_genus)] <- 0
rownames(list_all_genus) <- list_all_genus$Genus
list_all_genus$Genus <- NULL

# Create UpSet plot and save it to file

png(file = "Figures/upset_genus.png", height = 1500, width = 1200, res = 300, type = "cairo")
upset(list_all_genus, 
      main.bar.color = "#009E73", 
      sets.bar.color = "lightgray", 
      matrix.color = "black", 
      order.by = "freq", 
      decreasing = TRUE, 
      set_size.show = TRUE, 
      set_size.scale_max = 890,
      mainbar.y.label = "No. of genera", 
      sets.x.label = "No. of genera", 
      text.scale = c(1.3, 1.3, 1, 1, 1.3, 1.5))
dev.off()




### ---------------- comparison DP001/DP002 ---------------
# Create a data frame with the data
data <- data.frame(
  Marker = rep(c("COI", "18S", "ITS"), each = 2),
  DataPaper = rep(c("DP1", "DP2"), times = 3),
  ReadsPerSample = c(8292, 9392, 20688, 8098, 1185, 682)
)

# Plot the data
comparison <- ggplot(data, aes(x = Marker, y = ReadsPerSample, fill = DataPaper)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Marker", y = "Average Reads per Sample") +
  theme_minimal() +
  scale_fill_manual(values = c("#20A387FF", "#95D840FF"),
                    name = "Data Release",
                    labels = c("DP001", "DP002"))


ggsave(comparison,file="Figures/Comparison_reads_DP001_DP002.png",width=5,height=4,bg="transparent",dpi=600)


# Libraries
library(ggplot2)

# Data preparation
data <- data.frame(
  Marker = c("COI", "18S", "ITS"),
  DP001 = c(273.978836, 50.77248677, 12.0952381),
  DP002 = c(112.254902, 42.28378378, 12.9)
)

# Reshaping data for visualization
data_long <- reshape2::melt(data, id.vars = "Marker", variable.name = "DataRelease", value.name = "ASVs_OTUs")

# Set the order of markers
data_long$Marker <- factor(data_long$Marker, levels = c("COI", "18S", "ITS"))

# Custom colors
custom_colors <- c("DP1" = "#20A387FF", "DP2" = "#95D840FF")

# Bar plot to compare ASVs/OTUs per sample
comparison2 <- ggplot(data_long, aes(x = Marker, y = ASVs_OTUs, fill = DataRelease)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  labs(
    x = "Marker Gene",
    y = "Average Number of ASVs/OTUs per Sample",
    fill = "Data Release"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 12)
  )

ggsave(comparison2,file="Figures/Comparison_ASVsOTUs_DP001_DP002.png",width=5,height=4,bg="transparent",dpi=600)


