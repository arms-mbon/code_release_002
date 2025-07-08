# Our samples are named according to their Accession numbers ("ERRXXXX")
# This code is used to rename them with their correct Material Sample IDs
# This is STEP 2: this step comes after blank correction and before merging the different final tables

# Load necessary libraries
library(dplyr)


#----------------COI-------------------------

#Set working directory
setwd("~/Datapaper/2.sampleRenaming/COI")

# Read the extended tables into a list
file_paths_COI <- c(
  "extenedFinalTable_Apr21_COI_200924_postdecontam.csv",
  "extenedFinalTable_Sep20_COI_200924_postdecontam.csv",
  "extenedFinalTable_Aug23_COI_170924_postdecontam.csv"
)

# Load all tables into a list
extended_tables_COI <- lapply(file_paths_COI, function(path) {
  read.csv(path, row.names = 1)
})

# Read the mapping table (ERR number to Material Sample ID)
mapping_table_COI <- read.csv("mappingTable_COI.csv")
colnames(mapping_table_COI) <- c("NewName", "OldName")

# Create a named vector for replacement
name_replacements_COI <- setNames(mapping_table_COI$NewName, mapping_table_COI$OldName)

# Function to replace sample names in a table while preserving specific columns
replace_sample_names <- function(table, name_replacements) {
  # Identify columns to preserve based on their position
  #first_col <- colnames(table)[1]  # First column
  second_last_col <- colnames(table)[length(colnames(table)) - 1]  # Second to last column
  last_col <- colnames(table)[length(colnames(table))]  # Last column
  
  # Rename columns
  new_colnames <- ifelse(colnames(table) %in% names(name_replacements),
                         name_replacements[colnames(table)],
                         colnames(table))
  
  # Ensure specific columns remain unchanged
  #new_colnames[1] <- first_col
  new_colnames[length(new_colnames) - 1] <- second_last_col
  new_colnames[length(new_colnames)] <- last_col
  
  colnames(table) <- new_colnames
  return(table)
}

# Apply the renaming function to each extended table
extended_tables_COI <- lapply(extended_tables_COI, function(table) {
  replace_sample_names(table, name_replacements_COI)
})

# Save the modified tables to new CSV files
output_paths_COI <- c(
  "extenedFinalTable_Apr21_COI_200924_renamed.csv",
  "extenedFinalTable_Sep20_COI_200924_renamed.csv",
  "extenedFinalTable_Aug23_COI_170924_renamed.csv"
)

mapply(function(table, path) {
  write.csv(table, path)
}, extended_tables_COI, output_paths_COI)

# Then put the outputs csv file in COI folder to merge them in the next step "3. COI_merge_tables"

#----------------18S-------------------------

#Set working directory
setwd("~/Datapaper2.sampleRenaming/18S")

# Read the extended tables into a list
file_paths <- c(
  "extenedFinalTable_Apr21_18S_pr2_0623_postdecontam.csv",
  "extenedFinalTable_Sep20_18S_pr2_0623_postdecontam.csv",
  "extenedFinalTable_Aug23_18S_pr2_0624_postdecontam.csv"
)

# Load all tables into a list
extended_tables_18S <- lapply(file_paths, function(path) {
  read.csv(path, row.names = 1)
})

# Read the mapping table
mapping_table_18S <- read.csv("mappingTable_18S.csv")
colnames(mapping_table_18S) <- c("NewName", "OldName")

# Create a named vector for replacement
name_replacements_18S <- setNames(mapping_table_18S$NewName, mapping_table_18S$OldName)

# Function to replace sample names in a table
# Function to replace sample names in a table while preserving specific columns
replace_sample_names <- function(table, name_replacements) {
  # Identify columns to preserve based on their position
  #first_col <- colnames(table)[1]  # First column
  second_last_col <- colnames(table)[length(colnames(table)) - 1]  # Second to last column
  last_col <- colnames(table)[length(colnames(table))]  # Last column
  
  # Rename columns
  new_colnames <- ifelse(colnames(table) %in% names(name_replacements),
                         name_replacements[colnames(table)],
                         colnames(table))
  
  # Ensure specific columns remain unchanged
  #new_colnames[1] <- first_col
  new_colnames[length(new_colnames) - 1] <- second_last_col
  new_colnames[length(new_colnames)] <- last_col
  
  colnames(table) <- new_colnames
  return(table)
}


# Apply the renaming function to each extended table
extended_tables_18S <- lapply(extended_tables_18S, function(table) {
  replace_sample_names(table, name_replacements_18S)
})

# Save the modified tables to new CSV files
output_paths <- c(
  "extenedFinalTable_Apr21_18S_pr2_0623_renamed.csv",
  "extenedFinalTable_Sep20_18S_pr2_0623_renamed.csv",
  "extenedFinalTable_Aug23_18S_pr2_0624_renamed.csv"
)

mapply(function(table, path) {
  write.csv(table, path)
}, extended_tables_18S, output_paths)

# Then put the outputs csv file in 18S folder to merge them in the next step "3. 18S_merge_tables"

#----------------ITS-------------------------

#Set working directory
setwd("~/Datapaper/2.sampleRenaming/ITS")

# Read the extended tables into a list
file_paths <- c(
  "extenedFinalTable_Apr21_ITS_200924_postdecontam.csv",
  "extenedFinalTable_Sep20_ITS_210924_postdecontam.csv"
)

# Load all tables into a list
extended_tables_ITS <- lapply(file_paths, function(path) {
  read.csv(path, row.names = 1)
})

# Read the mapping table
mapping_table_ITS <- read.csv("mappingTable_ITS.csv")
colnames(mapping_table_ITS) <- c("NewName", "OldName")

# Create a named vector for replacement
name_replacements <- setNames(mapping_table_ITS$NewName, mapping_table_ITS$OldName)

# Function to replace sample names in a table while preserving specific columns
replace_sample_names <- function(table, name_replacements) {
  # Identify columns to preserve based on their position
  #first_col <- colnames(table)[1]  # First column
  second_last_col <- colnames(table)[length(colnames(table)) - 1]  # Second to last column
  last_col <- colnames(table)[length(colnames(table))]  # Last column
  
  # Rename columns
  new_colnames <- ifelse(colnames(table) %in% names(name_replacements),
                         name_replacements[colnames(table)],
                         colnames(table))
  
  # Ensure specific columns remain unchanged
  #new_colnames[1] <- first_col
  new_colnames[length(new_colnames) - 1] <- second_last_col
  new_colnames[length(new_colnames)] <- last_col
  
  colnames(table) <- new_colnames
  return(table)
}

# Apply the renaming function to each extended table
extended_tables_ITS <- lapply(extended_tables_ITS, function(table) {
  replace_sample_names(table, name_replacements)
})

# Save the modified tables to new CSV files
output_paths <- c(
  "extenedFinalTable_Apr21_ITS_200924_renamed.csv",
  "extenedFinalTable_Sep20_ITS_210924_renamed.csv"
)

mapply(function(table, path) {
  write.csv(table, path)
}, extended_tables_ITS, output_paths)

# 
# Copy the outputs csv file in ITS folder to merge them in the next step "3. ITS_merge_tables"