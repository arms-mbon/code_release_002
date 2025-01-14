# code_release_002

This repository contains all scripts and code used for the analyses presented in **ARMS-MBON's Second Data Paper**. For related datasets, please refer to:  
- [data_release_002](https://github.com/arms-mbon/data_release_002): Occurrence and event data.  
- [analysis_release_002](https://github.com/arms-mbon/analysis_release_002): Bioinformatics pipeline outputs.

---

## Table of Contents

1. [decontam](#1-decontam)
2. [rename_samples](#2-rename_samples)
3. [merge_tables](#3-merge_tables)
   - [COI](#coi_merge_tables)
   - [18S](#18s_merge_tables)
   - [ITS](#its_merge_tables)
4. [gene_analysis](#4-gene_analysis)

---

### 1. decontam.R

**Description:**  
This script performs blank curation as the first step of data processing. Using the *prevalence method* from the `decontam` R package, it identifies and removes potential contaminants from the dataset.

---

### 2. rename_samples.R

**Description:**  
The second step involves renaming sample identifiers. PEMA outputs use ENA accession numbers as sample names, which are replaced with their corresponding material sample IDs for clarity and consistency.

---

### 3. merge_tables

**Description:**  
This step merges data from the PEMA outputs, including read count tables, taxonomy assignments, and FASTA files, for each genetic marker. Separate scripts handle each marker:  

#### 3. COI_merge_tables.R 
Processes data for the COI gene.

#### 3. 18S_merge_tables.R  
Processes data for the 18S gene.

#### 3. ITS_merge_tables.R  
Processes data for the ITS gene.

---

### 4. gene_analysis.R

**Description:**  
The final step involves exploratory data analysis and visualization, including:  

- Curation of merged datasets.  
- Assessment of sequencing depth.  
- Visualization of recovered phyla and species.  
- Creation of an UpSet plot to show the overlap in species identified across marker gene datasets.  
- Comparisons between datasets from the first data paper (DP001) and the second data paper (DP002).

---

By providing this comprehensive code and documentation, we aim to ensure transparency and reproducibility for all analyses conducted in the ARMS-MBON project.
