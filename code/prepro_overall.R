# Builds individual phenotype file and generates necessary directories for snakemake workflow
# This file would need to be updated based on the measures that the investigator is taking to adjust the continuous measure
# In this example, IDSA_severity is split into multiple separate data frames for analysis grouped by elix weighted score
source("code/log_smk.R") #this assigns the log file for the run
#LIBRARIES ----
library(tidyverse)
library(doFuture)

#These two lines assign proper variables for using the cluster
doFuture::registerDoFuture()
future::plan(future::multicore, workers = snakemake@resources[["ncores"]])

#FILES ----
pheno <- readr::read_csv(file = snakemake@input[["metadata"]]) %>%
  rename("genome_id" ="Genome ID") %>%
  mutate(genome_id = gsub("$", "_", genome_id))
phenotype_cols <- snakemake@wildcards[["phenotype"]]

paths <- paste0(c("data/pheno/",
                  "results/",
                  "benchmarks/",
                  "figures/",
                  "log/",
                  "data/mikropml/"),
                phenotype_cols)
paths <- c("data/pheno", "data/mikropml", paths)

for(i in 1:length(paths)){
  
  if(dir.exists(paths[i]) == FALSE){
    
    dir.create(paths[i])
    
  }else{
    
    print(paste0("directory ", paths[i], " already exists"))
    
  }
}

new_dir3 <- paste0("results/", phenotype_cols, "/runs")

if(dir.exists(new_dir3) == FALSE){
  
  dir.create(new_dir3)
    
  }else{
    
    print(paste0("directory ", new_dir3, " already exists"))
    
  }

if(!all(sapply(paths, dir.exists))){
  stop("Not all directories created")
}

#UNADJUSTED FILE ----
unadjusted_out <- pheno %>%
  filter(grepl("Kleb", Organism_ID)) %>%
  select(genome_id,
         all_of(phenotype_cols))

log2_convert <- unadjusted_out %>%
  mutate(across(all_of(phenotype_cols), ~log2(.x)))
  
write_tsv(log2_convert,
          file = paste0("data/pheno/", phenotype_cols, "/log2.tsv"))
