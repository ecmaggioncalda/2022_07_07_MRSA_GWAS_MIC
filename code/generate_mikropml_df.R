#mikropml input file generation
source("code/log_smk.R") #this assigns the log file for the run
#LIBRARIES ----
library(tidyverse)
library(doFuture)

#These two lines assign proper variables for using the cluster
doFuture::registerDoFuture()
future::plan(future::multicore, workers = snakemake@resources[["ncores"]])

# FUNCTION ----
# return vector of non-variant sites
index_of_non_var_sites <- function(mat){
  
  index = as.numeric(c(which(rowSums(mat) == 0), which(rowSums(mat) == ncol(mat))))
  return(index)
  
}

#PHENO ----
paste0(snakemake@input[["pheno"]])
pheno <- readr::read_delim(file = snakemake@input[["pheno"]],
                           delim = "\t")

pheno_merge <- as.data.frame(pheno)
rownames(pheno_merge) <- pheno_merge$genome_id
pheno_merge <- pheno_merge[, -1, drop = FALSE]

#GENO ----
paste0(snakemake@params[['path']])
geno <- read.table(snakemake@params[['path']],
                   sep = "\t",
                   stringsAsFactors = FALSE,
                   header = TRUE,
                   row.names = 1)

print("Geno matrix has completed read in")

colnames(geno) <- gsub("$", "_", colnames(geno))

if(any(length(intersect(colnames(geno), rownames(pheno_merge))) != ncol(geno) | length(intersect(colnames(geno), rownames(pheno_merge))) != nrow(pheno_merge))){
  print("Mistmatch in geno/pheno contents: subsetting contents")
  geno_diff <- setdiff(colnames(geno), rownames(pheno_merge))
  pheno_diff <- setdiff(rownames(pheno_merge), colnames(geno))
  non_overlap <- c(geno_diff, pheno_diff)
  geno <- geno[, !colnames(geno) %in% geno_diff, drop = FALSE]
  pheno_merge <- pheno_merge[!rownames(pheno_merge) %in% pheno_diff, , drop = FALSE]
  print("The following strains were removed as non-overlapping contents")
  print(paste0(non_overlap))
  print("overwriting previous phenotype file")
  pheno_overwrite <- pheno_merge %>%
    rownames_to_column(var = "genome_id")
  write_tsv(pheno_overwrite,
            file = snakemake@input[["pheno"]])
}

# geno <- geno[,colnames(geno) %in% rownames(pheno_merge)]
rownames(geno) <- make.names(rownames(geno), unique = TRUE)

# REMOVE VARS
to_remove <- index_of_non_var_sites(geno)
if(length(to_remove) >= 1){
  geno_subset <- geno[-to_remove,]
  geno_merge <- t(geno_subset)
}else{
  geno_merge <- t(geno)
}

# geno_distinct <- distinct(geno_subset)
# temp_geno <- as.data.frame(geno_distinct)
# temp_geno <- as.data.frame(geno_subset)
# temp_geno <- temp_geno %>% 
#   rownames_to_column('variant')
# geno_merge <- geno_merge[rownames(geno_merge) %in% rownames(pheno_merge), ]

if(sum(rownames(pheno_merge) %in% rownames(geno_merge)) != length(rownames(pheno_merge))){
  stop("mismatch between pheno and geno contents")
}

index <- match(rownames(pheno_merge), rownames(geno_merge))

geno_ordered <- geno_merge[index, , drop = FALSE]

if(sum(rownames(pheno_merge) == rownames(geno_ordered)) != length(rownames(pheno_merge))){
  stop("mismatch between pheno and geno contents")
}

complete_frame <- cbind(pheno_merge,
                        geno_ordered)

print("complete frame generated with pheno:geno, export to csv for mikropml preprocessing")

#GENERATE FILES ----
#patient and genome factors
write_csv(complete_frame,
          file = snakemake@output[['file_name']],
          col_names = TRUE)
