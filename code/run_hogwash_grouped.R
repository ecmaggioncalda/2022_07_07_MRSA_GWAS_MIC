#Script to run hogwash grouped
source("code/log_smk.R")
#LIBRARIES ----
library(ape)
library(tidyverse)
library(hogwash)
library(doFuture)

#These two lines assign proper variables for using the cluster
doFuture::registerDoFuture()
future::plan(future::multicore, workers = snakemake@resources[["ncores"]])

#USE MIKROPML RDS ----
rds <- readRDS(snakemake@input[["rds"]])$dat_transformed

#PHENO ----
paste0(snakemake@input[["pheno"]])

pheno <- readr::read_delim(file = snakemake@input[["pheno"]],
                           delim = "\t")

#GENO ----
geno <- rds[, !grepl(snakemake@wildcards[['phenotype']], colnames(rds))]

#CLEAN GENO & PHENO ----
pheno <- as.data.frame(pheno)
rownames(pheno) <- pheno$genome_id

geno <- as.data.frame(geno)
rownames(geno) <- pheno$genome_id

pheno <- pheno[, -1, drop = FALSE]

# TREE ----
tree <- read.tree(snakemake@params[["tree"]])

tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(pheno)))

#ORDER DATA ----
reorder_pheno <- match(tree$tip.label, rownames(pheno))

pheno <- as.matrix(pheno[reorder_pheno, , drop = FALSE])

reorder_geno <- match(tree$tip.label, rownames(geno))

geno <- as.matrix(geno[reorder_geno,  , drop = FALSE])

if(all(rownames(pheno) == tree$tip.label) | all(rownames(geno) == tree$tip.label)){
  
  print("Geno/Pheno/Tree contents are in the correct order")
  
}else{
  
  stop("Mismatch in geno/pheno/tree contents")
  
}

'Number of samples and variants'
dim(geno)

#LOAD GENE KEY ----
gene_key <- snakemake@params[["gene_key"]]

#RUN HOGWASH ----
hogwash(pheno = pheno, 
        geno = geno, 
        tree = tree, 
        file_name = snakemake@params[["file_name"]],
        dir = snakemake@params[["dir"]],
        group_genotype_key = gene_key,
        grouping_method = "post-ar",
        test='continuous')