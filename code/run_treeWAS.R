#run treeWAS
source("code/log_smk.R")
#LIBRARIES ----
library(ape)
library(tidyverse)
library(treeWAS)
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
geno <- as.data.frame(geno)
rownames(geno) <- pheno$genome_id

pheno <- deframe(pheno)

# TREE ----
tree <- read.tree(snakemake@params[["tree"]])

tree <- drop.tip(tree, setdiff(tree$tip.label, names(pheno)))

#PRINT STATEMENTS ----
print("dimension of geno, nrows (variants), ncol (genomes)")
dim(geno)

if(sum(rownames(geno) %in% names(pheno)) != nrow(geno) | sum(rownames(geno) %in% tree$tip.label) != nrow(geno) | sum(names(pheno) %in% tree$tip.label) != nrow(geno)){
  stop("Mismatch in geno/pheno/tree contents")
}

# TREEWAS ----
treeWAS.out <- treeWAS(snps = geno,
                       phen = pheno,
                       tree = tree,
                       filename.plot = snakemake@output[["plot"]])

save(treeWAS.out, file = snakemake@output[["rdata"]])