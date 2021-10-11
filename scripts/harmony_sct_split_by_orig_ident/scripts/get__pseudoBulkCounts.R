library(Matrix.utils)
library(SingleCellExperiment)
source("scripts/packages.R")
# pick umap res 0.6 
Idents(object = harm.sct) <- "SCT_snn_res.0.6"
# extract pseudoBulkCounts
sce <- as.SingleCellExperiment(harm.sct)
assays(sce)

groups <- colData(sce)[, c("SCT_snn_res.0.6", "orig.ident")]
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2),`[`, 1)
pb <- split.data.frame(pb, 
  factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
  stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

pb.matrix = lapply(pb, function(x){as.matrix(x)})

save(pb,pb.matrix, file = "scripts/harmony_sct_split_by_orig_ident/sct__harmony__res__0.6__PseudoBulk_Counts.Rdata")

#----------------------------------------------------------------------------------------
# pick umap res 0.8
Idents(object = harm.sct) <- "SCT_snn_res.0.8"
# extract pseudoBulkCounts
sce <- as.SingleCellExperiment(harm.sct)
assays(sce)

groups <- colData(sce)[, c("SCT_snn_res.0.8", "orig.ident")]
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2),`[`, 1)
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

pb.matrix = lapply(pb, function(x){as.matrix(x)})

save(pb,pb.matrix, file = "scripts/harmony_sct_split_by_orig_ident/sct__harmony__res__0.8__PseudoBulk_Counts.Rdata")