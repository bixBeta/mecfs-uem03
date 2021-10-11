.libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.0")

library("Seurat")
library("glmGamPoi")
library("harmony")

sobj = readRDS("/home/rstudio/rstudio-docker-sessions/01__mecfs/data/post__cellCycleScoring.RDS")

sobj.list <- SplitObject(sobj, split.by="batch10x")

sobj.list <- lapply(X = sobj.list, FUN = function(x){SCTransform(object = x, method = "glmGamPoi", vars.to.regress = "percent.mt", return.only.var.genes = FALSE)})

var.features <- SelectIntegrationFeatures(object.list = sobj.list, nfeatures = 3000)

sobj.sct <- merge(x = sobj.list[[1]], y = sobj.list[2:length(sobj.list)], merge.data=TRUE)

VariableFeatures(sobj.sct) <- var.features

sobj.sct <- RunPCA(sobj.sct, verbose = FALSE, dims = 1:50)

sobj.sct <- RunHarmony(sobj.sct, assay.use="SCT", group.by.vars = "batch10x")

sobj.sct <- RunUMAP(sobj.sct, reduction = "harmony", dims = 1:50)

sobj.sct <- FindNeighbors(sobj.sct, reduction = "harmony", dims = 1:50) %>% FindClusters(resolution = c(0.4, 0.5, 0.6, 0.8))


saveRDS(sobj.sct, file = "sct__harmony__splitByBatch10x__ndim50__res0.4_0.5_0.6_0.8.RDS")

save.image("sct__harmony__split__by__batch10x.Rdata")
