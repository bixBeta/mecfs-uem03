
.libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.0")

library("Seurat")
library("glmGamPoi")

sobj.filtered = readRDS("/home/rstudio/rstudio-docker-sessions/01__mecfs/data/post__cellCycleScoring.RDS")

sobj.list <- SplitObject(sobj.filtered, split.by = "batch10x")

sobj.list

sobj.list <- lapply(X = sobj.list, FUN = function(x){SCTransform(object = x, method = "glmGamPoi", vars.to.regress = c("percent.mt","S.Score", "G2M.Score"))})
	
features <- SelectIntegrationFeatures(sobj.list, nfeatures = 2000) #may need to drop nfeatures
 
sobj.list <- PrepSCTIntegration(object.list = sobj.list, anchor.features = features)

save.image(file = "/home/rstudio/rstudio-docker-sessions/01__mecfs/data/03c__pre__integration__image.Rdata")

anchors <- FindIntegrationAnchors(object.list = sobj.list, normalization.method = "SCT",  anchor.features = features)

saveRDS(anchors,  file = "/home/rstudio/rstudio-docker-sessions/01__mecfs/data/03c__anchors__object.RDS")

sobj.SCT.CCA.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

save.image(file = "/home/rstudio/rstudio-docker-sessions/01__mecfs/data/03c__post__integration__image.Rdata")

saveRDS(sobj.SCT.CCA.integrated,  file = "/home/rstudio/rstudio-docker-sessions/01__mecfs/data/03c__sct__integrated__3000__features.RDS")
