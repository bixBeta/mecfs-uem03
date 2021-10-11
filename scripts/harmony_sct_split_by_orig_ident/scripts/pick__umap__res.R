library(Seurat)
library(patchwork)
library(ggplot2)

harm.sct <- readRDS("~/rstudio-docker-sessions/01__mecfs/scripts/harmony_sct_split_by_orig_ident/sct__harmony__ndim50__res0.4_0.5_0.6_0.8.RDS")
harm.sct.meta = harm.sct@meta.data

Idents(object = harm.sct) <- "SCT_snn_res.0.4"
umap.0.4 = DimPlot(harm.sct, reduction = "umap", label = TRUE, label.size = 6, raster = F, group.by = "SCT_snn_res.0.4")

Idents(object = harm.sct) <- "SCT_snn_res.0.5"
umap.0.5 = DimPlot(harm.sct, reduction = "umap", label = TRUE, label.size = 6, raster = F, group.by = "SCT_snn_res.0.5")

Idents(object = harm.sct) <- "SCT_snn_res.0.6"
umap.0.6 = DimPlot(harm.sct, reduction = "umap", label = TRUE, label.size = 6, raster = F, group.by = "SCT_snn_res.0.6")

Idents(object = harm.sct) <- "SCT_snn_res.0.8"
umap.0.8 = DimPlot(harm.sct, reduction = "umap", label = TRUE, label.size = 6, raster = F, group.by = "SCT_snn_res.0.8")

png(filename = "umap.sct.harm.png",width = 1920, height = 1080)
umap.0.4 + umap.0.5 + umap.0.6 + umap.0.8
dev.off()


DefaultAssay(harm.sct) <- "RNA"
DotPlot(harm.sct,
        features=c("CD3E","CD247","CD8A","CCR7","CD28","IFNG","GATA3","KLRB1","CTLA4","TRBC1","TRGC1","TRDC","KLRK1","NCR3","IRF4","CD4","LILRB1","SPI1","AHR","CD14","LGALS2","FCGR3A","CCL5","PF4","MS4A1"),
        cluster.idents=T, scale.by="size") + RotatedAxis() + ggtitle("RNA-assay-0.8")


DefaultAssay(harm.sct) <- "SCT"
DotPlot(harm.sct,
        features=c("CD3E","CD247","CD8A","CCR7","CD28","IFNG","GATA3","KLRB1","CTLA4","TRBC1","TRGC1","TRDC","KLRK1","NCR3","IRF4","CD4","LILRB1","SPI1","AHR","CD14","LGALS2","FCGR3A","CCL5","PF4","MS4A1"),
        cluster.idents=T, scale.by="size") + RotatedAxis() + ggtitle("SCT-assay-0.8")


harm.sct.meta = harm.sct@meta.data
