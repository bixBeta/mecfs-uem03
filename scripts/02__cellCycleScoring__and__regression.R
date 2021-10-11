sobj.filtered = readRDS("data/cfs__116__samples__filtered.RDS")
load("/workdir/docker-dump-uem03/global_data/cycle.rda")

sobj.filtered.meta = sobj.filtered@meta.data
sobj.filtered = NormalizeData(object = sobj.filtered)

sobj.filtered <- FindVariableFeatures(sobj.filtered, selection.method = "vst")
sobj.filtered <- ScaleData(sobj.filtered, features = rownames(sobj.filtered))
sobj.filtered <- RunPCA(sobj.filtered, features = VariableFeatures(sobj.filtered), ndims.print = 6:10, nfeatures.print = 10, npcs = 50)


# Score cells for cell cycle
sobj.filtered <- CellCycleScoring(sobj.filtered, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

