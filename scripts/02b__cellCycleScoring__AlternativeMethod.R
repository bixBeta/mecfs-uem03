sobj.filtered = readRDS("data/post__cellCycleScoring.RDS")
sobj.filtered.meta = sobj.filtered@meta.data

sobj.filtered$CC.Difference <- sobj.filtered$S.Score - sobj.filtered$G2M.Score
sobj.filtered.meta = sobj.filtered@meta.data

saveRDS(sobj.filtered, "data/alternativ__method__post__cellCycleScoring.RDS")
