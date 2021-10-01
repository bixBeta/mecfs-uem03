source("scripts/packages.R")

cfs.all = readRDS("/workdir/docker-dump-uem03/global_data/CFSall.RDS")
cfs.all.meta = cfs.all@meta.data

df = select(cfs.all.meta, orig.ident, percent.mt, batch10x)



png(filename = "figures/ggViolin__percentMT.png", width = 30000,height = 2000)
ggplot(df, aes(x= orig.ident, y=percent.mt, fill=orig.ident)) + 
  geom_violin() + geom_jitter(aes(colour = orig.ident)) + theme(axis.text=element_text(size=12),
     axis.text.x =element_text(size=30,face="bold")) + theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold")) + scale_y_continuous(breaks = seq(min(df$percent.mt), max(df$percent.mt), by =5)) +  geom_hline(yintercept = 30) 
dev.off()

df.drop.4759.4769 = df %>% filter(orig.ident != "R4759") %>% filter(orig.ident != "R4769")
cfs.samples.to.use = levels(as.factor(df.drop.4759.4769$orig.ident))
