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

sobj = subset(x = cfs.all, idents = cfs.samples.to.use)
saveRDS(sobj, file = "data/cfs_dropped_R4759_R4769_unfiltered.RDS")

sobj.meta = sobj@meta.data

png(filename = "figures/ggViolin__DroppedOutlierSamples__percentMT.png", width = 30000,height = 2000)

ggplot(sobj.meta, aes(x= orig.ident, y=percent.mt, fill=orig.ident)) + 
  geom_violin() + geom_jitter(aes(colour = orig.ident)) + theme(axis.text=element_text(size=12),
   axis.text.x =element_text(size=30,face="bold")) + theme(axis.text=element_text(size=12), 
  axis.title=element_text(size=14,face="bold")) + scale_y_continuous(breaks = seq(min(sobj.meta$percent.mt), max(sobj.meta$percent.mt), by =5)) +  geom_hline(yintercept = 30) 

dev.off()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Add number of genes per UMI for each cell to metadata
sobj$log10GenesPerUMI <- log10(sobj$nFeature_RNA) / log10(sobj$nCount_RNA)
sobj.meta = sobj@meta.data


# Visualize the number UMIs/transcripts per cell
sobj.meta %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)


# Visualize the distribution of genes detected per cell via histogram
sobj.meta %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via boxplot
sobj.meta %>% 
  ggplot(aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI

sobj.meta %>%
  ggplot(aes(x=log10GenesPerUMI, color = phenoGroup, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) + facet_wrap(~ENID)

