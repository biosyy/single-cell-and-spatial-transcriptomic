rm(list = ls())
options(stringsAsFactors = F)
options("repos" = c(CRAN="http://mirrors.ustc.edu.cn/CRAN/"))
Sys.setenv(LANGUAGE= "en")
library(tidyverse)
library(Seurat)
library(patchwork)
library(parallel)
library(ggthemes)
library(magick)
library(ggsignif)
library(AUCell)
load('../share_data/sce.combined.sct.Rda')
sce.fibro <- subset(sce.combined.sct,subset = cluster == 'Fibroblasts')
save(sce.fibro,file = 'sce.fibro.Rda')
load('sce.fibro.Rda')
DimPlot(sce.fibro)
####cell type anotation
DefaultAssay(sce.fibro) <- 'integrated'
sce.fibro <- sce.fibro %>%  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = 0.4)
sce.fibro$cluster <- as.character(sce.fibro$seurat_clusters)
sce.fibro$cluster[sce.fibro$cluster == '6'] <- '1'
sce.fibro@active.ident <- as.factor(sce.fibro$cluster)
DimPlot(sce.fibro)
DefaultAssay(sce.fibro) <- 'RNA'
markers <- FindAllMarkers(sce.fibro,logfc.threshold = 0.25,only.pos = T,min.pct = 0.1)
top30 <- markers %>% group_by(cluster) %>% top_n(30,wt = avg_log2FC)
FeaturePlot(sce.fibro,features = 'ACTG2')

p <- DimPlot(sce.fibro,label = T)
pdf('test.pdf',width = 7,height = 5,onefile = FALSE)
print(p)
dev.off()

FeaturePlot(sce.fibro,features = 'POSTN')
sce.fibro$cluster[sce.fibro$cluster == '0'] <- 'POSTN_fibro'
FeaturePlot(sce.fibro,features = 'CFD')
sce.fibro$cluster[sce.fibro$cluster == '1'] <- 'CFD_fibro'
FeaturePlot(sce.fibro,features = 'SORBS2')
sce.fibro$cluster[sce.fibro$cluster == '2'] <- 'SORBS2_fibro'
FeaturePlot(sce.fibro,features = 'CD36')
sce.fibro$cluster[sce.fibro$cluster == '3'] <- 'CD36_fibro'
FeaturePlot(sce.fibro,features = 'ACTG2')
sce.fibro$cluster[sce.fibro$cluster == '4'] <- 'ACTG2_fibro'
FeaturePlot(sce.fibro,features = 'GFRA3')
sce.fibro$cluster[sce.fibro$cluster == '5'] <- 'GFRA3_fibro'
FeaturePlot(sce.fibro,features = 'KIT')
sce.fibro$cluster[sce.fibro$cluster == '7'] <- 'KIT_fibro'
FeaturePlot(sce.fibro,features = 'PRKG1')
sce.fibro$cluster[sce.fibro$cluster == '8'] <- 'PRKG1_fibro'



DimPlot(sce.fibro,group.by = 'cluster',label = T)
id <- CellSelector(p)
sce.fibro <- sce.fibro[,id]
p <- DimPlot(sce.fibro,group.by = 'cluster',label = T)
pdf('./Fig/fibro_subcluster_umap.pdf',width = 7,height = 5,onefile = FALSE)
print(p)
dev.off()

p <- DimPlot(sce.fibro,group.by = 'sample',label = T)
pdf('./Fig/fibro_sample_umap.pdf',width = 7,height = 5,onefile = FALSE)
print(p)
dev.off()
sce.fibro$group <- data$group

p <- DimPlot(sce.fibro,group.by = 'group',label = T)
pdf('./Fig/fibro_group_umap.pdf',width = 7,height = 5,onefile = FALSE)
print(p)
dev.off()


###cell percent
library(ggthemes)
library(scales)
p
x<-ggplot_build(p)
info = data.frame(colour = x$data[[1]]$colour, group = x$data[[1]]$group)
info <- unique((arrange(info, group)))
cols <- as.character(info$colour)



data <- FetchData(sce.fibro,vars = c('sample','cluster'))
data$group <- if_else(data$sample%in%c('N1'),'normal',
                      if_else(data$sample%in%c('T1','T3'),'chemo','pri'))

data1 <- data %>% group_by(sample,cluster) %>% summarise(counts = n())
data1.perc <- data1 %>% group_by(sample) %>% summarise(cluster = cluster,percent = counts/sum(counts))

p <- ggplot(data1.perc,mapping = aes(sample,percent,fill=cluster))+
  geom_bar(stat='identity',position='fill') +
  scale_fill_manual("cluster", values = cols) +
  labs(x = 'Samples',y = 'Proportion') +
  theme_few()

pdf('./Fig/fibro_sample_percent.pdf',width = 7,height = 5,onefile = FALSE)
print(p)
dev.off()
p
data2 <- data %>% group_by(group,cluster) %>% summarise(counts = n())
data2.perc <- data2 %>% group_by(group) %>% summarise(cluster = cluster,percent = counts/sum(counts))

p <- ggplot(data2.perc,mapping = aes(cluster,percent,fill=group))+
  geom_bar(stat='identity',position='fill') +
  scale_fill_few() +
  labs(x = 'cluster',y = 'Proportion') +
  theme_few()
p
pdf('./Fig/fibro_group_percent.pdf',width = 7,height = 5,onefile = FALSE)
print(p)
dev.off()

####signature heatmap
FeaturePlot(sce.fibro,features = 'ZSWIM6')
featuers <- c('ACTG2','DES','SMTN','MYLK','CNN1',
              'RGS5','CD36','HIGD1B','ADGRF5','ARHGDIB',
              'CFD','PLA2G2A','C7','ADH1B','C3',
              'PLP1','GFRA3','CDH19','S100B','NRXN1',
              'KIT','IQGAP2','DPP10','ANXA3','EDN3',
              'POSTN','COL1A1','COL3A1','INHBA','COL5A2',
              'FTX','PRKG1','GLIS3','SDK1','ZSWIM6',
              'RERGL','SORBS2','LBH','RCAN2','SNCG')

exp <- sce.fibro@assays$RNA@data[featuers,] %>% as.matrix() %>% t()%>% as.data.frame()
exp$celltype <- sce.fibro$cluster
mean.exp <- exp %>% group_by(celltype) %>% summarise(across(.cols = everything(),.fns = mean)) %>% 
  column_to_rownames('celltype')
library(pheatmap)
library(RColorBrewer)
display.brewer.all() 
color <- brewer.pal(11, "RdBu")[2:10]

bk = unique(c(seq(-2,2, length=100)))
p <- pheatmap::pheatmap(mean.exp,cluster_rows = F,cluster_cols = F,scale = 'column',border =F,
                   breaks = bk,
                   color = colorRampPalette(rev(color))(100))
pdf('./Fig/fibro_subtypes_signatures_heatmap.pdf',width = 8,height = 4,onefile = FALSE)
print(p)
dev.off()



#GO Ontology
markers.sig <- markers %>% filter(avg_log2FC>0.5)
library(clusterProfiler)
library(org.Hs.eg.db)
library(viridis)
markers.sig$cluster %>% table()
gene.df <- markers.sig %>% dplyr::filter(cluster == 'SORBS2_fibro') %>% dplyr::select('gene')
gene.df <- bitr(gene.df$gene, fromType = "SYMBOL",
                toType ='ENTREZID',
                OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene          = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
ACTG2_fibro.GO <- ego %>% as.data.frame()
CD36_fibro.GO <- ego %>% as.data.frame()
CFD_fibro.GO <- ego %>% as.data.frame()
GFRA3_fibro.GO <- ego %>% as.data.frame()
KIT_fibro.GO <- ego %>% as.data.frame()
POSTN_fibro.GO <- ego %>% as.data.frame()
PRKG1_fibro.GO <- ego %>% as.data.frame()
SORBS2_fibro.GO <- ego %>% as.data.frame()

DimPlot(sce.fibro,group.by = 'cluster',label = T)
markers.sig$gene %>% unique() %>% length()

table(markers.sig$gene)

###mCAF:ACTG2_fibro,CD36_fibro,SORBS2_fibro
###iCAF:CFD_fibro,高表达补体因子，调节炎症。
###nuronCAF:GFRA3_fibro, KIT_fibro,PRKG1_fibro,
###classicalCAF: POSTN_fibro
###POSTN_fibro在癌组织中丰度很高，提示重要促癌类型。
FeaturePlot(sce.fibro,features = c('DCN', 'OGN', 'ACTB','PDGFRA','PDGFRB','DDR2','FAP','CAV1'))

FeaturePlot(sce.fibro,features = c('FAP'))
gene.df <- markers.sig %>% dplyr::select(cluster,gene)
gene.id <- bitr(gene.df$gene, fromType = "SYMBOL",
                toType ='ENTREZID',
                OrgDb = org.Hs.eg.db)
gene.df <- merge(gene.id,gene.df,by.x = 'SYMBOL',by.y = 'gene')

ego <- compareCluster(ENTREZID~cluster,data = gene.df,fun = 'enrichGO',
                      OrgDb = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
go.df <- as.data.frame(ego)







###mCAF
rownames(ACTG2_fibro.GO) <- seq(254)
ACTG2_fibro.GO.se <- ACTG2_fibro.GO[c(2,3,30,35,141,143),]
rownames(CD36_fibro.GO) <- seq(539)
CD36_fibro.GO.se <- CD36_fibro.GO[c(13,14,89,54,4,5),]

rownames(SORBS2_fibro.GO) <- seq(231)
SORBS2_fibro.GO.se <- SORBS2_fibro.GO[c(1,2,3,4,71,72),]
###iCAF
rownames(CFD_fibro.GO) <- seq(824)
CFD_fibro.GO.se <- CFD_fibro.GO[c(104,95,25,5,1,2),]
####nuronCAF
rownames(GFRA3_fibro.GO) <- seq(56)
GFRA3_fibro.GO.se <- GFRA3_fibro.GO[c(1,3,26,27),]
rownames(KIT_fibro.GO) <- seq(204)
KIT_fibro.GO.se <- KIT_fibro.GO[c(1,2,4,8,63,65),]
rownames(PRKG1_fibro.GO) <- seq(756)
PRKG1_fibro.GO.se <- PRKG1_fibro.GO[c(8,85,47,50,5,6),]
####classicaCAF
rownames(POSTN_fibro.GO) <- seq(399)
POSTN_fibro.GO.se <- POSTN_fibro.GO[c(26,28,16,17,1,2),]
combine.GO <- c(ACTG2_fibro.GO.se$Description,CD36_fibro.GO.se$Description,
                SORBS2_fibro.GO.se$Description,CFD_fibro.GO.se$Description,
                GFRA3_fibro.GO.se$Description,KIT_fibro.GO.se$Description,
                PRKG1_fibro.GO.se$Description,POSTN_fibro.GO.se$Description) %>% unique()
combine.GO %>% length()
data <- go.df %>% filter(Description%in%combine.GO)

order.item <- c('oxidative phosphorylation',
  'energy derivation by oxidation of organic compounds',
  'regulation of blood circulation',
  'vascular process in circulatory system',
  'muscle contraction',
  'muscle system process',
  'muscle tissue development',
  'muscle cell differentiation',
  'regulation of inflammatory response',
  'myeloid leukocyte migration',
  'complement activation',
  'axon development',
  'regulation of nervous system development',
  'regulation of neuron projection development',
  'epithelial cell migration',
  'epithelial cell development',
  'epithelial cell proliferation',
  'neutrophil chemotaxis',
  'regulation of angiogenesis',
  'regulation of vasculature development',
  'extracellular matrix organization',
  'extracellular structure organization'
  
)
order.cluster <- c('ACTG2_fibro','CD36_fibro','SORBS2_fibro','CFD_fibro',
                   'GFRA3_fibro','KIT_fibro',
                   'PRKG1_fibro','POSTN_fibro')
data$Description <- factor(data$Description,levels = order.item)
data$cluster <- factor(data$cluster,levels = order.cluster)
data$p.adjust
library(RColorBrewer)
display.brewer.all()
color1 <- brewer.pal(11, "RdBu")[-6]
data$pvalue1 <- -1*log10(data$pvalue)
data$pvalue1[data$pvalue1>6] <- 6
expression(data$pvalue1)

p = ggplot(data,aes(x=cluster,y=Description)) + 
  geom_point(aes(size=Count,color=pvalue1))+
  scale_color_gradient2(midpoint = 3)+
  #scale_size_continuous(range = c(4,12))+
  labs(
    color=expression(pvalue1),
    size="Gene number",
    #x="group"
    # y="Pathway name",
    # title="Pathway enrichment")
  )+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
p
pdf('./Fig/Fibro_GO_dotplot.pdf',width = 7,height = 5,onefile = FALSE)
print(p)
dev.off()

###Heatmap
data1 <- data %>% dplyr::select(cluster,Description,pvalue1)
data1.wider <- data1 %>% pivot_wider(names_from = Description,values_from = pvalue1) %>% column_to_rownames('cluster')

data1.wider <- data1.wider[order.cluster,order.item]

data1.wider[is.na(data1.wider)] <- 0
pheatmap::pheatmap(t(data1.wider),cluster_rows = F,cluster_cols = F,scale = 'none')


FeaturePlot(sce.fibro,features = 'NRXN3')

###空间共定位(未分析到阳性结果),考虑做免疫荧光炎症。
###NCAM1，NRXN1,NRXN3,SEMA3B,SEMA3C
###COL28A1,COL9A3.
genes <- c('NCAM1','NRXN1','NRXN3','SEMA3B','SEMA3C','COL28A1','COL9A3')
save(sce.fibro,file = 'sce.fibro.Rda')
sptial.T2 <- Load10X_Spatial(data.dir = '../ourdata/spatial/T2-1/',
                             filename = 'filtered_feature_bc_matrix.h5',
                             assay = 'Spatial',
                             slice = 'B2',
                             filter.matrix = T,
                             to.upper = F,
                             image = NULL)
sptial.T2 <- sptial.T2 %>% 
  SCTransform(assay = "Spatial", verbose = FALSE) %>% 
  RunPCA(assay = "SCT",verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = 0.5, verbose = FALSE)
SpatialDimPlot(sptial.T2, label = TRUE, label.size = 3) 

DefaultAssay(sptial.T2) <- 'Spatial'
sptial.T2 <- NormalizeData(sptial.T2)
exp <- sptial.T2@assays$Spatial@data[genes,] %>% as.matrix() %>% as.data.frame()
exp <- t(exp) %>% as.data.frame()
cor.m <- cor(exp)

SpatialFeaturePlot(sptial.T2,features =c('NCAM1','COL28A1') )



###clinical target expression in CAF
target.id <- c('DHRS3','VDR','HAS2','CTGF','FGFR1','LOXL2','PTK2','ROCK1','CXCR4','TGFB1','AGTR1','SHH')

sce.fibro$cluster2 <- sce.fibro$cluster
sce.fibro$cluster2[sce.fibro$cluster2%in%c('ACTG2_fibro','CD36_fibro','SORBS2_fibro')] <- 'mCAF'

sce.fibro$cluster2[sce.fibro$cluster2%in%c('CFD_fibro')] <- 'iCAF'

sce.fibro$cluster2[sce.fibro$cluster2%in%c('CFD_fibro')] <- 'iCAF'
#Neurotropic
sce.fibro$cluster2[sce.fibro$cluster2%in%c('GFRA3_fibro', 'KIT_fibro','PRKG1_fibro')] <- 'nCAF'
sce.fibro$cluster2[sce.fibro$cluster2%in%c('POSTN_fibro')] <- 'POSTN_CAF'

table(sce.fibro$cluster2)
sce.fibro@active.ident <- as.factor(sce.fibro$cluster2)
p <- DotPlot(sce.fibro,features = target.id)
data <- p$data
data$avg.exp.scaled
library(viridis)
color2 <- brewer.pal(9, "Blues")
data$features.plot <- factor(data$features.plot,levels = c(
  'FGFR1','DHRS3','AGTR1','PTK2','ROCK1','CXCR4','HAS2','LOXL2','TGFB1','SHH'
))
p = ggplot(data,aes(x=features.plot,y=id)) + 
  geom_point(aes(size=pct.exp,color=avg.exp.scaled))+
  scale_color_gradientn(colours = color2)+
  scale_size_continuous(range = c(0, 10),breaks = c(0,20,40,60))+
  labs(
    size="Percent Expressed",
    fill = 'Average Expression'
  )+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))
p
pdf('./Fig/ibro_clinical_target_dotplot.pdf',width = 7,height = 3.5,onefile = FALSE)
print(p)
dev.off()
?scale_size_continuous
library(RColorBrewer)
display.brewer.all()


####POSTN_fibro
sce.POSTN_fibro <- subset(sce.fibro,subset = cluster == 'POSTN_fibro')
sce.POSTN_fibro$group2 <- sce.POSTN_fibro$group
sce.POSTN_fibro$group2 <- if_else(sce.POSTN_fibro$group2 == 'normal','normal','tumor')
table(sce.POSTN_fibro$group2)
sce.POSTN_fibro@active.ident <- as.factor(sce.POSTN_fibro$group)
DefaultAssay(sce.POSTN_fibro) <- 'SCT'
sce.POSTN_fibro <- PrepSCTFindMarkers(sce.POSTN_fibro)
DEG.POSTN_fibro <- FindAllMarkers(sce.POSTN_fibro,logfc.threshold = 0.4,only.pos = T,min.pct = 0.1)
top30.DEG <- DEG.POSTN_fibro %>% group_by(cluster) %>% top_n(30,wt = avg_log2FC)
DimPlot(sce.POSTN_fibro)
FeaturePlot(sce.POSTN_fibro,features = 'OGN')

FeaturePlot(sce.combined.sct,features = 'CD34')
###chemo, FBLN1,CCDC80,SLC35F1,OGN,HSPB6,
###pri,WNT5A,CHN1,MMP1,ISG15,CST1
features <- c('FBLN1','CCDC80','SLC35F1','FKBP5','IGF1','WNT5A','CHN1','MMP1','ISG15','CST1')
data <- sce.POSTN_fibro@assays$SCT@data[features,] %>% as.matrix() %>%  t() %>% as.data.frame()
data$group <- sce.POSTN_fibro$group

exp <- data %>% group_by(group) %>% summarize(across(.cols = everything(),.fns = mean))
exp <- exp[c(1,3),] %>% column_to_rownames('group')
exp <- t(exp) %>% as.data.frame()

display.brewer.all()
color3 <- brewer.pal(9, "YlOrRd")[-1]
bk = c(seq(-1,1,by=0.01))

p1 <- pheatmap::pheatmap(exp,cluster_rows = F,cluster_cols = F,scale = 'column',
                   border_color = "white",
                   cellwidth = 50, cellheight = 50,
                   gaps_row = c(seq(9)), gaps_col = c(1,2),
                   color = colorRampPalette(colors = color3)(length(bk)),
                   breaks = bk)
p1
pdf('./Fig/Fibro.DEG_pri.vs_chemo.pdf',width = 4,height = 8,onefile = FALSE)
print(p1)
dev.off()

##免疫组化验证POSTN_fibro成纤维细胞在癌组织中高度富集。

##TCGA数据库验证POSTN_fibro在癌组织中显著富集

##
SpatialFeaturePlot(sptial.T2,features = 'MS4A1')
rm(sptial.T2)


load('../Fig4_analysis/Neu.sce.Rda')
Neu.sce$seurat_clusters %>% table()
Neu.markers <- FindAllMarkers(Neu.sce,logfc.threshold = 0.25,only.pos = T,min.pct = 0.1)
top30.Neu <- Neu.markers %>% group_by(cluster) %>% top_n(30,wt = avg_log2FC)
DimPlot(Neu.sce)
BiocManager::install('nnls')
library(nnls)
POSTN_fibro.markers <- c('POSTN','CTHRC1','INHBA','MMP14','COL5A1','COL5A2','FAP','THBS2')
Neu.CXCL8.markers <- c('FCGR3B','CSF3R','CCL4L2','CCL4','CCL3','CD83','C15orf48','PPIF')

macro.markers <- c('SPP1','CD68','CD163')


###成纤维细胞内皮细胞的转换。
FeaturePlot(sce.combined.sct,features = c('MS4A1'))



fibro.markers <- c('DCN','COL3A1','ACTA2','RGS5','MYH11')
epi.markers <- c('KRT18', 'EPCAM', 'CD24','KRT19','KRT17')
endo.markers <- c('PECAM1', 'CLDN5', 'CD34','CDH5','VWF')
B.markers <-  c('CD79A', 'JCHAIN','CD79B','MZB1','MS4A1')
T.markers <- c('CD3E', 'CD3D','CD8A','CD2','TRBC2')
myeloid.markers <- c('LYZ','C1QB','APOE','CD14','VCAN')
mast.markers <- c('TPSB2','TPSAB1','CPA3','TPSD1','CTSG')
Neutrophil.markers <- c('CSF3R','FCGR3B','G0S2','CXCL8','ACSL1')

FeaturePlot(Neu.sce,features = 'PPIF')
sce.combined.sct$cluster
DimPlot(sce.combined.sct,group.by = 'cluster')
sce.combined.sct$group2 <- as.character(sce.combined.sct$seurat_clusters)
sce.combined.sct$cluster[sce.combined.sct$group2 == '3'] <- 'Neutrophils'
sce.combined.sct@active.ident <- as.factor(sce.combined.sct$cluster)
markers <- FindAllMarkers(sce.combined.sct,logfc.threshold = 0.25,min.pct = 0.1,only.pos = T)
top30.main <- markers %>% group_by(cluster) %>% top_n(30,wt = avg_log2FC)
SpatialFeaturePlot(sptial.T4,features = c('DCN','CSF3R','CXCL8','ACSL1'))
DimPlot(sce.combined.sct,group.by = 'cluster')
FeaturePlot(sce.combined.sct,features = 'ACSL1')

######CXCL8_Neu and POSTN_fibro colization

###T1
sptial.T1 <- AddModuleScore(sptial.T1,features = list(markers = fibro.markers),name = 'fibro.score')
sptial.T1 <- AddModuleScore(sptial.T1,features = list(markers = myeloid.markers),name = 'myeloid.score')
sptial.T1 <- AddModuleScore(sptial.T1,features = list(markers = epi.markers),name = 'epi.score')
sptial.T1 <- AddModuleScore(sptial.T1,features = list(markers = endo.markers),name = 'endo.score')
sptial.T1 <- AddModuleScore(sptial.T1,features = list(markers = B.markers),name = 'Bcell.score')
sptial.T1 <- AddModuleScore(sptial.T1,features = list(markers = T.markers),name = 'Tcell.score')
sptial.T1 <- AddModuleScore(sptial.T1,features = list(markers = mast.markers),name = 'mast.score')
sptial.T1 <- AddModuleScore(sptial.T1,features = list(markers = Neutrophil.markers),name = 'neutrophil.score')




data <- FetchData(sptial.T1,vars = c('fibro.score1','myeloid.score1','epi.score1','endo.score1',
                                     'Bcell.score1','Tcell.score1','mast.score1','neutrophil.score1'))

data <- data[data$fibro.score1>0.1,]
cocorrence <- c()
for (i in colnames(data)) {
  tmp <- data.frame(fibro.score = data$fibro.score1,score = data[,i])
  x <- apply(tmp,1,function(x)sum(x>0.1))
  perc <- length(x[x == 2])/nrow(data)
  cocorrence <- c(cocorrence,perc)
}
names(cocorrence) <- colnames(data)
x
t1.coc <- cocorrence
t1.coc
####T2
sptial.T2 <- AddModuleScore(sptial.T2,features = list(markers = fibro.markers),name = 'fibro.score')
sptial.T2 <- AddModuleScore(sptial.T2,features = list(markers = myeloid.markers),name = 'myeloid.score')
sptial.T2 <- AddModuleScore(sptial.T2,features = list(markers = epi.markers),name = 'epi.score')
sptial.T2 <- AddModuleScore(sptial.T2,features = list(markers = endo.markers),name = 'endo.score')
sptial.T2 <- AddModuleScore(sptial.T2,features = list(markers = B.markers),name = 'Bcell.score')
sptial.T2 <- AddModuleScore(sptial.T2,features = list(markers = T.markers),name = 'Tcell.score')
sptial.T2 <- AddModuleScore(sptial.T2,features = list(markers = mast.markers),name = 'mast.score')
sptial.T2 <- AddModuleScore(sptial.T2,features = list(markers = Neutrophil.markers),name = 'neutrophil.score')

data <- FetchData(sptial.T2,vars = c('fibro.score1','myeloid.score1','epi.score1','endo.score1',
                                     'Bcell.score1','Tcell.score1','mast.score1','neutrophil.score1'))
data <- data[data$fibro.score1>0.1,]
cocorrence <- c()
for (i in colnames(data)) {
  tmp <- data.frame(fibro.score = data$fibro.score1,score = data[,i])
  x <- apply(tmp,1,function(x)sum(x>0.1))
  perc <- length(x[x == 2])/nrow(data)
  cocorrence <- c(cocorrence,perc)
}
names(cocorrence) <- colnames(data)

cocorrence
t2.coc <- cocorrence
t2.coc
t1.coc
###T4
sptial.T4 <- AddModuleScore(sptial.T4,features = list(markers = fibro.markers),name = 'fibro.score')
sptial.T4 <- AddModuleScore(sptial.T4,features = list(markers = myeloid.markers),name = 'myeloid.score')
sptial.T4 <- AddModuleScore(sptial.T4,features = list(markers = epi.markers),name = 'epi.score')
sptial.T4 <- AddModuleScore(sptial.T4,features = list(markers = endo.markers),name = 'endo.score')
sptial.T4 <- AddModuleScore(sptial.T4,features = list(markers = B.markers),name = 'Bcell.score')
sptial.T4 <- AddModuleScore(sptial.T4,features = list(markers = T.markers),name = 'Tcell.score')
sptial.T4 <- AddModuleScore(sptial.T4,features = list(markers = mast.markers),name = 'mast.score')
sptial.T4 <- AddModuleScore(sptial.T4,features = list(markers = Neutrophil.markers),name = 'neutrophil.score')

data <- FetchData(sptial.T4,vars = c('fibro.score1','myeloid.score1','epi.score1','endo.score1',
                                     'Bcell.score1','Tcell.score1','mast.score1','neutrophil.score1'))
data <- data[data$fibro.score1>0.1,]
cocorrence <- c()
for (i in colnames(data)) {
  tmp <- data.frame(fibro.score = data$fibro.score1,score = data[,i])
  x <- apply(tmp,1,function(x)sum(x>0.1))
  perc <- length(x[x == 2])/nrow(data)
  cocorrence <- c(cocorrence,perc)
}
names(cocorrence) <- colnames(data)

cocorrence
t4.coc <- cocorrence


t4.coc


coc.data <- data.frame(t1 = t1.coc,t2 = t2.coc,t4 = t4.coc)
coc.data <- coc.data[-1,]
###corcorence heatmap
library(RColorBrewer)
color3 <- brewer.pal(9, "YlOrRd")[-1]
bk = c(seq(0,1,by=0.01))

p1 <- pheatmap::pheatmap(coc.data,cluster_rows = F,cluster_cols = F,scale = 'none',
                         border_color = "white",display_numbers = T,number_color = 'black',
                         cellwidth = 50, cellheight = 50,
                         gaps_row = c(seq(7)), gaps_col = c(1,2,3),
                         color = colorRampPalette(colors = color3)(length(bk)),
                         breaks = bk)
pdf('./Fig/Fibro.corcorence.sp.heatmap.pdf',width = 5,height = 8,onefile = FALSE)
print(p1)
dev.off()

###cocorence HE visu
library(viridis)
sptial.T4$epi.score <- (sptial.T4$fibro.score1 + sptial.T4$epi.score1)/2
sptial.T4$endo.score <- (sptial.T4$fibro.score1 + sptial.T4$endo.score1)/2
sptial.T4$Bcell.score <- (sptial.T4$fibro.score1 + sptial.T4$Bcell.score1)/2
sptial.T4$Tcell.score <- (sptial.T4$fibro.score1 + sptial.T4$Tcell.score1)/2
sptial.T4$myeloid.score <- (sptial.T4$fibro.score1 + sptial.T4$myeloid.score1)/2
sptial.T4$neutrophil.score <- (sptial.T4$fibro.score1 + sptial.T4$neutrophil.score1)/2
sptial.T4$mast.score <- (sptial.T4$fibro.score1 + sptial.T4$mast.score1)/2

sp.features <- c('epi.score','endo.score','Bcell.score','Tcell.score','myeloid.score','neutrophil.score','mast.score')

p <- SpatialFeaturePlot(sptial.T4[,sptial.T4$mast.score1>0.1 & sptial.T4$fibro.score1>0.1],features = 'mast.score',max.cutoff = 1,)+scale_fill_viridis()
p
pdf('./Fig/spatialdim_fibro.score/fibro_mast_cocorence.T4.pdf',width = 5,height = 5,onefile = FALSE)
print(p)
dev.off()
###成纤维细胞亚群空间位置展示。
fibro.submarkers
top30.fibromarkers <- fibro.submarkers %>% group_by(cluster) %>% top_n(30,wt = avg_log2FC)
iCAF.markers <- c('PLA2G2A','CFD','C7','ADH1B','APOD','SERPINE2')
mCAF.markers <- c('MYH11','ACTA2','RGS5','TAGLN','ACTG2','PPP1R14A')
POSTN.markers <- c('POSTN','CTHRC1','INHBA','COL5A1','TMEM158','FAP')
nCAF.markers <- c('GFRA3','NRXN1','PLP1','DPP10','CDH19')
####T1
sce.T1 <- subset(sce.fibro,subset = sample == 'T1')

DefaultAssay(sce.T1) <- 'RNA'
allen_reference <- sce.T1 %>% SCTransform(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

DefaultAssay(sptial.T1) <- "SCT"
features <- c(iCAF.markers,mCAF.markers,POSTN.markers,nCAF.markers)
features <- features[features%in%rownames(allen_reference)]

anchors <- FindTransferAnchors(reference = allen_reference, query = sptial.T1, normalization.method = "SCT",
                               features = features,dims = 1:19)
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$cluster1, prediction.assay = TRUE,
                                  weight.reduction = sptial.T1[["pca"]], dims = 1:19)

sptial.T1[["predictions"]] <- predictions.assay
DefaultAssay(sptial.T1) <- "predictions"
table(sce.T1$cluster1)
library(viridis)
p <- SpatialFeaturePlot(sptial.T1, features = c("POSTN-fibro", "mCAF",'iCAF','nCAF'), pt.size.factor = 1.6, ncol = 4, crop = TRUE,
                   alpha = c(0.1,1))
p
pdf('./Fig/spatial_fibro_subcluster.T1.pdf',width =8,height = 4,onefile = FALSE )
print(p)
dev.off()
###T2
sce.T2 <- subset(sce.fibro,subset = sample == 'T2')

DefaultAssay(sce.T2) <- 'RNA'
allen_reference <- sce.T2 %>% SCTransform(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

DefaultAssay(sptial.T2) <- "SCT"
features <- c(iCAF.markers,mCAF.markers,POSTN.markers,nCAF.markers)
features <- features[features%in%rownames(allen_reference)]

anchors <- FindTransferAnchors(reference = allen_reference, query = sptial.T2, normalization.method = "SCT",
                               features = features,dims = 1:19)
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$cluster1, prediction.assay = TRUE,
                                  weight.reduction = sptial.T2[["pca"]], dims = 1:19)

sptial.T2[["predictions"]] <- predictions.assay
DefaultAssay(sptial.T2) <- "predictions"

library(viridis)
p2 <- SpatialFeaturePlot(sptial.T2, features = c("POSTN-fibro", "mCAF",'iCAF','nCAF'), pt.size.factor = 1.6, ncol = 4, crop = TRUE,
                         alpha = c(0.1,1))
p2
pdf('./Fig/spatial_fibro_subcluster.T2.pdf',width =8,height = 4,onefile = FALSE )
print(p)
dev.off()

###T4
sce.T4 <- subset(sce.fibro,subset = sample %in% c('N4','T4'))

DefaultAssay(sce.T4) <- 'RNA'
test <- FindAllMarkers(sce.T4,logfc.threshold = 0.25,only.pos = T,min.pct = 0.1)
test <- test %>% group_by(cluster) %>% top_n(30,wt = avg_log2FC)
allen_reference <- sce.T4 %>% SCTransform(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
FeaturePlot(sce.T4,features = 'KRT13')
sce.fibro$cluster1
sce.nCAF <- subset(sce.fibro,subset = cluster1 == 'nCAF')
FeaturePlot(sce.nCAF,features = 'NRXN1')
load('./sce.fibro.Rda')
















DefaultAssay(sptial.T4) <- "SCT"
features <- c(iCAF.markers,mCAF.markers,POSTN.markers,nCAF.markers)
features <- features[features%in%rownames(allen_reference)]
table(allen_reference$cluster1)
anchors <- FindTransferAnchors(reference = allen_reference, query = sptial.T4, normalization.method = "SCT",
                               features = features,dims = 1:19)
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$cluster1, prediction.assay = TRUE,
                                  weight.reduction = sptial.T4[["pca"]], dims = 1:19)

sptial.T4[["predictions"]] <- predictions.assay
DefaultAssay(sptial.T4) <- 'predictions'

library(viridis)
p4 <- SpatialFeaturePlot(sptial.T4, features = c("POSTN-fibro", "mCAF",'iCAF','nCAF'), pt.size.factor = 1.6, ncol = 4, crop = TRUE)
p4
pdf('./Fig/spatial_fibro_subcluster.T4.pdf',width =8,height = 4,onefile = FALSE )
print(p4)
dev.off()
DefaultAssay(sptial.T4) <- 'SCT'
SpatialFeaturePlot(sptial.T4,features = nCAF.markers)
sptial.T4@assays$Spatial


###成纤维细胞互作网络，空间定位。
###选取其中一个重点阐明。
###cellphoneDB


sce.combined.sct$cluster2 <- as.character(sce.combined.sct$seurat_clusters)
sce.combined.sct$cluster3 <- sce.combined.sct$cluster
sce.combined.sct$cluster3[sce.combined.sct$cluster2%in%c('22','10')] <- 'mCAF'
sce.combined.sct$cluster3[sce.combined.sct$cluster2%in%c('21')] <- 'nCAF'
sce.combined.sct$cluster3[sce.combined.sct$cluster2%in%c('11')] <- 'POSTN_CAF'
sce.combined.sct$cluster3[sce.combined.sct$cluster2%in%c('16')] <- 'iCAF'
DimPlot(sce.combined.sct,group.by = 'cluster3',label = T)
DefaultAssay(sce.combined.sct)
sce.combined.sct@active.ident <- as.factor(sce.combined.sct$cluster3)
BiocManager::install('SeuratDisk')
remotes::install_github("mojaveazure/seurat-disk")
####Seurat 转H5ad
library(SeuratDisk)

counts <- subset(sce.combined.sct,subset = sample!='N1')
counts <- CreateSeuratObject(counts@assays$RNA@data,assay = "RNA")
SaveH5Seurat(counts, filename = './nor_counts.RNA.h5Seurat',overwrite = T)
Convert('./nor_counts.RNA.h5Seurat', dest = "h5ad",overwrite = T)

metadata <- FetchData(sce.combined.sct,vars = c('cluster3','sample')) %>%
  dplyr::filter(sample!='N1') %>% dplyr::select(cluster3) %>% rownames_to_column('Cell')

colnames(metadata) <- c('Cell','cell_type')
write.table(metadata,file = 'metadata.txt',sep = '\t',row.names = F,quote = F)

####LR analysis
LR.pvalue <- read.table("D:\\jupyter-notebook\\ESCC\\cpdb\\results\\statistical_analysis_pvalues_method2.txt",
                        sep = '\t')
colnames(LR.pvalue) <- LR.pvalue[1,]
LR.pvalue <- LR.pvalue[-1,]
LR.pvalue <- cbind(LR.pvalue[,1:11],LR.pvalue[,grep('CAF',colnames(LR.pvalue))])

decov <- read.table("D:\\jupyter-notebook\\ESCC\\cpdb\\results\\statistical_analysis_deconvoluted_method2.txt",sep = '\t',header = T)
#####step1 填补LR_sig_long和RL_sig_long中基因名
combine_ID1 <- decov %>% select('gene_name','uniprot')
colnames(combine_ID1) <- c('genes','id')
combine_ID2 <- decov %>% select('gene_name','complex_name')
colnames(combine_ID2) <- c('genes','id')
combine_ID <- rbind(combine_ID1,combine_ID2)


partner_a_id <- LR.pvalue %>% select('partner_a') %>% pull(partner_a) %>% str_split(':',2) %>% 
  lapply(function(x){x[2]}) %>% unlist()
LR.pvalue$gene_a <-combine_ID %>% slice(match(partner_a_id,combine_ID$id)) %>% pull(genes)

partner_b_id <- LR.pvalue %>% select('partner_b') %>% pull(partner_b) %>% str_split(':',2) %>% 
  lapply(function(x){x[2]}) %>% unlist()
LR.pvalue$gene_b <-combine_ID %>% slice(match(partner_b_id,combine_ID$id)) %>% pull(genes)


####整合配体受体顺序
#condition 1, receptor a is true and receptor b is false
RL_pdata <- LR.pvalue %>% filter(receptor_a == 'True' & receptor_b == 'False') %>% 
  select(!matches('CAF.*CAF'))
RL_sig_long <- RL_pdata %>% pivot_longer(12:length(colnames(RL_pdata)), 
                                         names_to = "group", 
                                         values_to = "p_value") %>% 
  filter(p_value<0.05) %>% separate(col = group, into = c('R_cells','L_cells'),sep = "\\|") %>% 
  select(gene_b,gene_a,L_cells,R_cells,p_value)
colnames(RL_sig_long) <- c('ligands','receptors','L_cells','R_cells','p_value')

#condition 2, receptor b is true and receptor a is false
LR_pdata <- LR.pvalue %>% filter(receptor_a == 'False' & receptor_b == 'True') %>% 
  select(!matches('CAF.*CAF'))
LR_sig_long <- LR_pdata %>% pivot_longer(12:length(colnames(LR_pdata)), 
                                         names_to = "group", 
                                         values_to = "p_value") %>% 
  filter(p_value<0.05) %>% separate(col = group, into = c('L_cells','R_cells'),sep = "\\|") %>% 
  select(gene_a,gene_b,L_cells,R_cells,p_value)
colnames(LR_sig_long) <- c('ligands','receptors','L_cells','R_cells','p_value')

interaction_sig_long <- rbind(RL_sig_long,LR_sig_long) %>% unite(group,L_cells,R_cells, sep = "|",remove = FALSE)
save(interaction_sig_long,file = 'fibro_RL_interaction.Rda')
fibroL <- interaction_sig_long[grep('CAF',interaction_sig_long$L_cells),]

sce.fibro$cluster1 <- sce.fibro$cluster
table(sce.fibro$cluster1)
sce.fibro$cluster1[sce.fibro$cluster1%in%c('ACTG2_fibro','CD36_fibro','SORBS2_fibro')] <- 'mCAF'
sce.fibro$cluster1[sce.fibro$cluster1%in%c('GFRA3_fibro','KIT_fibro','PRKG1_fibro')] <- 'nCAF'
sce.fibro$cluster1[sce.fibro$cluster1%in%c('CFD_fibro')] <- 'iCAF'

sce.fibro@active.ident <- as.factor(sce.fibro$cluster1)


fibro.submarkers <- FindAllMarkers(sce.fibro,logfc.threshold = 0.25,min.pct = 0.1)
sce.combined.sct@active.ident <- as.factor(sce.combined.sct$cluster)
markers <- FindAllMarkers(sce.combined.sct,logfc.threshold = 0.25,min.pct = 0.1)
DimPlot()
top300.submarkers <- fibro.submarkers %>% group_by(cluster) %>% top_n(300,wt = avg_log2FC)
top300.submarkers <- top300.submarkers[top300.submarkers$gene%in%c(fibroL$ligands),]
top300.submarkers <- top300.submarkers[top300.submarkers$avg_log2FC>0.5,]
exp <- sce.fibro@assays$RNA@data[unique(top300.submarkers$gene),] %>% as.matrix() %>% t() %>% as.data.frame()


exp$cluster <- sce.fibro$cluster1
exp <- exp %>% group_by(cluster) %>% summarize(across(.cols = everything(),.fns = mean))
exp <- exp[c(2,3,4,1),]

exp[,-1] <- scale(exp[,-1],center = F)
order.ligands <- colnames(exp)[-1]
exp <- exp %>% pivot_longer(-1,names_to = "genes",values_to = "exp")
exp$genes <- factor(exp$genes,levels = order.ligands)
top300.submarkers$log10p <- -log10(top300.submarkers$p_val_adj)
display.brewer.all()
color1 <- brewer.pal(11, "RdBu")[-6]
color1 <- color1[2:10]
###ligand mean expression
receptor.data <- fibroL[fibroL$ligands%in%top300.submarkers$gene,]

order.id <- receptor.data %>% distinct(ligands,.keep_all = T)




exp$genes %>% unique() %>% length()

test <- fibro.submarkers[fibro.submarkers$gene%in%order.id$ligands,]
test$FDR <- if_else(test$p_val_adj<0.001,3,if_else(test$p_val_adj<0.01,2,
                                                   if_else(test$p_val_adj<0.1,1,NA)))

p1 <- ggplot()+geom_tile(data = exp,aes(genes,cluster,fill = exp))+scale_fill_gradientn(colors = rev(color1))+
  coord_equal()+theme_minimal()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  geom_point(data = test,aes(gene,cluster,size=FDR),color = 'white')+
  scale_size_continuous(range = c(1,3),breaks = c(1,2,3),labels = c('<0.1','<0.01','<0.001'))+
  theme(legend.key = element_rect(fill = 'grey'))

p1
FeaturePlot(sce.fibro.T,features = 'CLU')
pdf('./Fig/fibro_ligand_exp.heatmap.pdf',width = 20,height = 5,onefile = FALSE)
print(p1)
dev.off()
###receptor expression dotplot
receptor.id <- unique(order.id$receptors)

sce.combined.sct$cluster[sce.combined.sct$cluster2 == '3'] <- 'Neutrophils'
celltype <- unique(sce.combined.sct$cluster)[-7]
perc.data <- data.frame()
for (i in celltype) {
  tmp <- subset(sce.combined.sct,subset = cluster == i)
  tmp.exp <- tmp@assays$RNA@data[receptor.id,] %>% as.matrix() %>% as.data.frame()
  cell.expressed.perc <- apply(tmp.exp,1,function(x){
    perc <- sum(x>0)/ncol(tmp.exp)
  })
  perc.data <- rbind(perc.data,cell.expressed.perc)
}
rownames(perc.data) <- celltype
colnames(perc.data) <- receptor.id
perc.data <- perc.data %>% rownames_to_column('celltype')


receptor.exp <- sce.combined.sct@assays$RNA@data[receptor.id,] %>% as.matrix() %>% t() %>% as.data.frame()

receptor.exp$celltype <- sce.combined.sct$cluster

receptor.mean <- receptor.exp %>% group_by(celltype) %>% summarize(across(.cols = everything(),.fns = mean)) %>% slice(-4)
receptor.mean <- receptor.mean[match(perc.data$celltype,receptor.mean$celltype),]
receptor.mean <- receptor.mean %>% pivot_longer(-1,names_to = 'receptors',values_to = 'mean.exp')
perc.data <- perc.data %>% pivot_longer(-1,names_to = 'receptors',values_to = 'pct')

receptor.mean$pct <- perc.data$pct

receptor.mean <- unite(receptor.mean,'mergeid',receptors,celltype,sep = '_',remove = F)
receptor.mean <- receptor.mean[receptor.mean$mergeid%in%order.id$mergeid,]

for (i in seq(48)) {
  if (order.id$mergeid[i]%in%receptor.mean$mergeid) {
    order.id[i,8] <- receptor.mean[which(receptor.mean$mergeid == order.id$mergeid[i]),4]
    order.id[i,9] <- receptor.mean[which(receptor.mean$mergeid == order.id$mergeid[i]),5]
  }else{
    next
  }
}
order.id$receptor.id <- seq(48)
order.id$receptors_exp <- as.numeric(order.id$receptors_exp)
order.id$receptors_percent <- as.numeric(order.id$receptors_percent)


p <- ggplot(order.id,aes(x=receptor.id,y=R_cells)) + 
  geom_point(aes(color= receptors_exp,size =receptors_percent))+
  scale_color_viridis(option = 'viridis')+
  labs(
    size="Percent Expressed",
    fill = 'Average Expression'
  )+
  scale_size_continuous(range = c(1,8),breaks = c(0.2,0.4,0.6))+
  scale_x_continuous(name = "receptors", breaks = seq(1, 48, by = 1), 
                     labels = order.id$receptors,
                     expand=c(0,1),
                     guide = guide_axis(position = "top"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 0))


p
pdf('./Fig/receptors_exp.dotplot.pdf',width = 20,height = 4,onefile = FALSE)
print(p)
dev.off()

####other cell ligand expression
fibroR <- interaction_sig_long[grep('CAF',interaction_sig_long$R_cells),]
sce.combined.sct@active.ident <- as.factor(sce.combined.sct$cluster)
DimPlot(sce.combined.sct)
markers <- FindAllMarkers(sce.combined.sct,logfc.threshold = 0.25,min.pct = 0.1,only.pos = T)
table(markers$cluster)
top300.submarkers <- markers %>% group_by(cluster) %>% top_n(300,wt = avg_log2FC)
top300.submarkers <- top300.submarkers[top300.submarkers$gene%in%c(fibroR$ligands),]
top300.submarkers <- top300.submarkers[top300.submarkers$avg_log2FC>0.5,] %>% dplyr::filter(cluster!='Fibroblasts')

exp <- sce.combined.sct@assays$RNA@data[unique(top300.submarkers$gene),] %>% as.matrix() %>% t() %>% as.data.frame()


sce.combined.sct$cluster
exp$cluster <- sce.combined.sct$cluster
exp <- exp %>% group_by(cluster) %>% summarize(across(.cols = everything(),.fns = mean)) %>% slice(-4)

exp[,-1] <- scale(exp[,-1],center = F)
order.ligands <- colnames(exp)[-1]
exp <- exp %>% pivot_longer(-1,names_to = "genes",values_to = "exp")
exp$genes <- factor(exp$genes,levels = order.ligands)
top300.submarkers$log10p <- -log10(top300.submarkers$p_val_adj)



ligand.data <- fibroR[fibroR$ligands%in%top300.submarkers$gene,]

order.id <- ligand.data %>% distinct(ligands,.keep_all = T)

order.id <- order.id[match(order.ligands,order.id$ligands),]


exp$genes %>% unique() %>% length()

test <- top300.submarkers[top300.submarkers$gene%in%order.id$ligands,]
test$FDR <- if_else(test$p_val_adj<0.001,3,if_else(test$p_val_adj<0.01,2,
                                                   if_else(test$p_val_adj<0.1,1,NA)))

p1 <- ggplot()+geom_tile(data = exp,aes(genes,cluster,fill = exp))+scale_fill_gradientn(colors = rev(color1))+
  coord_equal()+theme_minimal()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  geom_point(data = test,aes(gene,cluster,size=FDR),color = 'white')+
  
  theme(legend.key = element_rect(fill = 'grey'))
p1
pdf('./Fig/othercell_ligands.meanexp.pdf',width = 20,height = 5,onefile = FALSE )

print(p1)
dev.off()

###fibro receptors expression
receptor.id <- unique(order.id$receptors)


celltype <- unique(sce.fibro$cluster1)
celltype
perc.data <- data.frame()
for (i in celltype) {
  tmp <- subset(sce.fibro,subset = cluster1 == i)
  tmp.exp <- tmp@assays$RNA@data[receptor.id,] %>% as.matrix() %>% as.data.frame()
  cell.expressed.perc <- apply(tmp.exp,1,function(x){
    perc <- sum(x>0)/ncol(tmp.exp)
  })
  perc.data <- rbind(perc.data,cell.expressed.perc)
}
rownames(perc.data) <- celltype
colnames(perc.data) <- receptor.id
perc.data <- perc.data %>% rownames_to_column('celltype')


receptor.exp <- sce.fibro@assays$RNA@data[receptor.id,] %>% as.matrix() %>% t() %>% as.data.frame()

receptor.exp$celltype <- sce.fibro$cluster1

receptor.mean <- receptor.exp %>% group_by(celltype) %>% summarize(across(.cols = everything(),.fns = mean))
receptor.mean <- receptor.mean[match(perc.data$celltype,receptor.mean$celltype),]
receptor.mean <- receptor.mean %>% pivot_longer(-1,names_to = 'receptors',values_to = 'mean.exp')
perc.data <- perc.data %>% pivot_longer(-1,names_to = 'receptors',values_to = 'pct')

receptor.mean$pct <- perc.data$pct
receptor.mean$celltype <- gsub('POSTN_fibro','POSTN_CAF',x = receptor.mean$celltype)

receptor.mean <- unite(receptor.mean,'mergeid',receptors,celltype,sep = '_',remove = F)
order.id <- unite(order.id,'mergeid',receptors,R_cells,sep = '_',remove = F)
receptor.mean <- receptor.mean[receptor.mean$mergeid%in%order.id$mergeid,]
table(receptor.mean$mergeid%in%order.id$mergeid)

order.id$receptors_exp <- rep(NA,42)
order.id$receptors_percent <- rep(NA,42)
for (i in seq(42)) {
  if (order.id$mergeid[i]%in%receptor.mean$mergeid) {
    order.id[i,8] <- receptor.mean[which(receptor.mean$mergeid == order.id$mergeid[i]),4]
    order.id[i,9] <- receptor.mean[which(receptor.mean$mergeid == order.id$mergeid[i]),5]
  }else{
    next
  }
}
order.id$receptor.id <- seq(42)
order.id$receptors_exp <- as.numeric(order.id$receptors_exp)
order.id$receptors_percent <- as.numeric(order.id$receptors_percent)


p <- ggplot(order.id,aes(x=receptor.id,y=R_cells)) + 
  geom_point(aes(color= receptors_exp,size =receptors_percent))+
  scale_color_viridis(option = 'viridis')+
  labs(
    size="Percent Expressed",
    fill = 'Average Expression'
  )+
  scale_size_continuous(range = c(1,8),breaks = c(0.2,0.4,0.6))+
  scale_x_continuous(name = "receptors", breaks = seq(1, 42, by = 1), 
                     labels = order.id$receptors,
                     expand=c(0,1),
                     guide = guide_axis(position = "top"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 0))


p
pdf('./Fig/fibro_receptors_exp.dotplot.pdf',width = 20,height = 4,onefile = FALSE)
print(p)
dev.off()




####巨噬细胞论述与成纤维细胞之间的关系
####思路1，共定位，空间临近，受体配体互作，阐述关键通路，促癌机制。
####空间转录特征入手，筛选关键因子。
####中性粒细胞简单阐述空间分布，重点结合成纤维细胞分析。
####免疫球蛋白在微环境中大量产生于非B细胞组织。重点有意思现象。
####成纤维细胞内皮细胞转换，有一类群成纤维细胞和内皮细胞同时表达AQP1基因

rm(list = ls())
############role
DimPlot(sce.fibro,group.by = 'cluster')
table(sce.fibro$cluster)
testdata <- FetchData(sce.fibro,vars = c('cluster','group','sample')) %>%
  dplyr::filter(sample%in%c('N1','N4','N5','T1','T4','T5'))
testdata$group <- if_else(testdata$sample%in%c('N1','N4','N5'),'W','T')
data <- testdata
data$cluster <- if_else(data$cluster%in%c('CFD_fibro'),'CFD_fibro','zther')

mydata <- xtabs(~cluster+group,data = data)
mydata
fisher.test(mydata)


#ACTG2_fibro   CD36_fibro    CFD_fibro  GFRA3_fibro    KIT_fibro  POSTN_fibro  PRKG1_fibro 
#SORBS2_fibro 
OR.POSTN_fibro <- c(3.100037,2.2e-16)
OR.ACTG2_fibro <- c(0.5691036,4.298e-10)
OR.CD36_fibro <- c(0.9465307,0.5241)
OR.CFD_fibro <- c(0.7237159,4.021e-07)
OR.GFRA3_fibro <- c(1.397906,0.0005874)
OR.KIT_fibro <- c(0.1667971,2.2e-16)
OR.PRKG1_fibro <- c(0.5620057,0.02376)
OR.SORBS2_fibro <- c(0.7960821,0.001683)




nOR.POSTN_fibro <- c(0.3225768,2.2e-16)
nOR.ACTG2_fibro <- c(1.757149,4.298e-10)
nOR.CD36_fibro <- c(1.05649,0.5241)
nOR.CFD_fibro <- c(1.381758,4.021e-07)
nOR.GFRA3_fibro <- c(0.7153555,0.0005874)
nOR.KIT_fibro <- c(5.995307,2.2e-16)
nOR.PRKG1_fibro <- c(1.779341,0.02376)
nOR.SORBS2_fibro <- c(1.256152,0.001683)


data <- rbind(OR.POSTN_fibro,OR.ACTG2_fibro,OR.CD36_fibro,
              OR.CFD_fibro,OR.GFRA3_fibro,OR.KIT_fibro,
              OR.PRKG1_fibro,OR.SORBS2_fibro)

data2 <- rbind(nOR.POSTN_fibro,nOR.ACTG2_fibro,nOR.CD36_fibro,
               nOR.CFD_fibro,nOR.GFRA3_fibro,nOR.KIT_fibro,
               nOR.PRKG1_fibro,nOR.SORBS2_fibro)
id <- rownames(data)
id <- c(rep(id,2))
id
data <- rbind(data,data2) %>% as.data.frame()
colnames(data) <- c('Role','p')
data$id <- id
data$group <- c(rep('tumor',8),rep('adjacent',8))
data$value <- if_else(data$Role<1,'1',
                     if_else(data$Role<1.5,'2',
                             if_else(data$Role<3,'3','4')))
data$stra <- if_else(data$Role<1,'±',
                     if_else(data$Role<1.5,'+',
                             if_else(data$Role<3,'++','+++')))
data$id <- factor(data$id,levels = rev(c('OR.POSTN_fibro','OR.GFRA3_fibro','OR.KIT_fibro',
                                     'OR.ACTG2_fibro','OR.PRKG1_fibro',
                                     'OR.SORBS2_fibro',
                                     'OR.CFD_fibro','OR.CD36_fibro')))
library(RColorBrewer)
color1 <- brewer.pal(9, "YlOrBr")[c(3,5,6,7)]
display.brewer.all()
library(ggplot2)
cols <- color1
p <- ggplot(data, aes(id,group)) +
  geom_tile(aes(fill =value)) +
  geom_text(aes(label=stra), color="black", size=4.5) + # 把星号添加进去
  scale_fill_manual(values = cols)+
  labs(x=NULL,y=NULL) + # 去掉横纵坐标标题
  theme(axis.text.x = element_text(size=8,angle = 30,hjust = 1,color = "black"),
        axis.text.y = element_text(size=8,color = "black"),
        axis.ticks.y = element_blank(),
        panel.background=element_blank()) # 做一些简单的美化
p

pdf('./Fig/fibro_role.pdf',width = 6,height = 2.5,onefile = FALSE)
print(p)
dev.off()










