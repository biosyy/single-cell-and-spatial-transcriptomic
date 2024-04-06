rm(list = ls())
library(Seurat)
library(tidyverse)
library(cols4all)
library(ggthemes)
c4a_gui()
load('../Fig2_analysis/sce.epi.filter.Rda')
load('../spatial_integration/Spatial.data.Rda')
SpatialDimPlot(sptial.T3,label = T)
p <- SpatialFeaturePlot(sptial.T1,features = 'CD8A')
x <- p$data

markres <- FindAllMarkers(sce.epi.filter,logfc.threshold = 0.25,
                          min.pct = 0.1,only.pos = T)

top20.epi.markers <- markres %>% group_by(cluster) %>%
  top_n(30,wt = avg_log2FC)



marker.sp <- FindAllMarkers(se.t1,logfc.threshold = 0.25,min.pct = 0.1,only.pos = T)
top20.sp.markers <- marker.sp %>% group_by(cluster) %>%
  top_n(30,wt = avg_log2FC)


epi.cl1.signatures <- c('SPRR2A','LCN2','KRT13','MUC21','SPRR2D')

epi.cl2.signatures <- c('KRT15','PTHLH','COL17A1','SULF2','SLC7A8','EPCAM')
epi.cl2.signatures <- c('MT2A','PTHLH','COL17A1')
DefaultAssay(sptial.T1) <- 'Spatial'
sptial.T1 <- NormalizeData(sptial.T1)
sptial.T1 <- ScaleData(sptial.T1)
sptial.T1 <- AddModuleScore(sptial.T1,features = list(epi.cl1.signatures),name = 'epi.cl1.score')
sptial.T5 <- AddModuleScore(sptial.T5,features = list(epi.cl2.signatures),name = 'epi.cl2.score')
SpatialFeaturePlot(sptial.T1,features = 'MUC21')
p2 <-  VlnPlot(sptial.T1,features = 'epi.cl1.score1',pt.size = 0,group.by = 'seurat_clusters',sort = F)+
  ylab("cl1.score")+
  theme(legend.position = 'none',axis.title.x = element_blank(),
        plot.title = element_blank())
p2
SpatialDimPlot(sptial.T5,label = T)

pdf('./Fig/T4.epi.cl1.score,vlnplot.pdf',width = 7,height = 5,onefile = FALSE)
print(p2)
dev.off()
p3 <-  VlnPlot(sptial.T1,features = 'epi.cl2.score1',pt.size = 0,group.by = 'seurat_clusters',sort = F)+
  ylab("cl2.score")+
  theme(legend.position = 'none',axis.title.x = element_blank(),
        plot.title = element_blank())
p3
SpatialDimPlot(sptial.T3,label = T,pt.size.factor = 0)
SpatialFeaturePlot(sptial.T3,features = 'EPCAM',pt.size.factor = 2)
pdf('./Fig/T4.epi.cl2.score,vlnplot.pdf',width = 7,height = 5,onefile = FALSE)
print(p3)
dev.off()
#########spatialplot Seuratcluster
p <- SpatialDimPlot(sptial.T5, group.by = "seurat_clusters",alpha = 1)+
  theme(plot.title = element_blank())
p
c4a_gui()
mycol <- c4a('wright25',9)
data <- p$data
str(data)
p
p <- ggplot(data,aes(imagecol,imagerow,color = seurat_clusters))
p1 <- p+geom_point(shape = 16)+scale_y_reverse()+
  scale_color_manual(values = mycol)+
  theme_map()
p1
p2 <- p1+theme(legend.position = 'none')
pdf('./Fig/spatial.seurat_cluster.T5.pdf',width = 5,height = 5,onefile = FALSE)

print(p2)
dev.off()
library(ggpubr)
library(grid)
legend.p <- as_ggplot(get_legend(p1))

pdf('./Fig/T5.spatial.cluster.legend.pdf',width = 5,height = 5,onefile = FALSE)

print(legend.p)
dev.off()


p1 <- SpatialDimPlot(sptial.T2, group.by = "seurat_clusters")+
  theme(plot.title = element_blank())
p1
pdf('./Fig/T2.HE.cluster.pdf',width = 5,height = 5,onefile = FALSE)

print(p1)
dev.off()

mycol <- c4a('Set1',8)[2:3]

sptial.T1$group.cl1 <- if_else(sptial.T1$seurat_clusters%in%c(0),'epi.cl1','other')
sptial.T1$group.cl2 <- if_else(sptial.T1$seurat_clusters%in%c(1),'epi.cl2','other')

p <- SpatialDimPlot(sptial.T2, group.by = 'group.cl1')
p
data <- p$data
p <- ggplot(data,aes(imagecol,imagerow,color = group.cl1))
p1 <- p+geom_point(shape = 16)+scale_y_reverse()+
  scale_color_manual(values = rev(c('lightgray',mycol[1])))+
  theme_map()+theme(legend.position = 'none')
p1
pdf('./Fig/T2.epi.cl1.spatial.pdf',width = 5,height = 5,onefile = FALSE)
print(p1)
dev.off()

p <- SpatialDimPlot(sptial.T2, group.by = 'group.cl2')
data <- p$data
p <- ggplot(data,aes(imagecol,imagerow,color = group.cl2))
p1 <- p+geom_point(shape = 16)+scale_y_reverse()+
  scale_color_manual(values = rev(c('lightgray',mycol[2])))+
  theme_map()+theme(legend.position = 'none')
p1
pdf('./Fig/T1.epi.cl2.spatial.pdf',width = 5,height = 5,onefile = FALSE)
print(p1)
dev.off()
library(RColorBrewer)
display.brewer.all() 
color1 <- brewer.pal(11, "RdYlBu")
p <- SpatialFeaturePlot(sptial.T2,features = 'epi.cl1.score1')
data <- p$data
p <- ggplot(data,aes(imagecol,imagerow,color = epi.cl1.score1))
p1 <- p+geom_point(shape = 16)+scale_y_reverse()+
  scale_color_gradientn(colours = rev(color1))+
  theme_map()
p1
p1 <- p1+theme(legend.position = 'none')
pdf('./Fig/T2.epi.cl1.spatial.score.pdf',width = 5,height = 5,onefile = FALSE)
print(p1)
dev.off()






p <- SpatialFeaturePlot(sptial.T2,features = 'epi.cl1.score1')
p
data <- p$data
p <- ggplot(data,aes(imagecol,imagerow,color = epi.cl2.score1))
p1 <- p+geom_point(shape = 16)+scale_y_reverse()+
  scale_color_viridis_c()+
  theme_map()+theme(legend.position = 'none')
p1
pdf('./Fig/T2.epi.cl2.spatial.score.pdf',width = 5,height = 5,onefile = FALSE)
print(p1)
dev.off()


p <- SpatialFeaturePlot(sptial.T1,features = 'CD24')
data <- p$data
p <- ggplot(data,aes(imagecol,imagerow,color = CD24))
p1 <- p+geom_point(shape = 16)+scale_y_reverse()+
  scale_color_gradientn(colours = rev(color1))+
  theme_map()+theme(legend.position = 'none')
p1

p <- SpatialFeaturePlot(sptial.T2,features = 'CD24')
data <- p$data
p <- ggplot(data,aes(imagecol,imagerow,color = CD24))
p2 <- p+geom_point(shape = 16)+scale_y_reverse()+
  scale_color_gradientn(colours = rev(color1))+
  theme_map()+theme(legend.position = 'none')
p2
pdf('Fig/T1.CD24.spatialfeatureplot.pdf',width = 5,height = 5,onefile = FALSE)
print(p1)
dev.off()
pdf('Fig/T2.CD24.spatialfeatureplot.pdf',width = 5,height = 5,onefile = FALSE)
print(p2)
dev.off()


p <- SpatialFeaturePlot(sptial.T1,features = 'COL17A1')
p
data <- p$data
p <- ggplot(data,aes(imagecol,imagerow,color = COL17A1))
p1 <- p+geom_point(shape = 16)+scale_y_reverse()+
  scale_color_gradientn(colours = rev(color1))+
  theme_map()+theme(legend.position = 'right')
p1

pdf('Fig/T1.COL17A1.spatialfeatureplot.pdf',width = 5,height = 5,onefile = FALSE)
print(p1)
dev.off()




#####

markres <- FindAllMarkers(sptial.T3,logfc.threshold = 0.25,min.pct = 0.1,only.pos = T)
top20 <- markres %>% group_by(cluster) %>% top_n(30,wt = avg_log2FC)

SpatialFeaturePlot(sptial.T1,features = epi.cl2.signatures)

VlnPlot(sptial.T1,features = 'epi.cl2.score1',pt.size = 0)


FeaturePlot(sce.epi.filter,features = epi.cl2.signatures)



####
load('../share_data/sce.combined.filter.Rda')
sce.combined.filter@active.ident <- as.factor(sce.combined.filter$celltype)
markres <- FindAllMarkers(sce.combined.filter,logfc.threshold = 0.25,min.pct = 0.1,only.pos = T)
top20.markers <- markres %>% group_by(cluster) %>% top_n(20,wt = avg_log2FC)


genes <- list(Bcells = c('IGHG4','IGHG3','IGHG1','JCHAIN','MS4A1'),
              CD4.Tcells = c('FOXP3','CD4','ICOS'),
              CD8.Tcells = c('IFNG','NKG7','CD8A','CD8B'),
              cDC =c('CD1C','S100B','FCER1A','CLEC9A','LAMP3','CLEC10A'),
              Endothelial = c('PLVAP','GNG11','AQP1','VWF','CLDN5'),
              Fibroblasts = c('COL3A1','COL1A1','COL1A2','POSTN','TAGLN'),
              mast = c('TPSB2','TPSAB1','CPA3','TPSD1','CTSG'),
              mono.macro = c('LYZ','APOE','C1QA','C1QB','C1QC'),
              Neutrophil = c('CSF3R','FCGR3B','G0S2','CXCL8','ACSL1'),
              NK.NKT = c('XCL1','TRDC','XCL2','CTSW','PRF1'),
              epi.cl1 = c('SPRR2A','LCN2','KRT13','MUC21','SPRR2D'),
              epi.cl2 = c('KRT15','PTHLH','COL17A1','SULF2','SLC7A8')
)

FeaturePlot(sce.combined.filter,features = 'TAGLN')
sptial.T1 <- AddModuleScore(sptial.T1,features = list(epi.cl1.signatures),name = 'epi.cl1.score')

for (i in seq(12)){
  sptial.T1 <- AddModuleScore(sptial.T1,features = genes[i],name = names(genes)[i])
}

data <- sptial.T1@meta.data[,10:21]

cor.data <- cor(data)
library(pheatmap)
bk = c(seq(-0.4,0.4,by=0.01))
pheatmap::pheatmap(cor.data,
                   breaks = bk)

SpatialDimPlot(sptial.T1,label = T)

sp.T1.markers <- FindAllMarkers(sptial.T1,logfc.threshold = 0.25,min.pct = 0.1,only.pos = T)
top20.T1.markers <- sp.T1.markers %>% group_by(cluster) %>% top_n(20,wt = avg_log2FC)

sptial.T2 <- NormalizeData(sptial.T2)
sptial.T2 <- ScaleData(sptial.T2)
genes <- c('MYH11','MYL9','DES','SYNM','COL11A1','COL1A2','COL3A1','POSTN',
           'IGHG3','IGHG1','JCHAIN','CD74','MS4A1',
           'KRT19','EPCAM','KRT8','COL17A1','KRT15','KRT16','KRT13',
           'SPP1','LYZ','APOC1','CD68','C1QA','C1QB')
DefaultAssay(sptial.T4) <- 'Spatial'
sptial.T4 <- NormalizeData(sptial.T4)
sptial.T4 <- ScaleData(sptial.T4)
DoHeatmap(sptial.T1,features = genes)
p <- SpatialDimPlot(sptial.T1,label = T)
pdf('testsp1.pdf',width = 7,height = 5,onefile = FALSE)
print(p)
dev.off()

###spatial T1 region signature expression heatmap
exp <- sptial.T1@assays$Spatial@data[genes,] %>% as.matrix() %>% t() %>% as.data.frame()
exp$cluster <- sptial.T1$seurat_clusters
mean.exp <- exp %>% group_by(cluster) %>% summarise(across(.cols = everything(),.fns = mean)) %>% 
  column_to_rownames('cluster')
mean.exp <- mean.exp %>% t() %>% as.data.frame()
bk = c(seq(-2,2,by=0.01))
order.id <- c(0,10,11,2,6,12,4,1,5,9,3,7,8,13)
order.id <- order.id + 1
data <- mean.exp[,order.id]
pheatmap::pheatmap(data,cluster_rows = F,cluster_cols = F,scale = 'row',
                   border = F,
                   color = colorRampPalette(colors = c("#1BB7E5","#000000","#FFFF00"))(length(bk)))

col_anno <- data.frame(cluster = c(rep('CAF',6),rep('Bcells',2),rep('Epi',5),rep('mye',1)))
rownames(col_anno) <- colnames(data)
col_anno$cluster <- as.factor(col_anno$cluster)
library(RColorBrewer)
color1 <- brewer.pal(9, "Set1")[1:4]
table(col_anno$cluster)
ann_colors = list(
  cluster = c(Bcells=color1[1], CAF=color1[2],
              Epi=color1[3],mye=color1[4]))

color2 <- brewer.pal(11, "RdBu")
color2 <- color2[-6]
color2 <- color2[2:10]
pt1 <- pheatmap::pheatmap(data,cluster_rows = F,cluster_cols = F,scale = 'row',
                          border = F,
                          annotation_col = col_anno,
                          annotation_colors = ann_colors,
                          color = colorRampPalette(colors = rev(color2))(length(bk)),
                          breaks = bk)
pdf('./Fig/T1.region.signature.htp1.pdf',width = 5,height = 7,onefile = FALSE)
print(pt1)
dev.off()

###spatial T1 region cluster
sptial.T1$region.cluster <- sptial.T1$seurat_clusters %>% as.character()
order.id <- c(0,10,11,2,6,12,4,1,5,9,3,7,8,13)
col_anno <- data.frame(cluster = c(rep('CAF',6),rep('Bcells',2),rep('Epi',5),rep('mye',1)))
sptial.T1$region.cluster[sptial.T1$region.cluster%in%c('0','10','11','2','6','12')] <- 'CAF'
sptial.T1$region.cluster[sptial.T1$region.cluster%in%c('4','1')] <- 'Bcells'
sptial.T1$region.cluster[sptial.T1$region.cluster%in%c('5','9','3','7','8','14')] <- 'Epi'
sptial.T1$region.cluster[sptial.T1$region.cluster%in%c('13')] <- 'mye'

p <- SpatialDimPlot(sptial.T1, group.by = "region.cluster",alpha = 1,label = T)+
  theme(plot.title = element_blank())
p
color1 <- brewer.pal(9, "Set1")[1:4]
ann_colors = list(
  cluster = c(Bcells=color1[1], CAF=color1[2],
              Epi=color1[3],mye=color1[4]))
data <- p$data

p <- ggplot(data,aes(imagecol,imagerow,color = region.cluster))
p1 <- p+geom_point(shape = 16)+scale_y_reverse()+
  scale_color_manual(values = c(Bcells=color1[1], CAF=color1[2],
                                Epi=color1[3],mye=color1[4]))+
  theme_map()
p1
p2 <- p1+theme(legend.position = 'none')
pdf('./Fig/spatial.region.cluster.T1.pdf',width = 5,height = 5,onefile = FALSE)

print(p2)
dev.off()

legend.p <- as_ggplot(get_legend(p1))

pdf('./Fig/spatial.region.cluster.legend.pdf',width = 5,height = 5,onefile = FALSE)

print(legend.p)
dev.off()


#####spatial T2 region signature expression heatmap
exp <- sptial.T2@assays$Spatial@data[genes,] %>% as.matrix() %>% t() %>% as.data.frame()
exp$cluster <- sptial.T2$seurat_clusters
mean.exp <- exp %>% group_by(cluster) %>% summarise(across(.cols = everything(),.fns = mean)) %>% 
  column_to_rownames('cluster')
mean.exp <- mean.exp %>% t() %>% as.data.frame()
bk = c(seq(-2,2,by=0.01))
order.id <- c(6,3,2,4,8,1,7,0,5)
order.id <- order.id + 1
data.T2 <- mean.exp[,order.id]
pheatmap::pheatmap(data.T2,cluster_rows = F,cluster_cols = F,scale = 'row',
                   border = F,
                   color = colorRampPalette(colors = c("#1BB7E5","#000000","#FFFF00"))(length(bk)))

col_anno <- data.frame(cluster = c(rep('CAF',2),rep('Bcells',3),rep('Epi',3),rep('mye',1)))
rownames(col_anno) <- colnames(data.T2)
col_anno$cluster <- as.factor(col_anno$cluster)
library(RColorBrewer)
color1 <- brewer.pal(9, "Set1")[1:4]
table(col_anno$cluster)
ann_colors = list(
  cluster = c(Bcells=color1[1], CAF=color1[2],
              Epi=color1[3],mye=color1[4]))

color2 <- brewer.pal(11, "RdBu")
color2 <- color2[-6]
color2 <- color2[2:10]
pt2 <- pheatmap::pheatmap(data.T2,cluster_rows = F,cluster_cols = F,scale = 'row',
                          border = F,
                          annotation_col = col_anno,
                          annotation_colors = ann_colors,
                          color = colorRampPalette(colors = rev(color2))(length(bk)),
                          breaks = bk)
pdf('./Fig/T2.region.signature.htp.pdf',width = 5,height = 7,onefile = FALSE)
print(pt2)
dev.off()

###spatial T2 region cluster
sptial.T2$region.cluster <- sptial.T2$seurat_clusters %>% as.character()
SpatialDimPlot(sptial.T2,label = T)
order.id <- c(6,3,2,4,8,1,7,0,5)
col_anno <- data.frame(cluster = c(rep('CAF',3),rep('Bcells',2),rep('Epi',3),rep('mye',1)))
sptial.T2$region.cluster[sptial.T2$region.cluster%in%c('6','3')] <- 'CAF'
sptial.T2$region.cluster[sptial.T2$region.cluster%in%c('2','4','8')] <- 'Bcells'
sptial.T2$region.cluster[sptial.T2$region.cluster%in%c('1','7','0')] <- 'Epi'
sptial.T2$region.cluster[sptial.T2$region.cluster%in%c('5')] <- 'mye'

p <- SpatialDimPlot(sptial.T2, group.by = "region.cluster",alpha = 1,label = T)+
  theme(plot.title = element_blank())
p
color1 <- brewer.pal(9, "Set1")[1:4]
data <- p$data

p <- ggplot(data,aes(imagecol,imagerow,color = region.cluster))
p1 <- p+geom_point(shape = 16)+scale_y_reverse()+
  scale_color_manual(values = c(Bcells=color1[1], CAF=color1[2],
                                Epi=color1[3],mye=color1[4]))+
  theme_map()
p1
p2 <- p1+theme(legend.position = 'none')
pdf('./Fig/spatial.region.cluster.T2.pdf',width = 5,height = 5,onefile = FALSE)

print(p2)
dev.off()

legend.p <- as_ggplot(get_legend(p1))

pdf('./Fig/spatial.region.cluster.legend.pdf',width = 5,height = 5,onefile = FALSE)

print(legend.p)
dev.off()


#####spatial T3 region signature expression heatmap
exp <- sptial.T3@assays$Spatial@data[genes,] %>% as.matrix() %>% t() %>% as.data.frame()
exp$cluster <- sptial.T3$seurat_clusters
mean.exp <- exp %>% group_by(cluster) %>% summarise(across(.cols = everything(),.fns = mean)) %>% 
  column_to_rownames('cluster')
mean.exp <- mean.exp %>% t() %>% as.data.frame()
bk = c(seq(-2,2,by=0.01))

order.id <- c(0,1,2,8,4,6,3,5,7)
order.id <- order.id + 1
data.T3 <- mean.exp[,order.id]
pheatmap::pheatmap(data.T3,cluster_rows = F,cluster_cols = F,scale = 'row',
                   border = F,
                   color = colorRampPalette(colors = c("#1BB7E5","#000000","#FFFF00"))(length(bk)))

col_anno <- data.frame(cluster = c(rep('CAF',6),rep('Bcells',1),rep('Epi',1),rep('mye',1)))
rownames(col_anno) <- colnames(data.T3)
col_anno$cluster <- as.factor(col_anno$cluster)
library(RColorBrewer)
color1 <- brewer.pal(9, "Set1")[1:4]
table(col_anno$cluster)
ann_colors = list(
  cluster = c(Bcells=color1[1], CAF=color1[2],
              Epi=color1[3],mye=color1[4]))

color2 <- brewer.pal(11, "RdBu")
color2 <- color2[-6]
color2 <- color2[2:10]
pt3 <- pheatmap::pheatmap(data.T3,cluster_rows = F,cluster_cols = F,scale = 'row',
                          border = F,
                          annotation_col = col_anno,
                          annotation_colors = ann_colors,
                          color = colorRampPalette(colors = rev(color2))(length(bk)),
                          breaks = bk)

pdf('./Fig/T3.region.signature.htp.pdf',width = 5,height = 7,onefile = FALSE)
print(pt3)
dev.off()

###spatial T3 region cluster
markres <- FindAllMarkers(sptial.T3,logfc.threshold = 0.25,min.pct = 0.1,only.pos = T)
top20 <- markres %>% group_by(cluster) %>% top_n(30,wt = avg_log2FC)
SpatialDimPlot(sptial.T3,label = T)
SpatialFeaturePlot(sptial.T3,features = 'MYLK')
sptial.T2$region.cluster <- sptial.T2$seurat_clusters %>% as.character()
order.id <- c(0,1,2,8,4,6,3,5,7)
sptial.T3$region.cluster <- sptial.T3$seurat_clusters %>% as.character()
sptial.T3$region.cluster[sptial.T3$region.cluster%in%c('0','1','2','8','4','6')] <- 'CAF'
sptial.T3$region.cluster[sptial.T3$region.cluster%in%c('3')] <- 'Bcells'
sptial.T3$region.cluster[sptial.T3$region.cluster%in%c('5')] <- 'Epi'
sptial.T3$region.cluster[sptial.T3$region.cluster%in%c('7')] <- 'mye'

p <- SpatialDimPlot(sptial.T3, group.by = "region.cluster",alpha = 1,label = T)+
  theme(plot.title = element_blank())
p
p <- SpatialDimPlot(sptial.T2, group.by = "seurat_clusters",alpha = 1,label = T)
p
color1 <- brewer.pal(9, "Set1")[1:4]
data <- p$data

p <- ggplot(data,aes(imagecol,imagerow,color = region.cluster))
p1 <- p+geom_point(shape = 16)+scale_y_reverse()+
  scale_color_manual(values = c(Bcells=color1[1], CAF=color1[2],
                                Epi=color1[3],mye=color1[4]))+
  theme_map()
p1
p2 <- p1+theme(legend.position = 'none')
pdf('./Fig/spatial.region.cluster.T3.pdf',width = 5,height = 5,onefile = FALSE)

print(p2)
dev.off()

legend.p <- as_ggplot(get_legend(p1))

pdf('./Fig/spatial.region.cluster.legend.pdf',width = 5,height = 5,onefile = FALSE)

print(legend.p)
dev.off()

######spatial T4 region signature expression heatmap
exp <- sptial.T4@assays$Spatial@data[genes,] %>% as.matrix() %>% t() %>% as.data.frame()
exp$cluster <- sptial.T4$seurat_clusters
mean.exp <- exp %>% group_by(cluster) %>% summarise(across(.cols = everything(),.fns = mean)) %>% 
  column_to_rownames('cluster')
mean.exp <- mean.exp %>% t() %>% as.data.frame()

bk = c(seq(-2,2,by=0.01))
order.id <- c(0,1,4,5,3,2,6)
order.id <- order.id + 1
data.T4 <- mean.exp[,order.id]
pheatmap::pheatmap(data.T4,cluster_rows = F,cluster_cols = F,scale = 'row',
                   border = F,
                   color = colorRampPalette(colors = c("#1BB7E5","#000000","#FFFF00"))(length(bk)))

col_anno <- data.frame(cluster = c(rep('CAF',3),rep('Bcells',2),rep('Epi',1),rep('mye',1)))
rownames(col_anno) <- colnames(data.T4)
col_anno$cluster <- as.factor(col_anno$cluster)
library(RColorBrewer)
color1 <- brewer.pal(9, "Set1")[1:4]
table(col_anno$cluster)
ann_colors = list(
  cluster = c(Bcells=color1[1], CAF=color1[2],
              Epi=color1[3],mye=color1[4]))

color2 <- brewer.pal(11, "RdBu")
color2 <- color2[-6]
color2 <- color2[2:10]
pt4 <- pheatmap::pheatmap(data.T4,cluster_rows = F,cluster_cols = F,scale = 'row',
                          border = F,
                          annotation_col = col_anno,
                          annotation_colors = ann_colors,
                          color = colorRampPalette(colors = rev(color2))(length(bk)),
                          breaks = bk)
pdf('./Fig/T4.region.signature.htp.pdf',width = 5,height = 7,onefile = FALSE)
print(pt4)
dev.off()

###spatial T4 region cluster
sptial.T4$region.cluster <- sptial.T4$seurat_clusters %>% as.character()
order.id <- c(0,1,4,5,3,2,6)
sptial.T4$region.cluster <- sptial.T4$seurat_clusters %>% as.character()
sptial.T4$region.cluster[sptial.T4$region.cluster%in%c('0','1','4')] <- 'CAF'
sptial.T4$region.cluster[sptial.T4$region.cluster%in%c('5','3')] <- 'Bcells'
sptial.T4$region.cluster[sptial.T4$region.cluster%in%c('2')] <- 'Epi'
sptial.T4$region.cluster[sptial.T4$region.cluster%in%c('6')] <- 'mye'

p <- SpatialDimPlot(sptial.T4, group.by = "region.cluster",alpha = 1,label = T)+
  theme(plot.title = element_blank())
p
color1 <- brewer.pal(9, "Set1")[1:4]
data <- p$data

p <- ggplot(data,aes(imagecol,imagerow,color = region.cluster))
p
p1 <- p+geom_point(shape = 16)+scale_y_reverse()+
  scale_color_manual(values = c(Bcells=color1[1], CAF=color1[2],
                                Epi=color1[3],mye=color1[4]))+
  theme_map()
p1
p2 <- p1+theme(legend.position = 'none')
p2
pdf('./Fig/spatial.region.cluster.T4.pdf',width = 5,height = 5,onefile = FALSE)

print(p2)
dev.off()


######spatial T5 region signature expression heatmap
DefaultAssay(sptial.T5) <- 'Spatial'
markres <- FindAllMarkers(sptial.T5,logfc.threshold = 0.25,min.pct = 0.1,only.pos = T)
top20 <- markres %>% dplyr::group_by(cluster) %>% top_n(20,wt = avg_log2FC)
order.id <- c(0,2,3,5,1,6,4,7)
exp <- sptial.T5@assays$Spatial@data[genes,] %>% as.matrix() %>% t() %>% as.data.frame()
exp$cluster <- sptial.T5$seurat_clusters
mean.exp <- exp %>% group_by(cluster) %>% summarise(across(.cols = everything(),.fns = mean)) %>% 
  column_to_rownames('cluster')
mean.exp <- mean.exp %>% t() %>% as.data.frame()

bk = c(seq(-2,2,by=0.01))
order.id <- c(0,3,5,2,4,1,6,7)
order.id <- order.id + 1
data.T5 <- mean.exp[,order.id]
pheatmap::pheatmap(data.T5,cluster_rows = F,cluster_cols = F,scale = 'row',
                   border = F,
                   color = colorRampPalette(colors = c("#1BB7E5","#000000","#FFFF00"))(length(bk)))

col_anno <- data.frame(cluster = c(rep('CAF',3),rep('Bcells',2),rep('Epi',2),rep('mye',1)))
rownames(col_anno) <- colnames(data.T5)
col_anno$cluster <- as.factor(col_anno$cluster)
library(RColorBrewer)
color1 <- brewer.pal(9, "Set1")[1:4]
table(col_anno$cluster)
ann_colors = list(
  cluster = c(Bcells=color1[1], CAF=color1[2],
              Epi=color1[3],mye=color1[4]))

color2 <- brewer.pal(11, "RdBu")
color2 <- color2[-6]
color2 <- color2[2:10]
pt5 <- pheatmap::pheatmap(data.T5,cluster_rows = F,cluster_cols = F,scale = 'row',
                          border = F,
                          annotation_col = col_anno,
                          annotation_colors = ann_colors,
                          color = colorRampPalette(colors = rev(color2))(length(bk)),
                          breaks = bk)
pdf('./Fig/T5.region.signature.htp.pdf',width = 5,height = 7,onefile = FALSE)
print(pt5)
dev.off()

###spatial T5 region cluster
sptial.T5$region.cluster <- sptial.T5$seurat_clusters %>% as.character()
order.id <- c(0,3,5,2,4,1,6,7)
sptial.T5$region.cluster <- sptial.T5$seurat_clusters %>% as.character()
sptial.T5$region.cluster[sptial.T5$region.cluster%in%c('0','3','5')] <- 'CAF'
sptial.T5$region.cluster[sptial.T5$region.cluster%in%c('2','4')] <- 'Bcells'
sptial.T5$region.cluster[sptial.T5$region.cluster%in%c('1','6')] <- 'Epi'
sptial.T5$region.cluster[sptial.T5$region.cluster%in%c('7')] <- 'mye'

p <- SpatialDimPlot(sptial.T5, group.by = "region.cluster",alpha = 1,label = T)+
  theme(plot.title = element_blank())

p <- SpatialDimPlot(sptial.T5, group.by = "seurat_clusters",alpha = 1,label = T)
p
color1 <- brewer.pal(9, "Set1")[1:4]
data <- p$data

p <- ggplot(data,aes(imagecol,imagerow,color = region.cluster))
p1 <- p+geom_point(shape = 16)+scale_y_reverse()+
  scale_color_manual(values = c(Bcells=color1[1], CAF=color1[2],
                                Epi=color1[3],mye=color1[4]))+
  theme_map()
p1
p2 <- p1+theme(legend.position = 'none')
pdf('./Fig/spatial.region.cluster.T5.pdf',width = 5,height = 5,onefile = FALSE)

print(p2)
dev.off()



####heatmap epi score across all spatial region
DefaultAssay(sptial.T1)
sptial.T1 <- NormalizeData(sptial.T1)
sptial.T1 <- ScaleData(sptial.T1)

sptial.T1$group <- as.character(sptial.T1$seurat_clusters)
sptial.T1$group[sptial.T1$group%in%c('2','6','12','11')] <- 'CAF'
sptial.T1$group[sptial.T1$group%in%c('1','4')] <- 'Bcells'
sptial.T1$group[sptial.T1$group%in%c('3','5','7','8','9')] <- 'Epi'
sptial.T1$group[sptial.T1$group%in%c('13')] <- 'mye'
sptial.T1$epi.cl1.score1
data <- FetchData(sptial.T1,vars = c('epi.cl1.score1','epi.cl2.score1','group'))
data <- data[data$group != '14',]

mean.dat <- data %>% group_by(group) %>% summarise(across(.cols = everything(),.fns = mean)) %>% dplyr::slice(c(-1,-2,-5))

mean.dat[,c(2,3)] <- mean.dat[,c(2,3)]+0.07
data.long.T1 <- mean.dat %>% pivot_longer(cols = c(-1),names_to = 'stra',values_to = 'score')
data.long.T1$lable <- if_else(data.long.T1$score>0.01,'+','-')

sptial.T2$group <- as.character(sptial.T2$seurat_clusters)
sptial.T2$group[sptial.T2$group%in%c('6','3')] <- 'CAF'
sptial.T2$group[sptial.T2$group%in%c('8','4','2')] <- 'Bcells'
sptial.T2$group[sptial.T2$group%in%c('0','1','7')] <- 'Epi'
sptial.T2$group[sptial.T2$group%in%c('5')] <- 'mye'
sptial.T2$epi.cl2.score1
data <- FetchData(sptial.T2,vars = c('epi.cl1.score1','epi.cl2.score1','group'))


mean.dat.T2 <- data %>% group_by(group) %>% summarise(across(.cols = everything(),.fns = mean))
mean.dat.T2 <- mean.dat.T2 %>% slice(c(-3))
data.long.T2 <- mean.dat.T2 %>% pivot_longer(cols = c(-1),names_to = 'stra',values_to = 'score')
data.long.T2$lable <- if_else(data.long.T2$score>0,'+','-')


sptial.T3$group <- as.character(sptial.T3$seurat_clusters)
sptial.T3$group[sptial.T3$group%in%c('0','1','2','4','6','8')] <- 'CAF'
sptial.T3$group[sptial.T3$group%in%c('3')] <- 'Bcells'
sptial.T3$group[sptial.T3$group%in%c('5')] <- 'Epi'
sptial.T3$group[sptial.T3$group%in%c('7')] <- 'mye'
data <- FetchData(sptial.T3,vars = c('epi.cl1.score1','epi.cl2.score1','group'))


mean.dat.T3 <- data %>% group_by(group) %>% summarise(across(.cols = everything(),.fns = mean))
mean.dat.T3 <- mean.dat.T3 %>% slice(c(-3))
data.long.T3 <- mean.dat.T3 %>% pivot_longer(cols = c(-1),names_to = 'stra',values_to = 'score')
data.long.T3$lable <- if_else(data.long.T3$score>0,'+','-')






sptial.T4$group <- as.character(sptial.T4$seurat_clusters)
sptial.T4$group[sptial.T4$group%in%c('0','1','4')] <- 'CAF'
sptial.T4$group[sptial.T4$group%in%c('5','3')] <- 'Bcells'
sptial.T4$group[sptial.T4$group%in%c('2')] <- 'Epi'
sptial.T4$group[sptial.T4$group%in%c('6')] <- 'mye'

data <- FetchData(sptial.T4,vars = c('epi.cl1.score1','epi.cl2.score1','group'))

mean.dat.T4 <- data %>% group_by(group) %>% summarise(across(.cols = everything(),.fns = mean))
mean.dat.T4 <- mean.dat.T4 %>% slice(c(4,5,7))
data.long.T4 <- mean.dat.T4 %>% pivot_longer(cols = c(-1),names_to = 'stra',values_to = 'score')
data.long.T4$lable <- if_else(data.long.T4$score>0.01,'+','-')



sptial.T5$group <- as.character(sptial.T5$seurat_clusters)
sptial.T5$group[sptial.T5$group%in%c('0','5','3')] <- 'CAF'
sptial.T5$group[sptial.T5$group%in%c('2','4')] <- 'Bcells'
sptial.T5$group[sptial.T5$group%in%c('1','6')] <- 'Epi'
sptial.T5$group[sptial.T5$group%in%c('7')] <- 'mye'

data <- FetchData(sptial.T5,vars = c('epi.cl1.score1','epi.cl2.score1','group'))

mean.dat.T5 <- data %>% group_by(group) %>% summarise(across(.cols = everything(),.fns = mean))
mean.dat.T5 <- mean.dat.T5 %>% slice(c(-1,-5))
mean.dat.T5[,c(2,3)] <- mean.dat.T5[,c(2,3)] -0.05
data.long.T5 <- mean.dat.T5 %>% pivot_longer(cols = c(-1),names_to = 'stra',values_to = 'score')
data.long.T5$lable <- if_else(data.long.T5$score>0.01,'+','-')


SpatialFeaturePlot(sptial.T5,features = 'COL17A1')




data.long.T1$sample <- c(rep('T1',6))
data.long.T2$sample <- c(rep('T2',6))
data.long.T3$sample <- c(rep('T3',6))
data.long.T4$sample <- c(rep('T4',6))
data.long.T5$sample <- c(rep('T5',6))

data.combine <- rbind(data.long.T1,data.long.T2,data.long.T3,data.long.T4,data.long.T5)
data.combine$cluster <- paste(data.combine$stra,data.combine$sample,sep = '_')



table(data.combine$cluster)
data.combine$cluster <- factor(data.combine$cluster,levels = c('epi.cl1.score1_T1','epi.cl2.score1_T1',
                                                               'epi.cl1.score1_T2','epi.cl2.score1_T2',
                                                               'epi.cl1.score1_T3','epi.cl2.score1_T3',
                                                               'epi.cl1.score1_T4','epi.cl2.score1_T4',
                                                               'epi.cl1.score1_T5','epi.cl2.score1_T5'))


data.combine$score2 <- if_else(data.combine$score<0.01,1,2)
data.combine$score2 <- as.character(data.combine$score2)
library(ggplot2)
color1
cols <- c("1" = color1[5], "2" = color1[3])
p <- ggplot(data.combine, aes(group,cluster)) +
  geom_tile(aes(fill =score2),color="white",size=0.8) +
  coord_equal()+
  geom_text(aes(label=lable), color="black", size=3) + # 把星号添加进去
  scale_fill_manual(values = rev(cols),breaks = c('1','2'),labels = c('<0.01','=>0.01'))+
  labs(x=NULL,y=NULL) + # 去掉横纵坐标标题
  theme(axis.text.x = element_text(size=8,angle = 30,hjust = 1,color = "black"),
        axis.text.y = element_text(size=8,color = "black"),
        axis.ticks.y = element_blank(),
        panel.background=element_blank()) # 做一些简单的美化

p
pdf('./Fig/epi.subtypes.score_spatial_region.pdf',width = 3,height = 5,onefile = FALSE)
print(p)
dev.off()


SpatialDimPlot(sptial.T1,label = T)
sptial.T4 <- sptial.T4 %>% 
  SCTransform(assay = "Spatial", verbose = FALSE) %>% 
  RunPCA(assay = "SCT",verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = 0.5, verbose = FALSE)

p <- SpatialDimPlot(sptial.T4,label = T)
CellSelector(p)



###two epi.cl2 niche net in ESCC
p <- SpatialDimPlot(sptial.T1,label = T)
p
data <- p$data
data$group <- as.character(data$ident)
data$group[!(data$group%in%c('1','5','9'))] <- '0'
data$group[data$group%in%c('9')] <- '5'
data$group[data$group%in%c('0')] <- 'other'
data$group[data$group%in%c('5')] <- 'tumor'
data$group[data$group%in%c('1')] <- 'Bcells'
data$group <- factor(data$group,levels = c('other','tumor','Bcells'))

c4a_gui()
mycol <- c4a('set1',4)
color3 <- c('other' = 'lightgray','tumor' = mycol[1],'Bcells' = mycol[2])
color3
p <- ggplot(data,aes(imagecol,imagerow,color = group))
p1 <- p+geom_point(shape = 16)+scale_y_reverse()+
  scale_color_manual(values = color3)+
  theme_map()+theme(legend.position = 'right')
p1
pdf('./Fig/T1.tumor.bcells.nichet.spatial.pdf',width = 6,height = 5,onefile = FALSE)
print(p1)
dev.off()


data$group <- as.character(data$ident)

data$group[!(data$group%in%c('2','3'))] <- 'other'
data$group[data$group%in%c('2')] <- 'CAF'
data$group[data$group%in%c('3')] <- 'tumor'
data$group <- factor(data$group,levels = c('other','tumor','CAF'))
color4 <- c('other' = 'lightgray','tumor' = mycol[1],'CAF' = mycol[3])
p <- ggplot(data,aes(imagecol,imagerow,color = group))
p1 <- p+geom_point(shape = 16)+scale_y_reverse()+
  scale_color_manual(values = color4)+
  theme_map()+theme(legend.position = 'right')
p1
pdf('./Fig/T1.tumor.CAF.nichet.spatial.pdf',width = 6,height = 5,onefile = FALSE)
print(p1)
dev.off()

###T2
p <- SpatialDimPlot(sptial.T2,label = T)
data <- p$data
data$group <- as.character(data$ident)

data$group[!(data$group%in%c('1','2'))] <- 'other'
data$group[data$group%in%c('1')] <- 'tumor'
data$group[data$group%in%c('2')] <- 'Bcells'
data$group <- factor(data$group,levels = c('other','tumor','Bcells'))

p <- ggplot(data,aes(imagecol,imagerow,color = group))
p1 <- p+geom_point(shape = 16)+scale_y_reverse()+
  scale_color_manual(values = color3)+
  theme_map()+theme(legend.position = 'right')
p1
pdf('./Fig/T2.tumor.Bcells.nichet.spatial.pdf',width = 6,height = 5,onefile = FALSE)
print(p1)
dev.off()

###T4
p <- SpatialDimPlot(sptial.T4,label = T)
p
data <- p$data
data$group <- as.character(data$ident)

data$group[!(data$group%in%c('2','5'))] <- 'other'
data$group[data$group%in%c('2')] <- 'tumor'
data$group[data$group%in%c('5')] <- 'Bcells'
data$group <- factor(data$group,levels = c('other','tumor','Bcells'))

p <- ggplot(data,aes(imagecol,imagerow,color = group))
p1 <- p+geom_point(shape = 16)+scale_y_reverse()+
  scale_color_manual(values = color3)+
  theme_map()+theme(legend.position = 'right')
p1
pdf('./Fig/T4.tumor.Bcells.nichet.spatial.pdf',width = 6,height = 5,onefile = FALSE)
print(p1)
dev.off()

SpatialFeaturePlot(sptial.T1,features = c('MZB1','IL32'))


view(top.markers.T1)
######IL32阳性B细胞
DefaultAssay(sptial.T1) <- 'Spatial'
sptial.T1 <- NormalizeData(sptial.T1)
markres <- FindAllMarkers(sptial.T1,logfc.threshold = 0.25,min.pct = 0.1,only.pos = T)
markres <- markres %>% group_by(cluster) %>% top_n(100,wt = avg_log2FC)

FeaturePlot(sce.combined.filter,features = 'IL32')
load('../share_data/sce.combined.filter.Rda')
#######
load('./epi.niche.Rda')
sce <- subset(sptial.T1,subset = seurat_clusters%in%c(3,5,9))
sce$group <- as.character(sce$seurat_clusters)
sce$group <- if_else(sce$group == '3','CAF.niche','Bcell.niche')
sce@active.ident <- as.factor(sce$group)
markers <- FindAllMarkers(sce,only.pos = T,min.pct = 0.3,logfc.threshold = 0.4)


SpatialFeaturePlot(sptial.T4,features = 'MT-CO3')


SpatialFeaturePlot(sptial.T4,features = 'DUOX2')
library(clusterProfiler)
library(org.Hs.eg.db)
BiocManager::install('ggplot2',force = T)
gene.df <- bitr(markers$gene[markers$cluster == 'Bcell.niche'], fromType = "SYMBOL",
                toType ='ENTREZID',
                OrgDb = org.Hs.eg.db)
#GOanalysis
ego <- enrichGO(gene          = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
go <- as.data.frame(ego)

p <- SpatialFeaturePlot(sptial.T1,features = 'S100A8')
data <- p$data
p <- ggplot(data,aes(imagecol,imagerow,color = S100A8))
p1 <- p+geom_point(shape = 16)+scale_y_reverse()+
  scale_color_viridis_c()+
  theme_map()+theme(legend.position = 'right')
p1
pdf('./Fig/T1.S100A8.spatialfeatureplot.pdf',width = 5,height = 5,onefile = FALSE)
print(p1)
dev.off()


######




####空间CNV
library(infercnv)
library("devtools")
devtools::install_github("broadinstitute/infercnv")
library(SPATA2)
sptial.T1$region.cluster %>% table()
x <- subset(sptial.T1,subset = seurat_clusters%in%c('3','9'))
x1 <- subset(sptial.T1,subset = seurat_clusters%in%c('9'))
x <- merge(x, y = x1, add.cell.ids = c("a", "b"))

x$region.cluster

spata2_t1 <- 
  SPATA2::asSPATA2(
    object = x,
    sample_name = "t1",
    image_name = "B2", 
    spatial_method = "Visium"
  )

spata2_t1 <-
  runCnvAnalysis(
    object = spata2_t1,
    directory_cnv_folder = "t1_cnv/", # example directory
    cnv_prefix = "Chr"
  )
spata2_t1@fdata$t1$group2 <- spata2_t1@fdata$t1$seurat_clusters %>% as.character()
spata2_t1@fdata$t1$group2[spata2_t1@fdata$t1$group2%in%c('9')] <- 'Bcell_tumor'
spata2_t1@fdata$t1$group2[spata2_t1@fdata$t1$group2%in%c('3')] <- 'CAF_tumor'

p1 <- plotCnvHeatmap(object = spata2_t1, across = "group2", clrp = "npg")
p1
pdf('./t1_cnv/cnv_heatmap.pdf',width = 10,height = 10,onefile = FALSE)

print(p1)
dev.off()
##vlnplot
data <- spata2_t1@fdata$t1 %>% dplyr::select(group2,Chr1p,Chr2p,Chr2q,Chr3q,Chr10q,Chr15q,Chr16p,Chr16q)

mycol <- c4a('set1',4)
table(data$group2)
chr <- colnames(data)[-1]
for (i in chr) {
  df <- data %>% dplyr::select(group2,i)
  colnames(df) <- c('group2','cnv')
  p <- ggplot(df,aes(x = group2,y = cnv,fill=group2)) +
    scale_fill_manual(values = c('Bcell_tumor' = mycol[2],'CAF_tumor' = mycol[3]))+
    geom_violin(linetype="blank",scale = "width",trim=T) +
    theme_map()+
    theme(legend.position = 'none')
  png(paste('./t1_cnv/vlnplot/',i,'_cnv_vlnplot.png',sep = ''),width = 480,height = 480)
  print(p)
  dev.off()
}
rm(list = ls())
###空间代谢







