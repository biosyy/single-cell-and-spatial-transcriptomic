rm(list = ls())
Sys.setenv(LANGUAGE = "en")
options("repos" = c(CRAN="http://mirrors.ustc.edu.cn/CRAN/"))
library(tidyverse)
library(ggplot2)
library(seriation)
library(ggpubr)
library(RColorBrewer)
library(magick)
library(ComplexHeatmap)
library(circlize)
library(ggthemes)
library(ggh4x)
library(scales)
library(Seurat)
load('../share_data/sce.combined.sct.Rda')
load('../share_data/marker.Rda')
FeaturePlot(sce.combined.sct,features = c('KRT18', 'EPCAM', 'CD24'))
sce.combined.sct@active.ident = as.factor(sce.combined.sct$cluster)
DefaultAssay(sce.combined.sct) <- "RNA"
DimPlot(sce.combined.sct)
sce.combined.sct <- ScaleData(sce.combined.sct)

###seurat cluster
pdf('seurat_cluster_umap.pdf',width = 7,height = 5,onefile = FALSE)
DimPlot(sce.combined.sct,repel = TRUE,label = T,pt.size = 0.5)+theme_few()
dev.off()
###QC
data <- FetchData(sce.combined.sct,vars = c('sample','percent.mt','nCount_RNA','nFeature_RNA'))
pdf('QC_percent.mt.pdf',width = 5,height = 5,onefile = FALSE)
ggviolin(data, x="sample", y="percent.mt", fill = "sample",
         
         palette = 'npg',
         
         add = "mean_se", add.params = list(fill="white"),
         ggtheme = theme_pander())
dev.off()
data$nCount_RNA <- log1p(data$nCount_RNA)
pdf('QC_log2nCount_RNA.pdf',width = 5,height = 5,onefile = FALSE)
ggviolin(data, x="sample", y="nCount_RNA", fill = "sample",
         
         palette = 'npg',
         
         add = "mean_se", add.params = list(fill="white"),
         ggtheme = theme_pander())
dev.off()
data$nFeature_RNA <- log1p(data$nFeature_RNA)
pdf('QC_log2nFeature_RNA.pdf',width = 5,height = 5,onefile = FALSE)
ggviolin(data, x="sample", y="nFeature_RNA", fill = "sample",
         
         palette = 'npg',
         
         add = "mean_se", add.params = list(fill="white"),
         ggtheme = theme_pander())
dev.off()

###featureplot
genes <- c('KRT18', 'EPCAM', 'CD24','CD3D', 'CD3E', 'CD8A','CD14','C1QA','ICAM1',
           'PECAM1', 'CLDN5', 'CD34','DCN', 'OGN', 'ACTB','CD79A', 'JCHAIN','IGHG1',
           'TPSB2','TPSAB1','PTPRC')
path <- './featurePlot/'
for (i in genes) {
  p <- FeaturePlot(sce.combined.sct,features = i,
                   max.cutoff = 3,pt.size = 0.1)+
    scale_color_viridis_c()+theme_few()
  name <- paste(i,'featurePlot.pdf',sep = '_')
  pdf(paste(path,name,sep = ''),width = 6,height = 5,onefile = FALSE)
  print(p)
  dev.off()
}
rm(sce.combined.sct)


dev.off()
?FeaturePlot

###main cluster umap
pdf('main_cluster_umap.pdf',width = 7,height = 5,onefile = FALSE)
DimPlot(sce.combined.sct,group.by = 'cluster',repel = TRUE,label = T,pt.size = 0.5)+theme_few()
dev.off()
table(sce.combined.sct$sample)
###umap labeled by sample
display.brewer.all()
color2 <- brewer.pal(6, "Set1")
p = DimPlot(sce.combined.sct,group.by = 'sample',repel = TRUE,
        pt.size = 0.8,cols = color2)+theme_few()
dev.off()
pdf('main_cluster_sample.pdf',width = 7,height = 5,onefile = FALSE)
print(p)
dev.off()

### cell counts
#cell number
df <- table(sce.combined.sct@meta.data$sample) %>% as.data.frame()
colnames(df) <- c('sample','cellnum')
df <- df[order(df$cellnum),]
display.brewer.all() #展示所有色板
brewer.pal.info
color1 <- brewer.pal(3, "BuPu") #取颜色
show_col(color1) #可视化

p <- ggplot(df,aes(reorder(sample,cellnum),cellnum))+
  geom_col(width = 0.5,fill = color1[2])+
  coord_flip()+
  theme_classic2()+
  scale_y_continuous(position = 'right')+
  guides(x = guide_axis_truncated(trunc_lower = 0, trunc_upper = 22053), y="axis_truncated")
print(p)
dev.off()
pdf('sample_cellnum.pdf',width = 7,height = 5,onefile = FALSE) 
print(p)
dev.off()

####cell percent in different sample
sample_cat <- c('T2','N1','T4','N4','T1','T3')
color2 <- brewer.pal(7, "Set1")
df1 <- sce.combined.sct@meta.data %>% select('sample','cluster')
df1 <- df1 %>% group_by(sample,cluster) %>% summarise(n = n())
df1 <- df1 %>% group_by(sample) %>% summarise(cluster = cluster,percent = n/sum(n))
df1$cluster <- factor(df1$cluster)#固定顺序
df1$cluster
df1$sample <- factor(df1$sample,levels = rev(sample_cat))#固定顺序
p = ggplot(df1, aes( x = sample,y= percent,fill = cluster))+
  geom_col(position = 'stack', width = 0.5)+
  coord_flip()+
  theme_classic2()+
  scale_y_continuous(position = 'right')+
  guides(x = "axis_truncated", y="axis_truncated")+
  scale_fill_manual(values = color2 )+
  theme(legend.position = 'none')
p

dev.off()
pdf('main_cluster_sample_percent.pdf',width = 7,height = 5,onefile = FALSE) 
print(p)
dev.off()
p = ggplot(df1, aes( x = sample,y= percent,fill = cluster))+
  geom_col(position = 'stack', width = 0.5)+
  coord_flip()+
  theme_classic2()+
  scale_y_continuous(position = 'right')+
  guides(x = "axis_truncated", y="axis_truncated")+
  scale_fill_manual(values = color2)
p
sce.combined.sct@assays$RNA@scale.data
load('../share_data/sce.marker.Rda')

##marker heatmap top 100
load('./htpdata.Rda')
library(circlize)
color = colorRampPalette(colors = c("#1BB7E5","#000000","#FFFF00"))(11)
col_fun = colorRamp2(seq(-2, 2, len = 11), color)

p <- Heatmap(data,cluster_rows = F,cluster_columns = F,raster_by_magick = TRUE,
             raster_magick_filter = 'Undefined',
             use_raster = TRUE,show_row_names = F,
             show_column_names = F,
             col = col_fun,raster_quality = 20,
             heatmap_legend_param = list(title = 'scale.data',
                                         at = c(-2,0,2),
                                         labels = c(-2,0,2)),
             top_annotation = ha,
             left_annotation = ha_row
)
pdf('./fig/top100genes.heatmap.pdf',width = 8,height = 8,onefile = FALSE)
print(p)
dev.off()

###kegg##
library(clusterProfiler)
library(org.Hs.eg.db)
load('./100.DEG.main.Rda')


Gene_ID <- bitr(genes100$gene, fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")

#构建文件并分析
data  <- merge(Gene_ID,genes100,by.x='SYMBOL',by.y='gene')
data_GO <- compareCluster(
  ENTREZID~cluster, 
  data=data, 
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

data_GO_sim <- simplify(data_GO, 
                        cutoff=0.7, 
                        by="p.adjust", 
                        select_fun=min)
data_GO_sim_fil <- data_GO_sim@compareClusterResult

go <- data_GO_sim_fil %>% dplyr::group_by(cluster) %>% top_n(5,wt = -log10(p.adjust))
go$Description <- factor(go$Description,levels = go$Description)
bk <- c("Epithelial", "Endothelial", "Fibroblasts",
        'T-cells','B-cells','Myeloids','Mast-cells','MKI67.immune.cells')
go$Cluster <- factor(go$Cluster,levels = bk)

p = ggplot(go,aes(x=Cluster,y=Description)) + 
  geom_point(aes(size=Count,color=-1*log10(p.adjust)))+
  scale_color_viridis(option = 'C')+
  scale_size_continuous(range = c(4,12))+
  labs(
    color=expression(-log[10](p.adjust)),
    size="Gene number",
    x="group"
    # y="Pathway name",
    # title="Pathway enrichment")
  )+
  coord_flip()+
  theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank())
pdf('./fig/top100DEG.GO.dotplot.pdf',width = 17,height = 9,onefile = FALSE)
print(p)
dev.off()
###vlnplot(top 5 marker)
vln.dat=FetchData(sce.combined.sct,c(anno_marker,"cluster"))
vln.dat$cluster <- as.factor(vln.dat$cluster)
levels(vln.dat$cluster)
vln.dat.melt=vln.dat %>% 
  reshape2::melt(anno_marker) %>%
  rename("Gene"="variable") %>%
  group_by(cluster,Gene) %>%
  mutate(fillcolor=mean(value))
#颜色设置
pal=colorRampPalette(RColorBrewer::brewer.pal(n = 9,name="YlOrRd"))(100)

p1 = ggplot(vln.dat.melt,aes(x=Gene,y=value,fill=fillcolor))+
  # 把小提琴图的外缘轮廓去除
  geom_violin(linetype="blank",scale = "width")+
  scale_fill_gradientn(colors=pal,name="Average Expression")+
  # 分面
  facet_grid(cluster~.)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5),
        strip.text.y = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        legend.position = "left"
        
  )


####ST data
library(cols4all)
c4a_gui()
mycol <- c4a('rainbow',7)
load('../Fig2_analysis/ST.data.Rda')
p <- SpatialDimPlot(sptial.T1,alpha = 0)+theme(legend.position = 'none')

p1 <- SpatialDimPlot(sptial.T1)+ scale_fill_discrete(name="Cluster")+scale_fill_manual(values = mycol)

p2 <- SpatialDimPlot(sptial.T2,alpha = 0)+theme(legend.position = 'none')
p3 <- SpatialDimPlot(sptial.T2)+ scale_fill_discrete(name="Cluster")+scale_fill_manual(values = mycol)
mycol
p4 <- SpatialDimPlot(sptial.T4,alpha = 0)+theme(legend.position = 'none')
p5 <- SpatialDimPlot(sptial.T4)+ scale_fill_discrete(name="Cluster")+scale_fill_manual(values = mycol)

library(patchwork)
p6 <- p+p1+p2+p3+p4+p5+p4+p5+p4+p5+patchwork::plot_layout(design = '
                                        02468
                                        13579')
p6
pdf('./fig/combine.spatialdimplot.pdf',width = 25,height = 10,onefile = FALSE)
print(p6)
dev.off()


rm(list = ls())












