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
library(SPATA2)
load('./sce.lym.Rda')



library(cols4all)
c4a_gui()
mycol <- c4a('set2',4)
sce.lym <- subset(sce.lym,subset = cluster != 'PDE3B_CD8T')
sce.lym$cluster1 <- sce.lym$cluster
table(sce.lym$cluster)
sce.lym$cluster1[sce.lym$cluster1 == 'CPEB4_plasma'] <- 'Bcells'
sce.lym$cluster1[sce.lym$cluster1 == 'JCHAIN_plasma'] <- 'Bcells'
sce.lym$cluster1[sce.lym$cluster1 == 'naive_Bcells'] <- 'Bcells'


p <- DimPlot(sce.lym,group.by = 'cluster1')+scale_color_manual(values = mycol)+
  scale_color_discrete(name="celltype",
                      breaks=c("Bcells", "CD4T", "CD8T",'MKI67.lym'),
                      labels=c("Bcells (n = 12,239)", "CD4T (n = 7,119)", "CD8T (n = 8,790)",'MKI67.lym (n = 1,043)'))

p
pdf('./fig/lym_main_celltype.pdf',width = 7,height = 5,onefile = FALSE)
print(p)
dev.off()

load('./sce.CD4T.Rda')
DimPlot(sce.CD4T.SCT)
table(sce.CD4T.SCT$cluster)
sce.CD4T.SCT <- subset(sce.CD4T.SCT,subset = cluster%in%c('CCR7_CD4T','FOXP3_CD4T','GZMA_CD4T'))
library(cols4all)
c4a_gui()
mycol <- c4a('bold',3)
library(scales)
show_col(mycol)
p1 <- DimPlot(sce.CD4T.SCT,group.by = 'cluster')+
  scale_color_discrete(name="celltype",
                       breaks=c("CCR7_CD4T", "FOXP3_CD4T", "GZMA_CD4T"),
                       labels=c("CCR7_CD4T (n = 3,068)", "FOXP3_CD4T (n = 2,836)", "GZMA_CD4T (n = 733)"),type = mycol)
p1
DimPlot(sce.CD4T.SCT,group.by = 'cluster')+scale_color_manual(values = mycol)

pdf('./fig/CD4T_subcelltype.pdf',width = 7,height = 5,onefile = FALSE)
print(p1)
dev.off()


load('./sce.CD8TandNK.Rda')
DimPlot(sce.CD8T.SCT)
table(sce.CD8T.SCT$cluster)
sce.CD8T.SCT <- subset(sce.CD8T.SCT,subset = cluster != '11')
library(cols4all)
c4a_gui()
mycol <- c4a('pastel',8)
library(scales)
show_col(mycol)
table(sce.CD8T.SCT$cluster)
p2 <- DimPlot(sce.CD8T.SCT,group.by = 'cluster')+
  scale_color_discrete(name="celltype",
                       breaks=c("AREG_NK", "CCL4_CD8T", "CCR7_CD8T",'CXCL13_CD8T','FCGR3A_NK','GNLY_CD8T','GZMK_CD8T','IFNG_CD8T'),
                       labels=c("AREG_NK (n = 585)", "CCL4_CD8T (n = 284)", "CCR7_CD8T (n = 472)",
                                'CXCL13_CD8T (n = 2,062)','FCGR3A_NK (n = 526)',
                                'GNLY_CD8T (n = 415)','GZMK_CD8T (n = 1,411)','IFNG_CD8T (n = 2,784)'),type = mycol)
p2
DimPlot(sce.CD8T.SCT,group.by = 'cluster')+scale_color_manual(values = mycol)

pdf('./fig/CD8T_subcelltype.pdf',width = 7,height = 5,onefile = FALSE)
print(p2)
dev.off()



sce.Bcell <- subset(sce.lym,subset = cluster%in%c('Bcells','naive_Bcells','CPEB4_plasma','JCHAIN_plasma'))
rm(sce.lym)

p <- DimPlot(sce.Bcell,group.by = 'cluster')
id <- CellSelector(p)
sce.Bcell <- sce.Bcell[,id]
p
FeaturePlot(sce.Bcell,features = 'LPP')

sce.Bcell@active.ident <- as.factor(sce.Bcell$cluster)

markes <- FindAllMarkers(sce.Bcell,logfc.threshold = 0.25,min.pct = 0.1,only.pos = T)
top30 <- markes %>% group_by(cluster) %>% top_n(30,wt = avg_log2FC)
sce.Bcell$cluster1 <- sce.Bcell$cluster
sce.Bcell$cluster1[sce.Bcell$cluster1 == 'Bcells'] <- 'imature_Bcells'

sce.Bcell$cluster1[sce.Bcell$cluster1 == 'naive_Bcells'] <- 'LPP_Bcells'

p <- DimPlot(sce.Bcell,group.by = 'seurat_clusters')
p
c4a_gui()
mycol <- c4a('vivid',4)
library(scales)
show_col(mycol)
table(sce.Bcell$cluster1)
p3 <- DimPlot(sce.Bcell,group.by = 'cluster1')+
  scale_color_discrete(name="celltype",
                       breaks=c("CPEB4_plasma", "imature_Bcells", "JCHAIN_plasma",'LPP_Bcells'),
                       labels=c("CPEB4_plasma (n = 4,470)", "imature_Bcells (n = 4,316)", "JCHAIN_plasma (n = 2,649)",
                                'LPP_Bcells (n = 788)'),type = mycol)
p3


pdf('./fig/Bcells_subcelltype.pdf',width = 7,height = 5,onefile = FALSE)
print(p3)
dev.off()
rm(list = ls())

load('./sce.lym.Rda')
DimPlot(sce.lym,group.by = 'cluster')
sce.pro <- subset(sce.lym,subset = cluster == 'MKI67.lym')
DefaultAssay(sce.pro) <- 'integrated'
rm(sce.lym)
sce.pro <- sce.pro %>% FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 0.3)
DimPlot(sce.pro)

sce.pro$cluster <- sce.pro$seurat_clusters %>% as.character()
sce.pro$cluster[sce.pro$cluster == '2'] <- 'MKI_CD8T'
sce.pro$cluster[sce.pro$cluster == '3'] <- 'MKI_CD4T'
sce.pro$cluster[sce.pro$cluster == '4'] <- 'MKI_LPP_Bcells'
sce.pro$cluster[sce.pro$cluster%in%c('0','1')] <- 'MKI_Imature_Bcells'
DimPlot(sce.pro,group.by = 'cluster')
table(sce.pro$cluster)
c4a_gui()
mycol <- c4a('safe',4)
p4 <- DimPlot(sce.pro,group.by = 'cluster')
p4
p4 <- DimPlot(sce.pro,group.by = 'cluster')+
  scale_color_discrete(name="celltype",
                       breaks=c("MKI_CD4T", "MKI_CD8T", "MKI_Imature_Bcells",'MKI_LPP_Bcells'),
                       labels=c("MKI_CD4T (n = 155)", "MKI_CD8T (n = 214)", "MKI_Imature_Bcells (n = 525)",
                                'MKI_LPP_Bcells (n = 149)'),type = mycol)
p4
pdf('./fig/MKI67_lym_subtypes.pdf',width = 7,height = 5,onefile = FALSE)
print(p4)
dev.off()

save(sce.pro,file = 'sce.pro.Rda')

####cellpercent (Bcells, CD8T, CD4T)
sample_cat <- c('N1','T4','N5','T1','T2','T3','N4','T5')
library(RColorBrewer)
library(ggpubr)
rm(sce.lym)
DimPlot(sce.CD8T.SCT,group.by = 'cluster')
sce.CD8T.SCT <- subset(sce.CD8T.SCT,cluster != '11')
table(sce.Bcell$cluster)
df1 <- sce.Bcell@meta.data %>% select('sample','cluster')
df1 <- df1 %>% group_by(cluster,sample) %>% summarise(n = n())
df1 <- df1 %>% group_by(cluster) %>% summarise(sample = sample,percent = n/sum(n))
df1$cluster <- factor(df1$cluster)#固定顺序
df1$cluster
df1$sample <- factor(df1$sample,levels = rev(sample_cat))#固定顺序
#c4a_gui()
mycol <- c4a('colorblind',8)

p = ggplot(df1, aes( x = cluster,y= percent,fill = sample))+
  geom_col(position = 'stack', width = 0.5,color = 'black')+
  coord_flip()+
  theme_classic2()+
  scale_y_continuous(position = 'right')+
  guides(x = "axis_truncated", y="axis_truncated")+
  scale_fill_manual(values = mycol)+labs(x = NULL)
p
library(cowplot)
library(ggpubr)
legends = as_ggplot(get_legend(p))

df2 <- sce.Bcell@meta.data %>% select('sample','cluster')
df2$group <- if_else(df2$sample%in%c('N1','T4','N5'),'adj','tumor')

df2 <- df2 %>% group_by(cluster,group) %>% summarise(n = n())
df2 <- df2 %>% group_by(cluster) %>% summarise(group = group,percent = n/sum(n))
mycol <- c4a('set1',2)
p1 = ggplot(df2, aes( x = cluster,y= percent,fill = group))+
  geom_col(position = 'stack', width = 0.5,color = 'black')+
  coord_flip()+
  theme_classic2()+
  scale_y_continuous(position = 'right')+
  guides(x = "axis_truncated", y="axis_truncated")+
  scale_fill_manual(values = mycol)+labs(x = NULL)
p1


df3 <- sce.Bcell@meta.data %>% select('sample','cluster')
df3 <- df3 %>% filter(sample%in%c('T1','T2','T3','N4','T5'))
df3$group <- if_else(df3$sample%in%c('T1','T3','T5'),'chemo','pri')

df3 <- df3 %>% group_by(cluster,group) %>% summarise(n = n())
df3 <- df3 %>% group_by(cluster) %>% summarise(group = group,percent = n/sum(n))
mycol <- c4a('set3',3)
p2 = ggplot(df3, aes( x = cluster,y= percent,fill = group))+
  geom_col(position = 'stack', width = 0.5,color = 'black')+
  coord_flip()+
  theme_classic2()+
  scale_y_continuous(position = 'right')+
  guides(x = "axis_truncated", y="axis_truncated")+
  scale_fill_manual(values = mycol)+labs(x = NULL)
p2 <- p2+ theme(axis.line.y = element_blank(),
           axis.text.y  = element_blank(),
           axis.ticks.y = element_blank())
p1 <-  p1+ theme(axis.line.y = element_blank(),
                 axis.text.y  = element_blank(),
                 axis.ticks.y = element_blank())

legends <- as_ggplot(get_legend(p))
legends1 <-   as_ggplot(get_legend(p1))
legends2 <-   as_ggplot(get_legend(p2))

legends <- legends+legends1+legends2+plot_layout(design = '12
                                                           13')
p <- p+theme(legend.position = 'none')
p1 <- p1+theme(legend.position = 'none')
p2 <- p2+theme(legend.position = 'none')
combineP <- p+p1+p2+legends+plot_layout(design = '1234')
combineP
pdf('./fig/Bcells_subtypes_percent.pdf',width = 15,height = 5,onefile = FALSE)
print(combineP)
dev.off()



rm(list = ls())

VlnPlot(sce.combined.sct,features = 'IGHG1',group.by = 'cluster')
sce.combined.sct$group <- if_else(sce.combined.sct$sample%in%c('T1','T3'),'chemo',if_else(
  sce.combined.sct$sample%in%c('T2','N4'),'pri','adj'
))
VlnPlot(sce.combined.sct,features = 'IGHG1',group.by = 'cluster',split.by = 'group')
load('../Fig2_analysis/GSE160269.epi.Rda')

DimPlot(adata2)
adata2 <- subset(adata2,subset = sample == 'P39T')
VlnPlot(adata2,features = 'IGHG1',group.by = 'cluster')

rm(list = ls())


