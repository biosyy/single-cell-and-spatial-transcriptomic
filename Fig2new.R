
options(stringsAsFactors = F)
options("repos" = c(CRAN="http://mirrors.ustc.edu.cn/CRAN/"))
Sys.setenv(LANGUAGE= "en")
library(tidyverse)
library(Seurat)
library(patchwork)
library(parallel)
library(ggthemes)
library(magick)
library(cols4all)
load('./sce.epi.Rda')
######FigA epi subtypes
p <- DimPlot(sce.epi)
epiid <- CellSelector(p)
sce.epi.filter <- sce.epi[,epiid]
sce.epi.filter <- subset(sce.epi.filter,subset = PTPRC < 1)
DefaultAssay(sce.epi.filter) <- 'integrated'
sce.epi.filter <-sce.epi.filter %>%  
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = 0.2)

DimPlot(sce.epi.filter,label = T)
sce.epi.filter$cluster <- sce.epi.filter$seurat_clusters %>% as.character()
sce.epi.filter$cluster <- if_else(sce.epi.filter$cluster%in%c('2','4','5'),'cl1','cl2')
c4a_gui()
mycol <- c4a('egypt',2)
p1 <- DimPlot(sce.epi.filter,pt.size = 0.5,group.by = 'cluster')+scale_color_manual(values = mycol)+theme_map()+
  theme(legend.position = 'right',plot.title = element_text(hjust = 0.5))
p1
mycol <- c4a('bold',8)
p2 <- DimPlot(sce.epi.filter,group.by = 'sample',pt.size = 0.5)+scale_color_manual(values = mycol)+theme_map()+
  theme(legend.position = 'right',plot.title = element_text(hjust = 0.5))
p1
p2
###epi development(利用integrated数据进行分析得到非常好得结果！！！！！！)
#rwa counts
library(monocle)
pd <- sce.epi.filter@meta.data
fd <- data.frame(gene_short_name = rownames(sce.epi.filter))
rownames(fd) <- rownames(sce.epi.filter)
pd <- new("AnnotatedDataFrame", data = pd)
fd <- new("AnnotatedDataFrame", data = fd)
sce <- newCellDataSet(sce.epi.filter@assays$RNA@counts,
                      phenoData = pd,
                      featureData = fd,
                      expressionFamily=negbinomial.size())
sce <- estimateSizeFactors(sce)
sce <- estimateDispersions(sce)
##Trajectory step 1: choose genes that define a cell's progress
disp_table <- dispersionTable(sce)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.2 &
                           dispersion_empirical >= 1 * dispersion_fit)$gene_id
sce <- setOrderingFilter(sce, ordering_genes)
plot_ordering_genes(sce)
###Trajectory step 2: reduce data dimensionality
sce <- reduceDimension(sce, max_components = 2,
                       method = 'DDRTree')

##Trajectory step 3: order cells along the trajectory
sce <- orderCells(sce)
p3 <- plot_cell_trajectory(sce, color_by = "cluster",show_tree = F,
                     show_branch_points = F)
p3

data <- p1$data
data$Pseudotime <- p3$data$Pseudotime
library(RColorBrewer)
display.brewer.all()
color1 <- brewer.pal(8, "YlOrBr")
p4 <- ggplot(data,aes(x= UMAP_1 , y = UMAP_2 ,color = rev(Pseudotime))) +
  geom_point(size = 0.5 , alpha =1)+
  scale_color_gradientn(colors = color1)+theme_map()+theme(legend.position = 'bottom')
p4
#integrated data slot
pd <- sce.epi.filter@meta.data
fd <- data.frame(gene_short_name = sce.epi.filter@assays$integrated@var.features)
rownames(fd) <- sce.epi.filter@assays$integrated@var.features
pd <- new("AnnotatedDataFrame", data = pd)
fd <- new("AnnotatedDataFrame", data = fd)
counts <- sce.epi.filter@assays$integrated@data %>% as.matrix()
counts[counts<0] <- 0
sce <- newCellDataSet(counts,
                      phenoData = pd,
                      featureData = fd)
sce <- estimateSizeFactors(sce)
sce <- estimateDispersions(sce)
##Trajectory step 1: choose genes that define a cell's progress
disp_table <- dispersionTable(sce)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.2 &
                           dispersion_empirical >= 1 * dispersion_fit)$gene_id
sce <- setOrderingFilter(sce, ordering_genes)
plot_ordering_genes(sce)
###Trajectory step 2: reduce data dimensionality
sce <- reduceDimension(sce, max_components = 2,
                       method = 'DDRTree')

##Trajectory step 3: order cells along the trajectory
sce <- orderCells(sce)
p3 <- plot_cell_trajectory(sce, color_by = "cluster",show_tree = F,
                           show_branch_points = F)
p3

data <- p1$data
data$Pseudotime <- p3$data$Pseudotime
library(RColorBrewer)
display.brewer.all()
color1 <- brewer.pal(8, "YlOrBr")
p4 <- ggplot(data,aes(x= UMAP_1 , y = UMAP_2 ,color = Pseudotime)) +
  geom_point(size = 0.5 , alpha =1)+
  scale_color_gradientn(colors = color1)+theme_map()+theme(legend.position = 'bottom')
p4


p5 <- p1+p2+p4+
  patchwork::plot_layout(design = "
                        AAB
                        AAC
                         ")#分成四份，横向排列
p5
pdf('./Fig/epi_subtypes_umap.pdf',width = 8,height = 5,onefile = FALSE)
print(p5)
dev.off()
###cell percent epi
data <- FetchData(sce.epi.filter,vars = c('sample','cluster'))

data <-data %>%  group_by(sample,cluster) %>% summarise(counts = n())

perc <- data %>% dplyr::group_by(sample) %>% summarise(cluster = cluster,percent = counts/sum(counts))
c4a_gui()
mycol <- c4a('egypt',2)
p = ggplot(perc, aes( x = sample,y= percent,fill = cluster))+
  geom_col(position = 'stack', width = 0.5)+
  scale_fill_manual(values = mycol)+
  theme_few()+
  scale_y_continuous(position = 'left')+
  theme(legend.position = 'none')
p
pdf('./Fig/epi_subtypes_percent.pdf',width = 5,height = 5,onefile = FALSE)
print(p)
dev.off()

###
data.t <- data %>% dplyr::filter(sample%in%c('T1','T2','T3','T5','N4'))
data.t$group <- if_else(data.t$sample%in%c('T1','T3','T5'),'chemo','pri')
data.t <- data.t %>% group_by(group,cluster) %>% summarise(counts = sum(counts))
perc.group <- data.t %>% dplyr::group_by(group) %>% summarise(cluster = cluster,percent = counts/sum(counts))


display.brewer.all()
c4a_gui()
mycol <- c4a('java',2)
p = ggplot(perc.group, aes( x = group,y= percent,fill = cluster))+
  geom_col(position = 'stack', width = 0.5)+
  theme_few()+
  scale_y_continuous(position = 'left')+
  scale_fill_manual(values = mycol)+
  theme(legend.position = 'right')
p
pdf('./Fig/epi_subtypes_cluster_percent.pdf',width = 3,height = 5,onefile =FALSE)
print(p)
dev.off()


top100.epi.markers <- epi.submarkers %>% group_by(cluster) %>% top_n(100,wt = avg_log2FC)
library(clusterProfiler)
library(org.Hs.eg.db)
top100.epi.markers$cluster %>% table()
genes <- top100.epi.markers %>% dplyr::filter(cluster == 'cl1') %>% pull(gene)
genes
gene.df <- bitr(genes, fromType = "SYMBOL",
                toType ='ENTREZID',
                OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene          = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
ego <- as.data.frame(ego)
go.cl1 <- ego

go.cl2 <- ego
FeaturePlot(sce.epi,features = 'CAV1')
DimPlot(sce.epi)
DefaultAssay(sce.epi) <- 'RNA'
sce.epi <-  sce.epi %>%  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = 0.2)
DimPlot(sce.epi,label = T)

DefaultAssay(sce.epi) <- 'RNA'
tmp <- subset(sce.epi,subset = cluster == 'cluster2')
load('./epi.subtypes.go.Rda')

go.dat <- rbind(go.cl1[c(1,2,3,5),],go.cl2[c(1,2,3,10),])
go.dat$group <- c(rep('cl1',4),rep('cl2',4))
go.dat$Description <- factor(go.dat$Description,levels = go.dat$Description)

library(viridis)
p = ggplot(go.dat,aes(x=group,y=Description)) + 
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
  theme_bw()
p
pdf('./Fig/Epi_GO_dotplot.pdf',width = 6,height = 5,onefile = FALSE)
print(p)
dev.off()
###some gene expression
sce.epi.filter@active.ident <- as.factor(sce.epi.filter$cluster)
DefaultAssay(sce.epi.filter) <- 'RNA'
epi.markers <- FindAllMarkers(sce.epi.filter,logfc.threshold = 0.25,min.pct = 0.1,only.pos = T)

FeaturePlot(sce.epi.filter,features = 'LCN2')
top100.epi.markers <- epi.markers %>% dplyr::group_by(cluster) %>% top_n(100,wt = avg_log2FC)
genes <- c('TP53','LCN2')
data <- sce.epi.filter@assays$RNA@data[genes,] %>% as.matrix() %>% t() %>% as.data.frame()
data <- scale(data,center = F) %>% as.data.frame()
data$celltype <- sce.epi.filter$cluster
data <- data[data$TP53>0,]

ggplot(data = data,aes(x = celltype,y = TP53,fill=celltype)) +
  geom_violin(linetype="blank",scale = "width",trim=F)+ylim(c(-1,4))

p <- DotPlot(sce.epi.filter,features = 'TP53')
VlnPlot(sce.epi.filter,features = 'TOP2A')
p




dat <- p$data
library(RColorBrewer)
display.brewer.all() #展示所有色板
color1 <- brewer.pal(9, "YlOrBr")[3:7]

p = ggplot(dat,aes(x=features.plot,y=id)) + 
  geom_point(aes(size=pct.exp,color= avg.exp.scaled))+
  scale_color_gradientn(colours = color1)+
  labs(
    size="Percent Expressed",
    fill = 'Average Expression'
  )+
  scale_size_continuous(range = c(6,12))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
pdf('./Fig/Epi_TP53_dotplot.pdf',width = 3,height = 4,onefile = FALSE)
print(p)
dev.off()
###参考文献：
###(KRT13) Kruppel like factor 4 promotes esophageal squamous cell carcinoma differentiation by upregulating keratin 13 expression
###(KRT16,S100A9，S100A8) Spatially and cell-type resolved quantitativeproteomic atlas of healthy human skin
###(MUC21,MUC4,CXCL17)
p <- DotPlot(sce.epi.filter,features = c('KRT16','KRT13','KRT6B','KRT6C','MUC21','MUC4','CXCL17','S100A9','S100A8'))
dat <- p$data

p = ggplot(dat,aes(x=features.plot,y=id)) + 
  geom_point(aes(size=pct.exp,color= avg.exp.scaled))+
  scale_color_gradientn(colours = color1)+
  labs(
    size="Percent Expressed",
    fill = 'Average Expression'
  )+
  scale_size_continuous(range = c(3,12))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
p
pdf('./Fig/Epi_signature_dotplot.pdf',width = 7,height = 3,onefile = FALSE)
print(p)
dev.off()


##halmarker score
BiocManager::install('cogena')
halmarker <- cogena::gmt2list('./h.all.v7.5.1.symbols.gmt')
DefaultAssay(sce.epi.filter)

lapply(halmarker,function(x){
  sce.epi.filter <- AddModuleScore(
    object = sce.epi.filter,
    features = x,
    ctrl = 5,
    name = names(x)
  )
})
for (i in seq(50)) {
  sce.epi.filter <- AddModuleScore(
    object = sce.epi.filter,
    features = halmarker[i],
    ctrl = 5,
    name = names(halmarker[i]))
}
x <- sce.epi.filter@meta.data
halmarker.score <- x[,14:ncol(x)]

library(stringr)

colnames(halmarker.score) <- stringr::str_sub(colnames(halmarker.score),end = -2)

halmarker.score$cluster <- sce.epi.filter$cluster
mean.score <- halmarker.score %>% dplyr::group_by(cluster) %>% 
  summarize(across(.cols = everything(),.fns = mean)) %>% column_to_rownames('cluster')

library(pheatmap)

pheatmap::pheatmap(t(mean.score),scale = 'none')


DimPlot(sce.epi.filter)



p <- VlnPlot(sce.epi.filter,features = 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION1')
p
dat <- p$data
p <- ggplot(dat, aes(x=ident, y=HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION1, fill=ident)) +
  geom_boxplot() +
  theme_classic()
pdf('./Fig/EMT.score.pdf',width = 5,height = 5,onefile = FALSE)
print(p)
dev.off()
####SCENIC
write.csv(t(as.matrix(sce.epi.filter@assays$RNA@counts)),file = "epi.counts.csv")


rm(list = ls())

###DEG in chemo vs pri
load('./sce.epi.filter.Rda')
DimPlot(sce.epi.filter)
sce.cl1 <- subset(sce.epi.filter,subset = cluster == 'cluster1' & sample != 'N1')
sce.cl1$group <- if_else(sce.cl1$sample%in%c('T1','T3'),'chemo','pri')
DimPlot(sce.cl1)
sce.cl1@active.ident <- as.factor(sce.cl1$group)
markers.cl1 <- FindAllMarkers(sce.cl1,logfc.threshold = 0.25,min.pct = 0.1,only.pos = T)
top100.cl1 <- markers.cl1 %>% group_by(cluster) %>% top_n(100,wt = avg_log2FC)

genes <- top100.cl1 %>% dplyr::filter(cluster == 'chemo') %>% pull(gene)
genes
library(clusterProfiler)
library(org.Hs.eg.db)
gene.df <- bitr(genes, fromType = "SYMBOL",
                toType ='ENTREZID',
                OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene          = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE)

chemo.up.cl1 <- as.data.frame(ego)
rownames(chemo.up.cl1) <- seq(nrow(chemo.up.cl1))
chemo.down.cl1 <- as.data.frame(ego)
rownames(chemo.down.cl1) <- seq(nrow(chemo.down.cl1))
chemo.up <- as.data.frame(ego)
rownames(chemo.up) <- seq(nrow(chemo.up))
chemo.down <- as.data.frame(ego)
rownames(chemo.down) <- seq(nrow(chemo.down))

chemo.up.cl1.go <- chemo.up.cl1[c(2,10,26),]
chemo.up.cl2.go <- chemo.up[c(1,4,13,86),]

chemo.down.cl1.go <- chemo.down.cl1[c(1,2,3,4),]

chemo.down.cl2.go <- chemo.down[c(1,2,3,4),]

####up genes
combine.up.go <- rbind(chemo.up.cl1.go,chemo.up.cl2.go)


chemo.up.cl1.go1 <- chemo.up.cl1[chemo.up.cl1$ID%in%combine.up.go$ID,]
chemo.up.cl1.go1 <- chemo.up.cl1.go1[match(combine.up.go$Description,chemo.up.cl1.go1$Description),]

chemo.up.cl2.go1 <- chemo.up[chemo.up$ID%in%combine.up.go$ID,]
chemo.up.cl2.go1 <- chemo.up.cl2.go1[match(combine.up.go$Description,chemo.up.cl2.go1$Description),]


combine.up.go <- rbind(chemo.up.cl1.go1,chemo.up.cl2.go1)
combine.up.go$cluster <- c(rep('cl1',7),rep('cl2',7))

combine.up.go$Description <- factor(combine.up.go$Description,levels = combine.up.go$Description[1:7])


####down genes
combine.down.go <- rbind(chemo.down.cl1.go,chemo.down.cl2.go)


chemo.down.cl1.go1 <- chemo.down.cl1[chemo.down.cl1$ID%in%combine.down.go$ID,]
chemo.down.cl1.go1 <- chemo.down.cl1.go1[match(combine.down.go$Description,chemo.down.cl1.go1$Description),]

chemo.down.cl2.go1 <- chemo.down[chemo.down$ID%in%combine.down.go$ID,]
chemo.down.cl2.go1 <- chemo.down.cl2.go1[match(combine.down.go$Description,chemo.down.cl2.go1$Description),]


combine.down.go <- rbind(chemo.down.cl1.go1,chemo.down.cl2.go1)
combine.down.go$cluster <- c(rep('cl1',8),rep('cl2',8))

combine.down.go$Description <- factor(combine.down.go$Description,levels = combine.down.go$Description[1:8])



###combine GO
combine.go <- rbind(combine.up.go,combine.down.go)
combine.go$group <- c(rep('up',nrow(combine.up.go)),rep('down',nrow(combine.down.go)))

combine.go$Description <- factor(combine.go$Description,levels =c(combine.up.go$Description[1:7],combine.down.go$Description[1:8]) )


combine.go$cus.p <- -1*log10(combine.go$p.adjust)
combine.go$cus.p[combine.go$cus.p>9] <- 9
p = ggplot(combine.go,aes(x=cluster,y=Description)) + 
  geom_point(aes(size=Count,color=cus.p))+
  scale_color_viridis(option = 'viridis')+
  scale_size_continuous(range = c(1,10))+
  labs(
    color=expression(-log[10](p.adjust)),
    size="Gene number",
    x="cluster"
  )+
  theme_classic()+
  geom_hline(yintercept = 7.5,linetype = "dashed", size = 1)+
  annotate('text',x=2,y=4,
           label=expression('up genes in chemo'),
           size=5,color='black')
p
pdf('./Fig/DEG.chemo_vs.pri.GO.dotplot.pdf',width = 7,height = 9,onefile = FALSE)
print(p)
dev.off()
save(combine.go,file = 'DEG.chemo_vs.pri.GO.Rda')

rm(list = ls())


####cell percent of MKI67 cells
data <- FetchData(sce.epi.filter,vars = c('MKI67','TOP2A','sample','cluster'))
data <- data %>% dplyr::filter(MKI67>0 |TOP2A >0 )
data$group <- if_else(data$sample%in%c('N1','N5','T4'),'adj','tumor')

data.count <-data %>%  group_by(group,cluster,sample) %>% summarise(counts = n())
data.count$counts <- data.count$counts/1269

##pie
data.pie <- as.data.frame( table(data$cluster))
label_value <- paste('(', round(data.pie$Freq/sum(data.pie$Freq) * 100, 1), '%)', sep = '')
label <- paste(data.pie$Var1, label_value, sep = '')
p <- ggplot(data = data.pie, mapping = aes(x = 'Counts', y = Freq, fill = Var1)) + 
  geom_bar(stat = 'identity', position = 'stack',width = 1,color = 'white')
p2 <- p+coord_polar(theta = 'y') + labs(x = '', y = '', title = '% Cycling Cells by Subpopulation')+
  scale_fill_manual(values = mycol)+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right",
        panel.grid=element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  geom_text(aes(y = data.pie$Freq/2 + c(0, cumsum(data.pie$Freq)[-length(data.pie$Freq)]), x = sum(data.pie$Freq)/1269, label = label))
##barplot
table(sce.epi.filter$cluster)
data.bar <- as.data.frame(table(sce.epi.filter$cluster))
data.bar$cycling <- data.pie$Freq
data.bar$cycling <-  data.bar$cycling / data.bar$Freq
data.bar$nocycling <- data.bar$Freq - data.bar$cycling
data.bar$nocycling <-  data.bar$nocycling /data.bar$Freq
data.bar <- data.bar %>% pivot_longer(cols = c(3,4),names_to = 'stra',values_to = 'percent')

c4a_gui()
mycol <- c4a('dark3',2)
library(ggpubr)
p1 = ggplot(data.bar, aes( x = Var1,y= percent,fill = stra))+
  geom_col(position = 'stack', width = 0.5,color = 'black')+
  scale_fill_manual(values = mycol)+
  theme_classic2()+
  scale_y_continuous(position = 'left',expand = c(0,0))+
  labs(x = '', y = '', title = '% Subpopulation in Cycling Cells')+
  theme(plot.title = element_text(hjust = 0.5))
p1



p1 <- ggplot(data.bar, aes( x = Var1,y=100 * percent,fill = stra,
                 stratum = stra, alluvium = stra))+
  scale_fill_manual(values = mycol)+
  geom_stratum(width = 0.7, color='black')+
  geom_alluvium(alpha = 0.3,
                width = 0.7,
                color='black',
                size = 0.3,
                curve_type = "linear")+
  scale_y_continuous(expand = c(0,0))+
  labs(x="Samples",y="Relative Abundance(%)",
       fill="group",
       title = '% Subpopulation in Cycling Cells')+
  theme_classic2()+
  theme(legend.position = 'right')



p3 <- p2+p1+patchwork::plot_layout(design ='AB' )

pdf('./Fig/cycling.cells.percent.pdf',width = 7,height = 4,onefile = FALSE)
print(p3)
dev.off()

#####
DefaultAssay(sce.epi.filter) <- 'RNA'
save(sce.epi.filter,file = 'sce.epi.filter.Rda') 
load('sce.epi.filter.Rda')
sce.epi.filter@active.ident <- as.factor(sce.epi.filter$cluster)
sce.epi.filter <- SCTransform(sce.epi.filter,variable.features.n = 10000)
epi.markers <- FindAllMarkers(sce.epi.filter,logfc.threshold = 0.4,min.pct = 0.1,only.pos = T)
top50.markers <- epi.markers %>% dplyr::group_by(cluster) %>% top_n(30,wt = avg_log2FC)
p <- DoHeatmap(sce.epi.filter,features = c(cl1.signatures,cl2.signatures))

cl1.signatures <- c('LCN2','SPRR2A','SPRR1A','SLPI','S100P','AGR2','WFDC2','SPRR2D','SPRR2E','S100A8','S100A9','S100A7')
cl2.signatures <- c('MT2A','COL17A1','PTHLH','CAV1','DST','LGALS1','STMN1','TGFBI','TP63','TUBA1B','FLNA','KRT15','MT1E')
signature <- c(cl1.signatures,cl2.signatures)

data <- sce.epi.filter@assays$RNA@data[signature,] %>% as.matrix() %>% as.data.frame()
data <- t(data) %>% as.data.frame()
data$cluster <- sce.epi.filter$cluster
data <- data %>% dplyr::arrange(cluster)


data.cl1 <- data %>% dplyr::filter(cluster == 'cl1')
data.cl1 <- data.cl1[sample(nrow(data.cl1),replace = F),]
data.cl2 <- data %>% dplyr::filter(cluster == 'cl2')
data.cl2 <- data.cl2[sample(nrow(data.cl2),replace = F),]


data <- rbind(data.cl1,data.cl2)



cluster <- data %>% dplyr::select('cluster')
data <- data[,1:ncol(data)-1]
data <- t(data) %>% as.data.frame()

library(pheatmap)
bk = c(seq(-2,2,by=0.01))
cluster$cluster <- as.factor(cluster$cluster)
ann_colors = list(
  cluster = c(cl1 = mycol[1],cl2 = mycol[2]))



p <- pheatmap(data,border = F,
         cluster_row = FALSE,
         cluster_cols = F,
         show_rownames=T,show_colnames=F,
         scale = "row",
         annotation_col = cluster,
         annotation_colors = ann_colors,
         color = colorRampPalette(colors = c("#1BB7E5","#000000","#FFFF00"))(length(bk)),
         breaks = bk)
pdf('./Fig/signature.heatmap.pdf',width = 5,height = 5,onefile = FALSE)
print(p)
dev.off()
library(SeuratDisk)
?SeuratDisk

?SeuratDisk
BiocManager::install('GSVA')
