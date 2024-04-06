###
rm(list = ls())
Sys.setenv(LANGUAGE = "en")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
library(tidyverse)
library(Seurat)
library(stringr)
library(ggpubr)
library(ggthemes)
library(data.table)
library(rio)
load('../Fig7_analysis/sce.fibro.SCT.Rda')
expression <- sce.fibro.SCT@assays$RNA@data %>% as.matrix()  %>% t() %>% as.data.frame()

metadata <- sce.fibro.SCT@meta.data
#####L_R预处理########################
pdata <- read.table('D:/jupyter-notebook/ESCC/cpdb/fibro/results/statistical_analysis_pvalues_method2.txt',sep = '\t',header = T)
colnames(pdata)[12:length(colnames(pdata))] <- colnames(pdata)[12:length(colnames(pdata))]  %>% str_replace('\\.','|')
decov <- read.table('D:/jupyter-notebook/ESCC/cpdb/fibro/results/statistical_analysis_deconvoluted_method2.txt',sep = '\t',header = T)
#####step1 填补LR_sig_long和RL_sig_long中基因名
combine_ID1 <- decov %>% dplyr::select('gene_name','uniprot')
colnames(combine_ID1) <- c('genes','id')
combine_ID2 <- decov %>% dplyr::select('gene_name','complex_name')
colnames(combine_ID2) <- c('genes','id')
combine_ID <- rbind(combine_ID1,combine_ID2)

partner_a_id <- pdata %>% dplyr::select('partner_a') %>% pull(partner_a) %>% str_split(':',2) %>% 
  lapply(function(x){x[2]}) %>% unlist()
pdata$gene_a <-combine_ID %>% dplyr::slice(match(partner_a_id,combine_ID$id)) %>% pull(genes)
partner_b_id <- pdata %>% dplyr::select('partner_b') %>% pull(partner_b) %>% str_split(':',2) %>% 
  lapply(function(x){x[2]}) %>% unlist()
pdata$gene_b <-combine_ID %>% dplyr::slice(match(partner_b_id,combine_ID$id)) %>% pull(genes)

####整合配体受体顺序
#condition 1, receptor a is true and receptor b is false
RL_pdata <- pdata %>% filter(receptor_a == 'True' & receptor_b == 'False') %>% 
  dplyr::select(!matches('cl.*cl')) %>% 
  dplyr::select(!matches('fibro.*fibro'))

RL_sig_long <- RL_pdata %>% pivot_longer(12:length(colnames(RL_pdata)), 
                                         names_to = "group", 
                                         values_to = "p_value") %>% 
  dplyr::filter(p_value<0.05) %>% tidyr::separate(col = group, into = c('R_cells','L_cells'),sep = "\\|") %>% 
  dplyr::select(gene_b,gene_a,L_cells,R_cells,p_value)
colnames(RL_sig_long) <- c('ligands','receptors','L_cells','R_cells','p_value')

#condition 2, receptor b is true and receptor a is false
LR_pdata <- pdata %>% filter(receptor_a == 'False' & receptor_b == 'True') %>% 
  dplyr::select(!matches('cl.*cl')) %>% 
  dplyr::select(!matches('fibro.*fibro'))
LR_sig_long <- LR_pdata %>% pivot_longer(12:length(colnames(LR_pdata)), 
                                         names_to = "group", 
                                         values_to = "p_value") %>% 
  filter(p_value<0.05) %>% separate(col = group, into = c('L_cells','R_cells'),sep = "\\|") %>% 
  dplyr::select(gene_a,gene_b,L_cells,R_cells,p_value)
colnames(LR_sig_long) <- c('ligands','receptors','L_cells','R_cells','p_value')
interaction_sig_long <- rbind(RL_sig_long,LR_sig_long) %>% unite(group,L_cells,R_cells, sep = "|",remove = FALSE)
####保存整合得interaction_sig_long
save(interaction_sig_long,file = 'Epi_fibro_interaction_sig_long.Rda')
load('Epi_fibro_interaction_sig_long.Rda')
#####受体配体计数counts################################
sig_counts <- interaction_sig_long %>% group_by(group) %>%
  summarise(counts = n()) %>%
  separate(col = group,into = c("cell_L", "cell_R"), sep = "\\|")


##Epi is sender cells
sender_epi_counts <- sig_counts %>% dplyr::slice(grep('cl',sig_counts$cell_L))
p1 = ggplot(sender_epi_counts, aes(x=cell_R, y = counts,fill = cell_L)) +
  geom_bar(width = .7,stat = 'identity',position = 'dodge') +
  scale_fill_excel_new()+scale_y_continuous(limits = c(0,40),expand =c(0,0))+
  theme_few()+theme(axis.text.x = element_text(angle = 45,hjust = 1))
p1
pdf('./fig8/sender_epi_LR_counts.pdf',width = 6,height = 4,onefile = FALSE)
print(p1)
dev.off()
##fibro is sender cells
sender_fibro_counts <- sig_counts %>% dplyr::slice(grep('fibro',sig_counts$cell_L))
p2 = ggplot(sender_fibro_counts, aes(x=cell_L, y = counts,fill = cell_R)) +
  geom_bar(width = .7,stat = 'identity',position = 'dodge') +
  scale_fill_excel_new()+scale_y_continuous(limits = c(0,30),expand =c(0,0))+
  theme_test()+theme(axis.text.x = element_text(angle = 45,hjust = 1))
p2
pdf('./fig8/sender_fibro_LR_counts.pdf',width = 6,height = 4,onefile = FALSE)
print(p2)
dev.off()




#######-------------------------------------------------------------###
#########L-R pairs heatmap landscape#########

##Epi1 is sender cells and fibro1 is receiver
epi_sender_sig_1 <- interaction_sig_long %>% filter(R_cells%in%c('CFD_fibro') & L_cells%in%c('cl2')) %>% 
  unite(interaction,ligands,receptors,remove = FALSE) %>% distinct(interaction,.keep_all = T)
epi_sender_sig_2 <- interaction_sig_long %>% filter(R_cells%in%c('POSTN_fibro') & L_cells%in%c('cl2')) %>% 
  unite(interaction,ligands,receptors,remove = FALSE) %>% distinct(interaction,.keep_all = T)
epi_sender_sig <- rbind(epi_sender_sig_1,epi_sender_sig_2) %>% distinct(interaction,.keep_all = T)
  
  
sce.epi.filter@active.ident <- as.factor(sce.epi.filter$cluster)
epi.markers <- FindAllMarkers(sce.epi.filter,logfc.threshold = 0.3,min.pct = 0.1,only.pos = T)
cl2.markers <- epi.markers %>% dplyr::filter(cluster == 'cl2',avg_log2FC >0.3)
epi_sender_sig <- epi_sender_sig %>% dplyr::filter(ligands%in%cl2.markers$gene)

FeaturePlot(sce.epi.filter,features = epi_sender_sig$ligands)

FeaturePlot(sce.fibro.SCT,features = 'CD46',max.cutoff = 2)


FeaturePlot(sce.epi.filter,features = 'TGFBR1')
FeaturePlot(sce.fibro.SCT,features = 'CD74')
table(sce.epi.filter$group)
test <- subset(sce.epi.filter,subset = group == 'tumor')

DefaultAssay(sce.epi.filter)
epi_sender_sig$p_value <-  if_else(epi_sender_sig$p_value<0.05 & epi_sender_sig$p_value >= 0.01,'1',
                                   if_else(epi_sender_sig$p_value >0.001 & epi_sender_sig$p_value <= 0.01,'2','3'))

library(RColorBrewer)
#display.brewer.all() #展示所有色板
color1 <- brewer.pal(9, "RdBu") #取颜色
library(scales)
show_col(color1) #可视化
show_col(color2) #可视化
brewer.pal.info #显示所有颜色信息
#constrct long p matrix
epi_sender_sig_weight <- epi_sender_sig %>% 
  select(ligands,receptors,p_value) %>% 
  pivot_wider(names_from = receptors, values_from = p_value ) %>% 
  column_to_rownames('ligands')

epi_sender_sig_long <- epi_sender_sig_weight %>% rownames_to_column('ligands') %>% 
  pivot_longer(2:length(colnames(.)),names_to = 'receptors',values_to = 'value')
epi_sender_sig_long[is.na(epi_sender_sig_long)] <- '0'
epi_sender_sig_long$value <- epi_sender_sig_long$value %>% as.character()
color = c(brewer.pal(9, "Greys")[1],brewer.pal(9, "YlOrRd")[c(7,9)])
###plot p_value heatmap
epi_fibro_p_vis <- epi_sender_sig_long %>% 
  ggplot(aes(receptors,ligands,fill=value))+
  theme_minimal()+
  geom_tile(color="black",size=0.5)+scale_fill_manual(values = color)+
  xlab('')+
  ylab('')+
  theme(
    axis.text.x.bottom = element_blank(),#x轴刻度文字大小
    axis.text.y.left = element_blank(), #y轴刻度文字大小
    axis.ticks = element_blank(), #不显示坐标轴刻度
    panel.border = element_blank()
  )

epi_fibro_p_vis

###配体表达量热
metadata <- read.table('./metadata.txt',header = T)
load('../Fig2_analysis/sce.epi.filter.Rda')
epi.exp <- sce.epi.filter@assays$RNA@data %>% as.matrix() %>% t() %>% as.data.frame()
fibro.exp <- sce.fibro.SCT@assays$RNA@data %>% as.matrix() %>% t() %>% as.data.frame()
ligand_exp <- epi.exp[,levels(factor(epi_sender_sig_long$ligands))]

ligand_mean_exp <- ligand_exp %>% rownames_to_column('Cell') %>% 
  inner_join(metadata,by = 'Cell') %>% group_by(cell_type) %>%
  select(-Cell) %>% summarise_each(funs = mean) %>% slice(grep('cl',.$cell_type)) %>% 
  column_to_rownames('cell_type')
ligand_mean_exp <- ligand_mean_exp %>% t() %>% as.data.frame() %>% arrange(desc(cl2)) %>% t() %>% as.data.frame()
order.ligands <- colnames(ligand_mean_exp)
ligand_mean_exp <- ligand_mean_exp %>%
  as.data.frame() %>% rownames_to_column('celltype') %>% 
  pivot_longer(2:length(colnames(.)),names_to = 'ligands',values_to = 'value')
ligand_mean_exp$ligands <- factor(ligand_mean_exp$ligands,levels = rev(order.ligands))
ligand_mean_exp$value[ligand_mean_exp$value>0.8] <- 0.8
display.brewer.all()
color1 <- brewer.pal(11, "RdBu")[3:9] %>% rev()##颜色倒序排列
ligand_epi_exp <- ligand_mean_exp %>% 
  ggplot(aes(celltype,ligands,fill=value))+
  geom_tile(color="black",size=0.5)+scale_fill_gradientn(colours = color1)+
  theme_test()+
  xlab('Sender')+
  theme(
    axis.text.x.bottom = element_text(size=4,angle = 45,
                                      hjust = 1,vjust= 1),#x轴刻度文字大小
    axis.text.y.left = element_text(size = 4,vjust= 1,hjust = 1), #y轴刻度文字大小
    axis.ticks = element_blank(), 
    panel.border = element_blank()
  )
ligand_epi_exp
###受体表达量热图
receptor_exp <- fibro.exp[,unique(epi_sender_sig_long$receptors)]
receptor_mean_exp <- receptor_exp %>% rownames_to_column('Cell') %>% 
  inner_join(metadata,by = 'Cell') %>% group_by(cell_type) %>%
  select(-Cell) %>% summarise_each(funs = mean) %>% slice(grep('fibro',.$cell_type)) %>% 
  column_to_rownames('cell_type')%>% scale() %>% as.data.frame()

order.receptor <- receptor_mean_exp %>% t() %>% as.data.frame() %>% arrange(desc(POSTN_fibro)) %>% 
  rownames_to_column('id') %>% pull('id')
order.receptor
test <- t(receptor_mean_exp) %>% as.data.frame()
receptor_mean_exp <- receptor_mean_exp  %>% rownames_to_column('celltype') %>% 
  pivot_longer(2:length(colnames(.)),names_to = 'receptors',values_to = 'value')
receptor_mean_exp$receptors <- factor(receptor_mean_exp$receptors,levels = order.receptor )

receptor_fibro_vis <-  receptor_mean_exp %>% 
  ggplot(aes(receptors,celltype,fill=value))+
  geom_tile(color="black",size=0.5)+scale_fill_gradientn(colours = color1)+
  theme_test()+
  ylab('receiver')+
  theme(
    axis.text.x.bottom = element_text(size=4,angle = 45,
                                      hjust = 1,vjust= 1.3),#x轴刻度文字大小
    axis.text.y.left = element_text(size = 4,vjust= 1,hjust = 1), #y轴刻度文字大小
    axis.ticks = element_blank(), #不显示坐标轴刻度
    panel.border = element_blank()
  )
receptor_fibro_vis
###plot order p_value heatmap
epi_sender_sig_long$ligands <- factor(epi_sender_sig_long$ligands,levels = rev(order.ligands))
epi_sender_sig_long$receptors <- factor(epi_sender_sig_long$receptors,levels = order.receptor)
epi_fibro_p_vis <- epi_sender_sig_long %>% 
  ggplot(aes(receptors,ligands,fill=value))+
  theme_minimal()+
  geom_tile(color="black",size=0.5)+scale_fill_manual(values = color)+
  xlab('')+
  ylab('')+
  theme(
    axis.text.x.bottom = element_blank(),#x轴刻度文字大小
    axis.text.y.left = element_blank(), #y轴刻度文字大小
    axis.ticks = element_blank(), #不显示坐标轴刻度
    panel.border = element_blank()
  )

epi_fibro_p_vis
####combineplot
library(cowplot)
combine_plot <- plot_grid(ligand_epi_exp+theme(legend.position = 'none'),
                          epi_fibro_p_vis+theme(legend.position = 'none'),
                          NULL,
                          receptor_fibro_vis+theme(legend.position = 'none'),
                          align = "hv",
                          nrow = 2,
                          rel_widths =c(2,11),
                          rel_heights = c(10,3))

combine_plot

legends = plot_grid(
  as_ggplot(get_legend(receptor_fibro_vis)),
  as_ggplot(get_legend(epi_fibro_p_vis)),
  ncol = 1,
  align = "hv")
dev.off()
pdf('./fig8/epi_sender_LR_landscape1.pdf',width = 8,height = 6,onefile = FALSE)

plot_grid(combine_plot, 
          legends, 
          rel_widths =c(20,1), nrow = 1)
dev.off()


##Fibro is sender cells and epi is receiver
fibro_sender_sig <- interaction_sig_long %>% filter(R_cells%in%c('cl2') & L_cells%in%c('CFD_fibro','POSTN_fibro')) %>% 
  unite(interaction,ligands,receptors,remove = FALSE) %>% distinct(interaction,.keep_all = T)

###配体表达量热图
ligand_exp <- fibro.exp[,unique(fibro_sender_sig$ligands)]

ligand_mean_exp <- ligand_exp %>% rownames_to_column('Cell') %>% 
  inner_join(metadata,by = 'Cell') %>% group_by(cell_type) %>%
  select(-Cell) %>% summarise_each(funs = mean) %>% slice(grep('fibro',.$cell_type)) %>% 
  column_to_rownames('cell_type')

ligand_mean_exp <- ligand_mean_exp %>% scale() %>% t() %>% as.data.frame() %>% arrange(desc(POSTN_fibro)) %>% 
  t() %>% as.data.frame()
order.ligands <- colnames(ligand_mean_exp)
ligand_mean_exp <- ligand_mean_exp %>% rownames_to_column('celltype') %>% 
  pivot_longer(2:length(colnames(.)),names_to = 'ligands',values_to = 'value')

ligand_mean_exp$ligands <- factor(ligand_mean_exp$ligands,levels = rev(order.ligands))
color1 <- rev(brewer.pal(9, "RdBu"))##颜色倒序排列
ligand_fibro_vis <- ligand_mean_exp %>% 
  ggplot(aes(celltype,ligands,fill=value))+
  geom_tile(color="black",size=0.5)+scale_fill_gradientn(colours = color1)+
  theme_test()+
  xlab('Sender')+
  theme(
    axis.text.x.bottom = element_text(size=4,angle = 45,
                                      hjust = 1,vjust= 1),#x轴刻度文字大小
    axis.text.y.left = element_text(size = 4,vjust= 1,hjust = 1), #y轴刻度文字大小
    axis.ticks = element_blank(), 
    panel.border = element_blank()
  )
ligand_fibro_vis

###受体表达量热图
receptor_exp <- epi.exp[,unique(fibro_sender_sig$receptors)]
receptor_mean_exp <- receptor_exp %>% rownames_to_column('Cell') %>% 
  inner_join(metadata,by = 'Cell') %>% group_by(cell_type) %>%
  select(-Cell) %>% summarise_each(funs = mean) %>% slice(grep('cl',.$cell_type)) %>% 
  column_to_rownames('cell_type')
receptor_mean_exp <- receptor_mean_exp %>% t() %>% as.data.frame() %>% arrange(desc(cl2)) %>% 
  t() %>% as.data.frame()
order.receptor <- colnames(receptor_mean_exp)
receptor_mean_exp <- receptor_mean_exp  %>% as.data.frame() %>% rownames_to_column('celltype') %>% 
  pivot_longer(2:length(colnames(.)),names_to = 'receptors',values_to = 'value')
receptor_mean_exp$value[receptor_mean_exp$value>0.4] <- 0.4
receptor_mean_exp$receptors <- factor(receptor_mean_exp$receptors,levels = order.receptor)
receptor_epi_vis <-  receptor_mean_exp %>% 
  ggplot(aes(receptors,celltype,fill=value))+
  geom_tile(color="black",size=0.5)+scale_fill_gradientn(colours = color1)+
  theme_test()+
  ylab('receiver')+
  theme(
    axis.text.x.bottom = element_text(size=8,angle = 45,
                                      hjust = 1,vjust=1.2),#x轴刻度文字大小
    axis.text.y.left = element_text(size = 8,vjust= 1,hjust = 1), #y轴刻度文字大小
    axis.ticks = element_blank(), #不显示坐标轴刻度
    panel.border = element_blank()
  )
receptor_epi_vis
FeaturePlot(sce.epi.filter,features = 'CD46')

#constrct long p matrix
fibro_sender_sig$p_value <-  if_else(fibro_sender_sig$p_value<0.05 & fibro_sender_sig$p_value >= 0.01,'1',
                                   if_else(fibro_sender_sig$p_value >0.001 & fibro_sender_sig$p_value <= 0.01,'2','3'))
fibro_sender_sig_weight <- fibro_sender_sig %>% 
  select(ligands,receptors,p_value) %>% 
  pivot_wider(names_from = receptors, values_from = p_value ) %>% 
  column_to_rownames('ligands')

fibro_sender_sig_long <- fibro_sender_sig_weight %>% rownames_to_column('ligands') %>% 
  pivot_longer(2:length(colnames(.)),names_to = 'receptors',values_to = 'value')
fibro_sender_sig_long[is.na(fibro_sender_sig_long)] <- '0'
fibro_sender_sig_long$value <- fibro_sender_sig_long$value %>% as.character()
color = c(brewer.pal(9, "Greys")[1],brewer.pal(9, "YlOrRd")[c(7,9)])
fibro_sender_sig_long$ligands <- factor(fibro_sender_sig_long$ligands,levels = order.ligands)
fibro_sender_sig_long$receptors <- factor(fibro_sender_sig_long$receptors,levels = order.receptor)
###plot p_value heatmap
fibro_epi_p_vis <- fibro_sender_sig_long %>% 
  ggplot(aes(receptors,ligands,fill=value))+
  theme_minimal()+
  geom_tile(color="black",size=0.5)+scale_fill_manual(values = color)+
  xlab('')+
  ylab('')+
  theme(
    axis.text.x.bottom = element_blank(),#x轴刻度文字大小
    axis.text.y.left = element_blank(), #y轴刻度文字大小
    axis.ticks = element_blank(), #不显示坐标轴刻度
    panel.border = element_blank()
  )

fibro_epi_p_vis


####combineplot

library(cowplot)
combine_plot1 <- plot_grid(ligand_fibro_vis+theme(legend.position = 'none'),
                          fibro_epi_p_vis+theme(legend.position = 'none'),
                          NULL,
                          receptor_epi_vis+theme(legend.position = 'none'),
                          align = "hv",
                          nrow = 2,
                          rel_widths =c(4,11),
                          rel_heights = c(10,3))
combine_plot1

legends1 = plot_grid(
  as_ggplot(get_legend(receptor_epi_vis)),
  as_ggplot(get_legend(fibro_epi_p_vis)),
  ncol = 1,
  align = "hv")
dev.off()


pdf('./fig8/fibro_sender_LR_landscape1.pdf',width = 10,height = 8,onefile = FALSE)
plot_grid(combine_plot1, 
          legends1, 
          rel_widths =c(20,1), nrow = 1)
dev.off()


####key interaction between CL2 and CFD fibro
####JAG1-NOTCH2
library(patchwork)
p <- FeaturePlot(sce.epi.filter,features = 'JAG1',max.cutoff = 2)
p

p1 <- FeaturePlot(sce.fibro.SCT,features = 'NOTCH2',max.cutoff = 1)
p2 <- p+p1+patchwork::plot_layout(design = 'AB')
pdf('./fig8/JAG1-NOTCH2_featurePlot.pdf',width = 10,height = 5,onefile = FALSE)
print(p2)
dev.off()













