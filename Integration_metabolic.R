rm(list = ls())
library(dplyr)
library(ggplot2)
library(Seurat)
library("rlist")
###坐标转化
file <- list.files('F:\\biodata\\MWY-23-1620-a_BIOINF-2327608\\20231025_BIOINF-2327608.data_result\\')
path <- "F:\\biodata\\ESCC\\ESCC_代谢\\MWY-23-1620-a_BIOINF-2327608\\20231025_BIOINF-2327608.data_result\\T1-2_coordinates_intensities.txt"

quantification <- read.table(file = path,header = T)
trans_path <- 'D:/jupyter-notebook/image_pre/transform_coodination_T1.txt'
trans_coordinate <- read.table(file = trans_path,header = T)

trans_coordinate$value <- rep(1,17952)

##读取spatail图片像素
coord_ST <- read.table("D:/localdata/Desktop/ESCC/ourdata/spatial/T1-2/spatial/tissue_positions.csv",sep = ',',header = T)
coord_ST$value <- rep(1,4992)
coord_ST <- coord_ST %>% filter(in_tissue == 1)
coord_ST$pxl_row_in_fullres <- coord_ST$pxl_row_in_fullres * 0.2
coord_ST$pxl_col_in_fullres <- coord_ST$pxl_col_in_fullres * 0.2


###validation coordinate
# p <- ggplot(trans_coordinate, aes(x, y)) +
#   geom_point(aes(fill = value),colour = "black") + 
#   scale_fill_gradient(low = "black",high = "black")+
#   theme_grey(base_size = 0) + 
#   labs(x = NULL,y = NULL) +
#   scale_y_discrete(expand = c(0, 0)) +
#   scale_y_discrete(expand = c(0, 0)) + 
#   theme(legend.position = "none",
#         axis.ticks = element_blank(), 
#         axis.text.x = element_blank(),
#         axis.text = element_blank(),
#         panel.background = element_blank())
# p
# 
# coo <- read.table("D:/localdata/Desktop/ESCC/ourdata/spatial/T1-2/spatial/tissue_positions.csv",sep = ',',header = T)
# coo$value <- rep(1,4992)
# coo <- coo %>% filter(in_tissue == 1)
# coo$pxl_row_in_fullres <- coo$pxl_row_in_fullres * 0.2
# coo$pxl_col_in_fullres <- coo$pxl_col_in_fullres * 0.2
# 
# 
# p1 <- ggplot(coo, aes(pxl_col_in_fullres, pxl_row_in_fullres)) +
#   geom_point(aes(fill = value),colour = "black") + 
#   scale_fill_gradient(low = "black",high = "black")+
#   theme_grey(base_size = 10) + 
#   labs(x = NULL,y = NULL) +
#   scale_x_discrete(expand = c(0, 0)) +
#   scale_y_discrete(expand = c(0, 0)) + 
#   theme(legend.position = "none",
#         axis.ticks = element_blank(), 
#         axis.text.x = element_text(size = 8,angle=30,hjust = 1),
#         axis.text = element_text(color='black'),
#         panel.background = element_blank())
# p1
# p


####caculate Metabolism intensity for a give spatial spot
##通过观察spot像素分布规律，同一行，2个spot之间的像素距离为7个像素。由先验知识可知，spot之间的距离为100um，
###则每个像素之间的距离约为14um，则每个spot的距离转换为像素距离约4个像素。
###像素距离过少会影响数据匹配结果，则适当放开匹配的像素距离。
###注意，每个空间HE切片的像素距离需要独立观察计算，因为该像素距离不具备普适性。
##T1
rownames(coord_ST) <- seq(3653)
index_list <- lapply(1:nrow(coord_ST),function(i){
  xi = coord_ST[i,6]
  yi = coord_ST[i,5]
  index = c()
  for (n in 1:nrow(trans_coordinate)) {
    xn = trans_coordinate[n,1]
    yn = trans_coordinate[n,2]
    if ((xn<xi+3.7) & (xn > xi-3.7) & (yn<yi+3.7) & (yn > yi-3.7)) {
      index = append(index,n)
    }else{next}
  }
  if (is.null(index)) {
    index = NA
  }else{index = index}
  return(index)
})

index_list <- test
rm(test)



##计算给定spot平均代谢强度
  avg_intensity_list <- lapply(index_list, function(x){
  if(length(x)<2 & sum(is.na(x))){
      avg_intensity = rep(0,1474)
    }else{
      avg_intensity <- quantification[x,5:ncol(quantification)]%>% t() %>% as.numeric()
    }
  if (length(x)>1) {
    dat <- quantification[x,5:ncol(quantification)]
    avg_intensity = apply(dat,2,mean) %>% as.numeric()
  }
    return(avg_intensity)
})

avg_intensity_df <- list.rbind(avg_intensity_list) %>% as.data.frame()
rownames(avg_intensity_df) <- coord_ST$barcode

save(avg_intensity_df,file = './transform_data/Spot_intensity_transform.Rda')



#T2
###坐标转化
path <- "F:\\biodata\\MWY-23-1620-a_BIOINF-2327608\\20231025_BIOINF-2327608.data_result\\T2-1_coordinates_intensities.txt"
quantification <- read.table(file = path,header = T)
trans_path <- 'D:/jupyter-notebook/image_pre/transform_coodination_T2.txt'
trans_coordinate <- read.table(file = trans_path,header = T)
##读取spatail图片像素
coord_ST <- read.table("D:/localdata/Desktop/ESCC/ourdata/spatial/T2-1/spatial/tissue_positions.csv",sep = ',',header = T)
coord_ST <- coord_ST %>% filter(in_tissue == 1)
##转换为600X600像素（根据json文件确定）
coord_ST$pxl_row_in_fullres <- coord_ST$pxl_row_in_fullres * 0.2
coord_ST$pxl_col_in_fullres <- coord_ST$pxl_col_in_fullres * 0.2
###对于给定的Spot，计算该spot所映射到的代谢空间坐标index 
rownames(coord_ST) <- seq(3790)

index_list <- lapply(1:nrow(coord_ST),function(i){
  xi = coord_ST[i,6]
  yi = coord_ST[i,5]
  index = c()
  for (n in 1:nrow(trans_coordinate)) {
    xn = trans_coordinate[n,1]
    yn = trans_coordinate[n,2]
    if ((xn<xi+3.7) & (xn > xi-3.7) & (yn<yi+3.7) & (yn > yi-3.7)) {
      index = append(index,n)
    }else{next}
  }
  if (is.null(index)) {
    index = NA
  }else{index = index}
  return(index)
})

##计算给定spot平均代谢强度
avg_intensity_list <- lapply(index_list, function(x){
  if(length(x)<2 & sum(is.na(x))){
    avg_intensity = rep(0,1474)
  }else{
    avg_intensity <- quantification[x,5:ncol(quantification)]%>% t() %>% as.numeric()
  }
  if (length(x)>1) {
    dat <- quantification[x,5:ncol(quantification)]
    avg_intensity = apply(dat,2,mean) %>% as.numeric()
  }
  return(avg_intensity)
})

avg_intensity_df <- list.rbind(avg_intensity_list) %>% as.data.frame()
rownames(avg_intensity_df) <- coord_ST$barcode

save(avg_intensity_df,file = './transform_data/Spot_intensity_transform_T2.Rda')



#T4
###坐标转化
path <- "F:\\biodata\\MWY-23-1620-a_BIOINF-2327608\\20231025_BIOINF-2327608.data_result\\N4-2_coordinates_intensities.txt"
quantification <- read.table(file = path,header = T)
trans_path <- 'D:/jupyter-notebook/image_pre/transform_coodination_T4.txt'
trans_coordinate <- read.table(file = trans_path,header = T)
##读取spatail图片像素
coord_ST <- read.table("D:/localdata/Desktop/ESCC/ourdata/spatial/N4-2/spatial/tissue_positions.csv",sep = ',',header = T)
coord_ST <- coord_ST %>% filter(in_tissue == 1)
##转换为600X600像素（根据json文件确定）
coord_ST$pxl_row_in_fullres <- coord_ST$pxl_row_in_fullres * 0.2
coord_ST$pxl_col_in_fullres <- coord_ST$pxl_col_in_fullres * 0.2
###对于给定的Spot，计算该spot所映射到的代谢空间坐标index 
rownames(coord_ST) <- seq(1710)

index_list <- lapply(1:nrow(coord_ST),function(i){
  xi = coord_ST[i,6]
  yi = coord_ST[i,5]
  index = c()
  for (n in 1:nrow(trans_coordinate)) {
    xn = trans_coordinate[n,1]
    yn = trans_coordinate[n,2]
    if ((xn<xi+3.7) & (xn > xi-3.7) & (yn<yi+3.7) & (yn > yi-3.7)) {
      index = append(index,n)
    }else{next}
  }
  if (is.null(index)) {
    index = NA
  }else{index = index}
  return(index)
})

##计算给定spot平均代谢强度
avg_intensity_list <- lapply(index_list, function(x){
  if(length(x)<2 & sum(is.na(x))){
    avg_intensity = rep(0,1474)
  }else{
    avg_intensity <- quantification[x,5:ncol(quantification)]%>% t() %>% as.numeric()
  }
  if (length(x)>1) {
    dat <- quantification[x,5:ncol(quantification)]
    avg_intensity = apply(dat,2,mean) %>% as.numeric()
  }
  return(avg_intensity)
})

avg_intensity_df <- list.rbind(avg_intensity_list) %>% as.data.frame()
rownames(avg_intensity_df) <- coord_ST$barcode

save(avg_intensity_df,file = './transform_data/Spot_intensity_transform_T4.Rda')



#T3
###坐标转化
path <- "F:\\biodata\\MWY-23-1620-a_BIOINF-2327608\\20231025_BIOINF-2327608.data_result\\T3_coordinates_intensities.txt"
quantification <- read.table(file = path,header = T)
trans_path <- 'D:/jupyter-notebook/image_pre/transform_coodination_T3.txt'
trans_coordinate <- read.table(file = trans_path,header = T)
##读取spatail图片像素
coord_ST <- read.table("F:/biodata/ESCC/spatial_preprocess/ESCC_T3/outs/spatial/tissue_positions.csv",sep = ',',header = T)
coord_ST <- coord_ST %>% filter(in_tissue == 1)
##转换为600X600像素（根据json文件确定）
coord_ST$pxl_row_in_fullres <- coord_ST$pxl_row_in_fullres * 0.2
coord_ST$pxl_col_in_fullres <- coord_ST$pxl_col_in_fullres * 0.2
###对于给定的Spot，计算该spot所映射到的代谢空间坐标index 
rownames(coord_ST) <- seq(2856)

index_list <- lapply(1:nrow(coord_ST),function(i){
  xi = coord_ST[i,6]
  yi = coord_ST[i,5]
  index = c()
  for (n in 1:nrow(trans_coordinate)) {
    xn = trans_coordinate[n,1]
    yn = trans_coordinate[n,2]
    if ((xn<xi+3.7) & (xn > xi-3.7) & (yn<yi+3.7) & (yn > yi-3.7)) {
      index = append(index,n)
    }else{next}
  }
  if (is.null(index)) {
    index = NA
  }else{index = index}
  return(index)
})

##计算给定spot平均代谢强度
avg_intensity_list <- lapply(index_list, function(x){
  if(length(x)<2 & sum(is.na(x))){
    avg_intensity = rep(0,1474)
  }else{
    avg_intensity <- quantification[x,5:ncol(quantification)]%>% t() %>% as.numeric()
  }
  if (length(x)>1) {
    dat <- quantification[x,5:ncol(quantification)]
    avg_intensity = apply(dat,2,mean) %>% as.numeric()
  }
  return(avg_intensity)
})

avg_intensity_df <- list.rbind(avg_intensity_list) %>% as.data.frame()
rownames(avg_intensity_df) <- coord_ST$barcode

save(avg_intensity_df,file = './transform_data/Spot_intensity_transform_T3.Rda')



#T5
###坐标转化
path <- "F:\\biodata\\MWY-23-1620-a_BIOINF-2327608\\20231025_BIOINF-2327608.data_result\\Zhao-Ca1_coordinates_intensities.txt"
quantification <- read.table(file = path,header = T)
trans_path <- 'D:/jupyter-notebook/image_pre/transform_coodination_T5.txt'
trans_coordinate <- read.table(file = trans_path,header = T)
##读取spatail图片像素
coord_ST <- read.table("F:/biodata/ESCC/spatial_preprocess/ESCC-zhao-Ca1/outs/spatial/tissue_positions.csv",sep = ',',header = T)
coord_ST <- coord_ST %>% filter(in_tissue == 1)
##转换为600X600像素（根据json文件确定）
coord_ST$pxl_row_in_fullres <- coord_ST$pxl_row_in_fullres * 0.19
coord_ST$pxl_col_in_fullres <- coord_ST$pxl_col_in_fullres * 0.19
###对于给定的Spot，计算该spot所映射到的代谢空间坐标index 
rownames(coord_ST) <- seq(2169)

index_list <- lapply(1:nrow(coord_ST),function(i){
  xi = coord_ST[i,6]
  yi = coord_ST[i,5]
  index = c()
  for (n in 1:nrow(trans_coordinate)) {
    xn = trans_coordinate[n,1]
    yn = trans_coordinate[n,2]
    if ((xn<xi+3.7) & (xn > xi-3.7) & (yn<yi+3.7) & (yn > yi-3.7)) {
      index = append(index,n)
    }else{next}
  }
  if (is.null(index)) {
    index = NA
  }else{index = index}
  return(index)
})

##计算给定spot平均代谢强度
avg_intensity_list <- lapply(index_list, function(x){
  if(length(x)<2 & sum(is.na(x))){
    avg_intensity = rep(0,1474)
  }else{
    avg_intensity <- quantification[x,5:ncol(quantification)]%>% t() %>% as.numeric()
  }
  if (length(x)>1) {
    dat <- quantification[x,5:ncol(quantification)]
    avg_intensity = apply(dat,2,mean) %>% as.numeric()
  }
  return(avg_intensity)
})

avg_intensity_df <- list.rbind(avg_intensity_list) %>% as.data.frame()
rownames(avg_intensity_df) <- coord_ST$barcode

save(avg_intensity_df,file = './transform_data/Spot_intensity_transform_T5.Rda')



