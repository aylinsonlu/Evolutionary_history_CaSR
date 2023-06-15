library(data.table)
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(ggstatsplot)
library(tidyverse)
library(gridExtra)
library(tidyverse)
library(hrbrthemes)
library(viridis)
set.seed(1234)

pdf("hp.pdf",width = 10,height = 10)

df <- fread("casr_artcile_ml_check_11_06_2023/heatmap_file.csv", data.table = T)
df
aa_names <- unique(df$substitution)
mat <- c()
for (i in 1:892){
  m1 <- df[df$position == i, ]
  m1 <- m1[match(aa_names, m1$substitution),]
  mat <- cbind(mat, m1$prediction)
}
colnames(mat) <- 1:dim(mat)[2]
rownames(mat) <- aa_names


col_labels <- colnames(mat)


for (i in 1:nrow(mat)) {
  for (j in 1:ncol(mat)) {
    if(mat[i,j] < 0.5) {mat[i,j] = "gain";next}
    if(mat[i,j] >= 0.5 & mat[i,j] < 1) {mat[i,j] = "loss";next}
    if(mat[i,j] == 1 ) {mat[i,j] = "neutral";next}
  }
}



df <- filter(df,df$position<=892)
tbl <- with(df, table(position, result))
as.data.frame(tbl)

df_bar <- as.data.frame(tbl)


ha1 = HeatmapAnnotation(
  gain = anno_barplot(filter(df_bar, result=="gain")$Freq,bar_width = 1, gp = gpar(col = NA, fill = "#4C123E"), border = FALSE, axis = FALSE),
  loss = anno_barplot(filter(df_bar, result=="loss")$Freq, bar_width = 1, gp = gpar(col = NA, fill = "#147083"), border = FALSE, axis = FALSE),
  neutral = anno_barplot(filter(df_bar, result=="neutral")$Freq, bar_width = 1, gp = gpar(col = NA, fill ="red"), border = FALSE, axis = FALSE)
)



df_u <- distinct(df, position, .keep_all = TRUE)
df_u$domain


ha = HeatmapAnnotation(foo = df_u$domain)

ht_list =  Heatmap(mat,col=c("#4C123E", "#147083", "#FFDA73"),
                   layer_fun  = function(j, i, x, y, width, height, fill) {
                     for (i in 1:nrow(mat)) {
                       for (j in 1:ncol(mat)) {
                         if(mat[i,j] < 0.5) {mat[i,j] = "gain";next}
                         if(mat[i,j] >= 0.5 & mat[i,j] < 1) {mat[i,j] = "loss";next}
                         if(mat[i,j] == 1 ) {mat[i,j] = "neutral";next}
                       }
                     }
                     
                   },
                   width = ncol(mat)*unit(0.1, "mm"), 
                   height = nrow(mat)*unit(5, "mm"),
                   cluster_rows = FALSE, cluster_columns = FALSE,
                   row_names_gp = gpar(fontsize =10),
                   cluster_column_slices=FALSE,
                   cluster_row_slices = FALSE, 
                   column_title_gp = gpar(fontfamily = "sans", fontsize = 10),
                   column_names_gp = gpar(fontsize = 0),name = "prediction")



draw(ht_list, show_heatmap_legend = TRUE)


c
dev.off()


pdf("gain_count.pdf",width = 10,height = 10)

df_counts <- fread("heatmap_counts.csv", data.table = T)


df_counts
gain_counts <- filter(df_counts,df_counts$result=="gain")
gain_counts




result_type <- unique(gain_counts$result)
result_type
mat_count <- c()
for (i in 1:length(unique(gain_counts$position))){
  m1 <- gain_counts[gain_counts$position == i, ]
  m1 <- m1[match(result_type, m1$result),]
  mat_count <- cbind(mat_count, m1$count)
}
mat_count

colnames(mat_count) <- 1:dim(mat_count)[2]
rownames(mat_count) <- result_type


col_labels <- colnames(mat_count)
col_labels
rownames

library(circlize)
col_fun = colorRamp2(c(1, 19), c("white","#4C123E"))
col_fun(seq(-3, 3))


ht_count_list =  Heatmap(mat_count, col=col_fun,
                         width = ncol(mat_count)*unit(0.1, "mm"), 
                         height = nrow(mat_count)*unit(5, "mm"),
                         cluster_rows = FALSE, cluster_columns = FALSE,
                         row_names_gp = gpar(fontsize =10),
                         cluster_column_slices=FALSE,
                         cluster_row_slices = FALSE, 
                         column_title_gp = gpar(fontfamily = "sans", fontsize = 10),
                         column_names_gp = gpar(fontsize = 0),name = "count",heatmap_legend_param = list(
                           direction = "horizontal"))
draw(ht_count_list)

dev.off()





pdf("loss_count.pdf",width = 10,height = 10)

df_counts <- fread("heatmap_counts.csv", data.table = T)


df_counts
loss_counts <- filter(df_counts,df_counts$result=="loss")
loss_counts




result_type <- unique(loss_counts$result)
result_type
mat_count <- c()
for (i in 1:length(unique(loss_counts$position))){
  m1 <- loss_counts[loss_counts$position == i, ]
  m1 <- m1[match(result_type, m1$result),]
  mat_count <- cbind(mat_count, m1$count)
}
mat_count

colnames(mat_count) <- 1:dim(mat_count)[2]
rownames(mat_count) <- result_type


col_labels <- colnames(mat_count)
col_labels
rownames


library(circlize)
col_fun = colorRamp2(c(1, 19), c("white","#147083"))
col_fun(seq(-3, 3))


ht_count_list =  Heatmap(mat_count, col=col_fun,
                         width = ncol(mat_count)*unit(0.1, "mm"), 
                         height = nrow(mat_count)*unit(5, "mm"),
                         cluster_rows = FALSE, cluster_columns = FALSE,
                         row_names_gp = gpar(fontsize =10),
                         cluster_column_slices=FALSE,
                         cluster_row_slices = FALSE, 
                         column_title_gp = gpar(fontfamily = "sans", fontsize = 10),
                         column_names_gp = gpar(fontsize = 0),name = "count",heatmap_legend_param = list(
                           direction = "horizontal"))
draw(ht_count_list)

dev.off()





###########neutral count




pdf("neutral_count.pdf",width = 10,height = 10)

df_counts <- fread("heatmap_counts.csv", data.table = T)


df_counts
neutral_counts <- filter(df_counts,df_counts$result=="neutral")
neutral_counts




result_type <- unique(neutral_counts$result)
result_type
mat_count <- c()
for (i in 1:length(unique(neutral_counts$position))){
  m1 <- neutral_counts[neutral_counts$position == i, ]
  m1 <- m1[match(result_type, m1$result),]
  mat_count <- cbind(mat_count, m1$count)
}
mat_count

colnames(mat_count) <- 1:dim(mat_count)[2]
rownames(mat_count) <- result_type

col_labels <- colnames(mat_count)


library(circlize)
col_fun = colorRamp2(c(1, 3), c("white","#FFDA73"))
col_fun(seq(-3, 3))


ht_count_list =  Heatmap(mat_count, col=col_fun,
                         width = ncol(mat_count)*unit(0.1, "mm"), 
                         height = nrow(mat_count)*unit(5, "mm"),
                         cluster_rows = FALSE, cluster_columns = FALSE,
                         row_names_gp = gpar(fontsize =10),
                         cluster_column_slices=FALSE,
                         cluster_row_slices = FALSE, 
                         column_title_gp = gpar(fontfamily = "sans", fontsize = 10),
                         column_names_gp = gpar(fontsize = 0),name = "count",heatmap_legend_param = list(
                           direction = "horizontal"))
draw(ht_count_list)

dev.off()

