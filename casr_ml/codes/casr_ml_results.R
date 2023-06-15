library(xgboost)
library(ggplot2)
library(dplyr)
library(AUC)
library(ggpubr)

data <- read.csv("/with_pr_50_split_results_depth_2/casr_xgb_train_test_auc.csv")
data
data2 <-read.csv("/with_pr_50_split_results_depth_2/casr_xgb_train_test_auc_2.csv")
data2
mean_train <- mean(data2$train_auc)
mean_train <- signif(mean_train,3)
mean_train
mean_val <- mean(data2$test_auc)
mean_val <- signif(mean_val,3)
mean_val






p<-data %>%
  ggplot( aes(x=reorder(type, -AUROC), y=AUROC, fill=type)) +
  geom_boxplot(notch = TRUE) +
  scale_fill_manual("dataset", values = c("train" = "#023E8A", "test" = "#8DB4F2"))+
  
  labs(x = "replication", y = "AUROC",colour = "")+
  theme(text = element_text(size = 20),plot.subtitle = element_text(size = 20),plot.title = element_text(size = 20),panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"))+
  #ggtitle("Basic boxplot") +
  xlab("")+annotate("text",size = 20/.pt, x = 1, y = 0.9, label = paste("AUC =", round(mean_train, 2)))+ annotate("text", size = 20/.pt,x = 2, y = .9, label = paste("AUC =", round(mean_val, 2)))
p<-p+ geom_hline(yintercept=0.5, linetype="dashed", color = "black")
p
ggsave("/Users/aylin/Desktop/casr_maxdepth2_50split_aucroc.pdf")



# ## Plot TRAIN-VAL AUCPR


pdf(file = "/Users/aylin/Desktop/AUCPR-50rep.pdf")
data <- read.csv("/with_pr_50_split_results_depth_2/casr_xgb_train_test_prauc.csv")
data
data2 <-read.csv("/with_pr_50_split_results_depth_2/casr_xgb_train_test_prauc_2.csv")
data2
mean_train <- mean(data2$train_pr_auc)
mean_train <- signif(mean_train,3)
mean_train
mean_val <- mean(data2$test_pr_auc)
mean_val <- signif(mean_val,3)
mean_val





h<-data %>%
  ggplot( aes(x=reorder(type, -AUPR), y=AUPR, fill=type)) +
  geom_boxplot(notch = TRUE) +
  scale_fill_manual("dataset", values = c("train" = "#023E8A", "test" = "#8DB4F2"))+
  
  labs(x = "replication", y = "AUPR",colour = "")+
  scale_y_continuous(breaks = c(0.7,0.8,0.9,1.0), limits = c(0.7,1))+
  theme(text = element_text(size = 20),plot.subtitle = element_text(size = 12),plot.title = element_text(size = 12),panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"))+
  #ggtitle("Basic boxplot") +
  xlab("")+annotate(size = 20/.pt,"text", x = 1, y = 0.9999, label = paste("AUPR =", round(mean_train, 2)))+ annotate(size = 20/.pt,"text", x = 2, y = 0.9999, label = paste("AUPR =", round(mean_val, 2)))

h<-h+ geom_hline(yintercept=243/(94+243), linetype="dashed", color = "black")

h
#scale_y_continuous(sec.axis = dup_axis(breaks = round(243/(94+243),2), labels =round(243/(94+243),2)))
ggsave("/Users/aylin/Desktop/casr_maxdepth2_50split_aupr.pdf")





