library(xgboost)
library(ggplot2)
library(dplyr)
library(AUC)
library(ggpubr)
#library(SHAPforxgboost) 
library(data.table)
library(shapviz)



cnames <- c("value", "rfvalue", "stdfvalue", "mean_value")


# 
mod <- xgb.load("results_final_feature/casr_xgb_replication_2001_xgb.model")
mod
X_train <- read.csv("final_csv_for_feature_selection/2001_train.csv")
X_train <- X_train[,-c(1,86)]
X_train

dataX <- as.matrix(X_train)
x <- shapviz(mod, X_pred = dataX)
sv_importance(x, kind = "beeswarm",show_other = FALSE,show_numbers = TRUE,viridis_args = list(begin = 0.25, end = 0.85, option = "mako"))+ theme_classic()
ggsave("shapley_values.pdf")
x$S
ggsave("shapley_values.pdf",width = 10,height = 10)

dev.off()
