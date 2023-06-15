library(xgboost)
library(ggplot2)
library(dplyr)
library(AUC)
library(ggpubr)

for (i in 1001:1050) {
  load(sprintf("with_pr_50_split_results_depth_2/casr_xgb_replication_%i_state.RData", i))
  write.csv(state$best_msg,
            file=sprintf("with_pr_50_split_results_depth_2/casr_xgb_replication_%i_state.csv", i))
  
  write.csv(c(state$params$eta,state$params$gamma,state$params$max_depth,
              
              state$params$subsample,state$params$colsample_bytree,
              state$params$min_child_weight,state$params$lambda,
              state$params$alpha,state$params$scale_pos_weight,
              state$params$objective,state$params$eval_metric),
            file=sprintf("with_pr_50_split_results_depth_2/casr_xgb_replication_%i_params.csv", i))
}
for (i in 1001:1050) {
  load(sprintf("with_pr_50_split_results_depth_2/casr_xgb_replication_%i_prauc.RData", i))
  write.csv(prauc$test_prauc$auc.integral,
            file=sprintf("with_pr_50_split_results_depth_2/casr_xgb_replication_%i_test_prauc.csv", i))
  write.csv(prauc$train_prauc$auc.integral,
            file=sprintf("with_pr_50_split_results_depth_2/casr_xgb_replication_%i_train_prauc.csv", i))
  
  
}
