library(AUC)
library(xgboost)
library(PRROC)
source("./grid_search_cv.R")
source("./varImp_cv.R")

args <- commandArgs(trailingOnly = TRUE)
replication <- as.numeric(args[[1]])

id <- 1000 + replication
result_path <- "./with_pr_50_split_results_depth_2"
set.seed(1903 * replication)

X_train <- read.csv(file = sprintf("./train_csv_files/%s_train.csv", id))
y_train <- 2 * (X_train$consequence == "loss-of-function") - 1
X_train <- X_train[,-c(1,86)]
#X_train <- X_train[,-c(which(apply(X_train, 2, sd) == 0))]

X_test <- read.csv(file = sprintf("./train_csv_files/%s_test.csv", id))
y_test <- 2 * (X_test$consequence == "loss-of-function") - 1
X_test <- X_test[,-c(1,86)]
#X_test <- X_test[,-c(which(apply(X_train, 2, sd) == 0))]


train_negative_indices <- which(y_train == -1)
train_positive_indices <- which(y_train == 1)

####Algorithm parameters
# eta <- 0.1
# nrounds <- 250

maxdepth <- 2
childw <- 1
gamma <- 0
subsample <- 0.8
colsample_bytree <- 0.8
lambda <- 1
alpha <- 0

epsilon <- 1e-3
fold_count <- 5
train_ratio <- 0.8
####


negative_allocation <- sample(rep(1:fold_count, ceiling(length(train_negative_indices) / fold_count)), length(train_negative_indices))
positive_allocation <- sample(rep(1:fold_count, ceiling(length(train_positive_indices) / fold_count)), length(train_positive_indices))

fold_indices_train <- list()
fold_indices_test <- list()
for (fold in 1:fold_count) {
  train_indices <- c(train_negative_indices[which(negative_allocation != fold)], train_positive_indices[which(positive_allocation != fold)])
  test_indices <- c(train_negative_indices[which(negative_allocation == fold)], train_positive_indices[which(positive_allocation == fold)])
  
  fold_indices_train[[fold]] <- train_indices
  fold_indices_test[[fold]] <- test_indices
}


###############feature select
# parameters <- list()
# parameters$eta <- 0.0005
# parameters$gamma <- gamma
# parameters$max_depth <- maxdepth
# parameters$subsample <- subsample
# parameters$colsample_bytree <- colsample_bytree
# parameters$min_child_weight <- childw
# parameters$lambda <- lambda
# parameters$alpha <- alpha
# # parameters$scale_pos_weight <- sum(y_train2 == -1) / sum(y_train2 == 1)
# parameters$booster <- "gbtree"
# 
# nrounds <- 200
# ####
# features <- colnames(X_train)
# model <- varImp_cv(features, X_train, y_train, fold_count, fold_indices_train, fold_indices_test,  parameters, nrounds, id)
# 
# mean_auc <- mean(model$auroc_matrix)
# varImp_mean <- rowMeans(model$feature_importance_mat)
# 
# while (1) {
#   thr <- 0.015
#   
#   X_tr_temp <- X_train[,((varImp_mean <= thr) == F)]
#   X_test_temp <- X_test[, ((varImp_mean <= thr) == F)]
#   features <- colnames(X_tr_temp)
#   temp_model <- varImp_cv(features, X_tr_temp, y_train, fold_count, fold_indices_train, fold_indices_test,  parameters, nrounds, id)
#   
#   temp_auc <- mean(temp_model$auroc_matrix)
#   temp_varImp <- rowMeans(temp_model$feature_importance_mat)
#   
#   while(temp_auc <= mean_auc) {
#     thr <- round(thr - 0.005)
#     
#     if (thr == 0) {
#       break
#     }
#     
#     X_tr_temp <- X_train[,((varImp_mean <= thr) == F)]
#     X_test_temp <- X_test[,((varImp_mean <= thr) == F)]
#     
#     temp_model <- varImp_cv(features, X_tr_temp, y_train, fold_count, fold_indices_train, fold_indices_test,  parameters, nrounds, id)
#     
#     temp_auc <- mean(temp_model$auroc_matrix)
#     temp_varImp <- rowMeans(temp_model$feature_importance_mat)
#   }
#   
#   if (temp_auc > mean_auc) {
#     X_train <- X_tr_temp
#     X_test <- X_test_temp
#     varImp_mean <- temp_varImp
#     mean_auc <- temp_auc
#   }else {
#     break
#   }
# }
# 
# dim(X_train)
# ###############

print("Eta and nrounds are being tuned")

#eta_set <- c(0.00001,0.0001, 0.001, 0.01, 0.1)

# eta_set <- c(1e-4, 1e-3, 1e-2)
# eta_set <- c(seq(0.0001, 0.01, 0.00005), seq(0.015, 0.2, 0.005))
# eta_set <- c(seq(0.001, 0.01, 0.001), seq(0.02, 0.5, 0.01))
# eta_set <- c(seq(0.0001, 0.01, 0.0001), seq(0.01, 0.1, 0.01))
# eta_set <- c(sapply(c(1:5), function(i){i*c(1e-4,1e-3,1e-2,1e-1)}))
eta_set <- c(seq(0.00001, 0.001, 0.00001))
# eta_set <- unique(c(seq(1e-4, 1e-3, 1e-4), seq(1e-3,1e-2, 1e-3), seq(1e-2, 1e-1, 1e-2)))
nrounds_set <- c(200) #n_estimator


nrow_mat <- length(eta_set)*length(nrounds_set)
auroc_matrix <- array(NA, dim = c(nrow_mat,7,fold_count), dimnames = list(c(1:nrow_mat), c("train_auc", "test_auc", "train_loss", "test_loss", "eta", "nrounds", "best_nrounds"), c(1:fold_count)))

for (fold in 1:fold_count) {
  # train_indices <- fold_indices_train[[fold]]
  # test_indices <- fold_indices_test[[fold]]
  # 
  # X_train2 <- as.matrix(X_train[train_indices,])
  # X_test2 <- as.matrix(X_train[test_indices,])
  # X_train2 <- scale(X_train2)
  # X_test2 <- (X_test2 - matrix(attr(X_train2, "scaled:center"), nrow = nrow(X_test2), ncol = ncol(X_test2), byrow = TRUE)) / matrix(attr(X_train2, "scaled:scale"), nrow = nrow(X_test2), ncol = ncol(X_test2), byrow = TRUE)
  # 
  # y_train2 <- y_train[train_indices]
  # y_test2 <- y_train[test_indices]
  
  X_train2 <- read.csv(sprintf("./train_csv_files/%s_%s_train.csv", id, fold))
  y_train2 <- 2 * (X_train2$consequence == "loss-of-function") - 1
  X_train2 <- X_train2[,-c(1,86)]
  X_train2 <- scale(X_train2)
  
  X_test2 <- read.csv(sprintf("./train_csv_files/%s_%s_val.csv", id, fold))
  y_test2 <- 2 * (X_test2$consequence == "loss-of-function") - 1
  X_test2 <- X_test2[,-c(1,86)]
  X_test2 <- (X_test2 - matrix(attr(X_train2, "scaled:center"), nrow = nrow(X_test2), ncol = ncol(X_test2), byrow = TRUE)) / matrix(attr(X_train2, "scaled:scale"), nrow = nrow(X_test2), ncol = ncol(X_test2), byrow = TRUE)
  
  i <- 0
  
  dtrain <- xgb.DMatrix(data = as.matrix(X_train2), label = 1*(y_train2==1))
  dtest <- xgb.DMatrix(data = as.matrix(X_test2), label = 1*(y_test2==1))
  watchlist <- list(train=dtrain, test=dtest)
  
  print(sprintf("running replication = %d, fold = %d", replication, fold))
  for (eta in eta_set) {
    for (nrounds in nrounds_set) {
      i <- i+1
      
      parameters <- list()
      parameters$eta <- eta
      parameters$gamma <- gamma
      parameters$max_depth <- maxdepth
      parameters$subsample <- subsample
      parameters$colsample_bytree <- colsample_bytree
      parameters$min_child_weight <- childw
      parameters$lambda <- lambda
      parameters$alpha <- alpha
      parameters$scale_pos_weight <- sum(y_train2 == -1) / sum(y_train2 == 1)
      parameters$booster <- "gbtree"
      
      state <- xgb.train(data = dtrain,  params = parameters,  watchlist = watchlist, early_stopping_rounds = 50,
                         nrounds = nrounds, objective = "binary:logistic", verbose = 1, set.seed(1), eval_metric = "auc", maximize = T)          
      prediction <- predict(state, as.matrix(X_test2))
      auroc <- auc(roc(prediction, as.factor(1 * (y_test2 == +1))))
      train_pred <- predict(state, as.matrix(X_train2))
      train_auroc <- auc(roc(train_pred, as.factor(1 * (y_train2 == +1))))
      train_loss <- -sum(log(train_pred)*(1*(y_train2==1)) + (1*(y_train2==-1) * log(1 - train_pred)))/length(y_train2)
      test_loss <- -sum(log(prediction)*(1*(y_test2==1)) + (1*(y_test2==-1) * log(1 - prediction)))/length(y_test2)
      
      auroc_matrix[i, ,fold] <- c(train_auroc, auroc, train_loss, test_loss, eta, nrounds, state$best_iteration)
    }
  }
}

mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"test_auc"]), ties.method = "last"), ]

eta_star <- best_auroc["eta"]
#nrounds_star <- ceiling(best_auroc["best_nrounds"])
nrounds_star <- c(200)
eta_set <- eta_star
depth_set <- c(2)
childw_set <- c(1:6)
gamma_set <- c(0)
subsample_set <- c(0.8)
colsample_bytree_set <- c(0.8)
lambda_set <- c(1)
alpha_set <- c(0)

print("Max depth and min child weight are being tuned")
auroc_matrix <- grid_search_cv(X_train, y_train, fold_count, id, fold_indices_train, fold_indices_test, eta_set, nrounds_star, gamma_set, subsample_set, colsample_bytree_set, childw_set, lambda_set, alpha_set, depth_set)  

mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]
maxdepth_star <- best_auroc["depth"]
childw_star <- best_auroc["childw"]

# if (maxdepth_star == depth_set[length(depth_set)]){
#   depth_set <- c((maxdepth_star ): (maxdepth_star + 2))
# }

if (childw_star == childw_set[length(childw_set)]){
  childw_set <- c((childw_star - 1): (childw_star + 2))
}

auroc_matrix <- grid_search_cv(X_train, y_train, fold_count, id, fold_indices_train, fold_indices_test, eta_set, nrounds_star, gamma_set, subsample_set, colsample_bytree_set, childw_set, lambda_set, alpha_set, depth_set)

mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]
maxdepth_star <- best_auroc["depth"]
childw_star <- best_auroc["childw"]

depth_set <- maxdepth_star
childw_set <- childw_star

########
########UPDATE nrounds 
# 
auroc_matrix <- grid_search_cv(X_train, y_train, fold_count, id, fold_indices_train, fold_indices_test, eta_set, nrounds_star, gamma_set, subsample_set, colsample_bytree_set, childw_set, lambda_set, alpha_set, depth_set)
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]
nrounds_star <- ceiling(best_auroc["best_nrounds"])

#######
#######


gamma_set <- seq(0,0.5,0.1)

print("gamma is being tuned")

auroc_matrix <- grid_search_cv(X_train, y_train, fold_count, id, fold_indices_train, fold_indices_test, eta_set, nrounds_star, gamma_set, subsample_set, colsample_bytree_set, childw_set, lambda_set, alpha_set, depth_set)  
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]
gamma_star <- best_auroc["gamma"]

if(gamma_star == gamma_set[1]){
  gamma_set <- c((gamma_star), (gamma_star + 0.05)) 
} else if (maxdepth_star == gamma_set[length(gamma_set)]){
  gamma_set <- seq((gamma_star - 0.05), (gamma_star + 0.3), 0.05) 
}else {
  gamma_set <- c((gamma_star - 0.05), gamma_star, (gamma_star + 0.05)) 
}

auroc_matrix <- grid_search_cv(X_train, y_train, fold_count, id, fold_indices_train, fold_indices_test, eta_set, nrounds_star, gamma_set, subsample_set, colsample_bytree_set, childw_set, lambda_set, alpha_set, depth_set)  
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]
gamma_star <- best_auroc["gamma"]
gamma_set <- gamma_star

maxdepth <- depth_set
childw <- childw_set
gamma <- gamma_set
subsample <- 0.8
colsample_bytree <- 0.8

colsample_bytree_set <- seq(0.5,1,0.05)
subsample_set <- seq(0.5,1,0.05)

print("Colsample by tree and subsample are being tuned")

auroc_matrix <- grid_search_cv(X_train, y_train, fold_count, id, fold_indices_train, fold_indices_test, eta_set, nrounds_star, gamma_set, subsample_set, colsample_bytree_set, childw_set, lambda_set, alpha_set, depth_set)  
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]
colsample_bytree_star <- best_auroc["colsample_bytree"]
subsample_star <- best_auroc["subsample"]

colsample_bytree_set <- colsample_bytree_star 
subsample_set <- subsample_star

#lambda_set <- c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1, 5, 10, 15, 20, 30, 50,100)
lambda_set <- c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100)

print("Lambda is being tuned")
auroc_matrix <- grid_search_cv(X_train, y_train, fold_count, id, fold_indices_train, fold_indices_test, eta_set, nrounds_star, gamma_set, subsample_set, colsample_bytree_set, childw_set, lambda_set, alpha_set, depth_set)  
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]
lambda_star <- best_auroc["lambda"]

if (lambda_star == 0) {
  lambda_set <- c(0, 1e-4*0.5)
} else if (lambda_star == 100) {
  lambda_set <- c(50, 100, 150, 200, 250)
} else(
  lambda_set <- c(lambda_star*0.2, lambda_star*0.5, lambda_star*0.75, lambda_star, lambda_star*2, lambda_star*5, lambda_star*7.5)
)

auroc_matrix <- grid_search_cv(X_train, y_train, fold_count, id, fold_indices_train, fold_indices_test, eta_set, nrounds_star, gamma_set, subsample_set, colsample_bytree_set, childw_set, lambda_set, alpha_set, depth_set)  
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]
lambda_star <- best_auroc["lambda"]
lambda_set <- lambda_star

alpha_set <- c(0, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100)
print("Alpha is being tuned")
auroc_matrix <- grid_search_cv(X_train, y_train, fold_count, id, fold_indices_train, fold_indices_test, eta_set, nrounds_star, gamma_set, subsample_set, colsample_bytree_set, childw_set, lambda_set, alpha_set, depth_set)  
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]
alpha_star <- best_auroc["alpha"]

if (alpha_star == 0) {
  alpha_set <- c(0, 1e-4*0.5)
} else if (alpha_star == 100) {
  alpha_set <- c(50, 100, 150, 200, 250)
} else(
  alpha_set <- c(alpha_star*0.2, alpha_star*0.5, alpha_star*0.75, alpha_star, alpha_star*2, alpha_star*5, alpha_star*7.5)
)

auroc_matrix <- grid_search_cv(X_train, y_train, fold_count, id, fold_indices_train, fold_indices_test, eta_set, nrounds_star, gamma_set, subsample_set, colsample_bytree_set, childw_set, lambda_set, alpha_set, depth_set)  
mean_auroc <- apply(auroc_matrix, c(1,2), mean)
best_auroc <- mean_auroc[max.col(t(mean_auroc[,"auc"]), ties.method = "last"), ]
alpha_star <- best_auroc["alpha"]
alpha_set <- alpha_star
nrounds_star <- ceiling(2*best_auroc["best_nrounds"])
print("best model is being trained")
parameters <- list()
parameters$eta <- eta_star
parameters$gamma <- gamma_star
parameters$max_depth <- maxdepth_star
parameters$subsample <- subsample_star
parameters$colsample_bytree <- colsample_bytree_star
parameters$min_child_weight <- childw_star
parameters$lambda <- lambda_star
parameters$alpha <- alpha_star
parameters$scale_pos_weight <- sum(y_train == -1) / sum(y_train == 1)
parameters$booster <- "gbtree"  

dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = 1*(y_train==1))
dtest <- xgb.DMatrix(data = as.matrix(X_test), label = 1*(y_test==1))
watchlist <- list(train=dtrain, test=dtest)

state <- xgb.train(data = dtrain, params = parameters, watchlist = watchlist, 
                   nrounds = nrounds_star,early_stopping_rounds = ceiling(nrounds_star * 0.2), objective = "binary:logistic", verbose = 1, eval_metric = "auc", maximize = T, set.seed(1))          



train_prediction <- predict(state, as.matrix(X_train))
train_auroc <- auc(roc(train_prediction, as.factor(1 * (y_train == 1))))
train_prauc <- pr.curve(scores.class0 = as.numeric(train_prediction), weights.class0 = (1 * (y_train == 1)), curve = T)

test_prediction <- predict(state, as.matrix(X_test))
test_auroc <- auc(roc(test_prediction, as.factor(1 * (y_test == 1))))
test_prauc <- pr.curve(scores.class0 = as.numeric(test_prediction), weights.class0 = (1 * (y_test == 1)), curve = T)


prediction <- list(train_prediction = train_prediction, test_prediction = test_prediction)

result <- list(train_auroc = train_auroc, test_auroc = test_auroc)
prauc <- list(train_prauc = train_prauc, test_prauc = test_prauc)

save("state", file = sprintf("%s/casr_xgb_replication_%d_state.RData", result_path, id))
save("prediction", file = sprintf("%s/casr_xgb_replication_%d_prediction.RData", result_path, id))
save("result", file = sprintf("%s/casr_xgb_replication_%d_result.RData", result_path, id))
save("prauc", file = sprintf("%s/casr_xgb_replication_%d_prauc.RData", result_path, id))
xgb.save(state,sprintf("%s/casr_xgb_replication_%d_xgb.model", result_path, id))
