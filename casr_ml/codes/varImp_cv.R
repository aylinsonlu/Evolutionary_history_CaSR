varImp_cv <- function(features, X_train, y_train, fold_count, fold_indices_train, fold_indices_test,  parameters, nrounds, id) {
  N <- dim(X_train)[2]
  
  feature_importance_mat <- array(0, dim = c(N,fold_count), dimnames = list(colnames(X_train), c(1:fold_count)))
  auroc_matrix <- array(NA, dim = c(1,fold_count), dimnames = list("auc", c(1:fold_count)))
  
  for (fold in 1:fold_count) {
    # train_indices <- fold_indices_train[[fold]]
    # test_indices <- fold_indices_test[[fold]]
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
    X_train2 <- X_train2[,features]
    X_train2 <- scale(X_train2)
    
    X_test2 <- read.csv(sprintf("./train_csv_files/%s_%s_val.csv", id, fold))
    y_test2 <- 2 * (X_test2$consequence == "loss-of-function") - 1
    X_test2 <- X_test2[,-c(1,86)]
    X_test2 <- X_test2[,features]
    
    X_test2 <- (X_test2 - matrix(attr(X_train2, "scaled:center"), nrow = nrow(X_test2), ncol = ncol(X_test2), byrow = TRUE)) / matrix(attr(X_train2, "scaled:scale"), nrow = nrow(X_test2), ncol = ncol(X_test2), byrow = TRUE)
    
    i <- 0
    
    dtrain <- xgb.DMatrix(data = as.matrix(X_train2), label = 1*(y_train2==1))
    dtest <- xgb.DMatrix(data = as.matrix(X_test2), label = 1*(y_test2==1))
    watchlist <- list(train=dtrain, test=dtest)
    
    # print(sprintf("running replication = %d, fold = %d", replication, fold))
    
    # parameters <- list()
    # parameters$eta <- eta
    # parameters$gamma <- gamma
    # parameters$max_depth <- maxdepth
    # parameters$subsample <- subsample
    # parameters$colsample_bytree <- colsample_bytree
    # parameters$min_child_weight <- childw
    # parameters$lambda <- lambda
    # parameters$alpha <- alpha
    parameters$scale_pos_weight <- sum(y_train2 == -1) / sum(y_train2 == 1)
    # parameters$booster <- "gbtree"
    
    state <- xgb.train(data = dtrain,  params = parameters, early_stopping_rounds = ceiling(nrounds*0.2), watchlist = watchlist, 
                       nrounds = nrounds, objective = "binary:logistic", verbose = 1, set.seed(1), eval_metric = "auc", maximize = T)          
    prediction <- predict(state, as.matrix(X_test2))
    auroc <- auc(roc(prediction, as.factor(1 * (y_test2 == +1))))
    train_pred <- predict(state, as.matrix(X_train2))
    train_auroc <- auc(roc(train_pred, as.factor(1 * (y_train2 == +1))))
    
    auroc_matrix[1,fold] <- auroc
    importance_matrix <- xgb.importance(model = state)
    feature_importance_mat[as.matrix(importance_matrix[,1]),fold] <-  unlist(importance_matrix[,2])
  }
  
  model <- list(auroc_matrix = auroc_matrix, feature_importance_mat = feature_importance_mat)
  
  return(model)
}

