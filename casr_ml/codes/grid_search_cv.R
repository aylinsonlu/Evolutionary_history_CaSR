grid_search_cv <- function(X, y, fold_count, id, fold_indices_train, fold_indices_test, eta_set, nrounds, gamma_set, subsample_set, colsample_bytree_set, childw_set, lambda_set, alpha_set, depth_set) {
  nrow_mat <- length(eta_set)*length(gamma_set)*length(depth_set)*length(subsample_set)*length(colsample_bytree_set)*length(childw_set)*length(lambda_set)*length(alpha_set)
  # auroc_matrix <- array(NA, dim = c(nrow_mat,9,fold_count), dimnames = list(c(1:nrow_mat), c("auc", "eta", "gamma", "depth", "subsample", "colsample_bytree", "childw", "lambda", "alpha"), c(1:fold_count)))
  auroc_matrix <- array(NA, dim = c(nrow_mat,11,fold_count), dimnames = list(c(1:nrow_mat), c("train_auc", "auc", "eta", "gamma", "depth", "subsample", "colsample_bytree", "childw", "lambda", "alpha", "best_nrounds"), c(1:fold_count)))
  
  for (fold in 1:fold_count) {
    # train_indices <- fold_indices_train[[fold]]
    # test_indices <- fold_indices_test[[fold]]
    # X_train <- as.matrix(X[train_indices,])
    # X_test <- as.matrix(X[test_indices,])
    # X_train <- scale(X_train)
    # X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
    # 
    # y_train <- y[train_indices]
    # y_test <- y[test_indices]
    
    X_train <- read.csv(sprintf("./train_csv_files/%s_%s_train.csv", id, fold))
    y_train <- 2 * (X_train$consequence == "loss-of-function") - 1
    X_train <- X_train[,-c(1,86)]
    X_train <- scale(X_train)
    
    X_test <- read.csv(sprintf("./train_csv_files/%s_%s_val.csv", id, fold))
    y_test <- 2 * (X_test$consequence == "loss-of-function") - 1
    X_test <- X_test[,-c(1,86)]
    X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
    
    dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = 1*(y_train==1))
    dtest <- xgb.DMatrix(data = as.matrix(X_test), label = 1*(y_test==1))
    watchlist <- list(train=dtrain, test=dtest)
    
    print(sprintf("running replication = %d, fold = %d", replication, fold))
    i <- 0
    for (eta in eta_set) {
      for (gamma in gamma_set) {
        for (maxdepth in depth_set) {
          for (subsample in subsample_set) {
            for (colsample_bytree in colsample_bytree_set) {
              for (childw in childw_set) {
                for (lambda in lambda_set) {
                  for (alpha in alpha_set) {
                    t1 <- Sys.time()
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
                    parameters$scale_pos_weight <- sum(y_train == -1) / sum(y_train == 1)
                    parameters$booster <- "gbtree"
                    # parameters$objective <- "binary:logistic"
                    # parameters$eval_metric <- "auc"
                    # parameters$verbose <- 1
                    #
                    state <- xgb.train(data = dtrain, params = parameters, watchlist = watchlist, early_stopping_rounds = ceiling(nrounds*0.2),
                                       nrounds = nrounds, objective = "binary:logistic", verbose = 1, eval_metric = "auc", maximize = T, set.seed(1))          
                    prediction <- predict(state, as.matrix(X_test))
                    train_pred <- predict(state, as.matrix(X_train))
                    train_auroc <- auc(roc(train_pred, as.factor(1 * (y_train == +1))))
                    auroc <- auc(roc(prediction, as.factor(1 * (y_test == +1))))
                    auroc_matrix[i, ,fold] <- c(train_auroc, auroc, eta, gamma, maxdepth, subsample, colsample_bytree, childw, lambda, alpha, state$niter)
                    # auroc_matrix[i, ,fold] <- c(train_auroc, auroc, eta, gamma, maxdepth, subsample, colsample_bytree, childw, lambda, alpha, state$best_iteration)
                    t2 <- Sys.time()
                    print(t2-t1)
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(auroc_matrix)
}