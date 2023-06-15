library(AUC)
library(xgboost)
library(PRROC)
# library(CORElearn)
library(CalibratR)

setwd("")

args <- commandArgs(trailingOnly = TRUE)
replication <- 1


id <- 2000 + replication
result_path <- "./uniprot_whole_test_cases"
set.seed(1903 * replication)

X_train <- read.csv(file = sprintf("./uniprot_whole_test_cases/%s_train.csv", id))
y_train <- 2 * (X_train$consequence == "loss-of-function") - 1



X_train <- X_train[,-c(1,86)]

X_train <- scale(X_train)

X_test <- read.csv("./uniprot_whole_test_cases/uniprot_final_unlabeled_test.csv")

X_test <- X_test[,-c(1,86)]
X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)


train_negative_indices <- which(y_train == -1)
train_positive_indices <- which(y_train == 1)

####Algorithm parameters

parameters <- list()
parameters$eta <- c(0.00064)
parameters$gamma <- c(0.55)
parameters$max_depth <- c(2)

parameters$subsample <- c(0.55)
parameters$colsample_bytree <- c(0.6)
parameters$min_child_weight <- c(2)
parameters$lambda <- c(1)
parameters$alpha <- c(0.1)
parameters$scale_pos_weight <- sum(y_train == -1) / sum(y_train == 1)
parameters$booster <- "gbtree"



dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = 1*(y_train==1))
dtest <- xgb.DMatrix(data = as.matrix(X_test))
dtest
watchlist <- list(train=dtrain, test=dtest)
dtrain
state <- xgb.train(data = dtrain, params = parameters,
                   nrounds = c(17), objective = "binary:logistic", verbose = 1, eval_metric = "auc", maximize = T, set.seed(1))

train_prediction <- predict(state, as.matrix(X_train))

test_prediction <- predict(state, as.matrix(X_test))

prediction <- list(train_prediction = train_prediction, test_prediction = test_prediction)
train_prediction
test_prediction
write.csv(test_prediction,
          file=("whole_test_predictions.csv"))
          
          write.csv(train_prediction,
                    file=("whole_train_predictions.csv"))
                    
test_pos <- list(c(15809,15221,7205,1089,1803,3382,15497,8242,12848,7964,2682,7298,14843,1099,15119,11164,733,2007,2540,2983,12888,12930,14911))

for (i in test_pos){
  print(2*(test_prediction[i] >= 0.5)-1)
}
