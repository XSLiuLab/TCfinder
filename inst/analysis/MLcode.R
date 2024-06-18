
### RF
library(dplyr)
library(data.table)
library(randomForest)
remove(list = ls())
setwd("~/project/mcIdentify/Revise_1/")
data <- fread("./new_data/GSE131928_new_data_score.csv",data.table = F) 
set.seed(123) 
split <- sample.split(data$type, SplitRatio = 0.8) 
train_data <- subset(data, split == TRUE) 
test_data <- subset(data, split == FALSE)

X_train <- train_data[, -1]
y_train <- as.factor(train_data[, 1])


ctrl <- trainControl(method = "cv", number = 5)
grid <- expand.grid(mtry = c(2, 4, 6)) 
rf_model <- train(x = X_train, y = y_train,
                  method = "rf",
                  trControl = ctrl,
                  tuneGrid = grid)

print(rf_model)

grid <- expand.grid(mtry = c(6)) 
modellist <- list()
for (ntree in c(100,200, 300)) {
  set.seed(123)
  fit <- train(x = X_train, y = y_train, method="rf", 
               metric="Accuracy", tuneGrid=grid, 
               trControl=ctrl, ntree=ntree)
  key <- toString(ntree)
  modellist[[key]] <- fit
}
results <- resamples(modellist)
summary(results)

model <- randomForest(x = X_train, y = y_train,mtry = 6,ntree = 200)
print(model)

x_test_data <- test_data[, -1]
y_test_data <- as.factor(test_data[, 1])
test_predictions <- predict(model, newdata = x_test_data)

confusion_matrix <- confusionMatrix(test_predictions, y_test_data)
accuracy <- confusion_matrix$overall["Accuracy"]
precision <- confusion_matrix$byClass["Pos Pred Value"]
recall <- confusion_matrix$byClass["Sensitivity"]
f1_score <- confusion_matrix$byClass["F1"]

print(confusion_matrix)
print(paste("Accuracy:", accuracy))
print(paste("Precision:", precision))
print(paste("Recall:", recall))
print(paste("F1 Score:", f1_score))






## SVM
library(e1071)
library(ggplot2)
library(caret)
remove(list = ls())
setwd("~/project/mcIdentify/Revise_1/")
data <- fread("./new_data/GSE131928_new_data_score.csv",data.table = F) 

set.seed(123)  
split <- sample.split(data$type, SplitRatio = 0.8)  
train_data <- subset(data, split == TRUE)  
test_data <- subset(data, split == FALSE) 
X_train <- train_data[, -1]
y_train <- as.factor(train_data[, 1])



param_grid <- expand.grid(
  sigma = c(0.1, 1, 10),
  C = c(0.1, 1, 10))
ctrl <- trainControl(method = "cv", number = 5, verboseIter = FALSE)

tuned_model <- train(
  x = X_train,
  y = y_train,
  method = "svmRadial",
  tuneGrid = param_grid,
  trControl = ctrl)
print(tuned_model)


svm_model <- svm(x = X_train, y = y_train,sigma = 0.1,C = 10)
train_predictions <- predict(svm_model, newdata = X_train)
table(y_train, train_predictions)

X_new_data <- test_data[, -1]
y_new_data <- as.factor(test_data[, 1])
test_predictions <- predict(svm_model, newdata = X_new_data)

table(test_predictions,y_new_data)
accuracy <- mean(test_predictions == y_new_data)
precision <- sum(test_predictions == "normal" & y_new_data == "normal") / sum(test_predictions == "normal")
recall <- sum(test_predictions == "normal" & y_new_data == "normal") / sum(y_new_data == "normal")
f1_score <- 2 * precision * recall / (precision + recall)
print(paste("accuracy:", accuracy))
print(paste("precision:", precision))
print(paste("recall:", recall))
print(paste("F1 score:", f1_score))




### xgboost
library(xgboost)
library(Matrix)
remove(list = ls())
setwd("~/project/mcIdentify/Revise_1/")
new_data <- fread("./new_data/GSE131928_new_data_score.csv",data.table = F)
new_data1 <- new_data %>% mutate(type = ifelse(type == "normal",0,1))

data <- fread("./new_data/GSE131928_new_data_score.csv",data.table = F)  
data <- data %>% mutate(type = ifelse(type == "normal",0,1))
set.seed(123) 
split <- sample.split(data$type, SplitRatio = 0.8) 
train_data <- subset(data, split == TRUE) 
test_data <- subset(data, split == FALSE)

X_train <- train_data[, -1]
y_train <- train_data[, 1]

ctrl <- trainControl(
  method = "cv", 
  number = 5,    
  verboseIter = FALSE)

param_grid <- expand.grid(
  nrounds = c(100, 200), 
  max_depth = c(3, 6), 
  eta = c(0.1), 
  gamma = c(0, 0.1),
  colsample_bytree = c(0.8),
  min_child_weight = c(1, 3),
  subsample = c(0.8))

xgb_model <- train(
  x = X_train,
  y = y_train,
  method = "xgbTree",
  trControl = ctrl,
  tuneGrid = param_grid)
print(xgb_model$bestTune)


dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
params <- list(objective = "binary:logistic", eval_metric = "logloss", eta = 0.1, max_depth = 3)
nrounds <- 100
xgb_model <- xgboost(params = params, data = dtrain, nrounds = nrounds)

train_predictions <- predict(xgb_model, newdata = dtrain)
train_predictions <- ifelse(train_predictions > 0.5,1,0)

confusion_matrix <- table(train_predictions,y_train)
accuracy <- mean(train_predictions == y_train)
precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
f1_score <- 2 * precision * recall / (precision + recall)

print(paste("accuracy:", accuracy))
print(paste("precision:", precision))
print(paste("recall:", recall))
print(paste("F1 score:", f1_score))

X_new_data1 <- new_data1[, -1]
y_new_data1 <- as.factor(new_data1[, 1])
dtest <- xgb.DMatrix(data = as.matrix(X_new_data1))
test_predictions <- predict(xgb_model, newdata = dtest)
test_predictions <- ifelse(test_predictions > 0.5,1,0)

confusion_matrix <- table(test_predictions,y_new_data1)
accuracy <- mean(test_predictions == y_new_data1)
precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
f1_score <- 2 * precision * recall / (precision + recall)

print(paste("accuracy:", accuracy))
print(paste("precision:", precision))
print(paste("recall:", recall))
print(paste("F1 score:", f1_score))



###LR
library(ggplot2)
library(dplyr)
library(caTools)
library(pROC)
library(caret)
remove(list = ls())
setwd("~/project/mcIdentify/Revise_1/")
data <- fread("./new_data/GSE131928_new_data_score.csv",data.table = F) 
data <- data %>% mutate(type = ifelse(type=="normal",1,0))
set.seed(123) 
split <- sample.split(data$type, SplitRatio = 0.7) 
train_data <- subset(data, split == TRUE)  
test_data <- subset(data, split == FALSE)  
model <- glm(type ~ ., data = train_data, family = gaussian)
summary(model)


predictions <- predict(model, newdata = test_data, type = "response")
threshold <- 0.5 
predicted_classes <- ifelse(predictions >= threshold, 1, 0)

confusion_matrix <- table(test_data$type, predicted_classes)
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
f1_score <- 2 * precision * recall / (precision + recall)

cat("Accuracy: ", accuracy, "\n")
cat("Precision: ", precision, "\n")
cat("Recall: ", recall, "\n")
cat("F1 Score: ", f1_score, "\n")


