data<-read.csv('data_final.csv', header = T)
data$Tumor.Class <- as.factor(data$Tumor.Class)

#Split the data into train and test data
smp_size <- floor(0.80 * nrow(data))
train_ind <- sample(seq_len(nrow(data)), size = smp_size)
train <- data[train_ind, ]
test <- data[-train_ind, ]
rf<-randomForest(Tumor.Class~.,data=train,importance=TRUE)
pred = predict(rf, newdata=test[,-1], type = "class")
#accuracy function with confusion matrix
accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}

#Confusion matrix
cm = table(test[,1], pred)
print(cm)
print(accuracy(cm))
var.imp1 <- data.frame(importance(rf, type=2))
var.imp1$Variables <- row.names(var.imp1)
varimp1 <- var.imp1[order(var.imp1$MeanDecreaseGini,decreasing = T),]
vars<-data.frame(varimp1$Variables)
vars<-transpose(vars)
vars<-vars[,c(1:1100)]
new_data<-data
new_data<-new_data[,(names(new_data) %in% vars)]
new_data<-cbind(data[,1],new_data)
setnames(new_data, "data[, 1]", "Tumor.Class")

new_data$Tumor.Class <- as.factor(new_data$Tumor.Class)

#Split the data into train and test data
smp_size <- floor(0.80 * nrow(new_data))
train_ind <- sample(seq_len(nrow(new_data)), size = smp_size)
train1 <- new_data[train_ind, ]
test1 <- new_data[-train_ind, ]

#Run random forest model and predict tumor class on test data
rf1 <- randomForest(Tumor.Class ~ .,data = train1)
pred1 = predict(rf1, newdata=test1[,-1], type = "class")


#accuracy function with confusion matrix
accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}

#Confusion matrix
cm1 = table(test1[,1], pred1)
print(cm1)
print(accuracy(cm1))

precision <- diag(cm1) / colSums(cm1)
recall <- diag(cm1) / rowSums(cm1)
f1 <-  ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
misclassification<-1-sum(diag(cm1))/sum(cm1)

#Predicting probabilities of test data for each class
predictions <- as.data.frame(predict(rf1, test1, type = "prob"))
#Taking maximum probability for the sample and assign to predict
predictions$predict <- names(predictions)[1:3][apply(predictions[,1:3], 1, which.max)]
#Actual tumor class
predictions$observed <- test1$Tumor.Class
head(predictions)

#Plotting ROC curve
roc.Ta <- roc(ifelse(predictions$observed=="Ta", "Ta", "non-Ta"), as.numeric(predictions$Ta))
roc.T1 <- roc(ifelse(predictions$observed=="T1", "T1", "non-T1"), as.numeric(predictions$T1))
roc.T2 <- roc(ifelse(predictions$observed=="T2+", "T2+", "non-T2+"), as.numeric(predictions$`T2+`))
plot(roc.Ta, col = "green",lty=1,lwd=4)
lines(roc.T2,col="red",lty=3,lwd=4)
lines(roc.T1, col = "blue",lty=2,lwd=4)
legend(x = "topright",legend = c("Ta", "T1","T2+"),lty = c(1, 2,3),col = c("green","blue","red"), lwd = 4)     

print(precision)
print(recall)
print(f1)
print(misclassification)
nrFolds <- 10
# generate array containing fold-number for each sample (row)
folds <- rep_len(1:nrFolds, nrow(new_data))
average_accuracy <- 0
# actual cross validation
for(k in 1:nrFolds) {
  # actual split of the data
  fold <- which(folds == k)
  train <- new_data[-fold,]
  test <- new_data[fold,]
  rf <- randomForest(Tumor.Class ~ .,data = train)
  pred = predict(rf, newdata=test[,-1], type = "class")
  cm = table(test[,1], pred)
  print(accuracy(cm))
  average_accuracy<- average_accuracy+accuracy(cm)
}
print(average_accuracy/1000)

bayes<-naiveBayes(Tumor.Class ~ .,data = train1)
pred2 = predict(bayes, newdata=test1[,-1], type = "class")
mean(pred2==new_data$Tumor.Class)

#accuracy function with confusion matrix
accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}

#Confusion matrix
cm2 = table(test1[,1], pred2)
print(cm2)
print(accuracy(cm2))

precision <- diag(cm2) / colSums(cm2)
recall <- diag(cm2) / rowSums(cm2)
f1 <-  ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
misclassification<-1-sum(diag(cm2))/sum(cm2)

#Predicting probabilities of test data for each class
predictions <- as.data.frame(predict(bayes, test1, type = "raw"))
#Taking maximum probability for the sample and assign to predict
predictions$predict <- names(predictions)[1:3][apply(predictions[,1:3], 1, which.max)]
#Actual tumor class
predictions$observed <- test1$Tumor.Class
head(predictions)

#Plotting ROC curve
roc.Ta <- roc(ifelse(predictions$observed=="Ta", "Ta", "non-Ta"), as.numeric(predictions$Ta))
roc.T1 <- roc(ifelse(predictions$observed=="T1", "T1", "non-T1"), as.numeric(predictions$T1))
roc.T2 <- roc(ifelse(predictions$observed=="T2+", "T2+", "non-T2+"), as.numeric(predictions$`T2+`))
plot(roc.Ta, col = "green",lty=1,lwd=4)
lines(roc.T2,col="red",lty=3,lwd=4)
lines(roc.T1, col = "blue",lty=2,lwd=4)
legend(x = "topright",legend = c("Ta", "T1","T2+"),lty = c(1, 2,3),col = c("green","blue","red"), lwd = 4)     

print(precision)
print(recall)
print(f1)
print(misclassification)

nrFolds <- 10
# generate array containing fold-number for each sample (row)
folds <- rep_len(1:nrFolds, nrow(new_data))
average_accuracy <- 0
# actual cross validation
for(k in 1:nrFolds) {
  # actual split of the data
  fold <- which(folds == k)
  train <- new_data[-fold,]
  test <- new_data[fold,]
  bayes2 <- naiveBayes(Tumor.Class ~ .,data = train)
  pred = predict(bayes2, newdata=test[,-1], type = "class")
  cm = table(test[,1], pred)
  print(accuracy(cm))
  average_accuracy<- average_accuracy+accuracy(cm)
}
print(average_accuracy/1000)
