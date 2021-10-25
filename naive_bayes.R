#Loading the data with Tumor sample class and gene expression data.
data<-read.csv('data_final.csv', header = T)

#Run randomforest classification algorithm on the original data.
library(ROCR)
library(dplyr)
library(caret)
library(mlbench)
library(e1071)
library(pROC)
library(corrplot)

#Make Tumor class column as factor
data$Tumor.Class <- as.factor(data$Tumor.Class)

#Split the data into train and test data
smp_size <- floor(0.80 * nrow(data))
train_ind <- sample(seq_len(nrow(data)), size = smp_size)
train_alldata <- data[train_ind, ]
test_alldata <- data[-train_ind, ]

#Run random forest model and predict tumor class on test data
bayes_alldata <- naiveBayes(Tumor.Class ~ .,data = train_alldata)
pred_alldata = predict(bayes_alldata, newdata=test_alldata[,-1], type = "class")
mean(pred_alldata==data$Tumor.Class)

#accuracy function with confusion matrix
accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}

#Confusion matrix
cm_alldata = table(test_alldata[,1], pred_alldata)
print(cm_alldata)
print(accuracy(cm_alldata))

precision <- diag(cm_alldata) / colSums(cm_alldata)
recall <- diag(cm_alldata) / rowSums(cm_alldata)
f1 <-  ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
misclassification<-1-sum(diag(cm_alldata))/sum(cm_alldata)

print(precision)
print(recall)
print(f1)
print(misclassification)


#Predicting probabilities of test data for each class
predictions <- as.data.frame(predict(bayes_alldata, test_alldata[,-1], type = "raw"))
#Taking maximum probability for the sample and assign to predict
predictions$predict <- names(predictions)[1:3][apply(predictions[,1:3], 1, which.max)]
#Actual tumor class
predictions$observed <- test_alldata$Tumor.Class
head(predictions)

#Plotting ROC curve
roc.Ta <- roc(ifelse(predictions$observed=="Ta", "Ta", "non-Ta"), as.numeric(predictions$Ta))
roc.T1 <- roc(ifelse(predictions$observed=="T1", "T1", "non-T1"), as.numeric(predictions$T1))
roc.T2 <- roc(ifelse(predictions$observed=="T2+", "T2+", "non-T2+"), as.numeric(predictions$`T2+`))
plot(roc.Ta, col = "green",lty=1,lwd=4)
lines(roc.T2,col="red",lty=3,lwd=4)
lines(roc.T1, col = "blue",lty=2,lwd=4)
legend(x = "topright",legend = c("Ta", "T1","T2+"),lty = c(1, 2,3),col = c("green","blue","red"), lwd = 4)     


#10-fold cross validation with random forest
nrFolds <- 10
average_accuracy <- 0
# generate array containing fold-number for each sample (row)
folds <- rep_len(1:nrFolds, nrow(data))
# actual cross validation
for(k in 1:nrFolds) {
  # actual split of the data
  fold <- which(folds == k)
  train <- data[-fold,]
  test <- data[fold,]
  nb <- naiveBayes(Tumor.Class ~ .,data = train)
  pred = predict(nb, newdata=test[,-1], type = "class")
  cm = table(test[,1], pred)
  print(accuracy(cm))
  average_accuracy <- average_accuracy+accuracy(cm)
  
}
print(average_accuracy/1000)

#Feature selection by removing highly correlated genes.
correlationMatrix <- cor(data[,2:7130])
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75)
reduced_data<-data[,-highlyCorrelated]
reduced_data<-cbind(data[,1],reduced_data)
reduced_data<-reduced_data[,-1]

#Make Tumor class column as factor
reduced_data$Tumor.Class <- as.factor(reduced_data$Tumor.Class)

#Split the data into train and test data
smp_size <- floor(0.80 * nrow(reduced_data))
train_ind <- sample(seq_len(nrow(reduced_data)), size = smp_size)
train_reduceddata <- reduced_data[train_ind, ]
test_reduceddata <- reduced_data[-train_ind, ]

#Run random forest model and predict tumor class on test data
bayes_reduceddata <- naiveBayes(Tumor.Class ~ .,data = train_reduceddata)
pred_reduceddata = predict(bayes_reduceddata, newdata=test_reduceddata[,-1], type = "class")

#Confusion matrix
cm_reduceddata = table(test_reduceddata[,1], pred_reduceddata)
print(cm_reduceddata)
print(accuracy(cm_reduceddata))

precision <- diag(cm_reduceddata) / col(cm_reduceddata)
recall <- diag(cm_reduceddata) / rowSums(cm_reduceddata)
f1 <-  ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
misclassification<-1-sum(diag(cm_reduceddata))/sum(cm_reduceddata)

print(precision)
print(recall)
print(f1)
print(misclassification)

#Predicting probabilities of test data for each class
predictions_r <- as.data.frame(predict(bayes_reduceddata, test_reduceddata, type = "raw"))
#Taking maximum probability for the sample and assign to predict
predictions_r$predict <- names(predictions_r)[1:3][apply(predictions_r[,1:3], 1, which.max)]
#Actual tumor class
predictions_r$observed <- test_reduceddata$Tumor.Class
head(predictions_r)

#Plotting ROC curve
roc.Ta <- roc(ifelse(predictions_r$observed=="Ta", "Ta", "non-Ta"), as.numeric(predictions_r$Ta))
roc.T1 <- roc(ifelse(predictions_r$observed=="T1", "T1", "non-T1"), as.numeric(predictions_r$T1))
roc.T2 <- roc(ifelse(predictions_r$observed=="T2+", "T2+", "non-T2+"), as.numeric(predictions_r$`T2+`))
plot(roc.Ta, col = "green",lty=1,lwd=4)
lines(roc.T2,col="red",lty=3,lwd=4)
lines(roc.T1, col = "blue",lty=2,lwd=4)
legend(x = "topright",legend = c("Ta", "T1","T2+"),lty = c(1, 2,3),col = c("green","blue","red"), lwd = 4)



#10-fold cross validation with random forest
nrFolds <- 10
# generate array containing fold-number for each sample (row)
folds <- rep_len(1:nrFolds, nrow(reduced_data))
average_accuracy<-0
# actual cross validation
for(k in 1:nrFolds) {
  # actual split of the data
  fold <- which(folds == k)
  train <- reduced_data[-fold,]
  test <- reduced_data[fold,]
  nb <- naiveBayes(Tumor.Class ~ .,data = train)
  pred = predict(nb, newdata=test[,-1], type = "class")
  cm = table(test[,1], pred)
  print(accuracy(cm))
  average_accuracy<-average_accuracy+accuracy(cm)
}
print(average_accuracy/1000)
