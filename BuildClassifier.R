library(data.table)
library(gdata)
library(randomForest)
library(verification)
library(ggplot2)
library(pROC)

#### load data ####
alldata = read.csv("alldata.csv", header = TRUE)
alldata.train = read.csv("alldata.train.csv", header = TRUE)
alldata.test = read.csv("alldata.test.csv", header = TRUE)

summary(alldata.train$element_type)
print(table(alldata.train$type))
summary(alldata.test$element_type)
print(table(alldata.test$type))
summary(alldata$element_type)
print(table(alldata$type))

############################################################################################################################

#### Cross-Valdidation for feature selection ####
# cros validation
set.seed(415)
rfcv.res = rfcv(alldata.train[,c(18, 20:29, 33:44, 48, 50:56, 60:67)], as.factor(alldata.train[,70]))

# plot
plot(rfcv.res$n.var, rfcv.res$error.cv, log="x",type="o", lwd=2,
     xlab="Number of Variables", ylab="Error Rate")
lines(rfcv.res$n.var, rfcv.res$error.cv)
title(main="Estimated Error Rate")

print(rfcv.res$error.cv)

#### Fit ####
#model fit with all predictors
set.seed(415)
fit.all <- randomForest(as.factor(type) ~ len34 + G34 + G34.len + Gfold34 + Gfold34.len + freq_MFE34 + diversity34 + avg_prob34 + avg_sl34_prob + mul_prob34 + mul_sl34_prob + len23 + end23 + G23 + G23.len + Gfold23 + Gfold23.len + freq_MFE23 + diversity23 + avg_prob23 + avg_sl23_prob + mul_prob23 + mul_sl23_prob + len + GP1 + GP1.len + overlap + ddG + fcG + ddG.SL + fcG.SL + Gfold23.G + Gfold23.G.1 + Gfold34.G + Gfold34.G.1 + dprob + dsl_prob + d_mul_prob + dsl_mul_prob, 
                        data=alldata.train)

#predictors with Gini importance value
import <- importance(fit.all, sort = TRUE)

#reordering of predictors by importance
import <- import[order(import, decreasing=TRUE),,drop = FALSE]

#15 predictors listed by Gini importance
import[1:15,,drop = FALSE]

#subset of 15 most important predictors
alldata.train.15 <- alldata.train[,rownames(import)[1:15]]

#random forest model fit with 24 most important predictors
set.seed(415)
fit.15 <- randomForest(as.factor(alldata.train$type) ~.,
                                data=alldata.train.15, importance=TRUE, proximity=TRUE)

#### Check ####
# Err rate plot (vs num of trees)
qplot(y=fit.15$err.rate[,1]) #+ coord_cartesian(xlim = c(1000,1500),ylim= c(0.16,0.18))

# importance plot
varImpPlot(fit.15)
importance(fit.15)

MDSplot(fit.15, as.factor(alldata.train$type), palette = c("red", "green"))
plot(outlier(fit.15), type="h", col=c("red", "blue")[as.factor(alldata.train$type)])
plot(fit.15, log="y")

fit.15$confusion

#### save ####
fit <- fit.15
save(fit, file = "Classifier.RData")
#load("Classifier.RData") => Change accordingly in the RFprediction.R file!!!

############################################################################################################################

#### k-fold cv ####
alldata.train.15$type = alldata.train$type
suff <- alldata.train.15[sample(nrow(alldata.train.15)), ]

k=3

n=floor(nrow(alldata.train.15)/k)
err.vect=0

for(i  in 1:k)
{
  s1 = ((i-1)* n+1)
  s2 = (i*n)
  subset = s1:s2
  
  cv.train = suff[-subset,]
  cv.test = suff[subset,]
  
  fit.CV <- randomForest(as.factor(type) ~., data=cv.train)

  Predictioncv = predict(fit.CV, cv.test, type="prob")
  
  cv.test$response <- ifelse(cv.test$type=="positive", 1, 0)
  
  err.vect[i] = roc.area(cv.test$response, Predictioncv[,2])$A
  roc.plot(cv.test$response, Predictioncv[,2])
  print(paste("AUC for fold", i, ":", err.vect[i]))
}
print(paste("Average AUC:", mean(err.vect)))
  

#### colored k-fold cv ####
alldata.train.15$type = alldata.train$type
suff <- alldata.train.15[sample(nrow(alldata.train.15)), ]

k=10
n=floor(nrow(alldata.train.15)/k)

err.vect=0

i=1
s1 = ((i-1)* n+1)
s2 = (i*n)
subset = s1:s2

cv.train = suff[-subset,]
cv.test = suff[subset,]

fit.CV <- randomForest(as.factor(type) ~., data=cv.train)

Predictioncv = predict(fit.CV, cv.test, type="prob")
cv.test$response <- ifelse(cv.test$type=="positive", 1, 0)

curr_roc = roc(cv.test$response, Predictioncv[,2])#, smooth = TRUE)
err.vect[i] = curr_roc$auc
plot.roc(curr_roc, legacy.axes=TRUE, print.auc = TRUE, print.auc.x=0.25)

for(i  in 2:k)
{
  s1 = ((i-1)* n+1)
  s2 = (i*n)
  subset = s1:s2
  
  cv.train = suff[-subset,]
  cv.test = suff[subset,]
  
  fit.CV <- randomForest(as.factor(type) ~., data=cv.train)
  
  Predictioncv = predict(fit.CV, cv.test, type="prob")
  cv.test$response <- ifelse(cv.test$type=="positive", 1, 0)
  
  curr_roc = roc(cv.test$response, Predictioncv[,2])#, smooth = TRUE)
  err.vect[i] = curr_roc$auc
  plot.roc(curr_roc, legacy.axes=TRUE, print.auc = TRUE, add = TRUE, col = i, print.auc.y=0.5-(i-1)*0.05, print.auc.x=0.25)
}
print(paste("Average AUC:", mean(err.vect)))
print(paste("STD", sd(err.vect)))

############################################################################################################################

#### predict test data ####
alldata.test <- alldata.test[sample(nrow(alldata.test)), ]         #shuffle
Prediction.test <- predict(fit.15, alldata.test, type="prob")

#### confusion mat
Prediction.test.res = as.factor(ifelse(Prediction.test[,2] >= 0.5, "positive", "negative"))
tab = table(Prediction.test.res, alldata.test$type)
tab
cat(paste(paste("sensitivity (TPR)", tab[2,2]/(tab[2,2]+tab[1,2])*100), 
          paste("specificity (TNR)", tab[1,1]/(tab[1,1]+tab[2,1])*100), 
          sep="\n"))

#plot
qplot(y = Prediction.test[,2], color = alldata.test$element_type) + 
  geom_point(size=6, aes(shape=alldata.test$type)) +
  scale_color_manual(values=c("chartreuse4", "chartreuse3", "chartreuse2", 
                              "tomato4", "darkred", "firebrick", "tomato3", "tomato2", "tomato1",
                              "sienna1",
                              "firebrick3")) +
  xlab("") + ylab("prediction") + theme(axis.title.y = element_text(size = rel(1.5))) +
  guides(color=guide_legend(title="Element Type"), shape=guide_legend(title=NULL)) +
  geom_hline(aes(yintercept=0.5))

# ROC
alldata.test$response <- ifelse(alldata.test$type=="positive", 1, 0)
roc.plot(alldata.test$response, Prediction.test[,2])
paste("AUC (test data)", roc.area(alldata.test$response, Prediction.test[,2])$A)

boxdata <- data.frame(Type = alldata.test$type, Prediction = Prediction.test[,2])
ggplot(data = boxdata, aes(Type, Prediction)) + geom_boxplot(outlier.shape = 1, notch = TRUE)

############################################################################################################################
