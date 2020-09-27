##=========================================================*
#**********************************************************#
#By: Emine Ozsahin
#Date: Feb 27, 2020 - Feb 28, 2020

#**********************************************************#
##=========================================================*

#loading the libraries

library(glmnet)
library(caret)
library(ggplot2)
library(tidyverse)
library(rpart)
library(class)
library(caret)
library(plotly)
library(lars)
library(GGally)
library(boot)
#install.packages("ggvis")
library(ggvis)
library(pROC)

setwd("/Users/emineozsahin/Documents/R_worksheets_winter2020/BINF6970_Stat")

##=========================================================*
##----------------------------------------------------*
## 	Problem 1 - i)----

# Examine the relationship between the response variable, mpg, and the predictor variable x.
##----------------------------------------------------*
##=========================================================*

Auto <- read.csv("auto_exam.csv", header=TRUE)
str(Auto)

par(mfrow=c(1,2))
hist(Auto$mpg)
hist(Auto$x)
par(mfrow=c(1,1))

# scatter plot
attach (Auto)

Auto %>% ggvis(~x, ~mpg, fill = ~x) %>% layer_points()

cor(Auto)

#ggpairs(Auto)
#pairs(Auto)
#boxplot(Auto)
#plot(x , mpg)

##=========================================================*
##----------------------------------------------------*
## 	Problem 1 - ii) ----

# Consider 9 linear regression models for mpg with increasing polynomial degrees of x. That is, Model-1 has polynomial degree 1; Model-2 has polynomial degrees 2, …, and Model-9 has polynomial degrees 9.

# Perform 5-fold cross-validation to choose the best parsimonious model among these 9 models. Repeat your cross-validation 100 times and use the overall average CV-error in your analysis.
##----------------------------------------------------*
##=========================================================*

# I generated 9 models. Each model has different lenght of polynomials. And I put them in a list called poly_list. Model-1 has one polynomial degree, Model-2 has two polynomial degree and so on up to 9.
# I also made a data frame to see and compare their AIC scores. 

poly_list <- list()
poly_AIC <- data.frame()

for (i in 1:9) {Model <- paste("Model", i, sep="-")
                poly_list[[Model]] <- glm(mpg ~ poly(x, i), data=Auto, family = gaussian)
                temp <- data.frame(Model= Model, AIC = poly_list[[i]]$aic)
                poly_AIC <- rbind(poly_AIC, temp)}

length(poly_list)
poly_list[[2]]

poly_AIC %>% ggvis(~ Model, ~ AIC, fill = ~ AIC, size := 300, opacity := 0.75) %>% layer_points() 

# I made a data frame of the results of cross validations of each Model found in poly_list. The first column indicate the Model name as explained above.

cv100_results <- data.frame()
cnt <- 0

repeat { 
  for (i in 1:length(poly_list)) { Model <- paste("Model", i, sep="-")
                                          cv.err.5 <- cv.glm(Auto, poly_list[[i]], K = 5)$delta
                                          temp<- data.frame (Model = Model, raw_cv_est_of_predic_error = cv.err.5[1])
                                          cv100_results <- rbind(cv100_results, temp)
                                          }

                                          cnt <- cnt+1
  if(cnt == 100) {
    break 
    }
  }

cv100_results[1:9,]
cv100_results[10:18,]

?cv.glm()

length(which(cv100_results$Model=="Model-1")) # 100
length(cv100_results[which(cv100_results$Model=="Model-1"), 2]) #100

cv100_results[which(cv100_results$Model=="Model-1"), 2]
mean(cv100_results[which(cv100_results$Model=="Model-1"), 2]) # 18.53931

# I put the avarage CV-ERROR for 100 times cross-validation of each models to a data frame to be able to compare easily. Following codes generate multiple data for each model for 100 times in the data frame because I did not add if statement. But later I removed the duplicated data. I know this is not the optimal solution especially for big data. I`ll try "if" statement in the for loop later. But for now it is okay because current data size is small.

df_mean_cv_error <- data.frame()

for (i in cv100_results$Model) {mean <- mean(cv100_results[which(cv100_results$Model==i), 2])
                                temp <- data.frame(Model = i, Mean_CV_ERROR = mean)
                                df_mean_cv_error <- rbind (df_mean_cv_error, temp) }

df_mean_cv_error <- df_mean_cv_error[!duplicated(df_mean_cv_error), ]

df_mean_cv_error


# Model-2 has the least error.


 ##=========================================================*
##----------------------------------------------------*
## 	Problem 1 - iii) ----

# Visualize how the averaged cross-validation error varies as a function of increasing polynomial degrees. Comment on the overall shape of the CV-error curve and identify the underlying factors driving this specific shape.
##----------------------------------------------------*
##=========================================================*
attach(df_mean_cv_error)

df_mean_cv_error %>% ggvis(~ Model, ~ Mean_CV_ERROR, fill = ~ Mean_CV_ERROR, size := 300, opacity := 0.75) %>% layer_points() 

# polynomial reduced the cross validation error.  But high polynomial degree may overfit to the data.

##=========================================================*
##----------------------------------------------------*
## 	Problem 1 - iv) ----

# Reinforce your argument above by visualizing the best and worst model fits to the auto_exam dataset. Comment on the goodness of model fits.
##----------------------------------------------------*
##=========================================================*

#First model is the worst which has one degree of polynomial meaning that there is only one variable which is x  and 

#second model is the best which has two degree of polynomials meaning that there are two variables x and x^2 in terms of CV-Error

poly_AIC

Model_worst <- poly_list[[1]]

Model_best <- poly_list[[2]]

par(mfrow=c(1,2))
plot(Model_worst, which = 1, main="The Worst Model")
plot(Model_best, which = 1, main = "The Best Model")

plot(Model_worst, which = 2, main="The Worst Model")
plot(Model_best, which = 2, main = "The Best Model")

plot(Model_worst, which = 3, main="The Worst Model")
plot(Model_best, which = 3, main = "The Best Model")

plot(Model_worst, which = 4, main="The Worst Model")
plot(Model_best, which = 4, main = "The Best Model")

plot(Model_worst, which = 5, main="The Worst Model")
plot(Model_best, which = 5, main = "The Best Model")

par(mfrow=c(1,1))

# missing part that I forgot to add the midterm
best <- glm(mpg ~ poly(x, 2), data=Auto, family = gaussian)

worst <- glm(mpg ~ poly(x, 1), data=Auto, family = gaussian)

ggplot(Auto, aes(x=x, y=mpg)) + geom_point(cex=2) + geom_smooth(method=lm, mapping=aes(y=predict(worst, Auto)))

ggplot(Auto, aes(x=x, y=mpg)) + geom_point(cex=2) + geom_smooth(method=lm, mapping=aes(y=predict(best, Auto)))

# polynomial regression of mpg and the quadratic effect of x^2 produced better model than the linear relation of mpg and x 

##=========================================================*
##----------------------------------------------------*
## 	Problem 2 Data Exploration, Wrangling and Splitting ----
##----------------------------------------------------*
##=========================================================*

# loding the data
diabets <- read.csv("diabetes.csv")

# Data exploration
str(diabets)
diabets$positive <- as.factor(ifelse(test = diabets$positive == 0, yes = "No", no = "Yes"))
str(diabets)
summary(diabets)
head(diabets)

# checked if there is any missing values
colSums(is.na(diabets))

# Checked the data in terms of the variable distributions based on positive
attach(diabets)
colnames(diabets)
xtabs(~ positive + age, data = diabets)
xtabs(~ positive + preg, data = diabets)
xtabs(~ positive + plas, data = diabets)
xtabs(~ positive + pres, data = diabets)
xtabs(~ positive + skin, data = diabets)
xtabs(~ positive + insu, data = diabets)
xtabs(~ positive + mass, data = diabets)

#Split the data for validation and training
set.seed(1245)
test.index <- sample.int(nrow(diabets), round(nrow(diabets)*.15),replace=FALSE)

head(diabets)

# Data converted to the format for glmnet() function
X <- model.matrix(~., data=diabets)[,-1]
head(X)
nrow(X) #768
str(diabets)
str(X)

X_training <- X [-test.index, -9]
dim(X_training) #653 8
head(X_training)

X_validation <- X [test.index, -9]
dim(X_validation) #115 8
head(X_validation)
  

class(X[-test.index, 9])
y_trainig <- as.factor(X[-test.index, 9])
length(y_trainig) #653
class(y_trainig)

y_validation <- as.factor(X[test.index, 9])
length(y_validation) #115

# Generated a list of 10 folds for cross validation to use for all models (Logistic Regression, KNN, Classification Trees)

length(diabets$positive[-test.index]) #653

#set.seed(1245)
set.seed(1244)
folds_for_caret <-  createFolds(diabets$positive[-test.index],k=10)

length(folds_for_caret[[1]]) #66

cntr <- trainControl(method = "cv", index=folds_for_caret)

##=========================================================*
##----------------------------------------------------*
## Problem 2 
## Regularized Logistic Regression with Elastic Net ----

# For the elastic net training, use seq(0.1, 0.9, 0.01) as the tuning grid for alpha. 
# Once you find the best lambda-alpha pair, use training data to find an optimal classification threshold as well.
##----------------------------------------------------*
##=========================================================*

set.seed(1245)
cv_enet <- train(X_training, y_trainig, method = "glmnet",  trControl = cntr, family="binomial", type.measure = "auc")

cv_enet$results

max_pair <- which.max(cv_enet$results[,3]) 

cv_enet$results[max_pair,]  # The best Accuracy is 0.7592418

#The alpha value of the model with the best accuracy value is 0.55.
cv_enet$results[max_pair,1] 

#The lambda value of the model with the best accuracy value is 0.0448592.
cv_enet$results[max_pair,2] 

# I generated a grid for different value of alpha and lamdas
cv_enet$results
set.seed(1245)
tuneGrid <- expand.grid(.alpha = seq(0.1, 0.9, 0.01),.lambda = seq(0.004, 0.05, 0.005))
unique(tuneGrid[,2])
cv_enet$results
length(tuneGrid[,1]) #810

set.seed(1245)
enet_tune <- train(
  X_training, y_trainig, method = "glmnet",
  trControl = cntr,
  tuneGrid = tuneGrid)

head(enet_tune$results)

#++++++++++++++
# I used one of alpha and lamda pair to build a model and make prediction to check if my codes to build all the Models in a for loop will work and predict correctly. Because I am going to use all alpha and lambda pairs to build Models for each pair. 

set.seed(1245)

head(X_training)
head(y_trainig)
control <- glmnet(X_training, y_trainig, family="binomial", alpha = enet_tune$results[1,1], lambda = enet_tune$results [1,2], type.measure = "auc")

pr_control <- predict(control, newx = X_training, type = "response")[,1]
pr_control[1]

#+++++++++++++++

# I used all alpha and lambda pairs to fit a Model and I put all Models (each represent a pair of alpha and lambda) in a list. I had 810 Models. 
set.seed(1245)
list.of.enet.fits <- list()
for (i in 1:length(enet_tune$results$alpha)) {fit <- paste("fit", i, sep="-")
list.of.enet.fits[[fit]] <- glmnet(X_training,
                                   y_trainig,
                                   family="binomial",
                                   alpha = enet_tune$results[i,1],
                                   lambda = enet_tune$results[i,2],
                                   type.measure = "auc")}

# I checked if the codes worked properly. 
length(enet_tune$results$alpha) #810
length(list.of.enet.fits) #810
names(list.of.enet.fits)
list.of.enet.fits$`fit-1`

# I used the Models in the list to predict the training response (in this data it called positive) by using the variables in training data set. 

prds_train_elastic <- lapply(list.of.enet.fits, predict, newx=X_training, type = "response") 

plot(prds_train_elastic$`fit-84`, y_trainig, main = "fit-84")

# I checked the predictions. My expect is to have a list of list. In the first list there must be 810 lists each represent a list of predicted y values by each Model. Therefore each list of 810 should have 653 predicted y values. 
dim(X_training) #653 8
length(y_trainig) #653 
length(prds_train_elastic) #810
length(prds_train_elastic[[1]]) #653

# I checked if the predictions are the same when I used one of Model separetly to fit the model and predict the response and they are same. Calculations are correct.
prds_train_elastic[[1]] == pr_control #TRUE
prds_train_elastic[[1]][1] 
pr_control[1] 

class(prds_train_elastic[[1]]) #matrix
class(pr_control) # numeric

# I generate a list of accuracy of the Models 1:810 and I also made a dataframe for Model names and accuracy scores to compare the models
auc_list_train_enet <- list()
results_enet_train <- data.frame()

for (i in 1:length(prds_train_elastic)) { fit.name <- paste("fit", i, sep="-")
  auc_list_train_enet[[fit.name]] <- roc(y_trainig, as.numeric(prds_train_elastic[[i]])) 
  temp <- data.frame(fit.name = fit.name, auc = auc_list_train_enet[[fit.name]]$auc)
  results_enet_train <- rbind(results_enet_train, temp)
  }

length(auc_list_train_enet) #810
auc_list_train_enet[1]
head(results_enet_train)

plot(results_enet_train$fit.name, results_enet_train$auc, xlab="810 Models", ylab="AUC")

plot(auc_list_train_enet$`fit-84`, which=1, main="Best Model (fit-84) on Training Data Set")

#tune grid

l2 <- enet_tune$results[,"lambda"]
a2 <- enet_tune$results[,"alpha"]
r1 <- results_enet_train$auc
plot_ly(x= a2, y= l2, z=r1, type="scatter3d", mode="markers", color=l2)

# I checked if the accuracies in the data frame are correct. I used the the First Model which I also fit and predict seperately as a control to confirm the calculated accuracies are correct. 
(roc(y_trainig,pr_control))$auc
auc_list_train_enet$`fit-1`$auc
results_enet_train[1, ]

# I found the highest accuracy in the data frame which is located at the 84th row. Therefore I looked at the columns of 84th row  
which.max(results_enet_train$auc) #84
results_enet_train[84, ] # fit-84, AUC = 0.8387446 # accuracy is better now

list.of.enet.fits[84] # lambda 0.019

# I made a data frame to select the best sensitivity and specificity based on train dataset of the best model based on auc score.
snsp_train_fit84 <- data.frame(sensitivity= auc_list_train_enet$`fit-84`$sensitivities, specificity = auc_list_train_enet$`fit-84`$specificities, cutoff = auc_list_train_enet$`fit-84`$thresholds)

dim(snsp_train_fit84) # 654 3

indx <- which.max(apply(snsp_train_fit84[,1:2],1,min))

plot(1-snsp_train_fit84[,1:2], main="fit84")

# Finally I selected the best cutoff of the Model fit742 based on the best specificity and sensitivity
snsp_train_fit84[indx, ]
cutoff <- snsp_train_fit84[indx, ][1,3]

## A function to compute Sensitivity and Specificity
sn.sp <- function(mat){
  sn <- mat[2,2]/sum(mat[2,])
  sp <- mat[1,1]/sum(mat[1,])
  return(unlist(list(sensitivity=sn, specificity=sp)))
}

table(y=y_trainig,yhat=as.numeric(prds_train_elastic[[1]]>cutoff))

sn.sp(table(y=y_trainig,yhat=as.numeric(prds_train_elastic[[1]]>cutoff)))

# But I tried all the theresholds belonged to the fit84 to find the best thereshold onthe predicted training respond values by fit84

names(prds_train_elastic[84])

snsp_df <- data.frame()
for (i in (auc_list_train_enet$`fit-84`$thresholds[2:653])) { 
    snsp <- sn.sp(table(y=y_trainig, yhat=as.numeric(prds_train_elastic[[84]]>i)))
    temp<- data.frame(sensitivity = snsp[1], specificity = snsp[2], cutoff = i)
    snsp_df <-rbind(snsp_df,temp)
  }

length(auc_list_train_enet$`fit-84`$thresholds[2:653])

head(snsp_df)

plot(snsp_df$sensitivity, snsp_df$cutoff, main = 'fit84')

#The best cutoff which produced the highest sensitivity and specificity is 0.3270323
head(snsp_df[, 1:2])
which.max(apply(snsp_df[, 1:2], 1, min)) #380

# sn and sp from the Model-84 (fit-84). This was calculated by the roc function
snsp_train_fit84[indx,]
# sensitivity specificity    cutoff
# 381   0.7633929   0.7622378 0.3270323

# sn and sp from confusion matrixes produced by the predictions of response variable this is calculated by the sn.sp function 
snsp_df[380, ]
# sensitivity specificity    cutoff
# sensitivity379   0.7633929   0.7622378 0.3270323

# I assigned the best thereshould as cutoff 
cutoff <- snsp_df[380, ][[3]]

# Confusion matrix of traning set 
table(y=y_trainig, yhat=as.numeric(prds_train_elastic[[84]]>cutoff))

# # Those two best cutoff are similar which means I guess roc and snsp fuction calculates the same thing.I investigated this. And yes they produced the same results.
sn.sp(table(y=y_trainig, yhat=as.numeric(prds_train_elastic[[84]]>cutoff)))
roc_84 <- roc(y_trainig, as.numeric(prds_train_elastic[[84]]>cutoff))
roc_84$sensitivities; roc_84$specificities

l2 <- enet_tune$results[,"lambda"]
a2 <- enet_tune$results[,"alpha"]
r1 <- results_enet_train$auc
plot_ly(x= a2, y= l2, z=r1, type="scatter3d", mode="markers", color=l2)

##+++++++++++++++++++++++++++++++++++
### Prediction on the Validation Set for the best Model which is fit84 from training test 
### Prediction Accuracy on validation data set
prds_test_logistic <- predict(list.of.enet.fits$`fit-84`,newx = X_validation, type = "response")[,1]

length(prds_test_logistic)
length(y_validation)
class(y_validation)

y <- as.numeric(y_validation)

snsp_validation_logistic <- sn.sp(table(y=y,yhat=as.numeric(prds_test_logistic>cutoff)))
snsp_validation_logistic
# sensitivity specificity 
# 0.7727273   0.8028169 

# AUC for the validation set
auc_test_logistic <- roc(y_validation, prds_test_logistic)
auc_test_logistic$auc # 0.8473

plot(auc_test_logistic, main="Logistic Regression on Validation Data Set")

## COnfusion Matrix
table(y=y,yhat=as.numeric(prds_test_logistic>cutoff))

##=========================================================*
##----------------------------------------------------*
## 	Problem 2 
##  KNN Classification ----
## Once you find the best lambda-alpha pair, use training data to find an optimal classification threshold as well.
##----------------------------------------------------*
##=========================================================*
## cntr variable that I used for trControl is the one that I assigned for cross validation with the certain folds generated by the caret package (Problem 2 Data Exploration, Wrangling and Splitting) to be ably to compare the models accurately.

set.seed(1244)
suppressWarnings(cvknn <- train(positive ~ ., data = diabets[-test.index,], method = "knn",
                                trControl = cntr, preProcess = c("center","scale"), 
                                tuneGrid = data.frame(k = 1:20)))

plot(cvknn)

cvresults <- cvknn$results
head(cvresults)

cvresults[which.max(cvresults[,2]),]
## The best fitting model has k=13 neighbors 

trainPrdsKNN <- predict(cvknn,newdata = diabets[-test.index,])
plot(roc(y_trainig,as.numeric(trainPrdsKNN)), main="KNN on Training Data Set")

#Prediction on the Validation Set
testPrdsKNN <- predict(cvknn,newdata = diabets[test.index,])

## COnfusion Matrix
table(diabets$positive[test.index], testPrdsKNN) ##  Do we have a useful model?

sn_knn <- sn.sp(table(diabets$positive[test.index], testPrdsKNN))[1]
sn_knn
sp_knn <- sn.sp(table(diabets$positive[test.index], testPrdsKNN))[2]
sp_knn 
## How about area under the ROC curve?
auc.test_knn <- roc(y_validation,as.numeric(testPrdsKNN))
plot(auc.test_knn, main="KNN on Test Data Set")
auc.test_knn$auc

##=========================================================*
##----------------------------------------------------*
## 	Problem 2 
##  Classification Tress ----
## Once you find the best lambda-alpha pair, use training data to find an optimal classification threshold as well.
##----------------------------------------------------*
##=========================================================*

## Let's use the caret package

set.seed(1244)

# Here I designed a new train control because the one that use for Logistic and KNN did not work here. But I used folds_for_caret variable which holds the folds that I used for previous methods in the train control therefore I used same folds as for previous methods.

train.control <- trainControl(method = "repeatedcv",
                              repeats = 3,## repeated three times # USE AUC
                              summaryFunction = twoClassSummary, 
                              classProbs = TRUE,
                              index=folds_for_caret)

system.time (rpartFit <- train(positive ~ ., data = diabets[-test.index,], 
                               method = "rpart2", 
                               tuneLength = 20,
                               trControl = train.control,
                               metric = "ROC"
))

rpartFit$results

fancyRpartPlot(rpartFit$finalModel)  

# I used the final model to predict the validation responce. 
testPrdsTree <- predict(rpartFit$finalModel, diabets[test.index,],type = "class")

#Confusion Matrix 
table(diabets$positive[test.index], testPrdsTree) 

hist(predict(rpartFit$finalModel, diabets[test.index, ],type = "prob")[,2], main = "Predicted Response Probabilities", xlab  ="Probabilities")

## The 0.5 threshold sounds reasonable?
sn_tree <- sn.sp(table(diabets$positive[test.index], testPrdsTree))[1]
sn_tree
sp_tree <- sn.sp(table(diabets$positive[test.index], testPrdsTree))[2]
sp_tree

sn_tree[[1]]

## How about area under the ROC curve?
auc.test_tree <- roc(diabets$positive[test.index],as.numeric(testPrdsTree))
auc.test_tree$auc
##=========================================================*
##----------------------------------------------------*
## 	Problem 2 Comparisons ----

# Report AUC, optimal sensitivity and specificity values under the three trained models for the validation data set (you should summaries them in a table).
##----------------------------------------------------*
##=========================================================*
#Confusion Matrixes of three models
table(diabets$positive[test.index], testPrdsKNN)
table(y=y,yhat=as.numeric(prds_test_logistic>cutoff))
table(diabets$positive[test.index], testPrdsTree) 

# A data frame to compare and plot three methods 
final_comparison <- data.frame("AUC" = c(auc_test_logistic$auc, auc.test_knn$auc, auc.test_tree$auc), 
                               "Sensitivity" = c(snsp_validation_logistic[1][[1]], sn_knn[[1]], sn_tree[[1]]), 
                               "Specificity" = c(snsp_validation_logistic[2][[1]], sp_knn[[1]], sp_tree[[1]]))

rownames(final_comparison) <- c("Logistic", "KNN", "Trees")

final_comparison
### Regularized Logistic Regression with elastic net ended up to be the best method based for this data set on AUC values and Sensitivity. 
##=========================================================*
##----------------------------------------------------*
## 	Problem 2 - VISUALIZATION  ----

# Plot validation set ROC curves for the three models.
##----------------------------------------------------*
##=========================================================*
snsp_test_logistic <- data.frame(sn=auc_test_logistic$sensitivities, sp=auc_test_logistic$specificities)

snsp_test_knn <- data.frame(sn=auc.test_knn$sensitivities, sp=auc.test_knn$specificities)

snsp_test_tree <- data.frame(sn=auc.test_tree$sensitivities, sp=auc.test_tree$specificities)

snsp_test_logistic$model <- factor("Logistic",c("Logistic","KNN", "Trees"))
snsp_test_knn$model <- factor("KNN",c("Logistic","KNN", "Trees"))
snsp_test_tree$model <- factor("Trees",c("Logistic","KNN", "Trees"))

roc_dat <- rbind(snsp_test_logistic, snsp_test_knn, snsp_test_tree)
head(roc_dat)
tail(roc_dat)

ggplot(data=roc_dat, aes(x=1-sp, y=sn, color=model)) +  geom_point()

## _____++++++++++++++++++++++++++++++++
# I trained the polynomial degree 2 for Auto data set
set.seed(1245)
test.index_Auto <- sample.int(nrow(Auto), round(nrow(Auto)*.15),replace=FALSE)
x_tr <- Auto[-test.index_Auto, ]

length(Auto$mpg) #390
dim(x_tr) #332

Model_best <- glm(mpg ~ poly(x, 2), data=x_tr, family = gaussian)

Model_best

mse <- function (i,y) {mean ( (y-i)^2 ) }

y_val <- Auto$mpg[test.index_Auto]
length(y_val) #58

x_val<- Auto[test.index_Auto, ]
length(x_val) #58

