fitpredictModel<-
function(pa,lab){
        datax1 <- data.frame(features = pa,labels = lab)
        folds <- createFolds(y = datax1$labels,k = 5)
          for (j in 1:5) {
              traindata <- datax1[-folds[[j]],]
              testdata <- datax1[folds[[j]],]
              model <- glm(labels ~ .,family = binomial(link = "logit"),data = traindata)
              pred <- predict(model,testdata,type = "response", interval = "confidence")
              modelroc <- roc(testdata$labels,pred)
              modelauc <- modelroc$auc
              auc<-auc+modelauc
              }
              return(auc/5)
}              
