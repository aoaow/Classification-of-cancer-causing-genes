library(caret) 
library(corrplot) 
library(plyr)

## load data
train <- read.csv("training.csv") 
test <- read.csv("test.csv")

## change column names
colnames(train) <- c("id", paste("x",1:97,sep=""),"class") 
colnames(test) <- c("id", paste("x",1:97,sep=""),"class")

## functions
mm_scale <- function(x){
  # normalizes values of a vector x such that all values range from 0 to 1 # params: x: the vector to be normalized
  (x-min(x))/(max(x)-min(x))
}

add_log_cols <- function(df,cols){
  # function that will create a log column from a column that has been normalized
  # params: df: the dataframe to add the columns to
  #         cols:
  # return: returns a vector containing the prediction (0, 1, or 2) for each gene
  log_df <- log(df[cols]+1) # scale by 1 for 0s then log
  colnames(log_df) <- paste("log", colnames(log_df), sep = "_") # add log_ prefix 
  df <- cbind(df,log_df)
  return(df) 
  }
WCA <- function(data,lev=NULL,model=NULL){
  # custom WCA metric to use in cross validation
  # in place of the default Accuracy metric
  # params: data: a dataframe containing the observed
  # and predicted values for class
  # (other params are boilerplate)
  # return: the WCA score for the model
  # check that data is correct format
  if(length(data$obs) != length(data$pred)){ 
    stop("length(data$obs)␣!=␣length(data$pred)")
  }
  # Scores for correct predictions
  score_vec <- c(1, 20, 20) # 1 NG, 20 OG, 20 TSG 
  gene_vec <- c("NG","OG", "TSG")
  # Running total of WCA score
  WCA <- 0
  for(i in 0:2){
    # check that at least one appearance of class i 
    obs_i <- data$pred==gene_vec[i+1]
    if(sum(obs_i) == 0){
      stop("sum(obs_i)␣==␣0␣(data␣mislabeled?)") 
      }
    # Calculate points scored for each class
    score_i <- score_vec[i + 1]*sum(data$obs[obs_i] == data$pred[obs_i]) 
    # Add to running total
    WCA <- WCA + score_i
  }
  names(WCA) = "WCA"
  return(WCA) 
  }
pred_LR <- function(ptable){
  # function that will classify genes based on custom thresholds
  # params: ptable: a dataframe object containing probabilities of each gene
  # type (NG, OG, TSG) as columns
  # return: returns a vector containing the prediction (0, 1, or 
  i = 0
  preds = list() # vector to hold predictions
  for(row in 1:nrow(ptable)){
    if(ptable[row,1] >= 0.9){ # if prob NG >= 0.90, classify as NG
      preds[row] = 0 }
    else{
      if(ptable[row,2] >= ptable[row,3]){ # classify OG if OG prob
        preds[row] = 1 }
      else{ # otherwise classify TSG 
        preds[row] = 2
        } 
      }
  }
  preds <- unlist(preds) return(preds)
}

## variable selection
# first narrow down to important variables, by choosing those with p values < 0.10
glm.sig_p <- glm(class~.-id, data = train[1:99]) 
summary(glm.sig_p)

# of the remaining, filter out those with > 0.7 correlation
vars <- c("x1","x2","x3","x4","x5","x6","x7","x8","x14","x17","x18","x19","x20",
          "x21","x22","x23","x24","x25","x29","x30","x35","x38","x42","x43",
          "x45","x46","x50","x51","x52","x55","x62","x69","x72","x74","x75",
          "x76","x78","x83","x87","x90","x93","x94","x95","x96","x97") 
corr <- data.frame(cor(train[vars]))
corr <- corr []
# only keep variables that do not have a correlation > 0.7 with any other vars 
corr[!rowSums(corr > 0.70 & corr != 1),]

## transformations
# normalizing predictor variables 
train[2:98] <- apply(train[2:98],2,mm_scale) 
test[2:98] <- apply(test[2:98],2,mm_scale)

# adding log-transformed cols to train and test
train <- add_log_cols(train ,paste("x",c(3,4,5,7,17,19,22,24,25,29,38,43,75,78,83,90,97), sep=""))
test <- add_log_cols(test,paste("x",c(3,4,5,7,17,19,22,24,25,29,38,43,75,78,83,90,97), sep=""))

# adding factor version of class
train$class_f <- as.factor(train$class)
train$class_f <- revalue(train$class_f, c("0"="NG", "1"="OG","2"="TSG")) 
test$class_f <- test$class

## model and predicted probabilites
train_control <- trainControl(method="cv", number = 10,
                              log_fit <- form =
                                summaryFunction = WCA , selectionFunction="best", classProbs = TRUE , savePredictions = TRUE)
log_fit <- train(
  class_f ~x1 +log_x3+log_x4+log_x7+log_x19+x20+x21+log_x22+log_x25+log_x38
  +x42+log_x43+x45+x51+x52+log_x75+log_x78+log_x83+log_x97, train ,
  data =
    trControl = train_control , method = "multinom", metric = "WCA",
  maximize = TRUE
)
p.table <- predict(log_fit,newdata = test,"prob")

# using classifier function to classify probabilities for test set
test_preds <- pred_LR(p.table)
# pair predictions with test set ids and write to csv
df <- data.frame("id" = test$id, "class" = test_preds)

write.csv(df,"submission.csv", row.names = FALSE)