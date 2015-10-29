
library(protr)
library(readxl)
library(caret)
library(Biostrings)
library(RWeka)
library(dplyr)
library(e1071)
library(randomForest)
library(nnet)
library(Rcpi)


library(readxl)
data <- read_excel("data.xlsx")
Activity <- data$Activity
compound <- data[, 3:309]
protein <- data[, 311: 937]
### remove zero Variance
zeroVar <- function(dat) {
  out <- lapply(dat, function(x) length(unique(x)))
  want <- which(!out > 1)
  unlist(want)
}

compound <- compound[, -zeroVar(compound)]
protein <- protein[, -zeroVar(protein)]
#### centering and scaling
compound <- compound[, -nearZeroVar(compound, allowParallel = TRUE)]
protein <- protein[, -nearZeroVar(protein, allowParallel = TRUE)]

compound <- apply(compound, MARGIN = 2, FUN = function(x) (x - min(x))/diff(range(x)))
protein <- apply(protein, MARGIN = 2, FUN = function(x) (x -min(x))/diff(range(x)))
#### Correltation Removed 0.7 percent

cor_removed <- function(x) {
  x <- x[, !apply(x, 2, function(x) length(unique(x)) ==1)]
  raw <- cor(x)
  raw_2 <- raw[1: ncol(raw), 1: ncol(raw)]
  high <- findCorrelation(raw_2, cutoff = 0.7)
  remove <- x[, -high]
  input <- remove
  return(input)
}

compound <- cor_removed(compound) ### 33
protein <- cor_removed(protein) ### 26


#pca_compound <- prcomp(compound, retx = TRUE, center = TRUE, scale. = TRUE)
#pca_protein <- prcomp(protein, retx = TRUE, center = TRUE, scale. = TRUE)
#### 
CxP <- getCPI(compound, protein, type = "tensorprod")
CxP <- as.data.frame(CxP)
dfcompound <- names(data.frame(compound[,1:33]))
dfprotein <- names(data.frame(protein[,1:26]))
compoundNamecross <- rep(dfcompound, each = 26)
proteinNamecross <- rep(dfprotein, times = 26)
label <- paste(compoundNamecross, proteinNamecross, sep="_")
colnames(CxP) <- label


PxP <- getCPI(protein, protein, type = "tensorprod")
proteinName2 <- rep(dfprotein, times = 26)
proteinName1 <- rep(dfprotein, each = 26)
label_protein <- paste(proteinName1, proteinName2, sep = "_")
colnames(PxP) <- label_protein
index <- seq(1, 676, by = 27)
protein_selfcross <- PxP[, -index]
transposedIndexed_protein <- t(protein_selfcross)
index1 <- which(duplicated(transposedIndexed_protein))
removed_duplicated_protein <- transposedIndexed_protein[-index1, ]
PxP <- t(removed_duplicated_protein)

CxC <- getCPI(compound, compound, type = "tensorprod")
compoundName2 <- rep(dfcompound, times = 33)
compoundName1 <- rep(dfcompound, each = 33)
label <- paste(compoundName1, compoundName2, sep = "_")
colnames(CxC) <- label
index3 <- seq(1, 1089, by = 34)
compound_selfcross <- CxC[, -index3]
transposedIndexed_compound <- t(compound_selfcross)
index4 <- which(duplicated(transposedIndexed_compound))
removed_compound <- transposedIndexed_compound[-index4, ]
compound_finalcrossterms <- t(removed_compound)
CxC <- compound_finalcrossterms

#compound
C <- compound
#protein
P <- protein
#CxP
#CxC
#PxP
C_P <- cbind(C, P)
C_P_CxP_data_block_scale <- cbind(C, P, CxP) * (1/sqrt(length(C)+length(P)+length(CxP)))
#A_B_AxB_data <- cbind(affinity, A_B_AxB_data_block_scale)
C_P_CxC_data_block_scale <- cbind(C, P, 
                                  CxC) * (1/sqrt(length(C)+length(P)+length(CxC)))
#A_B_AxA_data <- cbind(affinity, A_B_AxA_data_block_scale)
C_P_PxP_data_block_scale <- cbind(C, P,
                                  PxP) * (1/sqrt(length(C)+length(P)+length(PxP)))
#A_B_BxB_data <- cbind(affinity, A_B_BxB_data_block_scale)
C_P_CxP_CxC_data_block_scale <- cbind(C, P, CxP,
                                      CxC) * (1/sqrt(length(C)+length(P)+length(CxP)+length(CxC)))
#A_B_AxB_AxA_data <- cbind(affinity, A_B_AxB_AxA_data_block_scale)
C_P_CxP_PxP_data_block_scale <- cbind(C, P, CxP,
                                      PxP) * (1/sqrt(length(C)+length(P)+length(CxP)+length(PxP)))
#A_B_AxB_BxB_data <- cbind(affinity, A_B_AxB_BxB_data_block_scale)
C_P_CxC_PxP_data_block_scale <- cbind(C, P, CxC,
                                      PxP) * (1/sqrt(length(C)+length(P)+length(CxC)+length(PxP)))
#A_B_AxA_BxB_data <- cbind(affinity, A_B_AxA_BxB_data_block_scale)
C_P_CxP_CxC_PxP_data_block_scale <- cbind(C, P, CxP, CxC, PxP) * (1/sqrt(length(C)+length(P)+
                                                                           length(CxC)+length(PxP)))
#A_B_AxB_AxA_BxB_data <- cbind(affinity, A_B_AxB_AxA_BxB_data_block_scale)

C <- cbind(Activity, C)
P <- cbind(Activity, P)
CxP <- cbind(Activity, CxP)
CxC <- cbind(Activity, CxC)
PxP <- cbind(Activity, PxP)
C_P <- cbind(Activity, C_P)
C_P_CxP <- cbind(Activity, C_P_CxP_data_block_scale)
C_P_CxC <- cbind(Activity, C_P_CxC_data_block_scale)
C_P_PxP <- cbind(Activity, C_P_PxP_data_block_scale)
C_P_CxP_CxC <- cbind(Activity, C_P_CxP_CxC_data_block_scale)
C_P_CxP_PxP <- cbind(Activity, C_P_CxP_PxP_data_block_scale)
C_P_CxC_PxP <- cbind(Activity, C_P_CxC_PxP_data_block_scale)
C_P_CxP_CxC_PxP <- cbind(Activity, C_P_CxP_CxC_PxP_data_block_scale)


#### J48 testing
J48_training <- function(x) {
  high <- subset(x, affinity == "High")
  low <- subset(x, affinity =="Low")
  results <- list(100)
  for (i in 1:100) {
    train_high <- sample_n(high, size = 54)
    test_high <- sample_n(high, size = 14)
    train_low <- sample_n(low, size = 51)
    test_low <- sample_n(low, size = 13)
    train <- rbind(train_high, train_low)
    test <- rbind(test_high, test_low)
    model_train <- J48(affinity~., data = train)
    summary <- summary(model_train)
    confusionmatrix <- summary$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}

### mean and SD value

mean_and_sd <- function(x) {
  c(round(mean(x, na.rm = TRUE), digits = 4),
    round(sd(x, na.rm = TRUE), digits = 4))
}

J48_train <- function(x) {
  ok <- J48_training(x)
  results <- data.frame(ok)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  return(data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
}

#protein_A_data
#protein_B_data
#cross_terms_data
#AxA_data
#BxB_data
#A_B_data 
#A_B_AxB_data 
#A_B_AxA_data
#A_B_BxB_data
#A_B_AxB_AxA_data
#A_B_AxB_BxB_data 
#A_B_AxA_BxB_data
#A_B_AxB_AxA_BxB_data

#### results for training J48
protein_A_training <- J48_train(protein_A_data)
protein_B_training <- J48_train(protein_B_data)
protein_cross_terms_training <- J48_train(cross_terms_data)
protein_AxA_training <- J48_train(AxA_data)
protein_BxB_training <- J48_train(BxB_data)
protein_A_B_training <- J48_train(A_B_data)
protein_A_B_AxB_training <- J48_train(A_B_AxB_data)
protein_A_B_AxA_training <- J48_train(A_B_AxA_data)
protein_A_B_BxB_training <- J48_train(A_B_BxB_data)
protein_A_B_AxB_AxA_training <- J48_train(A_B_AxB_AxA_data)
protein_A_B_AxB_BxB_training <- J48_train(A_B_AxB_BxB_data)
protein_A_B_AxA_BxB_training <- J48_train(A_B_AxA_BxB_data)
protein_A_B_AxB_AxA_BxB_training <- J48_train(A_B_AxB_AxA_BxB_data)

### function for 10 fold cross validation

J48_10fold <- function(x) {
  high <- subset(x, affinity == "High")
  low <- subset(x, affinity =="Low")
  results <- list(100)
  for (i in 1:100) {
    train_high <- sample_n(high, size = 54)
    test_high <- sample_n(high, size = 14)
    train_low <- sample_n(low, size = 51)
    test_low <- sample_n(low, size = 13)
    train <- rbind(train_high, train_low)
    test <- rbind(test_high, test_low)
    model_train <- J48(affinity~., data = train)
    eval_j48 <- evaluate_Weka_classifier(model_train, numFolds = 10, complexity = FALSE, seed = 1, class = TRUE)
    confusionmatrix <- eval_j48$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}

J48_cross_validation <- function(x) {
  ok <- J48_10fold(x)
  results <- data.frame(ok)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  return(data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
}

#### results for 10 fold cross validation

protein_A_cross_validation <- J48_cross_validation(protein_A_data)
protein_B_cross_validation <- J48_cross_validation(protein_B_data)
protein_crossterms_cross_validation <- J48_cross_validation(cross_terms_data)
protein_AxA_cross_validation <- J48_cross_validation(AxA_data)
protein_BxB_cross_validation <- J48_cross_validation(BxB_data)
protein_A_B_cross_validation <- J48_cross_validation(A_B_data)
protein_A_B_AxB_cross_validation <- J48_cross_validation(A_B_AxB_data)
protein_A_B_AxA_cross_validation <- J48_cross_validation(A_B_AxA_data)
protein_A_B_BxB_cross_validation <- J48_cross_validation(A_B_BxB_data)
protein_A_B_AxB_AxA_cross_validation <- J48_cross_validation(A_B_AxB_AxA_data)
protein_A_B_AxB_BxB_cross_validation <- J48_cross_validation(A_B_AxB_BxB_data)
protein_A_B_AxA_BxB_cross_validation <- J48_cross_validation(A_B_AxA_BxB_data)
protein_A_B_AxB_AxA_BxB_cross_validation <- J48_cross_validation(A_B_AxB_AxA_BxB_data)



#### J48 modeling testing results
J48_testing <- function(x) {
  high <- subset(x, affinity == "High")
  low <- subset(x, affinity =="Low")
  results <- list(100)
  for (i in 1:100) {
    train_high <- sample_n(high, size = 54)
    test_high <- sample_n(high, size = 14)
    train_low <- sample_n(low, size = 51)
    test_low <- sample_n(low, size = 13)
    train <- rbind(train_high, train_low)
    test <- rbind(test_high, test_low)
    model_train <- J48(affinity~., data = train)
    eval_external <- evaluate_Weka_classifier(model_train, newdata = test, numFolds = 0, complexity = FALSE, seed = 1, class = TRUE)
    confusionmatrix <- eval_external$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}

J48_external <- function(x) {
  ok <- J48_testing(x)
  results <- data.frame(ok)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  return(data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
}

### results for testing set


protein_A_testing <- J48_external(protein_A_data)
protein_B_testing <- J48_external(protein_B_data)
protein_crossterms_testing <- J48_external(cross_terms_data)
protein_AxA_testing <- J48_external(AxA_data)
protein_BxB_testing <- J48_external(BxB_data)
protein_A_B_testing <- J48_external(A_B_data)
protein_A_B_AxB_testing <- J48_external(A_B_AxB_data)
protein_A_B_AxA_testing <- J48_external(A_B_AxA_data)
protein_A_B_BxB_testing <- J48_external(A_B_BxB_data)
protein_A_B_AxB_AxA_testing <- J48_external(A_B_AxB_AxA_data)
protein_A_B_AxB_BxB_testing <- J48_external(A_B_AxB_BxB_data)
protein_A_B_AxA_BxB_testing <- J48_external(A_B_AxA_BxB_data)
protein_A_B_AxB_AxA_BxB_testing <- J48_external(A_B_AxB_AxA_BxB_data)
