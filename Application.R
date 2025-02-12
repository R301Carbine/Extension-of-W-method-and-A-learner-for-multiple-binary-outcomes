#Real data application
################################
###set data
#install.packages("speff2trial")
source("MainFunctions.R")
library("speff2trial")
data("ACTG175")
library(MASS)
#set the control and test group
test<-ACTG175
test$arms[test$arms==0]<-(-1)
sample<-test[test$arms>=(-1)&test$arms<=1,]
sam<-sample[!is.na(sample$cd496),]

#set the outcomes
s1<-sam$cd420
s2<-sam$cd496
s3<-sam$cd820
obje<-cbind(s1,s2,s3)
dichotomize_column <- function(column) {
  median_value <- median(column)
  ifelse(column >= median_value, 1, 0)
}

# Apply the function to each column of the matrix
Y <- apply(obje, 2, dichotomize_column)
#Y <- scale(obje, center = TRUE)

#set the covariances
x1<-sam[,2:10]
x2<-sam[,12:14]
x3<-sam[16]
x4<-sam[,19]
x5<-sam[,24]
X_dag<-cbind(x1,x2,x3,x4,x5)
X_dag <- X_dag[,-9]
X_ori <- X_dag
X_mean <- colMeans(X_ori)
X_sd <- apply(X_ori, 2, sd)
X <- scale(X_dag, center = TRUE)
colnames(X)[13:14] <- c("cd40","cd80")
colnames(X_ori)[13:14] <- c("cd40","cd80")

#set the diagonal matrix T
Tr<-diag(sam$arms)

#set the initial value
prop_score <- rep(1/2,nrow(X))

n <- nrow(X)
p <- ncol(X)
m <- ncol(Y)
#create A
A <- matrix(0,n,n)
for (i in 1:n){
  A[i,i] <- (Tr[i,i]*prop_score[i]+(1-Tr[i,i])/2)^(-1)
}

M <- matrix(0,n,n)
for (i in 1:n){
  M[i,i] <- (Tr[i,i]+1)/2-prop_score[i]
}


###############################
#rank selection
library(caret)
K <- 3
threshold<- (10)^(-10)
nMultiStart <- 2
max_iter <- 2000
rank <- c(1,2,3)

set.seed(200)
folds <- createFolds(Y[,1], k = K)
testLoss <- matrix(0,K,max(rank))

for (r in rank){
  for (i in 1:K){
    test_indices <- folds[[i]]
    #y
    train_Y <- Y[-test_indices,]
    test_Y <- Y[test_indices,]
    #X
    train_X <- X[-test_indices, ]
    test_X <- X[test_indices, ]
    #Tr
    train_Tr <- Tr[-test_indices, -test_indices]
    test_Tr <- Tr[test_indices, test_indices]
    #A
    train_A <- A[-test_indices, -test_indices]
    test_A <- A[test_indices, test_indices]
    res_RRRWB <- RRR_WB(max_iter, nMultiStart, threshold, nrow(train_Y), ncol(train_Y), ncol(train_X), r, train_Tr, train_X, train_Y, train_A)
    V <- res_RRRWB$V
    W <- res_RRRWB$W
    
    #calculate loss (original loss L_W(V,W,X))
    testLoss[i,r] <- sum(rowSums(test_Y*(log(1+exp(-test_Tr%*%test_X%*%W%*%t(V)))))/diag(test_A))
  }
}

#estimate using best rank
bestrank <- which.min(colMeans(testLoss))
best_RRRWB <- RRR_WB(max_iter, nMultiStart, threshold, nrow(Y), ncol(Y), ncol(X), bestrank, Tr, X, Y, A)
V_best <- best_RRRWB$V
W_best <- best_RRRWB$W
print( W_best%*%t(V_best) )
