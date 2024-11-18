##W&A-learner package
###functions
##RRR_WB: proposed method for W-method
##RRR_AB: proposed method for A-learner
#If you are running RRR_WB or RRR_AB, you can obtain the estimates of V, W.
#Once you have these estimates, you can calculate the HTE using the inner product of V and W.

#max_iter: maximum number of iterations
#nMultiStart: number of multistart runs (to avoid local solutions)
#threshold: threshold for termination of updates
#n: sample size, m: number of outcome variables, p: number of covariates, r: number of rank for heterogeneous treatment effects
#Tr: diagonal matrix of treatment assignment, X: matrix of explanatory variables
#Y: matrix of multiple outcomes, prop_score: vector of propensity score

#A: diagonal matrix of weighting (created as follows:)
#A <- matrix(0,n,n)
#for (i in 1:n){
#  A[i,i] <- (Tr[i,i]*prop_score[i]+(1-Tr[i,i])/2)^(-1)
#}

################################
##Logistic RRR for W-learner
RRR_WB <- function(max_iter, nMultiStart, threshold, n, m, p, r, Tr, X, Y, A){
  obj <- matrix(0, max_iter, nMultiStart)
  stop_iter <- rep(0,nMultiStart)
  X_dag <- A^(1/2)%*%Tr%*%X
  for (o in 1:nMultiStart){
    set.seed(o)
    W <- matrix(0.1*rnorm(p * r), nrow = p, ncol = r)
    set.seed(o)
    V <- qr.Q(qr(matrix(0.1*rnorm(m * r), nrow = m, ncol = r)))
    for (i in 1:max_iter){
      #calculate the Z
      logit <- pmin(pmax(-Tr %*% X %*% W %*% t(V), -100), 100)  # avoid the overflow
      #Z <- Tr%*%X%*%W%*%t(V)+4*(Y*(exp(-Tr%*%X%*%W%*%t(V))/(1+exp(-Tr%*%X%*%W%*%t(V)))))
      Z <- -logit+4*(Y*(exp(logit)/(1+exp(logit))))
      #calculate the objective function
      obj[i,o] <- norm(A^(1/2)%*%Z-A^(1/2)%*%Tr%*%X%*%W%*%t(V),type="F")^2
      
      #update W (by gradient method)
      W <- solve(t(X)%*%A%*%X)%*%t(X)%*%Tr%*%A%*%Z%*%V
      #calculate the objective function
      obj[i,o] <- norm(A^(1/2)%*%Z-A^(1/2)%*%Tr%*%X%*%W%*%t(V),type="F")^2
      
      #update V (by singular value decomposition)
      decomp <- svd(2*t(Z)%*%A%*%Tr%*%X%*%W)
      V <- decomp$u%*%t(decomp$v)
      #calculate the objective function
      obj[i,o] <- norm(A^(1/2)%*%Z-A^(1/2)%*%Tr%*%X%*%W%*%t(V),type="F")^2
      if (i > 1){
        if (abs(obj[i,o] - obj[(i-1),o]) < threshold){
          #message(message(paste("Converged at iteration", i)))
          stop_iter[o] <- i
          break
        }else if (i == max_iter){
          #message(message(paste("Converged at iteration", max_iter)))
          stop_iter[o] <- max_iter
        }# end elif
      }#end if
    } #end for i
    if (o==1){
      W_res <- W
      V_res <- V
    }else if (o>1 && nMultiStart != 1){
      if (obj[(stop_iter[(o-1)]),(o-1)] > obj[(stop_iter[o]),o]){
        W_res <- W
        V_res <- V
      } #end if
    } #end if 
  } #end for r
  invisible(list(obj = obj, V = V_res, W = W_res))
}
##

##Logistic RRR for A-learner
RRR_AB <- function(max_iter, nMultiStart, threshold, n, m, p, r, X, Y, M, pena){
  obj <- matrix(0, max_iter, nMultiStart)
  loss <- as.vector(NULL)
  stop_iter <- rep(0,nMultiStart)
  X_dag <- M%*%X
  for (o in 1:nMultiStart){
    set.seed(o)
    W <- matrix(0.1*rnorm(p * r), nrow = p, ncol = r)
    set.seed(o)
    V <- qr.Q(qr(matrix(0.1*rnorm(m * r), nrow = m, ncol = r)))
    for (i in 1:max_iter){
      logit <- pmin(pmax(-M %*% X %*% W %*% t(V), -700), 700)  # avoid the overflow
      #calculate the Z
      #Z <- M%*%X%*%W%*%t(V)+4*(Y*(exp(-M%*%X%*%W%*%t(V))/(1+exp(-M%*%X%*%W%*%t(V)))))
      Z <- -logit+4*(Y*(exp(logit)/(1+exp(logit))))
      #calculate the objective function
      obj[i,o] <- norm(Z-M%*%X%*%W%*%t(V),type="F")^2
      
      #update W (by gradient method)
      W <- solve(t(X)%*%M^2%*%X)%*%t(X)%*%M%*%Z%*%V
      
      #calculate the objective function
      obj[i,o] <- norm(Z-M%*%X%*%W%*%t(V),type="F")^2
      
      #update V (by singular value decomposition)
      decomp <- svd(2*t(Z)%*%M%*%X%*%W)
      V <- decomp$u%*%t(decomp$v)
      #calculate the objective function
      obj[i,o] <- norm(Z-M%*%X%*%W%*%t(V),type="F")^2
      if (i > 1){
        if (abs(obj[i,o] - obj[(i-1),o]) < threshold){
          #message(message(paste("Converged at iteration", i)))
          stop_iter[o] <- i
          break
        }else if (i == max_iter){
          #message(message(paste("Converged at iteration", max_iter)))
          stop_iter[o] <- max_iter
        }# end elif
      }#end if
    } #end for i
    if (o==1){
      W_res <- W
      V_res <- V
    }else if (o>1 && nMultiStart != 1){
      if (obj[(stop_iter[(o-1)]),(o-1)] > obj[(stop_iter[o]),o]){
        W_res <- W
        V_res <- V
      } #end if
    } #end if 
  } #end for o
  invisible(list(obj = obj, V = V_res, W = W_res))
}
##
################################
##create the biplot of V and W (if r=2)
#For V and W obtained from RRR_WB or RRR_AB, assign appropriate row names and column names to each.
#Then, input W into the "W_ori" argument and V into the "V_ori" argument of the "Biplot" function to create a biplot.

Biplot <- function(W_ori, V_ori){
  normalize <- function(mat) {
    t(apply(mat, 1, function(row) row / sqrt(sum(row^2))))  # 各行ベクトルをそのノルムで割る
  }
  W <- normalize(W_ori)
  V <- normalize(V_ori)
  df_W <- data.frame(x = W[, 1], y = W[, 2], label = rownames(W))
  df_V <- data.frame(x = V[, 1], y = V[, 2], label = rownames(V))
  
  #create unit circle
  theta <- seq(0, 2*pi, length.out = 100)
  circle <- data.frame(x = cos(theta), y = sin(theta))
  
  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +  # x軸の点線
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  # y軸の点線
    
    geom_segment(data = df_W, aes(x = 0, y = 0, xend = x, yend = y), 
                 arrow = arrow(length = unit(0.3, "cm")), color = "black") +
    geom_text(data = df_W, aes(x = x, y = y, label = label), vjust = -0.5, hjust = 0.5, color = "black") +
    
    geom_path(data = circle, aes(x = x, y = y), linetype = "dashed", color = "gray") +
    
    geom_segment(data = df_V, aes(x = 0, y = 0, xend = x, yend = y), 
                 arrow = arrow(length = unit(0.3, "cm")), color = "red") +
    geom_text(data = df_V, aes(x = x, y = y, label = label), vjust = -0.5, hjust = 0.5, color = "red") +
    
    xlim(c(-1,1)) + ylim(c(-1,1)) +
    theme_minimal() +  
    labs(title = "Biplot of W's and V's Row Vectors", x = "Dim 1", y = "Dim 2")
}

#########################