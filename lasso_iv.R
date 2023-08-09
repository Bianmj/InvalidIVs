library(glmnet)
#Obtain the first-stage residuals and fitted value of X
#X is the exposure and G is the matrix for genetic variants
firststage <- function(X,G){
  n <- nrow(G)
  GGinv <- chol2inv(chol(t(cbind(1,G))%*%cbind(1,G)))
  alpha_hat <- GGinv%*%t(cbind(1,G))%*%X
  x.fitted <- cbind(1,G)%*%alpha_hat
  res <- X-x.fitted
  return(list(resid=res,xfit=x.fitted))
}

#Run the two-stage least squared regression and get the TSLS estimate,
#robust standard error and regular standard error
fit2sls<- function(Y,X,G){
  n <- nrow(G)
  qrg <- qr(cbind(1,G))
  x.fitted <- qr.fitted(qrg,X)
  res <- X-x.fitted
  qrxfit <- qr(cbind(1,x.fitted))
  beta_hat <- qr.coef(qrxfit,Y)
  mod2.res <- Y-cbind(1,X)%*%beta_hat
  x.fitted.2 <- t(cbind(1,x.fitted))%*%cbind(1,x.fitted)
  xxinv <-chol2inv(chol(x.fitted.2))
  inner = matrix(0,2,2)
  inner[1,1]=sum(mod2.res^2)
  inner[1,2]=inner[2,1]=sum(mod2.res^2*x.fitted)
  inner[2,2]=sum(mod2.res^2*x.fitted^2)
  var.robust <- xxinv%*%inner%*%xxinv
  var.norobust <- xxinv*sum(mod2.res^2)/(n-2)
  se.robust <- sqrt(var.robust[2,2])
  se.norobust <- sqrt(var.norobust[2,2])
  return(c(beta_hat=beta_hat[2],se.robust=se.robust, se.norobust= se.norobust))
  
  }


#Get the G-X and G-Y association estimates and their standard errors by linear regression
#X and Y are the exposure and outcome variables, and G is the matrix for genetic variants
#outputs: bx is a vector for the estimate of G-X association; by is a vector for the estimate of G-Y association;
#se_x is a vector for the standard error of bx; se_y is a vector the standard error of by 
est <- function(X,Y,G){
  M = ncol(G); n = nrow(G)
  bx = by=se_x=se_y=rep(0,M)
  for (i in 1:M){
    GGinv = chol2inv(chol(t(cbind(1,G[,i]))%*%cbind(1,G[,i])))
    beta_x = GGinv %*%t(cbind(1,G[,i]))%*%X
    resid_x = X-cbind(1,G[,i])%*%beta_x
    sigma2_x = sum(resid_x^2)/(n-M-1)
    bx[i]=beta_x[2]
    se_x[i]=sqrt(sigma2_x*GGinv[2,2])
    beta_y = GGinv %*%t(cbind(1,G[,i]))%*%Y
    resid_y = Y-cbind(1,G[,i])%*%beta_y
    sigma2_y = sum(resid_y^2)/(n-M-1)
    by[i]=beta_y[2]
    se_y[i]=sqrt(sigma2_y*GGinv[2,2])
  }
  cbind(bx,by,se_x,se_y)
}




#Apply the Lasso method by conditional on the first-stage residuals, exposure (X) and genetic variants (G) to find the invalid IVs
cv.my_alas <- function(Y,X,G){
  est_out <- est(X=X,Y=Y,G=G)
  se_y <- est_out[,4]; bx <- est_out[,1]; by <- est_out[,2]
  S = diag(se_y^-2)
  qrg <- qr(cbind(1,G))
  x.fitted <- qr.fitted(qrg,X)
  res <- X-x.fitted
  M <- ncol(G)
  n <- nrow(G)
  fs<- firststage(X,G)
  res <- fs$resid
  xfitted <- fs$xfit
  #Apply the lasso regression, since X and first-stage residuals(res) are not penalized,
  #so the penalty. factor is set to 0 for these two variables, and we only penalize the coefficient
  # of each genetic variant, which represents the potential directional pleiotropy effect on the outcome Y
  # The basic idea is to identify those invalid IVs with coefficients not equal to zero. 
  las_fit = glmnet(cbind(X,res,G),Y,penalty.factor=c(rep(0,2),rep(1,M)),intercept = FALSE)
  lamb = las_fit$lambda
  lamseq = sort(lamb)
  lamlen = length(lamseq)
  rse =c()
  # Applying the heterogeneity-stopping rule 
  # In glmnet, lambda values are arranged in a descending sequence. 
  # As a result, we analyze the shrinkage result from the final lambda, which holds the highest numerical value. 
  # The model corresponding to the highest lambda value will encompass the maximum count of valid IVs, 
  # while as lambda values decrease, the count of invalid IVs tends to increase.
  for (i in 1:lamlen){
    av = which(las_fit$beta[-c(1,2), (lamlen - i + 1)] == 0)
    # we will only consider the models with greater than 2 IVs.
    if (any(is.na(av)) | length(av)<=2) {next}
    mod = lm(S[av, av]^(1/2) %*% by[av] ~ S[av, av]^(1/2) %*% bx[av] - 1)
    ivw_est= mod$coefficients
    w= bx[av]^2/se_y[av]^2
    chival <- qchisq(0.95,df=length(av)-1)
    # rse is a matrix with four rows. The first row contains the Cochran'Q stats, and the second row is for the 
    # critical value from a Chi-squared distribution. The third row includes the number of valid IVs for
    # each lambda, and the fourth row is the associated value of lambda.
    rse=cbind(rse,c(sum(w*(by[av]/bx[av]-ivw_est)^2),chival, length(av),lamseq[i]))
  }
  nlam= ncol(rse)
  het = which(rse[1,] < rse[2,])
  if(length(het)==0){
    lam_pos = nlam
  } else{
    # We want to find the model with the greatest number of valid IVs.
    max.av<- max(rse[3,het])
    # Check whether multiple lambda values lead to the same highest count of valid IV.
    if(length(which(rse[3,het]==max.av))>1){
    # If there are multiple lambda values, then we select the lambda which has the smaller Cochran's Q stats
      lam_pos = which(lamseq==min(rse[4,het[which(rse[1,het]==min(rse[1,het[which(rse[3,het]==max.av)]]))]]))
    }
    # Else we just select the unique lambda with the largest number of valid IVs.
    else{lam_pos = which(lamseq==min(rse[4,het[which(rse[3,het]==max.av)]]))}
  }
  
  # Get the estimate of phi for each genetic variant. If phi is non-zero, then it is called "invalid".
  phi = las_fit$beta[-c(1,2), (lamlen - lam_pos + 1)]
  invalid_ind = which(phi != 0)
  
  # Return the index of invalid IVs
  invalid_ind
}


#Perform bootstrap procedure to get the bootstrap estimate of causal effect
#B is the number of bootstraps, and train is the number of observations for within-sample bootstrap
Boot.las <- function(X,Y,G,train,B){
  n = nrow(G)
  boot.est = rep(0,B)
  for (i in 1:B){
    shuffle <- sample(1:n, n, replace=TRUE)
    ind <- sample(1:n,train, replace=FALSE)
    insample=shuffle[ind];outsample=shuffle[-ind]
    invalid_ivs <- cv.my_alas(Y=Y[insample],X=X[insample],G=G[insample,])
    fit <- fit2sls(Y=Y[outsample],X=X[outsample],G=G[outsample,-invalid_ivs])
    boot.est[i]=fit[1]
  }
  boot.est
}
