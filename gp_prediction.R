gp.predict = function(newdata,gpmodel,list.noise,meanTr,typePred){
  y = matrix(gpmodel$y)
  X = gpmodel$X
  n = nrow(gpmodel$F)
  p = ncol(gpmodel$F)
  d = ncol(gpmodel$X)
  nms = gpmodel$inputNames   
  
  ## trend for training data:
  XTr = X[ , nms, drop = FALSE]
  tt = delete.response(terms(gpmodel))
  mt = model.frame(tt, data = data.frame(XTr))
  FTr = model.matrix(tt, data = mt)
  trend = FTr %*% gpmodel$betaHat
  
  ## trend for new data:
  if (meanTr=="linearAg"){
    FNew = model.matrix(~ newdata$x1 + as.factor(newdata$x3), newdata)
  } else {
    FNew = model.matrix(~ newdata$x1 + newdata$x2 + as.factor(newdata$x3), newdata)
    # linear in both Age and Year
  }
  trendNew = FNew %*% gpmodel$betaHat
  
  ## covariance:
  covXXnew = covMat(gpmodel$covariance,
                    X = gpmodel$X,
                    Xnew = newdata)
  K = covMat(gpmodel$covariance,
             X = gpmodel$X,
             checkNames = TRUE)
  pop.num = length(list.noise)
  freq.noise = table(X$x3)
  Sigma = diag(rep(list.noise,freq.noise))
  
  nn = nrow(newdata)
  SigmaNew = diag(rep(list.noise,each=nn/pop.num))
  covXnew = covMat(gpmodel$covariance,X=newdata)
  invCov = solve(K+Sigma)
  
  # posterior mean:
  condMean = trendNew + t(covXXnew)%*%invCov%*%(y-trend)
  # posterior covariance of the prediction:
  if (typePred=="fnoise"){
    condVar = covXnew + SigmaNew - t(covXXnew)%*%invCov%*%covXXnew
  } else {
    condVar = covXnew - t(covXXnew)%*%invCov%*%covXXnew
    # f noise-free
  }
  
  newdata$mean = condMean
  newdata$vr = diag(condVar)
  newdata$std = sqrt(newdata$vr)
  names(newdata) = c("age","year","popN","mean","vr","std")
  
  return(list(res = newdata, covMa = condVar))
}
