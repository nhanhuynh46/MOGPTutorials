covICM_Notched <- function(num_pop, rank){
  
  # function to list out parameters given no. of populations and rank 
  listnamesPar = function(num_pop,rank){
    a_loading = c()
    for (i in 1:rank){
      a_loading = c(a_loading,paste("a",1:num_pop,i,sep=""))
    }
    
    return(c("theta.ag","theta.yr",a_loading))
  }

  # function to compute the ICM kernel, given num_pop and rank
  myCov <- covMan(
    kernel = function(x1, x2, par) {
      
      # unique attributes in x1 and x2:
      xx1 <- unique(x1[,c(1,2)])
      xx2 <- unique(x2[,c(1,2)])
      
      npop <- length(unique(x1$x3))
      SS2 <- 0
      
      for (j in 1:2){
        Aj <- outer(xx1[, j], xx2[, j], "-")
        Aj2 <- Aj^2
        SS2 <- SS2 + Aj2 / (2*par[j]^2)
      }
      # covariance over the design space:
      K <- exp(-SS2)
      
      # cross-covariance:
      cross <- par[(j+1):length(par)]
      groups <- split(cross,ceiling(seq_along(cross)/npop))
      B <- Reduce("+",lapply(groups,function(x) x%*%t(x)))
      
      # kronecker product between B and K:
      Kstar <- kronecker(B,K)
      
      # eigenMapMatMult --> C++ backend to fasten matrix multiplication
      # final covariance for notched setup:
      if (identical(x1,x2)){
        kern <- tcrossprod(ST %*% Kstar, ST)
        # kern <- eigenMapMatMult(eigenMapMatMult(ST, Kstar), t(ST))
        kern <- kern + diag(1e-5, nrow(kern))
      } else {
        kern <- ST %*% Kstar
        # kern <- eigenMapMatMult(ST, Kstar)
      }
      
    },
    acceptMatrix = TRUE, # except matrices as arguments
    d = 3, # number of dimensions
    label = "myGauss",
    hasGrad = FALSE,
    parNames = listnamesPar(num_pop, rank), # names of the parameters
    parLower = do.call("<-",list(listnamesPar(num_pop, rank),c(5,5,rep(1e-6,num_pop*rank)))), # lower bound
    parUpper = do.call("<-",list(listnamesPar(num_pop, rank),c(30,30,rep(0.5,num_pop*rank)))), # upper bound
    par = rep(NA,2+rank*num_pop) # not given values of the parameters
  )
  
}
