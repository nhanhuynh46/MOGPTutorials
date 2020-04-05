covFull <- function(num_pop){
  
  # function to name parameters, given num_pop:
  listnamesParFull = function(num_pop){
    pop.list = c()
    pop.combination = combn(1:num_pop,2)
    for (l in 1:ncol(pop.combination)){
      current = paste("pop",pop.combination[,l][1],pop.combination[,l][2],sep="")
      pop.list = c(pop.list,current)
    }
    
    return(c("theta.ag","theta.yr",pop.list,"eta2"))
  }
  
  # create mapper:
  pop.combination = combn(1:num_pop,2)
  mapping = matrix(c(rep(NA,num_pop*num_pop)), 
                   nrow = num_pop, 
                   ncol = num_pop,
                   dimnames = list(1:num_pop,1:num_pop))
  diag(mapping) = 0
  for (i in 1:ncol(pop.combination)){
    mapping[pop.combination[1,i],pop.combination[2,i]] = i
    mapping[pop.combination[2,i],pop.combination[1,i]] = i
  }
  mapper <- function(X, Y) mapping[cbind(X, Y)]
  
  # function to compute the full-rank kernel, given num_pop and mapper:
  myCov <- covMan(
  kernel = function(x1, x2, par) {
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    d <- ncol(x1)
    npop <- length(unique(x1$x3))
    SS2 <- 0
    
    for (j in 1:2){
      Aj <- outer(x1[, j], x2[, j], "-")
      Aj2 <- Aj^2
      SS2 <- SS2 + Aj2 / (2*par[j]^2)
    }
    
    x1st = x1$x3; x2st = x2$x3;
    Aj = outer(x1st,x2st,FUN = mapper)
    for (j in 1:ncol(combn(1:npop,2))){
      Ajj = ifelse(Aj==j,1,0)
      SS2 <- SS2 + Ajj * par[j+2]
    }
    
    if (ncol(SS2)==nrow(SS2)){
      D6 <- exp(-SS2+diag(1e-5,nrow(SS2)))
    } else {
      D6 <- exp(-SS2)}
    kern <- par[2+ncol(combn(1:npop,2))+1] * D6
    
  },
  acceptMatrix = TRUE, # except matrices as arguments
  d = 3, # number of dimensions
  label = "myGauss",
  hasGrad = FALSE,
  parNames = listnamesParFull(num_pop), # names of the parameters
  parLower = do.call("<-",list(listnamesParFull(num_pop),c(5,5,rep(1e-6,ncol(pop.combination)),1e-6))), # lower bound
  parUpper = do.call("<-",list(listnamesParFull(num_pop),c(30,30,rep(0.25,ncol(pop.combination)),0.25))), # upper bound
  par = rep(NA,3+ncol(pop.combination)) # not given values of the parameters
)

}