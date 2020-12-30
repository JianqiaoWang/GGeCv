#' Functions for generating data and coeficient vectors
#'
#' ... Complete documentation ...
#'
#' @description  A class for generating regressor matrix 
#' @field n numeric. 
#' @field p numeric. 
#' @field Cov matrix. 
#'
#' @return a matrix of data
#' @export
#'
#' @examples
#' 
Regressor <- setRefClass("Regressor",
                         fields = list(n = "numeric",
                                       p = "numeric",
                                       Cov = "matrix"),
                         methods = list(
                           
                           normal = function(){
                             
                             library(mvtnorm)
                             
                             X <- rmvnorm(n, mean = rep(0,p), sigma = Cov)
                             
                             return(X)
                           },
                           uniform = function(size = 1){
                             
                             X <- matrix(runif(n*p, min = -size, max = size), nrow =n, ncol = p)
                             
                             return(X)
                             
                           },
                           
                           SNP = function(){
                             
                             prob = 0.1
                             
                             normal.cop = normalCopula(0, dim= p ,dispstr="ar1")
                             
                             u.1 <- rCopula(n, normal.cop)
                             
                             X <- array(qbinom(u.1,size=2,prob=prob), dim=c(n, P))
                             
                             return(X)
                           },
                           geno = function(geno){
                             
                             if(n < nrow(geno)){
                             
                             X <- geno[sample(1:nrow(geno), n, replace = T),]#
                             }else{
                               
                               X <- geno[sample(1:nrow(geno), n, replace = T),]#
                               
                             }
                             
                             #X <- array(qbinom(u.1,size=2,prob=prob), dim=c(n, P))
                             
                             return(X)
                           }
                           
                           )
)

#' @description  A class for generating coefficient vectors 
#'
#' @field MAG . 
#' @field p numeric. 
#'
#' @return a p*2 matrix
#' @export
#'
#' @examples
#' 

Coef <- setRefClass("Coefficient",
                    fields = list(p = "numeric", MAG = "numeric"),
                    methods = list(
                      sparse = function(s){
                        Beta = c(rep(1,s),rep(0,p-s)) * MAG
                        return(Beta)
                      },
                      normal = function(beta.abs){
                        if(beta.abs == 1){
                          Beta = abs(c(round(rnorm(p,0,1),4)))
                        }else{
                          Beta = (c(round(rnorm(p,0,1),4)))
                        }
                        return(Beta)
                      },
                      decreasing = function(s){
                        
                        if(s == 0){
                          
                          Beta = rep(0,  p)
                          
                        }else{
                        
                        Beta = rep(MAG,  p) * c((1:s)/s,rep(0, p-s))
                        
                        }
                        
                        return(Beta)
                        
                      },
                      sparse.uniform = function(s){
                        Beta = c(runif(s,0,1),rep(0,p-s)) * MAG
                        return(Beta)
                      },
                      zero = function(){
                        
                        Beta = rep(0, p)
                        
                        return(Beta)
                      },
                      cross = function(s){
                        
                        Beta = MAG*(-1)^abs(outer(1:s, 1, "-"))
                        
                        Beta = c(Beta, rep(0, p -s))
                        
                        return(Beta)
                      } 
                      
                    )
)

Response <- setRefClass("Response",
                        fields = list(regressor = "matrix",
                                      coef = "vector",
                                      n = "numeric",
                                      pr = "vector"),
                        methods = list(
                          logit = function(a = 0){
                            
                            z = regressor %*% coef        # linear combination with a bias
                            
                            prob = 1/(1+exp(-z))  
                            
                            pr <<- 1/(1+exp(-z))         # pass through an inv-logit function
                            
                            y = rbinom(n,1,prob)      # bernoulli response variable
                            
                            return(y)
                            
                          },
                          probit = function(){
                            
                            z = rep(0, n) + regressor %*% coef        # linear combination with a bias
                            
                            pr<<-pnorm(z)   # pass through an inv-logit function
                            
                            y = rbinom(n,1,pr)      # bernoulli response variable
                            
                            return(y)
                            
                          },
                          linear = function(){
                            
                            pr <<- regressor %*% coef   # pass through an inv-logit function
                            
                            sd.ep = 1
                            
                            response <- regressor %*% coef + rnorm(n,0,1)
                            
                            return(response)
                          },
                          quad = function(){
                            
                            #pr <<- (regressor %*% coef)^2   # pass through an inv-logit function
                            
                            pr <<- -5 +sin(regressor[,1] * regressor[,2] * pi) +  
                              4*(regressor[,3] - 0.5)^2 + 2*regressor[,5] + regressor[,6] # pass through an inv-logit function
                            
                            sd.ep = 1
                            
                            response <- pr + rnorm(n,0,sd.ep)
                            
                            return(response)
                          }
                        )
)


Correlation <- function(pr1, pr2, r){
  
  oddsratio = (pr1/pr2)*((1 - pr2)/(1 - pr1))
  
  temp.r = r
  
  if(temp.r > 0){
  
  r = pmin(r, 0.4*pmin(sqrt(oddsratio),sqrt(1/oddsratio)))
  }
  
  if(temp.r < 0){
    
    prod = (pr1 * pr2) /( (1 - pr1)*(1-pr2) )
    
    r = pmax(r, -0.4*pmin(sqrt(prod),sqrt(1/prod)))
    
  }
  
  p_1_1 <- pr1 * pr2 + r* sqrt(pr1 * pr2*(1 - pr1)*(1 - pr2))
  
  p_1_0 <- pr1 - p_1_1
  
  p_0_1 <- pr2 - p_1_1
  
  p_0_0 <- 1 - p_1_1 - p_1_0 - p_0_1
  
  prob <- cbind('(0,0)'= p_0_0, '(1,0)'= p_1_0, '(0,1)'= p_0_1, '(1,1)'= p_1_1)
  
  dict <- vector(mode="list", length=4)
  
  names(dict) <- c('(0,0)','(1,0)', '(0,1)', '(1,1)')
  
  dict[[1]] <- c(0,0) 
  
  dict[[2]] <- c(1,0) 
  
  dict[[3]] <- c(0,1)
  
  dict[[4]] <- c(1,1)
  
  result <- apply(prob,1, function(x){
    dict[[ sample(c('(0,0)',  '(1,0)', '(0,1)', '(1,1)' ),1, replace = F, x)]]
  }
 )
  result = t(result)
  return(result)
}


Correlation_linear <- function(pr1, pr2, r){
  
  sd.1 = 1/2 * sd(pr1)
  
  sd.2 = 1/2 * sd(pr2)
  
  if(sd.1 < 0.01){
    sd.1 = 1
  }
  
  if(sd.2 < 0.01){
    sd.2 = 1
  }  
  
  rmv
  
  result = t(result)
  
  return(result)
}