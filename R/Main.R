#' Function for estimating the Genetic Covariance
#'
#' @param Dataset the `Dataset` argument should be a n * (p + 2) input matrix and each row is an observation vector. First two columns are two interested outcomes where missing values are allowed.
#'  The last p columns are collected genotypes matrix for n individuals
#' @param type 1-character string giving the type of "split" for sample splitting and cross-fitting procedure is implemented ; "whole" for using whole samples for both fitting the data and construct the confidence interval.
#' @param intercept Should intercept be fitted (default=TRUE) or set to zero (FALSE)?
#' @param scale.lasso Should the scale lasso or glmnet be used? 
#' @param f Specified function form f(x_i beta) for outcome 1. 'linear' for linear model and 'logit' for logistic model of binary outcomes
#' @param g Specified function form g(x_i gamma) for outcome 2 . 'linear' for linear model and 'logit' for logistic model of binary outcomes
#' @param sd.epsilon,sd.v   prior knownlege on the standard deviation of error terms epsilon and v
#'
#' @return A vector of two elements: estimated genetic covariance and standard errors
#' @export
#'
#' @examples
#' 
#' 

GECV.EST.MAIN = function(Dataset, type = "split", 
                         intercept = T,scale.lasso = T, 
                         f = linear , g = linear,
                         sd.epsilon = NULL, sd.v = NULL
                         ){
  
  if(type == "whole"){
    
    Fitted.parameter = Fitted.model(Dataset, 
                                    intercept = intercept, 
                                    scale.lasso = scale.lasso,
                                    f = f,
                                    g= g,
                                    sd.epsilon = sd.epsilon,
                                    sd.v = sd.v)
    
    result = GECV.EST(Dataset, 
                      beta.hat = Fitted.parameter$beta.hat,
                      gamma.hat = Fitted.parameter$gamma.hat,
                      f = f,
                      g = g)
    
  #  print(result)
    
  }
  
  if(type == "split"){
    
    ##################################
    # Decide what's I.1 and I.2
    ##################################
    
    group = paste0(as.numeric(!is.na(Dataset$Y)), as.numeric(!is.na(Dataset$Z)))
    
    I.1 = vector()
    
    for(i in unique(group)){
      
      temp = which(group == i)
      
      I.1 = c(I.1, sample(temp, floor(length(temp)/2)) )
      
    }
    
    I.2 = setdiff(1:nrow(Dataset), I.1)
    
    ###################################
    
    Fitted.parameter = Fitted.model(Dataset[I.1,], 
                                    intercept = intercept, scale.lasso = scale.lasso,
                                    f = f, g = g,
                                    sd.epsilon = sd.epsilon,
                                    sd.v = sd.v)
    
    result.1 = GECV.EST(Dataset[I.2,], 
                        beta.hat = Fitted.parameter$beta.hat,
                        gamma.hat = Fitted.parameter$gamma.hat,
                        f = f,
                        g = g)
    
    Fitted.parameter = Fitted.model(Dataset[I.2,], intercept = intercept, 
                                    scale.lasso = scale.lasso,
                                    f = f, g = g,
                                    sd.epsilon = sd.epsilon,
                                    sd.v = sd.v)
    
    result.2 = GECV.EST(Dataset[I.1,], 
                        beta.hat = Fitted.parameter$beta.hat,
                        gamma.hat = Fitted.parameter$gamma.hat,
                        f = f,
                        g = g)
    
    result = c(Est = 0.5*(result.1$Est + result.2$Est),
               sd = 0.5*sqrt(result.1$sd^2 + result.2$sd^2),
               class.est.sd = 0.5*sqrt(result.1$class.est.sd^2 + result.2$class.est.sd^2) )
  }
  
  return(result)
  
}

#' Title
#'
#' @param Dataset 
#' @param beta.hat 
#' @param gamma.hat 
#' @param f 
#' @param g 
#' @param case.control 
#' @param p1 
#' @param p0 
#'
#' @return
#' @export
#'
#' @examples
#' 
GECV.EST = function(Dataset, beta.hat, gamma.hat, f , g, case.control = F, p1 = NULL, p0 = NULL){
  
  n.z = sum(!is.na(Dataset$Z))
  
  n.y = sum(!is.na(Dataset$Y))
  
  method = GeCv$new(X = as.matrix(Dataset[,-c(1,2)]), 
                    beta.hat = beta.hat, 
                    gamma.hat = gamma.hat, 
                    Y = Dataset$Y, 
                    Z = Dataset$Z, 
                    f = f, g = g)
  
  if(case.control){
    
    method$case.control(p1 = p1, p0 = p0, 
                        pi1 = sum(Dataset$Y == 1, na.rm = T)/ sum(!is.na(Dataset$Y)),
                        pi0 = sum(Dataset$Y == 0, na.rm = T)/ sum(!is.na(Dataset$Y ))  
    )
    
  }
  
  Estimate = method$estimator()
  
  Estimate = method$cov.est.hat()
  
  class.est.var = Estimate[2]
  
  Estimator = Estimate[1]
  
  Estimator.Var = max(class.est.var, 
                      max(var(Dataset$Y, na.rm = T)/n.z^(3/2) , var(Dataset$Z, na.rm = T)/n.y^(3/2)))
  
  res = data.frame(Est = Estimator, 
                   sd = sqrt(Estimator.Var), 
                   class.est.sd = sqrt(class.est.var) )
  
  return(res)
}

