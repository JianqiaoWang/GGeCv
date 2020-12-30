#' Title
#'
#' @param Dataset 
#' @param intercept 
#' @param scale.lasso 
#' @param f 
#' @param g 
#' @param sd.epsilon 
#' @param C 
#' @param sd.v 
#' @param cvnfolds 
#'
#' @return
#' @import glmnet
#' @export
#'
#' @examples
#' 
Fitted.model = function(Dataset, intercept = T, scale.lasso = T, 
                        f = linear, g= linear,
                        sd.epsilon = sd.epsilon,
                        C = 0.75,
                        sd.v = sd.v,
                        cvnfolds = 3){
  
  Threshold.Choice.y <- "lambda.min"
  
  Threshold.Choice.z <- "lambda.min"
  
  if( identical(f, expit) ){
    
    if(is.null(sd.epsilon)){
      
    lambda_vec = seq(0, 150, 1.5)/100 * sd(Dataset$Y[!is.na(Dataset$Y)]) * 
      sqrt(
        log(ncol(Dataset[,-c(1,2)]))/ (sum(!is.na(Dataset$Y)) * (cvnfolds -1)/cvnfolds)
          ) 
    
    fit1 <-cv.glmnet(x = as.matrix(Dataset[!is.na(Dataset$Y),-c(1,2)]),
                     y = Dataset$Y[!is.na(Dataset$Y)], 
                     alpha = 1,
                     lambda = lambda_vec,
                     intercept = intercept,
                     family="binomial", nfolds = cvnfolds) # estimate the model with that
    
    lambda.select = fit1$lambda.min * sqrt((cvnfolds -1)/cvnfolds  )
 #   fit2 <- cv.glmnet(x = W, y = Z, alpha = 1, family="binomial")
    
    fit1 <-glmnet(x = as.matrix(Dataset[!is.na(Dataset$Y),-c(1,2)]),
                     y = Dataset$Y[!is.na(Dataset$Y)], 
                     alpha = 1,
                     lambda = ((1:5)/4) * lambda.select,
                     intercept = intercept,
                     family="binomial") # estimate the model with that
    
   # Threshold.Choice.y <- "lambda.min"
    
   # Threshold.Choice.z <- "lambda.min"
    
    beta.hat = coef(fit1, s= lambda.select)
    
    }else{
      
         lambda1 = (1:4)/10*sd.epsilon * sqrt(log(ncol(Dataset[,-c(1,2)]))/sum(!is.na(Dataset$Y)))
        
         lambda2 = C * sd.epsilon* sqrt(log(ncol(Dataset[,-c(1,2)]))/sum(!is.na(Dataset$Y)))
       #  
         lambda_vec =  c(lambda1, lambda2)
      #  
        fit1 <-glmnet(x = as.matrix(Dataset[!is.na(Dataset$Y),-c(1,2)]),
                      y = Dataset$Y[!is.na(Dataset$Y)],
                      alpha = 1,
                      intercept = intercept,
                      lambda = lambda_vec,
                      family="binomial") # estimate the model with that
        #
        beta.hat = coef(fit1,s= lambda2)
      
    }
    
    #  resid = Dataset$Y - predict(fit1,
    #           newx = as.matrix(Dataset[,-c(1,2)]), 
    #           type = "response",
    #           s = "lambda.min")
    #  
    #   lambda1 = (1:4)/10*sd(resid) * sqrt(log(800)/sum(!is.na(Dataset$Y)))
    #  
    #   lambda2 = 0.9*sd(resid)* sqrt(log(800)/sum(!is.na(Dataset$Y)))
    # #  
    #   lambda_vec =  c(lambda1, lambda2)
    #  
    #   fit1 <-glmnet(x = as.matrix(Dataset[,-c(1,2)]),
    #                 y = Dataset$Y, 
    #                 alpha = 1,
    #                 intercept = intercept,
    #                 lambda = lambda_vec,
    #                 family="binomial") # estimate the model with that
    #   #
    #   beta.hat = coef(fit1,s= lambda2)


   # gamma.hat = coef(fit2, s=Threshold.Choice.z)
    
  }
  
  if(identical(g, expit)){
    
    #fit1 <-cv.glmnet(x = X, y = Y, alpha = 1, family="binomial") #estimate the model with that
    if(is.null(sd.v)){
      
      lambda_vec = seq(0, 150, 1.5)/100 * sd(Dataset$Z[!is.na(Dataset$Z)]) * 
        sqrt(
          log(ncol(Dataset[,-c(1,2)]))/ (sum(!is.na(Dataset$Z)) * (cvnfolds -1)/cvnfolds)
        ) 
      
    fit2 <- cv.glmnet(x = as.matrix(Dataset[!is.na(Dataset$Z),-c(1,2)]),
                      y = Dataset$Z[!is.na(Dataset$Z)], 
                      alpha = 1,
                      lambda = lambda_vec,
                      intercept = intercept,
                      family = "binomial",
                      nfolds = cvnfolds)    
    #Threshold.Choice.y <- "lambda.min"
    
    # Threshold.Choice.z <- "lambda.min"
    
    lambda.select = fit2$lambda.min * sqrt( (cvnfolds -1)/cvnfolds  )
    #   fit2 <- cv.glmnet(x = W, y = Z, alpha = 1, family="binomial")
    
    fit2 <- glmnet(x = as.matrix(Dataset[!is.na(Dataset$Z),-c(1,2)]),
                      y = Dataset$Z[!is.na(Dataset$Z)], 
                      alpha = 1,
                      lambda = ((1:5)/4) * lambda.select,
                      intercept = intercept,
                      family = "binomial") 
    
    
     gamma.hat = coef(fit2, s = lambda.select)
     
    }else{
      
      lambda1 = (1:4)/10*sd.v * sqrt(log(ncol(Dataset[,-c(1,2)]))/sum(!is.na(Dataset$Z)))
      
      lambda2 = C * sd.v* sqrt(log(ncol(Dataset[,-c(1,2)]))/sum(!is.na(Dataset$Z)))
      #  
      lambda_vec =  c(lambda1, lambda2)
      #  
      fit2 <-glmnet(x = as.matrix(Dataset[!is.na(Dataset$Z),-c(1,2)]),
                    y = Dataset$Z[!is.na(Dataset$Z)],
                    alpha = 1,
                    intercept = intercept,
                    lambda = lambda_vec,
                    family="binomial") # estimate the model with that
      #
      gamma.hat = coef(fit2,s= lambda2)
      
    }
    
    #beta.hat = coef(fit1, s= Threshold.Choice.y)
     
     #   resid = Dataset$Z - predict(fit2, 
     #                               newx = as.matrix(Dataset[,-c(1,2)]),
     #                               type = "response",
     #                               s = "lambda.min")
     #  # 
     #  # xres = t(as.matrix(Dataset[,-c(1,2)])) %*% resid/sum(!is.na(Dataset$Z))
     #  
     #   lambda1 = (1:4)/10*sd(resid) * sqrt(log(800)/sum(!is.na(Dataset$Z)))
     #   
     #   lambda2 = 1*sd(resid) * sqrt(log(800)/sum(!is.na(Dataset$Z)))
     # #  
     #   lambda_vec =  c(lambda1, lambda2)
     #   
     #   fit2 <-glmnet(x = as.matrix(Dataset[,-c(1,2)]),
     #                 y = Dataset$Z, 
     #                 alpha = 1,
     #                 intercept = intercept,
     #                 lambda = lambda_vec,
     #                 family="binomial") # estimate the model with that
     # #  #
     #   gamma.hat = coef(fit2,s= lambda2)
     #  
     # 
     # 
     # fit1 <-glmnet(x = as.matrix(Dataset[,-c(1,2)]),
     #               y = Dataset$Z, 
     #               alpha = 1,
     #               intercept = intercept,
     #               lambda = c(lambda1, 2*lambda1, lambda2, 2*lambda2),
     #               family="binomial") # estimate the model with that
     # 
     # gamma.hat = coef(fit2, s=lambda1)
     

  }
  
  
    
  
  if(identical(f, linear)  & (!scale.lasso) ){
    
    fit1 <-cv.glmnet(x = as.matrix(Dataset[,-c(1,2)]),
                     y = Dataset$Y, alpha = 1,
                     intercept = intercept ) # estimate the model with that

  #  fit2 <- cv.glmnet(x = as.matrix(Dataset[,-c(1,2)]),
  #                    y = Dataset$Z, alpha = 1,
  #                    intercept = intercept)

     beta.hat = coef(fit1, s= Threshold.Choice.y)

   #  gamma.hat = coef(fit2, s= Threshold.Choice.z)
      
  }
  
  if(identical(g, linear) & (!scale.lasso) ){
    
   # fit1 <-cv.glmnet(x = as.matrix(Dataset[,-c(1,2)]),
    #                 y = Dataset$Y, alpha = 1,
     #                intercept = intercept ) # estimate the model with that
    
      fit2 <- cv.glmnet(x = as.matrix(Dataset[,-c(1,2)]),
                        y = Dataset$Z, alpha = 1,
                        intercept = intercept)
    
    #beta.hat = coef(fit1, s= Threshold.Choice.y)
    
    gamma.hat = coef(fit2, s= Threshold.Choice.z)
    
  }
      
      if(identical(f, linear) & scale.lasso){
        
        # if(intercept){
        # 
        # Dataset$Y = Dataset$Y - mean(Dataset$Y, na.rm = T)
        # 
        # Dataset$Z = Dataset$Z - mean(Dataset$Z, na.rm = T)
        # }
        
        fit1 <- scalreg(X = as.matrix( cbind( rep(1, sum(!is.na(Dataset$Y))),
                                              Dataset[!is.na(Dataset$Y),-c(1,2)])),
                        y = Dataset$Y[!is.na(Dataset$Y)], LSE = T)
        
       # fit2 <- scalreg(X = as.matrix( cbind(rep(1, sum(!is.na(Dataset$Z))),
       #                                       Dataset[!is.na(Dataset$Z),-c(1,2)])), 
       #                 y = Dataset$Z[!is.na(Dataset$Z)], LSE = T)
        
        beta.hat = fit1$coefficients
        
        # gamma.hat =  fit2$coefficients
        
        #beta.hat = fit1$lse$coefficients
        #gamma.hat = fit2$lse$coefficients
      }
  
  if(identical(g, linear) & scale.lasso){
    
    # if(intercept){
    # 
    # Dataset$Y = Dataset$Y - mean(Dataset$Y, na.rm = T)
    # 
    # Dataset$Z = Dataset$Z - mean(Dataset$Z, na.rm = T)
    # }
    
    #fit1 <- scalreg(X = as.matrix( cbind( rep(1, sum(!is.na(Dataset$Y))),
    #                                      Dataset[!is.na(Dataset$Y),-c(1,2)])),
    #                y = Dataset$Y[!is.na(Dataset$Y)], LSE = T)
    
    fit2 <- scalreg(X = as.matrix( cbind(rep(1, sum(!is.na(Dataset$Z))),
                                         Dataset[!is.na(Dataset$Z),-c(1,2)])), 
                    y = Dataset$Z[!is.na(Dataset$Z)], LSE = T)
    
    # beta.hat = fit1$coefficients
    
    gamma.hat =  fit2$coefficients
    
    #beta.hat = fit1$lse$coefficients
    #gamma.hat = fit2$lse$coefficients
  }
     
  # hat.sigma = fit1$hsigma
  # 
  # hat.Q.beta = t(beta.hat) %*% ( t(as.matrix(Dataset[,-c(1,2)]) ) %*% 
  #                                  as.matrix(Dataset[,-c(1,2)] ) ) %*% 
  #                                  beta.hat/nrow(Dataset)
  # 
  # 
  # hat.Q.2.beta = var(
  #   (as.matrix(Dataset[,-c(1,2)] ) %*% beta.hat)^2 - as.numeric(hat.Q.beta)
  # )
  # 
  # var.1 = (4*hat.sigma*hat.Q.beta + hat.Q.2.beta)/nrow(Dataset)
  # 
  #print(sqrt(var.1))
  

  return(data.frame(beta.hat = as.vector(beta.hat), 
                    gamma.hat =  as.vector(gamma.hat)))
  
}