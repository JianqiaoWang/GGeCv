library(R6)
library(glmnet)
library(scalreg)

Fitted.model = function(Dataset, intercept = T, scale.lasso = T, 
                        f = linear, g= linear,
                        sd.epsilon = sd.epsilon,
                        C = 0.75,
                        sd.v = sd.v,
                        cvnfolds = 10){
  
  Threshold.Choice.y <- "lambda.min"
  
  Threshold.Choice.z <- "lambda.min"
  
  if( identical(f, expit) ){
    
    if(is.null(sd.epsilon)){
    
    fit1 <-  glmnet::cv.glmnet(x = as.matrix(Dataset[!is.na(Dataset$Y),-c(1,2)]),
                     y = Dataset$Y[!is.na(Dataset$Y)], 
                     alpha = 1,
                     intercept = intercept,
                     family="binomial", nfolds = cvnfolds) # estimate the model with that
    
 #   fit2 <- cv.glmnet(x = W, y = Z, alpha = 1, family="binomial")
    
    Threshold.Choice.y <- "lambda.min"
    
   # Threshold.Choice.z <- "lambda.min"
    
    beta.hat = coef(fit1, s= Threshold.Choice.y)
    
    }else{
      
         lambda1 = (1:4)/10*sd.epsilon * sqrt(log(ncol(Dataset[,-c(1,2)]))/sum(!is.na(Dataset$Y)))
        
         lambda2 = C * sd.epsilon* sqrt(log(ncol(Dataset[,-c(1,2)]))/sum(!is.na(Dataset$Y)))
       #  
         lambda_vec =  c(lambda1, lambda2)
      #  
        fit1 <- glmnet::glmnet(x = as.matrix(Dataset[!is.na(Dataset$Y),-c(1,2)]),
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
      
    
    fit2 <- cv.glmnet(x = as.matrix(Dataset[!is.na(Dataset$Z),-c(1,2)]),
                      y = Dataset$Z[!is.na(Dataset$Z)], 
                      alpha = 1,
                      intercept = intercept,
                      family = "binomial",
                      nfolds = cvnfolds)    
    #Threshold.Choice.y <- "lambda.min"
    
     Threshold.Choice.z <- "lambda.min"
     
     gamma.hat = coef(fit2, s=Threshold.Choice.z)
     
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
    
    fit1 <-cv.glmnet(x = as.matrix(Dataset[!is.na(Dataset$Y),-c(1,2)]),
                     y = Dataset$Y[!is.na(Dataset$Y)], 
                     alpha = 1,
                     intercept = intercept,
                     nfolds = cvnfolds) # estimate the model with that
    
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
    
  #  H = model.matrix( ~ -1 + ., Dataset[!is.na(Dataset$Z),-c(1,2)])
    
   # H = model.matrix( ~  ., Dataset[!is.na(Dataset$Z),-c(1,2)])
    
    
    fit2 <- cv.glmnet(x = as.matrix(Dataset[!is.na(Dataset$Z),-c(1,2)]),
                        y = Dataset$Z[!is.na(Dataset$Z)], 
                        intercept = intercept,
                        family="gaussian",
                        nfolds = cvnfolds) 
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
        
        beta.0 = mean(Dataset$Y[!is.na(Dataset$Y)], na.rm = T)
        
        fit1 <- scalreg(X = as.matrix(Dataset[!is.na(Dataset$Y),-c(1,2)]),
                        y = Dataset$Y[!is.na(Dataset$Y)] - beta.0, LSE = T)
        
       # fit2 <- scalreg(X = as.matrix( cbind(rep(1, sum(!is.na(Dataset$Z))),
       #                                       Dataset[!is.na(Dataset$Z),-c(1,2)])), 
       #                 y = Dataset$Z[!is.na(Dataset$Z)], LSE = T)
        
        beta.hat = c(beta.0, fit1$coefficients)
        
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
    
    gamma.0 = mean(Dataset$Z[!is.na(Dataset$Z)], na.rm = T)
    
    
    fit2 <- scalreg(X = as.matrix( Dataset[!is.na(Dataset$Z),-c(1,2)]), 
                    y = Dataset$Z[!is.na(Dataset$Z)] - gamma.0, LSE = T)
    
    # beta.hat = fit1$coefficients
    
    gamma.hat =  c(gamma.0, fit2$coefficients)
    
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