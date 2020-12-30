Print_Coverage = function(H){
  
  simuresult = as.data.frame(do.call(rbind, H))[-c(nsim + 1,nsim + 2),]
  
  true.para = c(H[[nsim + 1]])
  
  simuresult$lower <- simuresult[,"Est"] - simuresult[,"sd"] * z_alpha
  
  simuresult$upper <- simuresult[,"Est"] + simuresult[,"sd"] * z_alpha
  
  simuresult$lower.2 <- simuresult[,"Est"] - simuresult[,"class.est.sd"] * z_alpha
  
  simuresult$upper.2 <- simuresult[,"Est"] + simuresult[,"class.est.sd"] * z_alpha
  
  count <- 0
  
  count = sum(((true.para >= simuresult$lower) + (true.para <= simuresult$upper) == 2),na.rm = T)
  
  count.2 = sum(((true.para >= simuresult$lower.2) + (true.para <= simuresult$upper.2) == 2),na.rm = T)
  
  print(true.para)
  
  print(count/200)
  
  print(name)
  
  rmse.proposed = sqrt(sum((simuresult[,"Est"] - rep(true.para, length(simuresult[,"Est"] )))^2)/nsim) 
  
  rmse.plugin = rmse.proposed
  
  result.data = data.frame(coverage = count/nsim, 
                           coverprob.raw = count.2/nsim,
                           rmse.proposed = rmse.proposed, 
                           rmse.plugin = rmse.plugin, 
                           n = n, rho = rho, 
                           est.mean = mean(simuresult[,"Est"]),  
                           se = mean(simuresult[,"sd"]), 
                           type = type, s = sparsity )
  
  print(result.data)
  
  #count.2 = sum(((true.para[1] >= simuresult$lower.2) + (true.para[1] <= simuresult$upper.2) == 2),na.rm = T)
  
  return(c(count/nsim, count.2/nsim))
  
}


TRUEPARA = function(type){
  
  #print(Sigma[1:2,1:2])
  
  X <- Regressor(n = n, p = P, Cov = Sigma)$normal()
  
  #X <- Regressor(n = n, p = P, Cov = Sigma)$uniform(size = 1)
  
  #print(mean(X %*% beta.1))
  
  Response.1 <- Response(regressor = X, coef = beta.1, n = n)
  
  Response.2 <- Response(regressor = X, coef = gamma.1, n = n)
  
  if(type == "logit"){
    
    temp1 = Response.1$logit()
    
    temp2 =  Response.2$logit()
  }
  
  if(type == "linear"){
    temp1 = Response.1$linear()
    
    temp2 =  Response.2$linear()
  }
  if(type == "linear-logit"){
    temp1 = Response.1$linear()
    
    temp2 =  Response.2$logit()
  }
  
  if(type == "mis-linear-logit"){
    temp1 = Response.1$linear()
    
    temp2 =  Response.2$logit()
  }
  
  para1 = cov(Response.1$pr, Response.2$pr)
  
  #para2 = crossprod(Response.1$pr, Response.2$pr)/n
  
 # print(mean(X %*% beta.1))

  para2 = para1
  
  return(c(para1,para2))
}

#coef.choice = function(s){
  
# return(Coef(p = P,  MAG = coef)$decreasing(s = s))
  
#}
