
library(R6)

#GeCov class

GeCv <-  R6::R6Class("GeCov", list( X = "matrix",
                                beta = "vector",
                                gamma = "vector",
                                beta.hat = "matrix",
                                gamma.hat = "matrix",
                                Y = "vector",
                                Z = "vector",
                                f = "function",
                                g = "function",
                                inner = NULL,
                                Y.hat = NULL,
                                Z.hat = NULL,
                                eps.hat = NULL,
                                v.hat = NULL,
                                N = NULL,
                                N.Y = NULL,
                                N.Z = NULL,
                                Delta = NULL,
                                I.y = NULL,
                                I.z = NULL,
                                I = NULL,
                                mu.correct = NULL,
                                epsilon.mu.z.correct = NULL,
                                v.mu.y.correct = NULL,
                                weight = NULL,
                                initialize = function(Y, X,  beta.hat, f,  Z , W , gamma.hat, g
                                ){
                                  
                                  stopifnot(is.vector(as.vector(beta.hat)), 
                                            is.vector(as.vector(gamma.hat)),
                                            length(beta.hat) == length(gamma.hat))
                                  
                                  stopifnot(is.numeric(Y), is.numeric(Z))
                                  
                                  self$beta.hat = as.vector(beta.hat)
                                  
                                  self$gamma.hat = as.vector(gamma.hat)
                                  
                                  self$f = f
                                  
                                  self$g = g
                                  
                                  self$X = X
                                  
                                  self$Y = Y
                                  
                                  self$Z = Z
                                  
                                  self$weight = rep(1, length(self$Y))
                                  
                                },
                                case.control = function(p1,p0,pi1,pi0){
                                  
                                  w0 = p0/pi0 
                                  
                                  w1 = p1/pi1 
                                  
                                  self$weight = self$Y
                                  
                                  self$weight[is.na(self$Y)] = 1
                                  
                                  self$weight[which(self$Y == 1)] = w1
                                  
                                  self$weight[which(self$Y == 0)] = w0
                                  
                                  self$weight = self$weight #* 100
                                  
                                  self$beta.hat[1] = self$beta.hat[1] + log(w1/w0)
                                  
                                },
                                estimator =  function(){
                                  
                                  self$Y.hat <- self$f(self$X %*% self$beta.hat[-1] + self$beta.hat[1] )
                                  
                                  self$Z.hat <-  self$g(self$X %*% self$gamma.hat[-1] + self$gamma.hat[1] )
                                  
                                  self$eps.hat <- self$Y - self$Y.hat
                                  
                                  self$v.hat <- self$Z - self$Z.hat
                                  
                                  self$N = length(self$Y.hat)
                                  
                                  self$N.Y = sum(!is.na(self$Y))
                                  
                                  self$N.Z = sum(!is.na(self$Z))
                                
                                },
                                plug = function(){ # plug-in estimate
                                  
                                  mu.y.hat = mean(self$Y.hat) + sum(self$eps.hat)/self$N.Y
                                  
                                  mu.z.hat = mean(self$Z.hat) + sum(self$v.hat)/self$N.Z
                                  
                                  return(crossprod(self$Y.hat, self$Z.hat)/self$N - mu.y.hat*mu.z.hat  )
                                  
                                },
                                Inner.est = function(){ # estimate of covariance 
                                  
                                  self$eps.hat[is.na(self$eps.hat)] = 0
                                  
                                  self$v.hat[is.na(self$v.hat)] = 0
                                  
                                  Delta2 = self$weight * self$Y.hat * self$Z.hat /self$N + 
                                    self$weight * self$eps.hat * self$Z.hat /self$N.Y + 
                                    self$weight * self$v.hat * self$Y.hat/self$N.Z
                                  
                                  self$inner = sum(Delta2)
                                  
                                },
                                cov.est.hat = function(){
                                  
                                  ### cov estimation
                                  
                                  self$Inner.est()
                                  
                                  mu.y.hat = mean(self$weight * self$Y.hat) + sum( self$weight * self$eps.hat)/self$N.Y
                                  
                                  mu.z.hat = mean(self$weight * self$Z.hat) + sum( self$weight * self$v.hat)/self$N.Z
                                  
                                  estimator = self$inner - mu.y.hat * mu.z.hat
                                  
                                  ### cov variance 
                                  
                                  y.tilde = mean(self$Y, na.rm = T)
                                  
                                  z.tilde = mean(self$Z, na.rm = T)
                                  
                                  y.tilde = mu.y.hat
                                  
                                  z.tilde = mu.z.hat
                                  
                                  # self$Delta = self$weight * self$Y.hat * self$Z.hat /self$N + 
                                  #   self$weight * self$eps.hat * self$Z.hat /self$N.Y + 
                                  #   self$weight * self$v.hat * self$Y.hat/self$N.Z - 
                                  #   mu.y.hat * ( self$weight * self$Z.hat/self$N  + self$weight * self$v.hat/self$N.Z )-
                                  #   mu.z.hat * ( self$weight * self$Y.hat/ self$N + self$weight * self$eps.hat / self$N.Y)
                                  # 
                                  self$Delta = self$weight *(self$Y.hat - y.tilde)  * (self$Z.hat - z.tilde) /self$N + 
                                    self$weight * self$eps.hat * (self$Z.hat - z.tilde) /self$N.Y + 
                                    self$weight * self$v.hat * (self$Y.hat - y.tilde)/self$N.Z
                                  
                                  estimator.var = length(self$Delta)*var(self$Delta ) +  
                                    (sum(self$Delta) -  estimator)^2/self$N
                                  
                                 # plugest = self$plug() - mu.y*mu.z
                                  
                                  return(c(estimator, estimator.var))
                                },
                                cc.cov.est.hat = function(){
                                  
                                  ### cov estimation
                                  
                                  self$Inner.est()
                                  
                                  mu.y.hat = mean(self$weight * self$Y.hat) + sum( self$weight * self$eps.hat)/self$N.Y
                                  
                                  mu.z.hat = mean(self$weight * self$Z.hat) + sum( self$weight * self$v.hat)/self$N.Z
                                  
                                  y.tilde = mean(self$Y, na.rm = T)
                                  
                                  z.tilde = mean(self$Z, na.rm = T)
                                  
                                  estimator = self$inner - mu.y.hat  * mu.z.hat
                                  
                                  y.tilde = mu.y.hat
                                  
                                  z.tilde = mu.z.hat
                                  
                                  # self$Delta = self$weight * self$Y.hat * self$Z.hat /self$N + 
                                  #   self$weight * self$eps.hat * self$Z.hat /self$N.Y + 
                                  #   self$weight * self$v.hat * self$Y.hat/self$N.Z - 
                                  #   mu.y.hat * ( self$weight * self$Z.hat/self$N  + self$weight * self$v.hat/self$N.Z )-
                                  #   mu.z.hat * ( self$weight * self$Y.hat/ self$N + self$weight * self$eps.hat / self$N.Y)
                                  # 
                                  self$Delta = self$weight *(self$Y.hat - y.tilde)  * (self$Z.hat - z.tilde) /self$N + 
                                    self$weight * self$eps.hat * (self$Z.hat - z.tilde) /self$N.Y + 
                                    self$weight * self$v.hat * (self$Y.hat - y.tilde)/self$N.Z
                                  
                                
                                 # self$Delta = self$weight * self$Y.hat * self$Z.hat /self$N + 
                                 #   self$weight * self$eps.hat * self$Z.hat /self$N.Y + 
                                 #   self$weight * self$v.hat * self$Y.hat/self$N.Z
                                  
                                  
                                  
                                  estimator.var = 
                                    length(which(self$Y == 1))*var(self$Delta[which(self$Y == 1)])+
                                    length(which(self$Y == 0))*var(self$Delta[which(self$Y == 0)])+ 
                                    length(which(!is.na(self$Z) ))*var(self$Delta[!is.na(self$Z)])
                                    #(sum(self$Delta) -  estimator)^2/self$N
                                  
                                  # plugest = self$plug() - mu.y*mu.z
                                  
                                  return(c(estimator, estimator.var))
                                },
                                cov.est.tilde = function(){
                                  
                                  mu.y = mean(self$Y)
                                  
                                  mu.z = mean(self$Z)
                                  
                                  self$Delta = self$Y.hat.X.U * self$Z.hat.X.U/length(self$Y.hat.X.U) + 
                                    self$res.hat.epsilon * self$Z.hat.X.U/length(self$I.y) + 
                                    self$res.hat.v * self$Y.hat.X.U/length(self$I.z)
                                  
                                  Delta2 = (self$Y.hat.X.U - mu.y)  * (self$Z.hat.X.U - mu.z)/length(self$Y.hat.X.U) + 
                                    self$res.hat.epsilon * (self$Z.hat.X.U - mu.z) /length(self$I.y) + 
                                    self$res.hat.v * (self$Y.hat.X.U - mu.y)/length(self$I.z)
                                  
                                  estimator = sum(Delta2)
                                  
                                  estimator.var = sum((Delta2 - mean(Delta2) - self$mu.correct - self$v.mu.y.correct -  self$epsilon.mu.z.correct)^2)
                                  
                                 # estimator.var = sum((self$Delta - mean(self$Delta) - self$mu.correct - self$v.mu.y.correct -  self$epsilon.mu.z.correct)^2)
                                 
                                  return(c(estimator, estimator.var))
                             
                                },
                                cov.est.tilde.correct = function(){
                                  
                                  mu.y = mean(self$Y)
                                  
                                  mu.z = mean(self$Z)
                                  
                                  self$Delta = self$Y.hat.X.U * self$Z.hat.X.U/length(self$Y.hat.X.U) + 
                                    self$res.hat.epsilon * self$Z.hat.X.U/length(self$I.y) + 
                                    self$res.hat.v * self$Y.hat.X.U/length(self$I.z)
                                  
                                  estimator = sum(self$Delta) - mu.y*mu.z 
                                  
                                  #mu.correct = (Y.all *g.bar/self$N.Y + Z.all *f.bar/self$N.Z)

                                  #estimator.var = sum((self$Delta - mean(self$Delta) - mu.correct - estimator/length(self$I))^2)
                                  estimator.var = sum((self$Delta - mean(self$Delta) - self$mu.correct)^2)
                             
                                 # estimator.var = max(length(self$Delta)*var(self$Delta), 
                                 #                     log(self$N.Y)/self$N.Y^(1.5), 
                                  #                    log(self$N.Z)/self$N.Z^(1.5))
                                  
                                  return(c(estimator, estimator.var))
                                  
                                },
                                cov.est.hat.2 = function(){
                                
                                  mu.y = mean(self$Y.hat.X.U)
                                  
                                  mu.z = mean(self$Z.hat.X.U)
                                  
                                  
                                  self$res.hat.epsilon[self$I.y] = self$res.hat.epsilon[self$I.y] - mean(self$res.hat.epsilon[self$I.y])
                                  
                                  self$res.hat.v[self$I.z] = self$res.hat.v[self$I.z] - mean(self$res.hat.v[self$I.z]) 
                                  
                                  self$Delta = self$Y.hat.X.U * self$Z.hat.X.U/length(self$Y.hat.X.U) + 
                                    self$res.hat.epsilon * self$Z.hat.X.U/length(self$I.y) + 
                                    self$res.hat.v * self$Y.hat.X.U/length(self$I.z)
                                  
                                  Delta2 = (self$Y.hat.X.U - mu.y)  * (self$Z.hat.X.U - mu.z)/length(self$Y.hat.X.U) + 
                                    self$res.hat.epsilon * (self$Z.hat.X.U - mu.z) /length(self$I.y) + 
                                    self$res.hat.v * (self$Y.hat.X.U - mu.y)/length(self$I.z)
                                  
                                  estimator = sum(Delta2)
                                
                                  delta.1.hat = sum((self$Y.hat.X.U - mu.y)  * (self$Z.hat.X.U - mu.z)/length(self$Y.hat.X.U))/length(self$Y.hat.X.U) 
                                  
                                  delta.2.hat = sum(self$res.hat.epsilon * (self$Z.hat.X.U - mu.z) /length(self$I.y))/length(self$I.y)
                                  
                                  delta.3.hat = sum(self$res.hat.v * (self$Y.hat.X.U - mu.y)/length(self$I.z))/length(self$I.z)
                                  
                                  #######################################
                                  #
                                  # correction part
                                  #
                                  #########################################
                                  

                                  epsilon.mu.z.correct = rep(0, length = self$N.U)
                                  
                                  names(epsilon.mu.z.correct) = self$I
                                  
                                  epsilon.mu.z.correct[self$I.y] =  rep(delta.2.hat, length(self$I.y))
                                

                                  ###############################
                                  #
                                  # v.all
                                  #
                                  ################################
                                  
                                  v.mu.y.correct = rep(0, length = self$N.U)
                                  
                                  names(v.mu.y.correct) = self$I
                                  
                                  v.mu.y.correct[self$I.z] = rep(delta.3.hat, length(self$I.z))  
                                  
                                  #self$v.mu.y.correct = (v.all) *mu.y/self$N.Z + self$Y.hat.X.U/self$N.Y *(sum(self$res.hat.v)/length(self$I.z))
                           
                                  # estimator.var = sum((Delta2 - mean(Delta2) - self$mu.correct - self$v.mu.y.correct -  self$epsilon.mu.z.correct)^2)
                                  
                                  estimator.var = sum((Delta2 - delta.1.hat - epsilon.mu.z.correct  - v.mu.y.correct )^2)
                                  
                                #  plugest = self$plug() - mu.y*mu.z
                                  
                                  #estimator.var = length(Delta2)*var(Delta2)
                                
                                  
                                  # Temp = (self$Y.hat.X.U - mu.y.hat)*(self$Z.hat.X.U - mu.z.hat)/length(self$Y.hat.X.U) + 
                                  #   (self$Y.hat.X.U - mu.y.hat)*(mu.z.hat)/length(self$Y.hat.X.U) +
                                  #   (self$Z.hat.X.U - mu.z.hat)*(mu.y.hat)/length(self$Z.hat.X.U) -
                                  #   (self$Y.hat.X.U - mu.y.hat)*(self$I %in% self$I.y)*(mu.z.hat)/length(self$I.y)-
                                  #   (self$Z.hat.X.U - mu.z.hat)*(self$I %in% self$I.z)*(mu.y.hat)/length(self$I.z)+
                                  #   res.hat.epsilon * (self$Z.hat.X.U -mu.z.hat)/length(self$I.y) + 
                                  #   res.hat.v * (self$Y.hat.X.U - mu.y.hat)/length(self$I.z)
                                  # 
                                  # Temp2 = (self$Y.hat.X.U - mu.y.hat)*(self$Z.hat.X.U - mu.z.hat)/length(self$Y.hat.X.U) + 
                                  #   res.hat.epsilon * (self$Z.hat.X.U -mu.z.hat)/length(self$I.y) + 
                                  #   res.hat.v * (self$Y.hat.X.U - mu.y.hat)/length(self$I.z)
                                  
                                  return(c(estimator, estimator.var))
                                  
                                },
                                cov.est.hat.correct = function(){
                                  
                                  mu.y = mean(self$Y.hat)
                                  
                                  mu.z = mean(self$Z.hat)
                                  
                                  self$Delta = self$Y.hat.X.U * self$Z.hat.X.U/length(self$Y.hat.X.U) + 
                                    self$res.hat.epsilon * self$Z.hat.X.U/length(self$I.y) + 
                                    self$res.hat.v * self$Y.hat.X.U/length(self$I.z)
                                  
                                  estimator = sum(self$Delta) - mu.y*mu.z 
                                  
                                  estimator.var = max(length(self$Delta)*var(self$Delta), 
                                                      log(self$N.Y)/self$N.Y^(1.5), 
                                                      log(self$N.Z)/self$N.Z^(1.5))
                                  
                                  return(c(estimator, estimator.var))
                                  
                                },
                                cov.est.semi.hat.correct = function(){
                                  
                                  mu.y = mean(self$Y.hat.X.U) + sum(self$res.hat.epsilon)/self$N.Y
                                  
                                  mu.z = mean(self$Z.hat.X.U) + sum(self$res.hat.v)/self$N.Z
                                  
                                  self$Delta = self$Y.hat.X.U * self$Z.hat.X.U/length(self$Y.hat.X.U) + 
                                    self$res.hat.epsilon * self$Z.hat.X.U/length(self$I.y) + 
                                    self$res.hat.v * self$Y.hat.X.U/length(self$I.z)
                                  
                                  estimator = sum(self$Delta) - mu.y*mu.z 
                                  
                                  estimator.var = max(length(self$Delta)*var(self$Delta), 
                                                      log(self$N.Y)/self$N.Y^(1.5), 
                                                      log(self$N.Z)/self$N.Z^(1.5))
                                  
                                  return(c(estimator, estimator.var))
                                  
                                }
                                
)
)




