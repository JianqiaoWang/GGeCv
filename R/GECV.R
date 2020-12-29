
library(R6)

row.match <- function(x, table, nomatch=NA){
  if (class(table)=="matrix") table <- as.data.frame(table)
  if (is.null(dim(x))) x <- as.data.frame(matrix(x,nrow=1))
  cx <- do.call("paste",c(x[,,drop=FALSE],sep="\r"))
  ct <- do.call("paste",c(table[,,drop=FALSE],sep="\r"))
  match(cx,ct,nomatch=nomatch)
}

match_mat <- function(mat1, mat2){
  
  # match row of m1 to m2, and return corresponding index in m2
  
  #library(data.table)
  
  #mat1 = setkey(data.table(mat1))
  
  #mat2 = setkey(data.table(mat2))
  
  #index = mat2[mat1 , which=TRUE]
  
  #apply(mat1, 2, function(x){
  
  #  which( x %in% mat2)
  
  #})
  
  m3 = rbind(mat1, mat2)
  
  index = duplicated(m3)[-c(1:nrow(mat1))]
  
  return(index)
  
}

mat_union = function(m1, m2){
  
  return(unique(rbind(m1,m2)))
  
}

mat_intersect = function(m1, m2){
  
  m3 <- rbind(m1, m2)
  
  intersect = m3[duplicated(m3), , drop = FALSE]
  
  return(intersect)
}

mat_setdiff = function(m1, m2){
  
  m3 <- rbind(m1, m2)
  
  setdiff = m3[!duplicated(m3), , drop = FALSE]
  
  return(setdiff)
}

#GeCov class

GeCv <-  R6Class("GeCov", list( X = "matrix",
                                W = "matrix",
                                beta = "vector",
                                gamma = "vector",
                                beta.hat = "matrix",
                                gamma.hat = "matrix",
                                Y = "vector",
                                Z = "vector",
                                f = "function",
                                g = "function",
                                Y.hat.X = NULL,
                                Z.hat.W = NULL,
                                Y.hat.X.U = NULL,
                                Z.hat.X.U = NULL,
                                res.hat.epsilon = NULL,
                                res.hat.v = NULL,
                                N.U = NULL,
                                N.Y = NULL,
                                N.Z = NULL,
                                Delta = NULL,
                                I.y = NULL,
                                I.z = NULL,
                                I = NULL,
                                mu.correct = NULL,
                                epsilon.mu.z.correct = NULL,
                                v.mu.y.correct = NULL,
                                initialize = function(Y, X,  beta.hat, f,  Z , W , gamma.hat, g
                                ){
                                  
                                  stopifnot(is.vector(as.vector(beta.hat)), is.vector(as.vector(gamma.hat)), length(beta.hat) == length(gamma.hat))
                                  
                                  stopifnot(is.numeric(Y), is.numeric(Z))
                                  
                                  self$beta.hat = as.vector(beta.hat)
                                  
                                  self$gamma.hat = as.vector(gamma.hat)
                                  
                                  self$f = f
                                  
                                  self$g = g
                                  
                                  self$X = X
                                  
                                  self$W = W
                                  
                                  self$Y = Y
                                  
                                  self$Z = Z
                                  
                                  self$N.Z = length(Z)
                                  
                                  self$N.Y = length(Y)
                                },
                                estimator =  function(){
                                  
                                  self$I.y =  rownames(self$X)
                                  
                                  self$I.z =  rownames(self$W)
                                  
                                  self$Y.hat.X <- self$f(self$X %*% self$beta.hat[-1] + rep(self$beta.hat[1], nrow(self$X)))
                                  
                                  Y.hat.W <- self$f(self$W %*% self$beta.hat[-1] + rep(self$beta.hat[1], nrow(self$W)))
                                  
                                  res.hat.X <- self$Y - self$Y.hat.X
                                  
                                  self$Z.hat.W <-  self$g(self$W %*% self$gamma.hat[-1] + rep(self$gamma.hat[1], nrow(self$W)))
                                  
                                  Z.hat.X <-  self$g(self$X %*% self$gamma.hat[-1] + rep(self$gamma.hat[1], nrow(self$X)))
                                  
                                  res.hat.W <- self$Z - self$Z.hat.W
                       
                                  X.U = unique(rbind(self$X, self$W))
                                  
                                  self$I = rownames(X.U)
                                  
                                  self$N.U = nrow(X.U)
                                  
                                  #############################################
                                  #
                                  # Calculate each component
                                  #
                                  #############################################
                                  
                                  self$Y.hat.X.U = self$f(X.U %*% self$beta.hat[-1] + rep(self$beta.hat[1], nrow(X.U)))
                                  
                                  self$Z.hat.X.U =  self$g(X.U %*% self$gamma.hat[-1] + rep(self$gamma.hat[1], nrow(X.U)))
                                  
                                  names(self$Y.hat.X.U) = self$I
                                  
                                  names(self$Z.hat.X.U) = self$I
                                
                                  self$res.hat.epsilon = rep(0, length = self$N.U)
                                  
                                  names(self$res.hat.epsilon) = self$I
                                  
                                  self$res.hat.epsilon[self$I.y] = res.hat.X
                         
                                  self$res.hat.v = rep(0, length = self$N.U)
                                  
                                  names(self$res.hat.v) = self$I
                                  
                                  self$res.hat.v[self$I.z] = res.hat.W
                                  
                                  ######################################
                                  # 
                                  # Correct part
                                  #
                                  ##################################
                                  
                                   mu.y.tilde = mean(self$Y)
                                  # 
                                   mu.z.tilde = mean(self$Z)
                                  
                                  Y.all = rep(0, length = self$N.U)
                                  
                                  names(Y.all) = self$I
                                  
                                  Y.all[self$I.y] = self$Y - mean(self$Y)
                                  
                                  Z.all = rep(0, length = self$N.U)
                                  
                                  names(Z.all) = self$I
                                  
                                  Z.all[self$I.z] = self$Z - mean(self$Z)
                                  
                                  #self$mu.correct = (Y.all) *mu.z.tilde/self$N.Y + (Z.all) *mu.y.tilde/self$N.Z
                                  
                                  
                                  ############################
                                  #
                                  # 
                                  #
                                  ############################
                                  
                                  #epsilon.all = rep(0, length = self$N.U)
                                  
                                  #names(epsilon.all) = self$I
                                  
                                  #epsilon.all[self$I.y] =  self$res.hat.epsilon[self$I.y] -  sum(self$res.hat.epsilon)/length(self$I.y)
                                  
                                  
                                  #self$epsilon.mu.z.correct = (epsilon.all) *mean(self$Z)/self$N.Y + (Z.all) *(sum(self$res.hat.epsilon)/length(self$I.y))/self$N.Z
                                  
                                  ###############################
                                  #
                                  # v.all
                                  #
                                  ################################
                                  
                                  # v.all = rep(0, length = self$N.U)
                                  
                                  #names(v.all) = self$I
                                  
                                  #v.all[self$I.z] =  self$res.hat.v[self$I.z] -  sum(self$res.hat.v)/length(self$I.z)
                                  
                                  #self$v.mu.y.correct = (v.all) *mean(self$Y)/self$N.Z + (Y.all)/self$N.Y *(sum(self$res.hat.v)/length(self$I.z))
                                  
                                  
                                  
                                  # 
                                  # 
                                  # mu.y.tilde = mean(self$Y)
                                  # 
                                  # mu.z.tilde = mean(self$Z)
                                  # 
                                  # estimator = sum(self$Delta) - mu.y.tilde*mu.z.tilde
                                  # 
                                  
                                  #print("Estimator")
                                  
                                  #print(estimator)
                                  
                                  #mu.y.hat = mean(self$Y.hat.X.U)
                                  
                                  #mu.z.hat = mean(self$Y)
                                  
                                 # mu.y.hat = mean(self$Y.hat.X.U)
                                  
                                  #mu.z.hat = mean(self$Z)
                                  # 
                                  # 
                                  # print("Large variance")
                                  # 
                                  # print(length(Temp)*var(Temp))
                                  # 
                                  # print("small variance:")
                                  # 
                                  # print(length(Temp2)*var(Temp2))
                                  # 
                                  # estimator.var = sum((Temp)^2)
                                  # 
                                  # print(estimator.var)
                                  # 
                                 # return(c(estimator, estimator.var))
                                  
                                 # Estimator = crossprod(self$Y.hat.X.U, self$Z.hat.X.U)/length(self$Y.hat.X.U) +
                                 #  crossprod(res.hat.W, Y.hat.W)/length(Y.hat.W) +
                                 #    crossprod(res.hat.X, Z.hat.X)/length(Z.hat.X)
                                },
                                plug = function(){
                                  
                                  return(crossprod(self$Y.hat.X.U, self$Z.hat.X.U)/length(self$Y.hat.X.U))
                                  
                                },
                                Inner.est = function(mu.y = 0, mu.z = 0){
                                  
                                  self$Delta = (self$Y.hat.X.U - mu.y) * (self$Z.hat.X.U - mu.z) /length(self$Y.hat.X.U) + 
                                    self$res.hat.epsilon * (self$Z.hat.X.U - mu.z) /length(self$I.y) + 
                                    self$res.hat.v * (self$Y.hat.X.U - mu.y)/length(self$I.z)
                                  
                                  estimator = sum(self$Delta) 
                                  
                                  estimator.var = max(length(self$Delta)*var(self$Delta), 
                                                      log(self$N.Y)/self$N.Y^(1.5), 
                                                      log(self$N.Z)/self$N.Z^(1.5)) 
                                  
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
                                cov.est.hat = function(){
                                
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
                                cov.est.semi.hat = function(){
                                  
                                  mu.y = mean(self$Y.hat.X.U) + sum(self$res.hat.epsilon)/self$N.Y
                                  
                                  mu.z = mean(self$Z.hat.X.U) + sum(self$res.hat.v)/self$N.Z
                                  
                                  Delta2 = (self$Y.hat.X.U - mu.y)  * (self$Z.hat.X.U - mu.z)/length(self$Y.hat.X.U) + 
                                    self$res.hat.epsilon * (self$Z.hat.X.U - mu.z) /length(self$I.y) + 
                                    self$res.hat.v * (self$Y.hat.X.U - mu.y)/length(self$I.z)
                                  
                                  estimator = sum(Delta2) 
                           
                                  estimator.var = length(Delta2)*var(Delta2)
                                  
                                  plugest = self$plug() - mu.y*mu.z
                                  
                                  return(c(estimator, estimator.var))
                                },
                                cov.est.hat.correct = function(){
                                  
                                  mu.y = mean(self$Y.hat.X.U)
                                  
                                  mu.z = mean(self$Z.hat.X.U)
                                  
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

expit = function(x){
  exp(x)/(1 + exp(x))
}

linear = function(x){
  return(x)
}

expn = function(x){
  return(x)
}



