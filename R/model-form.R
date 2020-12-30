#' Define the potential model forms 
#'
expit = function(x){
  exp(x)/(1 + exp(x))
}

linear = function(x){
  return(x)
}
