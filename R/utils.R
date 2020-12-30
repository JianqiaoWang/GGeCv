#' @description  Some helper functions
#'

row.match <- function(x, table, nomatch=NA){
  
  if (class(table)=="matrix") table <- as.data.frame(table)
  if (is.null(dim(x))) x <- as.data.frame(matrix(x,nrow=1))
  cx <- do.call("paste",c(x[,,drop=FALSE],sep="\r"))
  ct <- do.call("paste",c(table[,,drop=FALSE],sep="\r"))
  match(cx,ct,nomatch=nomatch)
}

match_mat <- function(mat1, mat2){
  
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