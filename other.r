mergeSeuratClusters <- function(x, clusters = NA, newName="merged1"){
  
  old <- as.character(x@ident)
  ind <-which(x@ident %in% clusters)
  old[ind] <- newName
  old <- as.factor(old)
  names(old) <- names(x@ident)
  x@ident <- old
  x
  
}

relabelCluster <- function(x, cluster =NA, newName=NULL){
  
  
  old <- as.character(x@ident)
  ind <-which(x@ident == cluster)
  old[ind] <- newName
  old <- as.factor(old)
  names(old) <- names(x@ident)
  x@ident <- old
  x
}