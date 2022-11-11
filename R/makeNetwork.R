#' Make contact graph
#' 
#' This function will create a contact network based on observed degree distributions
#' 
#' @param N integer. Size of network/number of nodes.
#' @param nTries integer. It sometimes takes a few tries before a network can be made which matches the 
#' degree distribution.
#' @return igraph object
#' @export
makeNetwork <-
function(N,nTries=10){
  A = NULL
  for(it in 1:nTries){
    cat("\n")
    cat(paste0("Attempt number ",it))
    cat("\n")
    try({
      N_adjusted = N/(1-dnbinom(0,mu=12.5,size=1/(1 - 3/12.5)))
      degs = rnbinom(N_adjusted,mu=12.5,size=1/(1 - 3/12.5))
      degs = degs[which(degs > 0)]
      A = sample_degseq(degs,method="simple") %>% simplify()
    },silent=TRUE)
    if(!is.null(A)) break
  }
  return(A)
}
