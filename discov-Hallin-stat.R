library(clue)
library(adagio)
library(energy)

do <- function(p) { x = rnorm(p); return(x/sqrt(sum(x^2))) }

do.grid <- function(n,p){
  if (p == 1) {
    nR = n %/% 2; n0 = n %% 2
    g = matrix(c(-(nR:1)/(nR+1), rep(0,n0), (1:nR)/(nR+1)), n)
  } else {
    nR = floor(n^0.25); nS = n %/% nR; n0 = n %% nR
    g = rbind(((1:nR)/(nR+1)) %x% t(replicate(nS, do(p))), matrix(0,n0,p))
  }
  return(g)
}

Hallin.rank <- function(X, const=1e4, constC=1e8){
  n = nrow(X); p = ncol(X)
  X.Hallin = do.grid(n,p)
  if (p == 1){
    X.rank = X.Hallin[rank(X),]
  } else {
    X.DIST = matrix(NA,n,n)
    for (u in 1:n){
      for (v in 1:n){
        X.DIST[u,v] = sum((X[u,] - X.Hallin[v,])^2)
      }
    }
    if (const*max(X.DIST) < constC) {
      X.rank = X.Hallin[assignment(round(const*X.DIST))$perm,]
    } else {
      X.rank = X.Hallin[as.vector(solve_LSAP(X.DIST)),]
    }
  }
  return(X.rank)
}

discov.Hallin.stat <- function(X, Y, const=1e4, constC=1e8) {
  if (nrow(X) != nrow(Y)) {
    print("Error: sample sizes must agree!")
  } else {
    n = nrow(X)
    X.rank = Hallin.rank(X, const, constC)
    Y.rank = Hallin.rank(Y, const, constC)
    return(n*dcovU(X.rank, Y.rank))
  }
}