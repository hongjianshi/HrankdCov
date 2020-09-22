source("discov-Hallin-stat.R")

library(mvtnorm)
library(clue)
library(adagio)
library(energy)

n.l = c(216, 432, 864, 1728)
p.l = c(2, 3, 5, 7)
tau.l = c(0, 0.5, 0.9)
rho.l = seq(0,0.15,0.005)
repository = array(NA, dim = c(length(n.l),length(p.l),length(tau.l),length(rho.l),2,3))

for (n in n.l){
  for (p in p.l){
    for (tau in tau.l){
      for (rho in rho.l){
        M = diag(2*p)
        M[1,2] = M[2,1] = tau
        M[1,(p+1)] = M[(p+1),1] = rho
        
        X  = rmvnorm(n, sigma = M)
        X1 = X[,1:p]
        X2 = X[,(p+1):(2*p)]
        discov.Hallin = discov.Hallin.stat(X1, X2)
        discov.p      = dcov.test(X1, X2, R=n)$p
        discov.rank.p = dcov.test(apply(X1, 2, rank), apply(X2, 2, rank), R=n)$p
        repository[match(n,n.l),match(p,p.l),match(tau,tau.l),match(rho,rho.l),1, ] = 
          c(discov.Hallin, discov.rank.p, discov.p)
        
        X  = qt(pnorm(rmvnorm(n, sigma = M)), df=1)
        X1 = X[,1:p]
        X2 = X[,(p+1):(2*p)]
        discov.Hallin = discov.Hallin.stat(X1, X2)
        discov.p      = dcov.test(X1, X2, R=n)$p
        discov.rank.p = dcov.test(apply(X1, 2, rank), apply(X2, 2, rank), R=n)$p
        repository[match(n,n.l),match(p,p.l),match(tau,tau.l),match(rho,rho.l),2, ] = 
          c(discov.Hallin, discov.rank.p, discov.p)
      }
    }
  }
}