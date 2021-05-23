#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @author : Romain Graux, Amandine Evrard, Céline Everaert
# @date : 2021 May 22, 12:10:38
# @last modified : 2021 May 23, 09:29:46

require("coda")

# Data declaration
flanders_data = c(25,69, 65,106, 80,106,136,94,76,46)
wallonia_data = c(17,36, 47,58 , 47,53 ,59 ,54,33,21)

data = matrix(data=c(flanders_data, wallonia_data), nrow=2)
rownames(data) <- c("Flanders", "Wallonia") 
colnames(data) <- c("<1200", "[1200,1500)", "[1500,1800)", "[1800,2300)", "[2300,2700)", "[2700,3300)", "[3300,4000)", "[4000,4900)", "[4900,6000)", ">6000")

interval<-c(0,1200,1500,1800,2300,2700,3300,4000,4900,6000,Inf)

muprior<-3000
sigmaprior<-300

# [3](b)

lpost <-function(theta,freq){
  lambda = 1/(theta[1]*theta[2])
  kappa = 1/theta[2]
  #Likelihood
  highlow = function(low,high){
    pgamma(high,kappa,lambda)-pgamma(low,kappa,lambda)
  }
  likelihood = sum(sapply(1:length(freq),function(j){
    freq[j] *log(highlow(interval[j],interval[j+1]))
  }))
  logpost = dnorm(theta[1],muprior,sigmaprior,log=T) + dunif(theta[2],0,10,log=T)
  
  lpost = likelihood+logpost
  return(lpost)
}

# [5](a)

componentwise_metropolis <- function(n_run, theta, sd.prop, frequencies, burnin=0.1){
  m = length(sd.prop)
  n_accepted = c(0, 0)
  walk = matrix(theta, ncol=m, byrow=TRUE)
  
  for (i in 2:(n_run+1)){
    current_theta = walk[i-1,]
    
    thetas = matrix(rep(current_theta, m), ncol=m, byrow = TRUE)
    alphas = rep(0, m)
    
    for (j in 1:m){
      thetas[j, j] = thetas[j, j] + rnorm(1, 0, sd.prop[j])
      alphas[j] = min(1, exp(lpost(thetas[j,], frequencies) - lpost(current_theta, frequencies)))
    }
    
    is_accepted = runif(m) <= alphas
    
    walk = rbind(walk, current_theta)
    walk[i, is_accepted] = diag(thetas)[is_accepted]
    
    n_accepted = n_accepted + as.integer(is_accepted)
  }
  
  walk = tail(walk, -burnin*n_run)
  accepted_rate = n_accepted / n_run
  
  return(list(walk=walk, accepted_rate=accepted_rate))
}


init_thetas <- c(2900, 0.3)
sd.prop <- c(174, 0.0475)
metropolis <- componentwise_metropolis(12200, init_thetas, sd.prop, data[1,])

sprintf("Acceptance rate for mu    in Flanders : %.3f", metropolis$accepted_rate[1])
sprintf("Acceptance rate for sigma in Flanders : %.3f", metropolis$accepted_rate[2])


plot(metropolis$walk[,1], metropolis$walk[,2], xlab = "mu", ylab="phi", main=paste("Acceptance rate -> mu:", accepted_rate[1], "phi:", accepted_rate[2]))
plot(as.mcmc(metropolis$walk))
print(paste("Effective size :", effectiveSize(walk)))
