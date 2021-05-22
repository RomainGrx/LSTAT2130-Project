#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @author : Romain Graux, Amandine Evrard, Céline Everaert
# @date : 2021 May 22, 12:10:38
# @last modified : 2021 May 22, 12:11:31


# Data declaration
flanders_data = c(25,69, 65,106, 80,106,136,94,76,46)
wallonia_data = c(17,36, 47,58 , 47,53 ,59 ,54,33,21)

data = matrix(data=c(flanders_data, wallonia_data), nrow=2)
rownames(data) <- c("Flanders", "Wallonia") 
colnames(data) <- c("<1200", "[1200,1500)", "[1500,1800)", "[1800,2300)", "[2300,2700)", "[2700,3300)", "[3300,4000)", "[4000,4900)", "[4900,6000)", ">6000")

interval<-c(0,1200,1500,1800,23002700,3300,4000,4900,6000,Inf)

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

componentwise_metropolis <- function(n_run, theta, frequencies, factors, standard_deviations){
  m <- length(theta)
  
  n_accepted = rep(0, m)
  walk = matrix(theta, ncol=m) 
  
  sigmas = factors * standard_deviations
  
  for (i in 2:(n_run+1)) {
    current_theta = tail(walk, 1)
    
    thetas = matrix(rep(current_theta, m), ncol=m, byrow=T)
    alphas = rep(0, m)
    
    for (j in 1:m){
      thetas[j, j] = thetas[j, j] + runif(1, 0, sigmas[j])
      delta = exp(lpost(thetas[j,], frequencies) - lpost(current_theta, frequencies))
      alphas[j] <- min(1, delta)
    }
    
    is_accepted <- runif(m) <= alphas
    
    next_theta = current_theta
    next_theta[is_accepted] = diag(thetas)[is_accepted]
    
    walk <- rbind(walk, next_theta) 
    n_accepted <- n_accepted + as.integer(is_accepted)
    
  }
  
  accepted_rate = n_accepted / n_run
  
  return(list(walk=walk, accepted_rate=accepted_rate))
  }
