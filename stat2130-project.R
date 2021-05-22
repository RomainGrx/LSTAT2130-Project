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

interval = c(0, 1200, 1500, 1800, 23002700, 3300, 4000, 4900, 6000, Inf)

muprior = 3000
sigmaprior = 300

kappa <- function(phi) 1/phi
lambda <- function(phi, mu) 1/(phi * mu)


# [3](b)

lpost = function(theta, freq) {
  lambda = 1 / (theta[1] * theta[2])
  kappa = 1 / theta[2]
  #Likelihood
  highlow = function(low, high) {
    pgamma(high, kappa, lambda) - pgamma(low, kappa, lambda)
  }
  likelihood = sum(sapply(1:length(freq), function(j) {
    freq[j] * log(highlow(interval[j], interval[j + 1]))
  }))
  logpost = dnorm(theta[1], muprior, sigmaprior, log = T) + dunif(theta[2], 0, 10, log = T)
  
  lpost = likelihood + logpost
  return(lpost)
}
