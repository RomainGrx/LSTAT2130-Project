#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @author : Romain Graux, Amandine Evrard, Céline Everaert
# @date : 2021 May 22, 12:10:38
# @last modified : 2021 May 22, 12:11:31


interval = c(0, 1200, 1500, 1800, 23002700, 3300, 4000, 4900, 6000, Inf)
muprior = 3000
sigmaprior = 300

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
