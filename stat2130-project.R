#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @author : Romain Graux, Amandine Evrard, Céline Everaert
# @date : 2021 May 22, 12:10:38
# @last modified : 2021 May 23, 09:29:46

require("coda")

# Data declaration
flanders_data = c(25, 69, 65, 106, 80, 106, 136, 94, 76, 46)
wallonia_data = c(17, 36, 47, 58 , 47, 53 , 59 , 54, 33, 21)

data = matrix(data = c(flanders_data, wallonia_data),
              nrow = 2)
rownames(data) <- c("Flanders", "Wallonia")
colnames(data) <-
  c(
    "<1200",
    "[1200,1500)",
    "[1500,1800)",
    "[1800,2300)",
    "[2300,2700)",
    "[2700,3300)",
    "[3300,4000)",
    "[4000,4900)",
    "[4900,6000)",
    ">6000"
  )

interval <- c(0, 1200, 1500, 1800, 2300, 2700, 3300, 4000, 4900, 6000, 2**32)

muprior <- 3000
sigmaprior <- 306.12

# [3](b)

lpost <- function(theta, freq) {
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

# [5](a) Metropolis algorithm

componentwise_metropolis <-
  function(n_run,
           theta,
           sd.prop,
           frequencies,
           burnin = 0.1) {
    #' @param n_run. Number of samples in the chain.
    #' @param theta. The initial values for each component.
    #' @param sd.prop. The proposed standard deviations used in the normal distribution in order to get the next theta for each component.
    #' @param frequencies. The frequencies of the data.
    #' @param burnin. The ratio of first burn-in values to remove from the final chain.
    
    m = length(theta) # Number of variables
    n_accepted = c(0, 0)
    walk = matrix(theta, ncol = m, byrow = TRUE) # Contains the whole chain
    
    for (i in 2:(n_run + 1)) {
      current_theta = walk[i - 1, ] # Current theta at time t-1
      
      theta.props = matrix(rep(current_theta, m), ncol = m, byrow = TRUE) # Matrix with all component moves
      probs = rep(0, m) # Will contains the prob for each component
      
      # Compute for each component
      for (j in 1:m) {
        theta.props[j, j] = theta.props[j, j] + rnorm(1, 0, sd.prop[j])  # The proposed theta for component j with centered normal at the previous theta and particular sd
        probs[j] = min(1, exp(
          lpost(theta.props[j, ], frequencies) - lpost(current_theta, frequencies)
        )) # Get the prob for component j
      }
      
      is_accepted = runif(m) <= probs # Check if the probs are greater than random uniform
      
      walk = rbind(walk, current_theta) # Append the last theta on the walk
      
      walk[i, is_accepted] = diag(theta.props)[is_accepted] # Only move at the next theta if the probs are accepted
      
      n_accepted = n_accepted + as.integer(is_accepted) # Increment the number of accepted per component
    }
    
    walk = tail(walk,-burnin * n_run) # Remove the first burn-in values
    accepted_rate = n_accepted / n_run # Get the rate over all runs
    
    return(list(walk = walk, accepted_rate = accepted_rate))
  }


init_thetas <- c(2090, 0.4)
sd.prop <- c(500, 0.033)
metropolis <-
  componentwise_metropolis(25000, init_thetas, sd.prop, flanders_data)

sprintf("Acceptance rate for mu    in Flanders : %.3f",
        metropolis$accepted_rate[1])
sprintf("Acceptance rate for sigma in Flanders : %.3f",
        metropolis$accepted_rate[2])

plot(
  metropolis$walk[, 1],
  metropolis$walk[, 2],
  xlab = "mu",
  ylab = "phi",
  main = gettextf("Acceptance rate:: mu %.3f, phi %.3f", metropolis$accepted_rate[1], metropolis$accepted_rate[2])
)
plot(as.mcmc(metropolis$walk))
print(paste("Effective size :", effectiveSize(as.mcmc(metropolis$walk))))


# library(furrr)
# library(itertools)
# future::plan(multisession, workers=20)
#
# thetas.init = expand.grid(mu=seq(from=2500, to=3500, by=500), phi=seq(from=0.2, to=0.4, by=0.1))
# all_metropolis_runs = future_pmap(thetas.init, function(mu, phi){componentwise_metropolis(100, c(mu, phi), sd.prop, data[1,])})

# Gelman statistics
run1 = componentwise_metropolis(25000, c(2500, 0.3), sd.prop, flanders_data)
run2 = componentwise_metropolis(25000, c(2500, 0.5), sd.prop, flanders_data)
run3 = componentwise_metropolis(25000, c(3000, 0.3), sd.prop, flanders_data)
run4 = componentwise_metropolis(25000, c(3000, 0.5), sd.prop, flanders_data)
run5 = componentwise_metropolis(25000, c(3500, 0.3), sd.prop, flanders_data)
run6 = componentwise_metropolis(25000, c(3500, 0.5), sd.prop, flanders_data)

pmc = mcmc.list(
  as.mcmc(run1$walk),
  as.mcmc(run2$walk),
  as.mcmc(run3$walk),
  as.mcmc(run4$walk),
  as.mcmc(run5$walk),
  as.mcmc(run6$walk)
)

plot(pmc)

gelman.diag(pmc)
gelman.plot(pmc)

# Geweke diagnostic

for (obj in pmc){
  geweke.diag(obj)
  geweke.plot(obj)
}

# [5](c) Credible interval
# TODO : Need Laplace approximation bounds


#[6] JAGS implementation on Flanders data

require(R2WinBUGS)
require(runjags)
require(rjags)

model <- "model {
  kappa = 1 / phi
  lambda = 1 / (mu * phi) 
  
  mu ~ dnorm(3000, pow(306.12, -2))
  phi ~ dunif(0, 10)
  
  for (i in 1:10){
    pi[i] =  pgamma(x[i+1], kappa, lambda) - pgamma(x[i], kappa, lambda)
  }
  
  y ~ dmulti(pi, n)
}"

flanders_model = run.jags(model, c("mu", "phi"), burnin=2500, sample=25000, data=list(n=sum(flanders_data), x=interval, y=flanders_data), n.chains=5, inits=list(mu=init_thetas[1], phi=init_thetas[2]))
flanders_pcm = as.mcmc.list(flanders_model)

plot(flanders_pcm)

plot(
  x=as.double(unlist(flanders_model$mcmc[1][,1])),
  y=as.double(unlist(flanders_model$mcmc[1][,2])), 
  xlab = "mu",
  ylab = "phi",
)

# Gelman diagnostic
gelman.plot(flanders_pcm)
gelman.diag(flanders_pcm)

# Geweke diagnostic on first chain
geweke.plot(flanders_pcm[1])
geweke.diag(flanders_pcm[1])

write.jagsfile(flanders_model, "models/flanders.bug")

# [7] JAGS implementation on wallonia data

wallonia_model = run.jags(model, c("mu", "phi"), burnin=1000, sample=10000, data=list(n=sum(wallonia_data), x=interval, y=wallonia_data), n.chains=5, inits=list(mu=init_thetas[1], phi=init_thetas[2]))
wallonia_pcm = as.mcmc.list(wallonia_model)

plot(wallonia_pcm)

plot(
  x=as.double(unlist(wallonia_model$mcmc[1][,1])),
  y=as.double(unlist(wallonia_model$mcmc[1][,2])), 
  xlab = "mu",
  ylab = "phi",
)

# Gelman diagnostic
gelman.plot(wallonia_pcm)
gelman.diag(wallonia_pcm)

# Geweke diagnostic on first chain
geweke.plot(wallonia_pcm[1])
geweke.diag(wallonia_pcm[1])

write.jagsfile(wallonia_model, "models/wallonia.bug")
