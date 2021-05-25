#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @author : Romain Graux, Amandine Evrard, Céline Everaert
# @date : 2021 May 22, 12:10:38
# @last modified : 2021 May 23, 09:29:46

library(dae)
library(coda)
library(rjags)
library(mnormt)
library(ggplot2)
library(runjags)
library(EnvStats)
library(R2WinBUGS)

# Data declaration
muprior <- 3000
sigmaprior <- 306.12

flanders_frequencies = c(25, 69, 65, 106, 80, 106, 136, 94, 76, 46)
wallonia_frequencies = c(17, 36, 47, 58 , 47, 53 , 59 , 54, 33, 21)

interval <-
  c(0, 1200, 1500, 1800, 2300, 2700, 3300, 4000, 4900, 6000, 2 ** 32)

raw_flanders <-
  c(
    rep(600, 25),
    rep(1350, 69),
    rep(1650, 65),
    rep(2050, 106),
    rep(2500, 80),
    rep(3000, 106),
    rep(3650, 136),
    rep(4450, 94),
    rep(5450, 76),
    rep(7000, 46)
  )

data = matrix(data = c(flanders_frequencies, wallonia_frequencies),
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

# [3](b)

GammaFlandre <- egamma(raw_flanders)
kappaflandre <- GammaFlandre$parameters["shape"]
lambdaflandre <- 1 / GammaFlandre$parameters["scale"]
muflandre <- kappaflandre / lambdaflandre
phiflandre <- 1 / kappaflandre

highlow <- function(low, high, kappa, lambda) {
  pgamma(high, kappa, lambda) - pgamma(low, kappa, lambda)
}

lpost <- function(theta, freq) {
  lambda = 1 / (theta[1] * theta[2])
  kappa = 1 / theta[2]
  
  #Likelihood
  likelihood = sum(sapply(1:length(freq), function(j) {
    freq[j] * log(highlow(low = interval[j], high = interval[j + 1], kappa, lambda))
  }))
  logpost = dnorm(theta[1], muprior, sigmaprior, log = T) + dunif(theta[2], 0, 10, log = T)
  
  lpost = likelihood + logpost
  names(lpost) = "lpost"
  return(lpost)
}

lpost(c(muflandre, phiflandre), flanders_frequencies)

# [4]

laplace = function (mu, phi, freq) {
  ft = optim(
    list(mu = mu, phi = phi),
    lpost,
    control = list(fnscale = -1),
    hessian = T,
    freq = freq
  )
  params = ft$par
  cov = solve(-ft$hessian)
  echantillon = rmvnorm(params, cov, method = "choleski")
  list(echantillon = echantillon,
       params = params,
       cov = cov)
}

laplacefl = laplace(muprior, 0.01, flanders_frequencies)

posterior_laplace_flanders <-
  rnorm(50000, laplacefl$params["mu"], sqrt(laplacefl$cov[1, 1]))
HPD_laplace_flanders <-
  HPDinterval(as.mcmc(posterior_laplace_flanders), prob = 0.95)

plot(
  density(posterior_laplace_flanders),
  main = "Credible interval of Net Income in Flanders with Laplace approximation",
  xlab = parse(text = paste0('~ mu[1]'))
)
abline(v = HPD_laplace_flanders, col = 'orange')
legend(
  "topright",
  legend = c("HPD interval"),
  col = c("orange") ,
  lty = 1
)

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
    
    colnames(walk) = c("mu", "phi")
    return(list(walk = walk, accepted_rate = accepted_rate))
  }


mu1 <- parse(text = paste0("~mu[1]"))
phi1 <- parse(text = paste0("~phi[1]"))

# Run metropolis chain
init_thetas <- c(3000, 0.4)
sd.prop <- c(151, 0.033)
metropolis <-
  componentwise_metropolis(125000, init_thetas, sd.prop, flanders_frequencies)
metropolis_mcmc <- as.mcmc(metropolis$walk)

sprintf("Acceptance rate for mu    in Flanders : %.3f",
        metropolis$accepted_rate[1])
sprintf("Acceptance rate for sigma in Flanders : %.3f",
        metropolis$accepted_rate[2])

# Plot the chain
ggplot(as.data.frame(metropolis$walk), aes(x = mu, y = phi)) +
  ggtitle(
    gettextf(
      "Acceptance rate on Flanders data:: mu %.3f, phi %.3f",
      metropolis$accepted_rate[1],
      metropolis$accepted_rate[2]
    )
  ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()

plot(metropolis_mcmc)
print(paste("Effective size :", effectiveSize(as.mcmc(metropolis$walk))))

# Compute 6 chains from different initial values
run1 = componentwise_metropolis(25000, c(2500, 0.3), sd.prop, flanders_frequencies)
run2 = componentwise_metropolis(25000, c(2500, 0.5), sd.prop, flanders_frequencies)
run3 = componentwise_metropolis(25000, c(3000, 0.3), sd.prop, flanders_frequencies)
run4 = componentwise_metropolis(25000, c(3000, 0.5), sd.prop, flanders_frequencies)
run5 = componentwise_metropolis(25000, c(3500, 0.3), sd.prop, flanders_frequencies)
run6 = componentwise_metropolis(25000, c(3500, 0.5), sd.prop, flanders_frequencies)

# Group all chains
pmc = mcmc.list(
  as.mcmc(run1$walk),
  as.mcmc(run2$walk),
  as.mcmc(run3$walk),
  as.mcmc(run4$walk),
  as.mcmc(run5$walk),
  as.mcmc(run6$walk)
)

# Plot chains for both mu and phi
plot(pmc)

# Gelman diagnostic
gelman.plot(pmc)
gelman.diag(pmc)

# Geweke diagnostic
geweke.plot(metropolis_mcmc)
geweke.diag(metropolis_mcmc)

# [5](c) Credible interval
HPD_metropolis_flanders <- HPDinterval(metropolis_mcmc[, 1])
sprintf(
  "95%% Credible interval on metropolis chain : [%.2f, %.2f]",
  HPD_metropolis_flanders[1],
  HPD_metropolis_flanders[2]
)

# Plot the density of mu chain with credible interval
plot(density(metropolis_mcmc[, 1]),
     main = "Credible interval of Net Income in Flanders with metropolis output",
     xlab = parse(text = paste0('~ mu[1]')))
abline(v = HPD_laplace_flanders, col = 'orange', lty = 2)
abline(v = HPD_metropolis_flanders, col = 'purple', lty = 1)
legend(
  "topright",
  legend = c("Laplace HPD interval", "Metropolis HPD interval"),
  col = c("orange", "purple") ,
  lty = c(2, 1)
)


#[6] JAGS implementation on Flanders data
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

# Run jags with 25000 steps on 5 chains
flanders_model = run.jags(
  model,
  c("mu", "phi"),
  burnin = 2500,
  sample = 25000,
  data = list(
    n = sum(flanders_frequencies),
    x = interval,
    y = flanders_frequencies
  ),
  n.chains = 5,
  inits = list(mu = init_thetas[1], phi = init_thetas[2])
)
flanders_pcm = as.mcmc.list(flanders_model)
flanders_pcm_list <-
  list(mu = unlist(flanders_model$mcmc[][, 1]),
       phi = unlist(flanders_model$mcmc[][, 2]))

# Plot chains for both mu and phi
plot(flanders_pcm)

# PLot the density of all chains by hex
ggplot(as.data.frame(list(
  mu = unlist(flanders_model$mcmc[][, 1]),
  phi = unlist(flanders_model$mcmc[][, 2])
)), aes(x = mu, y = phi)) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()

# Gelman diagnostic
gelman.plot(flanders_pcm)
gelman.diag(flanders_pcm)

# Geweke diagnostic on first chain
geweke.plot(flanders_pcm[3])
geweke.diag(flanders_pcm[3])

# 95% credible interval of mu
HPD_jags_flanders <- unlist(HPDinterval(flanders_pcm[3][, 1]))
sprintf(
  "95%% Credible interval for Flanders on jags chain : [%.2f, %.2f]",
  HPD_jags_flanders[1],
  HPD_jags_flanders[2]
)

# Plot the density with credible interval
plot(density(unlist(flanders_pcm[3][, 1])),
     main = "Credible interval of Net Income in Flanders with jags output",
     xlab = parse(text = paste0('~ mu[1]')))
abline(v = unlist(HPD_jags_flanders),
       col = 'purple',
       lty = 1)
legend(
  "topright",
  legend = c("Laplace HPD interval"),
  col = c("purple") ,
  lty = c(1)
)

write.jagsfile(flanders_model, "models/flanders.bug")

# [7] JAGS implementation on wallonia data

# Run jags with 25000 steps on 5 chains
wallonia_model = run.jags(
  model,
  c("mu", "phi"),
  burnin = 2500,
  sample = 25000,
  data = list(
    n = sum(wallonia_frequencies),
    x = interval,
    y = wallonia_frequencies
  ),
  n.chains = 5,
  inits = list(mu = init_thetas[1], phi = init_thetas[2])
)
wallonia_pcm = as.mcmc.list(wallonia_model)
wallonia_pcm_list <-
  list(mu = unlist(wallonia_model$mcmc[][, 1]),
       phi = unlist(wallonia_model$mcmc[][, 2]))

# Plot chains for both mu and phi
plot(wallonia_pcm)

# PLot the density of all chains by hex
ggplot(as.data.frame(list(
  mu = unlist(wallonia_model$mcmc[][, 1]),
  phi = unlist(wallonia_model$mcmc[][, 2])
)), aes(x = mu, y = phi)) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()

# Gelman diagnostic
gelman.plot(wallonia_pcm)
gelman.diag(wallonia_pcm)

# Geweke diagnostic on first chain
geweke.plot(wallonia_pcm[5])
geweke.diag(wallonia_pcm[5])

# 95% credible interval of mu
HPD_jags_wallonia <- unlist(HPDinterval(wallonia_pcm[5][, 1]))
sprintf(
  "95%% Credible interval for wallonia on jags chain : [%.2f, %.2f]",
  HPD_jags_wallonia[1],
  HPD_jags_wallonia[2]
)

# Plot the density with credible interval
plot(density(unlist(wallonia_pcm[5][, 1])),
     main = "Credible interval of Net Income in wallonia with jags output",
     xlab = parse(text = paste0('~ mu[1]')))
abline(v = unlist(HPD_jags_wallonia),
       col = 'purple',
       lty = 1)
legend(
  "topright",
  legend = c("Laplace HPD interval"),
  col = c("purple") ,
  lty = c(1)
)

write.jagsfile(wallonia_model, "models/wallonia.bug")

# [8]

# Group all data
all_df <- data.frame(region = c(rep("Wallonia", 125000), rep("Flanders", 125000)),
                     mu = c(unlist(wallonia_pcm_list["mu"]), unlist(flanders_pcm_list["mu"])))

# Plot the density of density posterior for both regions
ggplot(all_df, aes(x = mu, fill = region)) +
  ggtitle("Density plot of Flanders and Wallonia mu chains") +
  geom_density(alpha = 0.4)

# Plot the difference in mean with credible and quantile intervals
diff_mean <-
  unlist(flanders_model$mcmc[][, 1]) - unlist(wallonia_model$mcmc[][, 1])
HDP_diff <- HPDinterval(as.mcmc(diff_mean))
plot(density(diff_mean),
     main = "Credible interval of the difference in Net Income between Flanders and Wallonia",
     xlab = parse(text = paste0('~ mu[1]', '-', '~ mu[2]')))
abline(v = quantile(diff_mean, probs = c(0.025, 0.975), col = 'purple'))
abline(v = HDP_diff, col = 'orange')
legend(
  "topright",
  legend = c("Quantile-based interval", "HPD interval"),
  col = c("purple", "orange") ,
  lty = 1:1
)
