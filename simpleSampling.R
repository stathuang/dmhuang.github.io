##############################
####### MCMC for catalytic priors
####### Dongming Huang (statdongminghuang@gmail.com)
####### date: 2020-07-15
##############################

library(mvtnorm)

## Metropolis-Hasting with proposals drawn from
## the independent normal distribution matching the mode and covariance
MH.ind.norm = function(sampSize, # output sample size
                       gapsize = 2, # gap between sub-samples
                       burnin = 10000, # number of steps for burn-in
                       mod, # mode; initial point
                       mat.cov, # covariance matrix for proposals
                       logp,    # function to evaluate log density
                       par      # parameters for logp
                       ) {
  p = length(mod)
  nprop = gapsize * sampSize + burnin
  mc = matrix(0, nprop, p)
  
  # proposals are normal vectors centered at mod with covariance matched
  normsamp = rnorm(p * nprop)
  normsamp[1:p] = 0
  eS = eigen(mat.cov)
  props = sweep(tcrossprod(
    matrix(normsamp, nrow = nprop, byrow = T),
    eS$vectors %*% diag(sqrt(pmax(eS$values, 0)))
  )
  , 2, mod, FUN = '+')
  
  dnormof = matrix(dnorm(normsamp, log = T), ncol = p, byrow = T)
  w = logp(props, par) - rowSums(dnormof)   # log importance ratio
  
  mc[1, ] = mod
  for (i in 2:(nprop)) {
    r = w[i] - w[i - 1]   # log MH ratio
    #calculate the accepted prob.
    if (log(runif(1)) < r) {
      #jump to the new one with the acceptance probabiity
      mc[i, ] = props[i, ]
    } else{
      mc[i, ] = mc[i - 1, ]
    }
  }
  mc = mc[-(1:burnin), ]
  return(mc[gapsize * (1:sampSize) - gapsize + 1, ])
}

library(mcmc)
# import Radford Neal's HMC code
path = getwd()
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source('RNHMC.R')
setwd(path)

# fn: call the HMC function to draw samples
MCHMC = function(N,    # number of samples
                 init, # initial point
                 eps = 0.01, # scales of coordinates
                 logp,       # fn to evaluate the log-density
                 logp.grad,  # fn to evaluate the gradient of log-density
                 par,        # parameters for the density functions
                 gapsize = 1,
                 L = 10,
                 strategy = strategy,
                 bMiniBatch = F,
                 mini.batch.size = 40,
                 first.pseudo = 31) {
  mc = matrix(0, N + 1, length(init))  # samples
  mc[1, ] = init
  mcprob = array(dim = N + 1)   # log of current probability
  accept = array(dim = N)
  
  oldp = logp(init, par)   #record the old likelihood
  mcprob[1] = oldp
  n = length(par$y)
  step.par = par
  
  if (strategy == 1) {
    scales = rep(1, N) # constant scale
  } else{
    scales = rbeta(N, 0.5, 0.5) # ranodm scale
  }
  for (t in 1:N)
  {
    if (bMiniBatch) {
      # use a subsample of pseudo data to compute the gradient
      ind = sample(first.pseudo:n, mini.batch.size)
      full.ind = c(1:(first.pseudo - 1), ind)
      step.par$x = par$x[full.ind, ]
      step.par$y = par$y[full.ind]
      # re-weight the pseudo sub-sample
      temp = par$weights[ind] * (n - first.pseudo + 1) / mini.batch.size
      step.par$weights = c(par$weights[1:(first.pseudo - 1)], temp)
    }
    newpoint = HMC(function(q) - logp(par = par, q), 
                   function(q) - logp.grad(par = step.par, q), 
                   scales[t] * eps, L, mc[t, ])
    mc[t + 1, ] = newpoint[[1]]
    mcprob[t + 1] = -newpoint[[2]]
    accept[t] = newpoint[[3]]
  }
  ind = 1 + seq(1, N, gapsize)
  return(list(
    samples = mc[ind, ],
    prob = mcprob[ind],
    accept = mean(accept)
  ))
}


# fn: interface for using HMC adaptively
############################################################
## init : initial point of the chain
## sd   : scales of updates in different coordinates
## logp, logp.grad, par  : log density, its gradient, and parameters
## BurnIn : number of iterations for burn-in
## gapsize: gap between sub-samples
## L : number of leapfrog steps
## strategy : 1 for constant scale, 2 for random scale
## bMiniBatch,mini.batch.size,first.pseudo: for catalytic prior only;
## whether to use a random subset (of size mini.batch.size) of
## pseudo sample to compute the gradient; if yes, first.pseudo is
## the first index of pseudo samples in the data set
############################################################
HMC.catalytic = function(sampSize,
                 init,
                 sd,
                 logp,
                 logp.grad,
                 par,
                 BurnIn = 5000,
                 gapsize = 10,
                 L = 10,
                 strategy = 1,
                 bMiniBatch = F,
                 mini.batch.size = 40,
                 first.pseudo = NULL) {
  eps = sd / L   # scales of coordinates
  
  # burn-in
  out <-
    MCHMC(
      logp = logp,
      init = init,
      N = BurnIn,
      eps = eps,
      logp.grad = logp.grad
      ,
      par = par,
      L = L,
      gapsize = 1,
      strategy = strategy,
      bMiniBatch,
      mini.batch.size,
      first.pseudo
    )
  N = sampSize * gapsize
  out2 <-
    MCHMC(
      logp = logp,
      init = out$samples[BurnIn, ],
      N = N,
      eps = eps,
      logp.grad = logp.grad
      ,
      par = par,
      L = L,
      gapsize = gapsize,
      strategy = strategy,
      bMiniBatch,
      mini.batch.size,
      first.pseudo
    )
  return(out2)
}
