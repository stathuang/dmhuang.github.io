##############################
####### Functions for Catalytic priors
####### Dongming Huang (statdongminghuang@gmail.com)
####### date: 2020-07-15
##############################


#Catalytic prior-----------
## Fn for generating covariates
resampling = function(x,  # original covariate matrix
                      size, # synethetic sample size
                      method = "default",
                      groups = NULL, # a list of vectors, each of which contains indices of a group
                      par = 1,
                      bMainEffect = F, # if T then only permute main effect; if F, all predictors will be permuted
                      dimen = NULL # the dimension of main effect
) {
  x = as.matrix(x)
  if (method == 'default') {
    if (bMainEffect) { # only permute the first dimen variables (the main effects)
      pseudox0 = apply(x[, 1:dimen], 2, sample, size = size, replace = T)
      pseudox.interact = t(apply(pseudox0, 1, combn, 2, prod))
      pseudox = cbind(pseudox0, pseudox.interact)
    } else{ #  default: resampling each coordinate independently
      pseudox = apply(x, 2, sample, size = size, replace = T) 
    }
  }else if (method == 'gaussian') { # sample from a normal dist. with matching moments
    pseudox = matrix(rnorm(size * ncol(x), sd = 1 * apply(x, 2, sd)), nrow = size)
  }else if (method == 'symmetry') {  # w.p. 1/2, flip around the mean
    pseudox = apply(x, 2, function(v)
      return(sample(c(v, 2 * mean(v) - v), size = size, replace = T))
    )
  }else if (method == 'flaten') { # w.p. 1/2, sample from uniform
    pseudox = apply(x, 2, function(v) {
      u = sample(v, size = size, replace = T)
      if (diff(range(scale(v))) > 5) {
        u[sample(size, size / 2)] = sample(unique(v), size = size / 2, replace =  T)
      }
      return(u)
    })
  }else if (method == 'unique') { # unform sample from the unique values 
    # -> avoid the issue of rare x
    pseudox = apply(x, 2, function(v)
      return(sample(
        unique(v), size = size, replace = T
      )))
  }  else if (method == 'extreme') { # only sample from max and min
    pseudox = apply(x, 2, function(v) 
      sample(c(max(v),min(v)), size, replace = T)
    )
  }else if (method == 'group_bootstrap') { # permuate each group together
    pseudox = matrix(nc = ncol(x), nr = size)
    n = nrow(x)
    for (ind in groups) {
      pseudox[, ind] = x[sample(1:n, size, rep = T), ind]
    }
  } else if (method == 'random_permute') { 
    # resample x, then randomly permute some coordinates (par)
    p = ncol(x)
    n = nrow(x)
    pseudox = x[sample(n, M, rep = T),]
    perm.ind = matrix(replicate(M, sample(p, par, rep = F)), ncol = par)
    for (i in 1:par) {
      pseudox[cbind(1:M, perm.ind[, i])] = x[cbind(sample(n, M, rep = T), perm.ind[, i])]
    }
  }
  else{
    print('Error: wrong method!') 
    return(NULL)
  }
  return(pseudox)
}
## fn: initialize a catalytic prior; generate pseudo data
InitCatalytic=function(M,x,y,p=ncol(x),pseudo.cov.dist,pseudo.model ){
  data.cata=list()
  n=length(y)
  
  data.cata$n=n
  data.cata$M=M
  data.cata$pseudo.cov.dist=pseudo.cov.dist
  data.cata$pseudo.model=pseudo.model
  
  x.star = resampling(x,M,method=pseudo.cov.dist)
  switch(pseudo.model,
         'equal'={
           y.star =rep(1/2,M)
         },
         'intercept'={
           y.star =rep((sum(y)+0.5)/(n+1),M)
         },
         'maineffect'={
           submodelmle=glm(y~x[,1:p],family=quasibinomial)$coefficients
           submodelmle[is.na(submodelmle)]=0
           fhi.star=x.star[,1:p]%*%submodelmle[-1] + submodelmle[1]
           
           if(bMixNull){
             y.star = (1/(1+exp(-fhi.star))+mean(y.obs))/2
           }else{
             fhi.star.new=linearization(fhi.star) # a method to avoid extreme values
             y.star = 1/(1+exp(-fhi.star.new)) 
           }
         }
  )
  
  
  data.cata$x.cata = cbind(1, rbind(x,x.star) )
  data.cata$y.cata = c(y,y.star)
  class(data.cata)='data.catalytic'
  return(data.cata)
}

## fn for fitting a logistic model with pseudo-sample and a given tau
CatalyticLogistic=function(data.cata,
                           tau, # prior weight
                           N.postSamp=5000,
                           bPostSamp=F, # whether to include posterior samples
                           bStan=F, # whether to compute posterior by Stan
                           N.chains.stan=5,
                           bMeasure = F, # whether to compute error
                           mu = NULL, 
                           fhi = NULL,
                           test.x = NULL,
                           test.y= NULL,
                           test.p= NULL,
                           criterion = "kl"  # error measurement
){
  
  Cata=list()
  Cata$selectedtau=tau
  n=data.cata$n
  M=data.cata$M
  
  den.par=list(x = data.cata$x.cata,y = data.cata$y.cata)
  den.par$weights=c(rep(1,n),rep(Cata$selectedtau/M,M))
  
  # compute posterior mode
  model.glm = glm(data.cata$y.cata ~ -1+data.cata$x.cata,weights = den.par$weights,family =  quasibinomial)
  Cata$betaHat=model.glm$coefficients
  
  
  # sampling from posterior
  if(bPostSamp | bMeasure){
    if(bStan){ # using Stan
      N.chains = N.chains.stan
      fit = stan(file='cata.stan', 
                 data=list(X=data.cata$x.cata[1:n,,drop=F],
                           Y=data.cata$y.cata[1:n],
                           n=n,M=M, 
                           Xstar=data.cata$x.cata[-(1:n),,drop=F],
                           Ystar=data.cata$y.cata[-(1:n)],
                           p=p+1,tau=8), 
                 warmup = BurnIn/N.chains, iter =(BurnIn+ N.postSamp) /N.chains,chains = N.chains)
      post.samp.cata=extract(fit)$beta
    }else{    # sampling using HMC
      intermed.cov = approxCov(model.glm,data.cata$x.cata,intercept = F)
      hmc.out = HMC.catalytic(   
        sampSize=N.postSamp,
        init=Cata$betaHat, 
        sd=sqrt( diag(intermed.cov) ),logp = logden,
        logp.grad = logden.grad,par=den.par,
        BurnIn = N.postSamp,gapsize = 1,L=5,strategy=2    
      )
      post.samp.cata=hmc.out$samples
    }
    
    if(bPostSamp){
      Cata$post.samp.cata=post.samp.cata
    }
  }
  
  # compute error; only if the true dist. is known
  if (bMeasure)
    Cata=c(Cata,Measurements.binary(  Cata$betaHat,data.cata$x.cata[1:n,-1,drop=F],data.cata$y.cata[1:n],mu,fhi,test.x,test.y,test.p,post.samp=post.samp.cata, criterion = criterion))
  
  return(Cata)
}


# fn for fitting a logistic regression with a full bayesian catalytic prior
CatalyticLogisticFullBayesian = function(data.cata,
                                         alpha=2, gammainv = 1, # hyperparameter alpha and 1/gamma
                                         N.postSamp = 5000,burnin = 5000, 
                                         Grid.tau,   
                                         init.tau=1,  # initial tau for Gibbs
                                         hmc.len = 5, hmc.steps = 5, # para. for HMC
                                         bOuputSamp=T, # whether to include posterior samples
                                         bTrun=T,  # whether to trim small values of tau 
                                         bVerb=T,  # whether to output progress
                                         bMeasure = F, # whether to compute errors
                                         mu,fhi,test.x,test.y,
                                         test.p,       # true test mean
                                         criterion = 'kl' # error measurement
) {
  Cata = list()
  n = data.cata$n
  M = data.cata$M
  x.obs = data.cata$x.cata[1:n,-1,drop=F]
  y.obs = data.cata$y.cata[1:n]
  x.star=data.cata$x.cata[-(1:n),,drop=F]
  y.star=data.cata$y.cata[-(1:n)]
  d = ncol(x.obs)
  
  
  ## find a lower bound for the tau-sequence
  if(bTrun){
    Plan.tau.local = Grid.tau
    post.cata.cov = list()
    idx.end = which.max(Plan.tau.local > 3)
    post.cata.mod = array(dim = c(idx.end, d + 1))
    md = array(data = 0, dim = idx.end)
    for (idx.tau in rev(1:idx.end))
    {
      tau.cata = Plan.tau.local[idx.tau]
      model.glm = glm(data.cata$y.cata ~ -1 + data.cata$x.cata,
                      weights = c(rep(1, n), rep(tau.cata / M, M)),
                      family =  quasibinomial)
      
      post.cata.mod[idx.tau,] = model.glm$coefficients#posterior mode
      W = model.glm$fitted.values * (1 - model.glm$fitted.values) * model.glm$prior.weights
      M1 = t(cbind(1, x.obs)) %*% diag(W[1:n]) %*% cbind(1, x.obs)
      M2 = t(x.star) %*% diag(W[-(1:n)]) %*% x.star
      post.cata.cov[[idx.tau]] = M1 + M2
      
      if (idx.tau < idx.end) {
        W = model.glm$fitted.values * (1 - model.glm$fitted.values) * model.glm$prior.weights
        M1 = t(cbind(1, x.obs)) %*% diag(W[1:n]) %*% cbind(1, x.obs)
        M2 = t(x.star) %*% diag(W[-(1:n)]) %*% x.star
        post.cata.cov[[idx.tau]] = M1 + M2
        u = post.cata.mod[idx.tau,] - post.cata.mod[idx.tau + 1,]
        md[idx.tau] = t(u) %*% post.cata.cov[[idx.tau + 1]] %*% u
      }
    }
    md = md[-idx.end]
    
    d2lmd = md#diff(log(md),differ=2)
    for (idx1 in rev(2:which.max(Plan.tau.local > 1))) {
      ref = d2lmd[(idx1):min(length(d2lmd), idx1 + 10)]
      threshold = ref[1] + 1.96 * sd(ref)
      if (d2lmd[idx1 - 1] > threshold)
        break
    }
    lowerb = idx1
    
    lowerbound = Plan.tau.local[lowerb]
    if (lowerb > 1) {
      Plan.tau.lim = Plan.tau.local[-(1:(lowerb - 1))]
    } else
      Plan.tau.lim = Plan.tau.local
    max.cov = post.cata.cov[[lowerb]]
    
    upperbound = max(Plan.tau.lim)
    Cata$bounds = c(lowerbound, upperbound)
    Grid.tau = Plan.tau.lim
    if(bVerb)print('Lower bound found.')
  }
  
  
  N.gibbs = N.postSamp + burnin   # number of Gibbs iteration
  gibbs.sample.tau = array(dim = N.gibbs)
  gibbs.accepted = array(dim = N.gibbs)
  gibbs.sample.beta = array(dim = c(N.gibbs, d + 1))
  
  # initialization
  gibbs.sample.tau[1] = init.tau
  model.glm = glm(data.cata$y.cata ~ -1 + data.cata$x.cata,
                  weights = c(rep(1, n), rep(gibbs.sample.tau[1] / M, M)),
                  family =  quasibinomial)
  init.cov = approxCov(model.glm, data.cata$x.cata, intercept = F)
  gibbs.sample.beta[1,] = model.glm$coefficients
  
  # kappa: the limit of C_tau
  null = glm(y.star ~ -1 + x.star, family = quasibinomial)
  kappa = logden(null$coefficients, list(x = x.star, y = y.star, weights =
                                           1 / M))
  
  ## preparation for Gibbs  ----------------------
  ## the space of tau is discretized into Grid.tau
  ## the HMC for updating beta is sensitive to scale parameters
  ## the idea is, if Grid.tau is not too many, we can find good scal parameters for HMC steps first
  
  
  Ftau = log(Grid.tau) * ( alpha + (d + 1) - 1) -
    Grid.tau * ( gammainv + kappa)
  
  ## this is very slow 
  proposal.scale = array(dim = c(length(Grid.tau), d + 1))
  proposal.beta = proposal.scale
  for (idx.tau in 1:length(Grid.tau)) {
    if (idx.tau %% 4 == 1) {
      model.glm = glm(data.cata$y.cata ~ -1 + data.cata$x.cata,
                      weights = c(rep(1, n), rep(Grid.tau[idx.tau] / M, M)),
                      family =  quasibinomial)
      intermed.cov = approxCov(model.glm, data.cata$x.cata, intercept = F)
      proposal.scale[idx.tau,] = sqrt(diag(intermed.cov))
      proposal.beta[idx.tau,] = model.glm$coefficients
    } else{
      proposal.scale[idx.tau,] =  proposal.scale[idx.tau - 1,]#diag(intermed.cov)
      proposal.beta[idx.tau,] =  proposal.beta[idx.tau - 1,]#model.glm$coefficients
    }
  }
  
  ##  Gibbs -----------------
  
  if(bVerb)print('Gibbs begins.')
  
  den.par = list(x = data.cata$x.cata, y = data.cata$y.cata)
  
  for (epo in 2:N.gibbs) {
    
    #update tau
    beta = gibbs.sample.beta[epo - 1,]
    theta = x.star %*% beta
    btheta = log(1 + exp(theta))
    btheta[is.infinite(btheta) == 1] = theta[is.infinite(btheta) == 1]
    vllh = theta * y.star - btheta
    cp.tau = Grid.tau * mean(vllh) + Ftau
    cp.tau = exp(cp.tau - max(cp.tau))
    sample.tau.id = sample(1:length(Grid.tau), 1, prob = cp.tau)
    gibbs.sample.tau[epo] = Grid.tau[sample.tau.id]
    
    #update beta
    den.par$weights = c(rep(1, n), rep(gibbs.sample.tau[epo] / M, M))
    hmc.out = HMC.catalytic(
      sampSize = hmc.len,
      init=gibbs.sample.beta[epo-1,] , 
      sd=proposal.scale[sample.tau.id,],
      logp = logden,
      logp.grad = logden.grad,
      par=den.par,
      BurnIn = 2,
      gapsize = 1,
      L = hmc.steps,
      strategy = 2   
    )
    gibbs.sample.beta[epo,] = hmc.out$samples[hmc.len,]
    gibbs.accepted[epo] = mean(hmc.out$accept)
    
  }
  if(bVerb)print('Gibbs finished.')
  
  post.samp.beta = gibbs.sample.beta[(N.gibbs - N.postSamp + 1):N.gibbs,]
  Cata$selectedtau = gibbs.sample.tau[(N.gibbs - N.postSamp + 1):N.gibbs]
  Cata$gibbs.accepted = gibbs.accepted
  
  if(bOuputSamp){
    Cata$betasample = post.samp.beta
  }
  Cata$betaHat = apply(post.samp.beta, 2, median)
  Cata$betaHat.trimmean = apply(post.samp.beta, 2, mean, trim = 0)
  
  # compute error; only if true dist. is known
  if(bMeasure){
    Cata$error.median=Measurements.binary(Cata$betaHat,x.obs,y.obs,mu,fhi,test.x,test.y,test.p,post.samp=post.samp.beta,criterion=criterion)
    Cata$error.trimmean=Measurements.binary(Cata$betaHat.trimmean,x.obs,y.obs,mu,fhi,test.x,test.y,test.p,post.samp=post.samp.beta,criterion=criterion)
  }
  
  return(Cata)
}

# Steinian estimate for risk of binary regression
EstimatePredErr.Stein = function(tau.seq,     # values of tau to be assessed
                                 data.cata,  # merged data set
                                 var.tau = 1  # tau_0 for estimating variance
) {
  
  n=data.cata$n
  M=data.cata$M
  
  # fit a model with tau_0 for variance estimate
  model.glm0 = glm(data.cata$y.cata ~ -1 + data.cata$x.cata,
                   weights = c(rep(1, n), rep(var.tau / M,M)),
                   family =  quasibinomial) # use for estimating Var_*
  w0 = with(model.glm0, fitted.values * (1 - fitted.values))
  
  est.stein=rep(NA, length(tau.seq))
  for (id.tau in 1:length(tau.seq)) {
    tau=tau.seq[id.tau]
    
    # fit a model with the given tau
    model.glm.tau = glm(data.cata$y.cata ~ -1 + data.cata$x.cata,
                        weights = c(rep(1, n), rep(tau / M, M)),
                        family =  quasibinomial)
    
    y.cata=data.cata$y.cata
    
    # compute the Steinian estiamte for covariance penalty
    cov.penalty.tau.Stein = 0
    for (i in 1:n) {
      y.cata[i] = 1 - y.cata[i]
      model.glm.tau.flipk = glm(y.cata ~ -1 + data.cata$x.cata,
                                weights = c(rep(1,n), rep(tau / M, M)),
                                family =  quasibinomial)
      y.cata[i] = 1 - y.cata[i]
      cov.penalty.tau.Stein = cov.penalty.tau.Stein + 
        w0[i] * (2 * y.cata[i] -1) * tcrossprod( 
          data.cata$x.cata[i, , drop = F], 
          (model.glm.tau$coefficient - model.glm.tau.flipk$coefficient)
        )
    }
    
    # compute negative log likelihood
    obs.par = list(x = data.cata$x.cata[1:n,,drop=F],
                   y = data.cata$y.cata[1:n],
                   weights = rep(1, n))
    neg.ll = -logden(model.glm.tau$coefficients, obs.par)
    
    est.stein[id.tau] = neg.ll + cov.penalty.tau.Stein
  }
  return(est.stein)
}


# Bootstrap estimate for risk of binary regression
EstimatePredErr.boot = function(tau.seq,    # values of tau to be assessed
                                data.cata,
                                var.tau = 1, # tau_0 for estimating variance
                                B = 100   # bootstrap replication
) {
  
  n =  data.cata$n
  M =  data.cata$M
  obs.par = list(x = data.cata$x.cata[1:n,,drop=F],
                 y = data.cata$y.cata[1:n],
                 weights = rep(1, n))
  y.star = data.cata$y.cata[-(1:n)]
  
  Y.boot     = array(dim = c(length(tau.seq), n, B))
  Fhihat.boot = array(dim = c(length(tau.seq), n, B))
  
  # fit a model for bootstrap sampling
  model.glm0 = glm( data.cata$y.cata ~ -1 + data.cata$x.cata,
                   weights = c(rep(1, n), rep(var.tau / M, M)),
                   family =  quasibinomial) 
  
  fhi.boot = obs.par$x %*% model.glm0$coefficients
  mu.boot = 1 / (1 + exp(-fhi.boot))
   
  for (b in 1:B) {
    # draw bootstrap sample
    y.boot = rbinom(n, 1, mu.boot)
    for (id.tau in 1:length(tau.seq)) {
      ## fit the model with a given tau using bootstrap samples
      tau = tau.seq[id.tau] 
      model.glm.tau.boot = glm(
        c(y.boot, y.star) ~ -1 + data.cata$x.cata,
        weights = c(rep(1, n), rep(tau / M, M)),
        family =  quasibinomial
      )
      
      Y.boot[id.tau, , b] = y.boot
      Fhihat.boot[id.tau, , b] = obs.par$x %*% model.glm.tau.boot$coefficients
      
    }
  }
  
  ## compute the negative loglikelihood and 
  ## the covariance penalty
  neg.ll = rep(NA, length(tau.seq))
  cov.penalty.tau.boot = rep(NA, length(tau.seq))
  est.boot = rep(NA, length(tau.seq))
  for (id.tau in 1:length(tau.seq)) {
    tau = tau.seq[id.tau]
    model.glm.tau = glm(data.cata$y.cata ~ -1 + data.cata$x.cata,
                        weights = c(rep(1, n), rep(tau / M, M)),
                        family =  quasibinomial)
    
    neg.ll[id.tau] = -logden(model.glm.tau$coefficients, obs.par)
    cov.penalty.tau.boot[id.tau] = 
      sum(
      sweep(Y.boot[id.tau, ,], 1, rowMeans(Y.boot[id.tau, ,]))  *  
                                         Fhihat.boot [id.tau, ,]
      ) / B
    est.boot[id.tau] = neg.ll[id.tau] + cov.penalty.tau.boot[id.tau]
  }
  
  
  return(est.boot)
}


# MLE
RunMLELogistic = function(x.obs,
                          y.obs,
                          N.postSamp=5000,
                          bMeasure = F,
                          mu,
                          fhi,
                          test.x,
                          test.y,
                          test.p,
                          criterion = "kl") {
  output = list()
  d = ncol(x.obs)
  model.pure.glm =   glm(y.obs ~ x.obs, family = quasibinomial)
  
  if (sum(is.na(model.pure.glm$coefficients)) == 0) {
    output$count = 1
    
    
    output$model = model.pure.glm
    output$betaHat =  model.pure.glm$coefficients
    output$betaHat[is.na(output$betaHat)] = 0
    post.cov = approxCov(model.pure.glm, sim.x)
    post.samp.mle = MH.ind.norm(
      N.postSamp,
      gapsize  = 2,
      burnin=10000,
      mod=output$betaHat,
      mat.cov=post.cov,
      logp = logden,
      par=list(
        x = cbind(1, x.obs),
        y = y.obs,
        weights = model.pure.glm$prior.weights
      )
    )
    
    # compute error; only if true dist. is known
    if (bMeasure)
      output = c(
        output,
        Measurements.binary(
          output$betaHat,x.obs,y.obs,
          mu, fhi,test.x,test.y,test.p,
          post.samp = post.samp.mle,
          criterion=criterion
        )
      )
    
  } else
    output$count = 0
  return(output)
}


# Cauchy prior
RunCauchy = function(x.obs,
                     y.obs,
                     N.postSamp=5000,
                     bSamp = F,
                     bMeasure=F,
                     mu,
                     fhi,
                     test.x,
                     test.y,
                     test.p,
                     criterion = "kl") {
  output = list()
  
  model.cauchy.glm = bayesglm(y.obs ~ x.obs, family = binomial)

  if (sum(is.na(model.cauchy.glm$coefficients)) == 0) {
    output$count = 1
    output$betaHat = model.cauchy.glm$coefficients
    
    
   
    post.samp.cauchy = sim(model.cauchy.glm, n.sims = N.postSamp)
    post.samp.cauchy = post.samp.cauchy@coef
    
    # compute error; only if true dist. is known
    if(bMeasure){
      output=c(output,Measurements.binary(output$betaHat,x.obs,y.obs,mu,fhi,test.x,test.y,test.p,post.samp = post.samp.cauchy,criterion=criterion))
    }
    
    if (bSamp)
      output$post.samp = post.samp.cauchy
  } else
    output$count = 0
  return(output)
}

# Experiment Setup------
AssignRegCoef = function(pattern,
                         p,
                         frac = 1,
                         seed = 1) {
  switch(
    pattern,
    multilevel_piecewiseconst_zero = {
      beta.true = c(rep(1, floor((p) / 2 + 1)), rep(-1, floor((p + 1) / 2)))
      beta.true = c(beta.true, rep(0, p * (p - 1) / 2))
      beta.true[1] = 0#0.2
    },
    multilevel_piecewiseconst_sparse_weak = {
      beta.true1 = c(rep(1, floor((p) / 2 + 1)), rep(-1, floor((p + 1) / 2)))
      beta.true2 = rep(1, (p - 1))
      beta.true2[2 * 1:floor((p - 1) / 2)] = -1
      beta.true2 = c(beta.true2, rep(0, (p - 2) * (p - 1) / 2))
      beta.true = c(beta.true1, beta.true2 * 0.8)
      beta.true[1] = 0#0.2
    },
    multilevel_piecewiseconst_dominate = {
      beta.true1 = c(rep(1, floor((p) / 2 + 1)), rep(-1, floor((p + 1) / 2)))
      beta.true2 = rep(1, (p - 1))
      beta.true2[2 * 1:floor((p - 1) / 2)] = -1
      beta.true2 = c(beta.true2, rep(0, (p - 2) * (p - 1) / 2))
      beta.true = c(beta.true1, beta.true2 * 1.0)
      beta.true[1] = 0#0.2
    },
    multilevel_piecewiseconst_diffuse = {
      beta.true1 = c(rep(1, floor((p) / 2 + 1)), rep(-1, floor((p + 1) / 2)))
      beta.true2 = rep(1, (p - 1) * p / 2)
      beta.true2[2 * 1:floor((p - 1) * p / 4)] = -1
      beta.true2 = beta.true2 * (1 / (1:((p - 1) * p / 2)))
      beta.true = c(beta.true1, beta.true2 * 0.4)
      beta.true[1] = 0#0.2
    },
    multilevel_mixture_mass = {
      main.frac  = 0.5
      indicator   = rep(0, p)
      indicator[sample(p, p * main.frac, rep = F)] = 1
      beta.true1 = sample(c(-1, 1), p, rep = T) *  indicator # rbinom(p,1,  main.frac)
      indicator   = rep(0, p * (p - 1) / 2)
      indicator[sample(p * (p - 1) / 2, p * (p - 1) / 2 * frac, rep = F)] =
        1
      beta.true2 = sample(c(-1, 1), p * (p - 1) / 2, rep = T) * indicator # rbinom(p*(p-1)/2,1,  frac)
      beta.true = c(0, beta.true1, beta.true2)
    },
    multilevel_mixture_norm = {
      main.frac  = 0.5
      indicator   = rep(0, p)
      indicator[sample(p, p * main.frac, rep = F)] = 1
      beta.true1 = rnorm(p) *  indicator # rbinom(p,1,  main.frac)
      indicator   = rep(0, p * (p - 1) / 2)
      indicator[sample(p * (p - 1) / 2, p * (p - 1) / 2 * frac, rep = F)] =
        1
      beta.true2 = rnorm(p * (p - 1) / 2) * indicator # rbinom(p*(p-1)/2,1,  frac)
      beta.true = c(0, beta.true1, beta.true2)
    },
    piecewiseconst = {
      beta.true = c(rep(1, floor((p) / 2 + 1)), rep(-1, floor((p + 1) / 2)))
      beta.true[1] = 0.2
    },
    alternatingconst = {
      beta.true = rep(1, (p + 1))
      beta.true[2 * 1:floor((p + 1) / 2)] = -1
      beta.true[1] = 0.2
    },
    piecewisei = {
      beta.true = c(rep(1, floor((p) / 2 + 1)), rep(-1, floor((p + 1) / 2)))
      beta.true = beta.true * (1 / (1:(p + 1)))
      beta.true[1] = 0.2
    },
    alternatingi = {
      beta.true = rep(1, (p + 1))
      beta.true[2 * 1:floor((p + 1) / 2)] = -1
      beta.true = beta.true * (1 / (1:(p + 1)))
      beta.true[1] = 0.2
    },
    piecewise = {
      beta.true = c(rep(1, floor((p) / 2 + 1)), rep(-1, floor((p + 1) / 2)))
      beta.true = beta.true * (1 / sqrt(1:(p + 1)))
      beta.true[1] = 0.2
    },
    alternating = {
      beta.true = rep(1, (p + 1))
      beta.true[2 * 1:floor((p + 1) / 2)] = -1
      beta.true = beta.true * (1 / sqrt(1:(p + 1)))
      beta.true[1] = 0.2
    },
    sparse = {
      beta.true = c(0.2, -1.1, -0.9, rep(0, p - 2))#too sparse
    },
    peak = {
      #  beta.true = c(rep(1,floor((p-1) / 2 )),rep(-1,floor((p) / 2)))
      # old:
      beta.true = c(rep(1, floor((p) / 2)), rep(-1, floor((p - 1) / 2)))
      beta.true = c(0.2, 4, beta.true)#too sparse
    },
    piecewiseisq = {
      beta.true = c(rep(1, floor((p) / 2 + 1)), rep(-1, floor((p + 1) / 2)))
      beta.true = beta.true * (1 / (1:(p + 1)) ^ 2)
      beta.true[1] = 0.2
    },
    alternatingisq = {
      beta.true = rep(1, (p + 1))
      beta.true[2 * 1:floor((p + 1) / 2)] = -1
      beta.true = beta.true * (1 / (1:(p + 1)) ^ 2)
      beta.true[1] = 0.2
    }  ,
    mixture_norm = {
      seed = 1
      indicator   = rep(0, p)
      indicator[sample(p, p * frac, rep = F)] = 1
      beta.true = c(0, rnorm(p) * indicator)
    },
    mixture_mass = {
      seed = 1
      indicator   = rep(0, p)
      indicator[sample(p, p * frac, rep = F)] = 1
      beta.true = c(0, sample(c(-1, 1), p, rep = T) * indicator)
    },
    random_fraction = {
      seed = 1
      s = floor((p) * frac)
      ind = sample(1 + 1:(p), s)
      beta.true = rep(0, (p + 1))
      beta.true[ind] = rnorm(s)
      beta.true[1] = 0.2
    },
    fixed = {
      beta.true = c(0.2, 1.1, -0.9, rep(0, p - 2))
    }
  )
  beta.true
}





#measurement of error in linear regression
MyPredErr = function(beta, beta.true, x) {
  mean((x %*% (beta - beta.true)) ^ 2 / 2)
}

# compute binomial deviance
CappedBinomialDeviance <- function(a, p, criterion = 'kl') {
  ### a is the ground true,p is the fitted
  if (length(a) !=  length(p))
    stop("Actual and Predicted need to be equal lengths!")
  p_capped <- pmin(0.999, p)
  p_capped <- pmax(0.001, p_capped)
  a_capped <- pmax(0.001, pmin(0.999, a))
  
  switch(
    criterion,
    "none" = {
      return(-mean(a * log(p_capped) + (1 - a) * log((1 - p_capped))))
    },
    "kl" = {
      return(-mean(a * log(p_capped / a_capped) + (1 - a) * log((1 - p_capped) /
                                                                  (1 - a_capped))))
      
    },
    "over-entropy" = {
      return(-mean(a * log(p_capped) + (1 - a) * log((1 - p_capped))) /
               -mean(a * log(a_capped) + (1 - a) * log((1 - a_capped))))
    }
  )
  
}
BinaryEntropy = function(p) {
  p_capped <- pmin(0.999, p)
  p_capped <- pmax(0.001, p_capped)
  mean(-p * log(p_capped) - (1 - p) * log(1 - p_capped))
}
KL = function(p, q)
{
  pcap = pmax(0.0001, pmin(0.9999, p))
  qcap = pmax(0.0001, pmin(0.9999, q))
  mean(q * log(qcap / pcap) + (1 - q) * log((1 - qcap) / pcap))
}

PredictiveInterval = function(samp,
                              x,
                              mu,
                              link = logit,
                              nominal = 0.025) {
  if (ncol(x) < ncol(samp))
    x = cbind(1, x)
  pred.mu = (tcrossprod(samp, x))
  lower = logit(apply(pred.mu, 2, function(v)
    quantile(v, nominal / 2)))
  upper = logit(apply(pred.mu, 2, function(v)
    quantile(v, 1 - nominal / 2)))
  result = list(
    coverage = mean((lower - mu) * (upper - mu) < 0),
    len = mean(upper - lower),
    loss = mean(upper - lower) + 2 / nominal * mean(pmax(0, lower -
                                                           mu, mu - upper))
  )
  result
}

# measurement of error in logistic regression
MyPredErrLogistic = function(beta,
                             test.x,
                             test.y,
                             test.p,
                             criterion = 'expdev',
                             ...) {
  theta = cbind(1, test.x) %*% beta
  pred.y = theta > 0
  pred.p = 1 / (1 + exp(-theta))
  if (criterion == 'expdev')
    return(CappedBinomialDeviance(test.p, pred.p, ...))
  if (criterion == 'kl')
    return(KL(test.p, pred.p))
}
Measurements.binary = function(betaHat,
                               x.obs,
                               y.obs,
                               mu,
                               fhi,
                               test.x,
                               test.y,
                               test.p,
                               post.samp = NULL,
                               criterion='kl') {
  measure = list()
  
  theta = cbind(1, x.obs) %*% betaHat
  
  # prediction
  measure$negloglike = mean(log(1 + exp(theta)) - theta * fhi)
  measure$pred.err = MyPredErrLogistic(betaHat, test.x, test.y, test.p, criterion =
                                         criterion)
  measure$pred.err.in = MyPredErrLogistic(betaHat, x.obs, y.obs, mu, criterion =
                                            criterion)
  
  
  measure$pred.err01 = MyPredErrBinary(betaHat, test.x, test.y, test.p, bProb = F)
  measure$pred.err01.in = MyPredErrBinary(betaHat, x.obs, y.obs, mu, bProb = F)
  
  measure$auc = MyPredErrBinary(betaHat,
                                test.x,
                                test.y,
                                test.p,
                                bProb = F,
                                bAUC = T)
  # interval estimate
  if (!is.null(post.samp)) {
    measure$pred.ci.in = PredictiveInterval(post.samp, x.obs, mu)
    measure$pred.ci = PredictiveInterval(post.samp, test.x, test.p)
    lb = apply(post.samp, 2, quantile, 0.025)
    ub = apply(post.samp, 2, quantile, 0.975)
    post.interval = cbind(lb, ub)
    measure$post.interval.loss = IntervalLoss(beta.true, lb, ub)
    containing = (beta.true - post.interval[, 1] > 0) * (beta.true - post.interval[, 2] < 0)
    measure$post.cover.rate  = containing
    measure$post.interval.len = post.interval[, 2] - post.interval[, 1]
  }
  
  # point estimate
  measure$mse = v.n.mse(betaHat, beta.true)
  return(measure)
}
MyROC = function(mu,
                 pred,
                 num = 20,
                 bPlot = F) {
  cutoffs = quantile(pred, probs = seq(0, num - 1) / num)
  tpfp = c()
  for (i in cutoffs) {
    posit = pred >= i
    tpr = sum(mu[posit]) / (sum(mu[posit]) + sum(mu[!posit]))
    fpr = sum(1 - mu[posit]) / (sum(1 - mu[posit]) + sum(1 - mu[!posit]))
    tpfp = rbind(tpfp, c(fpr, tpr))
  }
  if (bPlot)
    plot(tpfp)
  return(-sum ((c(tpfp[-1, 2], 0) + tpfp[, 2]) / 2 * diff(c(tpfp[, 1], 0))))
}
MyPredErrBinary = function(beta,
                           test.x,
                           test.y,
                           test.p,
                           bProb = F,
                           bAUC = F) {
  theta = cbind(1, test.x) %*% beta
  if (bAUC == T) {
    result = MyROC(test.p, theta)  # result=as.numeric(roc(test.p>0.5,pred=theta,plot=T,percent=T)$auc)
  } else{
    if (bProb) {
      pred.y = 1 / (1 + exp(-theta))
    } else{
      pred.y = theta > 0
    }
    result = mean(pred.y + test.p - 2 * test.p * pred.y)
  }
  return(result)
}

IntervalLoss = function(beta, lower, upper, nominal = 0.05) {
  mean(upper - lower + 1 / nominal * (abs(beta - upper) - upper + abs(lower -
                                                                        beta) + lower))
}

linearization = function(t, ss = 0.3, bb = 1 / ss * 4) {
  return(bb * (1 / (1 + exp(-ss * t)) - 0.5))
}

# sampling-----

## fn: generate the observed covariate
createX = function(cov.setting, n, p, par, bSecondOrder = F) {
  switch(
    cov.setting,
    uncor = {
      x = matrix(rnorm(n * p), nr = n)
    }, 
    acor = { 
      A = matrix(1:p ^ 2, p)
      A = abs(A - t(A)) / max(1, (p - 1))
      A = par ^ (A)
      x = matrix(rnorm(n * p), nr = n)
      x = x %*% chol(A)
    },  # AR(1) structure with rho_1=par
    intercor = {
      par = max(par, -1 / (p - 1))
      A = diag(p) * (1 - par) + par
      x = matrix(rnorm(n * p), nr = n)
      x = x %*% chol(A)
    },  # equal correlation (par)
    cor = {
      x = matrix(rnorm(n * p), nr = n)
      x = x %*% chol(par)
    },  # arbitrary covariance matrix (par)
    bern = {
      x = cbind(sapply(par, rbinom, n = n, size = 1))
    },  # Bernoulli with different probability (par)
    corbern = {
      eleidx = sample(2 ^ p, n, replace = T, par)
      basicx = factor.design(p)
      x = basicx[eleidx, ]
    }, # par contains the probability of each possible binary vector
    acorbern = {
      x = array(dim = c(n, p))
      x[, 1] = sample(c(1, -1), n, replace = T)
      if (p > 1) {
        for (i in 2:p) {
          x[, i] = rbinom(n, 1, prob = 1 / (1 + exp(-par * x[, i - 1]))) * 2 - 1
        }
      }
    },  # neighboring coordinates are linked by a logistic regression (coef.=par)
    acorbern_unbalance = {
      x = array(dim = c(n, p))
      x[, 1] = sample(c(1, -1), n, prob = c(2, 1), replace = T)  # a slight change
      if (p > 1) {
        for (i in 2:p) {
          x[, i] = rbinom(n, 1, prob = 1 / (1 + exp(-par * x[, i - 1]))) * 2 - 1
        }
      }
    }, 
    intercorbern = {
      par = max(par, -1 / (p - 1))
      A = diag(p) * (1 - par) + par
      x = matrix(rnorm(n * p), nr = n)
      x = x %*% chol(A)
      x = (x > 0) * 2 - 1
    },  # transform an intercorrelated normal vector into a binary vector
    intercormix = {
      par = max(par, -1 / (p - 1))
      A = diag(p) * (1 - par) + par
      x = matrix(rnorm(n * p), nr = n)
      x = x %*% chol(A)
      ind = 2 * (1:floor(p / 2))
      x[, ind] = (x[, ind] > 0) * 2 - 1
    },  # used in the paper; transform an intercorrelated normal vector into a mixed vector
    intercormix2 = {
      par = max(par, -1 / (p - 1))
      A = diag(p) * (1 - par) + par
      x = matrix(rnorm(n * p), nr = n)
      x = x %*% chol(A)
      ind = 1:floor(p / 2)   # a slight change 
      x[, ind] = (x[, ind] > 0) * 2 - 1
    },
    normal_factor = {
      A = matrix(1:p ^ 2, p)
      A = abs(A - t(A)) / max(1, (p - 1))
      A = par ^ (A)
      z = matrix(rnorm(n * p), nr = n)
      z = z %*% chol(A)
      
      mat.factors = data.frame(apply(
        z,
        2,
        cut,
        breaks = c(-Inf, qnorm(c(0.25, 0.5, 0.75)), Inf),  # the cut-offs
        labels = F
      ))
      for (j in 1:p) {
        mat.factors[, j] = as.factor(mat.factors[, j])
        contrasts(mat.factors[, j]) = contr.sum(4)
      }
      x = c()
      for (j in 1:p)
        x = cbind(x, model.matrix( ~ mat.factors[, j], contrasts.arg = contr.sum)[, -1])
    }  # transform an intercorrelated normal vector into a multi-factor vector
  ) 
  
  ## if need to include interactions
  if (bSecondOrder) {
    colnames(x) = paste("main", 1:p, sep = "")
    x.interact = t(apply(x, 1, combn, 2, prod))
    colnames(x.interact) <-
      paste("Inter.V", combn(1:p, 2, paste, collapse = "V"), sep = "")
    x = cbind(x, x.interact)
  }
  
  return(x)
}

## functions for calculating likelihood for logistic regression-----
logit = function(x) {
  1 / (1 + exp(-x))
}
## fn: log logistic regression density with weights
logden = function(beta, par)
{
  if (is.vector(beta) == T) {
    # for one vector
    beta[is.na(beta)] = 0 
    theta = par$x %*% beta
    id = theta > 700
    btheta = log(1 + exp(theta))
    btheta[id] = theta[id]
    return(sum(par$weights * (theta * par$y - btheta)))
  } else{
    # for multiple vectors
    theta = tcrossprod(par$x, beta)
    id = theta > 700
    btheta = log(1 + exp(theta))
    btheta[id] = theta[id]
    return(colSums(par$weights * (theta * par$y - btheta)))
  }
}
## fn: gradient of log density with weights
## CAUTION: this function is slow
logden.grad = function(beta, par)
{
  if (is.vector(beta) == T) {
    # one beta
    theta = par$x %*% beta
    eta = 1 / (1 + exp(-theta))
    return(c(c(par$weights * (par$y - eta)) %*% par$x))
    
  } else{
    # multiple beta
    theta =  tcrossprod(par$x, beta)# par$x %*% t(beta)
    eta = 1 / (1 + exp(-theta))
    return(crossprod(par$x, (par$weights * (par$y - eta))))
  }
}
# fn: second derivative of likelihood of logistic regression at a fitted model
SecondDerivative = function(model, inputx, intercept = T) {
  W = model$fitted.values * (1 - model$fitted.values) * model$prior.weights
  if (intercept == F)
  {
    M = t(inputx) %*% diag(W) %*% inputx
  } else
  {
    M = t(cbind(1, inputx)) %*% diag(W) %*% cbind(1, inputx)
  }
  return (M)
}
# approximate covariance matrix by inverting the information matrix
approxCov = function(model, inputx, intercept = T) {
  M = SecondDerivative(model, inputx, intercept)
  if (length(grep('Error', try(chol2inv(chol(M)), silent = T)
  )) == 1) {
    Res = chol2inv(chol(M + diag(ncol(M)) * 1e-10))
  } else
  {
    Res = chol2inv(chol(M))
  }
  return(Res)
}

## Operations----------
#create directory
MyCreateDir = function(path) {
  if (dir.exists(path))
    return
  string = path
  if (string[1] == '/')
    string = string[-1]
  i = 0
  l = regexpr('/', string)[1]
  while (l != -1) {
    #print(substr(path,0,i+l))#
    dir.create(substr(path, 0, i + l))
    string = substring(string, l + 1)
    i = i + l
    l = regexpr('/', string)[1]
  }
}
#string
MyExtract = function(element, location, begin) {
  as.numeric(sub(begin, '', element[location]))
}


##Calculation------
# fn: mean squared error
mse = function(a, b)
{
  if (is.vector(a))
    return(mean((a - b) ^ 2))
  return(rowMeans((a - b) ^ 2, na.rm = T))
}
# fn: normalized mse
n.mse = function(a, b){
  return(rowMeans(as.matrix((a - b) ^ 2 / (1 + b ^ 2)), na.rm = T))
}
# fn: normalized se
n.mse.se = function(a, b)
{
  return(sqrt(apply(a, 1, var, na.rm = T) / (1 + b ^ 2) / ncol(a)))
}
# fn: compute normalized mean square error
v.n.mse = function(a, b)
{
  return(mean(((a - b) ^ 2 / (1 + b ^ 2)), na.rm = T))
}
# fns: compute standard error for vector and matrix
sev=function(v){sqrt(var(v)/length(v))}
se=function(v){if(is.vector(v))return(sev(v))
  return(apply(v,2,sev))}
# fn: return all possible p-element set of 1:grid
factor.design = function(p, grid = 2) {
  if (p > 20)
    stop('p is too large')
  x = 1:grid[1]
  for (i in 2:p) {
    if (length(grid) > i) {
      y = c()
      for (j in 1:grid[i])
        y = rbind(y, cbind(x, j))
      x = y
    }
    else{
      y = c()
      #for(j in 1:grid[i])
      for (j in 1:2)
        y = rbind(y, cbind(x, j))
      x = y
    }
  }
  colnames(x) = paste('x', 1:p, sep = '')
  return(x)
}
# fn: calculate vector norm
vecnorm = function(x) {
  if (!is.matrix(x))
    return(sqrt(sum(x ^ 2)))
  return(sqrt(colSums(x ^ 2)))
}


### todo
SummarizeExp = function(Para.output,
                        idx.n,
                        MaxLoop,
                        bBinary = T,
                        except = NULL,
                        bSparse = F) {
  methodnames = names(Para.output[[1]][[1]])
  
  methodnames = setdiff(methodnames, except)
  m = length(methodnames)
  
  sum.pred.err = array(dim = c(m + 1, MaxLoop))
  sum.pred.err.in = array(dim = c(m + 1, MaxLoop))
  sum.pred.err01 = array(dim = c(m + 1, MaxLoop))
  sum.auc = array(dim = c(m + 1, MaxLoop))
  
  
  sum.pred.err.sparse = array(dim = c(m + 1, MaxLoop))
  sum.pred.err01.sparse = array(dim = c(m + 1, MaxLoop))
  
  
  sum.post.interval.loss = array(dim = c(m, MaxLoop))
  sum.pred.ci.in.coverage = array(dim = c(m, MaxLoop))
  sum.pred.ci.in.len = array(dim = c(m, MaxLoop))
  
  sum.pred.ci.coverage = array(dim = c(m, MaxLoop))
  sum.pred.ci.len = array(dim = c(m, MaxLoop))
  
  
  sum.count = rep(0, m)
  
  sum.post.cover.rate = array(dim = c(m, MaxLoop, d + 1))
  sum.betaHat = array(dim = c(m + 1, MaxLoop, d + 1))
  sum.post.interval.len = array(dim = c(m, MaxLoop, d + 1))
  sum.post.cover.rate = array(dim = c(m, MaxLoop, d + 1))
  
  sum.selectedtau = array(dim = c(MaxLoop, N.postSamp))
  hmc.accept = rep(NA, MaxLoop)
  
  for (j in 1:m) {
    for (i in 1:MaxLoop) {
      result = Para.output[[i]][[idx.n]][[methodnames[j]]]
      sum.count[j] <- sum.count[j] + result$count
      if (result$count == 1)
      {
        sum.betaHat[j, i, ] <- result$betaHat
        sum.post.cover.rate[j, i, ] = result$post.cover.rate
        sum.post.interval.len[j, i, ] = result$post.interval.len
        sum.post.interval.loss[j, i] = result$post.interval.loss [1]
        sum.pred.err[j, i] = result$pred.err
        sum.pred.err.in[j, i] = result$pred.err.in
        sum.pred.ci.in.coverage[j, i] = result$pred.ci.in$coverage
        sum.pred.ci.in.len[j, i] = result$pred.ci.in$len
        sum.pred.ci.coverage[j, i] = result$pred.ci$coverage
        sum.pred.ci.len[j, i] = result$pred.ci$len
        if (bBinary) {
          sum.pred.err01[j, i] = result$pred.err01
          sum.auc[j, i] = result$auc
        }
        if (bSparse) {
          sum.pred.err.sparse[j, i] = result$sparse$pred.err
          if (bBinary) {
            sum.pred.err01.sparse[j, i] = result$pred.err01
          }
        }
        if (methodnames[j] == 'Cata') {
          sum.pred.err[m + 1, i] = result$pred.err.trimmean
          sum.pred.err.in[m + 1, i] = result$pred.err.in.trimmean
          sum.betaHat[m + 1, i, ] = result$betaHat.trimmean
          sum.selectedtau[i, ] = result$selectedtau
          if (bBinary) {
            sum.pred.err01[m + 1, i] = result$pred.err01.trimmean
            sum.pred.err01[m + 1, i] = result$auc.trimmean
            hmc.accept[i] = sum(summary(result$gibbs.accepted[-1])[1:4] *
                                  10 ^ c(4, 2, 0, -2))
            
          }
          
        }
        
      }
    }
  }
  
  ######copy------
  names(sum.count) = methodnames
  
  row.names(sum.betaHat) = c(methodnames, 'Cata.trimmean')
  row.names(sum.pred.err) = c(methodnames, 'Cata.trimmean')
  row.names(sum.pred.err.in) = c(methodnames, 'Cata.trimmean')
  row.names(sum.pred.err01) = c(methodnames, 'Cata.trimmean')
  row.names(sum.auc) = c(methodnames, 'Cata.trimmean')
  
  
  if (bSparse) {
    row.names(sum.pred.err.sparse) = c(methodnames, 'Cata.trimmean')
    row.names(sum.pred.err01.sparse) = c(methodnames, 'Cata.trimmean')
  }
  
  
  row.names(sum.post.interval.loss) = methodnames
  row.names(sum.pred.ci.in.coverage) = methodnames
  row.names(sum.pred.ci.in.len) = methodnames
  row.names(sum.pred.ci.coverage) = methodnames
  row.names(sum.pred.ci.len) = methodnames
  
  row.names(sum.post.interval.len) = methodnames
  row.names(sum.post.cover.rate) = methodnames
  
  
  mse.data = apply(sum.betaHat, 1, function(mat) {
    apply(mat, 1, v.n.mse, beta.true)
  })
  result.mse = apply(sum.betaHat, 1, function(mat) {
    n.mse(t(mat), beta.true)
  })
  
  result.mse.se = apply(sum.betaHat, 1, n.mse.se0, beta.true)
  
  
  result.post.cover.rate = apply(sum.post.cover.rate, 1, function(mat)
    colMeans(mat))
  result.post.interval.len = apply(sum.post.interval.len, 1, function(mat)
    colMeans(mat))
  
  
  
  #####paste------
  
  list(
    result.mse = result.mse,
    result.mse.se = result.mse.se,
    overallmse = colMeans(result.mse),
    result.pred.err.mean = apply(sum.pred.err, 1, mean),
    result.pred.err.se = apply(sum.pred.err, 1, function(v) {
      sd(v) / sqrt(length(v))
    }),
    
    result.pred.err.in.mean = apply(sum.pred.err.in, 1, mean),
    result.pred.err.in.se = apply(sum.pred.err.in, 1, function(v) {
      sd(v) / sqrt(length(v))
    }),
    result.pred.err.mean.sparse = apply(sum.pred.err.sparse, 1, mean),
    
    result.pred.err01.mean = apply(sum.pred.err01, 1, mean),
    result.auc = apply(sum.auc, 1, mean),
    result.pred.err01.mean.sparse = apply(sum.pred.err01.sparse, 1, mean),
    
    result.post.cover.rate = result.post.cover.rate,
    result.post.interval.len = result.post.interval.len,
    
    
    hmc.accept = hmc.accept,
    selectedtau = sum.selectedtau,
    exist = sum.count,
    
    estimators = sum.betaHat,
    
    msedata = mse.data,
    pedata = t(sum.pred.err),
    peindata = t(sum.pred.err.in),
    pe01data = t(sum.pred.err01),
    aucdata = t(sum.auc),
    interval.len = sum.post.interval.len,
    interval.cover = sum.post.cover.rate,
    interval.loss = sum.post.interval.loss,
    pred.ci.in.coverage = sum.pred.ci.in.coverage,
    pred.ci.in.len = sum.pred.ci.in.len,
    pred.ci.coverage = sum.pred.ci.coverage,
    pred.ci.len = sum.pred.ci.len
    
  )
}
