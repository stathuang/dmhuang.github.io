##############################
####### Illustration of Catalytic priors for logistic regression
####### Dongming Huang (statdongminghuang@gmail.com)
####### date: 2020-07-15
####### Running the example in the paper of catalytic prior
##############################


##############################
## Meanings of the outputs
## pred.err : expected predictive binomial deviance
## ## pred.err.in : expected in-sample binomial deviance
## pred.err01: expected predictive classification error
## auc : AUC
## post.cover.rate: averaged coverage probability of 95% posterior interval
## post.interval.len: averaged length of 95% posterior interval
##############################

### load library
library(arm)
library(rstan) # can comment out if stan is not used
source("simpleSampling.R")
source("simpleHandy.R")

### set up environment  ----

output.path = 'output/'    # path for storing output
MyCreateDir(output.path)

MaxEpo = 2      # number of repeatation done on a single core; should change dependently
MaxParaLoop = 20  # number of cores used for a specific problem
ArrayID = 1       # job id; should be specified by computation platform

ParaID = 1 + ((ArrayID - 1) %/% MaxParaLoop)   # id of a problem, from 1 to 18
loop  = 1 + ((ArrayID - 1) %% MaxParaLoop)     # pointer of a job for this problem, from 1 to 20
print(c(ParaID, loop))



### experiment setup 

### Distributional parameters-----------

## Varying (pattern * 6, r * 3)
Pattern = c(
  'mixture_norm=0.25',
  'mixture_norm=0.5',
  'mixture_norm=0.75',
  'mixture_mass=0.25',
  'mixture_mass=0.5',
  'mixture_mass=0.75'
)
Grid.overlap = 0.1 * (1:3)   # factor r; oracle classification error

## Held fixed
bMulti = F                  # only main effects
Grid.covariate = c(
  'acorbern',
  'acor',
  'factors',
  'intercor',
  'intercorbern',
  'intercormix',
  'intercormix2'
)
cov.setting = Grid.covariate[7]

Grid.p = 16               # dimension
Grid.n = c(30, 200)        # sample sizes
Grid.coe = 0.5              # parameter for dependence among x's

TrueLink = function(fhi,
                    a = 100,
                    b = 1,
                    c = 0.2) {
  return(1 / (1 + exp(-fhi)))
}

TrueEta = function(x, beta, add.intercept = T) {
  if (add.intercept) {
    beta[1] + x %*% beta[-1]
  } else{
    x %*% beta
  }
}


N.test = 1000               # number of test samples

DevianAdjustments = c("none", "kl", "over-entropy")
criterion = DevianAdjustments[2]     # use KL diverence as error measurement

### Procedural parameters-----

Grid.M = c(400)
Grid.tau = seq(0, 60, 0.05)[-1]
pseudo.cov.dist = 'default'  # use independent resampling to generate pseudo x
pseudo.model = 'intercept'  # model used to generate pseudo y

## consider decreasing these parameter to speed up 
N.PostSamp = 1000         # number of posterior samples
BurnIn = 4000             # iterations of initalization for HMC
N.boot = 100              # number of Bootstrap replications

### Setting based on job id----

idx.pattern = 1 + ((ParaID - 1) %/% length(Grid.overlap))
idx.s =      1 + ((ParaID - 1) %% length(Grid.overlap))

pattern = Pattern[idx.pattern]
set.seed(ParaID)

p = Grid.p[1]
if (bMulti == T) {
  d = (p + 1) * p / 2
} else {
  d = p
}
print(p)
coe = Grid.coe[1]
M = Grid.M[1]   # for this example, M is fixed at 400

### main loop---------

Result = list()

for (epo in 1:MaxEpo)
{
  print(epo)
  seed.id = (loop - 1) * MaxEpo + epo
  set.seed(seed.id)
  
  # set up the true regression coefficients 
  beta.true0 = AssignRegCoef(strsplit(pattern, "=")[[1]][1],
                             p,
                             seed = 1,
                             frac = as.numeric(strsplit(pattern, "=")[[1]][2]))
  overlap = Grid.overlap[idx.s]
  x.ref = createX(cov.setting, max(Grid.n, 1000), p, coe, bSecondOrder = bMulti)
  scale.beta = 1
  # here, we keep resacling the true beta
  # until the oracle classification is close to overlap (r)
  repeat {
    beta.true = beta.true0 * scale.beta
    beta.true[1] = 0
    fhi.ref =  TrueEta(x.ref, beta.true, add.intercept = T)
    beta.true[1] = -mean(fhi.ref)
    fhi.ref =  TrueEta(x.ref, beta.true, add.intercept = T)
    prob.ref =  TrueLink(fhi.ref)
    oracle.risk = mean(1 - abs(1 - 2 * prob.ref)) / 2
    if (abs(oracle.risk - overlap) < 0.005) 
      break
    scale.beta = ifelse(oracle.risk > overlap, scale.beta * 1.25, scale.beta * 0.9) 
  }
  
  
  # generate real sample
  N = max(Grid.n)
  sim.x0 = createX(cov.setting, N, p, coe, bSecondOrder = bMulti)
  test.x = createX(cov.setting, N.test, p, coe, bSecondOrder = bMulti)
  
  fhi0 = TrueEta(sim.x0, beta.true, add.intercept = T)
  mu0 = TrueLink(fhi0)
  sim.y0 = rbinom(N, 1, mu0)
  
  # generate test sample
  test.fhi = TrueEta(test.x, beta.true, add.intercept = T)
  test.p = TrueLink(test.fhi)
  test.y = rbinom(length(test.p), 1, test.p)
  
  
  para.output = list()
  for (idx.n in 1:length(Grid.n)) {
    sim.n = Grid.n[idx.n]
    sim.x = sim.x0[1:sim.n, ]
    sim.y = sim.y0[1:sim.n]
    fhi = fhi0[1:sim.n]
    mu = mu0[1:sim.n]
    
    
    # MLE
    output = list()
    output$mle = RunMLELogistic(
      sim.x,
      sim.y,
      N.postSamp = N.PostSamp,
      bMeasure = T,
      mu = mu,
      fhi = fhi,
      test.x = test.x,
      test.y = test.y,
      test.p = test.p,
      criterion = criterion
    )
    output$mle$mse = v.n.mse(output$mle$betaHat, beta.true)
    
    # check if well-separated
    logistic = glm(sim.y ~ sim.x, family = binomial)
    output$completeseparated = sum(sim.y == as.integer(logistic$fitted.values > 0.5))
    
    # Cauchy prior
    
    set.seed(seed.id + 1)
    output$cauchy = RunCauchy(
      x.obs = sim.x,
      y.obs = sim.y,
      N.postSamp = N.PostSamp,
      bSamp = F,
      bMeasure = T,
      mu = mu,
      fhi = fhi,
      test.x = test.x,
      test.y = test.y,
      test.p = test.p,
      criterion = criterion
    )
    output$cauchy$mse = v.n.mse(output$cauchy$betaHat, beta.true)
    
    ### Catalytic priors -----
    set.seed(seed.id + 5)
    
    Cata = list()
    
    # generate pseudo sample
    data.cata = InitCatalytic(
      M = 400,
      x = sim.x,
      y = sim.y,
      p = p,
      pseudo.cov.dist = pseudo.cov.dist,
      pseudo.model = pseudo.model
    )
    
    
    
    ## select tau by Info-Criterion
    Grid.taus = (d + 1) * c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 2)
    
    ic.stein = EstimatePredErr.Stein(Grid.taus, data.cata, var.tau = (d + 1) / 4)
    ic.boot = EstimatePredErr.boot(Grid.taus,
                                   data.cata,
                                   var.tau = (d + 1) / 4,
                                   B = N.boot)
    
    
    tau.steinian = Grid.taus[which.min(ic.stein)]   # steinian
    tau.boot = Grid.taus[which.min(ic.boot)]  # bootstrap
    
    output$Cata.Steinian.quar = CatalyticLogistic(
      data.cata,
      tau.steinian,
      N.postSamp = N.PostSamp,
      bStan = F,
      bMeasure = T,
      fhi = fhi,
      test.x = test.x,
      test.y = test.y,
      test.p = test.p,
      criterion = criterion
    )
    output$Cata.boot = CatalyticLogistic(
      data.cata,
      tau.boot,
      N.postSamp = N.PostSamp,
      bStan = F,
      bMeasure = T,
      fhi = fhi,
      test.x = test.x,
      test.y = test.y,
      test.p = test.p,
      criterion = criterion
    )
    
    
    
    ## full bayesian
    output$fullbayesian = CatalyticLogisticFullBayesian(
      data.cata,
      N.postSamp = N.PostSamp,
      burnin = BurnIn,
      Grid.tau = Grid.tau,
      alpha = 2,
      gammainv = 1,
      bOuputSamp = F,
      bTrun = T,   # we remove the small values of tau here; useful for small n 
      bVerb = F,
      bMeasure = T,
      mu = mu,
      fhi = fhi,
      test.x = test.x,
      test.y = test.y,
      test.p = test.p,
      criterion = criterion
    )
    
    para.output[[idx.n]] = output
    cat('>')   
  }
  
  Result[[epo]] = para.output
  
  ## store
  save(Result,
       file = paste(output.path, ArrayID, '.RData', sep = ''))
}