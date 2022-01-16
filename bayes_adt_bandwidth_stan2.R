setwd("C:/od/Bayesian precision medicine trial")
library(MASS)
library(matlib)
library(doParallel)
library(rstan)
library(bayesplot)
library(mvtnorm)


### data generator ####
data_lg = function(p=2, a1=1, a2=1, n_sample=30, beta=seq(1,p+1), fam= "gaussian" ) {
  # a1 is the variance for the covariates
  # a2 is the variance for the error term
  
  covariate_mean = rep (100, p)
  covariate_cov = diag(a1, nrow = p, ncol = p)
  x =  MASS::mvrnorm (n_sample, mu= covariate_mean, Sigma = covariate_cov) # True Covariates matrix
  Z = cbind(rep(1, n_sample), x)     # True Design matrix
  if (fam == "gaussian") {
    Y = MASS::mvrnorm(1, mu= Z %*% beta , diag(a2, nrow = n_sample, ncol = n_sample) ) # outcome vector 
    #yd= mvtnorm::dmvnorm(x, covariate_mean, covariate_cov, log = FALSE)
    
  }
  # if (fam=="poisson"){
  #   Y = rpois(n_sample, lambda = exp (Z %*% beta))
  # }
  
  if (fam == "sin"){
    Y  = MASS::mvrnorm(1, mu= 2* sin ((Z %*% beta)*pi) , diag(a2, nrow = n_sample, ncol = n_sample) ) 
  }
  
  ## binary predictors
  
  # x=rbinom(n_sample,1, p=0.7)
  # x<-as.matrix(x)
  # Z = cbind(rep(1, n_sample), x)     # True Design matrix
  # Y = MASS::mvrnorm(1, mu= Z %*% beta , diag(a, nrow = n_sample, ncol = n_sample) ) # outcome vector 
  gdata<-list(n_sample = n_sample, #
              p=p,
              X=x,
              l0=1,
              Y=Y,
              Yf=Y
              # ,yd=yd
  )
  return(gdata)
}

set.seed(123)
p=2
d1.l<- data_lg(p=2, a1=1, a2=1, n_sample=100, beta=seq(1,p+1), fam= "gaussian" )
d2.sin<-data_lg(p=2, a1=1, a2=1, n_sample=100, beta=seq(1,p+1), fam= "sin" )

## data distribution
par(mfrow=c(2,2))
hist(d1.l$X[,1])
hist(d1.l$X[,2])
plot(d1.l$X[,1],d1.l$X[,2])
hist(d1.l$Y)

par(mfrow=c(2,2))
hist(d2.sin$X[,1])
hist(d2.sin$X[,2])
plot(d2.sin$X[,1],d2.sin$X[,2])
hist(d2.sin$Y)


##### Create an initial smoothing window for covariates #####
d1.l$alpha = 1
d1.l$beta = 1

## compute the posterior mean empirically by Stan
init.b=function(x, dat = d1.l, alpha0, beta0){
  # x is a vector of (dimension of p) the covariate pattern to be predicted
  # alpha is the shape parameter of inverse gamma dist
  # beta is the scale parameter of inverse gamma dist
  d1.l$alpha = alpha0
  d1.l$beta = beta0
  p = dat$p
  dat$x = x
  band.postm = rep(0, p)
  md<-stan_model(file = 'intial.window.stan')
  fit<-sampling(md, data=dat, cores = 1, warmup = 2000,iter = 2000+80000/4 ,
                  chains = 1,control = list(adapt_delta = 0.99)) # 1
  
  for (i in 1:p){
    band.postm [i] = get_posterior_mean(fit)[i]
  }
  return(band.postm)
}
# test 
init.b(x = c(100, 100), alpha0 = 0.001, beta0 = 0.001)

#d1.l.2 <- data_lg(p=2, a1=1, a2=1, n_sample=1000, beta=seq(1,p+1), fam= "gaussian" )
### Use the closed form to derive posterior mean ### 
postm.invg<-function(x, dat = d1.l, a0, b0){
  dat$alpha = a0
  dat$beta = b0
  p = dat$p
  X = dat$X
  n_sample = dat$n_sample
  nu.u = matrix(1, nrow = n_sample, ncol = p)  # the unit in the numerator
  de.u = matrix(1, nrow = n_sample, ncol = p)  # the unit in the denominator
  
  band.postm = rep(0, p)
  for (i in 1:p){
    for (j in 1:n_sample){
      nu.u[j,i] = (1 / (b0 * (X[j,i] - x[i])^2 +2 )    ) ^ a0
      de.u[j,i] = (1 / (b0 * (X[j,i] - x[i])^2 +2 )    ) ^ (a0 + 1/2)
    }
    band.postm [i] = (gamma(a0)/(sqrt(2*b0)*gamma(a0+0.5))) *
                    (sum(nu.u[,i])/ sum(de.u[,i]))
  }
  return(band.postm)
}

#test
postm.invg(x = c(100, 100), dat = d1.l, a0 = 0.5, b0 = 100^(0.2))  
## the range
r1<-c(100-4, 100+4)
r2<-c(100-4, 100+4)
n1<-which (d1.l$X[,1]>=r1[1] & d1.l$X[,1]<=r1[2])
n2<-which(d1.l$X[,2]>=r2[1] & d1.l$X[,2]<=r2[2])
is<-intersect(n1,n2)
d1.l.selected<-d1.l
d1.l.selected$X<-d1.l.selected$X[is,]
d1.l.selected$Y<-d1.l.selected$Y[is]
d1.l.selected$n_sample<-length(is)
## calling stan

# mg is bandwidth
mg1.1<-stan_model(file = 'global_bandwidth.stan')
#mg1.2<-stan_model(file = 'global_bandwidth2.stan') # Changed the log in the denominator as well, force no negative values in the simualted data

## fitg2 is the global bandwidth using all data 
fitg2<-sampling(mg1.1, data=d1.l, cores = 1, warmup = 2000,iter = 2000+10000/4 ,
                chains = 1,control = list(adapt_delta = 0.99)) # 1

## fitg3 is the local bandwidth using selected data
fitg3<-sampling(mg1.1, data=d1.l.selected, cores = 1, warmup = 2000,iter = 2000+20000/4 ,
               chains = 1,control = list(adapt_delta = 0.99)) # 0.1


#fit3<-vb(mg,iter = 100000, output_samples = 10000)
summary(fitg2)
save(fitg2, file = "fitg2.RData")
pos2<-as.array(fitg2)
mcmc_areas(pos2, pars= c("bandwidth[1]","bandwidth[2]" ) )
mcmc_hist(pos2, pars= c("bandwidth[1]","bandwidth[2]" ) )
mcmc_areas(pos2, pars=  "ysd" )
mcmc_trace(pos2, pars=  "ysd" )
mcmc_trace(pos2, pars= c("bandwidth[1]","bandwidth[2]" ) )


#f2<-rstan::vb(mg2.2, data=dat,tol_rel_obj= 0.00001, iter = 300000)

summary(fitg3)
pos3<-as.array(fitg3)
mcmc_areas(pos3, pars= c("bandwidth[1]","bandwidth[2]" ) )
mcmc_hist(pos3, pars= c("bandwidth[1]","bandwidth[2]" ) )
mcmc_areas(pos3, pars=  "ysd" )
mcmc_trace(pos3, pars=  "ysd" )
mcmc_trace(pos3, pars= c("bandwidth[1]","bandwidth[2]" ) )












