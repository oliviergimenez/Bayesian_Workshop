## ----setup, include = FALSE----------------------------------------------------------------------------------------
knitr::opts_chunk$set(cache = TRUE, 
                      echo = TRUE, 
                      message = FALSE, 
                      warning = FALSE,
                      fig.height=6, 
                      fig.width = 1.777777*6,
                      tidy = FALSE, 
                      comment = NA, 
                      highlight = TRUE, 
                      prompt = FALSE, 
                      crop = TRUE,
                      comment = "#>",
                      collapse = TRUE)
knitr::opts_knit$set(width = 60)
library(tidyverse)
library(reshape2)
theme_set(theme_light(base_size = 16))
make_latex_decorator <- function(output, otherwise) {
  function() {
      if (knitr:::is_latex_output()) output else otherwise
  }
}
insert_pause <- make_latex_decorator(". . .", "\n")
insert_slide_break <- make_latex_decorator("----", "\n")
insert_inc_bullet <- make_latex_decorator("> *", "*")
insert_html_math <- make_latex_decorator("", "$$")


## ---- out.width = '13cm',out.height='7cm',fig.align='center',echo=FALSE--------------------------------------------
knitr::include_graphics('img/brace_yourself.jpeg')


## ---- out.width = '13cm',out.height='7cm',fig.align='center',echo=FALSE--------------------------------------------
knitr::include_graphics('img/books.jpeg')


## ---- fig.align = 'center', echo = FALSE---------------------------------------------------------------------------
knitr::include_graphics('img/mccarthy.png')


## ---- fig.align = 'center', echo = FALSE---------------------------------------------------------------------------
knitr::include_graphics('img/kery.png')


## ---- fig.align = 'center', echo = FALSE---------------------------------------------------------------------------
knitr::include_graphics('img/kruschke.png')


## ---- fig.align = 'center', echo = FALSE---------------------------------------------------------------------------
knitr::include_graphics('img/mcelreath.png')


## ---- fig.align = 'center', echo = FALSE---------------------------------------------------------------------------
knitr::include_graphics('img/gelmanhill.png')


## ---- echo=FALSE, fig.align='center'-------------------------------------------------------------------------------
knitr::include_graphics('img/frequentists_vs_bayesians_2x.png')


## ----collapse=TRUE-------------------------------------------------------------------------------------------------
dbinom(x=3,size=10,prob=0.1)


## ----collapse=TRUE-------------------------------------------------------------------------------------------------
dbinom(x=3,size=10,prob=0.9)


## ----collapse=TRUE-------------------------------------------------------------------------------------------------
dbinom(x=3,size=10,prob=0.25)


## ----collapse=TRUE-------------------------------------------------------------------------------------------------
dbinom(x=3,size=10,prob=0.1)
dbinom(x=3,size=10,prob=0.9)
dbinom(x=3,size=10,prob=0.25)


## ----collapse=TRUE-------------------------------------------------------------------------------------------------
lik.fun <- function(parameter){
  ll <- dbinom(x=3, size=10, prob=parameter)
  return(ll)
}

lik.fun(0.3)

lik.fun(0.6)


## ---- echo=FALSE---------------------------------------------------------------------------------------------------
lik.fun <- function(parameter){
  ll <- dbinom(x=3, size=10, prob=parameter)
  return(ll)
}
p.grid = seq(0,1,by=0.01)
lik = rep(NA,length(p.grid))
for (i in 1:length(p.grid)){
  lik[i] <- lik.fun(p.grid[i])
}
plot(p.grid,lik,xlab='Probability of getting a successful breeding',ylab='Likelihood',type='l',lwd=3,cex.lab=1.5)
abline(v=0.3,lty=2,lwd=2,col='blue')


## ----collapse=TRUE-------------------------------------------------------------------------------------------------
lik.fun <- function(parameter) dbinom(x=3, size=10, prob=parameter)
# ?optimize
optimize(lik.fun,c(0,1),maximum=TRUE)


## ---- echo=FALSE---------------------------------------------------------------------------------------------------
lik.fun <- function(parameter) dbinom(x=3, size=10, prob=parameter)
plot(lik.fun,0,1,xlab="probability of success (p)",ylab="log-likelihood(p)",main="Binomial likelihood with 3 successes ot of 10 attempts",lwd=3,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
abline(v=0.3,h=0.26682,col='blue',lty=2,lwd=2)


## ------------------------------------------------------------------------------------------------------------------
# set seed for random numbers
set.seed(2020)
# simulate data from Normal distribution
n <- 100
height <- rnorm(n, mean=170, sd=10)


## ----echo=FALSE----------------------------------------------------------------------------------------------------
data.frame(height=height) %>%
  ggplot(aes(x=height))+ 
  geom_histogram(color="blue",fill="dodgerblue") + 
  labs(x = "Height", y = 'Density')


## ------------------------------------------------------------------------------------------------------------------
negloglik <- function(theta, data) {
  mu <- theta[1]
  sigma <- theta[2]
  x <- data
  -sum(dnorm(x, mean = mu, sd = sigma, log = TRUE))
}
negloglik(theta = c(150,1), height)


## ------------------------------------------------------------------------------------------------------------------
fit <- optim(par = c(1,1), fn = negloglik, data = height)
fit


## ---- echo = FALSE-------------------------------------------------------------------------------------------------
binwidth <- 1 
df <- data.frame(x = height) %>%
  ggplot(aes(x = height, mean = fit$par[1], sd = fit$par[2], binwidth = binwidth, n = n)) +
  geom_histogram(binwidth = binwidth, 
                 colour = "white", 
                 fill = "cornflowerblue", 
                 size = 0.1) +
  stat_function(fun = function(x) dnorm(x, mean = fit$par[1], sd = fit$par[2]) * n * binwidth,
    color = "darkred", size = 1)
df


## ---- echo=FALSE---------------------------------------------------------------------------------------------------
x <- seq(0, 1, length=200)
par(mfrow = c(2,3))
# distribution a posteriori beta
plot(x,dbeta(x, 1, 1),type='l',xlab='q',ylab='Density',main='beta(1,1)',lwd=3,col='red',ylim=c(0,1.5))
plot(x,dbeta(x, 2, 1),type='l',xlab='q',ylab='',main='beta(2,1)',lwd=3,col='red',ylim=c(0,2))
plot(x,dbeta(x, 1, 2),type='l',xlab='q',ylab='',main='beta(1,2)',lwd=3,col='red',ylim=c(0,2))
plot(x,dbeta(x, 2, 2),type='l',xlab='q',ylab='Density',main='beta(2,2)',lwd=3,col='red',ylim=c(0,1.5))
plot(x,dbeta(x, 10, 10),type='l',xlab='q',ylab='',main='beta(10,10)',lwd=3,col='red',ylim=c(0,3.5))
plot(x,dbeta(x, 0.8, 0.8),type='l',xlab='q',ylab='',main='beta(0.8,0.8)',lwd=3,col='red',ylim=c(0.5,2.5))


## ----echo=FALSE----------------------------------------------------------------------------------------------------
x <- seq(0, 1, length=200)
# distribution a posteriori beta
plot(x,dbeta(x, 20,39),type='l',xlab='',ylab='',main='',lwd=3,col='red')
# distribution a priori uniforme
points(x,dbeta(x, 1, 1),type='l',lwd=3)


## ----echo=FALSE----------------------------------------------------------------------------------------------------
x <- seq(0, 1, length=200)
# distribution a posteriori beta
plot(x,dbeta(x, 20,39),type='l',xlab='',ylab='',main='',lwd=3,col='red')
# distribution a priori uniforme
points(x,dbeta(x, 1, 1),type='l',lwd=3)
abline(v = 19/57, lwd = 3, lty = 2, col = 'blue')
text(x = 0.28, y = 0, 'MLE', col = 'blue')


## ----echo=FALSE----------------------------------------------------------------------------------------------------
x <- seq(0, 1, length=200)
# distribution a posteriori beta
plot(x,dbeta(x, .5+19,.5+57-19),type='l',xlab='',ylab='',main='',lwd=3,col='red')
# distribution a priori uniforme
points(x,dbeta(x, .5, .5),type='l',lwd=3)


## ----echo=FALSE----------------------------------------------------------------------------------------------------
x <- seq(0, 1, length=200)
# distribution a posteriori beta
plot(x,dbeta(x, 2+19,2+57-19),type='l',xlab='',ylab='',main='',lwd=3,col='red')
# distribution a priori uniforme
points(x,dbeta(x, 2, 2),type='l',lwd=3)


## ----echo=FALSE----------------------------------------------------------------------------------------------------
x <- seq(0, 1, length=200)
# distribution a posteriori beta
plot(x,dbeta(x, 20+19,1+57-19),type='l',xlab='',ylab='',main='',lwd=3,col='red')
# distribution a priori uniforme
points(x,dbeta(x, 20, 1),type='l',lwd=3)


## ---- out.width = '11cm',out.height='7cm',fig.align='center',echo=FALSE--------------------------------------------
knitr::include_graphics('img/falling_man.jpg')    


## ---- out.width = '10cm',out.height='6cm',fig.align='center',echo=FALSE--------------------------------------------
knitr::include_graphics('img/dipper.png')    


## ----include=FALSE-------------------------------------------------------------------------------------------------
# read in data
data <- as.matrix(read.table("dat/dipper.dat"))

# number of individuals 
n <- dim(data)[[1]] 

# number of capture occasions
K <- dim(data)[[2]] 

# compute the date of first capture
e <- NULL
for (i in 1:n){
	temp <- 1:K
	e <- c(e,min(temp[data[i,]==1]))
	}

# data
datax <- list(N=n,Years=K,obs=data,First=e)

# mark-recapture analysis for European Dippers
model <- 
paste("
model
{
for (i in 1:N){
	alive[i,First[i]] <- 1
	for (j in (First[i]+1):Years){
		alive[i,j] ~ dbern(alivep[i,j])
		alivep[i,j] <- surv * alive[i,j-1]
		obs[i,j] ~ dbern(sightp[i,j])
		sightp[i,j] <- resight * alive[i,j]
		}
	}
surv~dunif(0,1)
resight~dunif(0,1)
}
")
writeLines(model,"code/CJS.txt")

# In JAGS we have to give good initial values for the latent state alive. At all occasions when an individual was observed, its state is alive = 1 for sure. In addition, if an individual was not observed at an occasion, but was alive for sure, because it was observed before and thereafter (i.e. has a capture history of e.g. {101} or {10001}), then we know that the individual was alive at all of these occasions, and thus alive = 1. Therefore, we should provide initial values of alive = 1 at these positions as well. The following function provides such initial values from the observed capture histories (from Kery and Schaub book)

known.state.cjs <- function(ch){
   state <- ch
   for (i in 1:dim(ch)[1]){
      n1 <- min(which(ch[i,]==1))
      n2 <- max(which(ch[i,]==1))
      state[i,n1:n2] <- 1
      state[i,n1] <- NA
      }
   state[state==0] <- NA
   return(state)
   }

Xinit <- known.state.cjs(data)

# first list of inits
init1 <- list(surv=.1,resight=.1,alive=Xinit)
# second list of inits
init2 <- list(surv=.9,resight=.9,alive=Xinit)

# specify the parameters to be monitored
parameters <- c("resight","surv")

# load R2jags
library(R2jags)

# run the MCMC analysis WITHOUT PRIOR INFORMATION
CJS.sim <-jags(data=datax, inits=list(init1,init2), parameters,n.iter=1000,model.file="code/CJS.txt",n.chains=2,n.burnin=500)

# to see the numerical results
# CJS.sim
# traceplot(CJS.sim) # diagnostic de convergence

# keep 3 first years only
data = data[,1:3]
databis = NULL
for (i in 1:nrow(data)){
	# discard all non existing individuals i.e. those that were never captured
	# test whether there was at least 1 detection and keep this individual if it was the case
	if (sum(data[i,] == c(0,0,0))<3)  databis = rbind(databis,data[i,])
	}
data = databis

# number of individuals 
n <- dim(data)[[1]] 

# number of capture occasions
K <- dim(data)[[2]] 

# compute the date of first capture
e <- NULL
for (i in 1:n){
	temp <- 1:K
	e <- c(e,min(temp[data[i,]==1]))
	}

# data
datax <- list(N=n,Years=K,obs=data,First=e)

Xinit <- known.state.cjs(data)

# first list of inits
init1 <- list(surv=.1,resight=.1,alive=Xinit)
# second list of inits
init2 <- list(surv=.9,resight=.9,alive=Xinit)

# specify the parameters to be monitored
parameters <- c("resight","surv")

# run the MCMC analysis WITHOUT PRIOR INFORMATION
CJS.sim.wo.apriori <-jags(data=datax, inits=list(init1,init2), parameters,n.iter=1000,model.file="code/CJS.txt",n.chains=2,n.burnin=500)

# same model but with informative prior on survival 
model <- 
paste("
model
{
for (i in 1:N){
	alive[i,First[i]] <- 1
	for (j in (First[i]+1):Years){
		alive[i,j] ~ dbern(alivep[i,j])
		alivep[i,j] <- surv * alive[i,j-1]
		obs[i,j] ~ dbern(sightp[i,j])
		sightp[i,j] <- resight * alive[i,j]
		}
	}
surv~dnorm(0.57,187.6) # Norm(0.57,sd=0.073) ; precision = 1/var = 1/0.073^2
resight~dunif(0,1)
}
")
writeLines(model,"code/CJS2.txt")

CJS.sim.apriori <-jags(data=datax, inits=list(init1,init2), parameters,n.iter=1000,model.file="code/CJS2.txt",n.chains=2,n.burnin=500)


## ----echo=FALSE, message=FALSE, warning=FALSE----------------------------------------------------------------------
res = as.mcmc(CJS.sim.wo.apriori) 
res = rbind(res[[1]],res[[2]]) 
#head(res)

res2 = as.mcmc(CJS.sim.apriori) 
res2 = rbind(res2[[1]],res2[[2]]) 
#head(res2)

plot(density(res2[,'surv']),xlab='survival',ylab='probability density',col='red',lwd=4,main='',xlim=c(0.2,1))
lines(density(res[,'surv']),xlab='survival',ylab='probability density',col='blue',lwd=4,main='')
legend('topleft',lwd=2,legend=c('with prior info','without prior info'),col=c('red','blue'))


## ----echo = TRUE---------------------------------------------------------------------------------------------------
(alpha <- ( (1 - 0.57)/(0.073*0.073) - (1/0.57) )*0.57^2)
(beta <- alpha * ( (1/0.57) - 1))


## ----eval = FALSE, echo = TRUE-------------------------------------------------------------------------------------
## alpha <- ( (1 - 0.57)/(0.073*0.073) - (1/0.57) )*0.57^2
## beta <- alpha * ( (1/0.57) - 1)
## n <- 10000
## samp <- rbeta(n, alpha, beta)
## (mu <- mean(samp))
## (sigma <- sqrt(var(samp)))


## ----echo=1, fig.height=3, fig.width=3, fig.align='left'-----------------------------------------------------------
plot(density(rnorm(1000, 0, 1000)),   
     main="", xlab="Height (m)")


## ----echo=1, fig.height=3, fig.width=3, fig.align='left'-----------------------------------------------------------
plot(density(rnorm(1000, 2, 0.5)),   
      main="", xlab="Height (m)")


## ----echo=1, fig.height=3, fig.width=3, fig.align='left'-----------------------------------------------------------
plot(density(plogis(rnorm(1000,0,10)), 
from = 0, to = 1), main='', xlab='survival')


## ----echo=1, fig.height=3, fig.width=3, fig.align='left'-----------------------------------------------------------
plot(density(plogis(rnorm(1000,0,1.5)), 
from = 0, to = 1), main='', xlab='survival')


## ------------------------------------------------------------------------------------------------------------------
y <- 19 # nb of success
n <- 57 # nb of attempts


## ------------------------------------------------------------------------------------------------------------------
a <- 1; b <- 1; p <- seq(0,1,.002)
plot(p, dbeta(p,a,b), type='l', lwd=3)


## ------------------------------------------------------------------------------------------------------------------
numerator <- function(p) dbinom(y,n,p)*dbeta(p,a,b)


## ------------------------------------------------------------------------------------------------------------------
denominator <- integrate(numerator,0,1)$value


## ------------------------------------------------------------------------------------------------------------------
plot(p, numerator(p)/denominator,type="l", lwd=3, col="green", lty=2)


## ----eval = FALSE--------------------------------------------------------------------------------------------------
## lines(p, dbeta(p,y+a,n-y+b), col='darkred', lwd=3)


## ----echo = FALSE--------------------------------------------------------------------------------------------------
plot(p, numerator(p)/denominator,type="l", lwd=3, col="green", lty=2)
lines(p, dbeta(p,y+a,n-y+b), col='darkred', lwd=3)


## ----eval = FALSE--------------------------------------------------------------------------------------------------
## lines(p, dbeta(p,a,b), col='darkblue', lwd=3)


## ----echo = FALSE--------------------------------------------------------------------------------------------------
plot(p, numerator(p)/denominator,type="l", lwd=3, col="green", lty=2)
lines(p, dbeta(p,y+a,n-y+b), col='darkred', lwd=3)
lines(p, dbeta(p,a,b), col='darkblue', lwd=3)


## ---- out.width = '9cm',out.height='3cm',fig.align='center',echo=FALSE---------------------------------------------
knitr::include_graphics('img/metropolis.png')   


## ---- out.width = '11cm',out.height='7cm',fig.align='center',echo=FALSE--------------------------------------------
knitr::include_graphics('img/maniac.png')   


## ------------------------------------------------------------------------------------------------------------------
pd <- function(x){
  values <- c(5, 10, 4, 4, 20, 20, 12, 5)
  ifelse(x %in% 1:length(values), values[x], 0)
}
prob_dist <- data.frame(x = 1:8, prob = pd(1:8))
prob_dist


## ----echo=FALSE----------------------------------------------------------------------------------------------------
prob_dist %>%
  ggplot(aes(x = x, y = prob)) + 
  geom_col(width = 0.3) + 
  labs(x = 'x', y = 'Probability')


## ------------------------------------------------------------------------------------------------------------------
random_walk <- function(pd, start, num_steps){
  y <- rep(0, num_steps)
  current <- start
  for (j in 1:num_steps){
    candidate <- current + sample(c(-1, 1), 1)
    prob <- pd(candidate) / pd(current)
    if (runif(1) < prob) current <- candidate
    y[j] <- current
  }
  return(y)
}


## ------------------------------------------------------------------------------------------------------------------
out <- random_walk(pd, 4, 10000)
head(out)
tail(out)


## ----echo = FALSE--------------------------------------------------------------------------------------------------
sim <- data.frame(out) %>% 
  group_by(out) %>% 
  summarize(N = n(), prob = N / 10000) %>%
  select(out, prob) %>%
  add_column(Type = 'Simulated')

truth <- data.frame(out = 1:8, prob = pd(1:8)/sum(pd(1:8)), Type = 'Actual')

bind_rows(sim,truth) %>%
  ggplot() +
  aes(x = out, y = prob, fill = Type) +
  geom_col(width = 0.5, position = "dodge") + 
  labs(x = 'x', y = 'Probability')


## ---- out.width = '9cm',out.height='5cm',fig.align='center',echo=FALSE---------------------------------------------
knitr::include_graphics('img/plummer.png') 


## ---- out.width = '10cm',out.height='8cm',fig.align='center',echo=FALSE--------------------------------------------
knitr::include_graphics('img/stork_world.png')    


## ------------------------------------------------------------------------------------------------------------------
nbchicks <- c(151,105,73,107,113,87,77,108,118,122,112,120,122,89,69,71,
              53,41,53,31,35,14,18)

nbpairs <- c(173,164,103,113,122,112,98,121,132,136,133,137,145,117,90,80,
            67,54,58,39,42,23,23)

temp <- c(15.1,13.3,15.3,13.3,14.6,15.6,13.1,13.1,15.0,11.7,15.3,14.4,14.4,
         12.7,11.7,11.9,15.9,13.4,14.0,13.9,12.9,15.1,13.0)

rain <- c(67,52,88,61,32,36,72,43,92,32,86,28,57,55,66,26,28,96,48,90,86,
           78,87)

datax <- list(N = 23, nbchicks = nbchicks, nbpairs = nbpairs, 
              temp = (temp - mean(temp))/sd(temp), 
              rain = (rain - mean(rain))/sd(rain))


## ---- echo=TRUE, eval=FALSE----------------------------------------------------------------------------------------
## {
## # Likelihood
##   	for( i in 1 : N){
## 		nbchicks[i] ~ dbin(p[i],nbpairs[i])
## 		logit(p[i]) <- a + b.temp * temp[i] + b.rain * rain[i]
## 		}
## # ...


## ---- echo=TRUE, eval=FALSE----------------------------------------------------------------------------------------
## # Priors
## a ~ dnorm(0,0.001)
## b.temp ~ dnorm(0,0.001)
## b.rain ~ dnorm(0,0.001)
## }


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------
model <- 
paste("
model
{
	for( i in 1 : N) 
		{
		nbchicks[i] ~ dbin(p[i],nbpairs[i])
		logit(p[i]) <- a + b.temp * temp[i] + b.rain * rain[i]
		}
a ~ dnorm(0,0.001)
b.temp ~ dnorm(0,0.001)
b.rain ~ dnorm(0,0.001)
	}
")
writeLines(model,"code/logistic.txt")


## ---- message=FALSE, warning=FALSE, eval = FALSE-------------------------------------------------------------------
## logistic <- function() {
## 	for( i in 1 : N)
## 		{
## 		nbchicks[i] ~ dbin(p[i],nbpairs[i])
## 		logit(p[i]) <- a + b.temp * temp[i] + b.rain * rain[i]
## 		}
## 			
## # priors for regression parameters
## a ~ dnorm(0,0.001)
## b.temp ~ dnorm(0,0.001)
## b.rain ~ dnorm(0,0.001)
## 	}


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------
# list of lists of initial values (one for each MCMC chain)
init1 <- list(a = -0.5, b.temp = -0.5, b.rain = -0.5)
init2 <- list(a = 0.5, b.temp = 0.5, b.rain = 0.5)
inits <- list(init1,init2)

# specify parameters that need to be estimated
parameters <- c("a","b.temp","b.rain")

# specify nb iterations for burn-in and final inference 
nb.burnin <- 1000
nb.iterations <- 2000


## ----eval = FALSE--------------------------------------------------------------------------------------------------
## # load R2jags
## library(R2jags)
## # run Jags
## storks <- jags(data  = datax,
##                inits = inits,
##                parameters.to.save = parameters,
##                model.file = "code/logistic.txt",
##                # model.file = logistic, # if a function was written
##                n.chains = 2,
##                n.iter = nb.iterations,
##                n.burnin = nb.burnin)
## storks


## ----echo=FALSE, message=FALSE, warning=FALSE----------------------------------------------------------------------
library(R2jags)
storks <- jags(data  = datax,
               inits = inits,
               parameters.to.save = parameters,
               model.file = "code/logistic.txt",
               n.chains = 2,
               n.iter = nb.iterations,
               n.burnin = nb.burnin)
storks


## ---- echo=FALSE---------------------------------------------------------------------------------------------------
knitr::include_graphics('img/acf.png')   


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------
traceplot(storks,mfrow = c(1, 2), varname = c('b.rain','b.temp'), ask = FALSE)


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------
autocorr.plot(as.mcmc(storks),ask = FALSE) 


## ---- out.width = '11cm',out.height='7cm',fig.align='center',echo=FALSE--------------------------------------------
knitr::include_graphics('img/mcmc.png')   


## ------------------------------------------------------------------------------------------------------------------
storks


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------
res <- as.mcmc(storks) # convert outputs in a list
res <- rbind(res[[1]],res[[2]]) # put two MCMC lists on top of each other
head(res)


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------
# probability that the effect of rainfall is negative
mean(res[,'b.rain'] < 0)


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------
# probability that the effect of temperature is negative
mean(res[,'b.temp'] < 0)


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------
quantile(res[,'b.rain'],c(0.025,0.975))


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------
quantile(res[,'b.temp'],c(0.025,0.975))


## ----echo=FALSE, message=FALSE, warning=FALSE----------------------------------------------------------------------
par(mfrow=c(1,2))
plot(density(res[,'b.rain']),xlab="",ylab="", main="Rainfall",lwd=3)
abline(v=0,col='red',lwd=2)
plot(density(res[,'b.temp']),xlab="",ylab="", main="Temperature",lwd=3)
abline(v=0,col='red',lwd=2)


## ------------------------------------------------------------------------------------------------------------------
stupid_pd <- res[,'b.rain']^2 + cos(res[,'b.temp'])
head(stupid_pd)


## ------------------------------------------------------------------------------------------------------------------
plot(density(stupid_pd), xlab = '', main = '', lwd = 3)


## ----include=FALSE-------------------------------------------------------------------------------------------------
model <- 
paste("
model
{
	for( i in 1 : N) 
		{
		nbchicks[i] ~ dbin(p[i],nbpairs[i])
		logit(p[i]) <- a + b.temp * temp[i] + b.rain * rain[i]
		}
			
# priors for regression parameters
a ~ dnorm(0,0.001)
b.temp ~ dnorm(0,0.001)
b.rain ~ dnorm(0,0.001)
			
	}
")
writeLines(model,"code/logistic.txt")
init1 <- list(a = -0.5, b.temp = -0.5, b.rain = -0.5)
init2 <- list(a = 0.5, b.temp = 0.5, b.rain = 0.5)
inits <- list(init1,init2)
parameters <- c("a","b.temp","b.rain")
nb.burnin <- 1000
nb.iterations <- 2000
storks <- jags(data  = datax,
               inits = inits,
               parameters.to.save = parameters,
               model.file = "code/logistic.txt",
               n.chains = 2,
               n.iter = nb.iterations,
               n.burnin = nb.burnin)


## ------------------------------------------------------------------------------------------------------------------
# calculate wAIC with JAGS
# https://sourceforge.net/p/mcmc-jags/discussion/610036/thread/8211df61/#ea5c
samples <- jags.samples(storks$model,c("WAIC","deviance"), type = "mean", 
						n.iter = 2000,
						n.burnin = 1000,
						n.thin = 1)


## ------------------------------------------------------------------------------------------------------------------
samples$p_waic <- samples$WAIC
samples$waic <- samples$deviance + samples$p_waic
tmp <- sapply(samples, sum)
waic <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
waic


## ------------------------------------------------------------------------------------------------------------------
# model specification
model <- 
paste("
model
{
	for( i in 1 : N) 
		{
		nbchicks[i] ~ dbin(p[i],nbpairs[i])
		logit(p[i]) <- a + b * cov[i]
		}
			
# priors for regression parameters
a ~ dnorm(0,0.001)
b ~ dnorm(0,0.001)
	}
")
writeLines(model,"code/logtemp.txt")


## ------------------------------------------------------------------------------------------------------------------
# list of lists of initial values (one for each MCMC chain)
init1 <- list(a = -0.5, b = -0.5)
init2 <- list(a = 0.5, b = 0.5)
inits <- list(init1,init2)
# specify parameters that need to be estimated
parameters <- c("a","b")
# specify nb iterations for burn-in and final inference 
nb.burnin <- 1000
nb.iterations <- 2000
# read in data
datax <- list(N = 23, nbchicks = nbchicks, nbpairs = nbpairs, 
              cov = (temp - mean(temp))/sd(temp))


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------
# load R2jags to run Jags through R
storks_temp <- jags(data  = datax,
               inits = inits,
               parameters.to.save = parameters,
               model.file = "code/logtemp.txt",
               n.chains = 2,
               n.iter = nb.iterations,
               n.burnin = nb.burnin)


## ------------------------------------------------------------------------------------------------------------------
# compute WAIC
samples <- jags.samples(storks_temp$model,c("WAIC","deviance"), type = "mean", 
						n.iter = 2000,
						n.burnin = 1000,
						n.thin = 1)
samples$p_waic <- samples$WAIC
samples$waic <- samples$deviance + samples$p_waic
tmp <- sapply(samples, sum)
waic_temp <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)


## ----message=FALSE, warning=FALSE, paged.print=FALSE---------------------------------------------------------------
# read in data
datax <- list(N = 23, nbchicks = nbchicks, nbpairs = nbpairs, 
              cov = (rain - mean(rain))/sd(rain))


## ----message=FALSE, warning=FALSE, paged.print=FALSE---------------------------------------------------------------
# load R2jags to run Jags through R
storks_temp <- jags(data  = datax,
               inits = inits,
               parameters.to.save = parameters,
               model.file = "code/logtemp.txt",
               n.chains = 2,
               n.iter = nb.iterations,
               n.burnin = nb.burnin)


## ----message=FALSE, warning=FALSE, paged.print=FALSE---------------------------------------------------------------
# compute WAIC
samples <- jags.samples(storks_temp$model,c("WAIC","deviance"), type = "mean", 
						n.iter = 2000,
						n.burnin = 1000,
						n.thin = 1)
samples$p_waic <- samples$WAIC
samples$waic <- samples$deviance + samples$p_waic
tmp <- sapply(samples, sum)
waic_rain <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)


## ------------------------------------------------------------------------------------------------------------------
# model specification
model <- 
paste("
model
{
	for( i in 1 : N) 
		{
		nbchicks[i] ~ dbin(p[i],nbpairs[i])
		logit(p[i]) <- a
		}
			
# priors for regression parameters
a ~ dnorm(0,0.001)
	}
")
writeLines(model,"code/lognull.txt")


## ----message=FALSE, warning=FALSE, paged.print=FALSE---------------------------------------------------------------
# list of lists of initial values (one for each MCMC chain)
init1 <- list(a = -0.5)
init2 <- list(a = 0.5)
inits <- list(init1,init2)
# specify parameters that need to be estimated
parameters <- c("a")
# specify nb iterations for burn-in and final inference 
nb.burnin <- 1000
nb.iterations <- 2000
# read in data
datax <- list(N = 23, nbchicks = nbchicks, nbpairs = nbpairs)


## ----message=FALSE, warning=FALSE, paged.print=FALSE---------------------------------------------------------------
# load R2jags to run Jags through R
storks_temp <- jags(data  = datax,
               inits = inits,
               parameters.to.save = parameters,
               model.file = "code/lognull.txt",
               n.chains = 2,
               n.iter = nb.iterations,
               n.burnin = nb.burnin)


## ----message=FALSE, warning=FALSE, paged.print=FALSE---------------------------------------------------------------
# compute WAIC
samples <- jags.samples(storks_temp$model,c("WAIC","deviance"), type = "mean", 
						n.iter = 2000,
						n.burnin = 1000,
						n.thin = 1)
samples$p_waic <- samples$WAIC
samples$waic <- samples$deviance + samples$p_waic
tmp <- sapply(samples, sum)
waic_null <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)


## ------------------------------------------------------------------------------------------------------------------
data.frame(model = c('both_covariates', 'temp', 'rain', 'none'),
           waic = c(waic[1],waic_temp[1],waic_rain[1],waic_null[1]),
           p_waic = c(waic[2],waic_temp[2],waic_rain[2],waic_null[2])) %>%
  arrange(waic)


## ---- out.width = '8cm',out.height='6cm',fig.align='center',echo=FALSE---------------------------------------------
knitr::include_graphics('img/pdm.png')    


## ---- message=FALSE, warning=FALSE, include=FALSE------------------------------------------------------------------
# On lit le jeu de données à analyser et on le nettoie
VMG <- read.table("dat/VMG.csv", header=TRUE, dec= ".", sep =";")

# crée crée un vecteur contenant le nb de graines (en log)
y <- round(VMG$NGrTotest)

# crée un vecteur contenant la biomass
x <- VMG$Vm

# crée un vecteur contenant le nom des espèces
Sp <- VMG$Sp

# crée un vecteur contenant le numéro des espèces
species <- as.numeric(Sp)

# nombre d'espèces
nbspecies <- length(levels(Sp)) # ou bien length(unique(species))

# nombre de mesures
n <- length(y)

## ---- echo=FALSE, message=FALSE, warning=FALSE,fig.align='center',out.width = '12cm',out.height='7cm'--------------
# graphical representation
library(lattice)
dat <- data.frame(Biomass=x,nb_grain=y,Species=Sp)
xyplot(log(nb_grain) ~ Biomass | Species,data=dat,xlab='Biomass',ylab='Number of grains (log transformed)')


## ------------------------------------------------------------------------------------------------------------------
# read in data
VMG <- read.table("dat/VMG.csv", header=TRUE, dec= ".", sep =";")
# nb of seeds (log)
y <- VMG$NGrTotest
# biomass
x <- VMG$Vm
x <- (x - mean(x))/sd(x)
# species name
Sp <- VMG$Sp
# species label
species <- as.numeric(Sp)
# species name
nbspecies <- length(levels(Sp))
# total nb of measurements 
n <- length(y)


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------
model <- 
paste("
model{
for(i in 1:n){
	y[i] ~ dnorm(mu[i],tau.y)
	mu[i] <- a + b*x[i]
	}
tau.y <- 1/(sigma.y*sigma.y)
sigma.y ~ dunif(0,100)
a ~ dnorm(0,0.001)
b ~ dnorm(0,0.001)
}
")
writeLines(model,"code/completepooling.bug")


## ------------------------------------------------------------------------------------------------------------------
# data
allom.data <- list(y=y,n=n,x=x)

# initial values
init1 <- list(a=rnorm(1), b=rnorm(1),sigma.y=runif(1))
init2 <- list(a=rnorm(1), b=rnorm(1),sigma.y=runif(1))
inits <- list(init1,init2)

# parameters to be estimated
allom.parameters <- c("a", "b", "sigma.y")


## ------------------------------------------------------------------------------------------------------------------
allom.1 <- R2jags::jags(allom.data,inits,allom.parameters,
                        n.iter = 2500,model.file="code/completepooling.bug", 
                        n.chains = 2, n.burn = 1000)


## ------------------------------------------------------------------------------------------------------------------
allom.1


## ----echo=FALSE, message=FALSE, warning=FALSE----------------------------------------------------------------------
library(lattice)
xyplot(y ~ x | Sp,
       xlab = "Biomass", ylab = "Number of seeds",main="complete pooling",
       panel = function(x, y) {
           panel.xyplot(x, y)
           panel.abline(a=c(1991,1677),col='red',lwd=3)
       })


## ---- fig.align='center',echo=FALSE--------------------------------------------------------------------------------
knitr::include_graphics('img/varyingint.png')    


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------
model <- paste("
model {
  for (i in 1:n){
    y[i] ~ dnorm (mu[i], tau.y)
    mu[i] <- a[species[i]] + b *x[i]}
  tau.y <- pow(sigma.y, -2)
  sigma.y ~ dunif (0, 100)
  for (j in 1:nbspecies){ a[j] ~ dnorm (mu.a, tau.a)}
  mu.a ~ dnorm (0, 0.001)
  tau.a <- pow(sigma.a, -2)
  sigma.a ~ dunif (0, 100)
  b ~ dnorm (0, 0.001)    }")
writeLines(model,"code/varint.bug")


## ----message=FALSE, warning=FALSE, include=FALSE-------------------------------------------------------------------
allom.data <- list (n=n, nbspecies= nbspecies,x=x,y=y,species=species)

# on specifie le modele 
model <- 
paste("
model {
  for (i in 1:n){
    y[i] ~ dnorm (mu[i], tau.y)
    mu[i] <- a[species[i]] + b *x[i]
  }

  tau.y <- pow(sigma.y, -2)
  sigma.y ~ dunif (0, 100)

  for (j in 1:nbspecies){
    a[j] ~ dnorm (mu.a, tau.a)
  }
  
  mu.a ~ dnorm (0, 0.001)
  tau.a <- pow(sigma.a, -2)
  sigma.a ~ dunif (0, 100)

  b ~ dnorm (0, 0.001)

}
")
writeLines(model,"code/varint.bug")

init1 <- list(a=rnorm(nbspecies), b=rnorm(1), mu.a=rnorm(1),sigma.y=runif(1), sigma.a=runif(1))
init2 <- list(a=rnorm(nbspecies), b=rnorm(1), mu.a=rnorm(1),sigma.y=runif(1), sigma.a=runif(1))
inits<-list(init1,init2)
allom.parameters <- c ("a", "b", "mu.a","sigma.y", "sigma.a")
# run JAGS
allom.2 <- R2jags::jags(allom.data,inits,allom.parameters, n.iter = 2500,model.file="code/varint.bug", n.chains = 2, n.burn = 1000)
allom.2

## graph (correction BUG 2015)
acoef.sp <- allom.2$BUGSoutput$summary[1:33,1]
bcoef <- allom.2$BUGSoutput$summary[34,1]

# varying-intercept predicted values
yfit <- rep(0,length=n)
for (k in 1:n){yfit[k] <- acoef.sp[species[k]] + bcoef * x[k]}

# pooling model (no species effect) predicted values
ylinear <- rep(0,length=n)
for (k in 1:n){ylinear[k] <- 1991 + 1677 * x[k]}

## define function to fit observed and predicted values in species-specific panels 
panelfun2 <- 
  function(x, y, subscripts, ...){ 
           llines(x, lmhat[subscripts], type="p") # observed data
           llines(x, hat[subscripts], type="l", lty=1,col='green',lwd=3) # partial pooling
           llines(x, hat2[subscripts], type="l", lty=1,col='red',lwd=3) # no pooling 
} 

# assign observed and predicted values
lmhat <- y # observed data
hat <- yfit # partial pooling
hat2 <- ylinear # no pooling


## ----echo=FALSE, message=FALSE, warning=FALSE----------------------------------------------------------------------
# build a multipanel plot 
xyplot(y ~ x | Sp, panel=panelfun2,
       xlab="Biomass", 
       ylab="Number of seeds",      
       key = list(text = list(c("partial pooling", "no pooling")),
       lines = list(lwd = 3, col = c("green", "red"),
       type = c("l", "l"))))


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------
model <- paste("
model {
  for (i in 1:n){
    y[i] ~ dnorm (mu[i], tau.y)
    mu[i] <- a[species[i]] + b *x[i]}
  tau.y <- pow(sigma.y, -2)
  sigma.y ~ dunif(0, 100)
  for (j in 1:nbspecies){ a[j] ~ dnorm (0, 0.001)}
  b ~ dnorm (0,0.001)    }")
writeLines(model,"code/nopooling.bug")


## ----message=FALSE, warning=FALSE, include=FALSE-------------------------------------------------------------------
allom.data <- list (n=n, nbspecies= nbspecies,x=x,y=y,species=species)

# on specifie le modele 
model <- 
paste("
model {
  for (i in 1:n){
    y[i] ~ dnorm (mu[i], tau.y)
    mu[i] <- a[species[i]] + b *x[i]
  }

  tau.y <- pow(sigma.y, -2)
  sigma.y ~ dunif (0, 100)

  for (j in 1:nbspecies){
    a[j] ~ dnorm (0, 0.001)
  }
  
  b ~ dnorm (0, 0.001)

}
")
writeLines(model,"code/nopooling.bug")

init1 <- list(a=rnorm(nbspecies), b=rnorm(1), sigma.y=runif(1))
init2 <- list(a=rnorm(nbspecies), b=rnorm(1), sigma.y=runif(1))
inits<-list(init1,init2)
allom.parameters <- c ("a", "b","sigma.y")
# run JAGS
allom.3 <- R2jags::jags(allom.data,inits,allom.parameters, n.iter = 2500,model.file="code/nopooling.bug", n.chains = 2, n.burn = 1000)
allom.3

## graph (correction BUG 2015)
acoef.sp <- allom.3$BUGSoutput$summary[1:33,1]
bcoef <- allom.3$BUGSoutput$summary[34,1]

# fixed-effect predicted values
yfit2 <- rep(0,length=n)
for (k in 1:n){yfit2[k] <- acoef.sp[species[k]] + bcoef * x[k]}

# assign observed and predicted values
hat3 <- yfit2 # complete pooling

## define function to fit observed and predicted values in species-specific panels 
panelfun3 <- 
  function(x, y, subscripts, ...){ 
           llines(x, lmhat[subscripts], type="p") # observed data
           llines(x, hat[subscripts], type="l", lty=1,col='green',lwd=3) # partial pooling
           llines(x, hat2[subscripts], type="l", lty=1,col='red',lwd=3) # no pooling
           llines(x, hat3[subscripts], type="l", lty=1,col='blue',lwd=3) # complete pooling
}


## ----echo=FALSE, message=FALSE, warning=FALSE----------------------------------------------------------------------
# build a multipanel plot 
xyplot(y ~ x | Sp, panel=panelfun3,  
       xlab="Biomass", 
       ylab="Number of seeds",      
       key = list(text = list(c("partial pooling", "complete pooling", "no pooling")),
       lines = list(lwd = 3, col = c("green", "red", "blue"),
       type = c("l", "l", "l"))))


## ------------------------------------------------------------------------------------------------------------------
samples <- jags.samples(allom.1$model,c("WAIC","deviance"), type = "mean", 
                        n.iter = 2000, n.burnin = 1000, n.thin = 1)
samples$p_waic <- samples$WAIC; samples$waic <- samples$deviance + samples$p_waic
tmp <- sapply(samples, sum)
waic_completepooling <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)


## ------------------------------------------------------------------------------------------------------------------
samples <- jags.samples(allom.2$model,c("WAIC","deviance"), type = "mean", 
                        n.iter = 2000, n.burnin = 1000, n.thin = 1)
samples$p_waic <- samples$WAIC; samples$waic <- samples$deviance + samples$p_waic
tmp <- sapply(samples, sum)
waic_partialpooling <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)


## ------------------------------------------------------------------------------------------------------------------
samples <- jags.samples(allom.3$model,c("WAIC","deviance"), type = "mean", 
                        n.iter = 2000, n.burnin = 1000, n.thin = 1)
samples$p_waic <- samples$WAIC; samples$waic <- samples$deviance + samples$p_waic
tmp <- sapply(samples, sum)
waic_nopooling <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)


## ------------------------------------------------------------------------------------------------------------------
data.frame(model = c('no pooling', 'partial pooling', 'complete pooling'),
           waic = c(waic_nopooling[1],
                    waic_partialpooling[1],
                    waic_completepooling[1]),
           p_waic = c(waic_nopooling[2],
                      waic_partialpooling[2],
                      waic_completepooling[2])) %>%
  arrange(waic)


## ---- out.width = '13cm',out.height='5cm',fig.align='center',echo=FALSE--------------------------------------------
knitr::include_graphics('img/bayesian_evol.png')    


## ---- fig.align='center', echo=FALSE-------------------------------------------------------------------------------
knitr::include_graphics('img/whytwitter.png')    


## ------------------------------------------------------------------------------------------------------------------
transects <- 10
data <- NULL
for (tr in 1:transects){
  ref <- rnorm(1,0,.5) # random effect (intercept)
  t <- runif(1, 18,22) + runif(1,-.2,0.2)*1:20 # water temperature gradient
  ans <- exp(ref -14 + 1.8 * t - 0.045 * t^2) # Anemone gradient (expected response)
  an <- rpois(20, ans) # actual counts on 20 segments of the current transect
  data <- rbind(data,cbind(rep(tr, 20), t, an))
}


## ------------------------------------------------------------------------------------------------------------------
ref <- rnorm(1,0,.5) # random effect (intercept)
t <- runif(1, 18,22) + runif(1,-.2,0.2)*1:20 # water temperature gradient
plot(t,type='l')


## ------------------------------------------------------------------------------------------------------------------
ans <- exp(ref -14 + 1.8 * t - 0.045 * t^2) # Anemone gradient (expected response)
plot(t,log(ans),type='l')


## ------------------------------------------------------------------------------------------------------------------
data <- data.frame(Transect=data[,1],Temperature=data[,2],Anemones=data[,3])
plot(data$Temperature, data$Anemones)


## ------------------------------------------------------------------------------------------------------------------
data$Temp <- (data$Temperature - mean(data$Temperature))/sd(data$Temperature)
head(data)


## ------------------------------------------------------------------------------------------------------------------
model <- 
paste("
model {
  for (i in 1:n){
    count[i] ~ dpois(lambda[i])
    log(lambda[i]) <- a[transect[i]] + b[1] * x[i] + b[2] * pow(x[i],2)
  }
  for (j in 1:nbtransects){
    a[j] ~ dnorm (mu.a, tau.a)
  }
  mu.a ~ dnorm (0, 0.001)
  tau.a <- pow(sigma.a, -2)
  sigma.a ~ dunif (0, 100)
  b[1] ~ dnorm (0, 0.001)
  b[2] ~ dnorm (0, 0.001)
}
")
writeLines(model,"code/GLMMpoisson.bug")


## ------------------------------------------------------------------------------------------------------------------
dat <- list(n = nrow(data), 
            nbtransects = transects, 
            x = data$Temp, 
            count = data$Anemones, 
            transect = data$Transect)
init1 <- list(a=rnorm(transects), b=rnorm(2), mu.a=rnorm(1), sigma.a=runif(1))
init2 <- list(a=rnorm(transects), b=rnorm(2), mu.a=rnorm(1), sigma.a=runif(1))
inits <- list(init1,init2)
par <- c ("a", "b", "mu.a", "sigma.a")


## ------------------------------------------------------------------------------------------------------------------
fit <- jags(dat, inits, par, n.iter = 2500, model.file="code/GLMMpoisson.bug", n.chains = 2, n.burn = 1000)


## ------------------------------------------------------------------------------------------------------------------
fit


## ------------------------------------------------------------------------------------------------------------------
library(lme4)
fit_lme4 <- glmer(Anemones ~ Temp + I(Temp^2) + (1 | Transect), data=data, family=poisson)
fit_lme4


## ------------------------------------------------------------------------------------------------------------------
visreg::visreg(fit_lme4,xvar='Temp')

