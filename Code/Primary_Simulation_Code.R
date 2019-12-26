


# Simulation settings:
#  - primary and supplementary sources will have n = 20, and sigma = 4
#  - all scenarios will have 9 potential supplmentary sources
#  - Scen 1: all supplementary sources have similar sample means (~-2)
#  - Scen 2: one group of 2 is similar (~4), another group of 7 is similar (~-2)
#  - Scen 3: one group of 4 is similar (~4), another group of 5 is similar (~-2)

#  - True mean of the primary source will be varied on a grid from (-6.5,8.5)
#  - 1000 samples will be simulated at each value of the true mean
#  - Models will be compared: fully joing MEM, marginal iMEM with 2-7 supplementary sources in final model



root <- "C:/Users/rbrow_000/Documents/UofM Documents/Dissertation/Iterated-Multisource-Exchangeability-Models/"
simdata <- "Data/Simulation Results/"
code <- "Code/"


source(paste0(root,code,"iMEM_functions.R"))

library(lme4)
library(merTools)
library(lmerTest)

# set simulation settings

# sample mean scenarios
xbars1 <- rep(-2,9)
xbars2 <- c(4,rep(-2,7),4)
xbars3 <- c(4,4,rep(-2,5),4,4)
xbars4 <- c(-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6)
xbars5 <- c(3.9, -2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7, 4.1)
xbars6 <- c(3.9, 3.8, -2.2,-2.1,-2,-1.9,-1.8, 4.2, 4.1)
xbars7 <- c(-2.5,-1.5,1,1.2,1.3,2.3,2.8,3,3.8)
xbars8 <- seq(-2,4,length.out=9)

mus <- seq(-6.5,8.5,length.out=151) # true mean grid
n <- 20 # sample size for all sources
sigma <- 4 # sd for all sources
nsim <- 1000 # number of simulations


##############################################################################################
##############################################################################################
##############################################################################################

#### Simulations for marginal iMEMs, full MEM, and simple mean #######

fnums <- 2:7 # final group size values to be considered for marginal iMEMs


# run simulations for each supplmentary source scenario,
set.seed(515)
scen1simres <- sim_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars1,fnums=fnums,samp.sig=F,nsim=nsim)
save(scen1simres,file=paste0(root,simdata,"scen1res_v2b.rdata"))

set.seed(515)
scen2simres <- sim_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars2,fnums=fnums,samp.sig=F,nsim=nsim)
save(scen2simres,file=paste0(root,simdata,"scen2res_v2b.rdata"))

set.seed(515)
scen3simres <- sim_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars3,fnums=fnums,samp.sig=F,nsim=nsim)
save(scen3simres,file=paste0(root,simdata,"scen3res_v2b.rdata"))

set.seed(515)
scen4simres <- sim_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars4,fnums=fnums,samp.sig=F,nsim=nsim)
save(scen4simres,file=paste0(root,simdata,"scen4res_v2b.rdata"))

set.seed(515)
scen5simres <- sim_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars5,fnums=fnums,samp.sig=F,nsim=nsim)
save(scen5simres,file=paste0(root,simdata,"scen5res_v2b.rdata"))

set.seed(515)
scen6simres <- sim_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars6,fnums=fnums,samp.sig=F,nsim=nsim)
save(scen6simres,file=paste0(root,simdata,"scen6res_v2b.rdata"))

set.seed(515)
scen7simres <- sim_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars7,fnums=fnums,samp.sig=F,nsim=nsim)
save(scen7simres,file=paste0(root,simdata,"scen7res_v2b.rdata"))

set.seed(515)
scen8simres <- sim_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars8,fnums=fnums,samp.sig=F,nsim=nsim)
save(scen8simres,file=paste0(root,simdata,"scen8res_v2b.rdata"))


# this takes ~33 seconds, meaning the full simulation for one scenario should take ~13 hours
# system.time(sim_allmu(mus=rep(-2,20),sigma=4,n=20,xbars=xbars1,fnums=2:7,nsim=5))



##############################################################################################
##############################################################################################
##############################################################################################

#### Simulations for random effects models #####

# For the random effects models, we need to generate some sample datasets that have the requisite 
# sample means and standard deviations.
# We do this by simulating dataset, standardizing the data, then adding desired sample mean
# and multiplying by desired standard deviation
suppdatgen <- function( sampmeans, sampsd, n )
{
  ninds <- length(sampmeans)
  
  outdat <- NULL
  
  for(i in 1:ninds)
  {
    id <- paste0("supp",i)
    tmp <- rnorm(n)
    tmp <- sampsd*(tmp-mean(tmp))/sd(tmp)+sampmeans[i]
    tmp <- data.frame(id=id,outcm=tmp)
    outdat <- rbind(outdat,tmp)
  }
  
  return(outdat)
}


# Generate Data for All Scenarios 
set.seed(1234)
scen1dat <- suppdatgen(xbars1,sigma,n)
scen2dat <- suppdatgen(xbars2,sigma,n)
scen3dat <- suppdatgen(xbars3,sigma,n)
scen4dat <- suppdatgen(xbars4,sigma,n)
scen5dat <- suppdatgen(xbars5,sigma,n)
scen6dat <- suppdatgen(xbars6,sigma,n)
scen7dat <- suppdatgen(xbars7,sigma,n)
scen8dat <- suppdatgen(xbars8,sigma,n)


system.time(simre_allmu(mus=rep(-2,20),sigma=4,n=20,suppdat=scen1dat,nsim=1))


# Fit the Models 

# run random effects simulations for each supplmentary source scenario,
set.seed(515)
scen1reres <- simre_allmu(mus=mus,sigma=sigma,n=n,suppdat=scen1dat,nsim=nsim)
save(scen1reres,file=paste0(root,simdata,"scen1reres_v2.rdata"))

set.seed(515)
scen2reres <- simre_allmu(mus=mus,sigma=sigma,n=n,suppdat=scen2dat,nsim=nsim)
save(scen2reres,file=paste0(root,simdata,"scen2reres_v2.rdata"))

set.seed(515)
scen3reres <- simre_allmu(mus=mus,sigma=sigma,n=n,suppdat=scen3dat,nsim=nsim)
save(scen3reres,file=paste0(root,simdata,"scen3reres_v2.rdata"))

set.seed(515)
scen4reres <- simre_allmu(mus=mus,sigma=sigma,n=n,suppdat=scen4dat,nsim=nsim)
save(scen4reres,file=paste0(root,simdata,"scen4reres_v2.rdata"))

set.seed(515)
scen5reres <- simre_allmu(mus=mus,sigma=sigma,n=n,suppdat=scen5dat,nsim=nsim)
save(scen5reres,file=paste0(root,simdata,"scen5reres_v2.rdata"))

set.seed(515)
scen6reres <- simre_allmu(mus=mus,sigma=sigma,n=n,suppdat=scen6dat,nsim=nsim)
save(scen6reres,file=paste0(root,simdata,"scen6reres_v2.rdata"))

set.seed(515)
scen7reres <- simre_allmu(mus=mus,sigma=sigma,n=n,suppdat=scen7dat,nsim=nsim)
save(scen7reres,file=paste0(root,simdata,"scen7reres_v2.rdata"))

set.seed(515)
scen8reres <- simre_allmu(mus=mus,sigma=sigma,n=n,suppdat=scen8dat,nsim=nsim)
save(scen8reres,file=paste0(root,simdata,"scen8reres_v2.rdata"))


##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################


#### Simulations for random general iMEMs ####

# general iMEM parameters to be considered
us <- 2:4         # group sizes
rs <- c(0.5,0.75) # acceptance probabilities
qs <- 4           # final group size


system.time(simGIMEM_rand_allmu(mus=rep(-2,20),sigma=4,n=20,xbars=xbars1,us=us,rs=rs,qs=qs,samp.sig=F,nsim=5))


# run simulations for each supplmentary source scenario,
set.seed(515)
scen1grandres <- simGIMEM_rand_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars1,us=us,rs=rs,qs=qs,samp.sig=F,nsim=nsim)
save(scen1grandres,file=paste0(root,simdata,"scen1grandres_v2.rdata"))

set.seed(515)
scen2grandres <- simGIMEM_rand_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars2,us=us,rs=rs,qs=qs,samp.sig=F,nsim=nsim)
save(scen2grandres,file=paste0(root,simdata,"scen2grandres_v2.rdata"))

set.seed(515)
scen3grandres <- simGIMEM_rand_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars3,us=us,rs=rs,qs=qs,samp.sig=F,nsim=nsim)
save(scen3grandres,file=paste0(root,simdata,"scen3grandres_v2.rdata"))

set.seed(515)
scen4grandres <- simGIMEM_rand_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars4,us=us,rs=rs,qs=qs,samp.sig=F,nsim=nsim)
save(scen4grandres,file=paste0(root,simdata,"scen4grandres_v2.rdata"))

set.seed(515)
scen5grandres <- simGIMEM_rand_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars5,us=us,rs=rs,qs=qs,samp.sig=F,nsim=nsim)
save(scen5grandres,file=paste0(root,simdata,"scen5grandres_v2.rdata"))

set.seed(515)
scen6grandres <- simGIMEM_rand_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars6,us=us,rs=rs,qs=qs,samp.sig=F,nsim=nsim)
save(scen6grandres,file=paste0(root,simdata,"scen6grandres_v2.rdata"))

set.seed(515)
scen7grandres <- simGIMEM_rand_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars7,us=us,rs=rs,qs=qs,samp.sig=F,nsim=nsim)
save(scen7grandres,file=paste0(root,simdata,"scen7grandres_v2.rdata"))

set.seed(515)
scen8grandres <- simGIMEM_rand_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars8,us=us,rs=rs,qs=qs,samp.sig=F,nsim=nsim)
save(scen8grandres,file=paste0(root,simdata,"scen8grandres_v2.rdata"))


##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################

#### Simulations for marignal iMEMs with score threshold-based selection of q ####
ths <- c(0.05,0.1,0.15,0.2)

tst <- sim_MIMEMthresh_allmu(mus=rep(-8.5,2),sigma=4,n=20,xbars=xbars1,fnums=ths,samp.sig=F,nsim=2)
system.time(sim_MIMEMthresh_allmu(mus=rep(-2,20),sigma=4,n=20,xbars=xbars1,fnums=ths,samp.sig=F,nsim=5))


# run simulations for each supplmentary source scenario,
set.seed(515)
scen1margthresh <- sim_MIMEMthresh_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars1,fnums=ths,samp.sig=F,nsim=nsim)
save(scen1margthresh,file=paste0(root,simdata,"scen1margthresh_v2.rdata"))

set.seed(515)
scen2margthresh <- sim_MIMEMthresh_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars2,fnums=ths,samp.sig=F,nsim=nsim)
save(scen2margthresh,file=paste0(root,simdata,"scen2margthresh_v2.rdata"))

set.seed(515)
scen3margthresh <- sim_MIMEMthresh_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars3,fnums=ths,samp.sig=F,nsim=nsim)
save(scen3margthresh,file=paste0(root,simdata,"scen3margthresh_v2.rdata"))

set.seed(515)
scen4margthresh <- sim_MIMEMthresh_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars4,fnums=ths,samp.sig=F,nsim=nsim)
save(scen4margthresh,file=paste0(root,simdata,"scen4margthresh_v2.rdata"))

set.seed(515)
scen5margthresh <- sim_MIMEMthresh_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars5,fnums=ths,samp.sig=F,nsim=nsim)
save(scen5margthresh,file=paste0(root,simdata,"scen5margthresh_v2.rdata"))

set.seed(515)
scen6margthresh <- sim_MIMEMthresh_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars6,fnums=ths,samp.sig=F,nsim=nsim)
save(scen6margthresh,file=paste0(root,simdata,"scen6margthresh_v2.rdata"))

set.seed(515)
scen7margthresh <- sim_MIMEMthresh_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars7,fnums=ths,samp.sig=F,nsim=nsim)
save(scen7margthresh,file=paste0(root,simdata,"scen7margthresh_v2.rdata"))

set.seed(515)
scen8margthresh <- sim_MIMEMthresh_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars8,fnums=ths,samp.sig=F,nsim=nsim)
save(scen8margthresh,file=paste0(root,simdata,"scen8margthresh_v2.rdata"))

##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################

#### Simulations for marignal iMEMs with score threshold-based selection of q ####
qmaxs=c(7,9)
qmins=c(2,3,4)

tst <- simMIMEM_multq_allmu(mus=rep(-2,2),sigma=4,n=20,xbars=xbars1,qmaxs=qmaxs,qmins=qmins,samp.sig=F,nsim=2)
system.time(simMIMEM_multq_allmu(mus=rep(-2,10),sigma=4,n=20,xbars=xbars1,qmaxs=qmaxs,qmins=qmins,samp.sig=F,nsim=10))


# run simulations for each supplmentary source scenario,
set.seed(515)
scen1multq <- simMIMEM_multq_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars1,qmaxs=qmaxs,qmins=qmins,samp.sig=F,nsim=nsim)
save(scen1multq,file=paste0(root,simdata,"scen1multq_v2.rdata"))

set.seed(515)
scen2multq <- simMIMEM_multq_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars2,qmaxs=qmaxs,qmins=qmins,samp.sig=F,nsim=nsim)
save(scen2multq,file=paste0(root,simdata,"scen2multq_v2.rdata"))

set.seed(515)
scen3multq <- simMIMEM_multq_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars3,qmaxs=qmaxs,qmins=qmins,samp.sig=F,nsim=nsim)
save(scen3multq,file=paste0(root,simdata,"scen3multq_v2.rdata"))

set.seed(515)
scen4multq <- simMIMEM_multq_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars4,qmaxs=qmaxs,qmins=qmins,samp.sig=F,nsim=nsim)
save(scen4multq,file=paste0(root,simdata,"scen4multq_v2.rdata"))

set.seed(515)
scen5multq <- simMIMEM_multq_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars5,qmaxs=qmaxs,qmins=qmins,samp.sig=F,nsim=nsim)
save(scen5multq,file=paste0(root,simdata,"scen5multq_v2.rdata"))

set.seed(515)
scen6multq <- simMIMEM_multq_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars6,qmaxs=qmaxs,qmins=qmins,samp.sig=F,nsim=nsim)
save(scen6multq,file=paste0(root,simdata,"scen6multq_v2.rdata"))

set.seed(515)
scen7multq <- simMIMEM_multq_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars7,qmaxs=qmaxs,qmins=qmins,samp.sig=F,nsim=nsim)
save(scen7multq,file=paste0(root,simdata,"scen7multq_v2.rdata"))

set.seed(515)
scen8multq <- simMIMEM_multq_allmu(mus=mus,sigma=sigma,n=n,xbars=xbars8,qmaxs=qmaxs,qmins=qmins,samp.sig=F,nsim=nsim)
save(scen8multq,file=paste0(root,simdata,"scen8multq_v2.rdata"))






