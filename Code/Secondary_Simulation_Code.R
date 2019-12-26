

# root <- "C:/Users/roland/Documents/Roland-Thesis/MEMs/"
root <- "C:/Users/rbrow_000/Documents/UofM Documents/Dissertation/Roland-Thesis/MEMs/"
data <- "data/"
sim1 <- "code/Simulations/Simulation 1/"
source(paste0(root,"code/gaussian mem_v1.r"))
source(paste0(root,sim1,"simulation1 funx_v5.R"))

# sample means
xbars <- c( -3,-2.9,-2.8,-2.7,-2.6,rep(-2,5), 0,0.2,0.4,0.6,0.8,rep(1.3,10),3,3.3,3.6,3.9,rep(4.5,3)    )


mus <- seq(-6.5,8.5,length.out=151) # true mean grid
n <- 20 # sample size for all sources
sigma <- 4 # sd for all sources
nsim <- 200 # number of simulations
qs <- 8          # final group size




# process 1
set.seed(515)
simb_marg <- sim_MIMEM_allmu(mus=mus,sigma=sigma,n=20,xbars=xbars,fnums=qs,samp.sig=F,nsim=nsim) 
save(simb_marg,file=paste0(root,sim1,"results/simb_marg_v2.rdata"))

set.seed(515)
simb_gen_2to5 <- simGIMEM_allmu(mus=mus,sigma=sigma,n=20,xbars=xbars,us=2:5,rs=c(0.5,0.75),qs=qs,samp.sig=F,nsim=nsim) 
save(simb_gen_2to5,file=paste0(root,sim1,"results/simb_gen_2to5_v2.rdata"))

# process 2
set.seed(515)
simb_genrand_2to5 <- simGIMEM_rand_allmu(mus=mus,sigma=sigma,n=20,xbars=xbars,us=2:5,rs=c(0.5,0.75),qs=qs,samp.sig=F,nsim=nsim) 
save(simb_genrand_2to5,file=paste0(root,sim1,"results/simb_genrand_2to5_v2.rdata"))

# process 3
set.seed(515)
simb_gen_6 <- simGIMEM_allmu(mus=mus,sigma=sigma,n=20,xbars=xbars,us=6,rs=c(0.5,0.75),qs=qs,samp.sig=F,nsim=nsim) 
save(simb_gen_6,file=paste0(root,sim1,"results/simb_gen_6_v2.rdata"))

# process 4
set.seed(515)
simb_genrand_6 <- simGIMEM_rand_allmu(mus=mus,sigma=sigma,n=20,xbars=xbars,us=6,rs=c(0.5,0.75),qs=qs,samp.sig=F,nsim=nsim) 
save(simb_genrand_6,file=paste0(root,sim1,"results/simb_genrand_6_v2.rdata"))


##################################################################################


#### timing tests
mus <- rep(-2,5)

# process 1
system.time( sim_MIMEM_allmu(mus=mus,sigma=2,n=20,xbars=xbars,fnums=qs,nsim=2) )
system.time( simGIMEM_allmu(mus=mus,sigma=2,n=20,xbars=xbars,us=2:5,rs=c(0.5,0.75),qs=qs,nsim=2) )


# process 2
system.time( simGIMEM_rand_allmu(mus=mus,sigma=2,n=20,xbars=xbars,us=2:5,rs=c(0.5,0.75),qs=qs,nsim=2) )


# process 3
system.time( simGIMEM_allmu(mus=mus,sigma=2,n=20,xbars=xbars,us=6,rs=c(0.5,0.75),qs=qs,nsim=2) )

# process 4
system.time( simGIMEM_rand_allmu(mus=mus,sigma=2,n=20,xbars=xbars,us=6,rs=c(0.5,0.75),qs=qs,nsim=2) )





