

root <- "C:/Users/roland/Documents/iMEM Manuscript Data and Code/"
code <- "Code/"

source(paste0(root,code,"iMEM_functions.R"))
export <- T

###################################################################################################

### Here we're going to do a simulation to illustrate marginal iMEM asymptotics ###

calc.wts <- function(ssize,pmean,psigma,suppmeans,suppsigs,suppids,marg=T)
{
  ############################################################################
  #  Function to calculate iMEM weights as a function of sample size of all sources
  #   Inputs:
  #
  #   - ssize: sample size of all sources
  #   - pmean: primary source mean
  #   - psigma: primary source standard deviation
  #   - suppmeans: supplementary source means
  #   - suppsigs: supplementary source standard deviations (in same order as suppmeans)
  #   - suppids: supplementary source IDs (in same order as suppmeans)
  #   - marg: indicator variable for marginal or general algorithm
  #
  ############################################################################
  
  nsrc <- length(suppmeans)
  
  pdat <- c(pmean,psigma,ssize)
  suppsds <- suppsigs
  suppNs <- rep(ssize,length(suppmeans))
  
  if(marg) { fit.imem <- imem_marg_v2(pdat, suppmeans,suppsds, suppNs,suppids, prior='pi_e',final_grpsize=2) }
  if(!marg){ fit.imem <- imem_gen.rand_v2(pdat, suppmeans,suppsds, suppNs,suppids, prior='pi_e',grpsize=2,acc_prop=0.75,final_grpsize=2) }
  
  wts <- fit.imem$memlist[,c(1:(ncol(fit.imem$memlist)-4),ncol(fit.imem$memlist))]
  
  nms.zerowt <- suppids[!suppids %in% names(wts)]
  for(i in nms.zerowt){ wts[[i]] <- FALSE }
  
  wtmat <- merge(mods,wts,by=suppids,all.x=T)
  wtmat$postwts[is.na(wtmat$postwts)] <- 0
  wtmat <- wtmat[order(wtmat$modid),]
  wtmat <- wtmat[,c("modid","postwts")]
  
  out <- wtmat$postwts
  return(out)
}



### Set values for all function parameters ####
ssizes <- c(10:99, (10:99)*10,(10:99)*100,(10:99)*1000,(10:99)*10000,(10:100)*100000 )
psigma <- 4
pmean <- -2
suppsigs <- c(4,2,3,1)
suppmeans <- c(-2,-2,-2,4)
nsrc <- length(suppmeans)
suppids <- paste0("src",1:length(suppmeans))
mods <- expand.grid ( rep(list(c(FALSE,TRUE)), nsrc)); names(mods) <- paste0("src",1:nsrc) 
mods$modid <- 1:16


# calculate weights
wtmat.marg <- sapply(ssizes,calc.wts,pmean,psigma,suppmeans,suppsigs,suppids)


####### Plot of Weights by Sample Size ###############
cols <- c("blue","black","orange","black","orange","black","red",rep("black",9))
ltys <- c(2,1,3,1,3,1,4,rep(1,9))

if(export) { pdf( paste0(root,"Paper Graphics/Posterior Weights by Sample Size_v2.pdf")) }
plot("n",xlim=c(0,100000),ylim=c(0,1),xlab="Sample Sizes for all Sources",ylab="Posterior Model Weights")
for(i in 1:nrow(wtmat.marg)) { lines(ssizes,wtmat.marg[i,],col=cols[i],lty=ltys[i],lwd=1.5)}
legend(68000,0.9 ,c("Model 7","Models 3,5","Model 1","All Other Models"),col=cols[c(7,5,1,16)],lty=ltys[c(7,5,1,16)],lwd=1.5,bty="n")

text(40000,0.92,"M7",cex=0.8,col=cols[7])
text(20000,0.06,"M3",cex=0.8,col=cols[5])
text(28000,0.05,"M5",cex=0.8,col=cols[5])
text(4000,0.025,"M1",cex=0.8,col=cols[1])

if(export) { dev.off() }
