



############## FUNCTIONS FOR IMPLEMENTING STANDARD GAUSSIAN MEMS ########################## 

modpostmean <- function(modvec,prim, means, Vs)
{
  ########################################
  # calculate posterior mean for a single exchangeability model specified by modvec
  #
  #  modvec: logical vector specifying a particular exchangeability model
  #  prim: infomration on the primary source (format: c(mean,sd,N) )
  #  means: means for all supplementary sources (must be in same order as modvec)
  #  Vs: s^2/n for all supplmentary sources (must be in same order as modvec)
  #	
  ########################################
  
  # pull primary cohort information and calculate v = s^2/n
  m <- prim[1]; s <- prim[2]; n <- prim[3]
  v <- s^2/n
  
  # calculate posterior mean from Eq. 3.6 in paper
  pmean <- (m*prod(Vs^modvec) + sum( (v*means*modvec)/Vs )*prod(Vs^modvec) ) / 
    ( v*( sum(modvec/Vs)*prod(Vs^modvec) ) + prod(Vs^modvec) ) 
  
  return(pmean)
}

modpostvar <- function(modvec,prim, means, Vs)
{
  ########################################
  # calculate posterior variance for a single exchangeability model specified by modvec
  #
  #  modvec: logical vector specifying a particular exchangeability model
  #  prim: infomration on the primary source (format: c(mean,sd,N) )
  #  means: means for all supplementary sources (must be in same order as modvec)
  #  Vs: s^2/n for all supplmentary sources (must be in same order as modvec)
  #	
  ########################################
  
  # pull primary cohort information and calculate v = s^2/n
  m <- prim[1]; s <- prim[2]; n <- prim[3]
  v <- s^2/n
  
  # calculate posterior variance from Eq. 3.6 in paper
  pvar <- ( 1/v + sum( modvec/Vs ) )^-1	
  
  return( pvar )	
}



modmlik <- function(modvec,prim,means,Vs)
{
  ################################################################
  # calculate conditional marginal likelhiood for single exchangeability model specified by modvec 
  #
  #  modvec: logical vector specifying a particular exchangeability model
  #  prim: infomration on the primary source (format: c(mean,sd,N) )
  #  means: means for all supplementary sources (must be in same order as modvec)
  #  Vs: s^2/n for all supplmentary sources (must be in same order as modvec)
  #
  ################################################################
  
  # pull primary cohort information and calculate v = s^2/n
  m <- prim[1]; s <- prim[2]; n <- prim[3]
  v <- s^2/n
  
  # number of total sources
  nsrc <- length(means)
  
  #########
  
  # calculate the marginal likelihood given the exchangebility model specified by modvec (per Eq. 3.5 in paper)
  
  # calculate internal sum A needed for marginal likelihood computation
  As <- rep(NA,length(modvec))
  for( i in 1:length(modvec)) { As[i] <- sum( modvec[-i]/Vs[-i] ) }
  
  # calculate internal sum B needed for marginal likelihood computation
  Bs <- rep(NA,length(modvec))
  for(i in 1:length(modvec)) 
  {
    gtind <- (1:length(modvec)) > i
    
    modvec2 <- modvec[gtind]
    means2 <- means[gtind]
    Vs2 <- Vs[gtind]
    
    # calculate values for the sum within quantity B
    Cs <- rep(NA,length( (i+1):length(modvec) ))
    k <- 1
    for(j in (i+1):length(modvec)) 
    {
      Cs[k] <- sum( (modvec/Vs)[-c(i,j)] )
      k <- k+1
    }		
    
    Bs[i] <- sum( (modvec[i]*modvec2*(means[i]-means2)^2)/(Vs[i] + Vs2 + Vs[i]*Vs2*( 1/v + Cs) ))		
  }	
  
  mlik <- ( sqrt(2*pi)^(nsrc + 1 - sum(modvec) ) / sqrt( (1/v + sum(modvec/Vs))*prod( (1/Vs)^!modvec )) )*
    exp( -0.5*sum( (modvec*(m-means)^2)/(v+Vs+v*Vs*As ) + Bs ) )		
  
  return(mlik)	
}



modweights <- function(mmat,mliks,prim,means,Vs,prior)
{
  
  ################################################################
  # calculate weights for all exchangeability models specified in mmat 
  #
  #  mmat: logical matrix of all possible exchangeability models (rows: models, columns: sources)
  #  mliks: vector of conditional marginal likelhioods for each possible model (must be in same order as cols of mmat)
  #  prim: infomration on the primary source (format: c(mean,sd,N) )
  #  means: means for all supplementary sources (must be in same order as rows of mmat)
  #  Vs: s^2/n for all supplmentary sources (must be in same order as rows of mmat)
  #
  ################################################################
  
  # pull primary cohort information and calculate v = s^2/n
  m <- prim[1]; s <- prim[2]; n <- prim[3]
  v <- s^2/n
  
  # number of total sources
  nsrc <- length(means)
  
  #########
  
  # equal prior weights for pi_e
  if( prior=='pi_e') { modpriors <- 1 
  
  # for pi_n use Eq. 3.7 from paper
  }else if (prior=='pi_n') {
    
    # quantities needed for the source prior calcualtion
    As <- apply(mmat,1, function(modvec) { sum( modvec/Vs) } )
    Bs <- apply(mmat,1, function(modvec) { prod( (1/Vs)^(!modvec) ) } )
    Cs <- apply(mmat,1, function(modvec) { sum( modvec ) } )
    
    ### need to check if it's 1/s or 1/s^2, either typo in paper, 
    ### or mistake in Alex's code, currently set to match Alex's code
    Ds <- sqrt( (1/s + As)*Bs)/ ( sqrt(2*pi)^( nsrc + 1- Cs )) 
    
    # calculate source inclusion priors
    priors_1 <- apply( mmat,2,function(srcvec) { sum( srcvec*Ds )  } )
    priors_0 <- apply( mmat,2,function(srcvec) { sum( (!srcvec)*Ds) } )
    
    # calculate model inclusion priors from source inclusion priors
    modpriors <- apply(mmat,1,function(modvec){ prod(modvec*priors_1 + (!modvec)*priors_0) })
    
  } else { stop("Current support only for priors 'pi_e' and 'pi_n") }
  
  # calculate raw posterior weights
  rawwts <- modpriors*mliks
  
  # calculate scaled posterior weights
  wts <- rawwts/sum(rawwts)
  
  return(wts) 	
}


derivexp <- function(modvec,prim, means, Vs)
{
  ########################################
  # calculate derivative of exponential term in the marginal likelhiood
  #  for single exchangeability model specified by modvec
  #
  #  modvec: logical vector specifying a particular exchangeability model
  #  prim: infomration on the primary source (format: c(mean,sd,N) )
  #  means: means for all supplementary sources (must be in same order as modvec)
  #  Vs: s^2/n for all supplmentary sources (must be in same order as modvec)
  #	
  ########################################
  
  # pull primary cohort information and calculate v = s^2/n
  m <- prim[1]; s <- prim[2]; n <- prim[3]
  v <- s^2/n
  
  # calculate internal sum A 
  As <- rep(NA,length(modvec))
  for( i in 1:length(modvec)) { As[i] <- sum( modvec[-i]/Vs[-i] ) }
  
  # calcualte
  derexp <- sum ( (-modvec*(m-means))/(v+Vs+v*Vs*As ) )
  return( derexp )	
}



derivwts <- function(mods)
{
  ########################################
  # calculate derivative of weights for single exchangeability model specified by modvec
  #
  #  mods: the 
  #	
  ########################################
  
  num <- rep(NA,nrow(mods))
  for(i in 1:length(num) ) { num[i] <- sum( mods$postwts*(mods$expder- mods$expder[i]) )}
  return(num/mods$postwts)
  
}




derivpmean <- function(modvec,prim,Vs)
{
  # pull primary cohort information and calculate v = s^2/n
  m <- prim[1]; s <- prim[2]; n <- prim[3]
  v <- s^2/n
  
  return( prod(Vs^modvec)/( v*( sum(modvec/Vs)*prod(Vs^modvec) ) + prod(Vs^modvec) ) )
}





############################################################



mem_calc <- function( prim, means,sds, Ns, prior='pi_e' )
{
  ##########################################################################
  # This function carries out MEM calcualtions for the case of gaussian means for an arbitrary number
  # of supplmentary sources.  All calculations are closed form, but model space increases exponentially
  # with linearly increasing number of sources, so unsure how efficiently it will scale 
  # 
  #   prim: infomration on the primary source (format: c(mean,sd,N) )
  #   means: vector of means for supplementary sources
  #   sds: vector of standard deviations for supplementary sources
  #   Ns: vector of sample sizes for supplementary sources
  #   prior: type of prior
  #  
  # NOTE: means, sds, and Ns must all be the same length
  ############################################################################
  
  # if means, sds, and Ns are of different length, error out
  if( length(unique( c(length(means),length(sds),length(Ns)) )) > 1 ){ stop("'means','sds','Ns' must be of equal length")}
  
  # First, let's generate a matrix with all possible exchangebility models
  nsrc <- length(means)	
  mods <- expand.grid ( rep(list(c(FALSE,TRUE)), nsrc)); names(mods) <- paste0("src",1:nsrc) 
  mmat <- as.matrix ( mods )
  
  # calculate v = s^2/n for suppelementary sources
  Vs <- sds^2/Ns
  
  # Calculate model-specific posterior means
  mods$pmeans <- apply(mmat,1,modpostmean, prim, means,Vs)
  
  # Calculate model-specific posterior variances
  mods$pvars <- apply(mmat,1,modpostvar, prim, means,Vs)
  
  # calculate conditional marginal likelihood for each model
  mods$mlik <- apply(mmat,1,modmlik,prim,means,Vs)
  
  # calculate posterior weights for each model
  mods$postwts <- modweights(mmat,mods$mlik,prim,means,Vs,prior)
  
  ## calculate overall posterior mean ##
  pmean <- sum(mods$pmeans*mods$postwts)
  
  ### calcualte overall posterior variance ###
  mods$expder <- apply(mmat,1,derivexp,prim,means,Vs) # derivative of exponential part of marginal likelhiood
  mods$wtsder <- derivwts(mods) # derivative of weight calculation
  mods$pmeander <- apply(mmat,1,derivpmean,prim,Vs) # derivative of posterior mean
  mods$gder <- ( (mods$pmeander/mods$postwts) - (mods$pmeans*mods$wtsder) )/ (1/mods$postwts)^2
  gmat <- mods$gder%*%t(mods$gder)
  pvar <- rep(1,nrow(mods))%*%gmat%*%rep(1,nrow(mods))* ( (prim[2]^2)/prim[3])
  
  ### calculate effective historical sample sizes for individual models, and overall
  mods$ess <- (prim[2]^2/mods$pvars)-prim[3]
  ess_fin <- sum(mods$ess*mods$postwts)
  
  
  #############################
  
  # generate output list and return
  srcinfo <- cbind( prim, rbind(means,sds,Ns) )
  colnames(srcinfo) <- c("primary",paste0("src",1:(ncol(srcinfo)-1)))
  rownames(srcinfo) <- c("mean","sd","N")
  
  out <- list(postmean=pmean,postvar=pvar,ess=ess_fin,
              memlist=mods[,c(grep("src",names(mods),value=T),"pmeans","pvars","ess","postwts")],
              srcinfo=srcinfo)
  return(out)
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################



############## FUNCTIONS FOR IMPLEMENTING VARIOUS IMEM FLAVORS ########################## 


####### MARGINAL IMEM FUNCTIONS #######

imem_marg <- function(prim, means,sds, Ns, prior='pi_e',final_grpsize,qthresh)
{
  ####################################################################
  # function to fit marginal selection iterated MEM 
  #
  # prim: primary source information c(mean,sd,N)
  # means: supplementary  means
  # sds: supplementary sds
  # Ns: supplementary sample sizes
  # prior: prior desired, 'pi_e' or 'pi_n'
  # final_grpsize: q, the number of sources in the final MEM model, must be <= # supp. sources
  #
  #
  ###################################################################
  
  
  # final_grpsize <- 3
  # prim <- c(-2,4,20)
  # means <- xbars3
  # sds <- rep(sigma,length(means))
  # Ns <- rep(n,length(means))
  
  # first fit all marginal MEM models to be used in mem_calc()
  margmems <- lapply(1:length(means),function(x) { mem_calc(prim=prim,means=means[x],sds=sds[x],Ns=Ns[x],prior='pi_e' )} )
  
  # now extract the scores and order them
  margscores <- sapply(1:length(margmems), function(x){ margmems[[x]]$memlist$postwts[2] })
  ord <- rev(order(margscores))
  
  # if qthresh is specified
  if(final_grpsize=="thresh") {
    final_grpsize <- sum(margscores >= qthresh)
    if(final_grpsize > 12) { final_grpsize <- 12}
  }
  
  
  # print(final_grpsize)
  
  # grab the sources with q highest scores
  finmeans <- means[ord][1:final_grpsize]
  finsds <- sds[ord][1:final_grpsize]
  finNs <- Ns[ord][1:final_grpsize]
  finScores <- margscores[ord][1:final_grpsize]
  
  # fit the final mem model
  finmem <- mem_calc(prim=prim,means=finmeans,sds=finsds,Ns=finNs,prior=prior)
  finmem$srcinfo <- rbind(finmem$srcinfo,c(NA,finScores) )
  rownames(finmem$srcinfo)[4] <- "score"
  
  # if final_grpsize is zero
  if(final_grpsize==0){
    
    finmem$postmean <- finmem$memlist$pmeans[1]
    finmem$postvar <- finmem$memlist$pvars[1]
    finmem$ess <- 0
    finmem$memlist <- finmem$memlist[1,c("pmeans","pvars","ess","postwts")]
    finmem$memlist$postwts <- 1
    finmem$srcinfo <- as.matrix( finmem$srcinfo[,1],colnames="primary")
    colnames(finmem$srcinfo) <- "primary"
  }
  
  return(finmem)
}


imem_marg_v2 <- function(prim, means,sds, Ns, IDs,prior='pi_e',final_grpsize)
{
  ####################################################################
  # Same functionality as imem_marg(), but additionally retains the IDs of the sources
  #  selected into the final model.  Additional input variable is `IDs`, which are 
  #  identifiers for the supplementary sources
  #
  ###################################################################
  
  
  # final_grpsize <- 3
  # prim <- c(-2,4,20)
  # means <- xbars3
  # sds <- rep(sigma,length(means))
  # Ns <- rep(n,length(means))
  
  # first fit all marginal MEM models to be used in mem_calc()
  margmems <- lapply(1:length(means),function(x) { mem_calc(prim=prim,means=means[x],sds=sds[x],Ns=Ns[x],prior='pi_e' )} )
  
  # now extract the scores and order them
  margscores <- sapply(1:length(margmems), function(x){ margmems[[x]]$memlist$postwts[2] })
  ord <- rev(order(margscores))
  
  # grab the highest "final_grpsize" highest scores
  finmeans <- means[ord][1:final_grpsize]
  finsds <- sds[ord][1:final_grpsize]
  finNs <- Ns[ord][1:final_grpsize]
  finIDs <- IDs[ord][1:final_grpsize]
  
  # fit the final mem model
  finmem <- mem_calc(prim=prim,means=finmeans,sds=finsds,Ns=finNs,prior=prior)
  finmem$ids.final <- finIDs
  colnames(finmem$srcinfo) <- c("primary",finIDs)
  colnames(finmem$memlist)[1:(ncol(finmem$memlist)-4)] <- finIDs
  
  return(finmem)
}



imem_marg_multq <- function(prim, means,sds, Ns, prior='pi_e',qmax,qmin)
{
  ####################################################################
  # function to fit marginal selection iterated MEM for multiple q values
  #
  # prim: primary source information c(mean,sd,N)
  # means: supplementary  means
  # sds: supplementary sds
  # Ns: supplementary sample sizes
  # prior: prior desired, 'pi_e' or 'pi_n'
  # final_grpsize: q, the number of sources in the final MEM model, must be <= # supp. sources
  #
  #
  ###################################################################
  
  
  # final_grpsize <- 3
  # prim <- c(-2,4,20)
  # means <- xbars3
  # sds <- rep(sigma,length(means))
  # Ns <- rep(n,length(means))
  
  # first fit all marginal MEM models to be used in mem_calc()
  margmems <- lapply(1:length(means),function(x) { mem_calc(prim=prim,means=means[x],sds=sds[x],Ns=Ns[x],prior='pi_e' )} )
  
  # now extract the scores and order them
  margscores <- sapply(1:length(margmems), function(x){ margmems[[x]]$memlist$postwts[2] })
  ord <- rev(order(margscores))
  
  
  
  outlist <- vector(mod="list",length=qmax-qmin+1)
  j <- 1
  for(i in qmax:qmin)
  {
    # grab the sources with q highest scores
    finmeans <- means[ord][1:i]
    finsds <- sds[ord][1:i]
    finNs <- Ns[ord][1:i]
    finScores <- margscores[ord][1:i]
    
    # fit the final mem model
    finmem <- mem_calc(prim=prim,means=finmeans,sds=finsds,Ns=finNs,prior=prior)
    finmem$srcinfo <- rbind(finmem$srcinfo,c(NA,finScores) )
    rownames(finmem$srcinfo)[4] <- "score"
    outlist[[j]] <- finmem
    j <- j+1
  }
  
  names(outlist) <- paste0("q=",qmax:qmin)
  
  # posterior mean is average of all means
  meanavg <- mean( sapply(outlist,function(x) x$postmean) )
  
  # measurement of posterior variance is maximum posterior variance from all models
  varmax <- max( sapply(outlist,function(x) x$postvar) )
  
  
  out <- list(postmean=meanavg,postvar=varmax,modlist=outlist)
  
  return(out)

}

##########################################################################
##########################################################################
##########################################################################

####### GENERAL IMEM FUNCTIONS #######

imem_gen.internal <- function(grpids,prim,suppdat,prior=prior)
{
  ####################################################################
  #  Internal function for fitting group-specific MEMs in the general iMEM algorithm.  
  ###################################################################
  
  # prim <- prim
  # grpids <- grps[[1]]
  # suppdat <- suppdat
  # prior <- "pi_e"
  
  
  ###################
  # corrected to do the right labeling of ids on the last line
  
  tt <- suppdat$suppids %in% grpids
  suppmeans <- suppdat$means[tt]
  suppsds <- suppdat$sds[tt]
  suppNs <- suppdat$Ns[tt]
  
  out <- mem_calc( prim=prim,means=suppmeans,sds=suppsds,Ns=suppNs,prior=prior)
  out$suppids <- suppdat$suppids[tt]
  return(out)
}


scorecalc <- function(mem_genout){
  
  ####################################################################
  #  Internal function for for calculating similarity scores in general iMEM algorithm
  ###################################################################
  
  nsupp <- length(mem_genout$suppids)
  wtmat <- mem_genout$memlist[,c(1:nsupp,ncol(mem_genout$memlist))]
  scores <- sapply( names(wtmat)[1:nsupp],function(x) { sum(wtmat$postwts[wtmat[[x]]]) } )
  names(scores) <- mem_genout$suppids
  return(scores)
}


imem_gen <- function(prim,means,sds,Ns,suppids,order.start,prior='pi_e',grpsize=5,acc_prop=0.5,final_grpsize=11)
{
  ####################################################################
  # Function to fit ordered general iMEM algorithm
  #
  # prim: primary source information c(mean,sd,N)
  # means: supplementary  means
  # sds: supplementary sds
  # Ns: supplementary sample sizes
  # suppids: supplementary source IDs
  # order.start: an intial ordering for the supplementary sources
  # grpsize: the group size u
  # acc_prop: the acceptnace proprotion r
  # prior: prior desired, 'pi_e' or 'pi_n'
  # final_grpsize: q, the number of sources in the final MEM model, must be <= # supp. sources
  #
  #
  ###################################################################
  
  # grpsize <- 2
  # final_grpsize <- 4
  # prim <- c(-2,4,20)
  # means <- rep(xbars6)
  # sds <- rep(sigma,length(means))
  # Ns <- rep(n,length(means))
  # suppids <- paste0("id",1:length(sds))
  # order.start <- 1:length(sds)
  # acc_prop <- 0.5
  
  suppdat <- data.frame(means=means,sds=sds,Ns=Ns,suppids=suppids)
  
  
  # set the temporary supplementary ids variable that will change after each iteration
  suppids_tmp <- suppids[order.start]
  
  # iterate mem calculations on groups of size grpsize to determine the most important sources
  repeat{
    
    # calculate number of groups based on groupsize
    ngrp <- floor(length(suppids_tmp)/grpsize)
    
    # assign ids to groups 
    asnmnt <- sort( rep(1:ngrp, length.out=length(suppids_tmp) ))
    grps <- lapply(1:ngrp,function(x){ suppids_tmp[asnmnt==x] } )
    
    # fit mems to each group
    memlist <- lapply( grps,imem_gen.internal,prim,suppdat,prior)
    
    # extract source scores for each mem, involves summing weights over all models that include a source
    scores <-  unlist(lapply(memlist,scorecalc) )
    
    # order scores
    scores <- scores[rev(order(scores))]
    
    # for now, take the ids in the top half of scores for the next iteration
    nacc <- ceiling( length(scores)*(acc_prop) )
    if(nacc <= final_grpsize | grpsize==1) { break }
    
    suppids_tmp <- names(scores)[1:nacc]
  }
  
  suppids_tmp <- names(scores)[1:final_grpsize]
  
  # fit the final mem
  finmod <- imem_gen.internal(suppids_tmp,prim,suppdat,prior)
  return(finmod)
}



imem_gen.rand <- function(prim,means,sds,Ns,suppids,prior='pi_e',grpsize=5,acc_prop=0.5,final_grpsize=11)
{
  ####################################################################
  # Function to fit random general iMEM algorithm
  #
  # prim: primary source information c(mean,sd,N)
  # means: supplementary  means
  # sds: supplementary sds
  # Ns: supplementary sample sizes
  # suppids: supplementary source IDs
  # grpsize: the group size u
  # acc_prop: the acceptnace proprotion r
  # prior: prior desired, 'pi_e' or 'pi_n'
  # final_grpsize: q, the number of sources in the final MEM model, must be <= # supp. sources
  #
  #
  ###################################################################
  
  # ssize <- 20
  # sigma <- 4
  # pmean <- -2
  # suppmeans <- c(-2,-2,-2,4)
  # suppids <- paste0("src",1:length(suppmeans))
  # 
  # grpsize <- 2
  # final_grpsize <- 2
  # prim <- c(pmean,sigma,ssize)
  # means <- suppmeans
  # sds <- rep(sigma,length(suppmeans))
  # Ns <- rep(ssize,length(means))
  # # order.start <- 1:length(sds)
  # acc_prop <- 0.75
  # prior <- "pi_e"
  # 
  
  
  suppdat <- data.frame(means=means,sds=sds,Ns=Ns,suppids=suppids)
  
  
  # set the temporary supplementary ids variable that will change after each iteration
  suppids_tmp <- suppids
  
  # iterate mem calculations on groups of size grpsize to determine the most important sources
  repeat{
    
    # calculate number of groups based on groupsize
    ngrp <- floor(length(suppids_tmp)/grpsize)
    
    # assign ids to groups 
    suppids_smp <- sample(suppids_tmp,length(suppids_tmp),replace=F)
    asnmnt <- rep(1:ngrp, length.out=length(suppids_tmp) )
    grps <- lapply(1:ngrp,function(x){ suppids_smp[asnmnt==x] } )
    
    # fit mems to each group
    memlist <- lapply( grps,imem_gen.internal,prim,suppdat,prior)
    
    # extract source scores for each mem, involves summing weights over all models that include a source
    scores <-  unlist(lapply(memlist,scorecalc) )
    
    # order scores
    scores <- scores[rev(order(scores))]
    
    # for now, take the ids in the top half of scores for the next iteration
    ceil.val <- ceiling( length(scores)*(acc_prop) )
    nacc <- ifelse(ceil.val==length(suppids_smp),ceil.val-1,ceil.val)
    if(nacc <= final_grpsize | grpsize==1) { break }
    
    suppids_tmp <- names(scores)[1:nacc]
  }
  
  suppids_tmp <- names(scores)[1:final_grpsize]
  
  # fit the final mem
  finmod <- imem_gen.internal(suppids_tmp,prim,suppdat,prior)
  return(finmod)
}



imem_gen.rand_v2 <- function(prim,means,sds,Ns,suppids,prior='pi_e',grpsize=5,acc_prop=0.5,final_grpsize=11)
{
  ####################################################################
  # Same functionality as imem_gen.rand(), but additionally retains the supplemtnary source
  #  IDs that are included in the final model
  #
  #
  ###################################################################
  
  # ssize <- 20
  # sigma <- 4
  # pmean <- -2
  # suppmeans <- c(-2,-2,-2,4)
  # suppids <- paste0("src",1:length(suppmeans))
  # 
  # grpsize <- 2
  # final_grpsize <- 2
  # prim <- c(pmean,sigma,ssize)
  # means <- suppmeans
  # sds <- rep(sigma,length(suppmeans))
  # Ns <- rep(ssize,length(means))
  # # order.start <- 1:length(sds)
  # acc_prop <- 0.75
  # prior <- "pi_e"
  
  
  
  suppdat <- data.frame(means=means,sds=sds,Ns=Ns,suppids=suppids)
  
  
  # set the temporary supplementary ids variable that will change after each iteration
  suppids_tmp <- suppids
  
  # iterate mem calculations on groups of size grpsize to determine the most important sources
  repeat{
    
    # calculate number of groups based on groupsize
    ngrp <- floor(length(suppids_tmp)/grpsize)
    
    # assign ids to groups 
    suppids_smp <- sample(suppids_tmp,length(suppids_tmp),replace=F)
    asnmnt <- rep(1:ngrp, length.out=length(suppids_tmp) )
    grps <- lapply(1:ngrp,function(x){ suppids_smp[asnmnt==x] } )
    
    # fit mems to each group
    memlist <- lapply( grps,imem_gen.internal,prim,suppdat,prior)
    
    # extract source scores for each mem, involves summing weights over all models that include a source
    scores <-  unlist(lapply(memlist,scorecalc) )
    
    # order scores
    scores <- scores[rev(order(scores))]
    
    # for now, take the ids in the top half of scores for the next iteration
    ceil.val <- ceiling( length(scores)*(acc_prop) )
    nacc <- ifelse(ceil.val==length(suppids_smp),ceil.val-1,ceil.val)
    if(nacc <= final_grpsize | grpsize==1) { break }
    
    suppids_tmp <- names(scores)[1:nacc]
  }
  
  suppids_tmp <- names(scores)[1:final_grpsize]
  
  # fit the final mem
  finmod <- imem_gen.internal(suppids_tmp,prim,suppdat,prior)
  finmod$suppids <- as.character(finmod$suppids)
  colnames(finmod$srcinfo) <- c("primary",finmod$suppids)
  colnames(finmod$memlist)[1:(ncol(finmod$memlist)-4)] <- finmod$suppids
  return(finmod)
}

# extract posterior mean, posterior variacne, and effective supplemental sample size from a MEM object
memextr <- function(memobj,cred.lev=0.95) 
{ 
  # since this function uses the random number generator we must grab the Randomseed, and then reset it afterwards
  oldseed <- .GlobalEnv$.Random.seed
  
  # sample the multinomial distribution for the normal mixture (out of 10000000 draws)
  nsamp.mod <- c(rmultinom(1,100000,memobj$memlist$postwts))
  
  # sample the normal mixture
  post.dist <- unlist( apply( cbind(nsamp.mod,memobj$memlist$pmeans,memobj$memlist$pvars ) , 1, function(x) { rnorm(x[1],x[2],x[3]) } ))
  
  # reset the random number seed
  .GlobalEnv$.Random.seed <- oldseed
  
  # calculate the HPD interval
  hpdint <- boa.hpd(post.dist,alpha=1-cred.lev)
  upr <- hpdint[2]
  lwr <- hpdint[1]
  return(c(memobj$postmean,memobj$postvar,lwr,upr,memobj$ess)) 
  
}


# extract posterior mean, posterior variacne, and effective supplemental sample size from a multiple q iMEM object
memextr_multq <- function(multqobj,cred.lev=0.95) 
{ 
  # since this function uses the random number generator we must grab the Randomseed, and then reset it afterwards
  oldseed <- .GlobalEnv$.Random.seed
  
  # this calculates all 95% credible intervals for each iMEM object
  ints <- sapply( multqobj$modlist, 
                  function(memobj)
                  {
                    # sample the multinomial distribution for the normal mixture (out of 10000000 draws)
                    nsamp.mod <- c(rmultinom(1,100000,memobj$memlist$postwts))
                    
                    # sample the normal mixture
                    post.dist <- unlist( apply( cbind(nsamp.mod,memobj$memlist$pmeans,memobj$memlist$pvars ) , 1, 
                                                function(x) { rnorm(x[1],x[2],x[3]) } 
                    )
                    )
                    
                    # calculate the HPD interval
                    hpdint <- boa.hpd(post.dist,alpha=1-cred.lev)
                    upr <- hpdint[2]
                    lwr <- hpdint[1]
                    
                    return( c(lwr,upr) )
                  } 
  )
  
  # reset the random number seed
  .GlobalEnv$.Random.seed <- oldseed
  
  # final interval is the union of
  lwr <- min(ints[1,],na.rm=T)
  upr <- max(ints[2,],na.rm=T)
  
  return(c(multqobj$postmean,multqobj$postvar,lwr,upr,NA)) 
  
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################



boa.hpd <- function(x, alpha){
  ###Calculate HPD Intervals given vector of values and alpha
  n <- length(x)
  m <- max(1, ceiling(alpha * n))
  y <- sort(x)
  a <- y[1:m]
  b <- y[(n - m + 1):n]
  i <- order(b - a)[1]
  structure(c(a[i], b[i]), names = c("Lower Bound", "Upper Bound"))
}



####################### FUNCTIONS USED IN PRIMARY SIMULATION #################################

sim_oneiter <- function( mu, sigma, n, xbars, fnums=2:7,samp.sig=F)
{
  #######################################################################
  # function to run one simulation iteration for one supplementary source scenario, and one value of mu
  # simulates primary source data, fits simple mean, marginal iMEMs, and full MEM, and extracts the quantities desired:
  #  - effective supllemental sample size
  #  - posterior mean
  #  - posterior variacnce
  #  - 95% posterior credible intervals
  #
  #
  #######################################################################
  
  # final numbers of sources in marginal imem
  # fnums <- 2:7
  
  # simulate the primary source data from a normal distribution, calculate sample mean and sd
  dat <- rnorm(n=n,mean=mu,sd=sigma)
  xbarp <- mean(dat)
  if(samp.sig) { sdp <- sd(dat) } else { sdp <- sigma } # assume known or use sample sd
  np <- length(dat)
  
  # get data ready for input into MEM functions
  prim <- c(xbarp,sdp,np)
  suppmeans <- xbars
  suppsds <- rep(sigma,length(xbars))
  suppNs <- rep(n,length(xbars))
  
  # fit the models
  jmem <- mem_calc(prim=prim,means=suppmeans,sds=suppsds,Ns=suppNs,prior='pi_e')
  imems <- lapply(fnums,function(x){ imem_marg(prim=prim,means=suppmeans,sds=suppsds,Ns=suppNs,prior='pi_e',final_grpsize=x)} )
  names(imems) <- paste0("fsrc",fnums)
  
  # extract the relevant information
  nout <- c(xbarp,sdp/sqrt(np),xbarp-1.96*(sdp/sqrt(np)),xbarp+1.96*(sdp/sqrt(np)),0)   # simple mean
  jout <- memextr(jmem)             # fully joint model
  iout <- t(sapply(imems,memextr)) # marginal selection models
  out1 <- rbind(nout,iout,jout)
  out2 <- cbind(mu,c(0,fnums,-1),out1)
  
  return(out2)
}





# function to fit "nsim" runs of the simulation study for one value of mu and one supp source scenario
sim_onemu <- function(mu, sigma, n, xbars,fnums,samp.sig,nsim) { return ( replicate(nsim,sim_oneiter(mu,sigma,n,xbars,fnums,samp.sig),simplify="array") )}


# function to run simulation for all values of mu for "nsim" iterations each for one supp source scenario
sim_allmu <- function(mus,sigma,n,xbars,fnums,samp.sig,nsim) 
{ 
  out <- sapply(mus,sim_onemu,sigma,n,xbars,fnums,samp.sig,nsim,simplify="array" )
  dimnames(out) <- list( meth=c("simpmean",paste0("imem",fnums),"mem"),
                         quant=c("mu","model","postmean","postvar","cilo","ciup","esss"),
                         sim=paste0("sim",1:nsim),
                         truemean=paste0("mu",mus) 
  )
  
  return( out )
}


###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################

#### Simulation functions for fitting random effects models ######


# need to add functionality to draw comparisons between a random effects model and the MEMs.
# the problem is that we need a full dataset to run random effects models, not just a sample
# means and variances.  SO, I'll need to generate supplemental source datasets that satisfy the 
# sample mean and variance conditions for each scenario.  Then for each iteration of the simulation,
# I'll combine the supplemental datasets with the simulated primary dataset, and run the random effects
# model. 


re_oneiter <- function(mu,sigma,n,suppdat)
{
  primdat <- data.frame( id="prim",outcm=rnorm(n=n,mean=mu,sd=sigma) )
  primdat$id <- as.character(primdat$id)
  dat <- rbind(primdat,suppdat)
  
  # since predictInterval uses the random number generator we must grab the Randomseed, and then reset it afterwards
  oldseed <- .GlobalEnv$.Random.seed
  
  # grab the predicted value for the primary source
  reffmod <- lmer( outcm ~ 1 + (1|id) ,data=dat,na.action="na.exclude")
  ndat <- data.frame(id="prim")
  est <- predictInterval(reffmod,ndat,n.sims=100,type="linear.prediction") 
  
  # reset the random number seed
  .GlobalEnv$.Random.seed <- oldseed
  
  # format and return
  out1 <- rbind( c(est,NA),c(NA,NA,NA,NA) )
  out2 <- cbind(mu,c(-2,-3),out1)
  
  return(out2)
}




# function to fit "nsim" runs of the RE part of simulation study for one value of mu and one supp source scenario
simre_onemu <- function(mu, sigma, n, suppdat,nsim) { return ( replicate(nsim,re_oneiter(mu,sigma,n,suppdat),simplify="array") )}



# function to run simulation for all values of mu for "nsim" iterations each for one supp source scenario
simre_allmu <- function(mus,sigma,n,suppdat,nsim) 
{ 
  out <- sapply(mus,simre_onemu,sigma,n,suppdat,nsim,simplify="array" )
  dimnames(out) <- list( meth=c("repred","repop"),
                         quant=c("mu","model","postmean","cilo","ciup","esss"),
                         sim=paste0("sim",1:nsim),
                         truemean=paste0("mu",mus) 
  )
  
  return( out )
}





###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################

#### Simulation functions for fitting general iMEMs for various values of parameters ######


simGIMEM_oneiter <- function(mu,sigma,n,xbars,us=2:4,rs=c(0.5,0.75),qs=4,samp.sig=F)
{
  ######################################################################################################
  # function to run one simulation iteration for one supplementary source scneario, one value of mu
  # simulates primary source data, fits several general iMEMs with all combinations of specified parameters
  # and extracts 
  #  - effective supllemental sample size
  #  - posterior mean
  #  - posterior variacnce
  #
  # Function inputs:
  #   - mu: true mean of primary source
  #   - sigma: true standard deviation of primary source (as well as sample standard deviation of supp sources)
  #   - n: sample size of all sources
  #   - xbars: sample means of the supplementary sources
  #   - us: general iMEM group sizes 
  #   - rs: general iMEM acceptance ratios
  #   - qs: general iMEM final group size
  ################################################################################################
  
  # mu <- 2
  # sigma <- 4
  # n <- 20
  # xbars <- xbars6
  # us=2:4;rs=c(0.5,0.75);qs=4
  
  # simulate the primary source data from a normal distribution, calculate sample mean and sd
  dat <- rnorm(n=n,mean=mu,sd=sigma)
  xbarp <- mean(dat)
  if(samp.sig) { sdp <- sd(dat) } else { sdp <- sigma } # assume known or use sample sd
  np <- length(dat)
  
  # get data ready for input into MEM functions
  prim <- c(xbarp,sdp,np)
  suppmeans <- xbars
  suppsds <- rep(sigma,length(xbars))
  suppNs <- rep(n,length(xbars))
  suppids <- paste0("id",1:length(suppsds))
  
  
  # generate all combinations of general iMEM parameters to be considered
  comb <- expand.grid(us,rs,qs);names(comb) <- c("us","rs","qs")
  
  # fit all general imem models
  gimems <- apply(comb,1,function(x){ imem_gen(prim=prim,means=suppmeans,sds=suppsds,Ns=suppNs,
                                               suppids=suppids,order.start=1:length(suppmeans),prior='pi_e',
                                               grpsize=x[1],acc_prop = x[2],final_grpsize = x[3])  }  )
  
  
  out <- t(sapply(gimems,memextr)) # extract posterior mean, variance, and esss
  out2 <- cbind(mu,-4:(-4-nrow(comb)+1),out)
  return(out2)
}


# function to fit "nsim" runs of the simulation study for one value of mu and one supp source scenario
simGIMEM_onemu <- function(mu,sigma,n,xbars,us,rs,qs,samp.sig,nsim) 
{ return ( replicate(nsim,simGIMEM_oneiter(mu,sigma,n,xbars,us,rs,qs,samp.sig),simplify="array") ) }


# function to run simulation for all values of mu for "nsim" iterations each for one supp source scenario
simGIMEM_allmu <- function(mus,sigma,n,xbars,us,rs,qs,samp.sig,nsim) 
{
  
  comb <- expand.grid(us,rs,qs);names(comb) <- c("us","rs","qs")
  modnms <- paste0("gimem_u",comb$us,"_r",comb$rs,"_q",comb$qs)
  
  out <- sapply(mus,simGIMEM_onemu,sigma,n,xbars,us,rs,qs,samp.sig,nsim,simplify="array" )
  dimnames(out) <- list( meth=modnms,
                         quant=c("mu","model","postmean","postvar","cilo","ciup","esss"),
                         sim=paste0("sim",1:nsim),
                         truemean=paste0("mu",mus) 
  )
  
  return( out )
}





###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################

#### Simulation functions for fitting general iMEMs for various values of parameters ######


simGIMEM_rand_oneiter <- function(mu,sigma,n,xbars,us=2:4,rs=c(0.5,0.75),qs=4,samp.sig=F)
{
  ######################################################################################################
  # function to run one simulation iteration for one supplementary source scneario, one value of mu
  # simulates primary source data, fits several general random iMEMs with all combinations of specified parameters
  # and extracts 
  #  - effective supllemental sample size
  #  - posterior mean
  #  - posterior variacnce
  #
  # Function inputs:
  #   - mu: true mean of primary source
  #   - sigma: true standard deviation of primary source (as well as sample standard deviation of supp sources)
  #   - n: sample size of all sources
  #   - xbars: sample means of the supplementary sources
  #   - us: general iMEM group sizes 
  #   - rs: general iMEM acceptance ratios
  #   - qs: general iMEM final group size
  ################################################################################################
  
  # mu <- 2
  # sigma <- 4
  # n <- 20
  # xbars <- xbars6
  # us=2:4;rs=c(0.5,0.75);qs=4
  
  # simulate the primary source data from a normal distribution, calculate sample mean and sd
  dat <- rnorm(n=n,mean=mu,sd=sigma)
  xbarp <- mean(dat)
  if(samp.sig) { sdp <- sd(dat) } else { sdp <- sigma } # assume known or use sample sd
  np <- length(dat)
  
  # get data ready for input into MEM functions
  prim <- c(xbarp,sdp,np)
  suppmeans <- xbars
  suppsds <- rep(sigma,length(xbars))
  suppNs <- rep(n,length(xbars))
  suppids <- paste0("id",1:length(suppsds))
  
  
  # generate all combinations of general iMEM parameters to be considered
  comb <- expand.grid(us,rs,qs);names(comb) <- c("us","rs","qs")
  
  # since random imem uses the random number generator we must grab the Randomseed, and then reset it afterwards
  oldseed <- .GlobalEnv$.Random.seed
  
  # fit all general imem models
  gimems <- apply(comb,1,function(x){ imem_gen.rand(prim=prim,means=suppmeans,sds=suppsds,Ns=suppNs,
                                                    suppids=suppids,prior='pi_e',
                                                    grpsize=x[1],acc_prop = x[2],final_grpsize = x[3])  }  )
  
  # reset the random number seed
  .GlobalEnv$.Random.seed <- oldseed
  
  out <- t(sapply(gimems,memextr)) # extract posterior mean, variance, and esss
  out2 <- cbind(mu,-10:(-10-nrow(comb)+1),out)
  return(out2)
}


# function to fit "nsim" runs of the simulation study for one value of mu and one supp source scenario
simGIMEM_rand_onemu <- function(mu,sigma,n,xbars,us,rs,qs,samp.sig,nsim) 
{ return ( replicate(nsim,simGIMEM_rand_oneiter(mu,sigma,n,xbars,us,rs,qs,samp.sig),simplify="array") ) }


# function to run simulation for all values of mu for "nsim" iterations each for one supp source scenario
simGIMEM_rand_allmu <- function(mus,sigma,n,xbars,us,rs,qs,samp.sig,nsim) 
{
  
  comb <- expand.grid(us,rs,qs);names(comb) <- c("us","rs","qs")
  modnms <- paste0("gimemrand_u",comb$us,"_r",comb$rs,"_q",comb$qs)
  
  out <- sapply(mus,simGIMEM_rand_onemu,sigma,n,xbars,us,rs,qs,samp.sig,nsim,simplify="array" )
  dimnames(out) <- list( meth=modnms,
                         quant=c("mu","model","postmean","postvar","cilo","ciup","esss"),
                         sim=paste0("sim",1:nsim),
                         truemean=paste0("mu",mus) 
  )
  
  return( out )
}


###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################

#### Simulation functions for fitting marginal iMEMs only ######


sim_MIMEM_oneiter <- function( mu, sigma, n, xbars, fnums=2:7,samp.sig=F)
{
  #######################################################################
  # function to run one simulation iteration for one supplementary source scenario, and one value of mu
  # simulates primary source data, fits simple mean, marginal iMEMs, and full MEM, and extracts the quantities desired:
  #  - effective supllemental sample size
  #  - posterior mean
  #  - posterior variacnce
  #
  #
  #######################################################################
  
  # final numbers of sources in marginal imem
  # fnums <- 2:7
  
  # simulate the primary source data from a normal distribution, calculate sample mean and sd
  dat <- rnorm(n=n,mean=mu,sd=sigma)
  xbarp <- mean(dat)
  if(samp.sig) { sdp <- sd(dat) } else { sdp <- sigma } # assume known or use sample sd
  np <- length(dat)
  
  # get data ready for input into MEM functions
  prim <- c(xbarp,sdp,np)
  suppmeans <- xbars
  suppsds <- rep(sigma,length(xbars))
  suppNs <- rep(n,length(xbars))
  
  # fit the models
  imems <- lapply(fnums,function(x){ imem_marg(prim=prim,means=suppmeans,sds=suppsds,Ns=suppNs,prior='pi_e',final_grpsize=x)} )
  names(imems) <- paste0("fsrc",fnums)
  
  # extract the relevant information
  iout <- t(sapply(imems,memextr)) # marginal selection models
  out1 <- rbind(iout)
  out2 <- cbind(mu,fnums,out1)
  
  return(out2)
}




# function to fit "nsim" runs of the simulation study for one value of mu and one supp source scenario
sim_MIMEM_onemu <- function(mu, sigma, n, xbars,fnums,samp.sig,nsim) { return ( replicate(nsim,sim_MIMEM_oneiter(mu,sigma,n,xbars,fnums,samp.sig),simplify="array") )}


# function to run simulation for all values of mu for "nsim" iterations each for one supp source scenario
sim_MIMEM_allmu <- function(mus,sigma,n,xbars,fnums,samp.sig,nsim) 
{ 
  out <- sapply(mus,sim_MIMEM_onemu,sigma,n,xbars,fnums,samp.sig,nsim,simplify="array" )
  dimnames(out) <- list( meth=c(paste0("imem",fnums)),
                         quant=c("mu","model","postmean","postvar","cilo","ciup","esss"),
                         sim=paste0("sim",1:nsim),
                         truemean=paste0("mu",mus) 
  )
  
  return( out )
}


###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################

#### Simulation functions for fitting marginal iMEMs with q determined by thresholding similarity scores ######


sim_MIMEMthresh_oneiter <- function( mu, sigma, n, xbars, fnums=c(0.05,0.1,0.2,0.3,0.4,0.5),samp.sig=F)
{
  #######################################################################
  # function to run one simulation iteration for one supplementary source scenario, and one value of mu
  # simulates primary source data, fits simple mean, marginal iMEMs, and full MEM, and extracts the quantities desired:
  #  - effective supllemental sample size
  #  - posterior mean
  #  - posterior variacnce
  #
  #
  #######################################################################
  
  # final numbers of sources in marginal imem
  # fnums <- 2:7
  
  # simulate the primary source data from a normal distribution, calculate sample mean and sd
  dat <- rnorm(n=n,mean=mu,sd=sigma)
  xbarp <- mean(dat)
  if(samp.sig) { sdp <- sd(dat) } else { sdp <- sigma } # assume known or use sample sd
  np <- length(dat)
  
  # get data ready for input into MEM functions
  prim <- c(xbarp,sdp,np)
  suppmeans <- xbars
  suppsds <- rep(sigma,length(xbars))
  suppNs <- rep(n,length(xbars))
  
  # fit the models
  imems <- lapply(fnums,function(x){ imem_marg(prim=prim,means=suppmeans,sds=suppsds,Ns=suppNs,prior='pi_e',final_grpsize="thresh",qthresh=x)} )
  names(imems) <- paste0("thresh",fnums)
  
  # extract the relevant information
  iout <- t(sapply(imems,memextr)) # marginal selection models
  out1 <- rbind(iout)
  out2 <- cbind(mu,fnums,out1)
  
  return(out2)
}




# function to fit "nsim" runs of the simulation study for one value of mu and one supp source scenario
sim_MIMEMthresh_onemu <- function(mu, sigma, n, xbars,fnums,samp.sig,nsim) { return ( replicate(nsim,sim_MIMEMthresh_oneiter(mu,sigma,n,xbars,fnums,samp.sig),simplify="array") )}


# function to run simulation for all values of mu for "nsim" iterations each for one supp source scenario
sim_MIMEMthresh_allmu <- function(mus,sigma,n,xbars,fnums,samp.sig,nsim) 
{ 
  out <- sapply(mus,sim_MIMEMthresh_onemu,sigma,n,xbars,fnums,samp.sig,nsim,simplify="array" )
  dimnames(out) <- list( meth=c(paste0("imemthresh",fnums)),
                         quant=c("mu","model","postmean","postvar","cilo","ciup","esss"),
                         sim=paste0("sim",1:nsim),
                         truemean=paste0("mu",mus) 
  )
  
  return( out )
}




###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################

#### Simulation functions for fitting general iMEMs for various values of parameters ######


simMIMEM_multq_oneiter <- function(mu,sigma,n,xbars,qmaxs=c(7,9),qmins=c(2,3,4),samp.sig=F)
{
  ######################################################################################################
  # function to run one simulation iteration for one supplementary source scneario, one value of mu
  # simulates primary source data, fits several general random iMEMs with all combinations of specified parameters
  # and extracts 
  #  - effective supllemental sample size
  #  - posterior mean
  #  - posterior variacnce
  #
  # Function inputs:
  #   - mu: true mean of primary source
  #   - sigma: true standard deviation of primary source (as well as sample standard deviation of supp sources)
  #   - n: sample size of all sources
  #   - xbars: sample means of the supplementary sources
  #   - us: general iMEM group sizes 
  #   - rs: general iMEM acceptance ratios
  #   - qs: general iMEM final group size
  ################################################################################################
  
  # mu <- 2
  # sigma <- 4
  # n <- 20
  # xbars <- xbars6
  # us=2:4;rs=c(0.5,0.75);qs=4
  
  # simulate the primary source data from a normal distribution, calculate sample mean and sd
  dat <- rnorm(n=n,mean=mu,sd=sigma)
  xbarp <- mean(dat)
  if(samp.sig) { sdp <- sd(dat) } else { sdp <- sigma } # assume known or use sample sd
  np <- length(dat)
  
  # get data ready for input into MEM functions
  prim <- c(xbarp,sdp,np)
  suppmeans <- xbars
  suppsds <- rep(sigma,length(xbars))
  suppNs <- rep(n,length(xbars))
  suppids <- paste0("id",1:length(suppsds))
  
  
  # generate all combinations of general iMEM parameters to be considered
  comb <- expand.grid(qmaxs,qmins);names(comb) <- c("qmaxs","qmins")

  # fit all general imem models
  imems <- apply(comb,1,function(x){ imem_marg_multq(prim=prim,means=suppmeans,sds=suppsds,Ns=suppNs,
                                                    prior='pi_e',qmax=x[1],qmin=x[2])  }  )
  
  out <- t(sapply(imems,memextr_multq)) # extract posterior mean, variance, and esss
  out2 <- cbind(mu,-20:(-20-nrow(comb)+1),out)
  return(out2)
}


# function to fit "nsim" runs of the simulation study for one value of mu and one supp source scenario
simMIMEM_multq_onemu <- function(mu,sigma,n,xbars,qmaxs,qmins,samp.sig,nsim) 
{ return ( replicate(nsim,simMIMEM_multq_oneiter(mu,sigma,n,xbars,qmaxs,qmins,samp.sig),simplify="array") ) }


# function to run simulation for all values of mu for "nsim" iterations each for one supp source scenario
simMIMEM_multq_allmu <- function(mus,sigma,n,xbars,qmaxs,qmins,samp.sig,nsim) 
{
  
  comb <- expand.grid(qmaxs,qmins);names(comb) <- c("qmaxs","qmins")
  modnms <- paste0("imemmultq_",comb$qmins,"_",comb$qmaxs)
  
  out <- sapply(mus,simMIMEM_multq_onemu,sigma,n,xbars,qmaxs,qmins,samp.sig,nsim,simplify="array" )
  dimnames(out) <- list( meth=modnms,
                         quant=c("mu","model","postmean","postvar","cilo","ciup","esss"),
                         sim=paste0("sim",1:nsim),
                         truemean=paste0("mu",mus) 
  )
  
  return( out )
}

