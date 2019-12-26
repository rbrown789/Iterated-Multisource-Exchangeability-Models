


# root <- "C:/Users/roland/Documents/iMEM Manuscript Data and Code/"
root <- "C:/Users/rbrow_000/Documents/UofM Documents/Dissertation/Iterated-Multisource-Exchangeability-Models/"
root <- "C:/Users/rbrow/OneDrive/Documents/iMEM Manuscript Data and Code/"
code <- "Code/"
data <- "Data/"

# Load functions and iMEM data
source(paste0(root,code,"iMEM_functions.R"))
load(paste0(root,data,"imem_realdata.Rdata"))

library(ggplot2)
library(lme4)
library(merTools)
library(gridExtra)

umn <- imem_realdata

################################################################################


### Function to grab all trips/activities of a certain type for a certain i
typeid <- function(dat,id,type,subtype,varnms=names(dat))
{
  return(dat[dat$surveycode %in% id & dat$type %in% type & dat$primary_mode %in% subtype, varnms])
}


# function that given subtype, primary ID, and dataset, grabs all the supplementary indiviudals with
# "enough" data, maybe > 7 observations, and runs an iMEM
imem_dayn <- function(outname,type,subtype,primid,nfinsrc,qthresh=NULL,dat,sel="marg")
{
  # outname <- "happy"
  # type <- "TRIP"
  # subtype <- "BUS"
  # primid <- "2018"
  # primid <- "3096"
  # nfinsrc <- 12
  # dat <- umn
  
  # create primary source information vector
  prim_dat <- typeid(dat,primid,type,subtype,varnms=c("surveycode",outname))
  prim_out <- prim_dat[[outname]]
  prim_info <- c(mean(prim_out,na.rm=T),sd(prim_out,na.rm=T),sum(!is.na(prim_out))  )
  
  if(prim_info[2] ==0) { out <- list(postmean=prim_info[1],postvar=0,note="standard deviation 0 for primary source")
  } else {
  
  # find all the valid supplementary ids with >= 5 observations
  tmp <- typeid(dat,unique(dat$surveycode)[!unique(dat$surveycode) %in% primid],type,subtype,varnms=c("surveycode",outname))
  tmp <- tmp[!is.na(tmp[[outname]]),]
  tbl <- rev(sort(table(tmp$surveycode)))
  suppids <- names(tbl)[tbl >= 5]

  # create supplmentary information vectors
  suppmeans <- suppsds <- suppNs <- suppids2 <- c()
  suppdat_all <- NULL
  
  for(i in suppids)
  {
    supp_dat <- typeid(dat,i,type,subtype,varnms=c("surveycode",outname) )
    supp_out <- supp_dat[[outname]]
    
    tmp.sd <- sd(supp_out,na.rm=T)
    
    if(tmp.sd > 0) {
    suppmeans <- c(suppmeans,mean(supp_out,na.rm=T) )
    suppsds <- c(suppsds,tmp.sd)
    suppNs <- c(suppNs,sum(!is.na(supp_out)))
    suppids2 <- c(suppids2,i)
    suppdat_all <- rbind(suppdat_all,supp_dat)
    }
    
    
  }
  
  # apply MEM function to the primary and supplmentary information
  
  
  if(sel=="marg") { out <- imem_marg( prim=prim_info,means=suppmeans,sds=suppsds,Ns=suppNs,prior="pi_e",final_grpsize = nfinsrc,qthresh=qthresh)
  }  else if(sel=="gen"){
    out <- imem_gen.rand( prim=prim_info,means=suppmeans,sds=suppsds,Ns=suppNs,suppids=suppids2,prior="pi_e",grpsize=5,acc_prop=0.5,final_grpsize = nfinsrc)
  }
  
  }
  
  return(out)
}


remod_dayn <- function(outname,type,subtype,primid,dat)
{
  # outname <- "happy"
  # type <- "TRIP"
  # subtype <- "BUS"
  # primid <- "2018"
  # primid <- "3096"
  # nfinsrc <- 12
  # dat <- umn
  
  # create primary source information vector
  prim_dat <- typeid(dat,primid,type,subtype,varnms=c("surveycode",outname))
  
  
    
    # find all the valid supplementary ids with > 7 observations
    tmp <- typeid(dat,unique(dat$surveycode)[!unique(dat$surveycode) %in% primid],type,subtype,varnms=c("surveycode",outname))
    tmp <- tmp[!is.na(tmp[[outname]]),]
    tbl <- rev(sort(table(tmp$surveycode)))
    suppids <- names(tbl)[tbl >= 5]
    
    # create supplmentary information vectors
    suppdat_all <- NULL
    
    for(i in suppids)
    {
      supp_dat <- typeid(dat,i,type,subtype,varnms=c("surveycode",outname) )
      tmp.sd <- sd(supp_dat[[outname]],na.rm=T)
      if(tmp.sd > 0) { suppdat_all <- rbind(suppdat_all,supp_dat) }
    }
    
    # fit a random effects model
    redat <- rbind(prim_dat,suppdat_all)
    form <- as.formula(paste0(outname,"~ (1|surveycode)"))
    mod <- lmer(form,data=redat)
  

  
  
  return(mod)
}


# calculate posterior credible intervals from imem output
credint.imem <- function(imemobj,cred.lev=0.95)
{
 
  # sample the multinomial distribution for the normal mixture (out of 10000000 draws)
  nsamp.mod <- c(rmultinom(1,1000000,imemobj$memlist$postwts))
  
  # sample the normal mixture
  post.dist <- unlist( apply( cbind(nsamp.mod,imemobj$memlist$pmeans,imemobj$memlist$pvars ) , 1, function(x) { rnorm(x[1],x[2],x[3]) } ))
  
  # calculate the HPD interval
  hpdint <- boa.hpd(post.dist,alpha=1-cred.lev)
  pm <- imemobj$postmean
  pse <- sqrt(imemobj$postvar)
  upr <- hpdint[2]
  lwr <- hpdint[1]
  
  out <- c(pm,pse,lwr,upr)
  names(out) <- c("postmean","postsd","lwr","upr")
  
  return(out)
}



credint.smean <- function(imemobj,cred.lev=0.95)
{
  pm <- imemobj$memlist$pmeans[1]
  pse <- sqrt(imemobj$memlist$pvars[1])
  upr <- pm + pse*qnorm( 1-(1-cred.lev)/2 )
  lwr <- pm - pse*qnorm( 1-(1-cred.lev)/2 )
  
  out <- c(pm,pse,lwr,upr)
  names(out) <- c("smean","ssd","slwr","supr")
  
  return(out)
}


#################

# library(KScorrect)
# # 
# # fit the imem models
# emot <- "tired"; tmode <- "BUS"; id <- 3096
# imem12 <- imem_dayn(emot,type="TRIP",subtype=tmode,primid=id,nfinsrc=12,dat=umn)
# imem7 <-  imem_dayn(emot,type="TRIP",subtype=tmode,primid=id,nfinsrc=7,dat=umn)
# imem2 <-  imem_dayn(emot,type="TRIP",subtype=tmode,primid=id,nfinsrc=2,dat=umn)
# imemthresh <- imem_dayn(emot,type="TRIP",subtype=tmode,primid=id,nfinsrc="thresh",qthresh=0.2,dat=umn)
# imem12$srcinfo
# 
# x <- seq(1.8,3.2,length.out=200)
# d12 <- dmixnorm(x,imem12$memlist$pmeans,sqrt(imem12$memlist$pvars),imem12$memlist$postwts)
# d7 <- dmixnorm(x,imem7$memlist$pmeans,sqrt(imem7$memlist$pvars),imem7$memlist$postwts)
# d2 <- dmixnorm(x,imem2$memlist$pmeans,sqrt(imem2$memlist$pvars),imem2$memlist$postwts)
# 
# plot(x,d12,type="l",col="blue",lwd=2)
# lines(x,d7,col="red",lwd=2)
# lines(x,d2,lwd=2)
# abline(v=imem12$srcinfo[1,],lty=c(1,rep(2,ncol(imem12$srcinfo)-1)),
#        col=c("black",rep("grey60",ncol(imem12$srcinfo)-1)))

gendat.emotbymode <- function(primid,outcms)
{
  plotdatfin <- NULL
  for(i in outcms) 
  {
  

    # fit the imem models
    imem.car <- imem_dayn(i,type="TRIP",subtype="CAR",primid=primid,nfinsrc=12,dat=umn)
    imem.bus <- imem_dayn(i,type="TRIP",subtype="BUS",primid=primid,nfinsrc=12,dat=umn)
    imem.walk <- imem_dayn(i,type="TRIP",subtype="WALK",primid=primid,nfinsrc=12,dat=umn)
    
    # fit imem models with fewer supp sources
    imem.car7 <- imem_dayn(i,type="TRIP",subtype="CAR",primid=primid,nfinsrc=7,dat=umn)
    imem.bus7 <- imem_dayn(i,type="TRIP",subtype="BUS",primid=primid,nfinsrc=7,dat=umn)
    imem.walk7 <- imem_dayn(i,type="TRIP",subtype="WALK",primid=primid,nfinsrc=7,dat=umn)
    
    imem.car2 <- imem_dayn(i,type="TRIP",subtype="CAR",primid=primid,nfinsrc=2,dat=umn)
    imem.bus2 <- imem_dayn(i,type="TRIP",subtype="BUS",primid=primid,nfinsrc=2,dat=umn)
    imem.walk2 <- imem_dayn(i,type="TRIP",subtype="WALK",primid=primid,nfinsrc=2,dat=umn)
    
    # fit imem models with score-threshold determined q
    imem.car0.05 <- imem_dayn(i,type="TRIP",subtype="CAR",primid=primid,nfinsrc="thresh",qthresh=0.05,dat=umn)
    imem.bus0.05 <- imem_dayn(i,type="TRIP",subtype="BUS",primid=primid,nfinsrc="thresh",qthresh=0.05,dat=umn)
    imem.walk0.05 <- imem_dayn(i,type="TRIP",subtype="WALK",primid=primid,nfinsrc="thresh",qthresh=0.05,dat=umn)
    
    imem.car0.1 <- imem_dayn(i,type="TRIP",subtype="CAR",primid=primid,nfinsrc="thresh",qthresh=0.1,dat=umn)
    imem.bus0.1 <- imem_dayn(i,type="TRIP",subtype="BUS",primid=primid,nfinsrc="thresh",qthresh=0.1,dat=umn)
    imem.walk0.1 <- imem_dayn(i,type="TRIP",subtype="WALK",primid=primid,nfinsrc="thresh",qthresh=0.1,dat=umn)
    
    imem.car0.2 <- imem_dayn(i,type="TRIP",subtype="CAR",primid=primid,nfinsrc="thresh",qthresh=0.2,dat=umn)
    imem.bus0.2 <- imem_dayn(i,type="TRIP",subtype="BUS",primid=primid,nfinsrc="thresh",qthresh=0.2,dat=umn)
    imem.walk0.2 <- imem_dayn(i,type="TRIP",subtype="WALK",primid=primid,nfinsrc="thresh",qthresh=0.2,dat=umn)
    
    
    # fit the random effect models
    re.car <- remod_dayn(i,type="TRIP",subtype="CAR",primid=primid,dat=umn)
    re.bus <- remod_dayn(i,type="TRIP",subtype="BUS",primid=primid,dat=umn)
    re.walk <- remod_dayn(i,type="TRIP",subtype="WALK",primid=primid,dat=umn)
    
    lst <- list(imem.car,imem.bus,imem.walk)
    lst7 <- list(imem.car7,imem.bus7,imem.walk7)
    lst2 <- list(imem.car2,imem.bus2,imem.walk2)
    lst0.05 <- list(imem.car0.05,imem.bus0.05,imem.walk0.05)
    lst0.1 <- list(imem.car0.1,imem.bus0.1,imem.walk0.1)
    lst0.2 <- list(imem.car0.2,imem.bus0.2,imem.walk0.2)

    # generate plotting data with posterior means/intervals, simple means/intervals, esss
    plotdat <- data.frame( t(sapply(lst,credint.imem,cred.lev=0.95)))
    plotdat7 <- data.frame( t(sapply(lst7,credint.imem,cred.lev=0.95))); names(plotdat7) <- paste0(names(plotdat7),"7")
    plotdat2 <- data.frame( t(sapply(lst2,credint.imem,cred.lev=0.95))); names(plotdat2) <- paste0(names(plotdat2),"2")
    plotdat0.05 <- data.frame( t(sapply(lst0.05,credint.imem,cred.lev=0.95))); names(plotdat0.05) <- paste0(names(plotdat0.05),"0.05")
    plotdat0.1 <- data.frame( t(sapply(lst0.1,credint.imem,cred.lev=0.95))); names(plotdat0.1) <- paste0(names(plotdat0.1),"0.1")
    plotdat0.2 <- data.frame( t(sapply(lst0.2,credint.imem,cred.lev=0.95))); names(plotdat0.2) <- paste0(names(plotdat0.2),"0.2")
    
    smean <-  data.frame( t(sapply(lst,credint.smean,cred.lev=0.95)) )
    plotdat <- cbind(plotdat,plotdat7,plotdat2,plotdat0.05,plotdat0.1,plotdat0.2,smean)
    plotdat$esss <-  sapply(lst,function(x){x$ess})
    plotdat$esss7 <-  sapply(lst7,function(x){x$ess})
    plotdat$esss2 <-  sapply(lst2,function(x){x$ess})
    plotdat$esss0.05 <-  sapply(lst0.05,function(x){x$ess})
    plotdat$esss0.1 <-  sapply(lst0.1,function(x){x$ess})
    plotdat$esss0.2 <-  sapply(lst0.2,function(x){x$ess})
    plotdat$mode <- paste0( c("Car","Bus","Walk")," (N=", sapply(lst,function(x) { round(x$srcinfo[3,1],0) } ),")" )
    plotdat$outcm <- i
    
    # now for random effects models (use predictInterval, from the merTools package to do interval estimation)
    ndat <- data.frame(surveycode=primid)
    repdat <- rbind( predictInterval(re.car,ndat,n.sims=1000,type="linear.prediction"),
                     predictInterval(re.bus,ndat,n.sims=1000,type="linear.prediction"),
                     predictInterval(re.walk,ndat,n.sims=1000,type="linear.prediction") )
                     
    names(repdat) <- c("remean","relwr","reupr")
    plotdat <- cbind(plotdat,repdat)
    
    plotdatfin <- rbind(plotdatfin,plotdat)
  }
  
  plotdatfin$mode <- as.factor(plotdatfin$mode)
  plotdatfin$outcm <- as.factor(plotdatfin$outcm)
  plotdatfin$`Trip Mode` <- plotdatfin$mode
  plotdatfin$pctdiffmean <- abs(plotdat$postmean-plotdat$smean)/plotdat$smean
  plotdatfin$pctdecsd <- (plotdatfin$ssd-plotdatfin$postsd)/plotdatfin$ssd
  return(plotdatfin)
}


################################################################################################
################################################################################################
################################################################################################

#### Plotting functions for the primary real data figure #####

# original
plot.emotbymode <- function(plotdatfin,primid)
{
  shift <- 0.07
  dodge <- position_dodge(width = 0.9)
  limits <- aes(ymax = plotdatfin$upr, ymin = plotdatfin$lwr,x=as.numeric(outcm))
  limsmean <- aes(ymax=plotdatfin$supr,ymin=plotdatfin$slwr,x = as.numeric(outcm)+shift)
  limre <- aes(ymax=plotdatfin$reupr,ymin=plotdatfin$relwr,x=as.numeric(outcm)-shift)
  
  p <- ggplot(data = plotdatfin,aes(x = outcm, y = postmean,fill=`Trip Mode`)) 
  
  p <- p + geom_bar( stat = "identity", position = dodge) + 
    geom_point( aes(x = as.numeric(outcm)+shift, y = smean),stat = "identity", position = dodge,size=rel(1),shape=2) +
    geom_point( aes(x = as.numeric(outcm)-shift, y = remean),stat = "identity", position = dodge,size=rel(1),shape=3) +
    scale_y_continuous(name = "Intensity") + 
    coord_cartesian(ylim=c(0,6.2)) + 
    scale_fill_manual(values=c("grey70","steelblue","orange")) +
    
    geom_errorbar(limits, position = dodge, width = 0.15,size=0.7) +
    geom_errorbar(limsmean, position=dodge, width = 0.15,size=0.4,col="black",linetype=2) +
    geom_errorbar(limre, position=dodge, width = 0.15,size=0.4,col="black",linetype=3) +
    
    geom_text(aes(x=outcm,y=0.5,label=round(esss,0)),stat="identity",position=dodge,size=4.5) +
    
    
    
    
    guides(fill = guide_legend(override.aes = list(shape = NA) ) )     +
    
    theme(axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_text(size=rel(1.9)),
          axis.title.y=element_text(size = rel(1.9)),
          axis.text.x = element_text(size = rel(1.8)),
          axis.title.x= element_blank(),
          legend.position="bottom" ,
          legend.text=element_text(size=rel(1.3)),
          legend.title=element_text(size=rel(1.6)),
          plot.title =element_text(size=20, face='bold')
          
    )  +
    
    
    ggtitle(paste0("Emotion Intensity by Trip Mode: ID ",primid))
  
  print(p)
}




# increase size of points and font
plot.emotbymode2 <- function(plotdatfin,primid)
{
  shift <- 0.07
  dodge <- position_dodge(width = 0.9)
  limits <- aes(ymax = plotdatfin$upr, ymin = plotdatfin$lwr,x=as.numeric(outcm))
  limsmean <- aes(ymax=plotdatfin$supr,ymin=plotdatfin$slwr,x = as.numeric(outcm)+shift)
  limre <- aes(ymax=plotdatfin$reupr,ymin=plotdatfin$relwr,x=as.numeric(outcm)-shift)
  
  p <- ggplot(data = plotdatfin,aes(x = outcm, y = postmean,fill=`Trip Mode`)) 
  
  p <- p + geom_bar( stat = "identity", position = dodge) + 
    geom_point( aes(x = as.numeric(outcm)+shift, y = smean),stat = "identity", position = dodge,size=rel(3.5),shape=2) +
    geom_point( aes(x = as.numeric(outcm)-shift, y = remean),stat = "identity", position = dodge,size=rel(3.5),shape=4) +
    scale_y_continuous(name = "Intensity") + 
    coord_cartesian(ylim=c(0,6.2)) + 
    scale_fill_manual(values=c("grey70","steelblue","orange")) +
    
    geom_errorbar(limits, position = dodge, width = 0.15,size=0.7) +
    geom_errorbar(limsmean, position=dodge, width = 0.15,size=0.4,col="black",linetype=2) +
    geom_errorbar(limre, position=dodge, width = 0.15,size=0.4,col="black",linetype=3) +
    
    geom_text(aes(x=outcm,y=0.5,label=round(esss,0)),stat="identity",position=dodge,size=6) +
    
    
    
    
    guides(fill = guide_legend(override.aes = list(shape = NA) ) )     +
    
    theme(axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_text(size=rel(1.9)),
          axis.title.y=element_text(size = rel(1.9)),
          axis.text.x = element_text(size = rel(2.5)),
          axis.title.x= element_blank(),
          legend.position="bottom" ,
          legend.text=element_text(size=rel(1.6)),
          legend.title=element_text(size=rel(1.9)),
          plot.title =element_text(size=20, face='bold')
          
    )  +
    
    
    ggtitle(paste0("Emotion Intensity by Trip Mode: ID ",primid))
  
  print(p)
}



# black and white
plot.emotbymode_bw <- function(plotdatfinlst,primids)
{
  # plotdatfinlst <- pdatlst
  # primids <- primids
  # 
  shift <- 0.07
  dodge <- position_dodge(width = 0.9)
  
  ##########################################################################
  # first plot
  
  plotdatfin <- plotdatfinlst[[1]]
  primid <- primids[1]
  
  limits <- aes(ymax = plotdatfin$upr, ymin = plotdatfin$lwr,x=as.numeric(outcm))
  limsmean <- aes(ymax=plotdatfin$supr,ymin=plotdatfin$slwr,x = as.numeric(outcm)+shift)
  limre <- aes(ymax=plotdatfin$reupr,ymin=plotdatfin$relwr,x=as.numeric(outcm)-shift)
  
  p1 <- ggplot(data = plotdatfin,aes(x = outcm, y = postmean,fill=`Trip Mode`)) 
  
  p1 <- p1 + geom_bar( stat = "identity", position = dodge) + 
    geom_point( aes(x = as.numeric(outcm)+shift, y = smean),stat = "identity", position = dodge,size=rel(3.5),shape=2) +
    geom_point( aes(x = as.numeric(outcm)-shift, y = remean),stat = "identity", position = dodge,size=rel(3.5),shape=4) +
    scale_y_continuous(name = "Intensity") + 
    coord_cartesian(ylim=c(0,6.2)) + 
    scale_fill_manual(values=c("grey40","grey60","grey80")) +
    
    geom_errorbar(limits, position = dodge, width = 0.15,size=0.7) +
    geom_errorbar(limsmean, position=dodge, width = 0.15,size=0.4,col="black",linetype=2) +
    geom_errorbar(limre, position=dodge, width = 0.15,size=0.4,col="black",linetype=3) +
    
    geom_text(aes(x=outcm,y=0.5,label=round(esss,0)),stat="identity",position=dodge,size=6) +
    
    guides(fill = guide_legend(override.aes = list(shape = NA) ) )     +
    
    theme_bw() +
    
    theme(axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_text(size=rel(1.9)),
          axis.title.y=element_text(size = rel(1.9)),
          axis.text.x = element_text(size = rel(2.5)),
          axis.title.x= element_blank(),
          legend.position="bottom" ,
          legend.text=element_text(size=rel(1.6)),
          legend.title=element_text(size=rel(1.9)),
          plot.title =element_text(size=20, face='bold')
          
    )  +
    
    ggtitle(paste0("Emotion Intensity by Trip Mode: ID ",primid)) 
  
  
  ##########################################################################
  # second plot
  
  plotdatfin2 <- plotdatfinlst[[2]]
  primid2 <- primids[2]
  
  limits_2 <- aes(ymax = plotdatfin2$upr, ymin = plotdatfin2$lwr,x=as.numeric(outcm))
  limsmean_2 <- aes(ymax=plotdatfin2$supr,ymin=plotdatfin2$slwr,x = as.numeric(outcm)+shift)
  limre_2 <- aes(ymax=plotdatfin2$reupr,ymin=plotdatfin2$relwr,x=as.numeric(outcm)-shift)
  
  p2 <- ggplot(data = plotdatfin2,aes(x = outcm, y = postmean,fill=`Trip Mode`)) 
  
  p2 <- p2 + geom_bar( stat = "identity", position = dodge) + 
    geom_point( aes(x = as.numeric(outcm)+shift, y = smean),stat = "identity", position = dodge,size=rel(3.5),shape=2) +
    geom_point( aes(x = as.numeric(outcm)-shift, y = remean),stat = "identity", position = dodge,size=rel(3.5),shape=4) +
    scale_y_continuous(name = "Intensity") + 
    coord_cartesian(ylim=c(0,6.2)) + 
    scale_fill_manual(values=c("grey40","grey60","grey80")) +
    
    geom_errorbar(limits_2, position = dodge, width = 0.15,size=0.7) +
    geom_errorbar(limsmean_2, position=dodge, width = 0.15,size=0.4,col="black",linetype=2) +
    geom_errorbar(limre_2, position=dodge, width = 0.15,size=0.4,col="black",linetype=3) +
    
    geom_text(aes(x=outcm,y=0.5,label=round(esss,0)),stat="identity",position=dodge,size=6) +
    
    guides(fill = guide_legend(override.aes = list(shape = NA) ) )     +
    
    theme_bw() +
    
    theme(axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_text(size=rel(1.9)),
          axis.title.y=element_text(size = rel(1.9)),
          axis.text.x = element_text(size = rel(2.5)),
          axis.title.x= element_blank(),
          legend.position="bottom" ,
          legend.text=element_text(size=rel(1.6)),
          legend.title=element_text(size=rel(1.9)),
          plot.title =element_text(size=20, face='bold')
          
    )  +
    
    ggtitle(paste0("Emotion Intensity by Trip Mode: ID ",primid2)) 
  
  ########## mock plot for legend extraction #####################
  line_types <- c("Random Effects"=3,"Simple Mean"=2,"iMEM (q = 12)"=1)
  sizes <- c("Random Effects"=1,"Simple Mean"=1,"iMEM (q = 12)"=1)
  shps <- c("Random Effects"=4,"Simple Mean"=2,"iMEM (q = 12)"=1)
  legdat <- data.frame(y = 1:9,x = 9:1,grp=rep(c("Random Effects", "Simple Mean","iMEM (q = 12)"),3))
  legdat$grp <- factor(legdat$grp,levels(legdat$grp)[c(2,1,3)])
  
  legplot <- ggplot(data=legdat, aes(x=x,y=y,linetype=grp,size=grp,shape=grp)) + geom_line() + geom_point() +
    scale_linetype_manual(name="",values=line_types) +
    scale_size_manual(name="",values=sizes) +
    scale_shape_manual(name="",values=shps) +
    guides(shape=guide_legend(override.aes=list(size=3))) + theme_bw() + 
    theme(legend.text=element_text(size=rel(1.3)))
  
  # customize the legend (separte sizes for points and lines) by editing the low level grid objects. See:
  # https://stackoverflow.com/questions/25007324/can-ggplot2-control-point-size-and-line-size-lineweight-separately-in-one-lege
  build <- ggplot_build(legplot)
  gt <- ggplot_gtable(build)
  getgtable<-function(n) which(sapply(gt$grobs, `[[`, "name")==n)
  getgtable("guide-box")
  segs <- grepl("GRID.segments", sapply(gt$grobs[[getgtable("guide-box")]][[1]][[1]]$grobs, '[[', "name") )
  gt$grobs[[getgtable("guide-box")]][[1]][[1]]$grobs[segs]<-lapply(gt$grobs[[getgtable("guide-box")]][[1]][[1]]$grobs[segs], 
                                                                   function(x) {x$gp$lwd<-1; x})
  
  gt$grobs[[getgtable("guide-box")]][[1]][[1]]$grobs[which(segs)[2]][[1]]$gp$lwd <- 2
  
  pnts <- which(grepl("GRID.points", sapply(gt$grobs[[getgtable("guide-box")]][[1]][[1]]$grobs, '[[', "name") ))
  gt$grobs[[getgtable("guide-box")]][[1]][[1]]$grobs[pnts[2]][[1]]$pch <- ""
  
  # grid.draw(gt)
  
  leg <- get_legend(gt)
  
  ###########################################################################3
  
  twoplots <- plot_grid( p1,p2,ncol=1, rel_heights = c(1,1))
  outp <- ggdraw(twoplots) + draw_plot(leg,-0.325,0.41,1,1) + draw_plot(leg,0.38,-0.09,1,1)
  print(outp)

}

################################################################################################
################################################################################################
################################################################################################

### Plotting function for the supplementary figure comparing different # of supp sources ###


# black and white
plot.emotbymode_supp_bw <- function(plotdatfinlst,primids)
{
  # plotdatfinlst <- pdatlst
  # primids <- primids
  
  shift <- 0.07
  dodge <- position_dodge(width = 0.9)
  
  
  ##########################################################################
  # first plot
  
  plotdatfin <- plotdatfinlst[[1]]
  primid <- primids[1]
  
  limits <- aes(ymax = plotdatfin$upr, ymin = plotdatfin$lwr,x=as.numeric(outcm)+shift) # right: q=12
  limits7 <- aes(ymax=plotdatfin$upr7,ymin=plotdatfin$lwr7,x = as.numeric(outcm))      # middle: q=7
  limits2 <- aes(ymax=plotdatfin$upr2,ymin=plotdatfin$lwr2,x=as.numeric(outcm)-shift)  # left: q=2
  
  line_types <- c("iMEM (q = 2)"=3,"iMEM (q = 7)"=2,"iMEM (q = 12)"=1)

  
  p1 <- ggplot(data = plotdatfin,aes(x = outcm, y = postmean,fill=`Trip Mode`)) # bar height: q=12
  
  p1 <- p1 + geom_bar( stat = "identity", position = dodge) + 
    geom_point( aes(x = as.numeric(outcm), y = postmean7),stat = "identity", position = dodge,size=rel(3.5),shape=2) +
    geom_point( aes(x = as.numeric(outcm)-shift, y = postmean2),stat = "identity", position = dodge,size=rel(3.5),shape=4) +
    scale_y_continuous(name = "Intensity") + 
    coord_cartesian(ylim=c(0,6.2)) + 
    scale_fill_manual(values=c("grey40","grey60","grey80")) +
    scale_linetype_manual(name="Blah",values=line_types)+
    
    geom_errorbar(limits7, position = dodge, width = 0.15,size=0.4,linetype=2) +
    geom_errorbar(limits, position=dodge, width = 0.15,size=0.7,col="black") +
    geom_errorbar(limits2, position=dodge, width = 0.15,size=0.4,col="black",linetype=3) +
    
    # geom_text(aes(x=outcm,y=0.5,label=round(esss,0)),stat="identity",position=dodge,size=6) +
    
    guides(fill = guide_legend(override.aes = list(shape = NA) ) )     +
    
    theme_bw() +
    
    theme(axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_text(size=rel(1.9)),
          axis.title.y=element_text(size = rel(1.9)),
          axis.text.x = element_text(size = rel(2.5)),
          axis.title.x= element_blank(),
          legend.position="bottom" ,
          legend.text=element_text(size=rel(1.6)),
          legend.title=element_text(size=rel(1.9)),
          plot.title =element_text(size=20, face='bold')
          
    )  +
    
    ggtitle(paste0("Emotion Intensity by Trip Mode: ID ",primid)) 

  ##########################################################################
  # second plot
  
  plotdatfin2 <- plotdatfinlst[[2]]
  primid2 <- primids[2]
  
  limits_2 <- aes(ymax = plotdatfin2$upr, ymin = plotdatfin2$lwr,x=as.numeric(outcm)+shift) # right: q=12
  limits7_2 <- aes(ymax=plotdatfin2$upr7,ymin=plotdatfin2$lwr7,x = as.numeric(outcm))      # middle: q=7
  limits2_2 <- aes(ymax=plotdatfin2$upr2,ymin=plotdatfin2$lwr2,x=as.numeric(outcm)-shift)  # left: q=2
  
  line_types <- c("iMEM (q = 2)"=3,"iMEM (q = 7)"=2,"iMEM (q = 12)"=1)
  
  
  p2 <- ggplot(data = plotdatfin2,aes(x = outcm, y = postmean,fill=`Trip Mode`)) # bar height: q=12
  
  p2 <- p2 + geom_bar( stat = "identity", position = dodge) + 
    geom_point( aes(x = as.numeric(outcm), y = postmean7),stat = "identity", position = dodge,size=rel(3.5),shape=2) +
    geom_point( aes(x = as.numeric(outcm)-shift, y = postmean2),stat = "identity", position = dodge,size=rel(3.5),shape=4) +
    scale_y_continuous(name = "Intensity") + 
    coord_cartesian(ylim=c(0,6.2)) + 
    scale_fill_manual(values=c("grey40","grey60","grey80")) +
    scale_linetype_manual(name="Blah",values=line_types)+
    
    geom_errorbar(limits7_2, position = dodge, width = 0.15,size=0.4,linetype=2) +
    geom_errorbar(limits_2, position=dodge, width = 0.15,size=0.7,col="black") +
    geom_errorbar(limits2_2, position=dodge, width = 0.15,size=0.4,col="black",linetype=3) +
    
    # geom_text(aes(x=outcm,y=0.5,label=round(esss,0)),stat="identity",position=dodge,size=6) +
    
    guides(fill = guide_legend(override.aes = list(shape = NA) ) )     +
    
    theme_bw() +
    
    theme(axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_text(size=rel(1.9)),
          axis.title.y=element_text(size = rel(1.9)),
          axis.text.x = element_text(size = rel(2.5)),
          axis.title.x= element_blank(),
          legend.position="bottom" ,
          legend.text=element_text(size=rel(1.6)),
          legend.title=element_text(size=rel(1.9)),
          plot.title =element_text(size=20, face='bold')
          
    )  +
    
    ggtitle(paste0("Emotion Intensity by Trip Mode: ID ",primid2)) 
  
  
  ########## mock plot for legend extraction #####################
  line_types <- c("iMEM (q = 2)"=3,"iMEM (q = 7)"=2,"iMEM (q = 12)"=1)
  sizes <- c("iMEM (q = 2)"=1,"iMEM (q = 7)"=1,"iMEM (q = 12)"=1)
  shps <- c("iMEM (q = 2)"=4,"iMEM (q = 7)"=2,"iMEM (q = 12)"=1)
  legdat <- data.frame(y = 1:9,x = 9:1,grp=rep(c("iMEM (q = 2)", "iMEM (q = 7)","iMEM (q = 12)"),3))
  legdat$grp <- factor(legdat$grp,levels(legdat$grp)[c(2,3,1)])
  
  legplot <- ggplot(data=legdat, aes(x=x,y=y,linetype=grp,size=grp,shape=grp)) + geom_line() + geom_point() +
    scale_linetype_manual(name="",values=line_types) +
    scale_size_manual(name="",values=sizes) +
    scale_shape_manual(name="",values=shps) +
    guides(shape=guide_legend(override.aes=list(size=3))) + theme_bw() + 
    theme(legend.text=element_text(size=rel(1.3)))
  
  # customize the legend (separte sizes for points and lines) by editing the low level grid objects. See:
  # https://stackoverflow.com/questions/25007324/can-ggplot2-control-point-size-and-line-size-lineweight-separately-in-one-lege
  build <- ggplot_build(legplot)
  gt <- ggplot_gtable(build)
  getgtable<-function(n) which(sapply(gt$grobs, `[[`, "name")==n)
  getgtable("guide-box")
  segs <- grepl("GRID.segments", sapply(gt$grobs[[getgtable("guide-box")]][[1]][[1]]$grobs, '[[', "name") )
  gt$grobs[[getgtable("guide-box")]][[1]][[1]]$grobs[segs]<-lapply(gt$grobs[[getgtable("guide-box")]][[1]][[1]]$grobs[segs], 
                                              function(x) {x$gp$lwd<-1; x})
  
  gt$grobs[[getgtable("guide-box")]][[1]][[1]]$grobs[max(which(segs))][[1]]$gp$lwd <- 2
  
  pnts <- which(grepl("GRID.points", sapply(gt$grobs[[getgtable("guide-box")]][[1]][[1]]$grobs, '[[', "name") ))
  gt$grobs[[getgtable("guide-box")]][[1]][[1]]$grobs[max(pnts)][[1]]$pch <- ""
  
  # grid.draw(gt)
  
  leg <- get_legend(gt)
  
  ###########################################################################3
  
  twoplots <- plot_grid( p1,p2,ncol=1, rel_heights = c(1,1))
  outp <- ggdraw(twoplots) + draw_plot(leg,-0.35,0.41,1,1) + draw_plot(leg,0.40,-0.09,1,1)
  outp
  
  
  print(outp)
}


################################################################################################
################################################################################################
################################################################################################

### Plotting function for the supplementary figure comparing different threshold values ###


# black and white
plot.emotbymode_supp_bw2 <- function(plotdatfin,primid)
{
  shift <- 0.07
  dodge <- position_dodge(width = 0.9)
  limits0.05 <- aes(ymax = plotdatfin$upr0.05, ymin = plotdatfin$lwr0.05,x=as.numeric(outcm)+shift)
  limits0.1 <- aes(ymax=plotdatfin$upr0.1,ymin=plotdatfin$lwr0.1,x = as.numeric(outcm))
  limits0.2 <- aes(ymax=plotdatfin$upr0.2,ymin=plotdatfin$lwr0.2,x=as.numeric(outcm)-shift)
  
  p <- ggplot(data = plotdatfin,aes(x = outcm, y = postmean0.1,fill=`Trip Mode`)) 
  
  p <- p + geom_bar( stat = "identity", position = dodge) + 
    geom_point( aes(x = as.numeric(outcm)+shift, y = postmean0.05),stat = "identity", position = dodge,size=rel(3.5),shape=2) +
    geom_point( aes(x = as.numeric(outcm)-shift, y = postmean0.2),stat = "identity", position = dodge,size=rel(3.5),shape=4) +
    scale_y_continuous(name = "Intensity") + 
    coord_cartesian(ylim=c(0,6.2)) + 
    scale_fill_manual(values=c("grey40","grey60","grey80")) +

    geom_errorbar(limits0.1, position = dodge, width = 0.15,size=0.7) +
    geom_errorbar(limits0.2, position=dodge, width = 0.15,size=0.4,col="black",linetype=2) +
    geom_errorbar(limits0.05, position=dodge, width = 0.15,size=0.4,col="black",linetype=3) +
    
    # geom_text(aes(x=outcm,y=0.5,label=round(esss,0)),stat="identity",position=dodge,size=6) +
    
    guides(fill = guide_legend(override.aes = list(shape = NA) ) )     +
    
    theme_bw() +
    
    theme(axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_text(size=rel(1.9)),
          axis.title.y=element_text(size = rel(1.9)),
          axis.text.x = element_text(size = rel(2.5)),
          axis.title.x= element_blank(),
          legend.position="bottom" ,
          legend.text=element_text(size=rel(1.6)),
          legend.title=element_text(size=rel(1.9)),
          plot.title =element_text(size=20, face='bold')
          
    )  +
    
    ggtitle(paste0("Emotion Intensity by Trip Mode: ID ",primid)) 
  
  print(p)
}


################################################################################################
################################################################################################
################################################################################################

### Generate the plotting data ###
set.seed(123)
outcms <- c("happy","meaningful","stressful","tired")
pdat3096 <- gendat.emotbymode("3096",outcms)
pdat2018 <- gendat.emotbymode("2018",outcms)


### calculate average/max percent deviation from mean and percent dcrease in posterior standard deviation
mean(c(pdat3096$pctdecsd,pdat2018$pctdecsd))
max(c(pdat3096$pctdecsd,pdat2018$pctdecsd))
mean(c(pdat3096$pctdiffmean,pdat2018$pctdiffmean))
max(c(pdat3096$pctdiffmean,pdat2018$pctdiffmean))


pdatlst <- list(pdat2018,pdat3096)
primids <- c(2018,3096)


export <- T

### Make the primary plots ###
# if(export) { pdf(paste0(root,"iMEM_emotionbymode_3096_v4.pdf"),height=7,width=13) }
# plot.emotbymode_bw(pdat3096,"3096")
# if(export) {dev.off()}
# 
# if(export) { pdf(paste0(root,"iMEM_emotionbymode_2018_v4.pdf"),height=7,width=13) }
# plot.emotbymode_bw(pdat2018,"2018")
# if(export) {dev.off()}

### Manuscript Figure 5 ###
if(export) { pdf(paste0(root,"iMEM_emotionbymode_bothids_v4.pdf"),height=14/1.25,width=13/1.25,onefile=F) }
plot.emotbymode_bw(pdatlst,primids)
if(export) {dev.off()}

# 
# ### Make the supplementary plots ###
# if(export) { pdf(paste0(root,"iMEM_emotionbymode_3096_v4_supp.pdf"),height=7,width=13) }
# plot.emotbymode_supp_bw(pdat3096,"3096")
# if(export) {dev.off()}
# 
# if(export) { pdf(paste0(root,"iMEM_emotionbymode_2018_v4_supp.pdf"),height=7,width=13) }
# plot.emotbymode_supp_bw(pdat2018,"2018")
# if(export) {dev.off()}


### Supplementary Figure 8 ###
if(export) { pdf(paste0(root,"iMEM_emotionbymode_bothids_v4_supp.pdf"),height=14/1.25,width=13/1.25,onefile=F) }
plot.emotbymode_supp_bw(pdatlst,primids)
if(export) {dev.off()}


# ### Make the second set of supplementary plots ###
# if(export) { pdf(paste0(root,"iMEM_emotionbymode_3096_v4_thresh.pdf"),height=7,width=13) }
# plot.emotbymode_supp_bw2(pdat3096,"3096")
# if(export) {dev.off()}
# 
# if(export) { pdf(paste0(root,"iMEM_emotionbymode_2018_v4_thresh.pdf"),height=7,width=13) }
# plot.emotbymode_supp_bw2(pdat2018,"2018")
# if(export) {dev.off()}

