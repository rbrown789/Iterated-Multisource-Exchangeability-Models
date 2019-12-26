

# Simulation 1 Results

library(abind)
library(RColorBrewer)
root <- "C:/Users/rbrow_000/Documents/UofM Documents/Dissertation/Iterated-Multisource-Exchangeability-Models/"
root <- "C:/Users/rbrow/OneDrive/Documents/iMEM Manuscript Data and Code/"
simdata <- "Data/Simulation Results/"
code <- "Code/"

means.along <- function(a, i) {
  n <- length(dim(a))
  b <- aperm(a, c(seq_len(n)[-i], i))
  rowMeans(b, dims = n - 1)
}

###########################################################



load(paste0(root,simdata,"scen1res_v2b.rdata"))
load(paste0(root,simdata,"scen2res_v2b.rdata"))
load(paste0(root,simdata,"scen3res_v2b.rdata"))
load(paste0(root,simdata,"scen4res_v2b.rdata"))
load(paste0(root,simdata,"scen5res_v2b.rdata"))
load(paste0(root,simdata,"scen6res_v2b.rdata"))
load(paste0(root,simdata,"scen7res_v2b.rdata"))
load(paste0(root,simdata,"scen8res_v2b.rdata"))

############################################################



# combine MEM results into single 5D array for processing, then remove the separate result arrays
resar <- abind(scen1simres,scen2simres,along=5)
resar <- abind(resar,scen3simres,along=5)
resar <- abind(resar,scen4simres,along=5)
resar <- abind(resar,scen5simres,along=5)
resar <- abind(resar,scen6simres,along=5)
resar <- abind(resar,scen7simres,along=5)
resar <- abind(resar,scen8simres,along=5)
rm(list=paste0("scen",1:8,"simres"))

names(dimnames(resar)) <- c("mod","quant","sim","mu","scen")
dimnames(resar)[[5]] <- paste0("scen",1:8)

# interval length, and coverage indicator
intl <- resar[,"ciup",,,]-resar[,"cilo",,,]
cover <- 1*( resar[,"cilo",,,] <= resar[,"mu",,,] & resar[,"ciup",,,] >= resar[,"mu",,,])

# append to results array and drop unneeded arrays
resar <- aperm(resar,c(1,3,4,5,2))
resar <- abind(resar,intl,cover,along=5)
rm(list=c("intl","cover"))
dimnames(resar)[[5]] <- c("mu","model","postmean","postvar","cilo","ciup","ess","intl","cover")

# calculate coverage rate, bias, variance, mse
res <- means.along(resar[,,,,c("mu","model","postmean","ess","intl","cover")],2)
bias <- res[,,,"postmean"]-res[,,,"mu"]
var <- resar[,,,,"postmean"]
var <- apply(var, c(1,3,4), var)
mse <- bias^2 + var

res <- abind(res,bias,var,mse,along=4)
names(dimnames(res)) <- c("model","mu","scen","quant")
dimnames(res)[[4]] <- c("mu","model","postmean","ess","intl","cover","bias","var","mse")


##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################

# Now we're going to put the random effects models in the same format
load(paste0(root,simdata,"scen1reres_v2.rdata"))
load(paste0(root,simdata,"scen2reres_v2.rdata"))
load(paste0(root,simdata,"scen3reres_v2.rdata"))
load(paste0(root,simdata,"scen4reres_v2.rdata"))
load(paste0(root,simdata,"scen5reres_v2.rdata"))
load(paste0(root,simdata,"scen6reres_v2.rdata"))
load(paste0(root,simdata,"scen7reres_v2.rdata"))
load(paste0(root,simdata,"scen8reres_v2.rdata"))

dimm <- dim(scen1reres)
dimnms <- dimnames(scen1reres)

scen1reres <- array(  unlist(scen1reres),dim=dimm,dimnames=dimnms )
scen2reres <- array(  unlist(scen2reres),dim=dimm,dimnames=dimnms )
scen3reres <- array(  unlist(scen3reres),dim=dimm,dimnames=dimnms )
scen4reres <- array(  unlist(scen4reres),dim=dimm,dimnames=dimnms )
scen5reres <- array(  unlist(scen5reres),dim=dimm,dimnames=dimnms )
scen6reres <- array(  unlist(scen6reres),dim=dimm,dimnames=dimnms )
scen7reres <- array(  unlist(scen7reres),dim=dimm,dimnames=dimnms )
scen8reres <- array(  unlist(scen8reres),dim=dimm,dimnames=dimnms )



# combine random effects results into single 5D array for processing, then remove the separate result arrays
resar <- abind(scen1reres,scen2reres,along=5)
resar <- abind(resar,scen3reres,along=5)
resar <- abind(resar,scen4reres,along=5)
resar <- abind(resar,scen5reres,along=5)
resar <- abind(resar,scen6reres,along=5)
resar <- abind(resar,scen7reres,along=5)
resar <- abind(resar,scen8reres,along=5)
rm(list=paste0("scen",1:8,"reres"))

names(dimnames(resar)) <- c("mod","quant","sim","mu","scen")
dimnames(resar)[[5]] <- paste0("scen",1:8)

# calculate interval length and coverage indicator
intl <- resar[,"ciup",,,]-resar[,"cilo",,,]
cover <- 1*( resar[,"cilo",,,] <= resar[,"mu",,,] & resar[,"ciup",,,] >= resar[,"mu",,,])
pvar <- array(NA,dim=dim(cover),dimnames=dimnames(cover))

# append to main array and drop unneeded arrays, then rearrange to get in same format as the other array
resar <- aperm(resar,c(1,3,4,5,2))
resar <- abind(resar,pvar,intl,cover,along=5)
rm(list=c("intl","cover"))
dimnames(resar)[[5]] <- c("mu","model","postmean","cilo","ciup","ess","postvar","intl","cover")
resar <- resar[,,,,c("mu","model","postmean","postvar","ess","cilo","ciup","intl","cover")]

# calculate coverage rate, bias, variance, mse
res.re <- means.along(resar[,,,,c("mu","model","postmean","ess","intl","cover")],2)
bias <- res.re[,,,"postmean"]-res.re[,,,"mu"]
var <- resar[,,,,"postmean"]
var <- apply(var, c(1,3,4), var)
mse <- bias^2 + var

res.re <- abind(res.re,bias,var,mse,along=4)
names(dimnames(res.re)) <- c("model","mu","scen","quant")
dimnames(res.re)[[4]] <- c("mu","model","postmean","ess","intl","cover","bias","var","mse")

##################################################################################################
##################################################################################################
##################################################################################################


# Now we're going to put the general imem model results in the same format
load(paste0(root,simdata,"scen1grandres_v2.rdata"))
load(paste0(root,simdata,"scen2grandres_v2.rdata"))
load(paste0(root,simdata,"scen3grandres_v2.rdata"))
load(paste0(root,simdata,"scen4grandres_v2.rdata"))
load(paste0(root,simdata,"scen5grandres_v2.rdata"))
load(paste0(root,simdata,"scen6grandres_v2.rdata"))
load(paste0(root,simdata,"scen7grandres_v2.rdata"))
load(paste0(root,simdata,"scen8grandres_v2.rdata"))


# combine into single 5D array for processing, then remove the separate result arrays
resar <- abind(scen1grandres,scen2grandres,along=5)
resar <- abind(resar,scen3grandres,along=5)
resar <- abind(resar,scen4grandres,along=5)
resar <- abind(resar,scen5grandres,along=5)
resar <- abind(resar,scen6grandres,along=5)
resar <- abind(resar,scen7grandres,along=5)
resar <- abind(resar,scen8grandres,along=5)
rm(list=paste0("scen",1:8,"grandres"))

names(dimnames(resar)) <- c("mod","quant","sim","mu","scen")
dimnames(resar)[[5]] <- paste0("scen",1:8)

# interval length, and coverage indicator
intl <- resar[,"ciup",,,]-resar[,"cilo",,,]
cover <- 1*( resar[,"cilo",,,] <= resar[,"mu",,,] & resar[,"ciup",,,] >= resar[,"mu",,,])

# append to results array and drop unneeded arrays
resar <- aperm(resar,c(1,3,4,5,2))
resar <- abind(resar,intl,cover,along=5)
rm(list=c("intl","cover"))
dimnames(resar)[[5]] <- c("mu","model","postmean","postvar","cilo","ciup","ess","intl","cover")

# calculate coverage rate, bias, variance, mse
res.gr <- means.along(resar[,,,,c("mu","model","postmean","ess","intl","cover")],2)
bias <- res.gr[,,,"postmean"]-res.gr[,,,"mu"]
var <- resar[,,,,"postmean"]
var <- apply(var, c(1,3,4), var)
mse <- bias^2 + var

res.gr <- abind(res.gr,bias,var,mse,along=4)
names(dimnames(res.gr)) <- c("model","mu","scen","quant")
dimnames(res.gr)[[4]] <- c("mu","model","postmean","ess","intl","cover","bias","var","mse")


##################################################################################################
##################################################################################################
##################################################################################################


# Now we're going to put the general imem model results in the same format
load(paste0(root,simdata,"scen1margthresh_v2.rdata"))
load(paste0(root,simdata,"scen2margthresh_v2.rdata"))
load(paste0(root,simdata,"scen3margthresh_v2.rdata"))
load(paste0(root,simdata,"scen4margthresh_v2.rdata"))
load(paste0(root,simdata,"scen5margthresh_v2.rdata"))
load(paste0(root,simdata,"scen6margthresh_v2.rdata"))
load(paste0(root,simdata,"scen7margthresh_v2.rdata"))
load(paste0(root,simdata,"scen8margthresh_v2.rdata"))


# combine into single 5D array for processing, then remove the separate result arrays
resar <- abind(scen1margthresh,scen2margthresh,along=5)
resar <- abind(resar,scen3margthresh,along=5)
resar <- abind(resar,scen4margthresh,along=5)
resar <- abind(resar,scen5margthresh,along=5)
resar <- abind(resar,scen6margthresh,along=5)
resar <- abind(resar,scen7margthresh,along=5)
resar <- abind(resar,scen8margthresh,along=5)
rm(list=paste0("scen",1:8,"margthresh"))

names(dimnames(resar)) <- c("mod","quant","sim","mu","scen")
dimnames(resar)[[5]] <- paste0("scen",1:8)

# interval length, and coverage indicator
intl <- resar[,"ciup",,,]-resar[,"cilo",,,]
cover <- 1*( resar[,"cilo",,,] <= resar[,"mu",,,] & resar[,"ciup",,,] >= resar[,"mu",,,])

# append to results array and drop unneeded arrays
resar <- aperm(resar,c(1,3,4,5,2))
resar <- abind(resar,intl,cover,along=5)
rm(list=c("intl","cover"))
dimnames(resar)[[5]] <- c("mu","model","postmean","postvar","cilo","ciup","ess","intl","cover")

# calculate coverage rate, bias, variance, mse
res.thresh <- means.along(resar[,,,,c("mu","model","postmean","ess","intl","cover")],2)
bias <- res.thresh[,,,"postmean"]-res.thresh[,,,"mu"]
var <- resar[,,,,"postmean"]
var <- apply(var, c(1,3,4), var)
mse <- bias^2 + var

res.thresh <- abind(res.thresh,bias,var,mse,along=4)
names(dimnames(res.thresh)) <- c("model","mu","scen","quant")
dimnames(res.thresh)[[4]] <- c("mu","model","postmean","ess","intl","cover","bias","var","mse")

##################################################################################################
##################################################################################################
##################################################################################################


# Now we're going to put the general imem model results in the same format
load(paste0(root,simdata,"scen1multq_v2.rdata"))
load(paste0(root,simdata,"scen2multq_v2.rdata"))
load(paste0(root,simdata,"scen3multq_v2.rdata"))
load(paste0(root,simdata,"scen4multq_v2.rdata"))
load(paste0(root,simdata,"scen5multq_v2.rdata"))
load(paste0(root,simdata,"scen6multq_v2.rdata"))
load(paste0(root,simdata,"scen7multq_v2.rdata"))
load(paste0(root,simdata,"scen8multq_v2.rdata"))


# combine into single 5D array for processing, then remove the separate result arrays
resar <- abind(scen1multq,scen2multq,along=5)
resar <- abind(resar,scen3multq,along=5)
resar <- abind(resar,scen4multq,along=5)
resar <- abind(resar,scen5multq,along=5)
resar <- abind(resar,scen6multq,along=5)
resar <- abind(resar,scen7multq,along=5)
resar <- abind(resar,scen8multq,along=5)
rm(list=paste0("scen",1:8,"multq"))

names(dimnames(resar)) <- c("mod","quant","sim","mu","scen")
dimnames(resar)[[5]] <- paste0("scen",1:8)

# interval length, and coverage indicator
intl <- resar[,"ciup",,,]-resar[,"cilo",,,]
cover <- 1*( resar[,"cilo",,,] <= resar[,"mu",,,] & resar[,"ciup",,,] >= resar[,"mu",,,])

# append to results array and drop unneeded arrays
resar <- aperm(resar,c(1,3,4,5,2))
resar <- abind(resar,intl,cover,along=5)
rm(list=c("intl","cover"))
dimnames(resar)[[5]] <- c("mu","model","postmean","postvar","cilo","ciup","ess","intl","cover")

# calculate coverage rate, bias, variance, mse
res.multq <- means.along(resar[,,,,c("mu","model","postmean","ess","intl","cover")],2)
bias <- res.multq[,,,"postmean"]-res.multq[,,,"mu"]
var <- resar[,,,,"postmean"]
var <- apply(var, c(1,3,4), var)
mse <- bias^2 + var

res.multq <- abind(res.multq,bias,var,mse,along=4)
names(dimnames(res.multq)) <- c("model","mu","scen","quant")
dimnames(res.multq)[[4]] <- c("mu","model","postmean","ess","intl","cover","bias","var","mse")


################################

# combine MEM results random effect results into final array
res.fin <- abind(res,res.re,along=1)
res.fin <- abind(res.fin,res.gr,along=1)
res.fin <- abind(res.fin,res.thresh,along=1)
res.fin <- abind(res.fin,res.multq,along=1)

####################################################################################################

# convert to data frame for plotting
todf <- function(resarray)
{
  ##############################################################
  # convert all information on all quantitites into a dataframe
  # with columns for simulation settings 
  #
  #################################################################
  
  out <- NULL
  for(i in dimnames(resarray)[[2]])
  {
    for(j in dimnames(resarray)[[3]] )
    {
      tmp <- data.frame( resarray[,i,j,] )
      tmp$model <- row.names(tmp)
      tmp$scen <- j
      out <- rbind(out,tmp)
    }
  }
  return(out)
}

resdf <- todf(res.fin)

###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################

############## MAKING PLOTS #########################
ver <- 1
export <- T

##### plotting function ######


pfunc3 <- function(scen=1,quant,quantlab,mods,modlabs,cols,ltys,lwds=1,ylim,legloc="topright",legcex=1.2,dat)
{
  # scen <- 6
  # quant <- "cover"
  # dat <- rdfp
  # quantlab="Coverage Rate"
  # mods <- c("imem4","gimemrand_u2_r0.5_q4","gimemrand_u2_r0.75_q4","gimemrand_u3_r0.5_q4","gimemrand_u3_r0.75_q4",
  #           "gimemrand_u4_r0.5_q4","gimemrand_u4_r0.75_q4")
  # 
  # modlabs <- c("Marginal iMEM",paste0("Gen. iMEM (",c("u=2, r=0.5)","u=2, r=0.75)","u=3, r=0.5)","u=3, r=0.75)","u=4, r=0.5)","u=4, r=0.75)")))
  # 
  # cols <- c("darkviolet", brewer.pal(n=length(mods)-1,"RdBu"))
  # ltys <- lwds <- rep(1,length(mods))
  # ylim  <- c(0.8,1)
  
  
  ###################
  # set the sample mean locations
  if(scen == 1) { smlocs <- -2 }
  if(scen == 2) { smlocs <- c(-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6) }
  if(scen == 3) { smlocs <- c(3.9, -2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7, 4.1) }
  if(scen == 4) { smlocs <- c(3.9, 3.8, -2.2,-2.1,-2,-1.9,-1.8, 4.2, 4.1) }
  if(scen == 5) { smlocs <- c(-2.5,-1.5,1,1.2,1.3,2.3,2.8,3,3.8) }
  if(scen == 6) { smlocs <- seq(-2,4,length.out=9) }
  
  ###################
  
  sctxt <- paste0("scen",scen)
  tmp <- dat[dat$scen %in% sctxt,]
  
  xmin <- min(tmp$mu); xmax <- max(tmp$mu)

  plot("n",xlim=c(xmin,xmax),ylim=ylim,ylab=quantlab,xlab="mu",main=paste0("Scenario ",scen),cex.axis=1.2,cex.lab=1.2,cex.main=1.2)
  
  
  for(i in 1:length(mods))
  {
    x <- tmp$mu[tmp$mod%in% mods[i]]
    y <- tmp[[quant]][ tmp$mod%in% mods[i] ]
    lines(x,y,col=cols[i],lty=ltys[i],lwd=lwds[i])
  }
  
  abline(v=smlocs,col="grey",lty=2)
  legend(legloc,legend=modlabs,col=cols,lty=ltys,lwd=lwds,bty="n",cex=legcex)
  
  
  
  # if(quant=="bias") { abline(h=0,col="grey",lty=2) }
  # if(quant=="cover") { abline(h=0.95,col="grey",lty=2) }
}




# grab only scenarios 1,4,5,6,7,8, and rename to 1,2,3,4,5,6 #
rdfp <- resdf[!resdf$scen %in% c("scen2","scen3"),]
nums <- c(1,4,5,6,7,8)
for(i in 1:length(nums)) { rdfp$scen[rdfp$scen %in% paste0("scen",nums[i])] <- paste0("scen",i) }

# set supplementary source mean locations
xbars1 <- rep(-2,9)
xbars2 <- c(-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6)
xbars3 <- c(3.9, -2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7, 4.1)
xbars4 <- c(3.9, 3.8, -2.2,-2.1,-2,-1.9,-1.8, 4.2, 4.1)
xbars5 <- c(-2.5,-1.5,1,1.2,1.3,2.3,2.8,3,3.8)
xbars6 <- seq(-2,4,length.out=9)


############### MEM vs. iMEM ###################
mods <- c("imem2","imem4","imem7","mem")
modlabs <- c(paste0("Marg. iMEM (q=",c("2)","4)","7)")),"MEM")

cols <- brewer.pal(n=length(mods),"RdBu")
cols <- c("black","grey30","grey60","grey85")
ltys <- lwds <- rep(1,length(mods))



### MSE - Manuscript Figure 4 ###
if(export) { pdf( paste0(root,"memVimem_MSE_v",ver,".pdf"),height=10.5,width=10) }
par(mfrow=c(3,2))
sapply(1:6,pfunc3,quant="mse",quantlab="MSE",mods=mods,modlabs=modlabs,cols=cols,ltys=ltys,lwds=lwds,
       ylim=c(0,2),legloc="topright",dat=rdfp)
if(export){dev.off()}

### ESSS - Manuscript Figure 3 ###
if(export) { pdf( paste0(root,"memVimem_ESSS_v",ver,".pdf"),height=10.5,width=10) }
par(mfrow=c(3,2))
sapply(1:6,pfunc3,quant="ess",quantlab="ESSS",mods=mods,modlabs=modlabs,cols=cols,ltys=ltys,lwds=lwds,
       ylim=c(0,50),legloc="topright",dat=rdfp)
if(export){dev.off()}


### Interval Length - Supplementary Figure 3 ###
if(export) { pdf( paste0(root,"memVimem_intLength_v",ver,".pdf"),height=10.5,width=10) }
par(mfrow=c(3,2))
sapply(1:6,pfunc3,quant="intl",quantlab="Average Interval Length",mods=mods,modlabs=modlabs,cols=cols,ltys=ltys,lwds=lwds,
       ylim=c(1.5,3.5),legloc="bottomright",dat=rdfp)
if(export){dev.off()}


### Coverage ###
if(export) { pdf( paste0(root,"memVimem_coverage_v",ver,".pdf"),height=10.5,width=10) }
par(mfrow=c(3,2))
sapply(1:6,pfunc3,quant="cover",quantlab="Coverage Rate",mods=mods,modlabs=modlabs,cols=cols,ltys=ltys,lwds=lwds,
       ylim=c(0.55,1),legloc="bottomright",dat=rdfp)
if(export){dev.off()}



############### General vs. Marginal Selection ###################

mods <- c("imem4","gimemrand_u2_r0.5_q4","gimemrand_u2_r0.75_q4","gimemrand_u3_r0.5_q4","gimemrand_u3_r0.75_q4",
                    "gimemrand_u4_r0.5_q4","gimemrand_u4_r0.75_q4")

modlabs <- c("Marginal iMEM",paste0("Gen. iMEM (",c("u=2, r=0.5)","u=2, r=0.75)","u=3, r=0.5)","u=3, r=0.75)","u=4, r=0.5)","u=4, r=0.75)")))

cols <- c("darkviolet", brewer.pal(n=length(mods)-1,"RdBu"))
ltys <- lwds <- rep(1,length(mods))


### ESSS - Supplementary Figure 4 ###
if(export) { pdf( paste0(root,"Paper Graphics/gVm_ESSS_v",ver,".pdf"),height=10.5,width=10) }
par(mfrow=c(3,2))
sapply(1:6,pfunc3,quant="ess",quantlab="ESSS",mods=mods,modlabs=modlabs,cols=cols,ltys=ltys,lwds=lwds,
       ylim=c(0,27),legloc="topright",legcex=0.8,dat=rdfp)
if(export){dev.off()}





############### Simple Mean and Random Effects ###################
mods <- c("imem7","mem","repred","simpmean")
modlabs <- c("iMEM (q=7)","MEM","RE Prediction","Simple Mean")

cols <- c("black","steelblue","orange","red")
ltys <- lwds <- rep(1,length(mods))

### MSE - Supplementary Figure 6 ###
if(export) { pdf( paste0(root,"Paper Graphics/compmeth_MSE_v",ver,".pdf"),height=10.5,width=10) }
par(mfrow=c(3,2))
sapply(1:6,pfunc3,quant="mse",quantlab="MSE",mods=mods,modlabs=modlabs,cols=cols,ltys=ltys,lwds=lwds,
       ylim=c(0,5.5),legloc="topright",dat=rdfp)
if(export){dev.off()}

### Interval Length - Supplementary Figure 5 ###
if(export) { pdf( paste0(root,"Paper Graphics/compmeth_intlength_v",ver,".pdf"),height=10.5,width=10) }
par(mfrow=c(3,2))
sapply(1:6,pfunc3,quant="intl",quantlab="Average Interval Length",mods=mods,modlabs=modlabs,cols=cols,ltys=ltys,lwds=lwds,
       ylim=c(1.5,3.8),legloc="topright",dat=rdfp)
if(export){dev.off()}



############### Multiple q iMEM ###################
mods <- paste0("imemmultq_",c("2_7","2_9","3_7","3_9","4_7","4_9"))
modlabs <- paste0("iMEM (qmin = ",c(2,2,3,3,4,4),", qmax = ",c(7,9,7,9,7,9),")")

cols <- brewer.pal(n=length(mods),"RdBu")
# cols <- c("black","grey30","grey60","grey85")
ltys <- lwds <- rep(1,length(mods))


### Coverage ###
if(export) { pdf( paste0(root,"imemmultq_coverage_v",ver,".pdf"),height=10.5,width=10) }
par(mfrow=c(3,2))
sapply(1:6,pfunc3,quant="cover",quantlab="Coverage Rate",mods=mods,modlabs=modlabs,cols=cols,ltys=ltys,lwds=lwds,
       ylim=c(0.55,1),legloc="bottomright",dat=rdfp)
if(export){dev.off()}

### MSE ###
if(export) { pdf( paste0(root,"imemmultq_MSE_v",ver,".pdf"),height=10.5,width=10) }
par(mfrow=c(3,2))
sapply(1:6,pfunc3,quant="mse",quantlab="MSE",mods=mods,modlabs=modlabs,cols=cols,ltys=ltys,lwds=lwds,
       ylim=c(0,2),legloc="bottomright",dat=rdfp)
if(export){dev.off()}

### Interval Length ###
if(export) { pdf( paste0(root,"imemmultq_intLength_v",ver,".pdf"),height=10.5,width=10) }
par(mfrow=c(3,2))
sapply(1:6,pfunc3,quant="intl",quantlab="Average Interval Length",mods=mods,modlabs=modlabs,cols=cols,ltys=ltys,lwds=lwds,
       ylim=c(1.5,3.8),legloc="bottomright",dat=rdfp)
if(export){dev.off()}


############### Multiple q and Single q iMEM ###################
mods <- c( paste0("imemmultq_",c("2_9","4_7")),paste0("imem",c(2,4,7)),"mem" ) 
modlabs <- c(paste0("Marg. iMEM (q = ",c(2,4),"-",c(9,7),")"), paste0("Marg. iMEM (q=",c("2)","4)","7)")),"MEM")
cols <- c("black","grey60",c("black","grey30","grey60","grey85"))
ltys <- c( rep(1,2),rep(6,4) )
lwds <- c( rep(1.3,2),rep(0.78,4) )


### Coverage ###
if(export) { pdf( paste0(root,"imemmultsingq_coverage_v",ver,".pdf"),height=10.5,width=10) }
par(mfrow=c(3,2))
sapply(1:6,pfunc3,quant="cover",quantlab="Coverage Rate",mods=mods,modlabs=modlabs,cols=cols,ltys=ltys,lwds=lwds,
       ylim=c(0.55,1),legloc="bottomright",dat=rdfp)
if(export){dev.off()}

### MSE ###
if(export) { pdf( paste0(root,"imemmultsingq_MSE_v",ver,".pdf"),height=10.5,width=10) }
par(mfrow=c(3,2))
sapply(1:6,pfunc3,quant="mse",quantlab="MSE",mods=mods,modlabs=modlabs,cols=cols,ltys=ltys,lwds=lwds,
       ylim=c(0,2),legloc="topright",dat=rdfp)
if(export){dev.off()}

### Interval Length ###
if(export) { pdf( paste0(root,"imemmultsingq_intLength_v",ver,".pdf"),height=10.5,width=10) }
par(mfrow=c(3,2))
sapply(1:6,pfunc3,quant="intl",quantlab="Average Interval Length",mods=mods,modlabs=modlabs,cols=cols,ltys=ltys,lwds=lwds,
       ylim=c(1.5,3.8),legloc="bottomright",dat=rdfp)
if(export){dev.off()}


############### thresholded iMEM ###################
mods <- paste0("imemthresh",c(0.05,0.1,0.15,0.2))
modlabs <- paste0("iMEM (Score > ",c(.05,.1,.15,.2),")")

cols <- brewer.pal(n=length(mods),"RdBu")
cols <- c("black","grey30","grey60","grey85")
ltys <- lwds <- rep(1,length(mods))


### Coverage ###
if(export) { pdf( paste0(root,"imemthresh_coverage_v",ver,".pdf"),height=10.5,width=10) }
par(mfrow=c(3,2))
sapply(1:6,pfunc3,quant="cover",quantlab="Coverage Rate",mods=mods,modlabs=modlabs,cols=cols,ltys=ltys,lwds=lwds,
       ylim=c(0.55,1),legloc="bottomright",dat=rdfp)
if(export){dev.off()}

### MSE ###
if(export) { pdf( paste0(root,"imemthresh_MSE_v",ver,".pdf"),height=10.5,width=10) }
par(mfrow=c(3,2))
sapply(1:6,pfunc3,quant="mse",quantlab="MSE",mods=mods,modlabs=modlabs,cols=cols,ltys=ltys,lwds=lwds,
       ylim=c(0,2),legloc="bottomright",dat=rdfp)
if(export){dev.off()}

### ESSS ###
if(export) { pdf( paste0(root,"imemthresh_ESSS_v",ver,".pdf"),height=10.5,width=10) }
par(mfrow=c(3,2))
sapply(1:6,pfunc3,quant="ess",quantlab="ESSS",mods=mods,modlabs=modlabs,cols=cols,ltys=ltys,lwds=lwds,
       ylim=c(0,50),legloc="topright",dat=rdfp)
if(export){dev.off()}


### Interval Length ###
if(export) { pdf( paste0(root,"imemthresh_intLength_v",ver,".pdf"),height=10.5,width=10) }
par(mfrow=c(3,2))
sapply(1:6,pfunc3,quant="intl",quantlab="Average Interval Length",mods=mods,modlabs=modlabs,cols=cols,ltys=ltys,lwds=lwds,
       ylim=c(1.5,3.5),legloc="bottomright",dat=rdfp)
if(export){dev.off()}




#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################

########### SECONDARY SIMULATION PLOT ###############################

# Load data
load(paste0(root,simdata,"simb_gen_6_v2.rdata"))
load(paste0(root,simdata,"simb_genrand_6_v2.rdata"))
load(paste0(root,simdata,"simb_gen_2to5_v2.rdata"))
load(paste0(root,simdata,"simb_genrand_2to5_v2.rdata"))
load(paste0(root,simdata,"simb_marg_v2.rdata"))


# combine MEM results into single 4D array for processing, then remove the separate result arrays
resar <- abind(simb_marg,simb_gen_2to5,along=1)
resar <- abind(resar,simb_gen_6,along=1)
resar <- abind(resar,simb_genrand_2to5,along=1)
resar <- abind(resar,simb_genrand_6,along=1)

rm(list=c("simb_marg","simb_gen_2to5","simb_gen_6","simb_genrand_2to5","simb_genrand_6"))

names(dimnames(resar)) <- c("mod","quant","sim","mu")

# interval length, and coverage indicator
intl <- resar[,"ciup",,]-resar[,"cilo",,]
cover <- 1*( resar[,"cilo",,] <= resar[,"mu",,] & resar[,"ciup",,] >= resar[,"mu",,])


# append to results array and drop unneeded arrays
resar <- aperm(resar,c(1,3,4,2))
resar <- abind(resar,intl,cover,along=4)
rm(list=c("intl","cover"))
dimnames(resar)[[4]] <- c("mu","model","postmean","postvar","cilo","ciup","ess","intl","cover")

# calculate coverage rate, bias, variance, mse
res <- means.along(resar[,,,c("mu","model","postmean","ess","intl","cover")],2)
bias <- res[,,"postmean"]-res[,,"mu"]
var <- resar[,,,"postmean"]
var <- apply(var, c(1,3), var)
mse <- bias^2 + var

res <- abind(res,bias,var,mse,along=3)
names(dimnames(res)) <- c("model","mu","quant")
dimnames(res)[[3]] <- c("mu","model","postmean","ess","intl","cover","bias","var","mse")


# convert to data frame for plotting
todfb <- function(resarray)
{
  ##############################################################
  # convert all information on all quantitites into a dataframe
  # with columns for simulation settings 
  #
  #################################################################
  
  out <- NULL
  for(i in dimnames(resarray)[[2]])
  {
    tmp <- data.frame( resarray[,i,] )
    tmp$model <- row.names(tmp)
    out <- rbind(out,tmp)
  }
  return(out)
}
resdf <- todfb(res)

############################################################


pfunc4 <- function(quant,quantlab,mods,modlabs,cols,ltys,lwds,ylim,legloc="topright",dat)
{
  
  # quant <- "ess"
  # dat <- resdf
  # quantlab="ESSS"
  # mods <- c("imem8","gimemrand_u2_r0.5_q8","gimemrand_u2_r0.75_q8","gimemrand_u4_r0.5_q8","gimemrand_u4_r0.75_q8",
  #           "gimemrand_u6_r0.5_q8","gimemrand_u6_r0.75_q8")
  # 
  # mods <- c("imem8","gimem_u2_r0.5_q8","gimem_u2_r0.75_q8","gimem_u4_r0.5_q8","gimem_u4_r0.75_q8",
  #           "gimem_u6_r0.5_q8","gimem_u6_r0.75_q8")
  # 
  # modlabs <- c("Marginal iMEM",paste0("Gen. iMEM (",c("u=2, r=0.5)","u=2, r=0.75)","u=4, r=0.5)","u=4, r=0.75)","u=6, r=0.5)","u=6, r=0.75)")))
  # 
  # cols <- c("darkviolet", brewer.pal(n=length(mods)-1,"RdBu"))
  # ltys <- lwds <- rep(1,length(mods))
  # ylim  <- c(0,45)
  # legloc <- "topright"
  
  ###################
  # set the sample mean locations
  xbars <- c( -3,-2.9,-2.8,-2.7,-2.6,rep(-2,5), 0,0.2,0.4,0.6,0.8,rep(1.3,10),3,3.3,3.6,3.9,rep(4.5,3)    )
  xbarsmult <- xbars[duplicated(xbars)]
  
  ##################
  tmp <- dat[dat$mod %in% mods,]
  
  xmin <- min(tmp$mu); xmax <- max(tmp$mu)
  
  plot("n",xlim=c(xmin,9),ylim=ylim,ylab=quantlab,xlab="mu",main=quantlab)
  
  
  for(i in 1:length(mods))
  {
    x <- tmp$mu[tmp$mod%in% mods[i]]
    y <- tmp[[quant]][ tmp$mod%in% mods[i] ]
    lines(x,y,col=cols[i],lty=ltys[i])
  }
  
  legend(legloc,legend=modlabs,col=cols,lty=ltys,lwd=lwds,bty="n",cex=0.8)
  
  abline(v=xbars,col="grey",lty=2)
  abline(v=xbarsmult,col="grey",lty=1,lwd=1.5)
  
  for(i in unique(xbarsmult) ) 
  { 
    n <- sum(xbars %in% i)
    text <- paste0("x",n)
    yloc <- ylim[2]
    xloc <- i
    text(xloc,yloc,text)
  }
  
  # if(quant=="bias") { abline(h=0,col="grey",lty=2) }
  # if(quant=="cover") { abline(h=0.95,col="grey",lty=2) }
}


mods <- c("imem8","gimem_u2_r0.5_q8","gimem_u2_r0.75_q8","gimem_u4_r0.5_q8","gimem_u4_r0.75_q8",
          "gimem_u6_r0.5_q8","gimem_u6_r0.75_q8")

mods <- c("imem8","gimemrand_u2_r0.5_q8","gimemrand_u2_r0.75_q8","gimemrand_u4_r0.5_q8","gimemrand_u4_r0.75_q8",
          "gimemrand_u6_r0.5_q8","gimemrand_u6_r0.75_q8")

modlabs <- c("Marginal iMEM",paste0("Gen. iMEM (",c("u=2, r=0.5)","u=2, r=0.75)","u=4, r=0.5)",
                                                    "u=4, r=0.75)","u=6, r=0.5)","u=6, r=0.75)")))


# if(export) { pdf( paste0(root,"Paper Graphics/simb_results_",ver,".pdf"),height=10,width=13) }
# par(mfrow=c(2,2))
# 
# pfunc4(quant="ess",quantlab="Average ESSS",mods=mods,modlabs=modlabs
#        ,cols = c("darkviolet", brewer.pal(n=length(mods)-1,"RdBu"))
#        ,ltys = rep(1,length(mods)),lwds=rep(1,length(mods)),ylim=c(0,45) 
#        ,dat=resdf
# )
# 
# pfunc4(quant="intl",quantlab="Average Interval Length",mods=mods,modlabs=modlabs,legloc="bottomright"
#        ,cols = c("darkviolet", brewer.pal(n=length(mods)-1,"RdBu"))
#        ,ltys = rep(1,length(mods)),lwds=rep(1,length(mods)),ylim=c(1.5,4.4) 
#        ,dat=resdf
# )
# 
# 
# pfunc4(quant="cover",quantlab="Coverage Rate",mods=mods,modlabs=modlabs,legloc="bottomright"
#        ,cols = c("darkviolet", brewer.pal(n=length(mods)-1,"RdBu"))
#        ,ltys = rep(1,length(mods)),lwds=rep(1,length(mods)),ylim=c(0.3,0.97) 
#        ,dat=resdf
# )
# 
# pfunc4(quant="mse",quantlab="MSE",mods=mods,modlabs=modlabs,legloc="bottomright"
#        ,cols = c("darkviolet", brewer.pal(n=length(mods)-1,"RdBu"))
#        ,ltys = rep(1,length(mods)),lwds=rep(1,length(mods)),ylim=c(0,1.5) 
#        ,dat=resdf
# )
# 
# if(export) { dev.off() }
# 

#### Supplementary Figure 7 ####
if(export) { pdf( paste0(root,"simb_MSE_ESSS",ver,".pdf"),height=10,width=8) }
par(mfrow=c(2,1))

pfunc4(quant="ess",quantlab="Average ESSS",mods=mods,modlabs=modlabs
       ,cols = c("darkviolet", brewer.pal(n=length(mods)-1,"RdBu"))
       ,ltys = rep(1,length(mods)),lwds=rep(1,length(mods)),ylim=c(0,45) 
       ,dat=resdf
)


pfunc4(quant="mse",quantlab="MSE",mods=mods,modlabs=modlabs,legloc="bottomright"
       ,cols = c("darkviolet", brewer.pal(n=length(mods)-1,"RdBu"))
       ,ltys = rep(1,length(mods)),lwds=rep(1,length(mods)),ylim=c(0,1.5) 
       ,dat=resdf
)

if(export) { dev.off() }




