

library(lattice)
root <- "C:/Users/roland/Documents/iMEM Manuscript Data and Code/"
export <- T

###################################################################################################

###### DEFINE SOME FUNCTIONS TO QUIKCLY CACLULATE IMEM SIMILARITY SCORES #######

onesrcwt <- function( primvals, supp1vals)
{
	#######################################################################################
	#  Function to quickly calculate similarity score in the presence of only one supplementary source
  # 
	#  Source values specified as c(mean,v), where v is sd/n for that source
	#######################################################################################
	
	# primvals <- c(2,4/20)
	# supp1vals <- c(2,4/20)
	
	m <- primvals[1]; v <- primvals[2]
	m1 <- supp1vals[1]; v1 <- supp1vals[2]
	
	# marginal models
	p0 <- (2*pi)/sqrt(1/(v*v1))
	p1 <- (sqrt(2*pi)/sqrt(1/v + 1/v1))*exp(-0.5*((m-m1)^2/(v+v1)) )
	
	# inclusion weight for source 1
	return( p1/(p0+p1) )
}


twosrcwt_v2 <- function( supp2vals, primvals, supp1vals)
{
  #######################################################################################
  #  Function to quickly calculate similarity for Source 1  in the presence of two supplemetnary sources
  # 
  #  Source values specified as c(mean,v), where v is sd/n for that source
  #######################################################################################
  
  m <- primvals[1]; v <- primvals[2]
	m1 <- supp1vals[1]; v1 <- supp1vals[2]
	m2 <- supp2vals[1]; v2 <- supp2vals[2]
	
	# marginal models
	p0 <- (sqrt(2*pi)^3)/sqrt( 1/(v*v1*v2))
	p1 <- ((2*pi)/sqrt( (1/v + 1/v1)*(1/v2) )) * exp(-0.5*((m-m1)^2/(v+v1)) )
	p2 <- ((2*pi)/sqrt( (1/v + 1/v2)*(1/v1) )) * exp(-0.5*((m-m2)^2/(v+v2)) )
	p12 <- (sqrt(2*pi)/sqrt(1/v + 1/v1 + 1/v2)) * exp(-0.5*( ((m-m1)^2/(v+v1+v*v1*(1/v2))) + ((m-m2)^2/(v+v2+v*v2*(1/v1))) + ((m1-m2)^2/(v2+v1+v2*v1*(1/v))) ) )
	
	# inclusion weight for source 1
	return( (p1+p12)/ (p0+p1+p2+p12) )
}



###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

##### Setting the values and grids of values for primary source and two supplementary sources ####
pmean <- 0
s2ms <- s1ms <- seq( pmean - 3, pmean + 3, length.out=300)

pldf <- expand.grid(s1ms=s1ms,s2ms=s2ms)
pldf$s1wt <- NA

# calculate weights in presence of Source 2
for(i in 1:nrow(pldf))
{
  primvals <- c(pmean,0.2)
  s1vals <- c(pldf$s1ms[i],0.2)
  s2vals <- c(pldf$s2ms[i],0.01)
  
  pldf$s1wt[i] <- twosrcwt_v2(supp2vals=s2vals,primvals=primvals,supp1vals=s1vals)
}


# calculate marginal weights, and merge into joint weights dataframe
s1df <- data.frame(s1ms=s1ms,s1mwts=NA)
for(i in 1:nrow(s1df)) { s1df$s1mwts[i] <- onesrcwt(primvals=c(pmean,0.2),supp1vals=c(s1df$s1ms[i],0.2))}
pldf <- merge(pldf,s1df,by="s1ms",all.x=T)

pldf$wt_diff <- pldf$s1wt-pldf$s1mwts
pldf$jw_gt_mw <- 1*(pldf$wt_diff > 0)



##############################################################################

##### GENERATE THE HEATMAPS ###############


if(export) { pdf( paste0(root,"Paper Graphics/Inclusion Weights by Means_v2.pdf"),height=7,width=13) }

# Plotting raw joint weights
mn <- min(pldf$s1wt); mx <- max(pldf$s1wt)
col.l <- colorRampPalette(c('white','yellow','green', 'blue','purple'))
ncols <- (14*4)+1

colorBreaks1 <- seq(mn,mx,length.out=ncols)  
trellis.par.set(regions=list(col=col.l(ncols)))
p1 <- 					
  levelplot( s1wt ~ s1ms + s2ms, data=pldf , xlab="Source 1 Mean", ylab="Source 2 Mean",
             main=expression(paste("T"["1"])),
             scales=list(x=list(at=-3:3), y=list(at=-3:3 )),
             at=colorBreaks1,				
             panel=function(...) {
               arg <- list(...)
               panel.levelplot(...)
             })

# print(p)

print(p1, position=c(0, 0, 0.48, 1), more=TRUE)


# Plotting joint weights minus marginal weights
mn <- min(pldf$wt_diff); mx <- max(pldf$wt_diff)
col.l <- colorRampPalette(c('red','white','blue'))
ncols <- (14*4)+1
tol <- 1e-7   

colorBreaks1 <- c( seq(mn - 0.001, 0 - tol,length.out = (ncols+1)/2 ), # adding a small buffer on end      
                   seq(0 + tol, mx + 0.001,length.out = (ncols+1)/2))  

# pdf( paste0(root,"graphics/Inclusion Weights by Means_v1.pdf"))
trellis.par.set(regions=list(col=col.l(ncols)))
p2 <- 					
  levelplot( wt_diff ~ s1ms + s2ms, data=pldf , xlab="Source 1 Mean", ylab="Source 2 Mean",
             main= expression(paste("T"["1"],"-T"["1m"]) ),
             scales=list(x=list(at=-3:3), y=list(at=-3:3 )),
             at=colorBreaks1,				
             panel=function(...) {
               arg <- list(...)
               panel.levelplot(...)
             })

# print(p)
print(p2, position=c(0.52, 0, 1, 1))

if(export){ dev.off() }

