#!/usr/bin/Rscript --vanilla
# Function: Calculate indexes and draw heatmap for CSR data from HH2 in Davis 2015 plot designs.
# Time: 3/20/2015
# Author: Junli Zhang <zhjl86@gmail.com>
## changes
	# v2: add block order as input

## basic block layout for Davis 2015
#	Drought		|	Irr
#	3	4		|
#	2	6	8	|	10	12
#	1	5	7	|	9	11

# arguments: nrow in each block, ncolin each block, nblock, filename
# example
# hh2csr.R 15 24 AM_Dry.txt
#cat("\nYour three arguments should be: nrow, ncol, nblock, filename !!!\n\n")
cat("\nYour two arguments should be: filename, block_order !!!\n\n")

args <- commandArgs(trailingOnly = TRUE)
#nrow <- as.numeric(args[1])
#ncol <- as.numeric(args[2])
#nblock <- as.numeric(args[3])
#filename <- args[4]
filename <- args[1]
blocks <- as.numeric(unlist(strsplit(args[2], ","))) # block order during measurement
nrow = 8
ncol = 12
#nblock = 12
nblock = length(blocks)


dd <- read.delim(filename) #for normal file format
#dd <- read.delim(filename, header=FALSE, skip=3) # for T3 format
if (nrow*ncol*nblock != ncol(dd[-1])) cat("\nERROR: Plot number is not a whole rectangle !!!\n\n")

# remove .txt from file name
filename <- sub(".txt","",filename)

###### Calculate indices
window <- 2
#wl.need <- c(445, 472, 531, 550, 570, 670, 680, 780, 850, 880, 900, 920, 970)

wl.mean <- function(x, window=2){
	left <- x - window
	right <- x + window
	n1 <- findInterval(left,dd[,1])
	n2 <- findInterval(right,dd[,1])
	r.mean <- colMeans(dd[n1:n2,-1])
	return(r.mean)
}

R445 <- wl.mean(445)
R472 <- wl.mean(472) # Blue
R531 <- wl.mean(531)
R550 <- wl.mean(550)
R560 <- wl.mean(560)
R570 <- wl.mean(570)
R670 <- wl.mean(670)
R680 <- wl.mean(680)
R705 <- wl.mean(705)
R720 <- wl.mean(720)
R750 <- wl.mean(750)
R780 <- wl.mean(780)
R790 <- wl.mean(790)
R810 <- wl.mean(810)
R850 <- wl.mean(850)
R880 <- wl.mean(880)
R900 <- wl.mean(900)
R920 <- wl.mean(920)
R970 <- wl.mean(970)

# biomass, nitrogen, chlorophyll

SR <- R900/R680
PRI <- (R531-R570)/(R531+R570)
RNDVI <- (R780-R670)/(R780+R670)
GNDVI <- (R780-R550)/(R780+R550)
NDVI <- (R900-R680)/(R900+R680)
# EVI = 2.5*(NIR-Red)/(NIR+6*Red-7.5*Blue+1)
EVI1.1 <- 2.5*(R780-R670)/(R780+6*R670-7.5*R472+1)
EVI1.2 <- 2.5*(R810-R670)/(R810+6*R670-7.5*R472+1)
# EVI2 = 2.5*((NIR-Red)/(NIR+2.4*Red+1))
EVI2 <- 2.5*(R780-R670)/(R780+2.4*R670+1)
# chlorophyll
MSR <- (R750-R445)/(R705-R445)
NDRE <- (R790-R720)/(R790+R720)
# nitrogen index from Rice
NIO <- R810/R560



# water indieces
WI <- R970/R900
NWI1 <- (R970-R900)/(R970+R900)
NWI2 <- (R970-R850)/(R970+R850)
NWI3 <- (R970-R880)/(R970+R880)
NWI4 <- (R970-R920)/(R970+R920)

indices <- data.frame(SR, PRI, NDVI, RNDVI, GNDVI, EVI1.1, EVI1.2, EVI2, MSR, NDRE, NIO, WI, NWI1, NWI2, NWI3, NWI4)
names(indices) <- paste(names(indices),filename,sep="_") # add filename to each index
#write.table(indices, paste("Indices_",filename,".txt",sep=""),sep="\t",row.names=FALSE)
cat("\nIndices were calculated successfully!\n\n")

## two borders
longborder = matrix(NA, nrow=24, ncol=1) # long block
shortborder = matrix(NA, nrow=16, ncol=1) # short block

#longb = matrix(NA, nrow=24, ncol=25)
#shortb = matrix(NA, nrow=16, ncol=38)
blankb = matrix(NA, nrow=8, ncol=12)
###### heat plot to show field variation

library(lattice)
pdf(file=paste("heatmap_",filename,"_v2.pdf",sep=""), width=11, height=8.5)
for (i in 1:ncol(indices)){
    # get one trait and transform it to a matrix
    l1 <- split(indices[,i], rep(1:nblock, each=nrow*ncol))
    l2 <- lapply(l1, function(x){
    	m1 <- matrix(x,nrow=nrow) #default is by column
    	# transform to serpentene
    	for (j in 1:ncol) {if (j%%2 == 0) m1[,j] <- rev(m1[,j])}
    	return(m1)
    })
    
    # below code is necessary in case only part of blocks were measured due to reasons like bad weather
    l3 = rep(list(blankb), 12) # 12 blank blocks
    for (k in 1:nblock) l3[blocks[k]] = l2[k] # fill the blocks with values
    
    # arrange the blocks
    longb = cbind(do.call("rbind",l3[1:3]), longborder, do.call("rbind",l3[c(5,6,4)]))
    # the last plot in the first population is severly damaged, its indices are like soil's, so I remove it.
    longb[1,12] = NA
    shortb = cbind(do.call("rbind",l3[7:8]), shortborder, do.call("rbind",l3[9:10]), shortborder, do.call("rbind",l3[11:12]))

    
    # make heatmap
    #print(levelplot(t(mm),col.regions=terrain.colors,xlab="Col",ylab="Row",main=names(indices)[i]),cex=2)
    if (! all(is.na(longb))){
    	print(levelplot(t(longb),col.regions=terrain.colors,xlab="Col",ylab="Row",main=names(indices)[i]),cex=2)
    }
    if (! all(is.na(shortb))){
    	print(levelplot(t(shortb),col.regions=terrain.colors,xlab="Col",ylab="Row",main=names(indices)[i]),cex=2)
    }
}
dev.off()
cat("\nHeatmap PDF was created!\n")

###### chech index changes with time
wl <- c(472, 680, 780, 850, 880, 900, 920, 970)
pdf(file=paste("Time_change_",filename,"_v2.pdf",sep=""), width=11, height=8.5)
# Reflectance range
plot(0,0,xlim=range(dd[,1]),ylim=range(dd[-1]),type="n",ylab="Reflectance",xlab="Wavelength (nm)", xaxt="n")
axis(side=1,at=seq(ceiling(min(dd[,1])/50), floor(max(dd[,1])/50))*50)
for (i in 2:ncol(dd)) lines(dd[,1],dd[,i],type="l",col="blue")
for (i in wl) abline(v=i,lty=1,col="green")
# Time changes for each index
for (i in 1:ncol(indices)) plot(indices[,i], pch = i, col="red", xlab="Plot", ylab=names(indices)[i])
dev.off()
cat("\nTime change PDF was created!\n")


# output indices with correct plot order
b1 = 101:196
plot.list = vector("list", 12)
plot.v = c(1:8, 11:14)
for (i in 1:12) plot.list[[i]] = plot.v[i] * 100 + 1:96
plot.order = unlist(plot.list[blocks])
indices = data.frame(Plot=plot.order, indices)
write.table(indices, paste("Indices_",filename,"_with_plots.txt",sep=""),sep="\t",row.names=FALSE)

