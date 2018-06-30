#!/usr/bin/Rscript --vanilla
# Function: Calculate indexes and draw heatmap for CSR data from HH2.
# Time: 3/14/2014
# Author: Junli Zhang <zhjl86@gmail.com>

# arguments: nrow ncol filename
# example
# hh2csr.R 15 24 AM_Dry.txt
cat("\nYour three arguments should be: nrow, ncol, filename !!!\n\n")

args <- commandArgs(trailingOnly = TRUE)
nrow <- as.numeric(args[1])
ncol <- as.numeric(args[2])
filename <- args[3]
#dd <- read.delim(filename) #for normal file format
dd <- read.delim(filename, header=FALSE, skip=3) # for T3 format
if (nrow*ncol != ncol(dd[-1])) cat("\nERROR: Plot number is not a whole rectangle !!!\n\n")

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
write.table(indices, paste("Indices_",filename,".txt",sep=""),sep="\t",row.names=FALSE)

###### heat plot to show field variation
library(lattice)
pdf(file=paste("heatmap_",filename,".pdf",sep=""), width=ifelse(nrow < ncol,11,8.5), height=ifelse(nrow < ncol,8.5,11))
#pdf(file=paste("heatmap_",filename,".pdf",sep=""), width=27, height=27)
for (i in 1:ncol(indices)){
    # get one trait and transform it to a matrix
    m1 <- matrix(indices[,i],nrow=nrow)#default is by column
    # transform to serpentene
    for (j in 1:ncol) {if (j%%2 == 0) m1[,j] <- rev(m1[,j])}
    # make heatmap
    print(levelplot(t(m1),col.regions=terrain.colors,xlab="Col",ylab="Row",main=names(indices)[i]),cex=2)
}
dev.off()
cat("\nHeatmap PDF was created!\n")

###### chech index changes with time
wl <- c(472, 680, 780, 850, 880, 900, 920, 970)
pdf(file=paste("Time_change_",filename,".pdf",sep=""), width=11, height=8.5)
# Reflectance range
plot(0,0,xlim=range(dd[,1]),ylim=range(dd[-1]),type="n",ylab="Reflectance",xlab="Wavelength (nm)", xaxt="n")
axis(side=1,at=seq(ceiling(min(dd[,1])/50), floor(max(dd[,1])/50))*50)
for (i in 2:ncol(dd)) lines(dd[,1],dd[,i],type="l",col="blue")
for (i in wl) abline(v=i,lty=1,col="green")
# Time changes for each index
for (i in 1:ncol(indices)) plot(indices[,i], pch = i, col="red", xlab="Plot", ylab=names(indices)[i])
dev.off()
cat("\nTime change PDF was created!\n")
