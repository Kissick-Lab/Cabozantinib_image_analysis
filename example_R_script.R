###S20-23167-A6
cd4=read.csv('D:/Cabo_mIF_QuPath/QuPath/S20-23167-A6/CD4.csv')
cd8=read.csv('D:/Cabo_mIF_QuPath/QuPath/S20-23167-A6/CD8.csv')
mhc=read.csv('D:/Cabo_mIF_QuPath/QuPath/S20-23167-A6/MHC-II.csv')
dapi=read.csv('D:/Cabo_mIF_QuPath/QuPath/S20-23167-A6/DAPI.csv')

#merge to 1 data frame 
cells=rbind(cd8,cd4,mhc)
xlim=c(min(cells[6]),max(cells[6]))
ylim=c(-max(cells[7]),-min(cells[7]))
names(cells)

#notes to self about which column is which parameter (fluorophore=marker=column number)
#DAPI=14
#570=aSMA=20
#690=CD8=26
#480=MHC=32
#620=FOX=38
#780=CD4=44
#520=TCF=50

#plot CD4 vs CD8

plot(log(cells[,26]+1,10),log(cells[,44]+1,10),pch=1,cex=0.0001,col=rgb(0.4,0.4,0.4,0.2),xlim=c(0.0,2.6),ylim=c(0.0,2.6))

cells_CD4_neg_tf=apply(cells, 1, function(x) x[44]<2.4)

cells_CD4_neg=cells[cells_CD4_neg_tf,] 

plot(log(cells_CD4_neg[,26]+1,10),log(cells_CD4_neg[,44]+1,10),pch=1,cex=0.0001,col=rgb(0.4,0.4,0.4,0.2),xlim=c(0.0,2.6),ylim=c(0.0,2.6))

cd8_cells=cells[log(cells[26]+1,10)>1.25 & log(cells[44]+1,10)<1.6,]
cd4_cells=cells_CD4_neg[log(cells_CD4_neg[44]+1,10)>1.6 & log(cells_CD4_neg[26]+1,10)<1.25,]

pdf("/Users/BTVO/Desktop/K404_cd4xcd8_gateout.pdf")
plot(log(cells_CD4_neg[,26]+1,10),log(cells_CD4_neg[,44]+1,10),pch=1,cex=0.1,col=rgb(0.4,0.4,0.4,0.2),xlim=c(0.0,2.6),ylim=c(0.0,2.6))
par(new=T)
plot(log(cd8_cells[,26]+1,10),log(cd8_cells[,44]+1,10),pch=1,cex=0.1,col=rgb(1,0.4,0.4,0.2),xlim=c(0.0,2.6),ylim=c(0.0,2.6))
par(new=T)
plot(log(cd4_cells[,26]+1,10),log(cd4_cells[,44]+1,10),pch=1,cex=0.1,col=rgb(0.1,0.4,1,0.2),xlim=c(0.0,2.6),ylim=c(0.0,2.6))
abline(h=c(0,1,2),v=c(0,1,2),lwd=1,lty=3,col='grey50')
dev.off()

nrow(cd8_cells)
num_row_cd8=nrow(cd8_cells) 
num_row_DAPI=nrow(dapi)
percentage=(num_row_cd8/num_row_DAPI)*100
print(percentage)
# %CD8 = 16.23816
# FACS %CD8 = 

#CD8+ cells per mm2
xdim2=max(cd8_cells[6])
ydim2=max(cd8_cells[7])

#to convert to mm
xdim2_mm=(xdim2*0.4962)/1000
ydim2_mm=(ydim2*0.4962)/1000

#to get the area
area2=xdim2_mm*xdim2_mm
nCD8=nrow(cd8_cells)
amtCD8=nCD8/area2
amtCD8
#902.6362


#Number of CD4 positive cells in gated box
nrow(cd4_cells)
#  

#calculate CD4% of DAPI
num_row_cd4=nrow(cd4_cells) 
num_row_DAPI=nrow(dapi)
percentage=(num_row_cd4/num_row_DAPI)*100
print(percentage)
# 12.78741

#plot cd8 things to look at what you need to comp
plot(log(cd8_cells[,26]+1,10),log(cd8_cells[,44]+1,10),pch=1,cex=0.0001,col=rgb(1,0.4,0.4,0.2),xlim=c(0.0,2.6),ylim=c(0.0,2.6))
plot(log(cd8_cells[,26]+1,10),log(cd8_cells[,50]+1,10),pch=1,cex=0.0001,col=rgb(1,0.4,0.4,0.2),xlim=c(0.0,2.6),ylim=c(0.0,2.6))

#estimate your comps
cd8_cells_comp=cd8_cells
cd8_cells_comp[26]=cd8_cells_comp[26]-0.1*cd8_cells_comp[50]
cd8_cells_comp[50]=cd8_cells_comp[50]-0.3*cd8_cells_comp[26]

#check how you've comped with a plot 
plot(log(cd8_cells[,26]+1,10),log(cd8_cells[,50]+1,10),pch=1,cex=0.0001,col=rgb(1,0.4,0.4,0.2),xlim=c(0.0,2.6),ylim=c(0.0,2.6))
plot(log(cd8_cells_comp[,26]+1,10),log(cd8_cells_comp[,50]+1,10),pch=1,cex=0.0001,col=rgb(1,0.4,0.4,0.2),xlim=c(0.0,2.6),ylim=c(0.0,2.6))

#decide/delineate which ones are tcf+, then plot them 
tcf_cd8=cd8_cells_comp[cd8_cells_comp[50]>15,]
not_tcf_cd8=cd8_cells_comp[cd8_cells_comp[50]<20,]

#Plot CD8 and TCF1+
plot(log(cd8_cells_comp[,26]+1,10),log(cd8_cells_comp[,50]+1,10),pch=1,cex=0.0001,col=rgb(1.0,0.4,0.4,0.2),xlim=c(0.0,2.6),ylim=c(0.0,2.6))
par(new=T)
plot(log(tcf_cd8[,26]+1,10),log(tcf_cd8[,50]+1,10),pch=1,cex=0.0001,col=rgb(0.1,.8,0.2,0.8),xlim=c(0.0,2.6),ylim=c(0.0,2.6))
abline(h=c(0,1,2),v=c(0,1,2),lwd=1,lty=3,col='grey50')
dev.off()

#Number of TCF1+CD8+cells 
nrow(tcf_cd8)
#42200

#calculate %TCF1+ of CD8
num_row_tcf_cd8=nrow(tcf_cd8) 
num_row_cd8=nrow(cd8_cells)
percentage=(num_row_tcf_cd8/num_row_cd8)*100
print(percentage)
# 18.4204

#plot all CD8s and TCF1+ CD8s as immunomap -- with their locations 
plot(data.frame(cd8_cells_comp[6],-cd8_cells_comp[7]),pch=1,cex=0.0001,col=rgb(1.0,0.4,0.4,0.2),xlim=xlim,ylim=ylim)
par(new=T)
plot(data.frame(tcf_cd8[6],-tcf_cd8[7]),pch=1,cex=0.0001,col=rgb(0.1,.8,0.2,0.8),xlim=xlim,ylim=ylim)
par(new=T)
plot(data.frame(not_tcf_cd8[6],-not_tcf_cd8[7]),pch=1,cex=0.0001,col='blue3',xlim=xlim,ylim=ylim)
abline(v=c(5000,10000,15000),h=c(-5000,-15000),lwd=1,lty=3,col='grey50')
dev.off()

#CD4 and FoxP3
plot(log(cd4_cells[,44]+1,10),log(cd4_cells[,38]+1,10),pch=1,cex=0.0001,col=rgb(0.1,0.4,1,0.2),xlim=c(0.0,2.5),ylim=c(0.0,2.5))

cd4_cells_comp=cd4_cells

cd4_cells_comp[38]=cd4_cells_comp[38]-.5*cd4_cells_comp[26]
cd4_cells_comp[26]=cd4_cells_comp[26]-0.30*cd4_cells_comp[38]

plot(log(cd4_cells_comp[,44]+1,10),log(cd4_cells_comp[,38]+1,10),pch=1,cex=0.0001,col=rgb(0.1,0.4,1,0.2),xlim=c(0,2.5),ylim=c(0,2.5))

fox_cd4=cd4_cells_comp[cd4_cells_comp[38]>15,]
not_fox_cd4=cd4_cells_comp[cd4_cells_comp[38]<18,]

#Plot CD4 and FoxP3+
plot(log(cd4_cells_comp[,44]+1,10),log(cd4_cells_comp[,38]+1,10),pch=1,cex=0.0001,col=rgb(0.1,0.4,1,0.2),xlim=c(0.0,2.5),ylim=c(0,2.5))
par(new=T)
plot(log(fox_cd4[,44]+1,10),log(fox_cd4[,38]+1,10),pch=1,cex=0.0001,col='deeppink3',xlim=c(0.0,2.5),ylim=c(0,2.5))
abline(h=c(.5,1.5),v=c(.5,1,1.5,2,2.5),lwd=1,lty=3,col='grey50')
dev.off()

#calculate %FOXP3 of CD4
num_row_fox_cd4=nrow(fox_cd4) 
num_row_cd4=nrow(cd4_cells)
percentage=(num_row_fox_cd4/num_row_cd4)*100
print(percentage)
#10.71384

#Plot immunomap of CD4 and FoxP3+
plot(data.frame(cd4_cells_comp[6],-cd4_cells_comp[7]),pch=1,cex=0.001,col=rgb(0.1,0.4,1,0.2),xlim=xlim,ylim=ylim)
par(new=T)
plot(data.frame(fox_cd4[6],-fox_cd4[7]),pch=1,cex=0.001,col='deeppink3',xlim=xlim,ylim=ylim)
abline(v=c(5000,10000,15000),h=c(-5000,-15000),lwd=1,lty=3,col='grey50')
dev.off()

#CD4 and TCF1+
plot(log(cd4_cells[,44]+1,10),log(cd4_cells[,50]+1,10),pch=1,cex=0.0001,col=rgb(0.1,0.4,1,0.2),xlim=c(0.0,2.5),ylim=c(0.0,2.5))

cd4_cells_comp=cd4_cells

cd4_cells_comp[50]=cd4_cells_comp[50]-0.30*cd4_cells_comp[26]

plot(log(cd4_cells_comp[,44]+1,10),log(cd4_cells_comp[,50]+1,10),pch=1,cex=0.0001,col=rgb(0.1,0.4,1,0.2),xlim=c(0.0,2.5),ylim=c(0.0,2.5))

tcf_cd4=cd4_cells_comp[cd4_cells_comp[50]>20,]
not_tcf_cd4=cd4_cells_comp[cd4_cells_comp[50]<25,]

#Plot CD4 and TCF1+
plot(log(cd4_cells_comp[,44]+1,10),log(cd4_cells_comp[,50]+1,10),pch=1,cex=0.0001,col=rgb(0.1,0.4,1,0.2),xlim=c(0.0,2.5),ylim=c(0.0,2.5))
par(new=T)
plot(log(tcf_cd4[,44]+1,10),log(tcf_cd4[,50]+1,10),pch=1,cex=0.0001,col=rgb(0.1,.8,0.2,0.8),xlim=c(0.0,2.5),ylim=c(0.0,2.5))
abline(h=c(0,1,2),v=c(0,1,2),lwd=1,lty=3,col='grey50')
dev.off()

#calculate %TCF1+ of CD4
num_row_tcf_cd4=nrow(tcf_cd4) 
num_row_cd4=nrow(cd4_cells)
percentage=(num_row_tcf_cd4/num_row_cd4)*100
print(percentage)
# 49.67451

#Plot immunomap of CD4 and TCF1+
plot(data.frame(cd4_cells_comp[6],-cd4_cells_comp[7]),pch=1,cex=0.0001,col=rgb(0.1,0.4,1,0.2),xlim=xlim,ylim=ylim)
par(new=T)
plot(data.frame(tcf_cd4[6],-tcf_cd4[7]),pch=1,cex=0.0001,col=rgb(0.1,.8,0.2,0.8),xlim=xlim,ylim=ylim)
abline(v=c(5000,10000,15000,20000),h=c(-5000,-10000,-15000,-20000),lwd=1,lty=3,col='grey50')
dev.off()

#aSMA vs CD8
plot(log(cells[,26]+1,10),log(cells[,20]+1,10),pch=19,cex=0.1,col=rgb(0.4,0.4,0.4,0.2),xlim=c(0.0,2.6),ylim=c(0.0,2.6))

cells_asma_neg_tf=apply(cells, 1, function(x) x[20]<2.1)

cells_asma_neg=cells[cells_asma_neg_tf,] 

plot(log(cells_asma_neg[,26]+1,10),log(cells_asma_neg[,20]+1,10),pch=19,cex=0.1,col=rgb(0.4,0.4,0.4,0.2),xlim=c(0.0,2.6),ylim=c(0.0,2.6))

cd8_cells=cells[log(cells[26]+1,10)>1.65 & log(cells[20]+1,10)<1.5,]

pdf("/Users/BTVO/Desktop/K240_asmaxcd8_gateout.pdf")
plot(log(cells_asma_neg[,26]+1,10),log(cells_asma_neg[,20]+1,10),pch=19,cex=0.1,col=rgb(0.4,0.4,0.4,0.2),xlim=c(0.0,2.6),ylim=c(0.0,2.6))
par(new=T)
plot(log(cd8_cells[,26]+1,10),log(cd8_cells[,20]+1,10),pch=19,cex=0.1,col=rgb(1,0.4,0.4,0.2),xlim=c(0.0,2.6),ylim=c(0.0,2.6))
dev.off()

nrow(cd8_cells)
num_row_cd8=nrow(cd8_cells) 
num_row_DAPI=nrow(dapi)
percentage=(num_row_cd8/num_row_DAPI)*100
print(percentage)
# %CD8 = 4.94
# FACS %CD8 = 2.03

#MHC
plot(log(cells[,26]+1,10),log(cells[,32]+1,10),pch=1,cex=0.1,col=rgb(0.4,0.4,0.4,0.4),xlim=c(0,2.3),ylim=c(0.5,2.3))
abline(h=1.7)
plot(log(cells[,44]+1,10),log(cells[,32]+1,10),pch=1,cex=0.1,col=rgb(0.4,0.4,0.4,0.4),xlim=c(0,2.3),ylim=c(0.5,2.3))
abline(h=1.7)

mhc_cells=cells[log(cells[32]+1,10)>1.7 & log(cells[26]+1,10)<1.55,]

plot(log(cells[,26]+1,10),log(cells[,32]+1,10),pch=1,cex=0.1,col=rgb(0.4,0.4,0.4,0.4),xlim=c(0,2.3),ylim=c(0,2.3))
par(new=T)
plot(log(mhc_cells[,26]+1,10),log(mhc_cells[,32]+1,10),pch=1,cex=0.1,col='deeppink3',xlim=c(0,2.3),ylim=c(0,2.3))

plot(data.frame(mhc_cells[6],-mhc_cells[7]),pch=1,cex=0.001,col='deeppink3',xlim=xlim,ylim=ylim)
abline(v=c(2000,6000,10000,14000),h=c(-8000,-10000,-18000),lwd=1,lty=3,col='grey50')
dev.off()

#MHCII+ cells per mm2
xdim2=max(mhc_cells[6])
ydim2=max(mhc_cells[7])

#to convert to mm
xdim2_mm=(xdim2*0.4962)/1000
ydim2_mm=(ydim2*0.4962)/1000

#to get the area
area2=xdim2_mm*xdim2_mm
nMHC=nrow(mhc_cells)
amtMHC=nMHC/area2
amtMHC
#1181.928

#immunomap of mhc+ & tcf1+ cells 
plot(data.frame(mhc_cells[6],-mhc_cells[7]),pch=1,cex=0.001,col='deeppink3',xlim=xlim,ylim=ylim)
par(new=T)
plot(data.frame(tcf_cd8[6],-tcf_cd8[7]),pch=1,cex=0.0001,col=rgb(0.2,.7,0.2,0.8),xlim=xlim,ylim=ylim)
abline(v=c(7000,7500,8000,8500),h=c(-18000,-17500,-17000,-16500),lwd=1,lty=3,col='grey50')
dev.off()

#write MHC+ cells and tcf_cd8 cells to a separate csv for use in python
write.csv(mhc_cells, file = 'D:/Cabo_mIF_QuPath/QuPath/S20-23167-A6/mhc_cells.csv')
head(mhc_cells)

write.csv(tcf_cd8, file = 'D:/Cabo_mIF_QuPath/QuPath/S20-23167-A6/tcf_cd8_cells.csv')
head(tcf_cd8)

write.csv(cd8_cells, file = 'D:/Cabo_mIF_QuPath/QuPath/S20-23167-A6/cd8_cells.csv')
head(cd8_cells)

#####GO TO PYTHON HERE TO CALCULATE DENSITIES. RUN mhc densityquant first and then tcfdens, then come back here#####

####niche quant####
library(raster)
library(tibble)
library(tidyr)
install.packages('raster')
install.packages('tibble')
install.packages('tidyr')

#import density matrices from python for both mhc and tcf. they should always have the same dimensions. 
mhcdensity=read.csv('D:/Cabo_mIF_QuPath/QuPath/S20-23167-A6/mhc_density.csv', header = FALSE)
head(mhcdensity)

cd8density=read.csv('D:/Cabo_mIF_QuPath/QuPath/S20-23167-A6/cd8_density.csv', header = FALSE)
head(cd8density)

tcf_cd8_density=read.csv('D:/Cabo_mIF_QuPath/QuPath/S20-23167-A6/tcf_cd8_density.csv', header = FALSE)
head(tcf_cd8_density)

dim(tcf_cd8_density)
dim(cd8density)
dim(mhcdensity)

#####GO TO PYTHON HERE TO CALCULATE DENSITIES. RUN mhc densityquant first and then tcfdens, then come back here#####

####niche quant####
library(raster)
library(tibble)
library(tidyr)
install.packages('raster')
install.packages('tibble')
install.packages('tidyr')

#import density matrices from python for both mhc and tcf. they should always have the same dimensions. 
mhcdensity=read.csv('F:/Cabo_mIF_QuPath/QuPath/S20-23167-A6/mhc_density.csv', header = FALSE)
head(mhcdensity)

cd8density=read.csv('F:/Cabo_mIF_QuPath/QuPath/S20-23167-A6/cd8_density.csv', header = FALSE)
head(cd8density)

tcf_cd8_density=read.csv('F:/Cabo_mIF_QuPath/QuPath/S20-23167-A6/tcf_cd8_density.csv', header = FALSE)
head(tcf_cd8_density)

dim(tcf_cd8_density)
dim(cd8density)
dim(mhcdensity)

#change rows and columns values to match the dimensions of your density matrix above 
rows = 191
columns = 199
coord=sapply(1:columns, function(i) sapply(1:rows, function(j) paste(j,i,sep = ", ")))
dim(coord)
coord
typeof(coord)
coord=as.matrix(coord)
coordflat=as.vector(coord)
coordflat=as.data.frame(coordflat)
coordflatsep=separate(coordflat, col= 1, into = c("X", "Y"), sep = ",")
coordflatsepintx=as.integer(as.matrix(coordflatsep[,1]))
coordflatsepinty=as.integer(as.matrix(coordflatsep[,2]))
coordflatsepinty

mhcdensity=t(mhcdensity)
mhcdens=as.vector(mhcdensity)

tcf_cd8_density=t(tcf_cd8_density)
tcf_cd8_dens=as.vector(tcf_cd8_density)

#make it into one dataframe
dens = data.frame(X = coordflatsepintx, Y = coordflatsepinty, MHCDens = mhcdens, TCFdens = tcf_cd8_dens)

#write.csv(dens, file='/Users/careyjansen/Documents/Imaging/QuPathQuant/dens.csv')

#this is where we define what a niche is 
#can take some optimizing to see what is appropriate for each tumor/tissue type and for how stringent qupath cutoffs were. 
#i can help you decide this, sometimes trial and error is helpful. 
#brain mets will probably be much lower numbers than those seen here for kidney

nichedens =subset(dens, dens$MHCDens > 16 & dens$TCFdens > 4)
max(dens[,2])
nrow(nichedens)

#plot where the niches are 
#plot(dens[,1], (dens[,2]), cex=1,,pch=19, xlim=c(0,160), ylim=c(0,200))
#par(new=T)
plot(nichedens[,1], (-nichedens[,2]), cex=0.1, col='darkviolet', pch=19, xlim=c(0,180), ylim=c(-200,0))

#calculations
nicheproportion=nrow(nichedens) / nrow(dens)
nicheproportion
nichepercent=nicheproportion*100
nichepercent

hist(dens$MHCDens)

hist(log(dens$MHCDens+1,2))

hist(dens$TCFdens)

hist(log(dens$TCFdens+1,2))