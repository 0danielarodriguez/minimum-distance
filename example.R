source('functions.R')
library(circular)
library(CircStats)

# EXample
#the data
data1 <- rwrappedcauchy(100, mu=circular(0), rho=0.7, control.circular=list(units="degrees"))
alpha<-0.05
malpha<-1-alpha

MDEVM<-mindisthminconc(data1) # minimum distance estimator under von mises distribution
estmu<-MDEVM$mu
estconc<-MDEVM$conc
detection<-vmd2deteccion(data1,estmu)
distance<-detection$norma
corte<-cortando(estconc,malpha)
isout<-1*(distance<=corte)
isout # 0 is an outlier
porcout<-1-mean(isout)
porcout  # outliers proportions in de sample


