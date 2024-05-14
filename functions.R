
###
#FUNCTIONS

#M?nima distancia para un h m?nimo
################################################

#nucleo de epanechnikov
nucepan<-function(x)
{
  a=1.5*(1-x*x)
  kepan=a*(abs(x)<1)
  kepan
}

################################################

#me da el Kh
nucleo<-function(punto,xt,ventana)
{
  puntobis=cbind(cos(punto),sin(punto))
  xtbis=cbind(cos(xt),sin(xt))
  prodinter=as.vector((xtbis)%*%t(puntobis))
  arg=(1-prodinter)/(ventana*ventana)
  nucleovec=nucepan(arg)
  nucleovec
}

################################################

# suavizado para epanechnikov (fsombrero)
suavizado<-function(punto,x,ven)
{
  # funcion de la ventana
  intk=12/5
  lambda=intk*sqrt(2)
  cven=1/(lambda*ven)
  wei=nucleo(punto,x,ven)
  suav=cven*mean(wei,na.rm=T)
  suav
}

##############################################
library(CircStats)

argumento<-function(punto,cc,x,conc)
{
  mu=cc[1]
  ven=cc[2]
  arg=(suavizado(punto,x,ven)-dvm(punto,mu,conc))
  argum=arg*arg
  argum
}

##############################################

intargumento<-function(cc,xm,kk)
{
  mm=cc[1]
  hh=cc[2]
  grillas=seq(0,2*pi,length=1000)
  resul=grillas
  for(j in 1:1000)
  {
    resul[j]=argumento(grillas[j],cc,xm,kk)
  }
  res=mean(resul)
  res
}


###############################################

#estimador de m?nima distancia
mindisthminconc<-function(x)
{
  ini<-inicialesrob2(x)
  conc<-ini$conc
  muini=ini$mu
  ini=c(muini,1)
  res=optim(ini,method = "L-BFGS-B",intargumento,kk=conc,xm=x,lower = c(0,0.001), upper = c(2*pi,1.41))
  resultado<-list(mu=res$par[1],ventana=res$par[2],muini=muini,conc=conc)
  return(resultado)
}


################################################
inicialesrob2<-function(muestra)
{
  muini=median(muestra)
  muini=((muini<2*pi)*(muini>0)*muini)+((2*pi+muini)*(muini<0))+((muini-2*pi)*(muini>2*pi))
  mediana<-cbind(cos(muini),sin(muini))
  muestrat<-cbind(cos(muestra),sin(muestra))
  prodinter=median(as.vector((mediana)%*%t(muestrat)))
  #conc1<-log(2)/median(prodinter1)
  #conc<-log(2)/median(prodinter)
  #conc<-qnorm(0.75)/median(2*prodinter)
  #conc<-qnorm(0.75)/median(abs(muestra-muini))
  #conc<-conc*conc
  conc<-c2menos(prodinter)
  resultado<-list(mu=muini,conc=conc)
  return(resultado)
}

library(movMF)
#trabajokofunc2menos <- read.table("~/My Dropbox/daniela/trabajosfuturos/mindistesf/simulacion/circulo/trabajokofunc2menos.txt", quote="\"")
trabajokofunc2menos <- read.table("trabajokofunc2menos.txt", quote="\"")
c2menos<-function(tt)
{
  dimtr<-dim(trabajokofunc2menos)[1]
  dimtr1<-dimtr-1
  aa<-rep(0,dimtr)
  aa[1]<-trabajokofunc2menos[1,2]*(tt<=trabajokofunc2menos[1,1])
  aa[dimtr]<-trabajokofunc2menos[dimtr,2]*(tt>trabajokofunc2menos[dimtr,1])
  for(i in 2:dimtr1)
  {
    aa[i]<-trabajokofunc2menos[i,2]*(tt>trabajokofunc2menos[(i-1),1])*(tt<=trabajokofunc2menos[i,1])
  }
  sum(aa)
}
#Minima distancia para una Von Mises Fisher en dim 2 con k prefijados 

vmd2deteccion<-function(datos,mu)
{
  puntobis=cbind(cos(mu),sin(mu))
  xtbis=cbind(cos(datos),sin(datos))
  distancias=as.vector((xtbis)%*%t(puntobis))
  norma<-2-2*distancias
  resultado<-list(norma=norma)
  return(resultado)
}


cortando<-function(kk,poda)
{
  intfun2<-function(aa,kk,poda)
  {
    cte<-integrate(intfun,-1,1,k=kk)$value
    integrate(intfun,-1,aa,k=kk)$value-(cte*(1-poda))
  }
  y1mpoda<-uniroot(intfun2, lower = -0.9999, upper = 1,kk=kk,poda=poda)$root
  2-2*y1mpoda
  
}


intfun<-function(t,k)
{
  exp(k*t)*(1-t^2)^(-1/2)
}


