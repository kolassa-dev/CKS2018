#' Generate a single line of a table like Table 1 of Cohen, Kolassa, and Sackrowitz
#'
#' @param Mu A vector of sample means for simulation.
#' @param Cons A vector of constants for the generalized BIC.
#' @param Cstep A vector of constants for stepwise selection.
#' @param nsam A sample size for each of the means investigated.
#' @param Nrep Number of Monte Carlo samples for size and power calculations.
#' @param printme Logical flag indicating whether intermediate results are printed to the screen.
#' @return A line of a table of the format of Table 1.
#' @importFrom stats rnorm
#' @importFrom stats rchisq
#' @export
table1line<-function(Mu,Cons,Cstep,nsam,Nrep=10000,printme=FALSE){
   print("Entering table1line")
   print(Nrep)
   Numval <- length(Mu)
# Indicator for which null hypotheses are true.
   Nmu=(Mu==0)+0
# Both assessments of error will record Fwer and Fdr in first two positions, followed by proportion of time remainder are rejected.
   gbicerrors<-rep(0,Numval+2)
   names(gbicerrors)<-c("Fwer","Fdr",paste("p",1:Numval,sep=""))
   steperrors<-gbicerrors
   cstr<-makestepmats(Numval)
   for (Rep in 1:Nrep) {
      Cntrl <- rnorm(1)/sqrt(nsam)
      Treat <- rnorm(Numval)/sqrt(nsam)+Mu
      ssq<- rchisq(1,df= (Numval+1)*(nsam-1))/nsam
      FRgb<-gbic(Treat,Cntrl,ssq,nsam,Cons,cstr=cstr)
      gbicerrors<-updateerrors(gbicerrors,FRgb,Nmu)
# Now to get the step down procedure
      FRej<-steptest(Treat,Cntrl,ssq,nsam,Cstep)
      steperrors<-updateerrors(steperrors,FRej,Nmu)
   }
   gbicerrors <- gbicerrors/Nrep 
   steperrors <- steperrors/Nrep
   Powbic <- gbicerrors[-(1:2)]
   Powstp <- steperrors[-(1:2)]
   Avpowbic <- sum(Powbic*(1-Nmu))/sum(1-Nmu)
   Avpowstp <- sum(Powstp*(1-Nmu))/sum(1-Nmu)
   if(printme) print(c("-                   "),quote=FALSE)
   if(printme) print(c("-                   "),quote=FALSE)
   if(printme) print(c("Treatment means     ",Mu[9],Mu[10]),quote=FALSE)
   if(printme) print(c("BIC  FWER and FDR =    ",gbicerrors[1:2]),quote=FALSE)
   if(printme) print(c("Step FWER and FDR =    ",steperrors[1:2]),quote=FALSE)
#  print(c("BIC  power =    ",Powbic),quote=FALSE)
#  print(c("Step Power =    ",Powstp),quote=FALSE)
   if(printme) print(c("Avg. BIC  power =    ",Avpowbic),quote=FALSE)
   if(printme) print(c("Avg. Step  power =    ",Avpowstp),quote=FALSE)
   out<-c(Mu[(Numval-1):Numval],gbicerrors[1:2],Avpowbic,steperrors[1:2],Avpowstp)
   return(out)
}
