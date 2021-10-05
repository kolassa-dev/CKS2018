#' Generate a table like Table 1 of Cohen, Kolassa, and Sackrowitz
#'
#' @param Cons A vector of constants for the generalized BIC.
#' @param nsam A sample size for each of the means investigated.
#' @param gamma The constant related to stepwise constraints.
#' @param Nrep Number of Monte Carlo samples for size and power calculations.
#' @param printme Logical flag indicating whether intermediate results are printed.
#' @param setseed Logical flag indicating whether random seed will be.
#' @param BIC Integer. if not NA, number of categories to apply traditonal BIC.
#' @param maxp Integer. if not NA, maximum numbers of parameters in the model.
#' @param irange range for means for table construction.
#'
#' @return A table of the format of Table 1.
#' @examples
#' #Call to table1 regenerates Table 1.
#'\donttest{
#' tab1<-round(table1(setseed=TRUE,maxp=2)[c(6:4,13,12,10,14:16,19:24),],2)
#' tab2a<-table1(printme=TRUE,setseed=TRUE,BIC=10)
#' tab2b<-table1(printme=TRUE,setseed=TRUE,BIC=10,maxp=2)
#' tab2<-cbind(tab2a[,1:4],tab2b[,3:4])[c(5,4,12:10,14:17,19:21,23),]
#' }
#' #A version of Table 1 based on fewer samples.
#' tab1<-table1(printme=TRUE,setseed=TRUE,Nrep=2)
#' @export
table1<-function(Cons=c(rep(1000000,8),20.5,8,0.00),nsam=14,gamma=0.05, Nrep=5000,printme=FALSE,setseed=FALSE,BIC=NA,maxp=NA, irange=NULL){
   if(is.null(irange)) irange<-c(-2,4)
   if(setseed) set.seed(899884357)
   if(!is.na(BIC)){
      bicname<-"BIC"
      Cons<-rev((0:BIC)*log(nsam*(BIC+1)))
   }else{
      bicname<-"GBIC"
   }
   if(!is.na(maxp)){
      Cons[seq(length(Cons)-1-maxp)]<-100000
   }
   Numval<- length(Cons)-1
   Cons[Numval + 1] <- 0.00
   Cstep<-makeCstep(Numval,gamma,nsam,maxp)
   cat("Cstep",Cstep,"\n",file="temp")
#  print(Cstep)
#  print(CNstep)
   out<-NULL
   for (is in irange[1]:irange[2]) {
      for (ir in is:irange[2]) {
         Mu<-rep(0,Numval)
         Mu[Numval-1] <- .75*is
         Mu[Numval  ] <- .75*ir
         outline<-table1line(Mu,Cons,Cstep,nsam,Nrep,printme)
         out<-rbind(out,outline)
      }
   }
   if(printme) print (c("penalties =   ",Cons),quote=FALSE)
#  print (s)
#  print(Treat)
#  print (V1)
   dimnames(out)<-list(rep("",dim(out)[1]),
      c("Stand. Mean 1","Stand. Mean 2",
      paste(bicname,c("fwer","fdr","Av Pow"),sep=" "),
      paste("SD",c("fwer","fdr","Av Pow"),sep=" ")))
   return(out)
}
