#' Apply stepwise model selection.
#'
#' @param Treat A vector of treatment means.
#' @param Cntrl A single control mean
#' @param ssq Sum of squares within treatments and control, divided by the number of observations making up each mean.
#' @param nsam Number of independent observations making up the means.
#' @param Cstep Stepwise constants.
#' @return Models selected.
#' @examples
#' #Original single control multiple (3) treatments example of Dunnett (1955).
#' control<-c(55,47,48)
#' treatment<-array(c(55,55,50,64,49,44,64,52,41),c(3,3))
#' Cntrl<-mean(control)
#' Treat<-apply(treatment,1,mean)
#' ssq<-(sum((Cntrl-control)^2)+sum(apply(treatment,2,"-",Treat)^2))/length(control)
#' Cons<-c(1000000,20.5,8,0.00)
#' steptest(Treat,Cntrl,ssq,length(control),makeCstep(length(Treat),0.05,length(control),maxp=2))
#' @export
steptest<-function(Treat,Cntrl,ssq,nsam,Cstep){
   Numval<-length(Treat)
   Stat <- (nsam-1)*(Numval+1)*(Treat-Cntrl)^2/(2*ssq)
   RR <- sort(Stat,decreasing=T,index.return=T)
   FRej <- rep(0,times=Numval)
   for(j in 1:Numval) {
      if(RR$x[j] > Cstep[j]){
         FRej[RR$ix[j]] <- 1
      }else break
   }
   return(FRej)
}
