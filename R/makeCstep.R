#' Generate constants in stepwise constraints.
#'
#' @param Numval Number of parameters.
#' @param gamma The constant related to stepwise constraints.
#' @param nsam Number of independent observations making up the means.
#' @param maxp Maximum number of nonzero parameters to consider.
#' @return Stepwise constraints.
#' @export
#' @importFrom stats qt
makeCstep<-function(Numval,gamma,nsam,maxp=NA){
   if(is.na(maxp)) maxp<-Numval
   if(maxp>Numval) print("maxp too big")
   Cstep <- rep(0,Numval)
   for (j in 1:maxp) {
      Cstep[j] <- qt(1-(gamma/2)/(Numval+1-j), df= (Numval+1)*(nsam-1))^2
   }
   if(Numval>maxp){
      for (j in (maxp+1):Numval) {
         Cstep[j] <- 1000000
      }
   }
   return(Cstep)
}
