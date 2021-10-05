gbic<-function(Treat,Cntrl,ssq,nsam,Cons,cstr=NULL,verbose=FALSE){
#' Apply generalized BIC model selection.
#'
#' @param Treat  A vector of treatment means.
#' @param Cntrl  A single control mean
#' @param ssq  Sum of squares within treatments and control, divided by the number of observations making up each mean.
#' @param nsam Number of independent observations making up the means.
#' @param Cons Generalized BIC constants.
#' @param cstr  Work constants, regenerated if not given.
#' @param verbose Logical flag indicating whether intermediate results are printed.
#' @return Models selected.
#' @examples
#' #Original single control multiple (3) treatments example of Dunnett (1955).
#' control<-c(55,47,48)
#' treatment<-array(c(55,55,50,64,49,44,64,52,41),c(3,3))
#' Cntrl<-mean(control)
#' Treat<-apply(treatment,1,mean)
#' ssq<-(sum((Cntrl-control)^2)+sum(apply(treatment,2,"-",Treat)^2))/length(control)
#' Cons<-c(1000000,20.5,8,0.00)
#' gbic(Treat,Cntrl,ssq,length(control),Cons)
#' @importFrom stats var
#' @export
   Numval<-length(Treat)
   if(is.null(cstr)) cstr<-makestepmats(Numval)
   Numv2<-2^Numval
   LogD <- log(ssq)
   V2 <- LogD*nsam*(Numval+1) + Cons[1]#(22) from paper
   V1 <- V2
   Pos <- 1
   if(verbose) cat("i,LogD,V2,Pos",1,0,V1,Pos,"\n")
   for (i in 2:Numv2) {
      Tp <- cstr$Onev[i,]*Treat
#     cat(Tp)
      Tpp <- c(Cntrl,Tp[Tp!=0])
      LogD <- log(var(Tpp)*cstr$Onec[i] + ssq)
#     LogD <- log(var(Tpp)*cstr$Onec[i]/ssq + 1)
      V2 <- LogD*nsam*(Numval+1) + Cons[cstr$Onec[i]+1]#(22) from paper
#     V2 <- LogD + (Cons[cstr$Onec[i]+1])/(nsam*(Numval+1))
      if (V2<V1){
         V1 <- V2
         Pos <- i
      }
      if(verbose) cat("Cons",Cons[cstr$Onec[i]+1],"Onec[i]",cstr$Onec[i],"Tp",Tp,"i,LogD,V2,Pos",1,LogD,V2,Pos,"\n")
   }
   gbicout<-1-cstr$Onev[Pos,]
   return(gbicout)
}
