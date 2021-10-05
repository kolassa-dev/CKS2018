makestepmats<-function(Numval){
   Numv2<-2^Numval
   Onev <- array(as.integer(0),c(Numv2,Numval))
   for (i in 0:Numv2-1) {
      x <- i
      for (j in Numval:1) {
          y <- as.integer(x/2)
          One <- x-2*y
          Onev[i+1,j] <- One
          x <- y
      }
   }
   return(list(Onev=Onev,Onec=apply(Onev,1,sum)))
}
