#' Update counts of errors to track behavior of various tests, called within table1line
#'
#' @param errors A matrix of errors
#' @param FR output from gbic or steptest
#' @param Nmu number of comparisons.
#' @return the updated error matrix
#' @export
updateerrors<-function(errors,FR,Nmu){
      errors[-(1:2)] <- errors[-(1:2)]+FR
      T2 <- sum(Nmu*FR)
      if (T2 != 0 ) { 
         errors["Fwer"] <- errors["Fwer"]+1
         errors["Fdr"] <- errors["Fdr"] + T2/max(sum(FR),1)
      }
      return(errors)
}
#tab1<-table1(printme=TRUE,setseed=TRUE)
#tab2<-table1(printme=TRUE,setseed=TRUE,BIC=10)
