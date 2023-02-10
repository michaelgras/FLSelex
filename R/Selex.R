# {{{
# sel 
#
#' converts harvest() and catch.sel() selectivity at age Sa  
#'
#' @param stock Input FLStock object, FLStockLen or FLQuant.
#' @param nyears numbers of last years to compute selectivity
#' @param year option to specify year range, overwrites nyears
#' @param maturity if TRUE sel() extract mat_at_age or mat_at_length
#' @return FLQuant of Selectivity-at-age or Selectivity at length Se  
#' @export
sel <- function (stock, nyears=3, year=NULL, maturity=FALSE){

nyears = min(dims(stock)$maxyear-dims(stock)$minyear+1, nyears)  
  
if(is.null(year)) {
  yr <- (dims(stock)$"maxyear"-nyears+1):dims(stock)$"maxyear"
  } else {
    yr <- year  
}

if(class(stock)=="FLStock"){
  if(!maturity){flq <- harvest(stock)[,ac(yr)]} else {flq <- mat(stock)[,ac(yr)]}
} 

if(class(stock)=="FLStockLen"){
  if(!maturity){flq <- harvest(stock)[,ac(yr)]} else {flq <- mat(stock)[,ac(yr)]}
} 


if (class(stock)=="FLQuant") {
  flq <- stock[,ac(yr)]
}

Se <- apply(flq, 1, mean, na.rm=T) / max(apply(flq, 1, mean, na.rm=T), na.rm=T)

Se@units = "NA"

return(Se)
}


# {{{
# fit.VBGF 
#
#' estimates VBGF parameters based on a set of pairs age and length observed data
#'
#' @param age Observed aged data
#' @param length observed corresponding length
#' @return FLPar with the VBGF parameters  
#' @export

fit.VBGF <- function(age, length){
  dyn.load(dynlib("vonBertalanffy/cpp/VBGF"))
obj <- MakeADFun(data = list(length = length, age = age), parameters = list(Linf = 400, K = .3,T = 0 ,LogSigma = 0), DLL = "VBGF") 

optll <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr)

#VBGFpar <- FLPar(optll$par[1:3], params=c("Linf", "K", "T0")) # if iter is needed  FLPar(opt1$par[1:3], params=c("Linf", "K", "T0"), iter=20)
return(optll)
}