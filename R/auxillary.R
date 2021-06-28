#' Calculate ...
#'
#' Calculate ...
#'
#' @param XYZ.mean todo
#' @param XYZ.CV todo
#' @param YZ.mean todo
#' @param YZ.CV todo
#'
#' @details
#' todo
#'
#' @export
GetCVp <- function(XYZ.mean=19158000,XYZ.CV=0.11,YZ.mean=6853000,YZ.CV=0.16){
  sigma.logXYZ <- sqrt(log((XYZ.CV^2)+1))
  mu.logXYZ <- log(XYZ.mean)-((sigma.logXYZ^2)/2)
  XYZ<-rlnorm(10000,meanlog=mu.logXYZ,sdlog=sigma.logXYZ)

  sigma.logYZ <- sqrt(log((YZ.CV^2)+1))
  mu.logYZ <- log(YZ.mean)-((sigma.logYZ^2)/2)
  YZ<-rlnorm(10000,meanlog=mu.logYZ,sdlog=sigma.logYZ)

  p<-(XYZ-YZ)/XYZ
  cat("CVp =", sqrt(var(p))/mean(p))
}



#' Calculate ...
#'
#' Calculate ...
#'
#' @param x todo
#' @inheritParams TuneF3B
#'
#' @details
#' todo
#'
#' @export
ObjectiveF3B <- function(x,rates,resets,fishing){

  rates.tuning <- rates
  # fill the imports vector with the parameters being estimated
  imports.tuning <- c(amult=exp(x[1]), kmult=exp(x[1])*(exp(x[2])/(1+exp(x[2]))),rmult=exp(x[3]))
  resets.tuning <- resets
  fishing.tuning <- fishing

  tt <- MSEwithF3B(rates=rates.tuning,imports=imports.tuning,resets=resets.tuning,fishing=fishing.tuning)
  sumsqX <- sum((tt[1,2]-tt[(!is.na(tt[,12])&tt[,12]==1&tt[,13]==1&tt[,14]>1),2])^2)
  sumsqY <- sum((tt[1,3]-tt[(!is.na(tt[,12])&tt[,12]==1&tt[,13]==1&tt[,14]>1),3])^2)
  sumsqZ <- sum((tt[1,4]-tt[(!is.na(tt[,12])&tt[,12]==1&tt[,13]==1&tt[,14]>1),4])^2)
  sumsq <- sumsqX + sumsqY + sumsqZ # try to make as many biomasses on 1 Oct = initial biomass

  sumsq

}


#' Calculate ...
#'
#' Calculate ...
#'
#' @param rates todo
#' @param resets todo
#' @param fishing tooo
#'
#' @details
#' todo
#'
#' @return
#'
#' A list...
#'
#' @export
TuneF3B <- function(rates,resets,fishing){

  Rates<-rates
  # turn off process errors
  Resets <- resets
  Resets[c(2,3,5,7)] <- 1e-10
  # turn off fishing and observation errors
  Fishing <- fishing
  Fishing[4] <- 1e-10
  Fishing[7] <- 0


  # do the optimization
  steady <- optim(c(log(0.7),log((0.07/0.7)/(1-(0.07/0.7))),log(1)),ObjectiveF3B,rates=Rates,resets=Resets,fishing=Fishing)

  amult <- exp(steady$par[1])
  kmult <- exp(steady$par[1])*(exp(steady$par[2])/(1+exp(steady$par[2])))
  rmult <- exp(steady$par[3])

  # Annualized Biomass Import is multiplied by 12 for amult and kmult because these imports occur every month
  # whereas rmult only applies once per year
  cat("amult =", amult, "     Annualized Biomass Import =", amult*Resets[1]*12, fill=TRUE)
  cat("kmult =", kmult, "     Annualized Biomass Import =", kmult*Resets[1]*12, fill=TRUE)
  cat("rmult =", rmult, "     Annualized Biomass Import =", rmult*Resets[1]/12, fill=TRUE)

  # plot the outcome
  # time series should be flat
  # first set upt the parameter vectors
  Imports <- c(amult=exp(steady$par[1]), kmult=exp(steady$par[1])*(exp(steady$par[2])/(1+exp(steady$par[2]))), rmult=exp(steady$par[3]))
  outmat <- MSEwithF3B(rates=Rates,imports=Imports,resets=Resets,fishing=Fishing)
  plot.ss(outmat,pick.me=c("Y","Z","survey"))

  list(c(steady,Imports))

}


#' Calculate ...
#'
#' Calculate ...
#'
#' @param indata todo
#'
#' @details
#' todo
#'
#' @export
CalcH <- function(indata){
  denomX <- sum(!is.na(indata[,6])&indata[,6]>0)
  denomY <- sum(!is.na(indata[,7])&indata[,7]>0)
  denomZ <- sum(!is.na(indata[,8])&indata[,8]>0)

  numerX <- sum(!is.na(indata[,6])&indata[,6]>0&indata[,9]>=0.1)
  numerY <- sum(!is.na(indata[,7])&indata[,7]>0&indata[,10]>=0.1)
  numerZ <- sum(!is.na(indata[,8])&indata[,8]>0&indata[,11]>=0.1)

  c(HX = numerX/denomX, HY = numerY/denomY, HZ = numerZ/denomZ)
}


#' Calculate ...
#'
#' Calculate ...
#'
#' @inheritParams CalcH
#'
#' @details
#' todo
#'
#' @export
PropSurveyOofM <- function(indata){
  tt<-round(indata[!is.na(indata[,5]),5],0)
  tt<-nchar(as.character(tt))
  ttt<-0
  for(i in 2:length(tt)){
    if(tt[i]!=tt[i-1]){ttt<-ttt+1}
  }
  ttt/(length(tt)-1)
}
