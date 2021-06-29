#' A model for within-year biomass dynamics in three nested biomass pools - Z within Y within X
#'
#' A model for within-year biomass dynamics in three nested biomass pools - Z within Y within X
#'
#' @param t todo
#' @param state todo
#' @param parameters todo
#' @param survey todo
#' @param surveytime todo
#' @param catfuncX todo
#' @param catfuncY todo
#' @param catfuncZ todo
#'
#' @details
#' Run the model...
#'
#' @returns A list of ...
#'
#' @export
F3B <- function(t, state, parameters, survey, surveytime,
                catfuncX, catfuncY, catfuncZ){
  with(as.list(c(state,parameters)),{

    # get the interpolated catches
    CX <- catfuncX(t)
    CY <- catfuncY(t)
    CZ <- catfuncZ(t)

    # reduce catches if biomasses too small
    if(CX > X){CX <- X/2}
    if(CY > Y){CY <- Y/2}
    if(CZ > X){CZ <- Z/2}

    # rates of change
    dX <- a + ((b-c-d)*X) + (f*Y) - CX
    dY <- (d*X) + ((e-f-g)*Y) + (i*Z) - CY
    dZ <- (g*Y) + ((h-i)*Z) + k - CZ

    if(t==surveytime){
      junk <- survey
    } else {
      junk <- NA
    }

    # return the results
    list(c(dX, dY, dZ), c(survey=junk, CX=CX, CY=CY, CZ=CZ))
  })
}



#' An event function that tries to ensure biomasses >= 0
#'
#' An event function that tries to ensure biomasses >= 0
#'
#' @inheritParams F3B
#'
#' @export
PositiveF3B <- function(t, state, parameters, survey, surveytime,
                        catfuncX, catfuncY, catfuncZ){
  with(as.list(state),{
    if(X < 0){X <- 0}
    if(Y < 0){Y <- 0}
    if(Z < 0){Z <- 0}

    junk <- survey
    junk <- surveytime
    junk <- catfuncX(t)
    junk <- catfuncY(t)
    junk <- catfuncZ(t)

    return(state)
  })
}



#' Main function for running F3B
#'
#' Run F3B model, once or in Monte Carlo wrapper
#'
#' @param rates todo
#' @param imports todo
#' @param resets todo
#' @param fishing todo
#'
#' @details
#' todo
#'
#' See \url{https://pat-s.me/post/reproducibility-when-going-parallel/#reproducible-1}
#' for more about the reproducibility method used when running MC simulations in parallel.
#' Note that you will get different outputs when using the same \code{seed} value
#' when running the MC simulations in parallel vs not,
#' i.e., when \code{num.cores} is 1 vs otherwise.
#'
#' @return
#' \code{MSEwithF3B} returns a two-dimensional array
#'
#' \code{Mc.MSEwithF3B} returns a three-dimensional array;
#' the third dimension is the MC simulation runs
#'
#' @export
MSEwithF3B <- function(rates,imports,resets,fishing){

  # set the time sequence within years
  # nominal 4 weeks per month with time 1 = Oct 1; time 4 = January 1; time 7 = April 1; time 10 = July 1
  F3Btime <- seq(from=1,to=12.75,by=0.25)

  with(as.list(c(rates,imports,resets,fishing)),{

    nfishingweeks <- length(seq(from=startfishingtime,to=stopfishingtime,by=0.25))-1

    # be sure results from relative risk assessment sum to 1
    sumalpha <- sum(alphaX,alphaY,alphaZ)
    alphaX <- alphaX/sumalpha
    alphaY <- alphaY/sumalpha
    alphaZ <- alphaZ/sumalpha

    # wrap within-year dynamics inside some for() loops to simulate multiple years
    firstsurvey <- TRUE
    for(i in 1:nsurveys){
      for(j in 1:yrstillnextsurvey){

        if(firstsurvey){
          firstsurvey <- FALSE

          # set total biomass in subarea (XYZ)
          vv <- log(initXYZcv^2 + 1)
          mm <- log(initXYZmean)-(vv^2)/2
          XYZ <- rlnorm(1,mm,sqrt(vv))

          # set inputs from the Bellingshausen Sea (A) and Weddell Sea (K)
          # parameterized to be fractions of initial total biomass
          A <- (XYZ/12)*amult
          K <- (XYZ/12)*kmult

          # set mean recruitment
          # parameterized to be a fraction of initial total biomass
          meanR <- (XYZ/12)*rmult

        } else {

          # total biomass at start of year is function of biomass at end of previous year and random perturbation
          # need to subtract 1 from row index of outmat because added NAs to end
          Tend <- sum(outmat[dim(outmat)[1]-1,2:4])
          mm <- log(meanR)-(prosd^2)/2
          XYZ <- Tend + rlnorm(1,mm,prosd)

          proeps <- log(XYZ/Tend)

          # set inputs from the Bellingshausen Sea (A) and Weddell Sea (K)
          A <- (Tend*exp(proeps)/12)*amult
          K <- (Tend*exp(-proeps)/12)*kmult

        }

        # redistribute XYZ among the boxes and initialize state variables for the year
        mm<-pXinXYZmean
        vv<-(pXinXYZmean*pXinXYZcv)^2
        aa <- mm*(((mm*(1-mm))/vv)-1)
        bb <- (1-mm)*(((mm*(1-mm))/vv)-1)
        pXinXYZ<-rbeta(1,shape1=aa,shape2=bb)

        mm<-pYinYZmean
        vv<-(pYinYZmean*pYinYZcv)^2
        aa <- mm*(((mm*(1-mm))/vv)-1)
        bb <- (1-mm)*(((mm*(1-mm))/vv)-1)
        pYinYZ<-rbeta(1,shape1=aa,shape2=bb)

        X1 <- XYZ * pXinXYZ
        Y1 <- XYZ*(1-pXinXYZ)*pYinYZ
        Z1 <- XYZ*(1-pXinXYZ)*(1-pYinYZ)

        # set up for the ODEs
        F3Bstate <- c(X=X1, Y=Y1, Z=Z1)
        F3Bparameters <- c(a=A,rates,k=K)


        if(j==1){

          # quick run of within-year dynamics to estimate expected value of survey biomass
          # no catches removed until survey
          tt <- seq(from=1,to=surveytime,by=0.25)
          catchX <- approxfun(x=tt,y=rep(0,length(tt)),method="constant",rule=2)
          catchY <- approxfun(x=tt,y=rep(0,length(tt)),method="constant",rule=2)
          catchZ <- approxfun(x=tt,y=rep(0,length(tt)),method="constant",rule=2)
          ttt <- deSolve::ode(
            y=F3Bstate, times=tt, func=F3B, parms=F3Bparameters, method="ode45",
            events=list(func=PositiveF3B, time=tt), rtol=1e-12, atol=0,
            survey=NA, surveytime=surveytime,
            catfuncX=catchX, catfuncY=catchY, catfuncZ=catchZ)
          surveyYZmean <- sum(ttt[dim(ttt)[1],3:4])

          # survey estimate of biomass contaminated by observation error
          vv <- log(surveycv^2 + 1)
          mm <- log(surveyYZmean)-(vv^2)/2
          surveyYZ <- rlnorm(1,mm,sqrt(vv))

          #obseps <- rnorm(1,mean=0,sd=obssd)-(0.5*obssd^2)
          #surveyYZ <- surveyYZmean*exp(obseps)

          # expand survey biomass to subarea, compute catch limit, and apply risk assessment
          surveyXYZ <- surveyYZ/(1-pXinXYZmean)
          CLXYZ <- GYMgamma * surveyXYZ
          CLX <- CLXYZ * alphaX
          CLY <- CLXYZ * alphaY
          CLZ <- CLXYZ * alphaZ

          # make data sets used to establish interpolation functions that force catches
          # assume catches evenly distributed across fishing period
          weeklycatchX <- CLX/nfishingweeks
          weeklycatchY <- CLY/nfishingweeks
          weeklycatchZ <- CLZ/nfishingweeks

          catchX <- approxfun(x=F3Btime,y=ifelse(F3Btime<startfishingtime,0,ifelse(F3Btime>=stopfishingtime,0,weeklycatchX)),method="constant",rule=2)
          catchY <- approxfun(x=F3Btime,y=ifelse(F3Btime<startfishingtime,0,ifelse(F3Btime>=stopfishingtime,0,weeklycatchY)),method="constant",rule=2)
          catchZ <- approxfun(x=F3Btime,y=ifelse(F3Btime<startfishingtime,0,ifelse(F3Btime>=stopfishingtime,0,weeklycatchZ)),method="constant",rule=2)

        } else {

          surveyYZ <- NA

        }


        if(i==1 & j==1){
          outmat <- deSolve::ode(
            y=F3Bstate, times=F3Btime, func=F3B, parms=F3Bparameters, method="ode45",
            events=list(func=PositiveF3B, time=F3Btime), rtol=1e-12, atol=0,
            survey=surveyYZ, surveytime=surveytime,
            catfuncX=catchX, catfuncY=catchY, catfuncZ=catchZ)
        } else {
          tt <- deSolve::ode(
            y=F3Bstate,times=F3Btime, func=F3B,parms=F3Bparameters, method="ode45",
            events=list(func=PositiveF3B, time=F3Btime), rtol=1e-12, atol=0,
            survey=surveyYZ, surveytime=surveytime,
            catfuncX=catchX, catfuncY=catchY, catfuncZ=catchZ)
          # make the time stamp cumulative since time 1
          tt[,1] <- seq(from=outmat[dim(outmat)[1]-1,1]+0.25,length.out=48,by=0.25)
          outmat <- rbind(outmat,tt)
        }
        # add NAs at the end of each year so that time-series plots are a little nicer
        outmat <- rbind(outmat,rep(NA,8))
      }
    }

    # add effective harvest rate on X, Y and Z to the output
    EHR <- outmat[,6:8]/outmat[,2:4]
    dimnames(EHR)[[2]]<-c("EHRX","EHRY","EHRZ")

    outmat <- cbind(outmat,EHR)

    # add week, month, and year time stamps to output
    TS <- cbind(rep(c(rep(1:4,12),NA),nsurveys*yrstillnextsurvey),
                rep(c(rep(1:12,each=4),NA),nsurveys*yrstillnextsurvey),
                rep(1:(nsurveys*yrstillnextsurvey),each=(48+1)))
    TS[,3]<-ifelse(is.na(TS[,2]),NA,TS[,3])
    dimnames(TS)[[2]]<-c("week","month","year")

    outmat <- cbind(outmat,TS)

    return(outmat)

  })
}


#' @rdname MSEwithF3B
#'
#' @param nsims number of Monte Carlo simulations to run. The default is 3
#' @param suppress logical; if not running in parallel,
#'   should the function print which simulation it is currently running.
#'   Default is FALSE. Ignored if \code{num.cores} is not 1
#' @param num.cores Number of CPUs to over which to distribute computations.
#'   Defaults to \code{NULL}, which uses one fewer than the number of cores
#'   reported by \code{\link[parallel]{detectCores}}
#' @param seed numeric; if not NULL, then either passed to \code{seed} of
#'    \code{\link[base:Random]{set.seed}} or to \code{iseed} argument of
#'    \code{\link[parallel]{clusterSetRNGStream}},
#'    depending on \code{num.cores}. See Details for more information.
#'
#' @examples
#' MC.MSEwithF3B(
#'   nsims = 3,
#'   rates = c(b=-0.06/12, c=0.27/12, d=0.27/12, e=-0.06/12, f=0.27/12,
#'                g=0.27/12, h=-0.12/12, i=0.48/12),
#'   imports = c(amult = 0.09004751, kmult = 0.01793353, rmult = 1.340394),
#'   resets = c(initXYZmean=19158, initXYZcv=0.11, prosd=0.7,
#'                 pXinXYZmean=0.64, pXinXYZcv=0.11, pYinYZmean=0.63, pYinYZcv=0.15),
#'   fishing = c(nsurveys=10, yrstillnextsurvey=3, surveytime=4, surveycv=1e-10,
#'                  startfishingtime=7, stopfishingtime=11, GYMgamma=0,
#'                  alphaX=0.741, alphaY=0.174, alphaZ=0.085),
#'   num.cores = 1, seed = 1234
#' )
#'
#' @export
MC.MSEwithF3B <- function(nsims=3, rates, imports, resets, fishing,
                              suppress=FALSE, num.cores=NULL, seed=NULL) {

  # # set up dimensions of output array
  # outarray <- array(data=NA,dim=c(49*MC.fishing[1]*MC.fishing[2],14,nsims))

  # stopifnot(require(abind), require(deSolve), require(parallel))

  cl <- setupClusters(num.cores)
  outarray.list <- tryCatch({
    if (is.null(cl)) { # Don't parallelize if num.cores == 1
      if (!is.null(seed)) {
        on.exit(set.seed(NULL))
        set.seed(seed)
      }
      lapply(1:nsims, function(i) {
        if(!suppress) cat("Running Sim",i,"of",nsims,fill=TRUE)
        MSEwithF3B(rates=rates, imports=imports, resets=resets, fishing=fishing)
      })

    } else { # Run lapply using parLapply
      parallel::clusterExport(
        cl = cl, varlist = c("rates", "imports", "resets", "fishing"),
        envir = environment()
      )
      parallel::clusterEvalQ(cl, {
        library(deSolve)
      })

      # Set parallel seed and use parLapply (not parLapplyLB) for reproducibility
      # See: https://pat-s.me/post/reproducibility-when-going-parallel/
      #https://s3.amazonaws.com/assets.datacamp.com/production/course_6245/slides/chapter4.pdf
      if (!is.null(seed)) parallel::clusterSetRNGStream(cl = cl, seed)
      # parallel::parLapply(
      #   cl, 1:nsims, MSEwithF3B.wrap,
      #   rates=rates, imports=imports, resets=resets, fishing=fishing
      # )
      parallel::parLapply(
        cl, 1:nsims, function(i, ...) {MSEwithF3B(...)},
        rates=rates, imports=imports, resets=resets, fishing=fishing
      )
    }
  }, finally = if(!is.null(cl)) parallel::stopCluster(cl) else NULL)


  # for(i in 1:nsims){
  #   if(!suppress){cat("Running Sim",i,"of",nsims,fill=TRUE)}
  #   outarray[,,i]<-MSEwithF3B(rates=MC.rates,imports=MC.imports,resets=MC.resets,
  #                             fishing=MC.fishing)
  # }
  # outarray

  abind::abind(outarray.list, along = 3)
}
