# MSE with F3B (Flux thru 3 Biomass Pools)
# GEORGE WATTERS

# last edited - late Sept 2020

library(abind)
library(deSolve)
library(parallel)


# Parallel FUNCTIONS ##########################################################

# From swfscMisc: https://github.com/EricArcher/swfscMisc/blob/master/R/setupClusters.R
setupClusters <- function (num.cores = 1, max.cores = NULL) {
  if (is.null(num.cores)) num.cores <- parallel::detectCores() - 1
  if (is.null(max.cores)) max.cores <- parallel::detectCores() - 1
  if (is.na(max.cores)) max.cores <- 1
  if (max.cores < 1) max.cores <- 1
  num.cores <- min(num.cores, max.cores)
  
  cl.func <- ifelse(.Platform$OS.type == "windows", 
                    parallel::makePSOCKcluster, parallel::makeForkCluster)
  if (num.cores > 1) cl.func(num.cores) else NULL
}


MSEwithF3B.wrap <- function(i, ...)  MSEwithF3B(...)


# A wrapper to do many Monte Carlo simulations - the parallel version
MC.MSEwithF3B.par <- function(nsims=3,MC.rates,MC.imports,MC.resets,MC.fishing,
                              suppress=FALSE, num.cores=NULL, parallel.seed=NULL) {
  
  # # set up dimensions of output array
  # outarray <- array(data=NA,dim=c(49*MC.fishing[1]*MC.fishing[2],14,nsims))
  
  stopifnot(require(abind), require(deSolve), require(parallel))
  
  cl <- setupClusters(num.cores)
  outarray.list <- tryCatch({
    if (is.null(cl)) { # Don't parallelize if num.cores == 1
      lapply(1:nsims, function(i) {
        if(!suppress) cat("Running Sim",i,"of",nsims,fill=TRUE)
        MSEwithF3B(rates=MC.rates,imports=MC.imports,resets=MC.resets,fishing=MC.fishing)
      })
      
    } else { # Run lapply using parLapply
      parallel::clusterExport(
        cl = cl, varlist = c("MC.rates", "MC.imports", "MC.resets", "MC.fishing"),
        envir = environment()
      )
      parallel::clusterEvalQ(cl, {
        library(deSolve)
        source("../F3B/MSEwithF3Bv3_functions.R", echo = FALSE)
      })
      
      # Set parallel seed and use parLapply (not parLapplyLB) for reproducibility
      # See: https://pat-s.me/post/reproducibility-when-going-parallel/
      #https://s3.amazonaws.com/assets.datacamp.com/production/course_6245/slides/chapter4.pdf
      if (!is.null(parallel.seed)) parallel::clusterSetRNGStream(cl = cl, parallel.seed)
      parallel::parLapply( 
        cl, 1:nsims, MSEwithF3B.wrap, 
        rates=MC.rates, imports=MC.imports, resets=MC.resets, fishing=MC.fishing
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

