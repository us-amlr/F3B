#' Wrapper function
#'
#' Wrapper function for \code{\link{MSEwithF3B}}
#'
#' @param i index of the MC simulation; ignored
#' @param ... passed to \code{\link{MSEwithF3B}}
#'
#' @details
#' Wrapper for \code{\link{MSEwithF3B}}. The purpose of this function is to allow
#' \code{X} from \code{\link[parallel:clusterApply]{parLapply}}
#' (in this case, the index fo the current simulation) to be
#' passed to \code{\link{MSEwithF3B}}, and ignored
#'
#' @export
MSEwithF3B.wrap <- function(i, ...)  MSEwithF3B(...)


#' A wrapper to do many Monte Carlo simulations - the parallel version
#'
#' A wrapper to do many Monte Carlo simulations - the parallel version
#'
#' @param nsims number of Monte Carlo simulations to run. The default is 3
#' @param MC.rates todo
#' @param MC.imports todo
#' @param MC.resets todo
#' @param MC.fishing todo
#' @param suppress logical; if not running in parallel,
#'   should the function print which simulation it is currently running.
#'   Default is FALSE
#' @param num.cores Number of CPUs to over which to distribute computations.
#'   Defaults to \code{NULL}, which uses one fewer than the number of cores
#'   reported by \code{\link[parallel]{detectCores}}
#' @param parallel.seed Should be \code{NULL}; output of \code{\link{parallel::clusterSetRNGStream}}
#'
#' @details
#' todo
#' todo: likely rename inputs so they're the same as MSEwithF3B for consistency
#'
#' @return
#' A three dimensional array...
#'
#' @export
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
