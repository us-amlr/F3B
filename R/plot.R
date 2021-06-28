#' Plot F3B model outputs
#'
#' Plot F3B model outputs
#'
#' @param indata output of...
#' @param pick.me a character vector; one or more of "Y","Z", or "survey"
#'
#' @details
#' Plot the weekly biomass time series, along with survey estimates and effective harvest rates
#'
#' @export
plot.ss<-function(indata, pick.me=c("Y","Z","survey")){

  # set up the y-axes etc.
  axis.id<-NA
  col.id<-NA
  legend.text<-NA
  lty.id<-NA
  pch.id<-NA

  if("X" %in% pick.me){
    axis.id <- c(axis.id,2)
    col.id <- c(col.id,1)
    legend.text <- c(legend.text,"X - subarea")
    lty.id <- c(lty.id,1)
    pch.id <- c(pch.id,NA)
  }

  if("Y" %in% pick.me){
    axis.id <- c(axis.id,3)
    col.id <- c(col.id,2)
    legend.text <- c(legend.text,"Y - mesoscale area")
    lty.id <- c(lty.id,1)
    pch.id <- c(pch.id,NA)
  }

  if("Z" %in% pick.me){
    axis.id <- c(axis.id,4)
    col.id <- c(col.id,3)
    legend.text <- c(legend.text,"Z - fishing area")
    lty.id <- c(lty.id,1)
    pch.id <- c(pch.id,NA)
  }

  if("survey" %in% pick.me){
    col.id <- c(col.id,"blue")
    legend.text <- c(legend.text,"Y+Z - survey")
    lty.id <- c(lty.id,NA)
    pch.id <- c(pch.id,16)
  }

  axis.id<-axis.id[-1]
  col.id<-col.id[-1]
  legend.text<-legend.text[-1]
  lty.id<-lty.id[-1]
  pch.id<-pch.id[-1]

  if("survey" %in% pick.me){
    y.max <- max(indata[,c(axis.id,5)],na.rm=TRUE)*1.01
  } else {
    y.max <- max(indata[,axis.id],na.rm=TRUE)*1.01
  }
  # setup the layout
  layout(matrix(c(1,2),ncol=1,byrow=TRUE))

  # plot the weekly biomass time series
  par(mar=c(2,5,4,2))
  matplot(indata[,1],indata[,axis.id],type="l",xlab="",ylab="Biomass (k tons)", col=col.id, lty=3, lwd=0.3, ylim=c(0,y.max))
  # overlay the time series of biomass on Jan 1 every year
  if("X" %in% pick.me){
    lines(indata[!is.na(indata[,12])&indata[,12]==1&indata[,13]==4,1],indata[!is.na(indata[,12])&indata[,12]==1&indata[,13]==4,2],col=1,lwd=1.5)
  }
  if("Y" %in% pick.me){
    lines(indata[!is.na(indata[,12])&indata[,12]==1&indata[,13]==4,1],indata[!is.na(indata[,12])&indata[,12]==1&indata[,13]==4,3],col=2,lwd=1.5)
  }
  if("Z" %in% pick.me){
    lines(indata[!is.na(indata[,12])&indata[,12]==1&indata[,13]==4,1],indata[!is.na(indata[,12])&indata[,12]==1&indata[,13]==4,4],col=3,lwd=1.5)
  }
  # overlay the survey estimates
  if("survey" %in% pick.me){
    points(indata[,1],indata[,5],pch=16,col="blue")
  }
  # now the legend
  legend("topright", legend.text, col = col.id, lty = lty.id, pch=pch.id, bg="white")

  # plote effective harvest rates (EHR)
  par(mar=c(5,5,2,2))
  matplot(indata[,1],indata[,axis.id+7],type="l",xlab="Month",ylab="EHR", col=col.id, lty=1, lwd=0.3)
  abline(h=0.1,lty=3)
}
