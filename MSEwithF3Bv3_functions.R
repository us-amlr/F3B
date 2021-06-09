# MSE with F3B (Flux thru 3 Biomass Pools)
# GEORGE WATTERS

# last edited - late Sept 2020

library(deSolve)



# MAIN FUNCTIONS ##############################################################

# A model for within-year biomass dynamics in three nested biomass pools - Z within Y within X
F3B <- function(t,state,parameters,survey,surveytime,catfuncX,catfuncY,catfuncZ){
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

# An event function that tries to ensure biomasses >= 0
PositiveF3B <- function(t,state,parameters,survey,surveytime,catfuncX,catfuncY,catfuncZ){
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

# Piece everything together into the main function
MSEwithF3B <- function(rates,imports,resets,fishing){
  
  # set the time sequence within years
  F3Btime <- seq(from=1,to=12.75,by=0.25)     # nominal 4 weeks per month with time 1 = Oct 1; time 4 = January 1; time 7 = April 1; time 10 = July 1
  
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
          ttt <- ode(y=F3Bstate,times=tt,func=F3B,parms=F3Bparameters,method="ode45",events=list(func=PositiveF3B,time=tt),rtol=1e-12,atol=0,survey=NA,surveytime=surveytime,catfuncX=catchX,catfuncY=catchY,catfuncZ=catchZ)
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
          outmat <- ode(y=F3Bstate,times=F3Btime,func=F3B,parms=F3Bparameters,method="ode45",events=list(func=PositiveF3B,time=F3Btime),rtol=1e-12,atol=0,survey=surveyYZ,surveytime=surveytime,catfuncX=catchX,catfuncY=catchY,catfuncZ=catchZ)
        } else {
          tt <- ode(y=F3Bstate,times=F3Btime,func=F3B,parms=F3Bparameters,method="ode45",events=list(func=PositiveF3B,time=F3Btime),rtol=1e-12,atol=0,survey=surveyYZ,surveytime=surveytime,catfuncX=catchX,catfuncY=catchY,catfuncZ=catchZ)
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

# A wrapper to do many Monte Carlo simulations
MC.MSEwithF3B <- function(nsims=3,MC.rates,MC.imports,MC.resets,MC.fishing,suppress=FALSE){
  
  # set up dimensions of output array
  outarray <- array(data=NA,dim=c(49*MC.fishing[1]*MC.fishing[2],14,nsims))
  
  for(i in 1:nsims){
    if(!suppress){cat("Running Sim",i,"of",nsims,fill=TRUE)}
    outarray[,,i]<-MSEwithF3B(rates=MC.rates,imports=MC.imports,resets=MC.resets,fishing=MC.fishing)
  }
  
  outarray
}



# PLOTTING ###############################################################

plot.single.sim<-function(indata){
  
  layout(matrix(c(1,2),ncol=1,byrow=TRUE))
  
  par(mar=c(2,5,4,2))
  # plot the biomass time series for the mesoscale and fishing areas
  matplot(indata[,1],indata[,3:4],type="l",xlab="",ylab="Biomass (k tons)", col=2:3, lty=3, lwd=0.3)
  # plot the time series of biomass on Jan 1 every year
  lines(indata[!is.na(indata[,12])&indata[,12]==1&indata[,13]==4,1],indata[!is.na(indata[,12])&indata[,12]==1&indata[,13]==4,3],col=2,lwd=1.5)
  lines(indata[!is.na(indata[,12])&indata[,12]==1&indata[,13]==4,1],indata[!is.na(indata[,12])&indata[,12]==1&indata[,13]==4,4],col=3,lwd=1.5)
  # add the survey estimates
  points(indata[,1],indata[,5],pch=16,col="blue")
  legend("topright", c("Y - mesoscale area","Z - fishing area","Y+Z - survey"), col = c(2:3,"blue"), lty = c(1,1,NA), pch=c(NA,NA,16), bg="white")
  
  par(mar=c(5,5,2,2))
  # plot the effective harvest rates for the mesoscale and fishing areas
  matplot(indata[,1],indata[,10:11],type="l",xlab="Month",ylab="EHR", col=2:3, lty=1, lwd=0.3)
  abline(h=0.1,lty=3)
} # this function is deprecated -- see plot.ss()

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



# AUXILIARY FUNCTIONS ####################################################

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

GetRAconsumption <- function(){
  
  # Unless otherwise noted, code and data files used here provided by Vicky Warwick-Evans to WG-EMM e-group or to GW via email
  
  # NOTE - there are several errors produced that I have not diagnosed
  # Also, I'm not clear if the results are consistent with Vicky's paper
  # I will use them for now, but probably need to revisit this
  
  library(raster)
  library(rgdal)
  # GW - this also seems to be needed
  library(rgeos)
  
  #make a shapefile of the study area: the operational footprint for the krill fishery within CCAMLR Subarea 48.1 over the last 5 years 
  
  coords = matrix(c(	-52,	-63,
                     -50,-60,
                     -51,	-60,
                     -52,	-60,
                     -53,	-60,
                     -54,	-60,
                     -55,	-60,
                     -56,	-60,
                     -57,	-60,
                     -58,	-60,
                     -59,	-60,
                     -60,	-60,
                     -69,	-65,
                     -69,	-66,
                     -66,	-67,
                     -62,	-65,
                     -59,	-64,
                     -52,	-63), 
                  ncol = 2, byrow = TRUE)
  P2 = Polygon(coords)
  Ps1 = SpatialPolygons(list(Polygons(list(P2), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84"))
  
  #load in land shapefile
  # GW - edited the path for use on my local machine
  coast<-readOGR("ADD_Coastline_low_res_polygon.shp")
  
  #first calculate parameters for summer (October to March)
  #read in rasters of consumption for central place foragers (cpf) and for pelagic species (pel).
  # GW - edited paths for use on my local machine
  cpfmask1<-raster( "cpf_summerkgkrillperday.tif") #these data are in kg of krill per day 
  pelmask1<-raster("pelagics_summerkgkrillperday.tif")
  
  #read in SSMU shapefile
  # GW - edited the path for use on my local machine but commented out in favor of code in next paragraph
  #ssmu<-readOGR("ssmus.shp")
  
  # GW - this code bit from email sent to myself from Vicky on 17 Sept 2020
  # GW - edited the path for use on my local machine
  ssmu<-readOGR("ssmus.shp")
  e<-extent(-3000000,-2000000,1000000,2000000)
  ssm<-crop(ssmu,e)
  #plot(ssm,axes=T)
  ssmu<-spTransform(ssm,CRS=CRS("+proj=longlat +ellps=WGS84"))
  # GW - commented out following two line because they already appear below
  #plot(ssmu,axes=T)
  #ssmu<-crop(ssmu,Ps1)
  
  plot(ssmu,axes=T)
  
  #crop the ssmu shapefile to include only the areas in the study area
  ssmu<-crop(ssmu,Ps1)
  #crop the land shapefile to include only study area
  # GW - the following line is causing an error in rgeos
  # "TopologyException: Input geom 0 is invalid: Ring Self-intersection at or near point ..."
  # not sure if this is actually needed so am commenting it out for now
  #land<-crop(coast,Ps1)
  
  #create a data frame from data held in the ssmu shapefile
  da<-ssmu@data
  
  #now we need to extract the values of krill consumption within each ssmu 
  #first for cpf 
  lis<-extract(cpfmask1,ssmu[ssmu$ssmucode=="APW",])#first extract the value of all the raster cells within the ssmu
  li<-Reduce(`+`, lis)#extract the value for each cell from the list
  apw<-sum(na.omit(li)) #add the value for each cell together to get the total consumption within the ssmu
  lis<-extract(cpfmask1,ssmu[ssmu$ssmucode=="APDPW",])
  li<-Reduce(`+`, lis)
  apdpw<-sum(na.omit(li))
  lis<-extract(cpfmask1,ssmu[ssmu$ssmucode=="APBSW",])
  li<-Reduce(`+`, lis)
  apbsw<-sum(na.omit(li))
  lis<-extract(cpfmask1,ssmu[ssmu$ssmucode=="APDPE",])
  li<-Reduce(`+`, lis)
  apdpe<-sum(na.omit(li))
  lis<-extract(cpfmask1,ssmu[ssmu$ssmucode=="APBSE",])
  li<-Reduce(`+`, lis)
  apbse<-sum(na.omit(li))
  lis<-extract(cpfmask1,ssmu[ssmu$ssmucode=="APE",])
  li<-Reduce(`+`, lis)
  ape<-sum(na.omit(li))
  lis<-extract(cpfmask1,ssmu[ssmu$ssmucode=="APEI",])
  li<-Reduce(`+`, lis)
  apei<-sum(na.omit(li))
  
  #This is for all areas except the pelagic area. 
  #to extract the consumption in the pelagic area we first calculate the total from within all the other smmus
  shape<-apw+apdpw+apbsw+apdpe+apbse+ape+apei
  #then we take this off the total study area to calculate the amount from the pelagic area of the study area
  out<-cellStats(cpfmask1,sum)-shape 
  
  #now do the same for pelagic species
  lis<-extract(pelmask1,ssmu[ssmu$ssmucode=="APW",])
  li<-Reduce(`+`, lis)
  apwpel<-sum(na.omit(li))
  lis<-extract(pelmask1,ssmu[ssmu$ssmucode=="APDPW",])
  li<-Reduce(`+`, lis)
  apdpwpel<-sum(na.omit(li))
  lis<-extract(pelmask1,ssmu[ssmu$ssmucode=="APBSW",])
  li<-Reduce(`+`, lis)
  apbswpel<-sum(na.omit(li))
  lis<-extract(pelmask1,ssmu[ssmu$ssmucode=="APDPE",])
  li<-Reduce(`+`, lis)
  apdpepel<-sum(na.omit(li))
  lis<-extract(pelmask1,ssmu[ssmu$ssmucode=="APBSE",])
  li<-Reduce(`+`, lis)
  apbsepel<-sum(na.omit(li))
  lis<-extract(pelmask1,ssmu[ssmu$ssmucode=="APE",])
  li<-Reduce(`+`, lis)
  apepel<-sum(na.omit(li))
  lis<-extract(pelmask1,ssmu[ssmu$ssmucode=="APEI",])
  li<-Reduce(`+`, lis)
  apeipel<-sum(na.omit(li))
  
  
  shape1<-apwpel+apdpwpel+apbswpel+apdpepel+apbsepel+apepel+apeipel
  outpel<-cellStats(pelmask1,sum)-shape1 
  outpel
  
  
  kgperday_cpf<-c(apw,apdpw,apbsw,apdpe,apbse,ape,apei,out)
  kgperday_pelagics<-c(apwpel,apdpwpel,apbswpel,apdpepel,apbsepel,apepel,apeipel,outpel)
  
  #add the krill consumption within each ssmu to the data frame
  da<-as.data.frame(cbind(da,kgperday_cpf,kgperday_pelagics)) 
  da
  
  
  #the size of the area of each ssmu was already calculated and given in ssmu@data. However some of these ssmus were cut(the Antarctic Peninsula Pelagic Area, and Antarctic Peninsula East)
  # so calculate these again by taking the area of each ssmu from the size of the study area and divide it to make it into km2 
  areaout<-(area(Ps1)-area(ssmu[1,])-area(ssmu[2,])-area(ssmu[3,])-area(ssmu[4,])-area(ssmu[5,])-area(ssmu[6,])-area(ssmu[7,]))/1000000
  areaape<-(area(ssmu[ssmu$ssmucode=="APE",]))/1000000
  
  #so change these in the dataframe
  da$areakm2[da$ssmucode=="APE"]<-areaape
  da$areakm2[da$ssmucode=="APPA"]<-areaout
  
  
  
  # Now we calculate predator demand density in tonnes for all of summer.
  #we multiply by 182 to convert from day to summer, we /1000 to convert from kg to tonnes
  da$density_cpf<-da$kgperday_cpf/da$areakm2*182/1000
  da$density_pelagics<-da$kgperday_pelagics/da$areakm2*182/1000 #these are density of krill consumption (tonnes of krill per km2 all summer)
  
  # GW - the code provided by Vicky continues, but stop here, compute a couple things and pass back the output
  # this is all I need
  
  da$density_sum <- rowSums(da[,10:11])
  da$total_sum <- da$density_sum*da$areakm2
  
  da
  
}

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

CalcH <- function(indata){
  denomX <- sum(!is.na(indata[,6])&indata[,6]>0)
  denomY <- sum(!is.na(indata[,7])&indata[,7]>0)
  denomZ <- sum(!is.na(indata[,8])&indata[,8]>0)
  
  numerX <- sum(!is.na(indata[,6])&indata[,6]>0&indata[,9]>=0.1)
  numerY <- sum(!is.na(indata[,7])&indata[,7]>0&indata[,10]>=0.1)
  numerZ <- sum(!is.na(indata[,8])&indata[,8]>0&indata[,11]>=0.1)
  
  c(HX = numerX/denomX, HY = numerY/denomY, HZ = numerZ/denomZ)
}

PropSurveyOofM <- function(indata){
  tt<-round(indata[!is.na(indata[,5]),5],0)
  tt<-nchar(as.character(tt))
  ttt<-0
  for(i in 2:length(tt)){
    if(tt[i]!=tt[i-1]){ttt<-ttt+1}
  }
  ttt/(length(tt)-1)
}
