#' todo
#'
#' calculate predator demand density in tonnes for all of summer?
#'
#' @details
#' code and data files used here provided by Vicky Warwick-Evans to WG-EMM e-group or to GW via email
#' Todo: move paths to be function arguments
#'
#' @export
GetRAconsumption <- function(){

  # Unless otherwise noted, code and data files used here provided by Vicky Warwick-Evans to WG-EMM e-group or to GW via email

  # NOTE - there are several errors produced that I have not diagnosed
  # Also, I'm not clear if the results are consistent with Vicky's paper
  # I will use them for now, but probably need to revisit this

  # library(raster)
  # library(rgdal)
  # # GW - this also seems to be needed
  # library(rgeos)

  stopifnot(
    requireNamespace("raster"),
    requireNamespace("rgdal"),
    requireNamespace("rgeos")
  )

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
