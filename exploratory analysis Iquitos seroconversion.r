
require(MASS)

require(foreign)
require(maptools)
require(sp)
require(sf)

# for neighbours
# require(spdep)

# for choropleth
# require(GISTools)

# for polygon distances
# require(spatstat)

# require(raster)

if(Sys.info()[['user']]=="eidenale"){
   setwd("C:\\Users\\eidenale\\work\\project\\methods\\spatial CRT\\binary outcome\\Iquitos")
   # setwd("C:\\Users\\NEAL ALEXANDER\\work\\project\\methods\\spatial CRT\\binary outcome\\Iquitos")
   # aquí se leen los archivos (computador de Neal)
}else{
   # aquí se puede agregar donde leer los archivos para un usuario que no sea Neal
}

# CurtainBlocks<-readShapePoly("Curtain_blocks")
CurtainBlocks<-st_read("Curtain_blocks.shp")
summary(CurtainBlocks)
   plot(CurtainBlocks)
   plot(CurtainBlocks[,"Cortina_10"])
#  View(CurtainBlocks)
# Each row indicates intervention or control, but not cluster number.  

CurtainBlocks$consecutivo<-1:I(dim(CurtainBlocks)[1])
CurtainBlocksCentroids  <-st_centroid(CurtainBlocks)
CurtainBlocksCoordinates<-st_coordinates(CurtainBlocksCentroids)
# View(CurtainBlocksCoordinates)
   plot(CurtainBlocks[,"Cortina_10"])
   # plot(CurtainBlocksCoordinates[,"consecutivo"], add=T)

# par(mfrow=c(2,1))
   eqscplot(x=CurtainBlocksCoordinates[,"X"], y=CurtainBlocksCoordinates[,"Y"], type='n')
       text(x=CurtainBlocksCoordinates[,"X"], y=CurtainBlocksCoordinates[,"Y"], labels=CurtainBlocks$consecutivo)
   plot(CurtainBlocks[,"Cortina_10"])
# par(mfrow=c(1,1))

# it doesn't seem to respect the mfrow.  
# But the two plots can be compared to link up the cluster numbers from the map in the MS (file "cortinas_2019") to the polygons in this shapefile

   plot(CurtainBlocks[CurtainBlocks$consecutivo%in%c(49, 14),"Cortina_10"])

CurtainBlocks$cluster<-
   ifelse(CurtainBlocks$consecutivo%in%c(48)            , 20,
   ifelse(CurtainBlocks$consecutivo%in%c(45, 46, 47)    , 19,
   ifelse(CurtainBlocks$consecutivo%in%c(42, 43, 44)    , 18,
   ifelse(CurtainBlocks$consecutivo%in%c(39, 40, 41)    , 17,
   ifelse(CurtainBlocks$consecutivo%in%c(37, 38)        , 16,
   ifelse(CurtainBlocks$consecutivo%in%c(35, 36)        , 15,
   ifelse(CurtainBlocks$consecutivo%in%c(33, 34)        , 14,
   ifelse(CurtainBlocks$consecutivo%in%c(29, 30, 31, 32), 13,
   ifelse(CurtainBlocks$consecutivo%in%c(27, 28)        , 12,
   ifelse(CurtainBlocks$consecutivo%in%c(24, 25, 26)    , 11,
   ifelse(CurtainBlocks$consecutivo%in%c(22, 23)        , 10,
   ifelse(CurtainBlocks$consecutivo%in%c(20, 21)        ,  9,
   ifelse(CurtainBlocks$consecutivo%in%c(17, 18, 19)    ,  8,
   ifelse(CurtainBlocks$consecutivo%in%c(15, 16)        ,  7,
   ifelse(CurtainBlocks$consecutivo%in%c(12, 13, 14)    ,  6,
   ifelse(CurtainBlocks$consecutivo%in%c( 9, 10, 11)    ,  5,
   ifelse(CurtainBlocks$consecutivo%in%c( 7,  8)        ,  4,
   ifelse(CurtainBlocks$consecutivo%in%c( 5,  6)        ,  3,
   ifelse(CurtainBlocks$consecutivo%in%c( 3,  4)        ,  2,
   ifelse(CurtainBlocks$consecutivo%in%c( 1,  2)        ,  1, NA
   ))))))))))))))))))))

CurtainBlocks[is.na(CurtainBlocks$cluster), c("cluster", "consecutivo", "Cortina_10")]

# CurtainLots<-readShapePoly("Curtain_lots")
CurtainLots<-st_read("Curtain_lots.shp")
summary(CurtainLots)
# View(CurtainLots)
   plot(CurtainLots)
   plot(CurtainLots[,"LANDUSE"])

CurtainPoints<-st_read("Curtains_points.shp")
summary(CurtainPoints)
  names(CurtainPoints)
# View(CurtainPoints)
   plot(CurtainPoints)
  class(CurtainPoints)
   plot(CurtainPoints[,"LOC_CODE"])

# CurtainPointsCoords<-st_coordinates(CurtainPoints)

# View(CurtainPoints[,c("POINT_X","POINT_Y")])
# from summary(CurtainPoints), looks the coords should be metres (UTM)
head(CurtainPoints[,c("LOC_CODE", "POINT_X","POINT_Y", "Cluster_Co")])


if(Sys.info()[['user']]=="eidenale"){
   setwd("C:\\Users\\eidenale\\work\\project\\dengue\\iquitos\\data\\Esther\\spatial")
   # aquí se leen los archivos (computador de Neal)
}else{
   # aquí se puede agregar donde leer los archivos para un usuario que no sea Neal
}

seroconversion<-read.table(file="IquitosAtRiskSeroconversion.txt", sep = "\t", header=T, 
           na.strings = c("NA", ".", ""), stringsAsFactors = FALSE)
names(seroconversion)
 dim(seroconversion)

which(substr(names(seroconversion), 1, 2)=="X."  | names(seroconversion)=="X")
seroconversion<-seroconversion[,!(substr(names(seroconversion), 1, 2)=="X."  | names(seroconversion)=="X")]
     names(seroconversion)
sort(names(seroconversion))

table(                      seroconversion$cluster_co)
table(seroconversion$treat, seroconversion$cluster_co)
# These are the numbers at risk.  They agree with Table 2 of the main trial results MS.  

table(seroconversion$seroconvert_any)
# numbers seroconverting (to any serotype) by cluster
# This can't be compared directly to the MS because the MS has rates.
table(seroconversion$treat[as.logical(seroconversion$seroconvert_any)], seroconversion$cluster_co[as.logical(seroconversion$seroconvert_any)])

seroconversion$collection_date1Date<-as.Date(seroconversion$collection_date1, "%d/%m/%Y")
seroconversion$collection_date2Date<-as.Date(seroconversion$collection_date2, "%d/%m/%Y")
head(seroconversion[,as.logical(match(substr(names(seroconversion), 1, 15), "collection_date", nomatch=0))])

DaysAtRisk<-as.numeric(seroconversion$collection_date2Date - seroconversion$collection_date1Date)
summary(DaysAtRisk)

seroconversion$DaysAtRisk<-DaysAtRisk
rm(DaysAtRisk)

RiskYearsByCluster<-tapply(seroconversion$DaysAtRisk     , seroconversion$cluster_co, sum)/(365.25)
NConvertByCluster <-tapply(seroconversion$seroconvert_any, seroconversion$cluster_co, sum)
rbind(NConvertByCluster, RiskYearsByCluster, round(100*NConvertByCluster/RiskYearsByCluster, 1))
# these also agree with Table 2 of the main trial MS

# Now try to match up ‘location_code’ in the seroconversion file with ‘LOC_CODE’ in the ‘Curtains_points’ shapefile.
# View(seroconversion[,"location_code"])


# table(match(seroconversion$location_code, CurtainPoints$LOC_CODE, nomatch=0))

dim(seroconversion)
seroconversion<-merge(
   x=seroconversion, 
   y=as.data.frame(CurtainPoints)[,c("LOC_CODE", "POINT_X","POINT_Y")],
   by.x="location_code", by.y="LOC_CODE", all.x=T, all.y=F)
dim(seroconversion)
# View(seroconversion[,c("location_code", "seroconvert_any", "POINT_X", "POINT_Y")])
head(seroconversion[,c("location_code", "cluster_co", "treat", "seroconvert_any", "POINT_X", "POINT_Y")])

seroconversion[is.na(seroconversion$POINT_X),
   c("consent_id", "participant_id", "part_code", "location_code", "seroconvert_any", "POINT_X", "POINT_Y")]

# sort(unique(CurtainPoints$LOC_CODE))

# In CurtainPoints there is SCA054 (but not SCA054B) and SCA749 (but not suffix A, B or C)

# Plan:
#   for each person in 'seroconversion', 
#   look up their location code in 'CurtainPoints',
#   check that this point is within one of the polygons for that cluster in 'CurtainBlocks' (using 'st_within').

# seroconversion[1,c("location_code", "cluster_co", "treat", "seroconvert_any", "POINT_X", "POINT_Y")]

for(i in 1:I(dim(seroconversion)[1])){
           location_codeScalar<-seroconversion[i, "location_code"]
   SeroconversionClusterScalar<-seroconversion[i, "cluster_co"]
   CurtainPointsSubset<-CurtainPoints[CurtainPoints$LOC_CODE==location_codeScalar,]
   if(I(dim(CurtainPointsSubset)[1])>1){
      print(unlist(list("row matching more than one CurtainPoint"=i)))
   }else{
      WithinSubset<-st_within(CurtainPointsSubset , CurtainBlocks, sparse=F)
      nWithinSubset<-sum(WithinSubset)
      if(nWithinSubset==0){
         print(unlist(list("row outside all blocks"=i)))
      }
   else{
      # row is unproblematic (so far)
      CurtainPointsClusterScalar<-as.data.frame(CurtainPoints)[i,"Cluster_Co"]
      if(SeroconversionClusterScalar!=SeroconversionClusterScalar){
         print(unlist(list(
            "row in seroconversion"=i, 
            "seroconversion cluster"=SeroconversionClusterScalar, 
             "CurtainPoints cluster"= CurtainPointsClusterScalar)))
      }
   }
   }
   #   dim(CurtainPointsSubset)
   # class(CurtainPointsSubset)
}


# CurtainPointsSubset
# # CurtainPoints[1207,]
# # CurtainPoints[1207,"Cluster_Co"]
# 
# st_within(CurtainPointsSubset , CurtainBlocks, sparse=F)
# # st_within(CurtainPoints[1207,], CurtainBlocks, sparse=F)
# 
# CurtainBlocks[st_within(CurtainPointsSubset , CurtainBlocks, sparse=F),]



# also use 'CurtainLots'
# View(CurtainLots)
names(CurtainLots)
class(CurtainLots)

st_within(CurtainLots[1,], CurtainBlocks, sparse=F)

for(i in 1:I(dim(seroconversion)[1])){
           location_codeScalar<-seroconversion[i, "location_code"]
   SeroconversionClusterScalar<-seroconversion[i, "cluster_co"]
   CurtainPointsSubset<-CurtainPoints[CurtainPoints$LOC_CODE==location_codeScalar,]
   if(I(dim(CurtainPointsSubset)[1])>1){
      print(unlist(list("row matching more than one CurtainPoint"=i)))
   }else{
      WithinSubset<-st_within(CurtainPointsSubset , CurtainLots, sparse=F)
      nWithinSubset<-sum(WithinSubset)
      if(nWithinSubset==0){
         print(unlist(list("row outside all lots"=i)))
      }
   else{
      # row is unproblematic (so far)
      CurtainPointsClusterScalar<-as.data.frame(CurtainPoints)[i,"Cluster_Co"]
      if(SeroconversionClusterScalar!=SeroconversionClusterScalar){
         print(unlist(list(
            "row in seroconversion"=i, 
            "seroconversion cluster"=SeroconversionClusterScalar, 
             "CurtainPoints cluster"= CurtainPointsClusterScalar)))
      }
   }
   }
   #   dim(CurtainPointsSubset)
   # class(CurtainPointsSubset)
}

