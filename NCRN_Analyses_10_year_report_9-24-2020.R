#NCRN 11-year Forest Birds Report - Occupancy, Abundance, and BCI

#clear environment
rm(list=ls())

#load packages
library(ggplot2)
library(unmarked)

#setwd()
setwd("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN")

#load source functions
source("./Analyses/Source/FormatRawData_NCRN.R")
source("./Analyses/Source/CheckAOS_Codes.R")

#################################################################################################
#PART 1 - Format forest and grassland data

#read in forest data
forest.data<-read.csv("./Data/Bird_Data/qExport_R_UDel_Forest_2007_2017.csv",header=TRUE)
levels(forest.data$Survey_Type)

#check AOS_Codes
forest.data.1<-CheckAOS_Codes(DataIn=forest.data)

#remove spurious AOS Codes
forest.data.1<-subset(forest.data.1, c(AOU_Code != "WISP" & AOU_Code != "WPWI"))

#format data
formatted.data<-FormatRawData(DataIn = forest.data.1)

#read in grassland data
grassland.data<-read.csv("./Data/Bird_Data/qExport_NCRN_Grassland_Bird_Data.csv",header=TRUE)
levels(grassland.data$Survey_Type)

#check AOS_Codes
grassland.data.1<-CheckAOS_Codes(DataIn=grassland.data)

#Remove RODO and UNAH
grassland.data.1<-subset(grassland.data.1, c(AOU_Code != "RODO" & AOU_Code != "UNAH"))

#format data
FormatRawData(DataIn = grassland.data.1)

#################################################################################################
#################################################################################################
#################################################################################################
#Filtering data, and testing occupancy and abundance models

#clear environment (start fresh)
rm(list=ls())

#load packages
library(ggplot2)
library(unmarked)
library(reshape2)
library(AICcmodavg)
library(plyr)

#setwd()
setwd("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN")

#load source functions
source("./Analyses/Source/FormatRawData_NCRN.R")
source("./Analyses/Source/makePOccu_umf.R")
source("./Analyses/Source/runOccu.R")
source("./Analyses/Source/makePCount_umf.R")
source("./Analyses/Source/runPCount.R")
source("./Analyses/Source/makeGMMumf.R")
source("./Analyses/Source/runGmultMix.R")
source("./Analyses/Source/BCIappalachian.R")
source("./Analyses/Source/runBCI.R")
source("./Analyses/Source/summaryFunction.R")

#new source functions for running analyses across all years
source("./Analyses/Source/runOccuAllyears.R")
source("./Analyses/Source/runOccuAllyearsNetworkWide.R")
source("./Analyses/Source/runPCountAllyears.R")
source("./Analyses/Source/runPCountAllyearsNetworkWide.R")
source("./Analyses/Source/runPCountAllyearsNetworkWideSpecies.R")
source("./Analyses/Source/runPCountAllyearsSpecies.R")

#read in formatted data
data<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Data/Bird_Data/NCRN_Forest_Bird_Data_Formatted.csv",header=TRUE)

#species to remove
removeList<-c("UNBI","UNCH","UNHA","UNCR","UNDU","UNFL","UNOW","UNSP","UNSW","UNTH","UNWA","UNWO","UNWR")

data<-subset(data, ! AOU_Code %in% removeList)

#remove migrants 
migrantList<-c("Blackpoll Warbler", "Bay-breasted Warbler", "Canada Warbler", "Cape May Warbler", "Gray-cheeked Thrush", "Mourning Warbler", "Nashville Warbler", "Northern Waterthrush", "Swainson's Thrush", "Tennessee Warbler")

data<-subset(data, ! Common_Name %in% migrantList)


range(data$CountOfAOU_Code)
#################################################################################################
#filter out migratory species and those with fewere than 10 detections per year

#get all sampling locations within each year
allSamplingLocations<-unique(data[,c("Unit_Code","Plot_Name","Year")])

#get all species detection data
allSpecies<-data[,c("Common_Name","Unit_Code","Plot_Name","Year","CountOfAOU_Code")]

#get count of each species per point per year per park
allSpecies.agg<-aggregate(CountOfAOU_Code~Unit_Code+Plot_Name+Year+Common_Name, data=allSpecies, FUN="length")

#get sum of counts per species per year
allSpeciesPark.agg<-aggregate(CountOfAOU_Code~Unit_Code+Common_Name, data=allSpecies.agg, FUN="sum")

#filter removing any species
allSpeciesPark.sub<-subset(allSpeciesPark.agg, CountOfAOU_Code>50)
#create SpeciesPark column
allSpeciesPark.sub$SpeciesPark<-paste(allSpeciesPark.sub$Common_Name, allSpeciesPark.sub$Unit_Code, sep="+")

#make list to filter all data
speciesParkKeepList<-unique(allSpeciesPark.sub$SpeciesPark)

#make SpeciesPark column in data
data$SpeciesPark<-paste(data$Common_Name, data$Unit_Code, sep="+")

#filter data
data.filter<-subset(data, SpeciesPark %in% speciesParkKeepList)

#remove migrants ("Blackpoll Warbler")
data.filter.noMigrants<-subset(data.filter, Common_Name !="Blackpoll Warbler")

#save data without Migrants
write.csv(data.filter.noMigrants, file="/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Data/Bird_Data/NCRN_Forest_Bird_Data_Formatted_noMigrants.csv",row.names = FALSE)
############################################################################################################################

#function to get maxCountsPerVisit
data.sub<-subset(data, AOU_Code=="WOTH")
DataIn=data.sub

new.data.sub<-data.sub[,c("Unit_Code","Plot_Name","Year","Visit","CountOfAOU_Code","Common_Name")]

#get all points visited
all.points<-unique(data[,c("Unit_Code","Year","Plot_Name")])

#now aggregate data
new.data.agg<-aggregate(CountOfAOU_Code~Year+Unit_Code+Plot_Name, FUN="length", data=new.data.sub)

#convert >1 to 1 for occupancy
new.data.agg$CountOfAOU_Code<-ifelse(new.data.agg$CountOfAOU_Code>1,1,new.data.agg$CountOfAOU_Code)

#add zeros 
new.data.merge<-merge(all.points, new.data.agg, by=c("Year","Unit_Code","Plot_Name"), all.x=TRUE)

new.data.merge$CountOfAOU_Code<-ifelse(is.na(new.data.merge$CountOfAOU_Code),0, new.data.merge$CountOfAOU_Code)

#now get max counts per point across visits
max.counts<-aggregate(CountOfAOU_Code~Year+Unit_Code+Plot_Name, FUN="max", data=new.data.merge)


summary.dat<-summaryFunction(dataIn=max.counts, factor="Year",response="CountOfAOU_Code")

#shiny testing
all.coords<-unique(data[,c("Plot_Name","Lat_WGS84","Long_WGS84")])

points.agg<-aggregate(CountOfAOU_Code~Year+Unit_Code+Plot_Name, FUN="max", data=max.counts)

points.merge<-merge(points.agg,all.coords, by="Plot_Name",all.x=TRUE)

#################################################################################################
#test occu
umf<-makePOccuNetworkWide_umf(DataIn=data, SpeciesIn="AMGO")

mod<-occu(~1~Year, data=umf)

mod.predict<-unique(predict(mod, type="state",appendData=TRUE))

mf.gof.test(mod)

gof.test<-mb.gof.test(mod, nsim = 3, plot.hist = FALSE,report = NULL, parallel = FALSE)

gof.test$c.hat.est

source("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Analyses/Shiny/Bird_Data_Model_Explorer/NPS_Bird_Data_Model_Explorer/Source/mb.gof.test_source.R")

mb.gof.test.unmarkedFitOccu(mod,nsim=10,plot.hist=FALSE, report=NULL, parallel=TRUE)

mb.chisq(mod, print.table=TRUE)
#################################################################################################
#test Antietam AMCR

test.data<-subset(data, c(Unit_Code=="ANTI"))

umf.test<-makePOccu_umf(DataIn = test.data, SpeciesIn = "AMCR")

#################################################################################################
#Unit Testing - check months min and max dates of surveys
may.data<-subset(data, Month==5)
min(may.data$Day)

july.data<-subset(data, Month=7)
max(july.data$Day)

#################################################################################################
#remove flyovers
unique(data.filter.noMigrants$Flyover_Observed)
data.use<-subset(data.filter.noMigrants, Flyover_Observed==0)

#constrain to only 0-50m (double check this! Use 0-100 for grassland)
unique(data$Distance_id)
#data<-subset(data, Distance_id==1)

#################################################################################################
#test abundance estimation in ANTI for AMRO

ANTI.data<-subset(data, Unit_Code=="ANTI")

library(AICcmodavg)
library(plyr)
library(reshape2)
library(dplyr)
library(magrittr)
library(ggplot2)


#Fit N-mixture model

SpeciesName="American Goldfinch"

umf<-makePCount_umf(DataIn=data, SpeciesIn = "AMGO")

umf<-makePCountNetworkWide_umf(DataIn=data, SpeciesIn = "AMGO")


umf<-makePOccuNetworkWide_umf(DataIn=data, SpeciesIn="AMGO")

umf<-makePOccuNetworkWide_umf(DataIn=subset(data,Unit_Code=="ANTI"), SpeciesIn="AMGO")

#fit occu model
mod.test<-occu(~1~Year, data=umf)

gof.out<-mb.gof.test(mod.test, nsim = 10, plot.hist = FALSE,report = NULL, parallel = TRUE)

unique(predict(mod.test, type="state",append.data=FALSE))

umf@siteCovs$Plot_Name
range(umf@y)
sum(umf@y)

umf@obsCovs$Visit<-as.factor(umf@obsCovs$Visit)
umf@obsCovs$Observer<-as.factor(umf@obsCovs$Observer)

#fit model
mod<-pcount(~1~Year, data=umf, mixture="P")

#goodness-of-fit test
gof.out<-Nmix.gof.test(mod, nsim = 10, plot.hist = TRUE,
                       report = NULL, parallel = TRUE)
#gof.out
#get c-hat (measure of over or under dispersion of model)
c.hat<-gof.out$c.hat.est
#set c.hat.color
c.hat.color<-ifelse(c.hat<0.2, "red",ifelse(c.hat>4,"red","darkgreen"))

#get model-predicted esitmates
mod.predict<-unique(predict(mod.test, type="state",appendData=TRUE))

#get raw counts
max.count<-sapply(umf@y, max)

#get total detections
total.max.detections<-sum(max.count)
#make data.frame of max count per plot per year
data.raw<-data.frame(Year=umf@siteCovs$Year,Plot_Name=umf@siteCovs$Plot_Name, MaxCount=max.count)
#average max counts per year
data.agg<-aggregate(MaxCount~Year+Plot_Name, data=data.raw, FUN="max")
#summary function - mean and SE
data.summary<-summaryFunction(dataIn=data.agg, response="MaxCount",factor="Year")

#combine predicted and naive estimates

#get subset of model estimates
mod.predict.sub<-data.frame(unique(mod.predict[,c("Year","Predicted","SE")]))
names(mod.predict.sub)<-c("Year","mean","SE")

mod.predict.sub$Estimate<-"N-mixture Model"

#now get naive estimates
data.summary.sub<-data.summary[,c("Year","mean","SE")]
data.summary.sub$Estimate<-"Naive"

#merge naive and model-estimated abundance
data.combine<-rbind(data.summary.sub, mod.predict.sub)

#Make Estimate a factor
data.combine$Estimate<-as.factor(data.combine$Estimate)
levels(data.combine$Estimate)


#plot predicted and naive estimates
pred_naive_plot<-ggplot(data=data.combine, aes(x=Year, y=mean))+
  geom_point(aes(color=Estimate))+
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE, color=Estimate),width=0, alpha=0.7)+
  geom_line(aes(color=Estimate, group=Estimate))+
  scale_color_manual(values=c("royalblue1","gray50"))+
  theme(panel.background = element_rect(fill="white",color="black"),
        panel.border=element_rect(fill="transparent",color="black"),
        legend.position=c(0.15,0.83))+
  labs(y="Mean birds per point", color="Abundance\nestimate")+
  annotate("text", x = 9, y=max(data.combine$mean+data.combine$SE), label = paste("Total max detections: ",total.max.detections,sep=""))+
  annotate("text",x = 9, y=max(data.combine$mean+data.combine$SE)-1, label = paste ("GOF c-hat = ",round(c.hat,3),sep=""),color=c.hat.color)
pred_naive_plot

#################################################################################################
#################################################################################################
#run Occupancy All years pooled
occu.save<-runOccuAllyears(DataIn=data.use, modSel="AIC")

#save all output
write.csv(occu.save, file="./Analyses/Results/Occupancy/NCRN_occupancy_Forest_AIC_12-23-18.csv",row.names=FALSE)

##################################################################################################
#Occupancy running slow, use parallel processing

#load packages
library("doParallel")
library("parallel")
library("foreach")
library("dclone")

#set up cores
detectCores()
cl<-parallel::makeCluster(15, setup_strategy = "sequential")
#cl<-makeCluster(15)
registerDoParallel(cl)


#data.sub<-subset(data, Unit_Code !="PRWI")

#get list of units
unitList<-unique(sort(as.character(data.filter.noMigrants$Unit_Code)))

#parallel for loop
#foreach(i=1:length(unitList),.packages=c("reshape2","plyr","unmarked","AICcmodavg","doParallel","foreach","dclone")) %dopar% {

for(i in 1:length(unitList)){

  unitName<-unitList[i]
  
  unit.data<-subset(data.filter.noMigrants, Unit_Code==unitName)
 
  
  occu.save<-runOccuAllyears(DataIn=unit.data, modSel="AIC")
  
  #save all output
  write.csv(occu.save, file=paste("./Analyses/Results/Occupancy/", unitName,"_occupancy_Forest_POccu_AIC_0-100m_09-15-20.csv",sep=""),row.names=FALSE)
  
}

########################################################################################################################33
#run full model for each species on all NCRN data

#occupancy
ncrn.save<-runOccuAllyearsNetworkWide(DataIn=data.filter.noMigrants, modSel="AIC",networkName = "NCRN", nGOF=20)

write.csv(ncrn.save, file=paste("./Analyses/Results/Occupancy/", "NCRN","_occupancy_Forest_POccu_AIC_0-100m_NetworkWide_09-18-20.csv",sep=""),row.names=FALSE)

########################################################################################################################33
#for loop to read in each park and compile into one .csv

#get list of units
unitList<-unique(sort(as.character(data$Unit_Code)))

unitList<-c("NCRN",unitList)

all.park.results<-list()
for(i in 1:length(unitList)){
  
  unitName<-unitList[i]
  
  new.park.results<-read.csv(paste("./Analyses/Results/Occupancy/", unitName,"_occupancy_Forest_POccu_AIC_0-100m_06-24-20.csv",sep=""),header=TRUE)
  
  all.park.results<-rbind(all.park.results, new.park.results)
}


#save all output
write.csv(all.park.results, file="./Analyses/Results/Occupancy/NCRN_All_occupancy_Forest_POccu_AIC_0-100m_9-24-20.csv",row.names=FALSE)

##################################################################################################
#Abundance running slow, use parallel processing

#load packages
library("doParallel")
library("foreach")
library("dclone")

#set up cores
detectCores()
cl<-makeCluster(7)
registerDoParallel(cl)

#get list of units
unitList<-unique(sort(as.character(data.filter.noMigrants$Unit_Code)))
#unitList<-subset(unitList, unitList !="PRWI")

#parallel for loop
foreach(i=1:length(unitList),.packages=c("reshape2","plyr","unmarked","AICcmodavg","doParallel","foreach","dclone")) %dopar% {
  
  unitName<-unitList[i]
  
  unit.data<-subset(data.filter.noMigrants, Unit_Code==unitName)
  
  abun.save<-runPCountAllyears(DataIn=unit.data, modSel="AIC")
  
  #save all output
  write.csv(abun.save, file=paste("./Analyses/Results/Abundance/", unitName,"_abundance_Forest_PCount_AIC_0-100m_9-15-20.csv",sep=""),row.names=FALSE)
  
}

########################################################################################################################33
#abundance, use foreach, or will take forever

#load packages
library("doParallel")
library("parallel")
library("foreach")
library("dclone")


#stopCluster(cl)
detectCores()
cl<-parallel::makeCluster(5, setup_strategy = "sequential")
#cl<-makeCluster(15)
registerDoParallel(cl)

speciesList<-unique(data.filter.noMigrants$AOU_Code)

foreach(r=1:length(speciesList),.packages=c("unmarked","plyr","reshape2","AICcmodavg","doParallel","parallel","foreach","dclone")) %dopar% {
  
  speciesName<-as.character(speciesList[r])
  
  ncrn.save.abun<-runPCountAllyearsNetworkWideSpecies(DataIn=data.filter.noMigrants, modSel="AIC",networkName = "NCRN", speciesIn=speciesName, nGOF=5)
  
  write.csv(ncrn.save.abun, file=paste("./Analyses/Results/Abundance/NCRN_Species/", "NCRN_",speciesName,"_abundance_Forest_POccu_AIC_0-100m_NetworkWide_09-18-20.csv",sep=""),row.names=FALSE)
  
}

########################################################################################################################33
#for loop to compile individual species results for NCRN-wide

#first get list of files in the folder where they were saved
speciesFileList<-list.files("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Analyses/Results/Abundance/NCRN_Species",full.names = TRUE)

save.species.abun<-list()
for(i in 1:length(speciesFileList)){
  
  print(speciesFileList[i])
  
  new.file<-read.csv(speciesFileList[i],header=TRUE)
  
  save.species.abun<-rbind(save.species.abun, new.file)
  
}

#save as NCRN-wide estimates using previously-used naming convention
write.csv(save.species.abun, file="/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Analyses/Results/Abundance/NCRN_abundance_Forest_PCount_AIC_0-100m_06-28-20.csv",row.names = FALSE)

########################################################################################################################33
#get list of units
unitList<-unique(sort(as.character(data.filter.noMigrants$Unit_Code)))

#add "NCRN" to unitList

unitList<-c(unitList, "NCRN")

#for loop to read in each park and compile into one .csv
all.park.results<-list()
for(i in 1:length(unitList)){
  
  unitName<-unitList[i]
  print(unitName)
  
  new.park.results<-read.csv(paste("./Analyses/Results/Abundance/", unitName,"_abundance_Forest_PCount_AIC_0-100m_06-28-20.csv",sep=""),header=TRUE)
  
  all.park.results<-rbind(all.park.results, new.park.results)
  }


#save all output
write.csv(all.park.results, file="./Analyses/Results/Abundance/NCRN_All_abundance_Forest_PCount_AIC_0-100m_9-16-20.csv",row.names=FALSE)

####################################################################################
#runPcountAllyearsSpecies (species in parallel)

#load packages
library("doParallel")
library("foreach")
library("dclone")

#set up cores
detectCores()
cl<-makeCluster(7)
registerDoParallel(cl)


#pick Unit_Code to run
unitName<-"PRWI"

#get unit.data
unit.data<-subset(data, Unit_Code==unitName)

#get species list
speciesList<-sort(unique(as.character(unit.data$AOU_Code)))

foreach(r=1:length(speciesList),.packages=c("reshape2","plyr","unmarked","AICcmodavg","doParallel","foreach","dclone")) %dopar% {
  
  speciesName<-speciesList[r]
  
  runPCountAllyearsSpecies(DataIn=unit.data, SpeciesIn=speciesName, modSel="AIC")

}

#compile species results
species.save<-list()
for(i in 1:length(speciesList)){
  
  speciesName<-speciesList[i]
  
   new.sp.data<-read.csv(paste("./Analyses/Results/Abundance",unitName,paste(speciesName,"abundance_12-20-18.csv",sep="_"), sep="/"),header=TRUE)
  
  species.save<-rbind(species.save, new.sp.data)
  
}


write.csv(species.save, file=paste("./Analyses/Results/Abundance/", unitName,"_abundance_Forest_PCount_AIC_0-100m_12-12-18.csv",sep=""),row.names=FALSE)



#find out deal with NAs
species.NA<-subset(species.save, is.na(Predicted))
species.NA<-subset(species.NA, is.na(Metric))


NA.speciesList<-unique(sort(as.character(species.NA$AOU_Code)))


foreach(r=1:length(NA.speciesList),.packages=c("reshape2","plyr","unmarked","AICcmodavg","foreach","dclone","doParallel")) %dopar% {
  
  speciesName<-NA.speciesList[r]
  
  runPCountAllyearsSpecies(DataIn=unit.data, SpeciesIn=speciesName, modSel="AIC")
  
}


#################################################################################################
#Run Abundance models PCount (constrained to 0-50m, no flyovers, and 2 visits for Forest)
abun.save<-runPCount(DataIn=data, modSel="maxDet")

#save all output
write.csv(abun.save, file="./Analyses/Results/Abundance/NCRN_abundance_Forest_PCount_maxDet_0-100m.csv",row.names=FALSE)

#################################################################################################
#Run Abundance models PCountAllyears (constrained to 0-50m, no flyovers, and 2 visits for Forest)
abun.save<-runPCount(DataIn=data, modSel="AIC")

#save all output
write.csv(abun.save, file="./Analyses/Results/Abundance/NCRN_abundance_Forest_PCount_AIC_0-100m.csv",row.names=FALSE)


#################################################################################################
#Run Abundance models Gmultmix (constrained to 0-100m, no flyovers, and 3 visits for Grassland)
habitat.save<-runGmultMix(DataIn=data)

#save all output
write.csv(habitat.save, file=paste(getwd(), "Analyses","Abundance","Results","NCRN_abundance_results_Grassland_Gmultmix.csv",sep="/"),row.names=FALSE)

#################################################################################################
#Run BCI estimates (all data including flyovers)
bci.save<-runBCI(DataIn=data)

#save all output
write.csv(bci.save, file="./Analyses/Results/BCI/NCRN_BCI_results.csv",row.names=FALSE)

##########################################################################################
##########################################################################################
#PART 4 - Generate Tables and Figures for Report

#Glossary 1

#List of region and park codes within the National Capital Region Network

#create UnitNames.df
UnitNames.df<-unique(data[,c("Unit_Code","Unit_Name")])

row.names(UnitNames.df)<-NULL
#order by Unit_Name
UnitNames.df$Unit_Name<-trimws(as.character(UnitNames.df$Unit_Name))
UnitNames.df.order<-UnitNames.df[order(UnitNames.df$Unit_Name, decreasing=FALSE),]

#add NCRN
UnitNames.df.out<-rbind(data.frame(Unit_Code="NCRN",Unit_Name="National Capital Region Network"), UnitNames.df.order)
colnames(UnitNames.df.out)<-c("4-letter Code", "Name")
#save table
write.csv(UnitNames.df.out, file="./Analyses/Tables/Glossary_1_Park_names.csv", row.names=FALSE)

#################################################################################################
#Glossary 2

# List of all species detected within National Capital Region Network park units. Four-letter
# American Ornithological Society (AOS) codes are shown. Species of continental conservation
# concern as designated by Partners in Flight are shown in bold and gray-highlighted rows.

#read in PIF data
pif.data<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Data/PIF_data/PIF_2017_Global.csv")

pif.data.2<-pif.data[,c("Common.Name","Scientific.Name","Continental.Concern","IUCN.Red.List.2016")]
colnames(pif.data.2)[1]<-"Common_Name"

#species to remove
removeList<-c("UNBI","UNCH","UNHA","UNCR","UNDU","UNFL","UNOW","UNSP","UNSW","UNTH","UNWA","UNWO","UNWR")

#simplify data (get list of unique species)
unique.birds<-unique(data[,c("AOU_Code","Common_Name")])

unique.birds.2<-subset(unique.birds, !AOU_Code %in% removeList)

#now merge with PIF status
birds.merge.2<-merge(unique.birds.2, pif.data.2, by="Common_Name",all.x=TRUE)
birds.merge.2$Scientific.Name<-as.character(birds.merge.2$Scientific.Name)

#save as .csv
write.csv(birds.merge.2, file="./Analyses/Tables/Glosssary_2_Species_List.csv", row.names=FALSE)
#################################################################################################
#Figure 1: Study Area map

library(ggmap)
library(ggplot)
library(rgdal)
library(sp)
library(ggrepel)
library(sf)
library(grid)
#read in NPS shapefiles

park.shp<-readOGR("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/GIS/allparkboundariesshp.shp", stringsAsFactors=FALSE, verbose=FALSE)

#reproject shapefile to latlong
park.shp.proj<-spTransform(park.shp, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

plot(park.shp.proj, col="purple")

#fortify to create a d.f.
park.shp.df.1<-fortify(park.shp.proj)

park.shp.df.1.sub<-subset(park.shp.df.1, c(group != 1.1 & group != 14.1))

unique(park.shp.df.1$id)


#separate NACE into PISC/FOWA, NACE, and GREE
park.names<-data.frame(Unit_Code=c("ANTI","CATO","CHOH","GWMP","HAFE","MANA","MONO","GREE","NACE","PISC/FOWA","NACA","PRWI","ROCR","WOTR"), id=c("0","2","3","4","5","6","7","8","9","10","11","12","13","15"))

park.shp.merge<-merge(park.shp.df.1.sub, park.names, by="id",all.x=TRUE)


#get NACE
nace.shp<-subset(park.shp.merge, Unit_Code=="NACE")

ggplot(data=nace.shp,aes(x=long, y=lat,group=group))+
  geom_polygon(data=nace.shp, aes(fill=piece))+
  theme(legend.position="none")

#piscataway = piece 1
#greenbelt = piece 2
#FOWA = piece 7
ggplot(data=subset(nace.shp, piece=="2"),aes(x=long, y=lat,group=group))+
  geom_polygon(data=subset(nace.shp, piece=="2"), aes(fill=piece))+
  theme(legend.position="none")

#now subset shapefiles
pisca.shp<-subset(nace.shp, piece=="1")
fowa.shp<-subset(nace.shp, piece=="7")
gree.shp<-subset(nace.shp, c(piece=="2" & lat>38.97 & long< -76.88))
gree.nace.shp<-subset(nace.shp, c(piece=="2" & lat<38.97 & long> -76.88))


nace.shp.sub<-subset(nace.shp, c(piece!="1" & piece !="7" & piece != "2"))
nace.only.shp<-rbind(nace.shp.sub, gree.nace.shp)

ggplot(data=gree.nace.shp,aes(x=long, y=lat,group=group))+
  geom_polygon(data=gree.shp, aes(fill=piece))+
  theme(legend.position="none")

#combine pisc.fowa
pisc.fowa.shp<-rbind(pisca.shp, fowa.shp)

park.id<-unique(park.shp.merge[,c("Unit_Code","id")])

#fix ids
pisc.fowa.shp$id<-"10"
pisc.fowa.shp$Unit_Code<-"PISC/FOWA"
# pisc.fowa.shp$group<-ifelse(pisc.fowa.shp$group=="9.7","9.2",pisc.fowa.shp$group)

gree.shp$id<-"8"
gree.shp$Unit_Code<-"GREE"
gree.shp$group<-"8.1"

nace.only.shp$id<-"9"
nace.only.shp$Unit_Code<-"NACE"

#now remove NACE from all shapefiles
park.shp.merge.sub<-subset(park.shp.merge, c(Unit_Code !="GREE" & Unit_Code !="NACE" & Unit_Code !="PISC/FOWA"))


#combine with rbind
park.shp.merge.2<-rbind(rbind(rbind(park.shp.merge.sub, pisc.fowa.shp), gree.shp), nace.only.shp)

#reorder
park.shp.merge.order<-park.shp.merge.2[order(park.shp.merge.2$group, park.shp.merge.2$order),]
park.shp.merge.order$Unit_Code<-factor(park.shp.merge.order$Unit_Code, levels=sort(unique(park.shp.merge.order$Unit_Code)))

#get mean coords
meanLon<-mean(park.shp.merge.order$long)
meanLat<-mean(park.shp.merge.order$lat)

#get basemap

#register Google Maps API key
register_google(key="AIzaSyDerv2JQubJ_Dg8_RKBZCeyo5zqsM5USw0")

basemap<-get_map(location=c(lon=meanLon, lat=meanLat), maptype="terrain",zoom=9, color="bw")
ggmap(basemap)


#now aggregate to get mean lat/long for each Unit
park.long<-aggregate(long~Unit_Code+id, data=park.shp.merge.order, FUN="mean")
park.lat<-aggregate(lat~Unit_Code+id, data=park.shp.merge.order, FUN="mean")
names(park.lat)
#combine
park.coords.merge<-merge(park.long, park.lat, by=c("Unit_Code","id"),all=TRUE)

#move CHOH label a little
park.coords.merge$long<-ifelse(park.coords.merge$Unit_Code=="CHOH",-77.5, park.coords.merge$long)
park.coords.merge$lat<-ifelse(park.coords.merge$Unit_Code=="CHOH",39.2, park.coords.merge$lat)

#get mapExtent for scale bar
mapExtent<-attr(basemap, "bb")

#add polygons
NCRN_studyAreaMap<-ggmap(basemap, darken=c(0.4,"gray30"))+
  geom_map(data=subset(park.shp.merge.order, c(Unit_Code !="" & Unit_Code !="NACA")),map=park.shp.merge.order, aes(x=long, y=lat, map_id=id,fill=Unit_Code),alpha=0.8)+
  labs(x="Longitude",y="Latitude")+
  theme(
    axis.text=element_text(size=18),
    axis.title=element_text(size=22),
    panel.border=element_rect(fill="transparent", color="black"),
    legend.position="none"
    # legend.position = c(0.15,0.3),
    # legend.background = element_rect(fill=alpha("white",0.6))
  )+
  guides(fill=guide_legend(title="Park Unit"))+
  geom_label_repel(data=subset(park.coords.merge, c(Unit_Code !="" & Unit_Code !="NACA")) ,aes(x=long,y=lat, label=Unit_Code, color=Unit_Code), direction="both",fill=alpha(c("gray10"),0.8),show.legend=FALSE, label.padding=0.3, box.padding = 0.4, size=4, label.size=NA, segment.alpha=0, nudge_x=0.05, nudge_y=-0.006)+
  ggsn::scalebar(x.min=mapExtent$ll.lon, x.max=mapExtent$ur.lon, y.min=mapExtent$ll.lat, y.max=mapExtent$ur.lat,transform=TRUE, dist_unit="km",dist=20, st.dist=0.018, height = 0.01,anchor=c(x=-76.8,y=38.5))
NCRN_studyAreaMap

##################################################
#now produce map with inset

#get vector map of USA with State Borders
states <- ggplot2::map_data("state")

USmap <- ggplot(data = states,
                mapping = aes(x = long, y = lat,
                              group = group))

USmap.out<-USmap + geom_polygon(fill = "gray60", color = "white")+
  coord_map(xlim = c(-84, -67), ylim=c(25,47))+
  geom_point(aes(x=-77.4, y=39.0), color="red",shape=0, size=8, stroke=1.5)+
  theme(
    plot.background = element_blank(),
    panel.background=element_rect(fill=alpha("white",0.7),color="transparent"),
    panel.grid = element_blank(),
    panel.border=element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    axis.title = element_blank()
  )

USmap.out

library(cowplot)
#now inset eastern US map
gg_inset_map1 = ggdraw() +
  draw_plot(NCRN_studyAreaMap) +
  draw_plot(USmap.out, x = 0.14, y = -0.05, width = 0.25, height = 0.7)

#gg_inset_map1

#save plot (Figure 1)
# ggsave(gg_inset_map1, file="/Users/zach/Dropbox (ZachTeam)/Projects/eDNA_Brazil/Pilot_study/Results/Figures/Figure_1_eDNA_Study_Area_color.png",width=8, height=8, dpi=600)


###################################################
#add north arrow and save
png(filename="/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Analyses/Figures/NCRN_Study_Area_Map.png",width=11, height=11, units="in",res=600)
ggsn::north2(gg_inset_map1, x=0.8, y=0.8,symbol=12, scale=0.07)

dev.off()

#################################################################################################
#Table 1. Table to show how many visits to which points per year

#make GRTS.Year.Survey_Type column
data$Unit.Year<-paste(data$Unit_Code, data$Year,sep=".")

#levels(data$Survey_Type)
#subset data to inclued unique GRTS.Year rows
data.sub<-unique(data[,c("Unit_Code","Unit_Name","Year","Unit.Year","Plot_Name")])

#create UnitNames.df
UnitNames.df<-unique(data.sub[,c("Unit_Code","Unit_Name")])
row.names(UnitNames.df)<-NULL
#add NCRN
UnitNames.df<-rbind(UnitNames.df, data.frame(Unit_Code="NCRN",Unit_Name="National Capital Region Network"))

#make frequency table of num. visits per GRTS.Year.Survey_Type
table1<-aggregate(data.sub$Plot_Name ~ data.sub$Unit.Year, FUN="length")
colnames(table1)<-c("Unit.Year","Plot_Count")

#separate GRTS, Year, and Survey_Type
table1.covs<-read.table(text=as.character(table1$Unit.Year),sep=".")
colnames(table1.covs)<-c("Unit_Code","Year")
table1<-cbind(table1, table1.covs)

#melt and cast (pivot data)
table1.melt<-melt(table1, id.vars=c("Unit_Code","Year"), measure.vars="Plot_Count")
table1.cast<-cast(table1.melt, Unit_Code ~ Year, fun.aggregate ="max",fill=NA)
ncrn.sum<-data.frame(Unit_Code= "NCRN", t(colSums(table1.cast[2:12])))
colnames(ncrn.sum)[2:12]<-as.character(seq(from=2007, to=2017, by=1))

table1.comb<-rbind(ncrn.sum, table1.cast)

#add unit names
table1.merge<-merge(table1.comb, UnitNames.df, by="Unit_Code", all.x=TRUE,sort=FALSE)

#organize
table1.out<-cbind(table1.merge$Unit_Name, table1.merge$Unit_Code, table1.comb[,2:12])
colnames(table1.out)[1:2]<-c("Name","Code")

#save as .csv
write.csv(table1.out, file="./Analyses/Tables/Table_1_Plot_Count.csv", row.names = FALSE)

#################################################################################################
#Figure 2. Species Richness per park

birds.merge.2$PIF<-ifelse(birds.merge.2$Continental.Concern !="",1,0)
pif.merge<-merge(data, birds.merge.2, by="AOU_Code",all.x=TRUE)
pif.merge.df<-unique(pif.merge[,c("AOU_Code","Unit_Code","Continental.Concern","IUCN.Red.List.2016","PIF")])

#NCRN-wide species richness
ncrn.species.count<-data.frame(Unit_Code="NCRN", Count=length(unique(pif.merge.df$AOU_Code))) #157
ncrn.species.count$Type<-"Total"
ncrn.cc.count<-data.frame(Unit_Code="NCRN", Count=sum(birds.merge.2$PIF)) #27
ncrn.cc.count$Type<-"Continental Concern"
ncrn.species.df<-rbind(ncrn.species.count, ncrn.cc.count)


#get species count by park
park.species.count<-aggregate(data$AOU_Code~data$Unit_Code, FUN=function(x){length(unique(x))})
park.species.count$Type<-"Total"
colnames(park.species.count)<-c("Unit_Code","Count","Type")
park.cc.count<-aggregate(pif.merge.df$PIF~pif.merge.df$Unit_Code, FUN=sum)
park.cc.count$Type<-"Continental Concern"
colnames(park.cc.count)<-c("Unit_Code","Count","Type")

#combine
park.species.df<-rbind(park.species.count, park.cc.count)

#combine park and network-level counts
species.count.df<-rbind(ncrn.species.df, park.species.df)

#add park names
colnames(UnitNames.df.out)<-c("Unit_Code", "Unit_Name")
species.count.merge<-merge(species.count.df, UnitNames.df.out, by="Unit_Code")

#reorder by count

#get order of parks by total species (not cc)
species.count.sub<-subset(species.count.merge, Type=="Total")
species.count.sub.order<-species.count.sub[order(species.count.sub$Count, decreasing=FALSE),]

species.count.merge$Unit_Name<-factor(species.count.merge$Unit_Name, levels=unique(species.count.sub.order$Unit_Name))
species.count.merge$Type<-as.factor(species.count.merge$Type)
#create figure show counts on columns as labels

#create colors
mycolors=c("darkgreen","orange")

species.count.merge$Type<-factor(species.count.merge$Type, levels=c("Total","Continental Concern"))

spp.plot<-ggplot(data=species.count.merge)+
  # geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE),color="darkolivegreen", width=0)+
  geom_bar(data=species.count.merge, stat="identity",aes(x=Unit_Name, y=Count, fill=Type, group=Type))+
  theme(panel.background = element_rect(fill="transparent"),
        panel.border= element_rect(color="black", fill=NA),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))+
  theme(legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),legend.position=c(0.8,0.2))+
  scale_fill_manual(values=mycolors)+
  labs(x="Network/Park",y="Species Richness")+
  geom_text(aes(x=Unit_Name, y=Count, label=Count, group=Type),hjust = 1.3, 
            color="white", size = 3.5, fontface="bold",inherit.aes = TRUE)+
  coord_flip()
spp.plot

#save plot
ggsave(spp.plot, file="./Analyses/Figures/Figure_1_Species_Richness.png",width=8, height=7, dpi=600)

#################################################################################################
#Table 2. Count of unique species per park per year

#make GRTS.Year.Survey_Type column
data$Unit.Year<-paste(data$Unit_Code, data$Year,sep=".")

#levels(data$Survey_Type)
#subset data to inclued unique GRTS.Year rows
data.sub<-unique(data[,c("Unit_Code","Year","Unit.Year","AOU_Code")])

#make frequency table of num. visits per GRTS.Year.Survey_Type
table2<-aggregate(data.sub$AOU_Code ~ data.sub$Unit.Year, FUN="length")
colnames(table2)<-c("Unit.Year","CountOfAOU_Code")

#separate GRTS, Year, and Survey_Type
table2.covs<-read.table(text=as.character(table2$Unit.Year),sep=".")
colnames(table2.covs)<-c("Unit_Code","Year")
table2<-cbind(table2, table2.covs)

#melt and cast (pivot data)
table2.melt<-melt(table2, id.vars=c("Unit_Code","Year"), measure.vars="CountOfAOU_Code")
table2.cast<-cast(table2.melt, Unit_Code ~ Year, fun.aggregate ="max",fill=NA)


#get total unique species for NCRN per year
data.sub.2<-unique(data[,c("Year","AOU_Code")])

#make frequency table of num. visits per GRTS.Year.Survey_Type
table2.b<-aggregate(data.sub.2$AOU_Code ~ data.sub.2$Year, FUN="length")
colnames(table2.b)<-c("Year","CountOfAOU_Code")


ncrn.tot<-data.frame(Unit_Code= "NCRN", t(table2.b$CountOfAOU_Code))
colnames(ncrn.tot)[2:12]<-as.character(seq(from=2007, to=2017, by=1))

table2.comb<-rbind(ncrn.tot, table2.cast)

#add unit names
table2.merge<-merge(table2.comb, UnitNames.df.out, by="Unit_Code", all.x=TRUE,sort=FALSE)

#organize
table2.out<-cbind(table2.merge$Unit_Name, table2.merge$Unit_Code, table2.comb[,2:12])
colnames(table2.out)[1:2]<-c("Name","Code")

row.names(table2.out)<-NULL

#save .csv
write.csv(table2.out, file="./Analyses/Tables/Table_2_Count_Species_by_Park_Year.csv", row.names=FALSE)

#################################################################################################
#Table 3.List of PIF species

#read in PIF data
pif.data<-read.csv("./Data/PIF_data/PIF_2017_Global.csv")

#truncate pif.data
pif.data.2<-pif.data[,c("Common.Name","Scientific.Name","Continental.Concern","IUCN.Red.List.2016")]
colnames(pif.data.2)[1]<-"Common_Name"

#merge with data
data.merge<-merge(data, pif.data.2, by="Common_Name",all.x=TRUE)

#create Year.GRTS column
data.merge$Unit.Year<-paste(data.merge$Unit_Code, data.merge$Year, sep=".")

#subset to include only species of conservation concern
data.sub<-subset(data.merge, Continental.Concern !="")

#get unique list of PIF species
data.pif.out<-unique(data.sub[,c("Common_Name","AOU_Code","Scientific.Name","Continental.Concern","IUCN.Red.List.2016")])
length(row.names(data.pif.out)) #27 PIF species

#save for table in report
write.csv(data.pif.out, file="./Analyses/Tables/Table_3_PIF_species.csv",row.names=FALSE)

#################################################################################################
#Table 4. Proportion of monitoring locations with PIF species

#create UnitNames.df
UnitNames.df<-unique(data[,c("Unit_Code","Unit_Name")])

row.names(UnitNames.df)<-NULL
#order by Unit_Name
UnitNames.df$Unit_Name<-trimws(as.character(UnitNames.df$Unit_Name))
UnitNames.df.order<-UnitNames.df[order(UnitNames.df$Unit_Name, decreasing=FALSE),]

#add NCRN
UnitNames.df.out<-rbind(data.frame(Unit_Code="NCRN",Unit_Name="National Capital Region Network"), UnitNames.df.order)

#read in PIF data
pif.data<-read.csv("./Data/PIF_data/PIF_2017_Global.csv")

#truncate pif.data
pif.data.2<-pif.data[,c("Common.Name","Scientific.Name","Continental.Concern","IUCN.Red.List.2016")]
colnames(pif.data.2)[1]<-"Common_Name"

#merge with data
data.merge<-merge(data, pif.data.2, by="Common_Name",all.x=TRUE)

#create Year.GRTS column
data.merge$Unit.Year<-paste(data.merge$Unit_Code, data.merge$Year, sep=".")

#subset to include only species of conservation concern
data.sub<-subset(data.merge, Continental.Concern !="")

data.sub<-unique(data.sub[,c("Unit_Name","Unit_Code","Plot_Name","Unit.Year","Year")])
data.sub$Year<-as.factor(as.character(data.sub$Year))

#get Network-wide proportions
#create table (CC spp)
ncrn.cc.count<-aggregate(data.sub$Plot_Name ~ data.sub$Year, FUN=length)
colnames(ncrn.cc.count)<-c("Unit.Year","CountOfPoints")
ncrn.cc.count$Unit.Year<-paste("NCRN",ncrn.cc.count$Unit.Year, sep=".")

#create table (CC spp)
park.cc.count<-aggregate(data.sub$Plot_Name ~ data.sub$Unit.Year, FUN=length)
colnames(park.cc.count)<-c("Unit.Year","CountOfPoints")

#combine
spp.all.count<-rbind(ncrn.cc.count, park.cc.count)

#combine ncrn and parks
count.pts<-rbind(ncrn.cc.count, spp.cc.count)

#create table of total points in each park
points.all.count<-aggregate(data.merge$Plot_Name ~ data.merge$Unit.Year, FUN=function(x) length(unique(x)))
colnames(points.all.count)<-c("Unit.Year","CountOfPoints")

#add NCRN 
ncrn.df<-data.frame(Unit.Year=paste("NCRN",row.names(t(ncrn.sum)),sep="."), CountOfPoints=as.numeric(t(ncrn.sum)))[-1,]
points.all.count<-rbind(ncrn.df,points.all.count)

#species prop table
spp.prop.df<-merge(points.all.count, spp.all.count, by=c("Unit.Year"),all.x=TRUE)
spp.prop.df$Proportion<-round(spp.prop.df$CountOfPoints.y/spp.prop.df$CountOfPoints.x,2)

#get covs
table.covs<-read.table(text=as.character(spp.prop.df$Unit.Year), sep=".")
colnames(table.covs)<-c("Unit_Code","Year")
table.covs<-merge(table.covs, UnitNames.df.out, by="Unit_Code",all.x=TRUE,sort=FALSE)
colnames(table.covs)<-c("Unit_Code","Year", "Unit_Name")
table.covs<-table.covs[,c("Unit_Name","Unit_Code","Year")]

#create data frame
spp.prop.df.2<-data.frame(table.covs, Count_of_Plots=spp.prop.df$CountOfPoints.x , Proportion=spp.prop.df$Proportion)

#now melt and cast to get values of years across top of table
spp.prop.melt<-melt(spp.prop.df.2, id.vars=c("Unit_Name","Unit_Code","Year"), measure.vars=c("Proportion"))
spp.prop.cast<-cast(spp.prop.melt, Unit_Name+Unit_Code ~ Year, fun.aggregate = max, fill=0)

spp.prop.cast$Unit_Name<-as.character(spp.prop.cast$Unit_Name)

#fill NAs with 0s
spp.prop.cast[is.na(spp.prop.cast)] <- 0

#set column names
colnames(spp.prop.cast)[1:2]<-c("Name","4-letter Code")

#save as .csv
write.csv(spp.prop.cast, file="./Analyses/Tables/Table_4_Proportion_of_Points_CC.csv",row.names=FALSE)


#################################################################################################
#Figure 2. Detection Probability (NCRN) mean and SE

#read in formatted data
data<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Data/Bird_Data/NCRN_Forest_Bird_Data_Formatted.csv", header=TRUE)

common.names<-unique(data[,c("AOU_Code","Common_Name")])

#read in occupancy results
occu.results<-unique(read.csv(file="./Analyses/Results/Occupancy/NCRN_occupancy_Forest_maxDet_9-24-18.csv"))

#species to remove
removeList<-c("UNBI","UNCH","UNHA","UNCR","UNDU","UNFL","UNOW","UNSP","UNSW","UNTH","UNWA","UNWO","UNWR")

#remove unwanted species ()
occu.results<-subset(occu.results, !AOU_Code %in% removeList)

#retain needed columns
occu.det<-unique(occu.results[,c("Unit_Code","Year","AOU_Code","Modnames","OverallDet")])

#create Park.Species column
occu.det$Park.Species<-paste(occu.det$Unit_Code, occu.det$AOU_Code, sep=".")

#average across years
occu.summary<-summaryFunction(dataI=occu.det, factor="AOU_Code", response="OverallDet")

#add common namees
occu.merge<-merge(occu.summary, common.names, by=c("AOU_Code"))

#split occu.merge into 2 data sets

length(row.names(occu.merge))

occu.merge.1<-head(occu.merge,length(row.names(occu.merge))/2)

#get list of species in occu.merge.1
occu.merge.1.list<-as.character(occu.merge.1$Common_Name)

occu.merge.2<-subset(occu.merge, ! Common_Name  %in% occu.merge.1.list)



#plot
det.prob.plot.1<-ggplot(data=occu.merge.1, aes(x=Common_Name, y=mean))+
  coord_flip()+
  geom_bar(stat="identity", fill="darkgreen",color=I("gray40"),size=0.2)+
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), postion="dodge", color="darkgreen",width=0)+
  theme(panel.background=element_rect(fill="transparent"), panel.border=element_rect(fill=NA,color="black"),
        axis.text.x=element_text(angle=60, hjust=1,size=14),
        axis.text.y=element_text(size=14),
        axis.title = element_text(size=16))
  det.prob.plot.1

#plot
det.prob.plot.2<-ggplot(data=occu.merge.2, aes(x=Common_Name, y=mean))+
  coord_flip()+
  geom_bar(stat="identity", fill="darkgreen",color="white",size=0.055)+
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), postion="dodge", color="darkgreen",width=0)+
  theme(panel.background=element_rect(fill="transparent"), panel.border=element_rect(fill=NA,color="black"),
        axis.text.x=element_text(angle=0, hjust=1,size=14),
        axis.text.y=element_text(size=8),
        axis.title = element_text(size=16))
det.prob.plot.2


#########################################################################################
#Figure 5 - BCI figure

#clear environment (start fresh)
rm(list=ls())

#load packages
library(ggplot2)
library(unmarked)
library(reshape2)
library(AICcmodavg)
library(plyr)
library(magrittr)
library(dplyr)

#setwd()
setwd("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN")

#load source functions
source("./Analyses/Source/summaryFunction.R")

#read in formatted data
data<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Data/Bird_Data/NCRN_Forest_Bird_Data_Formatted_noMigrants.csv",header=TRUE)

#create UnitNames.df
UnitNames.df<-unique(data[,c("Unit_Code","Unit_Name")])

row.names(UnitNames.df)<-NULL
#order by Unit_Name
UnitNames.df$Unit_Name<-trimws(as.character(UnitNames.df$Unit_Name))
UnitNames.df.order<-UnitNames.df[order(UnitNames.df$Unit_Name, decreasing=FALSE),]

#add NCRN
UnitNames.df.out<-rbind(data.frame(Unit_Code="NCRN",Unit_Name="National Capital Region Network"), UnitNames.df.order)

#Get NCRN-wide annual BCI and linear trend

bci.data<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Analyses/Results/BCI/NCRN_BCI_results.csv", header=TRUE)

bci.summary.network<-summaryFunction(dataIn=bci.data, response="BCI", factor="yearName")
bci.summary.network$unitName<-"NCRN"
#convert year to integer
bci.summary.network$yearName<-as.integer(as.character(bci.summary.network$yearName))
bci.order<-bci.summary.network[order(bci.summary.network$mean),]

bci.lm<-lm(mean~yearName,data=bci.summary.network)
linear.trend<-round(bci.lm$coefficients[2],3)
trend.txt<-paste("Linear trend in BCI:",linear.trend,sep=" ")

bci.summary.network$Trend<-linear.trend

#Get Park-level annual mean BCI and trend

#combine annual mean BCI for network and parks
summary.network<-summaryFunction(dataIn=bci.data, response="BCI", factor="yearName")
summary.network<-cbind("NCRN",summary.network)
names(summary.network)[names(summary.network)=="yearName"]<-"Year"
colnames(summary.network)[1]<-c("Park")

#create Park.Year column
bci.data$Park.Year<-paste(bci.data$unitName, bci.data$yearName,sep=".")

#get mean annual BCI by park 
summary.park<-summaryFunction(dataIn=bci.data, response="BCI", factor="Park.Year")
covs<-read.table(text=as.character(summary.park$Park.Year), sep=".")
colnames(covs)<-c("Park","Year")
summary.park<-cbind(covs, summary.park)
summary.park$Park.Year<-NULL
head(summary.park)

#combine network and park mean annual BCI
summary.all<-rbind(summary.network, summary.park)


#get trend for each park and add to bci.data

#get list of parks
parkList<-c("NCRN",unique(sort(as.character(bci.data$unitName))))

trend.out<-list()
for(i in 1:length(parkList)){
  
  new.data<-summary.all
  new.sub<-subset(new.data, Park==parkList[i])
  
  row.length<-length(row.names(new.sub))
  
  #fit linear model and get slope (Beta coefficient)
  new.lm<-NULL
  ifelse(row.length > 1, {
    new.lm<- lm(mean~Year,data=new.sub)
    new.aov<-anova(new.lm)
    new.Pvalue<-new.aov$`Pr(>F)`[1]
    new.Pvalue<-ifelse(is.nan(new.Pvalue),100,new.Pvalue)
    linear.trend<-round(new.lm$coefficients[2],3)
    new.sub$Trend<-linear.trend
    #add categorical for slope (neg or pos)
    new.sub$TrendCat<-ifelse({new.sub$Trend > 0 & new.Pvalue < 0.1},"Positive",
                             ifelse({new.sub$Trend < 0 & new.Pvalue < 0.1},"Negative","Stable"))
    names(new.sub)[names(new.sub) =="Park"]<-"Unit_Code"
    
    
  }, {
    new.sub<-data.frame(Unit_Code=parkName, Year=yearList, mean=NA, SE=NA,AOU_Code=speciesName, Trend=NA, TrendCat="No Data")
  })
  
  
  unique(new.sub$TrendCat)
  
  trend.out<-rbind(trend.out,new.sub)
} 


#add Unit_Name
trend.merge<-merge(trend.out, UnitNames.df.out, by="Unit_Code",all=TRUE,sort=FALSE)

#choose line and fill colors
mycolors<-c("royalblue",I("gray30"),"red")
myfills<-c(alpha("royalblue",0.6),alpha(I("gray30"),0.6),alpha("red",0.6))

#reorder factor levels of TrendCat
trend.merge$TrendCat<-factor(trend.merge$TrendCat, levels=c("Positive", "Stable","Negative"))

#make year numeric
trend.merge$Year<-as.numeric(as.character(trend.merge$Year))


levels(trend.merge$TrendCat)

#plot bci mean annual BCI with linear trend
bci.trend.plot<-ggplot(data=trend.merge, aes(x=Year, y=mean))+
  #geom_path(aes(color=TrendCat, group=TrendCat), size=0.8, alpha=0.7)+
  geom_point(aes(color=TrendCat))+
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE, color=TrendCat), width=0)+
  #geom_smooth(method="lm",aes(color=TrendCat, fill=TrendCat), se = TRUE, size=0.7)+
  #scale_color_viridis(option="viridis", begin=0.2, end=1, direction=-1)+
  #scale_color_gradientn(colors=colfunc(10))+
  #scale_color_gradient2(low="red",mid="black", high="royalblue", midpoint=0, limits = c(-0.2, 0.2))+
  scale_color_manual(name="Trend",values=mycolors)+
  scale_fill_manual(name="Trend",values=myfills)+
  labs(x="Year", y="Mean bird community index (BCI)")+
  theme(panel.background=element_rect(fill="transparent"), panel.border=element_rect(fill=NA,color="black"),
        axis.text.x=element_text(angle=60, hjust=1,size=14),
        axis.text.y=element_text(size=14),
        axis.title = element_text(size=16))+
  theme(legend.justification = 'left', legend.position="top")+
  ggtitle("Trends in Annual Bird Community Index")+
  scale_x_continuous(breaks=c(2007, 2009, 2011, 2013, 2015, 2017))

#create list of trends at each park
trend.df<-unique(trend.out[,c("Park","Trend","TrendCat")])

#bci.trend.plot
bci.trend.plot.out<-bci.trend.plot+
  facet_wrap(~Unit_Name, labeller = labeller(Unit_Name=label_wrap_gen(10)))+
  theme(strip.background =element_rect(fill=I("gray95")))

# geom_text(
#   data    = trend.df,
#   mapping = aes(x = -Inf, y = -Inf, label = Trend, color=TrendCat),
#   hjust   = -0.3,
#   vjust   = -0.7,
#   show.legend = F
#)

bci.trend.plot.out

# #save bci plot
# ggsave(bci.trend.plot.out, file="./Analyses/Results/BCI/BCI_Trend_Fig_8-03-20.png", dpi=600, width=6.5, height=7.5)

##################################################################################################################
# make plots with background colors to show BCI integrity ranges

#create background data.frame
background <- data.frame( lower = c( -Inf, 40.1, 52.1, 60.1), 
                          upper = c( 40.1, 52.1, 60.1, Inf), 
                          Integrity = c("Low Integrity", "Medium Integrity","High Integrity","Highest Integrity"),
                          order=c(1,2,3,4))
#reorder Integrity factor
background.order<-background[order(background$lower)]
bg.levels <- background.order$Integrity
background.order$order <- factor(background.order$Integrity, levels = bg.levels)
#head(background.order)
####################################################################################################

#factor levels of of Network units
trend.merge$Unit_Name<-factor(trend.merge$Unit_Name, levels=unique(trend.merge$Unit_Name))

#plot bci mean annual BCI with linear trend
bci.trend.plot<-ggplot()+
  geom_rect(data = background.order, mapping= aes(ymin = lower, ymax = upper , xmin = -Inf , xmax = Inf, fill = order), alpha = 0.6)+
  #geom_path(aes(color=TrendCat, group=TrendCat), size=0.8, alpha=0.7)+
  # geom_point(data=trend.merge, aes(x=Year, y=mean, color=TrendCat))+
  # geom_errorbar(data=trend.merge, aes(x=Year, y=mean, ymin=mean-SE, ymax=mean+SE, color=TrendCat), width=0)+
  # geom_smooth(data=trend.merge, method="lm",aes(x=Year, y=mean, color=TrendCat), se = TRUE, size=0.7, show.legend=FALSE, alpha= 0.7)+
  geom_path(data=trend.merge, aes(x=Year, y=mean),color="darkgreen", size=0.8, alpha=0.7)+
  geom_point(data=trend.merge, aes(x=Year, y=mean), size=0.8, color="darkgreen")+
  geom_errorbar(data=trend.merge, aes(x=Year, y=mean, ymin=mean-SE, ymax=mean+SE),color="darkgreen", width=0)+
  #geom_smooth(data=subset(trend.merge, TrendCat !="Stable"), method="lm",aes(x=Year, y=mean), color="darkgreen", fill="darkgreen", se = TRUE, size=0.7, show.legend=FALSE, alpha= 0.4)+
  #scale_color_viridis(option="viridis", begin=0.2, end=1, direction=-1)+
  #scale_color_gradientn(colors=colfunc(10))+
  #scale_color_gradient2(low="red",mid="black", high="royalblue", midpoint=0, limits = c(-0.2, 0.2))+
  #scale_color_manual(name="Trend",values=mycolors )+
  #scale_fill_manual(name="Trend",values=myfills)+
  scale_fill_brewer(palette = "Greys")+
  labs(x="Year", y="Mean bird community index (BCI)", fill="Ecological Integrity")+
  theme(panel.background=element_rect(fill="transparent"), panel.border=element_rect(fill=NA,color="black"),
        axis.text.x=element_text(angle=60, hjust=1,size=12),
        axis.text.y=element_text(size=12),
        axis.title = element_text(size=14),
        strip.text.x=element_text(size=14))+
  #theme(legend.justification = 'left', legend.position="top")+
  theme(legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="black"),
        legend.position = c(0.7, 0.07), legend.box = "horizontal")+
  #guides(fill = guide_legend(title.position = "top"))+
  ylim(30,65)+
  #ggtitle("Trends in Annual Bird Community Index")+
  scale_x_continuous(breaks=c(2007, 2009, 2011, 2013, 2015, 2017))+
  guides(fill = guide_legend(reverse = TRUE))

#create list of trends at each park
trend.df<-unique(trend.merge[,c("Park","Trend","TrendCat")])

#bci.trend.plot
bci.trend.plot.out<-bci.trend.plot+
  facet_wrap(~Unit_Name, labeller = labeller(Unit_Name=label_wrap_gen(20)))+
  theme(strip.background =element_rect(fill=I("gray95")), strip.text.x=element_text(size=10))

# geom_text(
#   data    = trend.df,
#   mapping = aes(x = -Inf, y = -Inf, label = Trend, color=TrendCat),
#   hjust   = -0.3,
#   vjust   = -0.7,
#   show.legend = F
#)

bci.trend.plot.out

#save bci plot
ggsave(bci.trend.plot.out, file="./Analyses/Results/BCI/BCI_Trend_Fig_8-03-20.png", dpi=600, width=8.5, height=10)


#################################################################################################
######################################################################################################
#Appendix A - Park Unit 1-pg resource brief summaries

######################################################################################################
#Appendix B - Estimated Occupancy for each species per year in all parks (ADD Network-wide)

#read in formatted data
#data<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Data/Bird_Data/NCRN_Forest_Bird_Data_Formatted.csv", header=TRUE)
data<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Data/Bird_Data/NCRN_Forest_Bird_Data_Formatted_noMigrants.csv",header=TRUE)

#create UnitNames.df
UnitNames.df<-unique(data[,c("Unit_Code","Unit_Name")])

row.names(UnitNames.df)<-NULL
#order by Unit_Name
UnitNames.df$Unit_Name<-trimws(as.character(UnitNames.df$Unit_Name))
UnitNames.df.order<-UnitNames.df[order(UnitNames.df$Unit_Name, decreasing=FALSE),]

#add NCRN
UnitNames.df.out<-rbind(data.frame(Unit_Code="NCRN",Unit_Name="National Capital Region Network"), UnitNames.df.order)

common.names<-unique(data[,c("AOU_Code","Common_Name")])

#read in occupancy results
occu.results<-unique(read.csv(file="/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Analyses/Results/Occupancy/NCRN_occupancy_Forest_POccu_AIC_0-100m_06-24-20.csv"))

#species to remove
removeList<-c("UNBI","UNCH","UNHA","UNCR","UNDU","UNFL","UNOW","UNSP","UNSW","UNTH","UNWA","UNWO","UNWR")

occu.results<-subset(occu.results, !AOU_Code %in% removeList)

occu.sub.1<-subset(occu.results, Metric=="Occupancy")

#subset results to only keep models with GOF_chat > 0.2 and < 4 
occu.sub<-subset(occu.sub.1, c(GOF_chat > 0.2 & GOF_chat < 4))

length(unique(occu.sub$AOU_Code)) #75

#create Predicted mean(SE) column
occu.sub$pred.mean.SE<-paste(round(occu.sub$Predicted,2), " (", round(occu.sub$SE,2), ")",sep="")

#create Naive mean(SE) column
occu.sub$naive.mean.SE<-paste(round(occu.sub$Naive_mean,2), " (", round(occu.sub$Naive_SE,2), ")",sep="")

#now create a naive.pred column
occu.sub$naive.pred<-paste(occu.sub$naive.mean.SE, "; ",occu.sub$pred.mean.SE,sep="")

#round Predicted
#occu.sub$Predicted<-round(occu.sub$Predicted,2)

#add common names
occu.merge<-merge(occu.sub, common.names, by="AOU_Code",all.x=TRUE)

#add Unit Names
colnames(UnitNames.df.out)<-c("Unit_Code","Unit_Name")
occu.merge<-merge(occu.merge, UnitNames.df.out, by="Unit_Code",all.x=TRUE)

#create table
occu.df<-unique(occu.merge[,c("Common_Name","AOU_Code","Unit_Name","Unit_Code","Year","naive.pred")])

#create AOU.Unit.Year column
occu.df$AOU.Unit<-paste(occu.df$Common_Name, occu.df$AOU_Code, occu.df$Unit_Name, occu.df$Unit_Code, sep=".")

#melt and cast
occu.melt<-melt(occu.df, id.vars=c("AOU.Unit","Year"), measure.vars=c("naive.pred"))
occu.cast<-as.data.frame(dcast(occu.melt, AOU.Unit ~ Year, add.missing=FALSE, fun.aggregate=max))

occu.covs<-data.frame(as.character(occu.cast$AOU.Unit))
colnames(occu.covs)<-"AOU.Unit"

#get covs (Error: EOF within quoted stringnumber of items read is not a multiple of the number of columns)
covs<-data.frame(read.delim(text=as.character(occu.covs$AOU.Unit), sep=".", header=FALSE))
colnames(covs)<-c("Common_Name","AOU_Code","Unit_Name","Unit_Code")

#create data.frame
occu.table<-cbind(covs, occu.cast[,-1])

#remove Name columns
occu.table$Common_Name<-NULL
occu.table$Unit_Name<-NULL

#set column names
colnames(occu.table)[1:2]<-c("AOU_Code","Unit_Code")

#save occu.table
write.csv(occu.table, file="./Analyses/Tables/Appendix_C_Occupancy_by_species_per_park_07-01-20.csv")
#########################################################################################################
  #Occupancy Figures
  
  #clear environment (start fresh)
  rm(list=ls())
  
  #load packages
  library(ggplot2)
  library(unmarked)
  library(reshape2)
  library(AICcmodavg)
  library(plyr)
  library(magrittr)
  library(dplyr)
  
  #setwd()
  setwd("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN")
  
  #load source functions
  source("./Analyses/Source/summaryFunction.R")
  
occu.results<-unique(read.csv(file="./Analyses/Results/Occupancy/NCRN_occupancy_Forest_POccu_AIC_0-100m_9-16-20.csv",header=TRUE))

#read in formatted and filtered data
data<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Data/Bird_Data/NCRN_Forest_Bird_Data_Formatted_noMigrants.csv",header=TRUE)


#################################################################################################################  
  #Step 2) Filter data by including only species within parks where > 5 raw detections per park per year.
  
  # #get true list of birds to include in summary based on rule of minimum required detections to run Occupancy models
  # dataRaw<- read.csv("./Data/Bird_Data/NCRN_Forest_Bird_Data_Formatted.csv",header=TRUE)
  # 
  # #create species.unit column
  # dataRaw$Park.Species.Year<-paste(dataRaw$Unit_Code, dataRaw$AOU_Code, dataRaw$Year, sep=".")
  # 
  # #get frequency of counts by species
  # dataFreq<-as.data.frame(table(dataRaw$Park.Species.Year, dataRaw$CountOfAOU_Code))
  # colnames(dataFreq)<-c("Park.Species.Year","row.count","total")
  # range(dataFreq$total)
  # 
  # dataFreq.sub<-subset(dataFreq, total>5)
  # 
  # speciesKeepList<-sort(unique(as.character(dataFreq.sub$unit.habitat.species.year)))
  #length(speciesKeepList)
  
  #filter abun by speciesKeepList
  # occu.all.sub<-subset(occu.all, unit.habitat.species.year %in% speciesKeepList)
  # range(occu.all.sub$Predicted, na.rm=TRUE)
#################################################################################################################  
  
# names(occu.results)
# sort(unique(occu.results$AOU_Code))
# 
# #species to remove
# removeList<-c("UNBI","UNCH","UNHA","UNCR","UNDU","UNFL","UNOW","UNSP","UNSW","UNTH","UNWA","UNWO","UNWR")
# 
# occu.results.sub<-subset(occu.results, !AOU_Code %in% removeList)
# sort(unique(occu.results.sub$AOU_Code))

occu.sub.1<-subset(occu.results, Metric=="Occupancy")

#also look at detection here
occu.sub.det<-subset(occu.results, Metric=="Detection")


#subset results to only keep models with GOF_chat > 0.2 and < 4 
occu.sub<-subset(occu.sub.1, c(GOF_chat > 0.2 & GOF_chat < 4))

#Get Park-level annual mean Occupancy and trend

#list of species
speciesList<-sort(unique(occu.sub$AOU_Code))
#how many unique species?
length(speciesList) #75 of 157 total

#list of years
yearList<-seq(2007,2017,by=1)

#loop to get trend for each species in NCRN and park-level
species.save<-list()
for(i in 1:length(speciesList)){
  
  #subset data by species
  occu.data<-subset(occu.sub, AOU_Code==speciesList[i])
  
  speciesName<-unique(as.character(occu.data$AOU_Code))
  print(speciesName)
  
  #combine annual mean predicted Occupancy for network and parks
  occu.summary.network.pred<-summaryFunction(dataIn=occu.data, response="Predicted", factor="Year")
  occu.summary.network.pred<-cbind("NCRN",occu.summary.network.pred)
  colnames(occu.summary.network.pred)[1]<-c("Park")
  occu.summary.network.pred<-occu.summary.network.pred[,c("Park","Year","mean","SE")]
  #add Estimate column
  occu.summary.network.pred$Estimate<-"Predicted"
  

  #get mean annual predicted occupancy by park 
  occu.summary.park.pred<- unique(occu.data[,c("Unit_Code","Year","Predicted","SE")])
  colnames(occu.summary.park.pred)<-c("Park","Year","mean","SE")
  #add Estimate column
  occu.summary.park.pred$Estimate<-"Predicted"
  
  #get naive annual mean occupancy for network and parks
  occu.summary.network.naive<-summaryFunction(dataIn=occu.data, response="Naive_mean", factor="Year")
  occu.summary.network.naive<-cbind("NCRN",occu.summary.network.naive)
  colnames(occu.summary.network.naive)[1]<-c("Park")
  occu.summary.network.naive<-occu.summary.network.naive[,c("Park","Year","mean","SE")]
  #add Estimate column
  occu.summary.network.naive$Estimate<-"Naive"
  
  #get mean annual predicted occupancy by park 
  occu.summary.park.naive<- unique(occu.data[,c("Unit_Code","Year","Naive_mean","Naive_SE")])
  colnames(occu.summary.park.naive)<-c("Park","Year","mean","SE")
  #add Estimate column
  occu.summary.park.naive$Estimate<-"Naive"
  
  #combine network and park mean annual occupancy
  occu.summary.all<-rbind(rbind(rbind(occu.summary.network.pred, occu.summary.park.pred),occu.summary.network.naive),occu.summary.park.naive)
  
  #add speciesName
  occu.summary.all$AOU_Code<-speciesName
  
  #get list of parks
  parkList<-c("NCRN",unique(sort(as.character(occu.sub$Unit_Code))))
  
  #get estimateList
  estimateList<-unique(occu.summary.all$Estimate)
  
  estimate.out<-list()
  for(y in 1:length(estimateList)){
    
    estimate.data<-subset(occu.summary.all, Estimate==estimateList[y])
    print(estimateList[y])
    
      trend.out<-list()
      for(j in 1:length(parkList)){
        
        new.data<-estimate.data
        new.sub<-subset(new.data, Park==parkList[j])
        
        
        parkName<-as.character(parkList[j])
        print(parkName)
        
        ifelse(parkName =="NCRN",gof.sub<- data.frame(GOF_chat=NA),
               gof.sub<-subset(occu.sub, c(AOU_Code==speciesName & Unit_Code==parkName)))
        
        gof.value<-unique(gof.sub$GOF_chat)
        
        
        row.length<-length(row.names(new.sub))
        
        new.lm<-NULL
        ifelse(row.length > 1, {
          new.lm<- lm(mean~Year,data=new.sub)
          new.aov<-anova(new.lm)
          new.Pvalue<-new.aov$`Pr(>F)`[1]
          new.Pvalue<-ifelse(is.nan(new.Pvalue),100,new.Pvalue)
          linear.trend<-round(new.lm$coefficients[2],3)
          new.sub$Trend<-linear.trend
          #add categorical for slope (neg or pos)
          new.sub$TrendCat<-ifelse({new.sub$Trend > 0 & new.Pvalue < 0.1},"Positive",
                                   ifelse({new.sub$Trend < 0 & new.Pvalue < 0.1},"Negative","Stable"))
          new.sub$Rsq<-summary(new.lm)[[8]]
          #add GOF test results (c-hat)
          new.sub$GOF_chat<-gof.value
        }, {
          new.sub<-data.frame(Park=parkName, Year=yearList, mean=NA, SE=NA,Estimate=estimateList[y],AOU_Code=speciesName, Trend=NA, TrendCat="No Data", Rsq=NA, GOF_chat=NA)
        })
        
        trend.out<-rbind(trend.out,new.sub)
      } 
      
      estimate.out<-rbind(estimate.out, trend.out)
  }
  
  species.save<-rbind(species.save, estimate.out)
}
  
#get the percent declining for each Park

#create a Park.Trend column
percent.df<-subset(species.save, Estimate=="Predicted")
percent.df$Park.Trend<-paste(percent.df$Park, percent.df$TrendCat, sep=".")
percent.df$Park.Trend<-as.factor(as.character(percent.df$Park.Trend))

#remove Year and get unique rows
percent.df.2<-unique(percent.df[,c("Park","AOU_Code","Trend","TrendCat","Park.Trend")])

#get totals
percent.df.summary<-aggregate(percent.df.2$TrendCat~percent.df.2$Park.Trend, FUN="length")
colnames(percent.df.summary)<-c("Park.Trend","Count")
#get covs
percent.covs<-read.table(text=as.character(percent.df.summary$Park.Trend),sep=".")
#combine
percent.df.summary.2<-cbind(percent.covs, percent.df.summary)
colnames(percent.df.summary.2)<-c("Unit_Code","TrendCat","Park.Trend","Count")

#first remove no data
percent.df.summary.3<-subset(percent.df.summary.2, TrendCat !="No Data")

#get proportion from each

#get unitList
unitList<-unique(as.character(percent.df.summary.3$Unit_Code))

#get percents for each park (loop through)
save.pct<-list()
for(i in 1:length(unitList)){
  
  #subset data
  new.perc<-subset(percent.df.summary.3, Unit_Code==unitList[i])
  #get total
  total<-sum(as.integer(new.perc$Count))
  
  #add proportion and percent
  new.perc$Proportion<-new.perc$Count/total
  new.perc$Percent<-round(new.perc$Proportion*100,0)
  
  #check percent=100 with colSums
  colSums(new.perc[,4:6])
  
  #compile results
  save.pct<-rbind(save.pct, new.perc)
}

#subset just declines
positive.df<-subset(save.pct, TrendCat=="Positive")
stable.df<-subset(save.pct, TrendCat=="Stable")
decline.df<-subset(save.pct, TrendCat=="Negative")

#ggplot

#order by percent of species declining
decline.order<-decline.df[order(decline.df$Percent,decreasing=FALSE),]
decline.order$Unit_Code<-factor(decline.order$Unit_Code, levels=decline.order$Unit_Code)

#plot
decline.plot<-ggplot(data=decline.order, aes(x=Unit_Code, y=Percent))+
  coord_flip()+
  geom_bar(stat="identity")+
  lims(y=c(0,100))
decline.plot


#make PIF species bold (in facet titles)

#read in PIF data
pif.data<-read.csv("./Data/PIF_data/PIF_2017_Global.csv")

pif.data.2<-pif.data[,c("Common.Name","Scientific.Name","Continental.Concern","IUCN.Red.List.2016")]
colnames(pif.data.2)[1]<-"Common_Name"

#species to remove
removeList<-c("UNBI","UNCH","UNHA","UNCR","UNDU","UNFL","UNOW","UNSP","UNSW","UNTH","UNWA","UNWO","UNWR")

#simplify data (get list of unique species)
unique.birds<-unique(data[,c("AOU_Code","Common_Name")])

unique.birds.2<-subset(unique.birds, !AOU_Code %in% removeList)

#now merge with PIF status
birds.merge.2<-merge(unique.birds.2, pif.data.2, by="Common_Name",all.x=TRUE)
birds.merge.2$Scientific.Name<-as.character(birds.merge.2$Scientific.Name)

pif.ncrn<-subset(birds.merge.2, Continental.Concern !="")

#merge with species.save
species.merge<-merge(species.save, birds.merge.2, by="AOU_Code", all.x=TRUE)
species.merge$PIF<-ifelse(species.merge$Continental.Concern!="","PIF","Not_PIF")

#reorder factor levels of TrendCat
species.merge$TrendCat<-factor(species.merge$TrendCat, levels=c("Positive", "Stable","Negative","No Data"))
levels(species.merge$TrendCat)

#add SpeciesLabel column
species.merge$SpeciesLabel<-ifelse(species.merge$PIF=="PIF",
                        paste("bold(",species.merge$AOU_Code,")",sep=""), 
                        as.character(species.merge$AOU_Code))

species.merge$Common_Name_2<-ifelse(species.merge$PIF=="PIF", paste(species.merge$Common_Name,"*",sep=""), paste(species.merge$Common_Name))
species.merge$Common_Name_2<-as.factor(species.merge$Common_Name_2)

#remove mean with NAs
species.merge<-subset(species.merge, ! is.na(mean))

#Split Scientific name into genus and species
genus.species<-read.table(text=as.character(species.merge$Scientific.Name), sep=" ", fill=TRUE)
colnames(genus.species)<-c("Genus","Species")

#add columns to data
species.merge<-cbind(species.merge, genus.species)

#Add family and genus to data
aos.codes<-read.csv("./Data/AOS_Codes/SpeciesList_out.csv")

aos.codes.sub<-aos.codes[,c("ORDER","Family","AOU_Code")]
names(aos.codes.sub)[names(aos.codes.sub)=="ORDER"]<- "Order"

#merge
species.merge.2<-merge(species.merge, aos.codes.sub, by="AOU_Code", all.x=TRUE)

#modify facet labels dynamically
species.merge.2$SpeciesLabel <- factor(species.merge.2$SpeciesLabel, levels=unique(species.merge.2$SpeciesLabel), labels = unique(species.merge.2$SpeciesLabel))
unique(species.merge.2$SpeciesLabel)

species.merge.2<-subset(species.merge.2, ! is.na(Common_Name_2))

#add Unit Names
UnitNames.df<-data.frame(Park = c("NCRN","ANTI","CATO","CHOH","GREE","GWMP","HAFE","MANA","MONO","NACE","PISC_FOWA","PRWI","ROCR","WOTR"), Unit_Name = c("National Capital Region Network","Antietam National Battlefield", "Catoctin Mountain Park","Chesapeake and Ohio Canal National Historical Park","Greenbelt National Park","George Washington Memorial Parkway","Harpers Ferry National Historical Park","Manassas National Battlefield Park","Monocacy National Battlefield Park","National Capital Parks-East","Piscataway_Fort Washington","Prince William Forest Park","Rock Creek Park","Wolf Trap National Park for the Performing Arts"))

#######################################################################################################################
#Create trend figures for each Park (and NCRN)

#renamne "Chuck-will's-widow" to force wrapping 
levels(species.merge.2$Common_Name_2)[levels(species.merge.2$Common_Name_2)=="Chuck-will's-widow*"] <- "Chuck-will's- widow*"

for(i in 1:length(parkList)){
  
  print(parkList[i])
  
#subset data by Park
park.data<-subset(species.merge.2, Park==parkList[i])

park.data<-merge(park.data, UnitNames.df, by="Park", all.x=TRUE, sort=FALSE)

parkName<-unique(as.character(park.data$Park))

parkNameTitle<-unique(as.character(park.data$Unit_Name))

parkNameTitle<-ifelse(parkNameTitle=="Piscataway_Fort Washington","Piscataway/Fort Washington",
                      ifelse(parkNameTitle=="Chesapeake and Ohio Canal National Historical Park", "C&O Canal National Historical Park",parkNameTitle))

park.pct.positive<-ifelse(length(subset(positive.df, Unit_Code==parkList[i])[,6])==0,0,subset(positive.df, Unit_Code==parkList[i])[,6])
park.pct.stable<-ifelse(length(subset(stable.df, Unit_Code==parkList[i])[,6])==0,0, subset(stable.df, Unit_Code==parkList[i])[,6])
park.pct.decline<-ifelse(length(subset(decline.df, Unit_Code==parkList[i])[,6])==0,0,subset(decline.df, Unit_Code==parkList[i])[,6])

park.data$TrendCat<-factor(park.data$TrendCat, levels=c("Positive","Stable","Negative"))

#Assess GOF
gof.df<-unique(park.data[,c("Park","Common_Name_2","Estimate","Rsq","GOF_chat")])
gof.df<-subset(gof.df, Estimate=="Predicted")
gof.df %>%
  group_by(Common_Name_2)

#set up label
ifelse(parkName=="NCRN", gof.df$Label<- paste("Rsq =", round(gof.df$Rsq,4)),
       gof.df$Label<-paste("Rsq=", round(gof.df$Rsq,3),"\n","c-hat =", round(gof.df$GOF_chat,3)))


#pick colors for figures
# mycolors<-c("royalblue",I("gray30"),"red")
# myfills<-c(alpha("royalblue",0.5),alpha(I("gray30"),0.5),alpha("red",0.5))

mycolors<-c("gray60","darkgreen")
myfills<-c(alpha("gray60",0.5),alpha("darkgreen",0.5))

park.data$Estimate<-as.factor(park.data$Estimate)

#order by Year
park.data<-park.data[order(park.data$Year,decreasing = FALSE),]

#plot bci mean annual occupancy with linear trend
occu.trend.plot<-ggplot(data=park.data, aes(x=Year, y=mean))+
  #geom_line(aes(color=TrendCat),size=0.4,alpha=0.2)+
  #geom_path(aes(color=TrendCat, group=TrendCat),size=0.3, alpha=0.7)+
  geom_path(aes(color=Estimate),size=0.3, alpha=0.4)+
  #geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE, color=TrendCat), width=0, size=0.3)+
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE,color=Estimate), width=0, size=0.3)+
  #geom_point(aes(color=TrendCat),size=0.8)+
  geom_point(size=0.2,aes(color=Estimate))+
  #geom_smooth(method="lm",aes(color=TrendCat,fill=TrendCat), se = TRUE,size=0.7)+
  geom_smooth(data=subset(park.data, c(TrendCat !="Stable" & Estimate=="Predicted")),method="lm", se = TRUE,color="darkgreen", fill=alpha("darkgreen",0.5),size=0.5)+
  #scale_color_viridis(option="viridis", begin=0.2, end=1, direction=-1)+
  #scale_color_gradientn(colors=colfunc(10))+
  #scale_color_gradient2(low="red",mid="black", high="royalblue", midpoint=0, limits = c(-0.2, 0.2))+
  scale_color_manual(name="Trend",values=mycolors,drop=FALSE)+
  scale_fill_manual(name="Trend",values=myfills,drop=FALSE)+
  labs(x="Year", y="Probability of occupancy")+
  theme(panel.background=element_rect(fill="transparent"), panel.border=element_rect(fill=NA,color="black"),
        axis.text.x=element_text(angle=60, hjust = 1, size=10),plot.title=element_text(size=18),axis.text.y=element_text(size=10),
        strip.text.x=element_text(size=7, face="bold"))+
  scale_x_continuous(breaks=c(2007, 2012, 2017))+
  scale_y_continuous(breaks=seq(0,1,length.out=3), expand=c(0.1,0.1), limits=c(-50,100))+
  coord_cartesian(xlim=c(2007,2017), ylim=c(0,1))+
  #theme(legend.justification = 'left', legend.position="top")+
  theme(legend.position="none")+
  ggtitle(paste(parkNameTitle, " - Trends in Annual Occupancy\n",paste("(",park.pct.positive,"% increasing, ", park.pct.stable, "% stable, and ", park.pct.decline, "% declining)",sep=""),sep=""))
occu.trend.plot

#create new facet labels uisng 
new.lab <- as_labeller(park.data$SpeciesLabel, label_parsed)

#create list of trends at each park
trend.df<-unique(park.data[,c("Common_Name_2","Trend","TrendCat")])
colnames(trend.df)<-c("Common_Name_2","Trend","TrendCat")
length(row.names(trend.df))
levels(trend.df$TrendCat)
# trend.df$Trend<-as.character(trend.df$Trend)
# trend.df$TrendCat<-as.character(trend.df$TrendCat)



#occu.trend.plot
occu.trend.plot.out<-NULL
occu.trend.plot.out<-occu.trend.plot+
  facet_wrap(~Common_Name_2, labeller = labeller(Common_Name_2=label_wrap_gen(10)))+
  theme(strip.background =element_rect(fill=I("gray95")))
  #geom_text(data=gof.df,aes(x=2014,y=1,label=Label),size=2)
  # geom_text(
  #   data= trend.df,
  #   size=2,
  #   mapping = aes(x = -Inf, y = -Inf, label = Trend, color=TrendCat),
  #   hjust   = -0.2,
  #   vjust   = -0.7
  # )
#occu.trend.plot.out

#save plot
ggsave(occu.trend.plot.out, file=paste("./Analyses/Results/Occupancy/", paste(parkName, "Fig_Occupancy_Trends_9-16-20.png",sep="_"),sep="/"), width=12, height=15, dpi=600)


#For report, only show species with significant trends

#get species list that have significant trends
species.sig.data<-subset(park.data,c(TrendCat !="Stable" & Estimate=="Predicted"))
speciesSigList<-unique(species.sig.data$Common_Name_2)

gof.df.sub<-subset(gof.df, Common_Name_2 %in% speciesSigList)

park.data.sig<-subset(park.data, Common_Name_2 %in% speciesSigList)

occu.trend.plot.sig<-ggplot(data=park.data.sig, aes(x=Year, y=mean))+
  #geom_line(aes(color=TrendCat),size=0.4,alpha=0.2)+
  #geom_path(aes(color=TrendCat, group=TrendCat),size=0.3, alpha=0.7)+
  geom_path(aes(color=Estimate),size=0.5, alpha=0.4)+
  #geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE, color=TrendCat), width=0, size=0.3)+
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE,color=Estimate), width=0, size=0.3)+
  #geom_point(aes(color=TrendCat),size=0.8)+
  geom_point(size=0.6,aes(color=Estimate))+
  #geom_smooth(method="lm",aes(color=TrendCat,fill=TrendCat), se = TRUE,size=0.7)+
  geom_smooth(data=subset(park.data.sig, c(TrendCat !="Stable" & Estimate=="Predicted")),method="lm", se = TRUE,color="darkgreen", fill=alpha("darkgreen",0.5),size=0.7)+
  #scale_color_viridis(option="viridis", begin=0.2, end=1, direction=-1)+
  #scale_color_gradientn(colors=colfunc(10))+
  #scale_color_gradient2(low="red",mid="black", high="royalblue", midpoint=0, limits = c(-0.2, 0.2))+
  scale_color_manual(name="Trend",values=mycolors,drop=FALSE)+
  scale_fill_manual(name="Trend",values=myfills,drop=FALSE)+
  labs(x="Year", y="Probability of occupancy")+
  theme(panel.background=element_rect(fill="transparent"), panel.border=element_rect(fill=NA,color="black"),
        axis.text.x=element_text(angle=60, hjust = 1, size=10),plot.title=element_text(size=18),axis.text.y=element_text(size=10),
        strip.text.x=element_text(size=7, face="bold"))+
  scale_x_continuous(breaks=c(2007, 2012, 2017))+
  scale_y_continuous(breaks=seq(0,1,length.out=3), expand=c(0.1,0.1), limits=c(-50,100))+
  coord_cartesian(xlim=c(2007,2017), ylim=c(0,1))+
  #theme(legend.justification = 'left', legend.position="top")+
  theme(legend.position="none")+
  ggtitle(paste(parkNameTitle, " - Species with Significant Trends in Annual Occupancy\n",sep=""))
occu.trend.plot.sig

#create new facet labels uisng 
new.lab <- as_labeller(park.data$SpeciesLabel, label_parsed)

#create list of trends at each park
trend.df<-unique(park.data[,c("Common_Name_2","Trend","TrendCat")])
colnames(trend.df)<-c("Common_Name_2","Trend","TrendCat")
length(row.names(trend.df))
levels(trend.df$TrendCat)
# trend.df$Trend<-as.character(trend.df$Trend)
# trend.df$TrendCat<-as.character(trend.df$TrendCat)

#occu.trend.plot
occu.trend.plot.out.sig<-NULL
occu.trend.plot.out.sig<-try({
  occu.trend.plot.sig+
  facet_wrap(~Common_Name_2, labeller = labeller(Common_Name_2=label_wrap_gen(10)))+
  theme(strip.background =element_rect(fill=I("gray95")))
  #geom_text(data=gof.df.sub,aes(x=2015,y=1,label=Label),size=2)
})

tryCatch({
  ggsave(occu.trend.plot.out.sig, file=paste("./Analyses/Results/Occupancy/", paste(parkName, "Fig_Occupancy_Trends_Significant_9-16-20.png",sep="_"),sep="/"), width=12, height=10, dpi=600)
},error=function(cond2){
  ggsave(NULL, file=paste("./Analyses/Results/Occupancy/", paste(parkName, "Fig_Occupancy_Trends_Significant_9-16-20.png",sep="_"),sep="/"), width=12, height=10, dpi=600)
})

}
######################################################################################################################
######################################################################################################################
#Appendix C: NEW ABUNDANCE FIGS 8-3-2020

#clear environment (start fresh)
rm(list=ls())

#load packages
library(ggplot2)
library(unmarked)
library(reshape2)
library(AICcmodavg)
library(plyr)
library(magrittr)
library(dplyr)

#setwd()
setwd("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN")

#load source functions
source("./Analyses/Source/summaryFunction.R")

abun.results<-unique(read.csv(file="/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Analyses/Results/Abundance/NCRN_abundance_Forest_PCount_AIC_0-100m_9-16-20.csv",header=TRUE))

#read in formatted and filtered data
data<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Data/Bird_Data/NCRN_Forest_Bird_Data_Formatted_noMigrants.csv",header=TRUE)

#subset results
abun.sub.1<-subset(abun.results, Metric=="Abundance")


#subset results to only keep models with GOF_chat > 0.2 and < 4 
abun.sub<-subset(abun.sub.1, c(GOF_chat > 0.2 & GOF_chat < 4))

length(unique(abun.sub$AOU_Code))

#create density columns
abun.sub$Predicted.Density<-abun.sub$Predicted/3.14
abun.sub$Predicted.Density.SE<-abun.sub$SE/3.14

abun.sub$Naive.Density<-abun.sub$Naive_mean/3.14
abun.sub$Naive.Density.SE<-abun.sub$Naive_SE/3.14

#Get Park-level annual mean Occupancy and trend

#list of species
speciesList<-sort(unique(abun.sub$AOU_Code))
#how many unique species?
length(speciesList) #84 of 157 total

#list of years
yearList<-seq(2007,2017,by=1)

#loop to get trend for each species in NCRN and park-level
species.save<-list()
for(i in 1:length(speciesList)){
  
  #subset data by species
  abun.data<-subset(abun.sub, AOU_Code==speciesList[i])
  
  speciesName<-unique(as.character(abun.data$AOU_Code))
  print(speciesName)
  
  # #combine annual mean predicted Occupancy for network and parks
  # abun.summary.network.pred<-summaryFunction(dataIn=abun.data, response="Predicted.Density", factor="Year")
  # abun.summary.network.pred<-cbind("NCRN",abun.summary.network.pred)
  # colnames(abun.summary.network.pred)[1]<-c("Park")
  # abun.summary.network.pred<-abun.summary.network.pred[,c("Park","Year","mean","SE")]
  # #add Estimate column
  # abun.summary.network.pred$Estimate<-"Predicted"
  
  
  #get mean annual predicted occupancy by park 
  abun.summary.park.pred<- unique(abun.data[,c("Unit_Code","Year","Predicted.Density","Predicted.Density.SE")])
  colnames(abun.summary.park.pred)<-c("Park","Year","mean","SE")
  #add Estimate column
  abun.summary.park.pred$Estimate<-"Predicted"
  
  # #get naive annual mean occupancy for network and parks
  # abun.summary.network.naive<-summaryFunction(dataIn=abun.data, response="Naive.Density", factor="Year")
  # abun.summary.network.naive<-cbind("NCRN",abun.summary.network.naive)
  # colnames(abun.summary.network.naive)[1]<-c("Park")
  # abun.summary.network.naive<-abun.summary.network.naive[,c("Park","Year","mean","SE")]
  # #add Estimate column
  # abun.summary.network.naive$Estimate<-"Naive"
  
  #get mean annual predicted occupancy by park 
  abun.summary.park.naive<- unique(abun.data[,c("Unit_Code","Year","Naive.Density","Naive.Density.SE")])
  colnames(abun.summary.park.naive)<-c("Park","Year","mean","SE")
  #add Estimate column
  abun.summary.park.naive$Estimate<-"Naive"
  
  #combine network and park mean annual occupancy
  #abun.summary.all<-rbind(rbind(rbind(abun.summary.network.pred, abun.summary.park.pred),abun.summary.network.naive),abun.summary.park.naive)
  abun.summary.all<-rbind(abun.summary.park.pred,abun.summary.park.naive)
  
  #add speciesName
  abun.summary.all$AOU_Code<-speciesName
  
  #get list of parks
  parkList<-c("NCRN",unique(sort(as.character(abun.sub$Unit_Code))))
  
  #get estimateList
  estimateList<-unique(abun.summary.all$Estimate)
  
  estimate.out<-list()
  for(y in 1:length(estimateList)){
    
    estimate.data<-subset(abun.summary.all, Estimate==estimateList[y])
    print(estimateList[y])
    
    trend.out<-list()
    for(j in 1:length(parkList)){
      
      new.data<-estimate.data
      new.sub<-subset(new.data, Park==parkList[j])
      
      
      parkName<-as.character(parkList[j])
      print(parkName)
      
      ifelse(parkName =="NCRN",gof.sub<- data.frame(GOF_chat=NA),
             gof.sub<-subset(abun.sub, c(AOU_Code==speciesName & Unit_Code==parkName)))
      
      gof.value<-unique(gof.sub$GOF_chat)
      
      
      row.length<-length(row.names(new.sub))
      
      new.lm<-NULL
      ifelse(row.length > 1, {
        new.lm<- lm(mean~Year,data=new.sub)
        new.aov<-anova(new.lm)
        new.Pvalue<-new.aov$`Pr(>F)`[1]
        new.Pvalue<-ifelse(is.nan(new.Pvalue),100,new.Pvalue)
        linear.trend<-round(new.lm$coefficients[2],3)
        new.sub$Trend<-linear.trend
        #add categorical for slope (neg or pos)
        new.sub$TrendCat<-ifelse({new.sub$Trend > 0 & new.Pvalue < 0.1},"Positive",
                                 ifelse({new.sub$Trend < 0 & new.Pvalue < 0.1},"Negative","Stable"))
        new.sub$Rsq<-summary(new.lm)[[8]]
        #add GOF test results (c-hat)
        new.sub$GOF_chat<-gof.value
      }, {
        new.sub<-data.frame(Park=parkName, Year=yearList, mean=NA, SE=NA,Estimate=estimateList[y],AOU_Code=speciesName, Trend=NA, TrendCat="No Data", Rsq=NA, GOF_chat=NA)
      })
      
      trend.out<-rbind(trend.out,new.sub)
    } 
    
    estimate.out<-rbind(estimate.out, trend.out)
  }
  
  species.save<-rbind(species.save, estimate.out)
}

#get the percent declining for each Park

#create a Park.Trend column
percent.df<-subset(species.save, Estimate=="Predicted")
percent.df$Park.Trend<-paste(percent.df$Park, percent.df$TrendCat, sep=".")
percent.df$Park.Trend<-as.factor(as.character(percent.df$Park.Trend))

#remove Year and get unique rows
percent.df.2<-unique(percent.df[,c("Park","AOU_Code","Trend","TrendCat","Park.Trend")])

#get totals
percent.df.summary<-aggregate(percent.df.2$TrendCat~percent.df.2$Park.Trend, FUN="length")
colnames(percent.df.summary)<-c("Park.Trend","Count")
#get covs
percent.covs<-read.table(text=as.character(percent.df.summary$Park.Trend),sep=".")
#combine
percent.df.summary.2<-cbind(percent.covs, percent.df.summary)
colnames(percent.df.summary.2)<-c("Unit_Code","TrendCat","Park.Trend","Count")

#first remove no data
percent.df.summary.3<-subset(percent.df.summary.2, TrendCat !="No Data")

#get proportion from each

#get unitList
unitList<-unique(as.character(percent.df.summary.3$Unit_Code))

#get percents for each park (loop through)
save.pct<-list()
for(i in 1:length(unitList)){
  
  #subset data
  new.perc<-subset(percent.df.summary.3, Unit_Code==unitList[i])
  #get total
  total<-sum(as.integer(new.perc$Count))
  
  #add proportion and percent
  new.perc$Proportion<-new.perc$Count/total
  new.perc$Percent<-round(new.perc$Proportion*100,0)
  
  #check percent=100 with colSums
  colSums(new.perc[,4:6])
  
  #compile results
  save.pct<-rbind(save.pct, new.perc)
}

#subset just declines
positive.df<-subset(save.pct, TrendCat=="Positive")
stable.df<-subset(save.pct, TrendCat=="Stable")
decline.df<-subset(save.pct, TrendCat=="Negative")

#ggplot

#order by percent of species declining
decline.order<-decline.df[order(decline.df$Percent,decreasing=FALSE),]
decline.order$Unit_Code<-factor(decline.order$Unit_Code, levels=decline.order$Unit_Code)

#plot
decline.plot<-ggplot(data=decline.order, aes(x=Unit_Code, y=Percent))+
  coord_flip()+
  geom_bar(stat="identity")+
  lims(y=c(0,100))
decline.plot


#make PIF species bold (in facet titles)

#read in PIF data
pif.data<-read.csv("./Data/PIF_data/PIF_2017_Global.csv")

pif.data.2<-pif.data[,c("Common.Name","Scientific.Name","Continental.Concern","IUCN.Red.List.2016")]
colnames(pif.data.2)[1]<-"Common_Name"

#species to remove
removeList<-c("UNBI","UNCH","UNHA","UNCR","UNDU","UNFL","UNOW","UNSP","UNSW","UNTH","UNWA","UNWO","UNWR")

#simplify data (get list of unique species)
unique.birds<-unique(data[,c("AOU_Code","Common_Name")])

unique.birds.2<-subset(unique.birds, !AOU_Code %in% removeList)

#now merge with PIF status
birds.merge.2<-merge(unique.birds.2, pif.data.2, by="Common_Name",all.x=TRUE)
birds.merge.2$Scientific.Name<-as.character(birds.merge.2$Scientific.Name)

pif.ncrn<-subset(birds.merge.2, Continental.Concern !="")

#merge with species.save
species.merge<-merge(species.save, birds.merge.2, by="AOU_Code", all.x=TRUE)
species.merge$PIF<-ifelse(species.merge$Continental.Concern!="","PIF","Not_PIF")

#reorder factor levels of TrendCat
species.merge$TrendCat<-factor(species.merge$TrendCat, levels=c("Positive", "Stable","Negative","No Data"))
levels(species.merge$TrendCat)

#add SpeciesLabel column
species.merge$SpeciesLabel<-ifelse(species.merge$PIF=="PIF",
                                   paste("bold(",species.merge$AOU_Code,")",sep=""), 
                                   as.character(species.merge$AOU_Code))

species.merge$Common_Name_2<-ifelse(species.merge$PIF=="PIF", paste(species.merge$Common_Name,"*",sep=""), paste(species.merge$Common_Name))
species.merge$Common_Name_2<-as.factor(species.merge$Common_Name_2)

#remove mean with NAs
species.merge<-subset(species.merge, ! is.na(mean))

#Split Scientific name into genus and species
genus.species<-read.table(text=as.character(species.merge$Scientific.Name), sep=" ",fill=TRUE)
colnames(genus.species)<-c("Genus","Species")

#add columns to data
species.merge<-cbind(species.merge, genus.species)

#Add family and genus to data
aos.codes<-read.csv("./Data/AOS_Codes/SpeciesList_out.csv")

aos.codes.sub<-aos.codes[,c("ORDER","Family","AOU_Code")]
names(aos.codes.sub)[names(aos.codes.sub)=="ORDER"]<- "Order"

#merge
species.merge.2<-merge(species.merge, aos.codes.sub, by="AOU_Code", all.x=TRUE)

#modify facet labels dynamically
species.merge.2$SpeciesLabel <- factor(species.merge.2$SpeciesLabel, levels=unique(species.merge.2$SpeciesLabel), labels = unique(species.merge.2$SpeciesLabel))
unique(species.merge.2$SpeciesLabel)

species.merge.2<-subset(species.merge.2, ! is.na(Common_Name_2))

#add Unit Names
UnitNames.df<-data.frame(Park = c("NCRN","ANTI","CATO","CHOH","GREE","GWMP","HAFE","MANA","MONO","NACE","PISC_FOWA","PRWI","ROCR","WOTR"), Unit_Name = c("National Capital Region Network","Antietam National Battlefield", "Catoctin Mountain Park","Chesapeake and Ohio Canal National Historical Park","Greenbelt National Park","George Washington Memorial Parkway","Harpers Ferry National Historical Park","Manassas National Battlefield Park","Monocacy National Battlefield Park","National Capital Parks-East","Piscataway_Fort Washington","Prince William Forest Park","Rock Creek Park","Wolf Trap National Park for the Performing Arts"))

#######################################################################################################################
#Create trend figures for each Park (and NCRN)

#renamne "Chuck-will's-widow" to force wrapping 
levels(species.merge.2$Common_Name_2)[levels(species.merge.2$Common_Name_2)=="Chuck-will's-widow*"] <- "Chuck-will's- widow*"

for(i in 1:length(parkList)){
  
  print(parkList[i])
  
  #subset data by Park
  park.data<-subset(species.merge.2, Park==parkList[i])
  
  unique(park.data$Year)
  
  park.data<-merge(park.data, UnitNames.df, by="Park", all.x=TRUE, sort=FALSE)
  
  parkName<-unique(as.character(park.data$Park))
  
  parkNameTitle<-unique(as.character(park.data$Unit_Name))
  
  parkNameTitle<-ifelse(parkNameTitle=="Piscataway_Fort Washington","Piscataway/Fort Washington",
                        ifelse(parkNameTitle=="Chesapeake and Ohio Canal National Historical Park", "C&O Canal National Historical Park",parkNameTitle))
  
  park.pct.positive<-ifelse(length(subset(positive.df, Unit_Code==parkList[i])[,6])==0,0,subset(positive.df, Unit_Code==parkList[i])[,6])
  park.pct.stable<-ifelse(length(subset(stable.df, Unit_Code==parkList[i])[,6])==0,0, subset(stable.df, Unit_Code==parkList[i])[,6])
  park.pct.decline<-ifelse(length(subset(decline.df, Unit_Code==parkList[i])[,6])==0,0,subset(decline.df, Unit_Code==parkList[i])[,6])
  
  park.data$TrendCat<-factor(park.data$TrendCat, levels=c("Positive","Stable","Negative"))
  
  #Assess GOF
  gof.df<-unique(park.data[,c("Park","Common_Name_2","Estimate","Rsq","GOF_chat")])
  gof.df<-subset(gof.df, Estimate=="Predicted")
  gof.df %>%
    group_by(Common_Name_2)
  
  #set up label
  ifelse(parkName=="NCRN", gof.df$Label<- paste("Rsq =", round(gof.df$Rsq,4)),
         gof.df$Label<-paste("Rsq=", round(gof.df$Rsq,3),"\n","c-hat =", round(gof.df$GOF_chat,3)))
  
  
  #pick colors for figures
  # mycolors<-c("royalblue",I("gray30"),"red")
  # myfills<-c(alpha("royalblue",0.5),alpha(I("gray30"),0.5),alpha("red",0.5))
  
  mycolors<-c("gray60","darkgreen")
  myfills<-c(alpha("gray60",0.5),alpha("darkgreen",0.5))
  
  park.data$Estimate<-as.factor(park.data$Estimate)
  
  #order by Year
  park.data<-park.data[order(park.data$Year,decreasing = FALSE),]
  
  #plot bci mean annual occupancy with linear trend
  abun.trend.plot<-ggplot(data=park.data, aes(x=Year, y=mean))+
    #geom_line(aes(color=TrendCat),size=0.4,alpha=0.2)+
    #geom_path(aes(color=TrendCat, group=TrendCat),size=0.3, alpha=0.7)+
    geom_path(aes(color=Estimate),size=0.3, alpha=0.4)+
    #geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE, color=TrendCat), width=0, size=0.3)+
    geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE,color=Estimate), width=0, size=0.3)+
    #geom_point(aes(color=TrendCat),size=0.8)+
    geom_point(size=0.2,aes(color=Estimate))+
    #geom_smooth(method="lm",aes(color=TrendCat,fill=TrendCat), se = TRUE,size=0.7)+
    geom_smooth(data=subset(park.data, c(TrendCat !="Stable" & Estimate=="Predicted")),method="lm", se = TRUE,color="darkgreen", fill=alpha("darkgreen",0.5),size=0.5)+
    #scale_color_viridis(option="viridis", begin=0.2, end=1, direction=-1)+
    #scale_color_gradientn(colors=colfunc(10))+
    #scale_color_gradient2(low="red",mid="black", high="royalblue", midpoint=0, limits = c(-0.2, 0.2))+
    scale_color_manual(name="Trend",values=mycolors,drop=FALSE)+
    scale_fill_manual(name="Trend",values=myfills,drop=FALSE)+
    labs(x="Year", y="Density (birds per ha)")+
    theme(panel.background=element_rect(fill="transparent"), panel.border=element_rect(fill=NA,color="black"),
          axis.text.x=element_text(angle=60, hjust = 1, size=10),plot.title=element_text(size=18),axis.text.y=element_text(size=10),
          strip.text.x=element_text(size=7, face="bold"))+
    scale_x_continuous(breaks=c(2007, 2012, 2017))+
    #scale_y_continuous(expand=c(0,0.1))+
    coord_cartesian(xlim=c(2007,2017),ylim=c(0,NA))+
    #theme(legend.justification = 'left', legend.position="top")+
    theme(legend.position="none")+
    ggtitle(paste(parkNameTitle, " - Trends in Annual Density\n",paste("(",park.pct.positive,"% increasing, ", park.pct.stable, "% stable, and ", park.pct.decline, "% declining)",sep=""),sep=""))
  abun.trend.plot
  
  #create new facet labels uisng 
  new.lab <- as_labeller(park.data$SpeciesLabel, label_parsed)
  
  #create list of trends at each park
  trend.df<-unique(park.data[,c("Common_Name_2","Trend","TrendCat")])
  colnames(trend.df)<-c("Common_Name_2","Trend","TrendCat")
  length(row.names(trend.df))
  levels(trend.df$TrendCat)
  # trend.df$Trend<-as.character(trend.df$Trend)
  # trend.df$TrendCat<-as.character(trend.df$TrendCat)
  
  
  
  #abun.trend.plot
  abun.trend.plot.out<-NULL
  abun.trend.plot.out<-abun.trend.plot+
    facet_wrap(~Common_Name_2, labeller = labeller(Common_Name_2=label_wrap_gen(10)),scales="free_y")+
    theme(strip.background =element_rect(fill=I("gray95")))
    #geom_text(data=gof.df,aes(x=2015.5,y=1,label=Label),size=2)
  # geom_text(
  #   data= trend.df,
  #   size=2,
  #   mapping = aes(x = -Inf, y = -Inf, label = Trend, color=TrendCat),
  #   hjust   = -0.2,
  #   vjust   = -0.7
  # )
  #abun.trend.plot.out
  
  #save plot
  ggsave(abun.trend.plot.out, file=paste("./Analyses/Results/Abundance/", paste(parkName, "Fig_Abundance_Trends_9-24-20.png",sep="_"),sep="/"), width=12, height=15, dpi=600)
  
  
  #For report, only show species with significant trends
  
  #get species list that have significant trends
  species.sig.data<-subset(park.data,c(TrendCat !="Stable" & Estimate=="Predicted"))
  speciesSigList<-unique(species.sig.data$Common_Name_2)
  
  gof.df.sub<-subset(gof.df, Common_Name_2 %in% speciesSigList)
  
  park.data.sig<-subset(park.data, Common_Name_2 %in% speciesSigList)
  
  abun.trend.plot.sig<-ggplot(data=park.data.sig, aes(x=Year, y=mean))+
    #geom_line(aes(color=TrendCat),size=0.4,alpha=0.2)+
    #geom_path(aes(color=TrendCat, group=TrendCat),size=0.3, alpha=0.7)+
    geom_path(aes(color=Estimate),size=0.5, alpha=0.4)+
    #geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE, color=TrendCat), width=0, size=0.3)+
    geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE,color=Estimate), width=0, size=0.3)+
    #geom_point(aes(color=TrendCat),size=0.8)+
    geom_point(size=0.6,aes(color=Estimate))+
    #geom_smooth(method="lm",aes(color=TrendCat,fill=TrendCat), se = TRUE,size=0.7)+
    geom_smooth(data=subset(park.data.sig, c(TrendCat !="Stable" & Estimate=="Predicted")),method="lm", se = TRUE,color="darkgreen", fill=alpha("darkgreen",0.5),size=0.7)+
    #scale_color_viridis(option="viridis", begin=0.2, end=1, direction=-1)+
    #scale_color_gradientn(colors=colfunc(10))+
    #scale_color_gradient2(low="red",mid="black", high="royalblue", midpoint=0, limits = c(-0.2, 0.2))+
    scale_color_manual(name="Trend",values=mycolors,drop=FALSE)+
    scale_fill_manual(name="Trend",values=myfills,drop=FALSE)+
    labs(x="Year", y="Density (birds per ha)")+
    theme(panel.background=element_rect(fill="transparent"), panel.border=element_rect(fill=NA,color="black"),
          axis.text.x=element_text(angle=60, hjust = 1, size=10),plot.title=element_text(size=18),axis.text.y=element_text(size=10),
          strip.text.x=element_text(size=7, face="bold"))+
    scale_x_continuous(breaks=c(2007, 2012, 2017))+
    #scale_y_continuous(breaks=seq(0,1,length.out=3), expand=c(0.1,0.1))+
    coord_cartesian(xlim=c(2007,2017))+
    #theme(legend.justification = 'left', legend.position="top")+
    theme(legend.position="none")+
    ggtitle(paste(parkNameTitle, " - Species with Significant Trends in Annual Density\n",sep=""))
  abun.trend.plot.sig
  
  #create new facet labels uisng 
  new.lab <- as_labeller(park.data$SpeciesLabel, label_parsed)
  
  #create list of trends at each park
  trend.df<-unique(park.data[,c("Common_Name_2","Trend","TrendCat")])
  colnames(trend.df)<-c("Common_Name_2","Trend","TrendCat")
  length(row.names(trend.df))
  levels(trend.df$TrendCat)
  # trend.df$Trend<-as.character(trend.df$Trend)
  # trend.df$TrendCat<-as.character(trend.df$TrendCat)
  
  #abun.trend.plot
  abun.trend.plot.out.sig<-NULL
  abun.trend.plot.out.sig<-try({
    abun.trend.plot.sig+
      facet_wrap(~Common_Name_2, labeller = labeller(Common_Name_2=label_wrap_gen(10)),scales="free_y")+
      theme(strip.background =element_rect(fill=I("gray95")))
      #geom_text(data=gof.df.sub,aes(x=2015.5,y=1,label=Label),size=2)
  })
  
  tryCatch({
    ggsave(abun.trend.plot.out.sig, file=paste("./Analyses/Results/Abundance/", paste(parkName, "Fig_Abundance_Trends_Significant_9-24-20.png",sep="_"),sep="/"), width=12, height=10, dpi=600)
  },error=function(cond2){
    ggsave(NULL, file=paste("./Analyses/Results/Abundance/", paste(parkName, "Fig_Abundance_Trends_Significant_9-24-20.png",sep="_"),sep="/"), width=12, height=10, dpi=600)
  })
  
}

########################################################################################################
#Appendix D - Table of top models used in occu and abun estimates

#clear environment
rm(list=ls())

#read in data

#occupancy results
occu.results<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Analyses/Results/Occupancy/NCRN_All_occupancy_Forest_POccu_AIC_0-100m_9-24-20.csv",header=TRUE)

occu.sub<-subset(occu.results, Metric=="Occupancy")


#get unique models
occu.unique<-unique(occu.sub[,c("Unit_Code","AOU_Code","Modnames","Cum.Wt","GOF_chat")])
colnames(occu.unique)<-c("Unit_Code","AOU_Code","Occu_Model","Occu_Cum_Wt","Occu_GOF_chat")

#round values
occu.unique$Occu_Cum_Wt<-round(occu.unique$Occu_Cum_Wt,2)
occu.unique$Occu_GOF_chat<-round(occu.unique$Occu_GOF_chat,2)


#abundance results
abun.results<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Analyses/Results/Abundance/NCRN_All_abundance_Forest_PCount_AIC_0-100m_9-16-20.csv",header=TRUE)

abun.sub<-subset(abun.results, Metric=="Abundance")

#get unique models
abun.unique<-unique(abun.sub[,c("Unit_Code","AOU_Code","Modnames","Cum.Wt","GOF_chat")])
colnames(abun.unique)<-c("Unit_Code","AOU_Code","Abun_Model","Abun_Cum_Wt","Abun_GOF_chat")

#round values
abun.unique$Abun_Cum_Wt<-round(abun.unique$Abun_Cum_Wt,2)
abun.unique$Abun_GOF_chat<-round(abun.unique$Abun_GOF_chat,2)

#merge Occupancy and Abundance results
merge.out<-merge(occu.unique, abun.unique, by=c("Unit_Code","AOU_Code"),all=TRUE)

#export table
write.csv(merge.out, file="/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Analyses/Results/Model_Tables/NCRN_model_tables_9_24_20.csv",row.names=FALSE)




##########################################################################################################
#Appendix E - Count of unique species per plot

#read in formatted data
data<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Data/Bird_Data/NCRN_Forest_Bird_Data_Formatted.csv", header=TRUE)

#species to remove
removeList<-c("UNBI","UNCH","UNHA","UNCR","UNDU","UNFL","UNOW","UNSP","UNSW","UNTH","UNWA","UNWO","UNWR")

data<-subset(data, !AOU_Code %in% removeList)

#create Year.Plot column
data$Unit.Plot.Year<-paste(data$Unit_Name, data$Unit_Code, data$Plot_Name, data$Year,    sep=".")

#reduce dataframe
data.reduced<-unique(data[,c("Unit_Code","AOU_Code","Survey_Type","Admin_Unit_Code","Plot_Name","Year","Unit.Plot.Year")])

#create table (All spp)
spp.all.count<-aggregate(data.reduced$AOU_Code ~ data.reduced$Unit.Plot.Year, FUN=length)
colnames(spp.all.count)<-c("Unit.Plot.Year","CountOfSpecies")

#get covs
table.covs<-read.table(text=as.character(spp.all.count$Unit.Plot.Year), sep=".")
colnames(table.covs)<-c("Unit_Name","Unit_Code","Plot_Name","Year")

#create data frame
spp.count.df<-data.frame(table.covs, Total_Species=spp.all.count$CountOfSpecies)

#now melt and cast to get values of years across top of table
spp.count.melt<-melt(spp.count.df, id.vars=c("Unit_Name","Unit_Code","Plot_Name","Year"), measure.vars=c("Total_Species"))
spp.count.cast<-cast(spp.count.melt, Unit_Name+Unit_Code+Plot_Name ~ Year, fun.aggregate = max, fill=NA)

spp.count.cast$Unit_Name<-as.character(spp.count.cast$Unit_Name)

#order by Park Unit
spp.count.order<-spp.count.cast[order(spp.count.cast$Unit_Code, spp.count.cast$Plot_Name),]
row.names(spp.count.order)<-NULL

#set column names
colnames(spp.count.order)[1:3]<-c("Unit_Name","Unit_Code", "Plot_Name")

#write.csv
write.csv(spp.count.order, file="./Analyses/Tables/Appendix_B_Spp_Count_per_Point.csv",row.names=FALSE)

#########################################################################################################
library(scales)
library(viridis)
library(plotrix)

#Make a plot

#create format for ggplot (use all data)
yearList<-colnames(spp.count.order[,-1:-3])
park.spp.long<-list()
for(i in 1:length(yearList)){
  new.stack<-spp.count.order[,c("Unit_Name","Unit_Code","Plot_Name",as.character(yearList[i]))]
  new.stack$Year<-yearList[i]
  names(new.stack)[names(new.stack)==yearList[i]]<-"Species_Count"
  park.spp.long<-rbind(park.spp.long, new.stack)
}


#make plot with geom_tile
ggplot(data=park.spp.long, aes(x=Year, y=Plot_Name, fill = Species_Count)) + 
  geom_tile() +
  scale_fill_viridis() +
  #geom_text(aes(label = round(Species_Count, 0)), size = 2)+
  facet_wrap(~Unit_Code, nrow=3, scales = "free")



#add Unit Names
UnitNames.df<-data.frame(Unit_Code = c("NCRN","ANTI","CATO","CHOH","GREE","GWMP","HAFE","MANA","MONO","NACE","PISC_FOWA","PRWI","ROCR","WOTR"), Unit_Name = c("National Capital Region Network","Antietam National Battlefield", "Catoctin Mountain Park","Chesapeake and Ohio Canal National Historical Park","Greenbelt National Park","George Washington Memorial Parkway","Harpers Ferry National Historical Park","Manassas National Battlefield Park","Monocacy National Battlefield Park","National Capital Parks-East","Piscataway_Fort Washington","Prince William Forest Park","Rock Creek Park","Wolf Trap National Park for the Performing Arts"))

#########################################################################################################
#unitList without NCRN
unitList<- c("ANTI","CATO","CHOH","GREE","GWMP","HAFE","MANA","MONO","NACE", "PISC_FOWA", "PRWI","ROCR","WOTR")     
for(i in 1:length(unitList)){
  
  #get park.name
  park.name<-as.character(unitList[i])
  
  #subset data
  park.spp.sub<-subset(park.spp.long, Unit_Code==park.name)
  
  park.data<-merge(park.spp.sub, UnitNames.df, by="Unit_Code",all.x=TRUE, sort=FALSE)
  park.data$Unit_Name.x<-NULL
  names(park.data)[names(park.data)=="Unit_Name.y"]<-"Unit_Name"
  
  parkName<-unique(as.character(park.data$Unit_Code))
  
  parkNameTitle<-unique(as.character(park.data$Unit_Name))
  
  parkNameTitle<-ifelse(parkNameTitle=="Piscataway_Fort Washington","Piscataway/Fort Washington",
                        ifelse(parkNameTitle=="Chesapeake and Ohio Canal National Historical Park", "C&O Canal National Historical Park",parkNameTitle))
  
  
  
  #get length of rows/2
  split.length<-length(unique(park.data$Plot_Name))/2
  firstHalfList<-as.character(unique(park.data$Plot_Name)[1:split.length])
  secondHalfList<-setdiff(as.character(unique(park.data$Plot_Name) ), as.character(firstHalfList))
  
  if(parkName=="PRWI" | parkName=="CHOH")
  {
    park.data.A<-subset(park.data, Plot_Name %in% firstHalfList)
    
    spp.count.plot<-ggplot(data=park.data.A, aes(x=Year, y=Plot_Name, fill = Species_Count)) + 
      geom_tile(color="black") +
      theme(panel.background = element_rect(fill="transparent"),
            panel.border= element_rect(color="white", fill=NA),
            axis.text.x = element_text(size=12),
            axis.text.y = element_text(size=12),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14),
            title = element_text(size=16),
            legend.key.width=unit(2,"lines"),
            legend.key.height=unit(3,"lines"))+
      scale_fill_viridis(limits=c(0,33)) +
      labs(x="Year", y="Monitoring Location", fill="Number\nof Species")+
      geom_text(aes(label = round(Species_Count, 0)), size = 5, color="white")+
      ggtitle(paste(parkNameTitle, " (Part 1) - Number of species detetected",sep=""))
    
    
    #save plot
    ggsave(spp.count.plot, file=paste("./Analyses/Figures/Species_Richness/",paste(parkName,"part.1.species.count.plot.png",sep="_"),sep=""),width=10.5, height=13, dpi=600)
    
    park.data.B<-subset(park.data, Plot_Name %in% secondHalfList)
    
    spp.count.plot<-ggplot(data=park.data.B, aes(x=Year, y=Plot_Name, fill = Species_Count)) + 
      geom_tile(color="black") +
      theme(panel.background = element_rect(fill="transparent"),
            panel.border= element_rect(color="white", fill=NA),
            axis.text.x = element_text(size=12),
            axis.text.y = element_text(size=12),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14),
            title = element_text(size=16),
            legend.key.width=unit(2,"lines"),
            legend.key.height=unit(3,"lines"))+
      scale_fill_viridis(limits=c(0,33)) +
      labs(x="Year", y="Monitoring Location", fill="Number\nof Species")+
      geom_text(aes(label = round(Species_Count, 0)), size = 5, color="white")+
      ggtitle(paste(parkNameTitle, " (Part 2) - Number of species detetected",sep=""))
    
    ggsave(spp.count.plot, file=paste("./Analyses/Figures/Species_Richness/",paste(parkName,"part.2.species.count.plot.png",sep="_"),sep=""),width=10.5, height=13, dpi=600)
  } else {
    
    
    spp.count.plot<-ggplot(data=park.data, aes(x=Year, y=Plot_Name, fill = Species_Count)) + 
      geom_tile(color="black") +
      theme(panel.background = element_rect(fill="transparent"),
            panel.border= element_rect(color="white", fill=NA),
            axis.text.x = element_text(size=12),
            axis.text.y = element_text(size=12),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14),
            title = element_text(size=16),
            legend.key.width=unit(2,"lines"),
            legend.key.height=unit(3,"lines"))+
      scale_fill_viridis(limits=c(0,33)) +
      labs(x="Year", y="Monitoring Location", fill="Number\nof Species")+
      geom_text(aes(label = round(Species_Count, 0)), size = 5, color="white")+
      ggtitle(paste(parkNameTitle, " - Number of species detetected",sep=""))
    
    ggsave(spp.count.plot, file=paste("./Analyses/Figures/Species_Richness/",paste(parkName,"species.count.plot.png",sep="_"),sep=""),width=10.5, height=13, dpi=600)
  }
}

#########################################################################################################
#Occupancy 10-year means (not in report)

#create species.unit column
occu.sub$Species.Unit<-paste(occu.sub$AOU_Code, occu.sub$Unit_Code, sep="_")

occu.summary<-summaryFunction(dataIn = occu.sub, factor="Species.Unit",response="Predicted")

#separate species and unit
occu.covs<-read.table(text=as.character(occu.summary$Species.Unit),sep="_")
colnames(occu.covs)<-c("AOU_Code","Unit_Code")

occu.summary<-cbind(occu.covs, occu.summary)

#create plot
occuPlot.mean<-ggplot(data=occu.summary, aes(x=Unit_Code, y=mean))+
  #geom_area(aes(fill=AOU_Code),alpha=0.4)+
  #geom_line(stat="identity",aes(color=AOU_Code))+
  #geom_smooth(aes(color=AOU_Code),method="lm",alpha=0.4,size=0.5,se=FALSE)+
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE,color=Unit_Code),width=0)+
  geom_point(aes(color=Unit_Code),alpha=0.8)+
  labs(y="Mean Occupancy", "Park Unit")+
  theme(axis.text.x=element_text(angle=90))+
  ylim(c(0,1))

occuplot.mean.out<-occuPlot.mean+facet_wrap(~AOU_Code,ncol=8)+theme(legend.position="right")

print(occuplot.mean.out)

#########################################################################################################
#Abundance (10-yr means; (not in report))
abun.results<-read.csv(file=paste(getwd(),"Analyses/Abundance/Results","NCRN_abundance_results_Forest.csv",sep="/"))
names(abun.results)
sort(unique(abun.results$AOU_Code))

#remove any estimates above 20 inddividuals
abun.results.reduced<-subset(abun.results, Predicted<20)

#species to remove
removeList<-c("UNBI","UNCH","UNHA","UNCR","UNFL","UNOW","UNSP","UNTH","UNWA","UNWO","UNWR")

abun.results.sub<-subset(abun.results.reduced, !AOU_Code %in% removeList)
sort(unique(abun.results.sub$AOU_Code))


abun.sub<-subset(abun.results.sub, Metric=="Abundance")

#get unitList
unitList<-sort(as.character(unique(abun.sub$Unit_Code)))

#add size of Unit
# Unit_Area.df<-data.frame(Unit_Code=c("ANTI", "CATO", "CHOH", "GWMP", "HAFE", "MANA", "MONO", "NACE", "PRWI", "ROCR", "WOTR"),Area_ha=c(1315, 2490, 7788, 3198, 965, 2064, 667, 4378, 7518, 1100, 53))
# 
# Unit_Area.df<-Unit_Area.df[order(Unit_Area.df$Area_ha,decreasing=TRUE),]
# Unit_Area.df$Unit_Code
# 
# abun.merge<-merge(abun.sub, Unit_Area.df, by="Unit_Code",all.x=TRUE)
# 
# #sort factor levels by Park Unit area (ha)
# abun.merge$Unit_Code<-factor(abun.merge$Unit_Code, levels=as.character(Unit_Area.df$Unit_Code))

#area under the curve plots
unitList<-levels(abun.merge$Unit_Code)
for(i in 1:length(unitList)){
  
  new.unit<-as.character(unitList[i])
  
  abun.sub.unit<-subset(abun.sub, Unit_Code==new.unit)
  
  #make plot
  #create plot
  abunPlot.1<-ggplot(data=abun.sub.unit, aes(x=Year, y=Predicted))+
    geom_area(aes(fill=AOU_Code),alpha=0.4)+
    geom_line(stat="identity",aes(color=AOU_Code))+
    #geom_point(aes(color=AOU_Code, fill=AOU_Code))+
    scale_x_continuous(breaks=c(2007,2009,2011,2013, 2015, 2017),labels=c("2007","2009","2011","2013", "2015", "2017"))+
    labs(y="Predicted Abundance", "Year")+
    ggtitle(new.unit)
  
  abunplot.out<-abunPlot.1+facet_wrap(~AOU_Code,ncol=5)+theme(legend.position="none")
  
  print(abunplot.out)
}


#regressiojn plots
unitList<-levels(abun.merge.forest$Unit_Code)
for(i in 1:length(unitList)){
  
  new.unit<-as.character(unitList[2])
  
  abun.merge.unit<-subset(abun.sub, Unit_Code==new.unit)
  
  #make plot
  #create plot
  abunPlot.trend.1<-ggplot(data=abun.merge.unit, aes(x=Year, y=Predicted))+
    geom_area(aes(fill=AOU_Code),alpha=0.4)+
    geom_line(stat="identity",aes(color=AOU_Code))+
    #geom_smooth(aes(color=AOU_Code),method="lm",alpha=0.4,size=0.5,se=FALSE)+
    geom_errorbar(aes(ymin=Predicted-SE, ymax=Predicted+SE,color=AOU_Code),width=0)+
    geom_point(aes(color=AOU_Code),alpha=0.8)+
    scale_x_continuous(breaks=c(2007,2009,2011,2013,2015, 2017),labels=c("2007","2009","2011","2013", "2015", "2017"))+
    labs(y="Predicted Abundance", "Year")+
    theme(axis.text.x=element_text(angle=90))+
    ylim(c(0,20))+
    ggtitle(new.unit)
  
  abunplot.trend.out<-abunPlot.trend.1+facet_wrap(~AOU_Code,ncol=8)+theme(legend.position="none")
  
  print(abunplot.trend.out)
}
range(abun.sub$Predicted,na.rm=TRUE)

##########################################################################################
#Abundance 10-year means (not in report)

#create species.unit column
abun.sub$Species.Unit<-paste(abun.sub$AOU_Code, abun.sub$Unit_Code, sep="_")

abun.summary<-summaryFunction(dataIn = abun.sub, factor="Species.Unit",response="Predicted")

#separate species and unit
abun.covs<-read.table(text=as.character(abun.summary$Species.Unit),sep="_")
colnames(abun.covs)<-c("AOU_Code","Unit_Code")

#add species names and families back to data.frame
ncrn.data<- read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Data/BCI_Guild_Assignments/BirdGuildAssignments_Appalachian.csv", header=TRUE)

ncrn.species.data<-ncrn.data[,c("Scientific_Name","Common_Name","AOU_Code")]

all.species.data<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Data/BCI_Guild_Assignments/NACC_list_species.csv", header=TRUE)
names(all.species.data)[names(all.species.data)=="species"]<-"Scientific_Name"

bird.taxa<-all.species.data[,c("Scientific_Name", "order","family","genus")]

spp.data.merge<-merge(species.data, bird.taxa, by="Scientific_Name",all.x=TRUE)
head(species.data)
head(bird.taxa)


abun.summary<-cbind(abun.covs, abun.summary)

#create plot
abunPlot.mean<-ggplot(data=abun.summary, aes(x=Unit_Code, y=mean))+
  #geom_area(aes(fill=AOU_Code),alpha=0.4)+
  #geom_line(stat="identity",aes(color=AOU_Code))+
  #geom_smooth(aes(color=AOU_Code),method="lm",alpha=0.4,size=0.5,se=FALSE)+
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE,color=Unit_Code),width=0)+
  geom_point(aes(color=Unit_Code),alpha=0.8)+
  labs(y="Mean Abundance", "Park Unit")+
  theme(axis.text.x=element_text(angle=90))
  #ylim(c(0,))

abunplot.mean.out<-abunPlot.mean+facet_wrap(~AOU_Code,ncol=8)+theme(legend.position="right")

print(abunplot.mean.out)

################################################################################################################
#swicth levels of AOU_Code
abun.summary$AOU_Code<-factor(abun.summary$AOU_Code, levels=rev(sort(unique(abun.summary$AOU_Code))))

#plot with species on y axis
#create plot
abunPlot.mean.2<-ggplot(data=abun.summary, aes(x=AOU_Code, y=mean))+
  coord_flip()+
  #geom_area(aes(fill=AOU_Code),alpha=0.4)+
  #geom_line(stat="identity",aes(color=AOU_Code))+
  #geom_smooth(aes(color=AOU_Code),method="lm",alpha=0.4,size=0.5,se=FALSE)+
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE,color=AOU_Code),width=0)+
  geom_bar(stat="identity",aes(fill=AOU_Code),alpha=0.8)+
  labs(y="Mean Abundance", "Park Unit")+
  theme(axis.text.x=element_text(angle=0))
#ylim(c(0,))

abunplot.mean.out.2<-abunPlot.mean.2+facet_wrap(~Unit_Code,ncol=11)+theme(legend.position="none")

print(abunplot.mean.out.2)

#############################################################################################################
#Older abundance trend figs (no longer in report)
#read in formatted data
data<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Data/Bird_Data/NCRN_Forest_Bird_Data_Formatted.csv", header=TRUE)

common.names<-unique(data[,c("AOU_Code","Common_Name")])

#read in abunpancy results
abun.results<-unique(read.csv(file="/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Analyses/Abundance/Results/NCRN_abundance_results_Forest_PCount.csv"))

#convert all ETTI to TUTI and Eastern Tufted Titmouse to Tufted Titmouse
abun.results$AOU_Code<-gsub("ETTI","TUTI",abun.results$AOU_Code)

#species to remove
removeList<-c("UNBI","UNCH","UNHA","UNCR","UNDU","UNFL","UNOW","UNSP","UNSW","UNTH","UNWA","UNWO","UNWR")

abun.results<-subset(abun.results, !AOU_Code %in% removeList)


abun.sub<-subset(abun.results, Metric=="Abundance")

#replace any values > 30 with NA
abun.sub$Predicted[abun.sub$Predicted > 15] <- NA
abun.sub$SE[abun.sub$SE>10]<-NA
abun.sub$SE[is.na(abun.sub$Predicted)] <- NA


#add Unit Names
UnitNames.df<-data.frame(Unit_Code = c("NCRN","ANTI","CATO","CHOH","GWMP","HAFE","MANA","MONO","NACE","PRWI","ROCR","WOTR"), Unit_Name = c("National Capital Region Network","Antietam National Battlefield Park", "Catoctin Mountain Park","Chesapeake and Ohio Canal National Historical Park","George Washington Memorial Parkway","Harpers Ferry National Historical Park","Manassas National Battlefield Park","Monocacy National Battlefield Park","National Capital Parks-East","Prince William Forest Park","Rock Creek Park","Wolf Trap National Park for the Performing Arts"))

abun.merge<-merge(abun.sub, UnitNames.df, by="Unit_Code",all.x=TRUE)
abun.merge$Unit_Name<-as.character(abun.merge$Unit_Name)

#get unitList
unitList<-sort(as.character(unique(abun.merge$Unit_Name)))

#area under the curve plots
for(i in 1:length(unitList)){
  
  new.unit<-as.character(unitList[i])
  
  #subset by Unit
  abun.sub.unit<-subset(abun.merge, Unit_Name==new.unit)
  
  #get species list
  sppList<-unique(sort(as.character(abun.sub.unit$AOU_Code)))
  
  #get count of unique species
  sppCount<-length(sppList)
  n.rows<-ceiling(sppCount/5)
  
  ifelse(n.rows < 15,
         {
           #make plot
           # abunPlot.1<-ggplot(data=abun.sub.unit, aes(x=Year, y=Predicted))+
           #   geom_area(aes(fill=AOU_Code),alpha=0.4)+
           #   #geom_errorbar(aes(ymin=Predicted-SE, ymax=Predicted+SE, color=AOU_Code), width=0)+
           #   geom_line(stat="identity",aes(color=AOU_Code))+
           #   #geom_point(aes(color=AOU_Code, fill=AOU_Code))+
           # scale_x_continuous(breaks=c(2007,2009,2011,2013, 2015, 2017),labels=c("2007","2009","2011","2013", "2015", "2017"))+
           #     scale_y_continuous(breaks=c(0,1),labels=c("0","1"))+
           #   theme(axis.text.x=element_text(size=8, angle=90),
           #         axis.text.y=element_text(size=8))+
           #   labs(y="Estimated abundance", "Year")+
           #   ggtitle(new.unit)
           # 
           # abunplot.out<-abunPlot.1+facet_wrap(~AOU_Code, ncol=5, nrow=as.integer(n.rows))+theme(legend.position="none")
           # print(abunplot.out)
           
           
           #create plot
           abunPlot.trend.1<-ggplot(data=abun.sub.unit, aes(x=Year, y=Predicted))+
             #geom_area(aes(fill=AOU_Code),alpha=0.4)+
             #geom_line(stat="identity",aes(color=AOU_Code))+
             geom_smooth(aes(fill=AOU_Code),color="white",method="lm",alpha=0.4,size=0.5)+
             geom_errorbar(aes(ymin=Predicted-SE, ymax=Predicted+SE,color=AOU_Code),width=0)+
             geom_point(aes(color=AOU_Code),alpha=0.8)+
             scale_x_continuous(breaks=c(2007,2009,2011,2013,2015, 2017),labels=c("2007","2009","2011","2013", "2015", "2017"))+
             labs(y="Predicted Abundance", "Year")+
             theme(axis.text.x=element_text(angle=90),
                   axis.text.y=element_text(size=8))+
             #ylim(c(0,30))+
             ggtitle(new.unit)
           
           abunplot.trend.out<-abunPlot.trend.1+facet_wrap(~AOU_Code,ncol=5, nrow=n.rows,scales="free_y")+theme(legend.position="none")
           print(abunplot.trend.out)
           
           
         },
         {
           #get species list
           sppList<-unique(sort(as.character(abun.sub.unit$AOU_Code)))
           
           #get count of unique species
           sppCount<-length(sppList)
           
           #divide sppList into 2 fragments
           sppList.1<-sppList[c(seq(1,ceiling(sppCount/2)))]
           #use setdiff to get other part of sppList
           sppList.2<-setdiff(sppList, sppList.1)
           
           #plot first part
           
           #subset data to only inclue sppList.1
           abun.part.1.data<-subset(abun.sub.unit, AOU_Code %in% sppList.1)
           
           #get count of unique species
           sppCount<-length(unique(abun.part.1.data$AOU_Code))
           
           n.rows<-ceiling(sppCount/5)
           
           #make plot
           abunPlot.trend.1<-ggplot(data=abun.part.1.data, aes(x=Year, y=Predicted))+
             #geom_area(aes(fill=AOU_Code),alpha=0.4)+
             #geom_line(stat="identity",aes(color=AOU_Code))+
             geom_smooth(aes(fill=AOU_Code),color="white",method="lm",alpha=0.4,size=0.5)+
             geom_errorbar(aes(ymin=Predicted-SE, ymax=Predicted+SE,color=AOU_Code),width=0)+
             geom_point(aes(color=AOU_Code),alpha=0.8)+
             scale_x_continuous(breaks=c(2007,2009,2011,2013,2015, 2017),labels=c("2007","2009","2011","2013", "2015", "2017"))+
             labs(y="Predicted Abundance", "Year")+
             theme(axis.text.x=element_text(angle=90),
                   axis.text.y=element_text(size=8))+
             #ylim(c(0,30))+
             ggtitle(paste(new.unit, "(part 1)",sep=" "))
           
           abunplot.trend.out.1<-abunPlot.trend.1+facet_wrap(~AOU_Code,ncol=5, nrow=n.rows,scales="free_y")+theme(legend.position="none")
           print(abunplot.trend.out.1)
           
           cat("\n\n\\pagebreak\n")
           
           
           #subset data to only inclue sppList.2
           abun.part.2.data<-subset(abun.sub.unit, AOU_Code %in% sppList.2)
           
           #get count of unique species
           sppCount<-length(unique(abun.part.2.data$AOU_Code))
           #get number of needed rows, if there are 5 columns
           n.rows<-ceiling(sppCount/5)
           
           abunPlot.trend.1<-ggplot(data=abun.part.2.data, aes(x=Year, y=Predicted))+
             #geom_area(aes(fill=AOU_Code),alpha=0.4)+
             #geom_line(stat="identity",aes(color=AOU_Code))+
             geom_smooth(aes(fill=AOU_Code),color="white",method="lm",alpha=0.4,size=0.5)+
             geom_errorbar(aes(ymin=Predicted-SE, ymax=Predicted+SE,color=AOU_Code),width=0)+
             geom_point(aes(color=AOU_Code),alpha=0.8)+
             scale_x_continuous(breaks=c(2007,2009,2011,2013,2015, 2017),labels=c("2007","2009","2011","2013", "2015", "2017"))+
             labs(y="Predicted Abundance", "Year")+
             theme(axis.text.x=element_text(angle=90),
                   axis.text.y=element_text(size=8))+
             #ylim(c(0,30))+
             ggtitle(paste(new.unit, "(part 2)",sep=" "))
           
           abunplot.trend.out.2<-abunPlot.trend.1+facet_wrap(~AOU_Code,ncol=5, nrow=n.rows,scales="free_y")+theme(legend.position="none")
           print(abunplot.trend.out.2)
           
           
         })
  
#############################################################################################################
  #try with sorting species by slope (positive to negative) and color by family (get from .csv file)
  
  #get slopes
  sppList<-unique(as.character(abun.part.2.data$AOU_Code))
  
  mod.1<-lm(Predicted ~ Year+AOU_Code, data=abun.part.2.data)
  
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  beta.df<-as.data.frame(mod.1$coefficients)
  beta.df$AOU_Code<- substrRight(row.names(beta.df),4)
  row.names(beta.df)<-NULL
  colnames(beta.df)[1]<-"Slope"
  
  #remove Year and Intercept
  beta.df<-subset(beta.df, c(AOU_Code != "Year" & AOU_Code != "ept)"))
  
  #order by slope
  beta.order<-beta.df[order(beta.df$Slope,decreasing=FALSE),]
  
  sppOrder<-as.character(beta.order$AOU_Code)
  
  #now try with loop
  
  slopes.out<-list()
  for(i in 1:length(sppList)){
    #subset by species
    new.sp<-subset(abun.part.2.data, AOU_Code==sppList[i])
    #fit model
    new.mod<-NULL
    try(new.mod<-lm(Predicted ~ Year, data=new.sp))
    
    try(new.slope<-new.mod$coefficients[2])
    
    try(new.df<-data.frame(AOU_Code=sppList[i], Slope=new.slope))
    
    slopes.out<-rbind(slopes.out, new.df)
  }
  
  slopes.df<-unique(slopes.out)
  
  #order slopes
  slopes.order<-slopes.df[order(slopes.df$Slope,decreasing=FALSE),]
  
  sppOrder<-as.character(slopes.order$AOU_Code)
  
  
  #reorder AOU code in abun.part.2.data
  abun.sub.unit$AOU_Code<-factor(abun.sub.unit$AOU_Code, levels=rev(sppOrder))
  

  levels(abun.sub.unit.1$AOU_Code)
  abunPlot.trend.1<-ggplot(data=abun.sub.unit.1, aes(x=Year, y=Predicted))+
    #geom_area(aes(fill=AOU_Code),alpha=0.4)+
    #geom_line(stat="identity",aes(color=AOU_Code))+
    geom_smooth(aes(fill=AOU_Code),color="white",method="lm",alpha=0.4,size=0.5)+
    geom_errorbar(aes(ymin=Predicted-SE, ymax=Predicted+SE,color=AOU_Code),width=0)+
    geom_point(aes(color=AOU_Code),alpha=0.8)+
    scale_x_continuous(breaks=c(2007,2009,2011,2013,2015, 2017),labels=c("2007","2009","2011","2013", "2015", "2017"))+
    labs(y="Predicted Abundance", "Year")+
    theme(axis.text.x=element_text(angle=90),
          axis.text.y=element_text(size=8))+
    #ylim(c(0,30))+
    ggtitle(new.unit)
  
  abunplot.trend.out<-abunPlot.trend.1+facet_wrap(~AOU_Code,ncol=5, nrow=n.rows,scales="free_y")+theme(legend.position="none")
  print(abunplot.trend.out)
  
  #reorder AOU code in abun.part.2.data
  abun.part.2.data$AOU_Code<-factor(abun.part.2.data$AOU_Code, levels=rev(sppOrder))
  
  #now plot
  abunPlot.trend.1<-ggplot(data=abun.part.2.data, aes(x=Year, y=Predicted))+
    #geom_area(aes(fill=AOU_Code),alpha=0.4)+
    #geom_line(stat="identity",aes(color=AOU_Code))+
    geom_smooth(aes(fill=AOU_Code),color="white",method="lm",alpha=0.4,size=0.5)+
    geom_errorbar(aes(ymin=Predicted-SE, ymax=Predicted+SE,color=AOU_Code),width=0)+
    geom_point(aes(color=AOU_Code),alpha=0.8)+
    scale_x_continuous(breaks=c(2007,2009,2011,2013,2015, 2017),labels=c("2007","2009","2011","2013", "2015", "2017"))+
    labs(y="Predicted Abundance", "Year")+
    theme(axis.text.x=element_text(angle=90),
          axis.text.y=element_text(size=8))+
    #ylim(c(0,30))+
    ggtitle(paste(new.unit, "(part 2)",sep=" "))
  
  abunplot.trend.out.2<-abunPlot.trend.1+facet_wrap(~AOU_Code,ncol=5, nrow=n.rows,scales="free_y")+theme(legend.position="none")
  print(abunplot.trend.out.2)
  
#########################################################################################################
#Appendix A
#Park Summary Resource briefs

  #load libraries
  library(plyr)
  library(sp)
  library(rgdal)
  library(maptools)
  library(ggmap)
  library(ggplot2)
  library(pander)
  library(ggsn)
  library(png)
  library(grid)
  library(gridExtra)
  library(xtable)
  library(scales)
  library(spatstat)
  library(raster)
  library(akima)
  library(xtable)
  library(dismo)
  
  source("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Analyses/Source/summaryFunction.R")
  
  #subchunkify function (for dynamically modifying plot dimensions from within R chunk)
  subchunkify <- function(g, fig_height=7, fig_width=5) {
    g_deparsed <- paste0(deparse(
      function() {g}
    ), collapse = '')
    
    sub_chunk <- paste0("
                        `","``{r sub_chunk_", floor(runif(1) * 10000), ", fig.height=", fig_height, ", fig.width=", fig_width, ", echo=FALSE}",
                        "\n(", 
                        g_deparsed
                        , ")()",
                        "\n`","``
                        ")
    
    cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
  }
  
  
  #define padding around figures
  if(is_latex_output()) {
    plot_default <- knit_hooks$get("plot")
    knit_hooks$set(plot = function(x, options) { 
      x <- c(plot_default(x, options), "\\vspace{8pt}")
    })
  }
  
  #cat("\n\\newgeometry{top=0.1in,left=0.9in,right=0.9in,bottom=0.2in}\n")
  
  
  #############################################
  #read in raw data
  data<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Data/Bird_Data/NCRN_Forest_Bird_Data_Formatted.csv", header=TRUE)
  
  #############################################
  #species to remove
  removeList<-c("UNBI","UNCH","UNHA","UNCR","UNDU","UNFL","UNOW","UNSP","UNSW","UNTH","UNWA","UNWO","UNWR")
  
  data<-subset(data, !AOU_Code %in% removeList)
  
  #############################################
  #get species richness per point
  data$Park.Plot<-paste(trimws(data$Unit_Name), trimws(data$Unit_Code), trimws(data$Plot_Name), sep=".")
  
  #create table (All spp)
  spp.all.count<-aggregate(data$AOU_Code ~ data$Park.Plot, FUN=function(x){x=length(unique(x))})
  colnames(spp.all.count)<-c("Park.Plot","Count_of_Species")
  
  #get covs
  covs<-data.frame(read.delim(text=as.character(spp.all.count$Park.Plot), sep=".",header=FALSE))
  colnames(covs)<-c("Unit_Name","Unit_Code","Plot_Name")
  sp.richness<-cbind(covs, spp.all.count)
  #############################################
  #Get density (birds/ha) estimates
  
  abun.all<-unique(read.csv(file="/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Analyses/Results/Abundance/NCRN_abundance_Forest_PCount_AIC_0-100m_06-28-20.csv",header=TRUE))
  
  #create Park.Species.Year
  abun.all$Park.Species.Year<-paste(abun.all$Unit_Code, abun.all$AOU_Code, abun.all$Year,sep=".")
  #################################################################################################################  
  #Step 2) Filter data by including only species within parks where > 5 raw detections per park per year.
  
  #get true list of birds to include in summary based on rule of minimum required detections to run abundance models
  dataRaw<- read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Data/Bird_Data/NCRN_Forest_Bird_Data_Formatted.csv",header=TRUE)
  
  #create species.unit column
  dataRaw$Park.Species.Year<-paste(dataRaw$Unit_Code, dataRaw$AOU_Code, dataRaw$Year, sep=".")
  
  #get frequency of counts by species
  dataFreq<-as.data.frame(table(dataRaw$Park.Species.Year, dataRaw$CountOfAOU_Code))
  colnames(dataFreq)<-c("Park.Species.Year","row.count","total")
  
  dataFreq.sub<-subset(dataFreq, total>5)
  
  speciesKeepList<-sort(unique(as.character(dataFreq.sub$Park.Species.Year)))
  #length(speciesKeepList)
  
  #filter abun by speciesKeepList
  #abun.all.sub<-subset(abun.all, Park.Species.Year %in% speciesKeepList)
  # range(abun.all.sub$Predicted, na.rm=TRUE)
  #################################################################################################################  
  
  abun.results<-abun.all
  
  #species to remove
  removeList<-c("UNBI","UNCH","UNHA","UNCR","UNDU","UNFL","UNOW","UNSP","UNSW","UNTH","UNWA","UNWO","UNWR")
  
  abun.results.sub<-subset(abun.results, !AOU_Code %in% removeList)
  
  abun.sub<-subset(abun.results.sub, Metric=="Abundance")
  
  #filter abundance results
  #abun.sub<-subset(abun.sub, Predicted<40)
  
  #create density column
  abun.sub$Density<-abun.sub$Predicted/3.14
  abun.sub$Density.SE<-abun.sub$SE/3.14
  
  #Get Park-level annual mean Abundance and trend
  
  #list of species
  speciesList<-sort(unique(abun.sub$AOU_Code))
  #how many unique species?
  #length(speciesList) #157
  #list of years
  yearList<-seq(2007,2017,by=1)
  
  #loop to get trend for each species in NCRN and park-level
  species.save<-list()
  for(i in 1:length(speciesList)){
    
    #subset data by species
    abun.data<-subset(abun.sub, AOU_Code==speciesList[i])
    
    speciesName<-unique(as.character(abun.data$AOU_Code))
    
    #combine annual mean abundance for network and parks
    abun.summary.network<-summaryFunction(dataIn=abun.data, response="Density", factor="Year")
    abun.summary.network<-cbind("NCRN",abun.summary.network)
    colnames(abun.summary.network)[1]<-c("Park")
    abun.summary.network<-abun.summary.network[,c("Park","Year","mean","SE")]
    
    #get mean annual BCI by park 
    abun.summary.park<- unique(abun.data[,c("Unit_Code","Year","Density","Density.SE")])
    colnames(abun.summary.park)<-c("Park","Year","mean","SE")
    
    #combine network and park mean annual BCI
    abun.summary.all<-rbind(abun.summary.network, abun.summary.park)
    
    #add speciesName
    abun.summary.all$AOU_Code<-speciesName
    
    #get list of parks
    parkList<-c("NCRN",unique(sort(as.character(abun.sub$Unit_Code))))
    
    trend.out<-list()
    for(j in 1:length(parkList)){
      
      new.data<-abun.summary.all
      new.sub<-subset(new.data, Park==parkList[j])
      
      parkName<-as.character(parkList[j])
      
      row.length<-length(row.names(new.sub))
      
      new.lm<-NULL
      ifelse(row.length > 1, {
        new.lm<- lm(mean~Year,data=new.sub)
        linear.trend<-round(new.lm$coefficients[2],3)
        new.sub$Trend<-linear.trend
        #add categorical for slope (neg or pos)
        new.sub$TrendCat<-ifelse(new.sub$Trend > 0.05,"Positive",
                                 ifelse(new.sub$Trend < -0.05,"Negative","Stable"))
        
        
      }, {
        new.sub<-data.frame(Park=parkName, Year=yearList, mean=NA, SE=NA,AOU_Code=speciesName, Trend=NA, TrendCat="No Data")
      })
      
      #add to data
      trend.out<-rbind(trend.out,new.sub)
    } 
    
    species.save<-rbind(species.save, trend.out)
  }
  
  #make PIF species bold (in facet titles)
  
  #read in PIF data
  pif.data<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Data/PIF_data/PIF_2017_Global.csv")
  
  pif.data.2<-pif.data[,c("Common.Name","Scientific.Name","Continental.Concern","IUCN.Red.List.2016")]
  colnames(pif.data.2)[1]<-"Common_Name"
  
  #species to remove
  removeList<-c("UNBI","UNCH","UNHA","UNCR","UNDU","UNFL","UNOW","UNSP","UNSW","UNTH","UNWA","UNWO","UNWR")
  
  #simplify data (get list of unique species)
  unique.birds<-unique(dataRaw[,c("AOU_Code","Common_Name")])
  
  unique.birds.2<-subset(unique.birds, !AOU_Code %in% removeList)
  
  #now merge with PIF status
  birds.merge.2<-merge(unique.birds.2, pif.data.2, by="Common_Name",all.x=TRUE)
  birds.merge.2$Scientific.Name<-as.character(birds.merge.2$Scientific.Name)
  
  pif.ncrn<-subset(birds.merge.2, Continental.Concern !="")
  
  #merge with species.save
  species.merge<-merge(species.save, birds.merge.2, by="AOU_Code", all.x=TRUE)
  species.merge$PIF<-ifelse(species.merge$Continental.Concern!="","PIF","Not_PIF")
  
  #reorder factor levels of TrendCat
  species.merge$TrendCat<-factor(species.merge$TrendCat, levels=c("Positive", "Stable","Negative","No Data"))
  
  
  #add SpeciesLabel column
  species.merge$SpeciesLabel<-ifelse(species.merge$PIF=="PIF",
                                     paste("bold(",species.merge$AOU_Code,")",sep=""), 
                                     as.character(species.merge$AOU_Code))
  
  
  #remove mean with NAs
  species.merge<-subset(species.merge, ! is.na(mean))
  
  #Split Scientific name into genus and species
  genus.species<-read.table(text=as.character(species.merge$Scientific.Name), sep=" ")
  colnames(genus.species)<-c("Genus","Species")
  
  #add columns to data
  species.merge<-cbind(species.merge, genus.species)
  
  #Add family and genus to data
  aos.codes<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Data/AOS_Codes/SpeciesList_out.csv")
  
  aos.codes.sub<-aos.codes[,c("ORDER","Family","AOU_Code")]
  names(aos.codes.sub)[names(aos.codes.sub)=="ORDER"]<- "Order"
  
  #merge
  species.merge.2<-merge(species.merge, aos.codes.sub, by="AOU_Code", all.x=TRUE)
  
  #modify facet labels dynamically
  species.merge.2$SpeciesLabel <- factor(species.merge.2$SpeciesLabel, levels=unique(species.merge.2$SpeciesLabel), labels = unique(species.merge.2$SpeciesLabel))
  
  #############################################
  #calculate BCI
  
  #read in formatted data
  bci.data<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Analyses/Results/BCI/NCRN_BCI_results.csv", header=TRUE)
  
  #now create plot of BCI by year and park, and add linear trendlines
  bci.data$Unit_Code<-bci.data$unitName
  bci.data$Year<-as.factor(as.character(bci.data$yearName))
  
  #add Unit Names
  UnitNames.df<-data.frame(Unit_Code = c("NCRN","ANTI","CATO","CHOH","GREE","GWMP","HAFE","MANA","MONO","NACE","PRWI","PISCA_FOWA","ROCR","WOTR"), Unit_Name = c("National Capital Region Network","Antietam National Battlefield", "Catoctin Mountain Park","Chesapeake and Ohio Canal National Historical Park","Greenbelt Park","George Washington Memorial Parkway","Harpers Ferry National Historical Park","Manassas National Battlefield Park","Monocacy National Battlefield Park","National Capital Parks-East","Piscataway/Fort Washington","Prince William Forest Park","Rock Creek Park","Wolf Trap National Park for the Performing Arts"))
  
  bci.merge<-merge(bci.data, UnitNames.df, by="Unit_Code",all.x=TRUE)
  
  #create Unit.Year column
  bci.merge$Park.Year<-paste(bci.merge$Unit_Name, bci.merge$Unit_Code, bci.merge$Year, sep=".")
  
  bci.summary<-summaryFunction(dataIn=bci.merge, factor="Park.Year", response="BCI")
  
  covs<-read.table(text=as.character(bci.summary$Park.Year), sep=".")
  colnames(covs)<-c("Unit_Name","Unit_Code","Year")
  bci.summary.2<-cbind(covs, bci.summary)
  names(bci.summary.2)[names(bci.summary.2)=="mean"]<-"BCI"
  bci.summary.2$Year<-as.factor(as.character(bci.summary.2$Year))
  
  #############################################
  #Intro paragraph
  
  #get parkList
  #parkList<-sort(unique(as.character(data$Unit_Code)))
  parkList<-c("ANTI","CATO","CHOH","GWMP","GREE","HAFE","MANA","MONO","NACE","PISC_FOWA","PRWI","ROCR","WOTR")
  #parkList<-c("ANTI","CATO")
  
  #read in park shapefiles
  #NCRN.shp<-readOGR(dsn="/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/GIS/allparkboundariesshp.shp",verbose=FALSE)
  
  #rename PISC to PISC_FOWA ALPHACODE
  #levels(NCRN.shp$ALPHACODE)[levels(NCRN.shp$ALPHACODE)=="PISC"] <- "PISC_FOWA"
  
  #for loop to loop through parks
  for(i in 1:length(parkList)){
    
    header.filepath<-paste("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Reports/Markdown/10-year_report_park_summaries/Markdown",paste(parkList[i], "_header.png",sep=""),sep="/")
    
    #load header img
    cat("![](",header.filepath,")")
    
    #switch to 2-column format
    cat("\n\\btwocol\n")
    
    #subset data by park 
    park.data<-subset(data, Unit_Code==parkList[i])
    
    #get park.code
    park.code<-as.character(unique(park.data$Unit_Code))
    
    #get full unitName
    parkNameFull<-switch(park.code,
                         "ANTI"="Antietam National Battlefield",
                         "CATO"="Catoctin Mountain Park",
                         "CHOH"="Chesapeake & Ohio Canal National Historic Park",
                         "GREE"="Greenbelt Park",
                         "GWMP"="George Washington Memorial Parkway",
                         "HAFE" ="Harper's Ferry National Historic Park",
                         "MANA"="Mannassas National Battlefield Park",
                         "MONO"="Monocacy National Battlefield Park",
                         "NACE"="National Capital Parks-East",
                         "PISC_FOWA"="Piscataway Park/Fort Washington",
                         "PRWI"="Prince William Forest Park",
                         "ROCR"="Rock Creek Park",
                         "WOTR"="Wolf Trap Center for the Performing Arts")
    
    #make Unit-based zoom list
    parkZoom<-switch(park.code,
                     "ANTI"=14,
                     "CATO"=13,
                     "CHOH"= 9,
                     "GREE"= 15,
                     "GWMP"=12,
                     "HAFE"=13,
                     "MANA"=13,
                     "MONO"=13,
                     "NACE"=13,
                     "PISC_FOWA"=13,
                     "PRWI"=13,
                     "ROCR"=13,
                     "WOTR"=15)
    
    #get num.species
    num.species<-length(unique(park.data$AOU_Code))
    
    #get unique.Plot
    unique.plot<-sort(unique(park.data$Plot_Name))
    
    #get list of points with coordinates 
    park.points<-unique(park.data[,c("Plot_Name","Lat_WGS84","Long_WGS84","UTM_X_Coord","UTM_Y_Coord","UTM_Zone")])
    
    count.points<-length(row.names(park.points))
    
    xy.utm<-park.points[,c("UTM_X_Coord","UTM_Y_Coord")]
    park.sp <- SpatialPoints(xy.utm, proj4string=CRS("+proj=utm +zone=18 +datum=WGS84") ) 
    xy.latlong <- spTransform(park.sp, CRS("+proj=longlat +datum=WGS84"))
    
    park.points.2<-data.frame(park.points$Plot_Name, xy.latlong@coords)
    colnames(park.points.2)<-c("Plot_Name","Longitude","Latitude")
    
    #add species richness column
    park.points.merge<-merge(park.points.2, sp.richness, by="Plot_Name",all.x=TRUE)
    
    #####################################
    #BCI
    
    #subset by park
    bci.summary.sub<-subset(bci.summary.2, Unit_Code==park.code)
    bci.summary.sub$Year<-as.integer(as.character(bci.summary.sub$Year))
    
    #fit linear model and get slope (Beta coefficient)
    new.lm<-lm(BCI~Year,data=bci.summary.sub)
    linear.trend<-round(new.lm$coefficients[2],3)
    
    #add to data
    trend.df<-data.frame(Trend=linear.trend)
    row.names(trend.df)<-NULL
    
    #add categorical for slope (neg or pos)
    trend.df$TrendCat<-ifelse(trend.df$Trend > 0.05,"positive",
                              ifelse(trend.df$Trend < -0.05,"negative","stable"))
    
    #bci.trend (number)
    bci.trend<- trend.df$Trend
    
    #text (positive, stable, negative)
    bci.trend.txt<- trend.df$TrendCat
    ######################################
    
    #add text for header
    cat("\n\n")
    
    #add summary paragraph
    cat("We collected data at",count.points," monitoring locations (Dawson and Efford 2013; Fig. 1) which were each visited twice per year. During each 10-min visit, all individual birds were identified to species by a skilled observer within a 100-m radius circular plot and tallied. Using species detections, we computed a Bird Community Index (BCI) that links bird community assemblages with ecological condition based on vegetation structure and diversity (O'Connell 2003). We used the 'unmarked' (Fiske et al. 2011) package in R (R Core Team 2017) to estimate the probability of site occupancy (MacKenzie et al. 2002) and density (birds per ha; Chandler et al. 2011) for ", num.species, " unique species within ", parkNameFull, " between 2007-2017. Our models accounted for imperfect detection by including Observer, Ordinal Day, Time, and Wind covariates that are known to effect detection probability (Royle and Nichols 2003). We found a ", bci.trend.txt," trend of ", bci.trend," for BCI within ", parkNameFull,"(Fig. 2). Annual density estimates and linear model-based trends for the top-10 most abundant species are shown (Fig. 3).")
    
    # For more detailed information see: Ladin, Z. S., E. Tymkiw, S.Roberts, and W. G. Shriver. 2018. Forest Bird Monitoring in the National Capital Region Network: 2007 - 2017. Natural Resource Technical Report NPS//NCRN/NRDS-2018/***. National Park Service, Fort Collins, Colorado, USA."
    
    
    #######################################################################
    #Fig. 1. Park map with species richness
    #import shapefile
    #park.shp<-subset(NCRN.shp, NCRN.shp@data$ALPHACODE==parkList[i])
    #view
    #plot(park.shp,col="darkgreen")
    
    #look at current projection
    #proj4string(park.shp)
    
    #convert from UTM to lat/long
    #park.shp.proj <- spTransform(park.shp, CRS("+proj=longlat +datum=WGS84"))
    #proj4string(park.shp.proj)
    
    # convert to dataframe
    #park.map.df <- fortify(park.shp.proj)
    
    #set bounding box from coords (from shapefile)
    # min.x<-min(park.map.df$long, na.rm=TRUE)-0.5
    # max.x<-max(park.map.df$long, na.rm=TRUE)+0.5
    # min.y<-min(park.map.df$lat, na.rm=TRUE)-0.5
    # max.y<-max(park.map.df$lat, na.rm=TRUE)+0.5
    
    #set bounding box from coords (park.points.2)
    min.x<-min(park.points.2$Longitude, na.rm=TRUE)-0.5
    max.x<-max(park.points.2$Longitude, na.rm=TRUE)+0.5
    min.y<-min(park.points.2$Latitude, na.rm=TRUE)-0.5
    max.y<-max(park.points.2$Latitude, na.rm=TRUE)+0.5
    
    
    bbox<-c(left=min.x, bottom=min.y,right=max.x,top=max.y)
    
    #get polygon centroid (from shapefile)
    # mean.x<-mean(park.map.df$long)
    # mean.y=mean(park.map.df$lat)
    
    #get polygon centroid (from points)
    mean.x<-mean(park.points.2$Longitude)
    mean.y=mean(park.points.2$Latitude)
    
    
    map.center<-c(mean.x, mean.y)
    
    map.extent<-extent( min.x , max.x , min.y  , max.y )
    
    #get terrain map
    #park.map<-get_map(location=map.center, maptype="terrain",zoom=parkZoom)
    
    #get map with get_google to have finer control over map output
    park.map <- get_googlemap(center = map.center, zoom = parkZoom, maptype="terrain",
                              style = 'feature:administrative.country|element:labels|visibility:off')
    
    ifelse(park.code=="ANTI", {legend.coords<-c(0.005, 0.15)}, {legend.coords<-c(0.75,0.85)})
    
    park.map.out <- ggmap(park.map)+
      #geom_polygon(data=park.map.df, aes(x=long, y=lat, group=group),fill="darkgreen",color="black",alpha=0.6)+
      geom_point(data=park.points.merge, aes(x=Longitude,y=Latitude, size=Count_of_Species),color="red",alpha=0.9)+
      scale_size_continuous(range = c(0.5, 2))+
      #geom_jitter(width = 0.3, height = 0.3)+
      theme(panel.border=element_rect(fill="transparent",color="black"))+
      theme(panel.background=element_rect(fill='white',color="black"))+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      theme(axis.line=element_line(color="black"))+
      theme(axis.text.x = element_text(size=10, color="black"),
            axis.text.y = element_text(size=10, color="black"),
            axis.title.x = element_text(size=11, hjust=0.5, vjust=1.9),
            axis.title.y = element_text(angle = 90, vjust=1.2, size=11),
            title=element_text(size=8))+
      theme(legend.position = legend.coords,
            legend.justification = "left",
            legend.title=element_text(size=10),
            legend.background = element_rect(alpha("white", 0.7)))+
      theme(legend.key.size =  unit(0.05, "in"))+ 
      theme(legend.text = element_text(size=6))+ 
      #geom_text(data=park.points.merge, aes(x=Longitude,y=Latitude,label=Count_of_Species),color="white", hjust=0.5, vjust=0.4, size=2)+
      labs(x="Longitude",y="Latitude", size=paste("Species\nRichness"))+
      ggtitle("Figure 1. Map showing monitoring locations\nand species richness.")
    
    #print(park.map.out)
    
    subchunkify(park.map.out,4,4)
    
    ############################################################################################################
    #Figure 2 (BCI)
    
    # make plots with background colors to show BCI integrity ranges
    
    #create background data.frame
    background <- data.frame( lower = c( -Inf, 40.1, 52.1, 60.1), 
                              upper = c( 40.1, 52.1, 60.1, Inf), 
                              Integrity = c("Low Integrity", "Medium Integrity","High Integrity","Highest Integrity"),
                              order=c(1,2,3,4))
    #reorder Integrity factor
    background.order<-background[order(background$lower)]
    bg.levels <- background.order$Integrity
    background.order$order <- factor(background.order$Integrity, levels = bg.levels)
    #head(background.order)
    
    ##################################################################################################################
    
    bci.plot<-ggplot()+
      geom_rect(data = background.order, mapping= aes(ymin = lower, ymax = upper , xmin = -Inf , xmax = Inf, 
                                                      fill = order), alpha = 0.6)+
      geom_errorbar(data=bci.summary.sub, aes(x=Year, y=BCI, ymin=BCI-SE, ymax=BCI+SE, color=BCI), width=0, size=0.35, color="darkgreen")+     
      #geom_line(data=bci.summary.sub,aes(x=Year, y=BCI,group=1),color="darkgreen")+
      #geom_point(data=BCIdata,size=4)+
      #geom_smooth(data=bci.summary.sub, aes(x=Year, y=BCI), method="lm",se=FALSE, color="darkgreen")+
      geom_point(data=bci.summary.sub, aes(x=Year, y=BCI), color="darkgreen",shape=19,size=2)+
      theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
      theme(panel.background=element_rect(fill='white',color="black"))+
      theme(axis.text.x = element_text(angle = 55, hjust = 0.5,vjust=0.5, size=10, color="black"), 
            axis.text.y = element_text(size=10, color="black"),
            axis.title.x = element_text(size=11, hjust=0.5, vjust=0.2),
            axis.title.y = element_text(angle = 90, vjust=1.2, size=11),
            title=element_text(size=8))+
      theme(axis.line=element_line(color="black"))+
      theme(legend.position="none")+
      scale_fill_brewer(palette="Greys")+
      #scale_color_gradient2(midpoint=mid,low="red",mid="yellow", high="darkgreen")+
      # guides(fill = guide_legend(nrow=2, byrow=TRUE,reverse=TRUE, title="Ecological Integrity", title.position="top",title.hjust =0.5))+
      labs(x="Year", y="BCI score")+
      ylim(30,65)+
      scale_x_continuous(breaks=c(2007, 2009, 2011, 2013, 2015, 2017))+
      ggtitle("Figure 2. Trend in Bird Community Index")+
      annotate("text", x = 2016.5, y = 39.5, label = "Low Integrity")+
      annotate("text", x = 2016.5, y = 51.5, label = "Medium Integrity")+
      annotate("text", x = 2016.5, y = 59.5, label = "High Integrity")+
      annotate("text", x = 2016.5, y = 64.5, label = "Highest Integrity")
    
    bci.plot
    #print(bci.plot)
    
    subchunkify(bci.plot, 2.5,4)
    ############################################################################################################
    
    
    #Figure 3. Top ten species trends
    
    #create d.f. of common.names
    common.names<-unique(data[,c("AOU_Code","Common_Name")])
    
    
    #######################################################################################################################
    #Create trend figures for each Park (and NCRN)
    
    #subset data by Park
    park.data<-subset(species.merge.2, Park==parkList[i])
    
    #get 11-year mean to sort species by
    park.summary<-summaryFunction(dataIn=park.data, response="mean", factor="AOU_Code")
    
    #order park.summary (constrain to top-10)
    park.summary.order<-head(park.summary[order(park.summary$mean,decreasing=TRUE),],10)
    
    #get top-10 list
    top.list<-as.character(park.summary.order$AOU_Code)
    
    #subset park.data jby top.list
    park.data<-subset(park.data, AOU_Code %in% top.list)
    
    parkName<-unique(as.character(park.data$Park))
    
    park.data$TrendCat<-factor(park.data$TrendCat, levels=c("Positive","Stable","Negative"))
    
    #pick colors for figures
    trend.length<-length(unique(park.data$TrendCat))
    ifelse(trend.length<3,
           {
             mycolors<-c("royalblue","red")
           },
           mycolors<-c("royalblue",I("gray50"),"red"))
    
    
    #plot bci mean annual BCI with linear trend
    abun.trend.plot<-ggplot(data=park.data, aes(x=Year, y=mean))+
      #geom_line(aes(color=TrendCat),size=0.4,alpha=0.2)+
      geom_path(aes(color=TrendCat, group=TrendCat),size=0.3, alpha=0.7)+
      #geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE, color=TrendCat), width=0)+
      #geom_smooth(method="lm",aes(color=TrendCat), se = FALSE,size=0.7)+
      #scale_color_viridis(option="viridis", begin=0.2, end=1, direction=-1)+
      #scale_color_gradientn(colors=colfunc(10))+
      #scale_color_gradient2(low="red",mid="black", high="royalblue", midpoint=0, limits = c(-0.2, 0.2))+
      scale_color_manual(name="Trend",values=mycolors,drop=FALSE)+
      labs(x="Year", y="Estimated density (birds per ha)")+
      theme(panel.background=element_rect(fill="transparent"), panel.border=element_rect(fill=NA,color="black"),
            axis.text.x = element_text(angle=55,size=10, color="black", vjust=0.5),
            axis.text.y = element_text(size=10, color="black"),
            axis.title.x = element_text(size=11, hjust=0.5, vjust=2.2),
            axis.title.y = element_text(angle = 90, vjust=1.2, size=11),
            title=element_text(size=8))+
      scale_x_continuous(breaks=c(2007, 2012, 2017))+
      #scale_y_continuous(breaks=c(0, 15, 30))+
      theme(legend.position="top")+
      ggtitle(paste("Figure 3. Trends in Annual Population Density",sep=""))
    #abun.trend.plot
    
    #create new facet labels uisng 
    new.lab <- as_labeller(park.data$SpeciesLabel, label_parsed)
    
    #create list of trends at each park
    trend.df<-unique(park.data[,c("SpeciesLabel","Trend","TrendCat")])
    colnames(trend.df)<-c("SpeciesLabel","Trend","TrendCat")
    length(row.names(trend.df))
    levels(trend.df$TrendCat)
    # trend.df$Trend<-as.character(trend.df$Trend)
    # trend.df$TrendCat<-as.character(trend.df$TrendCat)
    
    # v.just<-switch(parkName,
    #                "NCRN"=-2.7,-9)
    
    
    #abun.trend.plot
    abun.trend.plot.out<-NULL
    abun.trend.plot.out<-abun.trend.plot+
      geom_text(
        data= trend.df,
        size=2.5,
        mapping = aes(x = -Inf, y = -Inf, label = Trend, color=TrendCat),
        hjust   = -0.2,
        vjust   = -1.2
      )+
      facet_wrap(~SpeciesLabel, labeller = new.lab, ncol=2)+
      theme(strip.background =element_rect(fill=I("gray95")))
    
    
    #print(abun.trend.plot.out)
    #use subchunkify() function (object, height, width)
    subchunkify(abun.trend.plot.out, 5, 4)
    ############################################################################################################
    
    
    
    #add text for header
    #cat(paste("Information","\n\n",sep=""))
    
    #add summary paragraph
    cat("\n\n")
    
    cat("References")
    
    cat("\n\n")
    
    cat("Chandler, R. B., J. A. Royle, and D. I. King. 2011. Inference about density and temporary emigration in
        unmarked populations. Ecology 92:1429–1435.") 
    
    cat("\n\n")
    
    cat("Dawson, D. K., and M. G. Efford. 2013. Protocol for monitoring forest-nesting birds in national park service
        parks. National Park Service.")
    
    cat("\n\n")
    
    cat("Fiske, I., R. Chandler, and others. 2011. Unmarked: An r package for fitting hierarchical models of wildlife
        occurrence and abundance. Journal of Statistical Software 43:1–23.") 
    
    cat("\n\n")
    
    cat("MacKenzie, D. I., J. D. Nichols, G. B. Lachman, S. Droege, J. Andrew Royle, and C. A. Langtimm. 2002.
        Estimating site occupancy rates when detection probabilities are less than one. Ecology 83:2248–2255.") 
    
    cat("\n\n")
    
    cat("O'Connell, T., R. Brooks, M. Lanzone, and J. Bishop. 2003. A bird community index for the mid-atlantic
        piedmont and coastal plain, final report to the usgs-patuxent wildlife research center. penn state cooperative
        wetlands center, university park, p 44. Report.")
    
    cat("\n\n")
    
    cat("R Core Team. 2017. R: A language and environment for statistical computing. R Foundation for Statistical
        Computing, Vienna, Austria.")
    
    cat("\n\n")
    
    cat("Royle, J. A., and J. D. Nichols. 2003. Estimating abundance from repeated presence-absence data or point
        counts. Ecology 84:777–790.") 
    
    
    cat("\n\\etwocol\n")
    
    cat("\n\\clearpage\n")
    











#####################################################################################################################