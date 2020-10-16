#make Generalized Multinomial Model unmarkedFrame
makeGMMumf<-function(DataIn,SpeciesIn){
  require(unmarked)

  #read in data
  data.new<-DataIn

  #get parkName for figures
  parkName<-sort(unique(as.character(data.new$Admin_Unit_Code)))
  
  #create habitat type
  habitatType=unique(data.new$Survey_Type)
  
  #create columns
  data.new$Year.Unit.Plot.Visit<-paste(data.new$Year, data.new$Unit_Code, data.new$Plot_Name, data.new$ Visit,sep="_")
  
  #subset data to remove visits 0 and >1 from data
  ifelse(habitatType=="Forest",
         {data.new<-subset(data.new, c(Visit != 0 & Visit !=3))
         },
         {data.new<-subset(data.new, c(Visit != 0))
         })
  
  #create Visit.Interval column
  data.new$Visit.Interval<-paste(data.new$Visit, data.new$Interval, sep="_")
  
  #create Year.Unit.Plot column
  data.new$Year.Unit.Plot<-paste(data.new$Year, data.new$Unit_Code, data.new$Plot_Name,sep="_")
  
  #get detection data and covs with melt and cast (this is where you can select which Minutes (intervals) to use, here Min1-Min5)
  data.count.melt <- melt(data.new, id=c( "Unit_Name", "Unit_Code","Plot_Name","Year", "Month","Day","Ord.Day","Date","Visit","Lat_WGS84","Long_WGS84", "Start_Time.dec","Interval", "AOU_Code","Sex_ID","Distance_id","Observer","Visit.Interval","Year.Unit.Plot"),
                          measure=c("CountOfAOU_Code"), na.rm=FALSE)

  data.count.melt$Visit.Interval<-as.factor(data.count.melt$Visit.Interval)
  
  habitatType<-as.character(unique(data.new$Survey_Type))
  
  
  #cast count data (one visit at first) add try()
  ifelse(habitatType=="Forest",
         {count.data <- cast(data.count.melt, Year.Unit.Plot ~ Visit.Interval, add.missing=TRUE, fun.aggregate=sum, subset=c(AOU_Code==SpeciesIn))
         colnames(count.data)<-c("Year.Unit.Plot","y.1","y.2","y.3", "y.4", "y.5", "y.6", "y.7", "y.8")},
         {count.data <- cast(data.count.melt, Year.Unit.Plot ~ Visit.Interval, add.missing=TRUE, fun.aggregate=sum, subset=c(AOU_Code==SpeciesIn))
         colnames(count.data)<-c("Year.Unit.Plot","y.1","y.2","y.3", "y.4", "y.5", "y.6", "y.7", "y.8","y.9","y.10","y.11","y.12")}
         )
  
  #get some site covariates from new.data
  head(data.new)
  data.new.covs<-data.new[,c("Year.Unit.Plot","Year","Unit_Code","Plot_Name","Visit", "Long_WGS84","Lat_WGS84","Date","Month","Ord.Day","Start_Time.dec","Temperature","Wind_Speed","Observer")]

  data.new.covs.unique<-unique(data.new.covs)
  
  all.sampled<-data.frame(Year.Unit.Plot=unique(data.new.covs[,c("Year.Unit.Plot")]))
  
  #merge with count data
  count.data.merge<-merge(count.data, all.sampled, by="Year.Unit.Plot",all.y=TRUE)
  count.data.merge[is.na(count.data.merge)]<-0
  
  
  ifelse(habitatType=="Forest",
         {count.data.out <- data.frame(count.data.merge,"1","2", "1", "2","3","4","1","2","3","4")
         colnames(count.data.out)<-c("Year.Unit.Plot","y.1","y.2","y.3", "y.4", "y.5", "y.6", "y.7", "y.8", "Visit_1", "Visit_2","v1.i1", "v1.i2","v1.i3","v1.i4","v2.i1","v2.i2","v2.i3","v2.i4")},
         {count.data.out <- data.frame(count.data.merge,"1","2","3", "1", "2","3","4","1","2","3","4","1","2","3","4")
         colnames(count.data.out)<-c("Year.Unit.Plot","y.1","y.2","y.3", "y.4", "y.5", "y.6", "y.7", "y.8","y.9","y.10","y.11","y.12", "Visit_1", "Visit_2","Visit_3","v1.i1", "v1.i2","v1.i3","v1.i4","v2.i1","v2.i2","v2.i3","v2.i4","v3.i1","v3.i2","v3.i3","v3.i4")}
  )
  
  #format count data
  covs.melt<- melt(data.new.covs.unique,id=c("Year.Unit.Plot","Year","Unit_Code","Plot_Name","Visit"), measure=c("Ord.Day","Month","Start_Time.dec","Temperature","Wind_Speed","Observer"))
  head(covs.melt)
  
  #cast data
  all.site.covs<-cast(covs.melt, Year.Unit.Plot ~ Visit, add.missing=TRUE,fun.aggregate="length")
  head(all.site.covs)
  

  covs.day<-cast(covs.melt, Year.Unit.Plot ~ Visit, add.missing=TRUE, fun.aggregate="max", subset=(variable=="Ord.Day"))
  covs.month<-cast(covs.melt, Year.Unit.Plot ~ Visit, add.missing=TRUE, fun.aggregate="max", subset=(variable=="Month"))
  covs.time<-cast(covs.melt, Year.Unit.Plot ~ Visit, add.missing=TRUE,fun.aggregate="max", subset=(variable=="Start_Time.dec"))
  covs.temp<-cast(covs.melt, Year.Unit.Plot ~ Visit, add.missing=TRUE,fun.aggregate="max", subset=(variable=="Temperature"))
  covs.wind<-cast(covs.melt, Year.Unit.Plot ~ Visit, add.missing=TRUE,fun.aggregate="max", subset=(variable=="Wind_Speed"))
  covs.observer<-cast(covs.melt, Year.Unit.Plot ~ Visit, add.missing=TRUE,fun.aggregate="max", subset=(variable=="Observer"))

  #create dataframe of covariates
  covs.combine<-cbind(covs.month, covs.day[,-1], covs.time[,-1], covs.temp[,-1], covs.wind[,-1], covs.observer[,-1])
  colnames(covs.combine)<-c("Year.Unit.Plot","Month_1","Month_2","Ord.Day_1","Ord.Day_2","Time_1","Time_2","Temp_1","Temp_2","Wind_1","Wind_2","Obs_1","Obs_2")
  
  ifelse(habitatType=="Forest",
         {covs.combine<-cbind(covs.month, covs.day[,-1], covs.time[,-1], covs.temp[,-1], covs.wind[,-1], covs.observer[,-1])
         colnames(covs.combine)<-c("Year.Unit.Plot","Month_1","Month_2","Ord.Day_1","Ord.Day_2","Time_1","Time_2","Temp_1","Temp_2","Wind_1","Wind_2","Obs_1","Obs_2")},
         {covs.combine<-cbind(covs.month, covs.day[,-1], covs.time[,-1], covs.temp[,-1], covs.wind[,-1], covs.observer[,-1])
         colnames(covs.combine)<-c("Year.Unit.Plot","Month_1","Month_2","Month_3","Ord.Day_1","Ord.Day_2","Ord.Day_3","Time_1","Time_2","Time_3","Temp_1","Temp_2","Temp_3","Wind_1","Wind_2","Wind_3","Obs_1","Obs_2","Obs_3")}
  )
  

  #combine count data with detection covariates
  bird.gmult.data <- merge(count.data.merge, covs.combine, by="Year.Unit.Plot", all.y=TRUE)
  head(bird.gmult.data)
  
  #add Visit
  ifelse(habitatType=="Forest",
         {bird.gmult.data$Visit_1<-1
         bird.gmult.data$Visit_2<-2},
         {bird.gmult.data$Visit_1<-1
         bird.gmult.data$Visit_2<-2
         bird.gmult.data$Visit_3<-3})
  
  

  #get nYears, nSites
  site.covs<-read.table(text=as.character(bird.gmult.data$Year.Unit.Plot),colClasses="character",sep="_")
  colnames(site.covs)<-c("Year","Unit_Code","Plot_Name")
  count.data.out<-cbind(site.covs,bird.gmult.data)
  count.data.out$Plot_Name<-as.factor(as.character(count.data.out$Plot_Name))
  count.data.out$Year<-as.factor(as.character(count.data.out$Year))

  #number of survey locations
  nSites<-length(unique(count.data.out$Plot_Name))

  #number of visits
  ifelse(habitatType=="Forest",
         nVisits<-as.integer(max(unique(as.integer(count.data.out$Visit_2)),na.rm=TRUE)),
         nVisits<-as.integer(max(unique(as.integer(count.data.out$Visit_3)),na.rm=TRUE))
         )

  #number of intervals
  nIntervals<-4


  bird.data<-count.data.out
  head(bird.data)

  #########################################################################################
  #fill in NAs in count data, where needed
  ifelse(habitatType=="Forest",
         {
           bird.data$y.1<-ifelse(is.na(bird.data$Month_1),"NA",bird.data$y.1)
           bird.data$y.2<-ifelse(is.na(bird.data$Month_1),"NA",bird.data$y.2)
           bird.data$y.3<-ifelse(is.na(bird.data$Month_1),"NA",bird.data$y.3)
           bird.data$y.4<-ifelse(is.na(bird.data$Month_1),"NA",bird.data$y.4)
           
           bird.data$y.5<-ifelse(is.na(bird.data$Month_2),"NA",bird.data$y.5)
           bird.data$y.6<-ifelse(is.na(bird.data$Month_2),"NA",bird.data$y.6)
           bird.data$y.7<-ifelse(is.na(bird.data$Month_2),"NA",bird.data$y.7)
           bird.data$y.8<-ifelse(is.na(bird.data$Month_2),"NA",bird.data$y.8)
           
           bird.data$Visit_1<-ifelse(is.na(bird.data$Month_1),"NA",bird.data$Visit_1)
           bird.data$Visit_2<-ifelse(is.na(bird.data$Month_2),"NA",bird.data$Visit_2)
         },
         {
           bird.data$y.1<-ifelse(is.na(bird.data$Month_1),"NA",bird.data$y.1)
           bird.data$y.2<-ifelse(is.na(bird.data$Month_1),"NA",bird.data$y.2)
           bird.data$y.3<-ifelse(is.na(bird.data$Month_1),"NA",bird.data$y.3)
           bird.data$y.4<-ifelse(is.na(bird.data$Month_1),"NA",bird.data$y.4)
           
           bird.data$y.5<-ifelse(is.na(bird.data$Month_2),"NA",bird.data$y.5)
           bird.data$y.6<-ifelse(is.na(bird.data$Month_2),"NA",bird.data$y.6)
           bird.data$y.7<-ifelse(is.na(bird.data$Month_2),"NA",bird.data$y.7)
           bird.data$y.8<-ifelse(is.na(bird.data$Month_2),"NA",bird.data$y.8)
           
           bird.data$y.9<-ifelse(is.na(bird.data$Month_3),"NA",bird.data$y.9)
           bird.data$y.10<-ifelse(is.na(bird.data$Month_3),"NA",bird.data$y.10)
           bird.data$y.11<-ifelse(is.na(bird.data$Month_3),"NA",bird.data$y.11)
           bird.data$y.12<-ifelse(is.na(bird.data$Month_3),"NA",bird.data$y.12)
           
           bird.data$Visit_1<-ifelse(is.na(bird.data$Month_1),"NA",bird.data$Visit_1)
           bird.data$Visit_2<-ifelse(is.na(bird.data$Month_2),"NA",bird.data$Visit_2)
           bird.data$Visit_3<-ifelse(is.na(bird.data$Month_3),"NA",bird.data$Visit_3)
           
         }
         )
  
  ##########################################################################################
  #convert variables to correct class
  ifelse(habitatType=="Forest",
    {bird.data$y.1<-as.integer(bird.data$y.1)
    bird.data$y.2<-as.integer(bird.data$y.2)
    bird.data$y.3<-as.integer(bird.data$y.3)
    bird.data$y.4<-as.integer(bird.data$y.4)
    bird.data$y.5<-as.integer(bird.data$y.5)
    bird.data$y.6<-as.integer(bird.data$y.6)
    bird.data$y.7<-as.integer(bird.data$y.7)
    bird.data$y.8<-as.integer(bird.data$y.8)
    },
    {bird.data$y.1<-as.integer(bird.data$y.1)
    bird.data$y.2<-as.integer(bird.data$y.2)
    bird.data$y.3<-as.integer(bird.data$y.3)
    bird.data$y.4<-as.integer(bird.data$y.4)
    bird.data$y.5<-as.integer(bird.data$y.5)
    bird.data$y.6<-as.integer(bird.data$y.6)
    bird.data$y.7<-as.integer(bird.data$y.7)
    bird.data$y.8<-as.integer(bird.data$y.8)
    bird.data$y.9<-as.integer(bird.data$y.9)
    bird.data$y.10<-as.integer(bird.data$y.10)
    bird.data$y.11<-as.integer(bird.data$y.11)
    bird.data$y.12<-as.integer(bird.data$y.12)
    }
  )
  
  
  ifelse(habitatType=="Forest",
         {
         bird.data$Month_1<-as.factor(as.character(bird.data$Month_1))
         bird.data$Month_2<-as.factor(as.character(bird.data$Month_2))
         bird.data$Ord.Day_1<-as.integer(as.character(bird.data$Ord.Day_1))
         bird.data$Ord.Day_2<-as.integer(as.character(bird.data$Ord.Day_2))
         bird.data$Time_1<-as.numeric(as.character(bird.data$Time_1))
         bird.data$Time_2<-as.numeric(as.character(bird.data$Time_2))
         bird.data$Temp_1<-as.numeric(as.character(bird.data$Temp_1))
         bird.data$Temp_2<-as.numeric(as.character(bird.data$Temp_2))
         bird.data$Wind_1<-as.integer(as.character(bird.data$Wind_1))
         bird.data$Wind_2<-as.integer(as.character(bird.data$Wind_2))
         bird.data$Obs_1<-as.factor(as.character(bird.data$Obs_1))
         bird.data$Obs_2<-as.factor(as.character(bird.data$Obs_2))
         bird.data$Visit_1<-as.factor(as.character(bird.data$Visit_1))
         bird.data$Visit_2<-as.factor(as.character(bird.data$Visit_2))
         },
         {
           bird.data$Month_1<-as.factor(as.character(bird.data$Month_1))
           bird.data$Month_2<-as.factor(as.character(bird.data$Month_2))
           bird.data$Month_3<-as.factor(as.character(bird.data$Month_3))
           bird.data$Ord.Day_1<-as.integer(as.character(bird.data$Ord.Day_1))
           bird.data$Ord.Day_2<-as.integer(as.character(bird.data$Ord.Day_2))
           bird.data$Ord.Day_3<-as.integer(as.character(bird.data$Ord.Day_3))
           bird.data$Time_1<-as.numeric(as.character(bird.data$Time_1))
           bird.data$Time_2<-as.numeric(as.character(bird.data$Time_2))
           bird.data$Time_3<-as.numeric(as.character(bird.data$Time_3))
           bird.data$Temp_1<-as.numeric(as.character(bird.data$Temp_1))
           bird.data$Temp_2<-as.numeric(as.character(bird.data$Temp_2))
           bird.data$Temp_3<-as.numeric(as.character(bird.data$Temp_3))
           bird.data$Wind_1<-as.integer(as.character(bird.data$Wind_1))
           bird.data$Wind_2<-as.integer(as.character(bird.data$Wind_2))
           bird.data$Wind_3<-as.integer(as.character(bird.data$Wind_3))
           bird.data$Obs_1<-as.factor(as.character(bird.data$Obs_1))
           bird.data$Obs_2<-as.factor(as.character(bird.data$Obs_2))
           bird.data$Obs_3<-as.factor(as.character(bird.data$Obs_3))
           bird.data$Visit_1<-as.factor(as.character(bird.data$Visit_1))
           bird.data$Visit_2<-as.factor(as.character(bird.data$Visit_2))
           bird.data$Visit_3<-as.factor(as.character(bird.data$Visit_3))
         }
  )
  
  ##########################################################################################
  
  #define response variable
  ifelse(habitatType=="Forest",
         y<-data.frame(bird.data$y.1,bird.data$y.2,bird.data$y.3,bird.data$y.4,
                       bird.data$y.5,bird.data$y.6,bird.data$y.7,bird.data$y.8),
         y<-data.frame(bird.data$y.1,bird.data$y.2,bird.data$y.3,bird.data$y.4,
                       bird.data$y.5,bird.data$y.6,bird.data$y.7,bird.data$y.8,
                       bird.data$y.9,bird.data$y.10,bird.data$y.11,bird.data$y.12)
         )
  
  #define site covariates - FOR MULTIPLE YEARS
  site.covs<- data.frame(Year=as.factor(bird.data$Year), Unit_Code=as.factor(bird.data$Unit_Code), 
                         Plot_Name=as.factor(bird.data$Plot_Name))
  
  #create yearlySiteCovs matrices (Abun and Availability)
  ifelse(habitatType=="Forest",
         {
           Ord.Day=data.frame(Ord.Day1=bird.data$Ord.Day_1, Ord.Day2=bird.data$Ord.Day_2)
           Visit=data.frame(Visit1=bird.data$Visit_1, Visit2=bird.data$Visit_2)
           Time=data.frame(Time1=bird.data$Time_1, Time2=bird.data$Time_2)
         },
         {
           Ord.Day=data.frame(Ord.Day1=bird.data$Ord.Day_1, Ord.Day2=bird.data$Ord.Day_2, Ord.Day3=bird.data$Ord.Day_3)
           Visit=data.frame(Visit1=bird.data$Visit_1, Visit2=bird.data$Visit_2, Visit3=bird.data$Visit_3)
           Time=data.frame(Time1=bird.data$Time_1, Time2=bird.data$Time_2, Time3=bird.data$Time_3)
         }
         )

  #define abundance covariates (yearly site covs)
  yearly.site.covs<-list(Visit=Visit, Ord.Day=Ord.Day, Time=Time)
  
  #create obsCovs matrices (Detection Prob)
  # Ord.Day.obs=data.frame(bird.data$Ord.Day_1, bird.data$Ord.Day_1,bird.data$Ord.Day_1,bird.data$Ord.Day_1,
  #                    bird.data$Ord.Day_2, bird.data$Ord.Day_2,bird.data$Ord.Day_2,bird.data$Ord.Day_2)
  # Visit.obs=data.frame(bird.data$Visit_1, bird.data$Visit_1,bird.data$Visit_1,bird.data$Visit_1,
  #                  bird.data$Visit_2,bird.data$Visit_2,bird.data$Visit_2,bird.data$Visit_2)
  # Time.obs=data.frame(bird.data$Time_1,bird.data$Time_1,bird.data$Time_1,bird.data$Time_1, 
  #                 bird.data$Time_2,bird.data$Time_2,bird.data$Time_2,bird.data$Time_2)
  # Temp.obs=data.frame(bird.data$Temp_1, bird.data$Temp_1,bird.data$Temp_1,bird.data$Temp_1,
  #                 bird.data$Temp_2,bird.data$Temp_2,bird.data$Temp_2,bird.data$Temp_2)
  # Wind.obs=data.frame(bird.data$Wind_1, bird.data$Wind_1,bird.data$Wind_1,bird.data$Wind_1,
  #                 bird.data$Wind_2,bird.data$Wind_2,bird.data$Wind_2,bird.data$Wind_2 ) 
  # Observer.obs=data.frame(bird.data$Obs_1, bird.data$Obs_1,bird.data$Obs_1,bird.data$Obs_1,
  #                     bird.data$Obs_2,bird.data$Obs_2,bird.data$Obs_2,bird.data$Obs_2)
  
  #create obsCovs matrices (Detection Prob)
  ifelse(habitatType=="Forest",
         {
           Ord.Day.obs=data.frame(bird.data$Ord.Day_1, bird.data$Ord.Day_1,bird.data$Ord.Day_1,bird.data$Ord.Day_1,
                                  bird.data$Ord.Day_2, bird.data$Ord.Day_2,bird.data$Ord.Day_2,bird.data$Ord.Day_2)
           Visit.obs=data.frame(bird.data$Visit_1, bird.data$Visit_1,bird.data$Visit_1,bird.data$Visit_1,
                                bird.data$Visit_2,bird.data$Visit_2,bird.data$Visit_2,bird.data$Visit_2)
           Time.obs=data.frame(bird.data$Time_1,bird.data$Time_1,bird.data$Time_1,bird.data$Time_1, 
                               bird.data$Time_2,bird.data$Time_2,bird.data$Time_2,bird.data$Time_2)
           Temp.obs=data.frame(bird.data$Temp_1, bird.data$Temp_1,bird.data$Temp_1,bird.data$Temp_1,
                               bird.data$Temp_2,bird.data$Temp_2,bird.data$Temp_2,bird.data$Temp_2)
           Wind.obs=data.frame(bird.data$Wind_1, bird.data$Wind_1,bird.data$Wind_1,bird.data$Wind_1,
                               bird.data$Wind_2,bird.data$Wind_2,bird.data$Wind_2,bird.data$Wind_2 ) 
           Observer.obs=data.frame(bird.data$Obs_1, bird.data$Obs_1,bird.data$Obs_1,bird.data$Obs_1,
                                   bird.data$Obs_2,bird.data$Obs_2,bird.data$Obs_2,bird.data$Obs_2)
         },
         {
           Ord.Day.obs=data.frame(bird.data$Ord.Day_1, bird.data$Ord.Day_1,bird.data$Ord.Day_1,bird.data$Ord.Day_1,
                                  bird.data$Ord.Day_2, bird.data$Ord.Day_2,bird.data$Ord.Day_2,bird.data$Ord.Day_2,
                                  bird.data$Ord.Day_3, bird.data$Ord.Day_3,bird.data$Ord.Day_3,bird.data$Ord.Day_3)
           Visit.obs=data.frame(bird.data$Visit_1, bird.data$Visit_1,bird.data$Visit_1,bird.data$Visit_1,
                                bird.data$Visit_2,bird.data$Visit_2,bird.data$Visit_2,bird.data$Visit_2,
                                bird.data$Visit_3,bird.data$Visit_3,bird.data$Visit_3,bird.data$Visit_3)
           Time.obs=data.frame(bird.data$Time_1,bird.data$Time_1,bird.data$Time_1,bird.data$Time_1, 
                               bird.data$Time_2,bird.data$Time_2,bird.data$Time_2,bird.data$Time_2,
                               bird.data$Time_3,bird.data$Time_3,bird.data$Time_3,bird.data$Time_3)
           Temp.obs=data.frame(bird.data$Temp_1, bird.data$Temp_1,bird.data$Temp_1,bird.data$Temp_1,
                               bird.data$Temp_2,bird.data$Temp_2,bird.data$Temp_2,bird.data$Temp_2,
                               bird.data$Temp_3,bird.data$Temp_3,bird.data$Temp_3,bird.data$Temp_3)
           Wind.obs=data.frame(bird.data$Wind_1, bird.data$Wind_1,bird.data$Wind_1,bird.data$Wind_1,
                               bird.data$Wind_2,bird.data$Wind_2,bird.data$Wind_2,bird.data$Wind_2,
                               bird.data$Wind_3,bird.data$Wind_3,bird.data$Wind_3,bird.data$Wind_3) 
           Observer.obs=data.frame(bird.data$Obs_1, bird.data$Obs_1,bird.data$Obs_1,bird.data$Obs_1,
                                   bird.data$Obs_2,bird.data$Obs_2,bird.data$Obs_2,bird.data$Obs_2,
                                   bird.data$Obs_3,bird.data$Obs_3,bird.data$Obs_3,bird.data$Obs_3)
         }
         
  )
  
  
  
  obs.covs<-list(Visit.obs=Visit.obs, Ord.Day.obs=Ord.Day.obs, Time.obs=Time.obs, Temp=Temp.obs, Wind=Wind.obs, Observer=Observer.obs)
  
  o2y <- diag(ncol(y))
  o2y[upper.tri(o2y)] <- 1
  o2y
  
  
  #set up unmarked FRAME GMM
         umf<-tryCatch({unmarkedFrameGMM(y=y, 
                               siteCovs = site.covs,
                               yearlySiteCovs= yearly.site.covs,
                               obsCovs= obs.covs,
                               numPrimary=2, 
                               obsToY=o2y,
                               type="removal")}
                       ,error=function(cond2){
                         cond2=unmarkedFrameGMM(y=y, 
                                                siteCovs = site.covs,
                                                yearlySiteCovs= yearly.site.covs,
                                                obsCovs= obs.covs,
                                                numPrimary=3, 
                                                obsToY=o2y,
                                                type="removal")
                         cond2
                       })

         
  return(umf)
}

#End