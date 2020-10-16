#make Generalized Multinomial Model unmarkedFrame
makePOccu_umf<-function(DataIn,SpeciesIn){
  require(unmarked)
  require(plyr)
  require(reshape2)

  #read in data
  data.new<-DataIn

  #get habitatType
  habitatType<-unique(as.character(data.new$Survey_Type))
  
  #get parkName for figures
  parkName<-sort(unique(as.character(data.new$Unit_Code)))
  
  #create columns
  data.new$Year.Unit.Plot.Visit<-paste(data.new$Year, data.new$Unit_Code, data.new$Plot_Name, data.new$ Visit,sep="+")
  
  #subset data to remove visits 0 and >1 from data
  ifelse(habitatType=="Forest",
  {data.new<-subset(data.new, c(Visit != 0 & Visit !=3 & Visit !=4))
    },
  {data.new<-subset(data.new, c(Visit != 0))
    })
  
  #create Visit.Interval column
  data.new$Visit.Interval<-paste(data.new$Visit, data.new$Interval, sep="_")
  
  #create Year.Unit.Plot column
  data.new$Year.Unit.Plot<-paste(data.new$Year, data.new$Unit_Code, data.new$Plot_Name,sep="+")
  
  #get detection data and covs with melt and cast (this is where you can select which Minutes (intervals) to use, here Min1-Min5)
  data.count.melt <- melt(data.new, id=c( "Unit_Name", "Unit_Code","Plot_Name","Year", "Month","Day","Ord.Day","Date","Visit","Lat_WGS84","Long_WGS84", "Start_Time.dec","Interval", "AOU_Code","Sex_ID","Distance_id","Observer","Visit.Interval","Year.Unit.Plot"),
                          measure=c("CountOfAOU_Code"), na.rm=FALSE)

  data.count.melt$Visit<-as.factor(data.count.melt$Visit)
  data.count.melt$Visit.Interval<-as.factor(data.count.melt$Visit.Interval)
  data.count.melt$Year.Unit.Plot<-as.factor(data.count.melt$Year.Unit.Plot)
  
  #cast count data (one visit at first) add try()
  ifelse(habitatType=="Forest",
         {count.data <- dcast(data.count.melt, Year.Unit.Plot ~ Visit, drop=FALSE, fun.aggregate=sum, subset=.(AOU_Code==SpeciesIn))
         colnames(count.data)<-c("Year.Unit.Plot","y.1","y.2")
         },
         {count.data <- dcast(data.count.melt, Year.Unit.Plot ~ Visit, drop=FALSE, fun.aggregate=sum, subset=.(AOU_Code==SpeciesIn))
         colnames(count.data)<-c("Year.Unit.Plot","y.1","y.2","y.3")}
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
         {count.data.out <- data.frame(count.data.merge,"1","2")
         colnames(count.data.out)<-c("Year.Unit.Plot","y.1","y.2", "Visit_1", "Visit_2")},
         {count.data.out <- data.frame(count.data.merge,"1","2","3")
         colnames(count.data.out)<-c("Year.Unit.Plot","y.1","y.2","y.3", "Visit_1", "Visit_2","Visit_3")}
  )
  
  #format count data
  covs.melt<- melt(data.new.covs.unique,id=c("Year.Unit.Plot","Year","Unit_Code","Plot_Name","Visit"), measure=c("Ord.Day","Month","Start_Time.dec","Temperature","Wind_Speed","Observer"))
  head(covs.melt)
  
  #cast data
  all.site.covs<-dcast(covs.melt, Year.Unit.Plot ~ Visit, drop=FALSE,fun.aggregate=length)
  head(all.site.covs)
  

  covs.day<-dcast(covs.melt, Year.Unit.Plot ~ Visit, drop=FALSE, fun.aggregate=max, subset=.(variable=="Ord.Day"))
  covs.month<-dcast(covs.melt, Year.Unit.Plot ~ Visit, drop=FALSE, fun.aggregate=max, subset=.(variable=="Month"))
  covs.time<-dcast(covs.melt, Year.Unit.Plot ~ Visit, drop=FALSE,fun.aggregate=max, subset=.(variable=="Start_Time.dec"))
  covs.temp<-dcast(covs.melt, Year.Unit.Plot ~ Visit, drop=FALSE,fun.aggregate=max, subset=.(variable=="Temperature"))
  covs.wind<-dcast(covs.melt, Year.Unit.Plot ~ Visit, drop=FALSE,fun.aggregate=max, subset=.(variable=="Wind_Speed"))
  covs.observer<-dcast(covs.melt, Year.Unit.Plot ~ Visit, drop=FALSE,fun.aggregate=max, subset=.(variable=="Observer"))

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
  site.covs<-read.table(text=as.character(bird.gmult.data$Year.Unit.Plot),colClasses="character",sep="+")
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
  #nIntervals<-2


  bird.data<-count.data.out
  head(bird.data)

  #########################################################################################
  #fill in NAs in count data, where needed
  ifelse(habitatType=="Forest",
         {
           bird.data$y.1<-ifelse(is.na(bird.data$Month_1),"NA",bird.data$y.1)
           bird.data$y.2<-ifelse(is.na(bird.data$Month_2),"NA",bird.data$y.2)

           bird.data$Visit_1<-ifelse(is.na(bird.data$Month_1),"NA",bird.data$Visit_1)
           bird.data$Visit_2<-ifelse(is.na(bird.data$Month_2),"NA",bird.data$Visit_2)
         },
         {
           bird.data$y.1<-ifelse(is.na(bird.data$Month_1),"NA",bird.data$y.1)
           bird.data$y.2<-ifelse(is.na(bird.data$Month_2),"NA",bird.data$y.2)
           bird.data$y.3<-ifelse(is.na(bird.data$Month_3),"NA",bird.data$y.3)

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
    },
    {bird.data$y.1<-as.integer(bird.data$y.1)
    bird.data$y.2<-as.integer(bird.data$y.2)
    bird.data$y.3<-as.integer(bird.data$y.3)
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
         y<-data.frame(bird.data$y.1,bird.data$y.2),
         y<-data.frame(bird.data$y.1,bird.data$y.2,bird.data$y.3)
         )
  
  #convert all values > 1 to 1
  y[y>1]<-1
  

  #define site covariates - FOR MULTIPLE YEARS
  site.covs<- data.frame(Year=as.factor(bird.data$Year), Unit_Code=as.factor(bird.data$Unit_Code), 
                         Plot_Name=as.factor(bird.data$Plot_Name))
  
  #create yearlySiteCovs matrices (Abun and Availability)
  ifelse(habitatType=="Forest",
         {
           Ord.Day=data.frame(bird.data$Ord.Day_1, bird.data$Ord.Day_2)
           Visit=data.frame(bird.data$Visit_1, bird.data$Visit_2)
           Time=data.frame(bird.data$Time_1, bird.data$Time_2)
           Temp=data.frame(bird.data$Temp_1, bird.data$Temp_2)
           Wind=data.frame(bird.data$Wind_1, bird.data$Wind_2) 
           Observer=data.frame(bird.data$Obs_1, bird.data$Obs_2)
         },
         {
           Ord.Day=data.frame(bird.data$Ord.Day_1, bird.data$Ord.Day_2, bird.data$Ord.Day_3)
           Visit=data.frame(bird.data$Visit_1, bird.data$Visit_2, bird.data$Visit_3)
           Time=data.frame(bird.data$Time_1, bird.data$Time_2, bird.data$Time_3)
           Temp=data.frame(bird.data$Temp_1, bird.data$Temp_2, bird.data$Temp_3)
           Wind=data.frame(bird.data$Wind_1, bird.data$Wind_2, bird.data$Wind_3) 
           Observer=data.frame(bird.data$Obs_1, bird.data$Obs_2, bird.data$Obs_3)
         }
         )

  #define abundance covariates (yearly site covs)
  obs.covs<-list(Visit=Visit, Ord.Day=Ord.Day, Time=Time, Temp=Temp, Wind=Wind, Observer=Observer)
  

  
  #set up unmarked FRAME GMM
         umf<-tryCatch({unmarkedFrameOccu(y=y, 
                               siteCovs = site.covs,
                               obsCovs= obs.covs)}
                       ,error=function(cond2){
                         cond2=unmarkedFramePCount(y=y, 
                                                siteCovs = site.covs,
                                                obsCovs= obs.covs)
                         cond2
                       })
         
         umf@obsCovs$Visit<-as.factor(umf@obsCovs$Visit)
         umf@obsCovs$Observer<-as.factor(umf@obsCovs$Observer)
         

         
  return(umf)
}

#End