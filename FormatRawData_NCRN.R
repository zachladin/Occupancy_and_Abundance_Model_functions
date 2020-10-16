FormatRawData_NCRN<-function(DataIn){
  
  data<-DataIn
  
  habitatType<-unique(as.character(data$Survey_Type))
  
  sort(unique(data$AOU_Code))
  
  #add in common names
  #read in BCI
  bci.data<-read.csv("/Users/zach/Dropbox (ZachTeam)/Projects/NPS/NCRN/Data/BCI_Guild_Assignments/BirdGuildAssignments_Appalachian.csv")
  bci.data$Common_Name<-gsub("Eastern Tufted Titmouse","Tufted Titmouse", bci.data$Common_Name)
  
  common.names<-unique(bci.data[,c("AOU_Code","Common_Name")])
  
  #add common names from BCI
  data<-merge(data, common.names, by=c("AOU_Code"),all.x=TRUE)
  
  #convert back to factor
  data$AOU_Code<-as.factor(as.character(data$AOU_Code))
  
  #format variables
  data$GRTS_Order<-as.factor(as.character(data$GRTS_Order))
  
  #add long name of parks
  unitCodeList<-sort(unique(as.character(data$Unit_Code)))
  
  #make data.frame with Unit_Code and Unit_Name
  ifelse(habitatType=="Forest",
         {unitCode.df<-data.frame(Unit_Code = unitCodeList, Unit_Name=c("Antietam National Battlefield", "Catoctin Mountain Park","Chesapeake and Ohio Canal National Historical Park","Greenbelt Park","George Washington Memorial Parkway","Harpers Ferry National Historical Park", " Manassas National Battlefield Park","Monocacy National Battlefield","National Capital Parks-East","Piscataway-Fort Washington", "Prince William Forest Park", "Rock Creek Park", "Wolf Trap National Park for the Performing Arts"))
         },
           {unitCode.df<-data.frame(Unit_Code = unitCodeList, Unit_Name=c("Antietam National Battlefield","Harpers Ferry National Historical Park", " Manassas National Battlefield Park","Monocacy National Battlefield"))
         })
  
  data<-merge(data, unitCode.df, by="Unit_Code",all.x=TRUE)
  
  #create Year, Month, Day, and Ordinal Day columns
  date.data<-data$Date
  ifelse(habitatType=="Forest",
         date.data<-as.character(date.data),
         date.data<-trimws(gsub("0:00:00","",as.character(date.data)))
         )
  
  date.string<-read.table(text=as.character(date.data), sep="/",colClasses = "character")
  colnames(date.string)<-c("Month","Day","Year")
  date.string$Year<-as.factor(paste("20",date.string$Year,sep=""))
  date.string$new.date<-as.character(paste(date.string$Year,date.string$Month,date.string$Day,sep="/"))
  date.string$Ord.Day <- as.integer(format(as.Date(date.string$new.date, format = "%Y/%m/%d"), "%j"))
  
  head(date.string)
  data<-cbind(data,date.string)
  
  #remove dates from Start_Time and End_Time
  data$Start_Time<-trimws(gsub("12/30/1899","",data$Start_Time))
  data$End_Time<-trimws(gsub("12/30/1899","",data$End_Time))
  
  
  time=data$Start_Time
  #convert Start_Time to decimal (numeric)
  timeConvertDecimal<-function(time){
    
    time.data<-as.character(time)
    time.data<-substr(time.data,1,4)
    
    time.new<-round(sapply(strsplit(time.data,":"),
                           function(x) {
                             x <- as.numeric(x)
                             x[1]+x[2]/60
                           }),3)
    return(time.new)
  }
  
  data$Start_Time.dec<-timeConvertDecimal(time=data$Start_Time)
  data$End_Time.dec<-timeConvertDecimal(time=data$End_Time)

  #scale numeric covariates
  data$Time.scale<-scale(data$Start_Time.dec)
  data$Temp.scale<-scale(as.numeric(data$Temperature))
  data$Humidity.scale<-scale(as.numeric(data$Humidity))
  
  #still need observer
  
  #create single Observer column
  ifelse(habitatType=="Forest",
         data$Observer<-paste(data$FirstOfLast_Name, data$FirstOfFirst_Name,sep="_"), 
         data$Observer<-as.character(data$ObserverLastName) 
         )
  data$Observer<-as.factor(data$Observer)
  
  #save as .csv file
  write.csv(data, file=paste("./Data/Bird Data/","NCRN_",habitatType,"_Bird_Data_Formatted.csv",sep=""),row.names=FALSE)

  message("Raw data has been formatted successfully and saved")
  return(data)
}
