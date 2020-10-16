runBCI<-function(DataIn){
  
    data<-DataIn
  
    #subset data by Survey_Type (Forest only)
    new.data<-subset(data, Survey_Type=="Forest")
    
    #create unitList
    unitList<-as.character(unique(new.data$Unit_Code))
    
    #now unit-level
    unit.save<-list()
    for(j in 1:length(unitList)){
      
      
      #get Unit data
      unit.data<-subset(new.data, Unit_Code==unitList[j])
      
      #create yearList
      yearList<-sort(as.character(unique(unit.data$Year)))
      
      year.save<-list()
      for(k in 1:length(yearList)){
        
        
        yearName<-yearList[k]
        
        year.data<-subset(unit.data, Year==yearName)
        
        unitName<-as.character(unique(year.data$Unit_Code))
        
        bci.out<-BCIappalachian(DataIn=year.data, unitName=unitName, yearName=yearName)
  
        print(paste("Saving BCI results for", unitName, yearName,sep=" "))
        year.save<-rbind(year.save, bci.out)
      }
      print(paste("Saving BCI results for", unitName,sep=" "))
      unit.save<-rbind(unit.save, year.save)
    }
  
  return(unit.save)
}