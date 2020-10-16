
#fucntion to read in 3 files and combine into one file for a given year
combineAnnualDataFiles<-function(Dir){
  
  #get list of files in folder
  newFileList<- list.files(Dir, full.names=TRUE)
  
  #read in new data files
  new.file1<-read.csv(newFileList[1])
  new.file2<-read.csv(newFileList[2])
  new.file3<-read.csv(newFileList[3])
  
  #merge newfile1 and newfile3 using Unit_Date_Point_Visit column
  merge.file.1<-merge(new.file1, new.file3, by=c("Admin_Unit_Code","EventDate","Point_Name", "Visit"))
  
  #now merge plot location data
  merge.file.2<-merge(merge.file.1, new.file2, by=c("Admin_Unit_Code","Point_Name"))
  
  return(merge.file.2)  
  
}
