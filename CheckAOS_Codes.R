#CheckAOS_Codes

#Function to check species names against names in data, update any names as needed, and generate table of changes that were made

CheckAOS_Codes<-function(DataIn){
  
  require(foreign)
  require(reshape2)
  require(plyr)
  require(dplyr)
  require(stringr)
  
  new.data<-DataIn
  
  #read in up-to-date AOU_Codes
  species.codes<-read.dbf(file=paste(getwd(), "Data","AOS_Codes","LIST17.dbf",sep="/"))
  species.codes<-species.codes[,c("COMMONNAME","SCINAME","SPEC")]
  colnames(species.codes)<-c("English_Common_Name","Scientific_Name","AOU_Code")
  species.codes$Scientific_Name<-trimws(as.character(species.codes$Scientific_Name))
  
  #get official list of all BBS species info from url
  raw_species = readLines(paste(getwd(),"Data","AOS_Codes","SpeciesList.txt",sep="/"))
  pasted_species = paste0(raw_species, collapse = "\n")
  
  # Convert to UTF-8 and then save.
  # The encoding wasn't actually latin1, it was something R couldn't handle.  But
  # on my machine at least, iconv seemed to accept it.
  write(
    iconv(gsub("\n([:alnum:])", "\1", pasted_species), from = "latin1", to = "UTF-8"), 
    file =paste(getwd(),"Data","AOS_Codes","SpeciesList_1.txt",sep="/")
  )
  
  dashes = readLines(paste(getwd(),"Data","AOS_Codes","SpeciesList_1.txt",sep="/"))[9]
  species_colnames = strsplit(
    readLines(paste(getwd(),"Data","AOS_Codes","SpeciesList_1.txt",sep="/"))[8],
    " +"
  )[[1]]
  
  species_list = read.fwf(
    paste(getwd(),"Data","AOS_Codes","SpeciesList_1.txt",sep="/"), 
    widths = nchar(strsplit(dashes, " ")[[1]]) + 1, 
    header = FALSE, 
    skip = 9,
    encoding = "UTF-8",
    stringsAsFactors = FALSE
  )
  colnames(species_list) = species_colnames
  
  #add Scientific_Name column
  species_list$Scientific_Name<-trimws(paste(trimws(species_list$Genus), trimws(species_list$Species),sep=" "))
  
  #merge species_list with read.dbf
  species.merge<-merge(species_list, species.codes, by="Scientific_Name",all.x=TRUE)
  
  write.csv(species.merge, paste(getwd(),"Data","AOS_Codes","SpeciesList_out.csv",sep="/"),row.names=FALSE)  
  
  #get list of unique AOU_Codes in data
  uniqueAOU<-data.frame(AOU_Code=sort(unique(new.data$AOU_Code)))
  
  #remove blanks and any "unidentified" codes
  uniqueAOU<-subset(uniqueAOU, AOU_Code !="")
  
  removeList<-c("UNBI","UNCH","UNHA","UNCR","UNDU","UNFL","UNOW","UNSP","UNSW","UNTH","UNWA","UNWO","UNWR")
  
  uniqueAOU<-subset(uniqueAOU, ! AOU_Code %in% removeList)
  
  uniqueAOU.merge<-merge(uniqueAOU, species.merge, by="AOU_Code",all.x=TRUE)
  
  #find rows that did not match up (with NAs)
  NA.species<-subset(uniqueAOU.merge, is.na(Scientific_Name))
  #for looping over
  NA.species.list<-as.character(NA.species$AOU_Code)
  #for printing message
  NA.species.list.print<-list(as.character(NA.species$AOU_Code))
  
  message("There was a problem with the following species: ", NA.species.list.print)
  message("For updated codes, visit https://www.birdpop.org/docs/misc/Alpha_codes_tax.pdf.")
  
  new.codes <- readline("Enter new codes using uppercase letters and delimit by commas with no spaces): ")
  
  new.codes.2<-t(read.table(text=as.character(new.codes), sep=","))

  #Reconcile changes in species names
  for(i in 1:length(new.codes.2)){
    temp.code<-new.codes.2[i]
    new.data$AOU_Code<-gsub(as.character(NA.species.list[[i]]), as.character(temp.code), new.data$AOU_Code)
  }
  
  #save table (.csv file) for inclusion in report
  
  #create file for report materials
  dir.create(paste(getwd(), "QAP_report",sep="/"),showWarnings = FALSE)
  
  #create .csv file of changes made
  spp.report<-data.frame(Old_Species_Code=unlist(NA.species.list), Updated_Species_Code=unlist(new.codes.2), Date=Sys.Date())
  
  row.names(spp.report)<-NULL
  
  spp.report.2<-spp.report
  spp.report.2$Old_Species_Code<-as.character(spp.report.2$Old_Species_Code)
  spp.report.2$Updated_Species_Code<-as.character(spp.report.2$Updated_Species_Code)
  
  spp.report.2<- spp.report.2 %>% rowwise() %>% mutate(Modified=!identical(Old_Species_Code,Updated_Species_Code))
  
  write.csv(spp.report.2, file=paste(getwd(),"QAP_report","spp_report.csv",sep="/"),row.names=FALSE)
  
  message("4-letter species AOS code(s) have been updated, and report has been saved.")
  
  return(new.data)
  }
