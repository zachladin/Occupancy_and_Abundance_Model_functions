#CheckAOS_Codes

#Function to check species names against names in data, update any names as needed, and generate table of changes that were made

CheckAOS_Codes<-function(DataIn){
  
  require(RCurl)
  require(curl)
  require(utils)
  require(foreign)
  require(reshape2)
  require(plyr)
  require(dplyr)
  require(stringr)
  
  Dir=paste(getwd(),"Data",sep="/")
  
  #create directory for saving updated species info
  message("Creating new folder to save updated species info.")
  dir.create(paste(Dir,"AOS_Codes",sep="/"))
  
  #download .zipped DBF file of up-to-date bird codes
  message("Downloading latest AOS codes from The Institute of Bird Populations (https://www.birdpop.org)")
  URL <- "https://www.birdpop.org/docs/misc/List18.zip"
  download.file(url=URL, destfile=paste(Dir,"AOS_Codes","List18.zip",sep="/"))
  
  #unzip AOS codes
  zipF<-paste(Dir,"AOS_Codes","List18.zip",sep="/")
  outDir<-paste(Dir,"AOS_Codes",sep="/")
  unzip(zipF,exdir=outDir)  # unzip your file 
  
  #read in up-to-date AOU_Codes
  species.codes<-read.dbf(file=paste(Dir,"AOS_Codes","LIST18.DBF",sep="/"))
  species.codes<-species.codes[,c("COMMONNAME","SCINAME","SPEC")]
  colnames(species.codes)<-c("English_Common_Name","Scientific_Name","AOU_Code")
  species.codes$Scientific_Name<-trimws(as.character(species.codes$Scientific_Name))
  species.codes$English_Common_Name<-trimws(as.character(species.codes$English_Common_Name))
  
  #Now add Order, Family, and Genus names
  message("Downloading latest list species detected in Breeding Bird Survey data from USGS (ftp://ftpext.usgs.gov)")

  BBSurl="ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/SpeciesList.txt"
  download.file(url=BBSurl, destfile=paste(Dir,"AOS_Codes","SpeciesList.txt",sep="/"))
  
  #get official list of all BBS species info from url
  raw_species = readLines(paste(Dir,"AOS_Codes","SpeciesList.txt",sep="/"))
  pasted_species = paste0(raw_species, collapse = "\n")
  
  # Convert to UTF-8 and then save.
  # The encoding wasn't actually latin1, it was something R couldn't handle.  But
  # on my machine at least, iconv seemed to accept it.
  write(
    iconv(gsub("\n([:alnum:])", "\1", pasted_species), from = "latin1", to = "UTF-8"), 
    file =paste(Dir,"AOS_Codes","SpeciesList_1.txt",sep="/")
  )
  
  dashes = readLines(paste(Dir,"AOS_Codes","SpeciesList_1.txt",sep="/"))[9]
  species_colnames = strsplit(readLines(paste(Dir,"AOS_Codes","SpeciesList_1.txt",sep="/"))[8]," +" )[[1]]
  
  species_list = read.fwf(
    paste(Dir,"AOS_Codes","SpeciesList_1.txt",sep="/"), 
    widths = nchar(strsplit(dashes, " ")[[1]]) + 1, 
    header = FALSE, 
    skip = 9,
    encoding = "UTF-8",
    stringsAsFactors = FALSE
  )
  colnames(species_list) = species_colnames
  
  #remove whitespace
  species_list$English_Common_Name<-trimws(as.character(species_list$English_Common_Name))

  #add Scientific_Name column
  species_list$Scientific_Name<-trimws(paste(trimws(species_list$Genus), trimws(species_list$Species),sep=" "))
  
  #merge species_list with read.dbf
  species.merge<-merge(species_list, species.codes, by="English_Common_Name",all.x=TRUE)
  
  #clean up species.merge
  species.merge$Scientific_Name.y<-NULL
  names(species.merge)[names(species.merge)=="Scientific_Name.x"]<-"Scientific_Name"
  names(species.merge)[names(species.merge)=="English_Common_Name"]<-"Common_Name"
  names(species.merge)[names(species.merge)=="ORDER"]<-"Order"


  #download PIF watchlist data
  message("Downloading latest Partners in Flight watchlist information from (http://pif.birdconservancy.org)")
  PIFurl ="http://pif.birdconservancy.org/ACAD/ajax/download.aspx?list=Glo&regions=US"
  download.file(url=PIFurl, destfile=paste(Dir,"AOS_Codes","PIF_watchlist.csv",sep="/"))
  
  #read in pif.data
  pif.data<-read.csv(paste(Dir,"AOS_Codes", "PIF_watchlist.csv",sep="/"),quote = "", 
                     row.names = NULL,stringsAsFactors = FALSE)

  pif.data.2<-as.data.frame(sapply(pif.data, function(x) gsub("\"", "", x)))
  pif.data.2$row.names<-NULL

  #fix column names
  pif.names<-names(pif.data.2)[-1]
  names(pif.data.2)<-pif.names
  
  #get desired columns
  pif.data.3<-pif.data.2[,c("common_name","continental_importance","iucn_red_list_2016")]
  colnames(pif.data.3)<-c("Common_Name","Continental_Concern","IUCN_Red_List_2016")
  pif.data.3$Common_Name<-trimws(as.character(pif.data.3$Common_Name))
  
  #merge pif.data with species.merge
  species.out<-merge(species.merge, pif.data.3, by="Common_Name",all.x=TRUE)
  
  #save species info to .csv file
  write.csv(species.out, paste(Dir,"AOS_Codes","SpeciesList_out.csv",sep="/"),row.names=FALSE)  

  
  data<-NCRN@Birds
  
  
  
  #simplify data (get list of unique species)
  unique.birds<-unique(data[,c("AOU_Code","Common_Name")])
  
  #get list of unique AOU_Codes in data
  uniqueAOU<-data.frame(AOU_Code=sort(unique(unique.birds$AOU_Code)))
  
  #remove blanks and any "unidentified" codes
  uniqueAOU<-subset(uniqueAOU, AOU_Code !="")
  
  #species to remove
  removeList<-c("UNBI","UNCH","UNHA","UNCR","UNDU","UNFL","UNOW","UNSP","UNSW","UNTH","UNWA","UNWO","UNWR")
  
  uniqueAOU<-subset(uniqueAOU, ! AOU_Code %in% removeList)
  
  uniqueAOU.merge<-merge(uniqueAOU, species.out, by="AOU_Code",all.x=TRUE)
  
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
