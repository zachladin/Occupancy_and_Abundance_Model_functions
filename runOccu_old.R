runOccu<-function(DataIn){
  
  data<-DataIn
  
  #create habitatList
  habitatList<-unique(as.character(data$Survey_Type))
  
  habitat.save<-list()
  for(i in 1:length(habitatList)){
    
    #set habitatType
    habitatType<-as.character(habitatList[i])
    
    #subset data by Survey_Type
    new.data<-subset(data, Survey_Type==habitatType)
    
    #create unitList
    unitList<-as.character(unique(new.data$Unit_Code))
    
    #now unit-level
    unit.save<-list()
    for(j in 1:length(unitList)){
      
      unitName<-as.character(unitList[j])
      
      #get Unit data
      unit.data<-subset(new.data, Unit_Code==unitName)
      
      #create yearList
      yearList<-sort(as.character(unique(unit.data$Year)))
      
      year.save<-list()
      for(k in 1:length(yearList)){
        
        yearName<-yearList[k]
        
        year.data<-subset(unit.data, Year==yearName)
        
        #get speciesList
        speciesList<-sort(as.character(unique(year.data$AOU_Code)))

        species.save<-list()
        for(r in 1:length(speciesList)){
          
          speciesName<-speciesList[r]
          
          #create unmarkedFramePOccu_umf
          #make umf
          new.umf<-tryCatch({
            makePOccu_umf(DataIn=year.data, SpeciesIn=speciesName)
          },error=function(cond2){
            cond2=NULL
          })
          
          #Fit models
          #set all models to equal NULL to clear previous results
          mod1=mod2=mod3=mod4=mod5=mod6=mod7=NULL
          
          #use convert.units=0.01 (0.01*0.01) to convert from m^2 to ha
          
          #model 1 - Null
          mod1<-tryCatch({occu(~1 ~1, data=new.umf)
          }, error= function(out){
            out=NULL
          })
          
          #model 2 - Observer (det)
          mod2<-tryCatch({occu(~Observer ~1, data=new.umf)
          }, error=function(out){
            out=NULL
          })
          
          
          #model 3 - Observer+Ord.Day (det)
          mod3<-tryCatch({occu(~Observer+Ord.Day ~1, data=new.umf)
          }, error=function(out){
            out=NULL
          })
          
          #model 4 - Observer+Ord.Day+Time (det) 
          mod4<-tryCatch({occu(~Observer+Ord.Day+Time ~ 1, data=new.umf)
          }, error=function(out){
            out=NULL
          })
          
          
          #model 5 - Observer+Ord.Day+Time+Temp (det)
          mod5<-tryCatch({occu(~Observer+Ord.Day+Time+Temp ~ 1, data=new.umf)
          }, error=function(out){
            out=NULL
          })
          
          #model 6 - Observer+Ord.Day+Time+Temp+Wind (det)
          mod6<-tryCatch({occu(~Observer+Ord.Day+Time+Temp+Wind ~ 1, data=new.umf)
          }, error=function(out){
            out=NULL
          })
          
          #######################################################################################################################
          #choose model names
          Modnames <- c("Null","Observer","Observer_Ord.Day", "Observer_Ord.Day_Time", "Observer_Ord.Day_Time_Temp","Observer_Ord.Day_Time_Temp_Wind") 
          
          #crazy error-handling section!
          aic.summary.df.1<-tryCatch({aictab(cand.set = list(mod1,mod2,mod3,mod4,mod5,mod6), modnames = unlist(Modnames))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod2,mod3,mod4,mod5), modnames = unlist(Modnames[-6]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod2,mod3,mod4,mod6), modnames = unlist(Modnames[-5]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod2,mod3,mod5,mod6), modnames = unlist(Modnames[-4]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod2,mod4,mod5,mod6), modnames = unlist(Modnames[-3]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod3,mod4,mod5,mod6), modnames = unlist(Modnames[-2]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod2,mod3,mod4,mod5,mod6), modnames = unlist(Modnames[-1]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod2,mod3,mod4), modnames = unlist(Modnames[-c(5, 6)]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod2,mod3,mod5), modnames = unlist(Modnames[-c(4, 6)]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod2,mod3,mod6), modnames = unlist(Modnames[-c(4, 5)]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod2,mod4,mod5), modnames = unlist(Modnames[-c(3, 6)]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod2,mod4,mod6), modnames = unlist(Modnames[-c(3, 5)]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod2,mod5,mod6), modnames = unlist(Modnames[-c(3, 4)]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod3,mod4,mod5), modnames = unlist(Modnames[-c(2, 6)]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod3,mod4,mod6), modnames = unlist(Modnames[-c(2, 5)]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod3,mod5,mod6), modnames = unlist(Modnames[-c(2, 4)]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod4,mod5,mod6), modnames = unlist(Modnames[-c(2, 3)]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod2,mod3,mod4,mod5), modnames = unlist(Modnames[-c(1, 6)]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod2,mod3,mod4,mod6), modnames = unlist(Modnames[-c(1, 5)]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod2,mod3,mod5,mod6), modnames = unlist(Modnames[-c(1, 4)]))
          },error = function(out){
            out=NULL
          })})})})})})})})})})})})})})})})})})})})
          
             
                 aic.summary.df.2<-tryCatch({aictab(cand.set = list(mod2,mod4,mod5,mod6), modnames = unlist(Modnames[-c(1, 3)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod3,mod4,mod5,mod6), modnames = unlist(Modnames[-c(1, 2)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod2,mod3), modnames = unlist(Modnames[-c(4, 5, 6)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod2,mod4), modnames = unlist(Modnames[-c(3, 5, 6)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod2,mod5), modnames = unlist(Modnames[-c(3, 4, 6)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod2,mod6), modnames = unlist(Modnames[-c(3, 4, 5)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod3,mod4), modnames = unlist(Modnames[-c(2, 5, 6)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod3,mod5), modnames = unlist(Modnames[-c(2, 4, 6)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod3,mod6), modnames = unlist(Modnames[-c(2, 4, 5)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod4,mod5), modnames = unlist(Modnames[-c(2, 3, 6)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod1,mod4,mod6), modnames = unlist(Modnames[-c(2, 3, 5)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod2, mod3,mod4), modnames = unlist(Modnames[-c(1, 5, 6)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod2, mod3,mod5), modnames = unlist(Modnames[-c(1, 4, 6)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod2, mod3,mod6), modnames = unlist(Modnames[-c(1, 4, 5)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod2, mod4,mod5), modnames = unlist(Modnames[-c(1, 3, 6)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod2, mod4,mod6), modnames = unlist(Modnames[-c(1, 3, 5)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod2, mod5,mod6), modnames = unlist(Modnames[-c(1, 3, 4)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod3, mod4,mod5), modnames = unlist(Modnames[-c(1, 2, 6)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod3, mod4,mod6), modnames = unlist(Modnames[-c(1, 2, 5)]))
                 },error = function(out){
                   out=NULL
                 })})})})})})})})})})})})})})})})})})})
                

                 aic.summary.df.3<-tryCatch({aictab(cand.set = list(mod3, mod5,mod6), modnames = unlist(Modnames[-c(1, 2, 4)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod4, mod5,mod6), modnames = unlist(Modnames[-c(1, 2, 3)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod1, mod2), modnames = unlist(Modnames[-c(3, 4, 5, 6)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod1, mod3), modnames = unlist(Modnames[-c(2, 4, 5, 6)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod1, mod4), modnames = unlist(Modnames[-c(2, 3, 5, 6)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod1, mod5), modnames = unlist(Modnames[-c(2, 3, 4, 6)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod1, mod6), modnames = unlist(Modnames[-c(2, 3, 4, 5)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod2, mod3), modnames = unlist(Modnames[-c(1, 4, 5, 6)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod2, mod4), modnames = unlist(Modnames[-c(1, 3, 5, 6)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod2, mod5), modnames = unlist(Modnames[-c(1, 3, 4, 6)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod2, mod6), modnames = unlist(Modnames[-c(1, 3, 4, 5)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod3, mod4), modnames = unlist(Modnames[-c(1, 2, 5, 6)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod3, mod5), modnames = unlist(Modnames[-c(1, 2, 4, 6)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod3, mod6), modnames = unlist(Modnames[-c(1, 2, 4, 5)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod4, mod5), modnames = unlist(Modnames[-c(1, 2, 3, 6)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod4, mod6), modnames = unlist(Modnames[-c(1, 2, 3, 5)]))
                 },error = function(out){out=tryCatch({aictab(cand.set = list(mod5, mod6), modnames = unlist(Modnames[-c(1, 2, 3, 4)]))
                 },error = function(out){
                   out=NULL
                 })})})})})})})})})})})})})})})})})
                 
                 
                 
                 ifelse(!is.null(aic.summary.df.1), aic.summary.df.4<-aic.summary.df.1,
                                          ifelse(!is.null(aic.summary.df.2), aic.summary.df.4<-aic.summary.df.2,
                                                 ifelse(!is.null(aic.summary.df.3), aic.summary.df.4<-aic.summary.df.3,
                                                        aic.summary.df.4<-data.frame("Modnames"=NA, "K"=NA, "AICc"=NA, "Delta_AICc"=NA,"ModelLik"=NA, 
                                                                   "AICcWt"=NA, "LL"=NA,"Cum.Wt"=NA))))
                                          
                                          
          #select top performing model
          #take model with lowest AIC score
          mod.keep<-aic.summary.df.4[1,]
          
          #make data.frame of stuff to add to model output to keep track of things
          species.data<-subset(year.data, AOU_Code==speciesName)
          keep.info<-unique(species.data[,c("Year","Unit_Code","AOU_Code")])[1,]
          
          #add keep.info to model output
          mod.keep.2<-cbind(keep.info, mod.keep)
          
          #get model name
          mod.name<-as.character(mod.keep.2$Modnames)
          
          Modnames
          
          tryCatch({
            ifelse(habitatType=="Forest",
                 
          {       
          #get occupancy estimates from top model
          occu.out<-NULL
          ifelse(mod.name=="Null",occu.out<-predict(mod1,type="state",appendData=TRUE),
                 ifelse(mod.name=="Observer", occu.out<-predict(mod2,type="state",appendData=TRUE),
                        ifelse(mod.name=="Observer_Ord.Day", occu.out<-predict(mod3,type="state",appendData=TRUE),
                               ifelse(mod.name=="Observer_Ord.Day_Time", occu.out<-predict(mod4,type="state",appendData=TRUE),
                                      ifelse(mod.name=="Observer_Ord.Day_Time_Temp", occu.out<-predict(mod5,type="state",appendData=TRUE),
                                             ifelse(mod.name=="Observer_Ord.Day_Time_Temp_Wind", occu.out<-predict(mod6,type="state",appendData=TRUE),
                                                    occu.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,Year=yearName,Unit_Code=unitName,Plot_Name=NA,Visit.1=NA,Visit.2=NA,Ord.Day.1=NA,Ord.Day.2=NA,Time.1=NA,Time.2=NA, Temp.1=NA,Temp.2=NA, Wind.1=NA,Wind.2=NA,Observer.1=NA,Observer.2=NA)
                                             ))))))
          

          #get p(detection) estimates from top model
          det.out<-NULL
          ifelse(mod.name=="Null",det.out<-predict(mod1,type="det",appendData=TRUE),
                 ifelse(mod.name=="Observer", det.out<-predict(mod2,type="det",appendData=TRUE),
                        ifelse(mod.name=="Observer_Ord.Day", det.out<-predict(mod3,type="det",appendData=TRUE),
                               ifelse(mod.name=="Observer_Ord.Day_Time", det.out<-predict(mod4,type="det",appendData=TRUE),
                                      ifelse(mod.name=="Observer_Ord.Day_Time_Temp", det.out<-predict(mod5,type="det",appendData=TRUE),
                                             ifelse(mod.name=="Observer_Ord.Day_Time_Temp_Wind", det.out<-predict(mod6,type="det",appendData=TRUE),
                                                    det.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,Year=yearName,Unit_Code=unitName,Plot_Name=NA,Visit.1=NA,Visit.2=NA,Ord.Day.1=NA,Ord.Day.2=NA,Time.1=NA,Time.2=NA, Temp.1=NA,Temp.2=NA, Wind.1=NA,Wind.2=NA,Observer.1=NA,Observer.2=NA)
                                             ))))))
          
          
          #get model formula
          formula.out<-NULL
          ifelse(mod.name=="Null",formula.out<-paste(deparse(mod1@call)[1],sep=""),
                 ifelse(mod.name=="Observer", formula.out<-paste(deparse(mod2@call)[1],sep=""),
                        ifelse(mod.name=="Observer_Ord.Day", formula.out<-paste(deparse(mod3@call)[1],sep=""),
                               ifelse(mod.name=="Observer_Ord.Day_Time", formula.out<-paste(deparse(mod4@call)[1],sep=""),
                                      ifelse(mod.name=="Observer_Ord.Day_Time_Temp", formula.out<-paste(deparse(mod5@call)[1],sep=""),
                                             ifelse(mod.name=="Observer_Ord.Day_Time_Temp_Wind", formula.out<-paste(deparse(mod6@call)[1],sep="") ,
                                                    formula.out<-NA
                                             ))))))
          
          
          
          
          
          
          #condense results
          #get model output
          occu.df<-unique(occu.out[,-5:-6])
          occu.df$Metric<-"Occupancy"
          occu.df$AOU_Code<-speciesName
          
          det<-unique(det.out[,-5:-6])
          det$Metric<-"Detection"
          det$AOU_Code<-speciesName
          
          #combine output
          results.comb<-rbind(occu.df, det)
          
          #add model
          results.comb$Model.Formula<-formula.out
          },
          
          {
          #get occupancy estimates from top model
          occu.out<-NULL
          ifelse(mod.name=="Null",occu.out<-predict(mod1,type="state",appendData=TRUE),
                 ifelse(mod.name=="Observer", occu.out<-predict(mod2,type="state",appendData=TRUE),
                        ifelse(mod.name=="Observer_Ord.Day", occu.out<-predict(mod3,type="state",appendData=TRUE),
                               ifelse(mod.name=="Observer_Ord.Day_Time", occu.out<-predict(mod4,type="state",appendData=TRUE),
                                      ifelse(mod.name=="Observer_Ord.Day_Time_Temp", occu.out<-predict(mod5,type="state",appendData=TRUE),
                                             ifelse(mod.name=="Observer_Ord.Day_Time_Temp_Wind", occu.out<-predict(mod6,type="state",appendData=TRUE),
                                                    occu.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,y.3=NA,Year=yearName,Unit_Code=unitName,Plot_Name=NA,Visit.1=NA,Visit.2=NA,Visit.3=NA,Ord.Day.1=NA,Ord.Day.2=NA,Ord.Day.3=NA,Time.1=NA,Time.2=NA,Time.3=NA, Temp.1=NA,Temp.2=NA,Temp.3=NA, Wind.1=NA,Wind.2=NA,Wind.3=NA,Observer.1=NA,Observer.2=NA,Observer.3=NA)
                                             ))))))
          
          
          #get p(detection) estimates from top model
          det.out<-NULL
          ifelse(mod.name=="Null",det.out<-predict(mod1,type="det",appendData=TRUE),
                 ifelse(mod.name=="Observer", det.out<-predict(mod2,type="det",appendData=TRUE),
                        ifelse(mod.name=="Observer_Ord.Day", det.out<-predict(mod3,type="det",appendData=TRUE),
                               ifelse(mod.name=="Observer_Ord.Day_Time", det.out<-predict(mod4,type="det",appendData=TRUE),
                                      ifelse(mod.name=="Observer_Ord.Day_Time_Temp", det.out<-predict(mod5,type="det",appendData=TRUE),
                                             ifelse(mod.name=="Observer_Ord.Day_Time_Temp_Wind", det.out<-predict(mod6,type="det",appendData=TRUE),
                                                    det.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,y.3=NA,Year=yearName,Unit_Code=unitName,Plot_Name=NA,Visit.1=NA,Visit.2=NA,Visit.3=NA,Ord.Day.1=NA,Ord.Day.2=NA,Ord.Day.3=NA,Time.1=NA,Time.2=NA,Time.3=NA, Temp.1=NA,Temp.2=NA,Temp.3=NA, Wind.1=NA,Wind.2=NA,Wind.3=NA,Observer.1=NA,Observer.2=NA,Observer.3=NA)
                                             ))))))
          
          
          #get model formula
          formula.out<-NULL
          ifelse(mod.name=="Null",formula.out<-paste(deparse(mod1@call)[1],sep=""),
                 ifelse(mod.name=="Observer", formula.out<-paste(deparse(mod2@call)[1],sep=""),
                        ifelse(mod.name=="Observer_Ord.Day", formula.out<-paste(deparse(mod3@call)[1],sep=""),
                               ifelse(mod.name=="Observer_Ord.Day_Time", formula.out<-paste(deparse(mod4@call)[1],sep=""),
                                      ifelse(mod.name=="Observer_Ord.Day_Time_Temp", formula.out<-paste(deparse(mod5@call)[1],sep=""),
                                             ifelse(mod.name=="Observer_Ord.Day_Time_Temp_Wind", formula.out<-paste(deparse(mod6@call)[1],sep="") ,
                                                    formula.out<-NA
                                             ))))))
          
          
          
          
          
          
          #condense results
          #get model output
          occu.df<-unique(occu.out[,-5:-7])
          occu.df$Metric<-"Occupancy"
          occu.df$AOU_Code<-speciesName
          
          det<-unique(det.out[,-5:-7])
          det$Metric<-"Detection"
          det$AOU_Code<-speciesName
          
          #combine output
          results.comb<-rbind(occu.df, det)
          
          #add model
          results.comb$Model.Formula<-formula.out
        })
          },error=function(cond2){
            cond2=NULL
          })
          
          #reorganize df
          results.out<-tryCatch({
            data.frame(Survey_Type=habitatType,Unit_Code=unitName, Year=yearName, AOU_Code=speciesName, Metric=results.comb$Metric, Predicted=results.comb$Predicted, SE=results.comb$SE, lower=results.comb$lower, upper=results.comb$upper)
          },error=function(cond2){
            cond2=data.frame(Survey_Type=habitatType,Unit_Code=unitName, Year=yearName, AOU_Code=speciesName, Metric=NA, Predicted=NA, SE=NA, lower=NA, upper=NA)
          })
            
          #get unique rows
          results.out.unique<-unique(results.out)
          
          
          #add AICc info back in 
          results.out.unique<-data.frame(results.out.unique, mod.keep.2[-1:-3])
          
          
          #fill park.output list with each iteration through park loop
          print(paste("Saving results for",habitatType, unitName, yearName, speciesName,sep=" "))
          species.save<-rbind(species.save, results.out.unique)
        }
        
        
        #fill species.output list with each iteration through the species loop
        print(paste("Saving results for",habitatType, unitName, yearName,sep=" "))
        year.save<-rbind(year.save, species.save)
        
      }
      
      #fill year.output list with each iteration throught the year loop
      print(paste("Saving results for",habitatType, unitName,sep=" "))
      unit.save<-rbind(unit.save, year.save)
      
    }
    
    #fill year.output list with each iteration throught the year loop
    print(paste("Saving results for",habitatType,sep=" "))
    habitat.save<-rbind(habitat.save, unit.save)
    
  }
  
  
  return(habitat.save)
}