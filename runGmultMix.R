runGmultMix<-function(DataIn){
  
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
          
          #create unmarkedFrameGMM
          #make umf
          #make umf
          new.umf<-tryCatch({
            makeGMMumf(DataIn=year.data, SpeciesIn=speciesName)
          },error=function(cond2){
            cond2=NULL
          })
          
          
          #Fit models
          #set all models to equal NULL to clear previous results
          mod1=mod2=mod3=mod4=mod5=mod6=mod7=NULL
          
          #use convert.units=0.01 (0.01*0.01) to convert from m^2 to ha
          
          #model 1 - Null, Negative Binomial
          mod1<-tryCatch({gmultmix(~1, ~1, ~1, data=new.umf, mixture="NB")
          }, error= function(out){
            out=NULL
          })
          
          #model 2 - Null, Poisson 
          mod2<-tryCatch({gmultmix(~1, ~1, ~1, data=new.umf, mixture="P")
          }, error=function(out){
            out=NULL
          })
          
          
          #model 3 - Observer (det), Negative Binomial
          mod3<-tryCatch({gmultmix(~1, ~1, ~Observer, data=new.umf, mixture="NB")
          }, error=function(out){
            out=NULL
          })
          
          #model 4 - Observer (det), Poisson
          mod4<-tryCatch({gmultmix(~1, ~1, ~Observer, data=new.umf, mixture="P")
          }, error=function(out){
            out=NULL
          })
          
          
          #model 5 - Observer (det) and Ordinal Day (avail), Negative Binomial
          mod5<-tryCatch({gmultmix(~1, ~Ord.Day, ~Observer, data=new.umf, mixture="NB")
          }, error=function(out){
            out=NULL
          })
          
          #model 6 - Observer (det) and Ordinal Day (avail), Poisson
          mod6<-tryCatch({gmultmix(~1, ~Ord.Day, ~Observer, data=new.umf, mixture="P")
          }, error=function(out){
            out=NULL
          })
          
          # #########################################
          # #get all possible combinations of mods1-6
          # comb6<-t(combn(x=c(1,2,3,4,5,6), m=6, FUN = NULL, simplify = TRUE))
          # comb6.out<-cbind(paste("mod",comb6[,1],sep=""), paste("mod",comb6[,2],sep=""), paste("mod",comb6[,3],sep=""), paste("mod",comb6[,4],sep=""),paste("mod",comb6[,5],sep=""),paste("mod",comb6[,6],sep=""))
          # modlist.1<-list(comb6.out)
          # 
          # comb5<-t(combn(x=c(1,2,3,4,5,6), m=5, FUN = NULL, simplify = TRUE))
          # comb5.out<-cbind(paste("mod",comb5[,1],sep=""), paste("mod",comb5[,2],sep=""), paste("mod",comb5[,3],sep=""), paste("mod",comb5[,4],sep=""),paste("mod",comb5[,5],sep=""))
          # 
          # comb4<-t(combn(x=c(1,2,3,4,5,6), m=4, FUN = NULL, simplify = TRUE))
          # comb4.out<-cbind(paste("mod",comb4[,1],sep=""), paste("mod",comb4[,2],sep=""), paste("mod",comb4[,3],sep=""), paste("mod",comb4[,4],sep=""))
          # 
          # comb3<-t(combn(x=c(1,2,3,4,5,6), m=3, FUN = NULL, simplify = TRUE))
          # comb3.out<-cbind(paste("mod",comb3[,1],sep=""), paste("mod",comb3[,2],sep=""), paste("mod",comb3[,3],sep=""))
          # 
          # comb2<-t(combn(x=c(1,2,3,4,5,6), m=2, FUN = NULL, simplify = TRUE))
          # comb2.out<-cbind(paste("mod",comb2[,1],sep=""), paste("mod",comb2[,2],sep=""))
          # 
          # comb1<-t(combn(x=c(1,2,3,4,5,6), m=1, FUN = NULL, simplify = TRUE))
          # comb1.out<-paste("mod",comb1,sep="")
          # #########################################
          
          ############################################################################################################################
          #choose model names
          Modnames <- c("Null_NB","Null_P","Observer_NB", "Observer_P", "Observer.Day_NB","Observer.Day_P") 
          
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
          #get abundance estimates from top model
          abun.out<-NULL
          ifelse(mod.name=="Null_NB",abun.out<-predict(mod1,type="lambda",appendData=TRUE),
                 ifelse(mod.name=="Null_P", abun.out<-predict(mod2,type="lambda",appendData=TRUE),
                        ifelse(mod.name=="Observer_NB", abun.out<-predict(mod3,type="lambda",appendData=TRUE),
                               ifelse(mod.name=="Observer_P", abun.out<-predict(mod4,type="lambda",appendData=TRUE),
                                      ifelse(mod.name=="Observer.Day_NB", abun.out<-predict(mod5,type="lambda",appendData=TRUE),
                                             ifelse(mod.name=="Observer.Day_P", abun.out<-predict(mod6,type="lambda",appendData=TRUE),
                                                    abun.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,y.3=NA,y.4=NA,y.5=NA,y.6=NA,           y.7=NA,y.8=NA,Year=yearName,Unit_Code=unitName,Plot_Name=NA,Visit.obs.1=NA,Visit.obs.2=NA,Visit.obs.3=NA,Visit.obs.4=NA,Visit.obs.5=NA,   Visit.obs.6=NA,Visit.obs.7=NA,Visit.obs.8=NA,Ord.Day.obs.1=NA,Ord.Day.obs.2=NA,Ord.Day.obs.3=NA,Ord.Day.obs.4=NA,Ord.Day.obs.5=NA,Ord.Day.obs.6=NA, Ord.Day.obs.7=NA,Ord.Day.obs.8=NA,Time.obs.1=NA,Time.obs.2=NA,Time.obs.3=NA,Time.obs.4=NA,Time.obs.5=NA,Time.obs.6=NA,Time.obs.7=NA,Time.obs.8=NA,    Temp.1=NA,Temp.2=NA, Temp.3=NA,Temp.4=NA,Temp.5=NA,Temp.6=NA, Temp.7=NA,Temp.8=NA, Wind.1=NA,Wind.2=NA,Wind.3=NA,Wind.4,Wind.5=NA,Wind.6=NA,Wind.7=NA, Wind.8=NA,Observer.1=NA,Observer.2=NA,Observer.3=NA,Observer.4=NA,Observer.5=NA,Observer.6=NA,Observer.7=NA,Observer.8=NA )
                                             ))))))
          
          #get p(availability) estimates from top model
          avail.out<-NULL
          ifelse(mod.name=="Null_NB",avail.out<-predict(mod1,type="phi",appendData=TRUE),
                 ifelse(mod.name=="Null_P", avail.out<-predict(mod2,type="phi",appendData=TRUE),
                        ifelse(mod.name=="Observer_NB", avail.out<-predict(mod3,type="phi",appendData=TRUE),
                               ifelse(mod.name=="Observer_P", avail.out<-predict(mod4,type="phi",appendData=TRUE),
                                      ifelse(mod.name=="Observer.Day_NB", avail.out<-predict(mod5,type="phi",appendData=TRUE),
                                             ifelse(mod.name=="Observer.Day_P", avail.out<-predict(mod6,type="phi",appendData=TRUE),
                                                    avail.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,y.3=NA,y.4=NA,y.5=NA,y.6=NA,y.7=NA,y.8=NA,Year=yearName,Unit_Code=unitName,Plot_Name=NA,Visit.obs.1=NA,Visit.obs.2=NA,Visit.obs.3=NA,Visit.obs.4=NA,Visit.obs.5=NA,Visit.obs.6=NA,Visit.obs.7=NA,Visit.obs.8=NA,Ord.Day.obs.1=NA,Ord.Day.obs.2=NA,Ord.Day.obs.3=NA,Ord.Day.obs.4=NA,Ord.Day.obs.5=NA,Ord.Day.obs.6=NA, Ord.Day.obs.7=NA,Ord.Day.obs.8=NA,Time.obs.1=NA,Time.obs.2=NA,Time.obs.3=NA,Time.obs.4=NA,Time.obs.5=NA,Time.obs.6=NA,Time.obs.7=NA,Time.obs.8=NA,Temp.1=NA,Temp.2=NA, Temp.3=NA,Temp.4=NA,Temp.5=NA,Temp.6=NA, Temp.7=NA,Temp.8=NA, Wind.1=NA,Wind.2=NA,Wind.3=NA,Wind.4,Wind.5=NA,Wind.6=NA,Wind.7=NA, Wind.8=NA,Observer.1=NA,Observer.2=NA,Observer.3=NA,Observer.4=NA,Observer.5=NA,Observer.6=NA,Observer.7=NA,Observer.8=NA )
                                             ))))))
          
          #get p(availability) estimates from top model
          det.out<-NULL
          ifelse(mod.name=="Null_NB",det.out<-predict(mod1,type="det",appendData=TRUE),
                 ifelse(mod.name=="Null_P", det.out<-predict(mod2,type="det",appendData=TRUE),
                        ifelse(mod.name=="Observer_NB", det.out<-predict(mod3,type="det",appendData=TRUE),
                               ifelse(mod.name=="Observer_P", det.out<-predict(mod4,type="det",appendData=TRUE),
                                      ifelse(mod.name=="Observer.Day_NB", det.out<-predict(mod5,type="det",appendData=TRUE),
                                             ifelse(mod.name=="Observer.Day_P", det.out<-predict(mod6,type="det",appendData=TRUE),
                                                    det.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,y.3=NA,y.4=NA,y.5=NA,y.6=NA,y.7=NA,y.8=NA,Year=yearName,Unit_Code=unitName,Plot_Name=NA,Visit.obs.1=NA,Visit.obs.2=NA,Visit.obs.3=NA,Visit.obs.4=NA,Visit.obs.5=NA,Visit.obs.6=NA,Visit.obs.7=NA,Visit.obs.8=NA,Ord.Day.obs.1=NA,Ord.Day.obs.2=NA,Ord.Day.obs.3=NA,Ord.Day.obs.4=NA,Ord.Day.obs.5=NA,Ord.Day.obs.6=NA, Ord.Day.obs.7=NA,Ord.Day.obs.8=NA,Time.obs.1=NA,Time.obs.2=NA,Time.obs.3=NA,Time.obs.4=NA,Time.obs.5=NA,Time.obs.6=NA,Time.obs.7=NA,Time.obs.8=NA,Temp.1=NA,Temp.2=NA, Temp.3=NA,Temp.4=NA,Temp.5=NA,Temp.6=NA, Temp.7=NA,Temp.8=NA, Wind.1=NA,Wind.2=NA,Wind.3=NA,Wind.4,Wind.5=NA,Wind.6=NA,Wind.7=NA, Wind.8=NA,Observer.1=NA,Observer.2=NA,Observer.3=NA,Observer.4=NA,Observer.5=NA,Observer.6=NA,Observer.7=NA,Observer.8=NA )
                                             ))))))
          
          
          #get model formula
          formula.out<-NULL
          ifelse(mod.name=="Null_NB",formula.out<-paste(deparse(mod1@call)[1],trimws(deparse(mod1@call)[2]),sep=""),
                 ifelse(mod.name=="Null_P", formula.out<-paste(deparse(mod2@call)[1],trimws(deparse(mod2@call)[2]),sep=""),
                        ifelse(mod.name=="Observer_NB", formula.out<-paste(deparse(mod3@call)[1],trimws(deparse(mod3@call)[2]),sep=""),
                               ifelse(mod.name=="Observer_P", formula.out<-paste(deparse(mod4@call)[1],trimws(deparse(mod4@call)[2]),sep=""),
                                      ifelse(mod.name=="Observer.Day_NB", formula.out<-paste(deparse(mod5@call)[1],trimws(deparse(mod5@call)[2]),sep=""),
                                             ifelse(mod.name=="Observer.Day_P", formula.out<-paste(deparse(mod6@call)[1],trimws(deparse(mod6@call)[2]),sep="") ,
                                                    formula.out<-NA
                                             ))))))
          
          
          
          
          
          
          #condense results
          #get model output
          abun<-unique(abun.out[,-5:-12])
          abun$Metric<-"Abundance"
          abun$AOU_Code<-speciesName
          
          
          avail<-unique(avail.out[,-5:-12])
          avail$Metric<-"Availability"
          avail$AOU_Code<-speciesName
          
          
          det<-unique(det.out[,-5:-12])
          det$Metric<-"Detection"
          det$AOU_Code<-speciesName
          
          #combine output
          results.comb<-rbind(abun, avail, det)
          
          #add model
          results.comb$Model.Formula<-formula.out
          },
          
          {
            #get abundance estimates from top model
            abun.out<-NULL
            ifelse(mod.name=="Null_NB",abun.out<-predict(mod1,type="lambda",appendData=TRUE),
                   ifelse(mod.name=="Null_P", abun.out<-predict(mod2,type="lambda",appendData=TRUE),
                          ifelse(mod.name=="Observer_NB", abun.out<-predict(mod3,type="lambda",appendData=TRUE),
                                 ifelse(mod.name=="Observer_P", abun.out<-predict(mod4,type="lambda",appendData=TRUE),
                                        ifelse(mod.name=="Observer.Day_NB", abun.out<-predict(mod5,type="lambda",appendData=TRUE),
                                               ifelse(mod.name=="Observer.Day_P", abun.out<-predict(mod6,type="lambda",appendData=TRUE),
                                                      abun.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,y.3=NA,y.4=NA,y.5=NA,y.6=NA,y.7=NA,y.8=NA,y.9=NA,y.10=NA,y.11=NA,y.12=NA,Year=yearName,Unit_Code=unitName,Plot_Name=NA,Visit.obs.1=NA,Visit.obs.2=NA,Visit.obs.3=NA,Visit.obs.4=NA,Visit.obs.5=NA,   Visit.obs.6=NA,Visit.obs.7=NA,Visit.obs.8=NA,Visit.obs.9=NA,Visit.obs.10=NA,Visit.obs.11=NA,Visit.obs.12=NA,Ord.Day.obs.1=NA,Ord.Day.obs.2=NA,Ord.Day.obs.3=NA,Ord.Day.obs.4=NA,Ord.Day.obs.5=NA,Ord.Day.obs.6=NA, Ord.Day.obs.7=NA,Ord.Day.obs.8=NA,Ord.Day.obs.9=NA,Ord.Day.obs.10=NA,Ord.Day.obs.11=NA,Ord.Day.obs.12=NA,Time.obs.1=NA,Time.obs.2=NA,Time.obs.3=NA,Time.obs.4=NA,Time.obs.5=NA,Time.obs.6=NA,Time.obs.7=NA,Time.obs.8=NA,Time.obs.9=NA,Time.obs.10=NA,Time.obs.11=NA,Time.obs.12=NA,Temp.1=NA,Temp.2=NA, Temp.3=NA,Temp.4=NA,Temp.5=NA,Temp.6=NA, Temp.7=NA,Temp.8=NA,Temp.9=NA,Temp.10=NA,Temp.11=NA,Temp.12=NA, Wind.1=NA,Wind.2=NA,Wind.3=NA,Wind.4,Wind.5=NA,Wind.6=NA,Wind.7=NA, Wind.8=NA, Wind.9=NA, Wind.10=NA, Wind.11=NA, Wind.12=NA,Observer.1=NA,Observer.2=NA,Observer.3=NA,Observer.4=NA,Observer.5=NA,Observer.6=NA,Observer.7=NA,Observer.8=NA,Observer.9=NA,Observer.10=NA,Observer.11=NA,Observer.12=NA )
                                               ))))))
            
            #get p(availability) estimates from top model
            avail.out<-NULL
            ifelse(mod.name=="Null_NB",avail.out<-predict(mod1,type="phi",appendData=TRUE),
                   ifelse(mod.name=="Null_P", avail.out<-predict(mod2,type="phi",appendData=TRUE),
                          ifelse(mod.name=="Observer_NB", avail.out<-predict(mod3,type="phi",appendData=TRUE),
                                 ifelse(mod.name=="Observer_P", avail.out<-predict(mod4,type="phi",appendData=TRUE),
                                        ifelse(mod.name=="Observer.Day_NB", avail.out<-predict(mod5,type="phi",appendData=TRUE),
                                               ifelse(mod.name=="Observer.Day_P", avail.out<-predict(mod6,type="phi",appendData=TRUE),
                                                      avail.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,y.3=NA,y.4=NA,y.5=NA,y.6=NA,y.7=NA,y.8=NA,y.9=NA,y.10=NA,y.11=NA,y.12=NA,Year=yearName,Unit_Code=unitName,Plot_Name=NA,Visit.obs.1=NA,Visit.obs.2=NA,Visit.obs.3=NA,Visit.obs.4=NA,Visit.obs.5=NA,   Visit.obs.6=NA,Visit.obs.7=NA,Visit.obs.8=NA,Visit.obs.9=NA,Visit.obs.10=NA,Visit.obs.11=NA,Visit.obs.12=NA,Ord.Day.obs.1=NA,Ord.Day.obs.2=NA,Ord.Day.obs.3=NA,Ord.Day.obs.4=NA,Ord.Day.obs.5=NA,Ord.Day.obs.6=NA, Ord.Day.obs.7=NA,Ord.Day.obs.8=NA,Ord.Day.obs.9=NA,Ord.Day.obs.10=NA,Ord.Day.obs.11=NA,Ord.Day.obs.12=NA,Time.obs.1=NA,Time.obs.2=NA,Time.obs.3=NA,Time.obs.4=NA,Time.obs.5=NA,Time.obs.6=NA,Time.obs.7=NA,Time.obs.8=NA,Time.obs.9=NA,Time.obs.10=NA,Time.obs.11=NA,Time.obs.12=NA,Temp.1=NA,Temp.2=NA, Temp.3=NA,Temp.4=NA,Temp.5=NA,Temp.6=NA, Temp.7=NA,Temp.8=NA,Temp.9=NA,Temp.10=NA,Temp.11=NA,Temp.12=NA, Wind.1=NA,Wind.2=NA,Wind.3=NA,Wind.4,Wind.5=NA,Wind.6=NA,Wind.7=NA, Wind.8=NA, Wind.9=NA, Wind.10=NA, Wind.11=NA, Wind.12=NA,Observer.1=NA,Observer.2=NA,Observer.3=NA,Observer.4=NA,Observer.5=NA,Observer.6=NA,Observer.7=NA,Observer.8=NA,Observer.9=NA,Observer.10=NA,Observer.11=NA,Observer.12=NA)
                                               ))))))
            
            #get p(detection) estimates from top model
            det.out<-NULL
            ifelse(mod.name=="Null_NB",det.out<-predict(mod1,type="det",appendData=TRUE),
                   ifelse(mod.name=="Null_P", det.out<-predict(mod2,type="det",appendData=TRUE),
                          ifelse(mod.name=="Observer_NB", det.out<-predict(mod3,type="det",appendData=TRUE),
                                 ifelse(mod.name=="Observer_P", det.out<-predict(mod4,type="det",appendData=TRUE),
                                        ifelse(mod.name=="Observer.Day_NB", det.out<-predict(mod5,type="det",appendData=TRUE),
                                               ifelse(mod.name=="Observer.Day_P", det.out<-predict(mod6,type="det",appendData=TRUE),
                                                      det.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,y.3=NA,y.4=NA,y.5=NA,y.6=NA,y.7=NA,y.8=NA,y.9=NA,y.10=NA,y.11=NA,y.12=NA,Year=yearName,Unit_Code=unitName,Plot_Name=NA,Visit.obs.1=NA,Visit.obs.2=NA,Visit.obs.3=NA,Visit.obs.4=NA,Visit.obs.5=NA,   Visit.obs.6=NA,Visit.obs.7=NA,Visit.obs.8=NA,Visit.obs.9=NA,Visit.obs.10=NA,Visit.obs.11=NA,Visit.obs.12=NA,Ord.Day.obs.1=NA,Ord.Day.obs.2=NA,Ord.Day.obs.3=NA,Ord.Day.obs.4=NA,Ord.Day.obs.5=NA,Ord.Day.obs.6=NA, Ord.Day.obs.7=NA,Ord.Day.obs.8=NA,Ord.Day.obs.9=NA,Ord.Day.obs.10=NA,Ord.Day.obs.11=NA,Ord.Day.obs.12=NA,Time.obs.1=NA,Time.obs.2=NA,Time.obs.3=NA,Time.obs.4=NA,Time.obs.5=NA,Time.obs.6=NA,Time.obs.7=NA,Time.obs.8=NA,Time.obs.9=NA,Time.obs.10=NA,Time.obs.11=NA,Time.obs.12=NA,Temp.1=NA,Temp.2=NA, Temp.3=NA,Temp.4=NA,Temp.5=NA,Temp.6=NA, Temp.7=NA,Temp.8=NA,Temp.9=NA,Temp.10=NA,Temp.11=NA,Temp.12=NA, Wind.1=NA,Wind.2=NA,Wind.3=NA,Wind.4,Wind.5=NA,Wind.6=NA,Wind.7=NA, Wind.8=NA, Wind.9=NA, Wind.10=NA, Wind.11=NA, Wind.12=NA,Observer.1=NA,Observer.2=NA,Observer.3=NA,Observer.4=NA,Observer.5=NA,Observer.6=NA,Observer.7=NA,Observer.8=NA,Observer.9=NA,Observer.10=NA,Observer.11=NA,Observer.12=NAs)
                                               ))))))
            
            
            #get model formula
            formula.out<-NULL
            ifelse(mod.name=="Null_NB",formula.out<-paste(deparse(mod1@call)[1],trimws(deparse(mod1@call)[2]),sep=""),
                   ifelse(mod.name=="Null_P", formula.out<-paste(deparse(mod2@call)[1],trimws(deparse(mod2@call)[2]),sep=""),
                          ifelse(mod.name=="Observer_NB", formula.out<-paste(deparse(mod3@call)[1],trimws(deparse(mod3@call)[2]),sep=""),
                                 ifelse(mod.name=="Observer_P", formula.out<-paste(deparse(mod4@call)[1],trimws(deparse(mod4@call)[2]),sep=""),
                                        ifelse(mod.name=="Observer.Day_NB", formula.out<-paste(deparse(mod5@call)[1],trimws(deparse(mod5@call)[2]),sep=""),
                                               ifelse(mod.name=="Observer.Day_P", formula.out<-paste(deparse(mod6@call)[1],trimws(deparse(mod6@call)[2]),sep="") ,
                                                      formula.out<-NA
                                               ))))))
            
            
            #condense results
            #get model output
            abun<-unique(abun.out[,-5:-16])
            abun$Metric<-"Abundance"
            abun$AOU_Code<-speciesName
            
            avail<-unique(avail.out[,-5:-16])
            avail$Metric<-"Availability"
            avail$AOU_Code<-speciesName
            
            det<-unique(det.out[,-5:-16])
            det$Metric<-"Detection"
            det$AOU_Code<-speciesName
            
            #combine output
            results.comb<-rbind(abun, avail, det)
            
            #add model
            results.comb$Model.Formula<-formula.out
            }
          )
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