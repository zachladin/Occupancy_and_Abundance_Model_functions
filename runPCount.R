runPCount<-function(DataIn, modSel){
  
  data<-DataIn
  
  #change Survey_Type to Habitat
  names(data)[names(data)=="Survey_Type"]<-"Habitat"
  
  #create habitatList
  habitatList<-unique(as.character(data$Habitat))
  
  habitat.save<-list()
  for(i in 1:length(habitatList)){
    
    #set habitatType
    habitatType<-as.character(habitatList[i])
    
    #subset data by habitat
    new.data<-subset(data, Habitat==habitatType)
    
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
          
          #create unmarkedFramePCount_umf
          #make umf
          new.umf<-tryCatch({
            suppressWarnings(makePCount_umf(DataIn=year.data, SpeciesIn=speciesName))
          },error=function(cond2){
            cond2=NULL
          })
          
          #Fit models
          #set all models to equal NULL to clear previous results
          mod1=mod2=mod3=mod4=mod5=mod6=mod7=NULL
          
          #use convert.units=0.01 (0.01*0.01) to convert from m^2 to ha
          
          #model 1 - Null - P
          mod1<-tryCatch({pcount(~1 ~1, data=new.umf, mixture="P")
          }, error= function(out){
            out=NULL
          })
          
          #model 2 - Null - NB (det)
          mod2<-tryCatch({pcount(~1 ~1, data=new.umf, mixture="NB")
          }, error=function(out){
            out=NULL
          })
          
          
          #model 3 - Observer - P (det)
          mod3<-tryCatch({pcount(~Observer ~1, data=new.umf, mixture="P")
          }, error=function(out){
            out=NULL
          })
          
          #model 4 - Observer - NB (det) 
          mod4<-tryCatch({pcount(~Observer ~ 1, data=new.umf, mixture="NB")
          }, error=function(out){
            out=NULL
          })
          
          
          #model 5 - Observer+Ord.Day - P
          mod5<-tryCatch({pcount(~Observer+Ord.Day ~ 1, data=new.umf, mixture="P")
          }, error=function(out){
            out=NULL
          })
          
          #model 6 - Observer+Ord.Day - NB 
          mod6<-tryCatch({pcount(~Observer+Ord.Day ~ 1, data=new.umf, mixture="NB")
          }, error=function(out){
            out=NULL
          })
          
          ####################################################################################################################
          #get detection probablilities from each model and sort to maximize detection
          
          try({
            mod.1.detect.A<-unique(suppressWarnings(predict(mod1, type="det")))
            mod.1.detect.B<-mod.1.detect.A$Predicted[1]
            
            mod.2.detect.A<-unique(suppressWarnings(predict(mod2, type="det")))
            mod.2.detect.A$One_minus_pred<-1-mod.2.detect.A$Predicted
            mod.2.detect.B<-1-prod(mod.2.detect.A$One_minus_pred)
            
            mod.3.detect.A<-unique(suppressWarnings(predict(mod3, type="det")))
            mod.3.detect.A$One_minus_pred<-1-mod.3.detect.A$Predicted
            mod.3.detect.B<-1-prod(mod.3.detect.A$One_minus_pred)
            
            mod.4.detect.A<-unique(suppressWarnings(predict(mod4, type="det")))
            mod.4.detect.A$One_minus_pred<-1-mod.4.detect.A$Predicted
            mod.4.detect.B<-1-prod(mod.4.detect.A$One_minus_pred)
            
            mod.5.detect.A<-unique(suppressWarnings(predict(mod5, type="det")))
            mod.5.detect.A$One_minus_pred<-1-mod.5.detect.A$Predicted
            mod.5.detect.B<-1-prod(mod.5.detect.A$One_minus_pred)
            
            mod.6.detect.A<-unique(suppressWarnings(predict(mod6, type="det")))
            mod.6.detect.A$One_minus_pred<-1-mod.6.detect.A$Predicted
            mod.6.detect.B<-1-prod(mod.6.detect.A$One_minus_pred)
            
          },silent=TRUE)
          
          #get detection prob mod1
          mod.1.det<-tryCatch({
            mod.1.det<-data.frame(ModNames="Null_P",Det=mod.1.detect.B)
          },error=function(cond2){
            cond2=data.frame(ModNames="Null_P",Det=NA)
          })
          
          #get detection prob mod2
          mod.2.det<-tryCatch({
            mod.2.det<-data.frame(ModNames="Null_NB",Det=mod.2.detect.B)
          },error=function(cond2){
            cond2=data.frame(ModNames="Null_NB",Det=NA)
          })
          
          #get detection prob mod3
          mod.3.det<-tryCatch({
            mod.3.det<-data.frame(ModNames="Observer_P",Det=mod.3.detect.B)
          },error=function(cond2){
            cond2=data.frame(ModNames="Observer_P",Det=NA)
          })
          
          #get detection prob mod4
          mod.4.det<-tryCatch({
            mod.4.det<-data.frame(ModNames="Observer_NB",Det=mod.4.detect.B)
          },error=function(cond2){
            cond2=data.frame(ModNames="Observer_NB",Det=NA)
          })
          
          #get detection prob mod5
          mod.5.det<-tryCatch({
            mod.5.det<-data.frame(ModNames="Observer.Day_P",Det=mod.5.detect.B)
          },error=function(cond2){
            cond2=data.frame(ModNames="Observer.Day_P",Det=NA)
          })
          
          #get detection prob mod6
          mod.6.det<-tryCatch({
            mod.6.det<-data.frame(ModNames="Observer.Day_NB",Det=mod.6.detect.B)
          },error=function(cond2){
            cond2=data.frame(ModNames="Observer.Day_NB",Det=NA)
          })
          
          #compile detection results
          mod.det.results<-rbind(mod.1.det, mod.2.det, mod.3.det, mod.4.det, mod.5.det, mod.6.det)
          
          #remove NAs
          mod.det.results.noNA<-na.omit(mod.det.results)
          
          #order by det high to low
          mod.det.results.order<-mod.det.results.noNA[order(mod.det.results.noNA$Det, decreasing=TRUE),]
          
          #mod to keep
          mod.max.det<-as.character(mod.det.results.order$ModNames[1])
          
          #######################################################################################################################
          #choose model names
          Modnames <- c("Null_P","Null_NB","Observer_P", "Observer_NB",  "Observer.Day_P","Observer.Day_NB") 
          
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
          
          aic.summary.df.4<-tryCatch({aictab(cand.set = list(mod1), modnames = unlist(Modnames[-c(2, 3, 4, 5, 6)]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod2), modnames = unlist(Modnames[-c(1, 3, 4, 5, 6)]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod3), modnames = unlist(Modnames[-c(1, 2, 4, 5, 6)]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod4), modnames = unlist(Modnames[-c(1, 2, 3, 5, 6)]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod5), modnames = unlist(Modnames[-c(1, 2, 3, 4, 6)]))
          },error = function(out){out=tryCatch({aictab(cand.set = list(mod6), modnames = unlist(Modnames[-c(1, 2, 3, 4, 5)]))
          },error = function(out){ 
            out=NULL
          })})})})})})
          
          
          ifelse(!is.null(aic.summary.df.1), aic.summary.df.5<-aic.summary.df.1,
                 ifelse(!is.null(aic.summary.df.2), aic.summary.df.5<-aic.summary.df.2,
                        ifelse(!is.null(aic.summary.df.3), aic.summary.df.5<-aic.summary.df.3,
                               ifelse(!is.null(aic.summary.df.4), aic.summary.df.5<-aic.summary.df.4,
                                      aic.summary.df.5<-data.frame("Modnames"=NA, "K"=NA, "AICc"=NA, "Delta_AICc"=NA,"ModelLik"=NA, 
                                                                   "AICcWt"=NA, "LL"=NA,"Cum.Wt"=NA)))))
          
          
          #select top performing model
          #take model with lowest AIC score
          mod.keep<-aic.summary.df.5[1,]
          
          #make data.frame of stuff to add to model output to keep track of things
          species.data<-subset(year.data, AOU_Code==speciesName)
          keep.info<-unique(species.data[,c("Year","Unit_Code","AOU_Code")])[1,]
          
          #add keep.info to model output
          mod.keep.2<-cbind(keep.info, mod.keep)
          
          #Let "modSel" function parameter determine whether to use maxDet or AIC to choose model predicted estimates
          mod.name<-switch(modSel,
                           "maxDet"=mod.max.det,
                           "AIC"=as.character(mod.keep.2$Modnames),
                           mod.max.det)
          
          Modnames
          
          tryCatch({
            ifelse(habitatType=="Forest",
                 
          {       
          #get abundance estimates from top model
          abun.out<-NULL
          ifelse(mod.name=="Null_P",abun.out<-predict(mod1,type="state",appendData=TRUE),
                 ifelse(mod.name=="Null_NB", abun.out<-predict(mod2,type="state",appendData=TRUE),
                        ifelse(mod.name=="Observer_P", abun.out<-predict(mod3,type="state",appendData=TRUE),
                               ifelse(mod.name=="Observer_NB", abun.out<-predict(mod4,type="state",appendData=TRUE),
                                      ifelse(mod.name=="Observer.Day_P", abun.out<-predict(mod5,type="state",appendData=TRUE),
                                             ifelse(mod.name=="Observer.Day_NB", abun.out<-predict(mod6,type="state",appendData=TRUE),
                                                    abun.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,Year=yearName,Unit_Code=unitName,GRTS_Order=NA,Visit.1=NA,Visit.2=NA,Ord.Day.1=NA,Ord.Day.2=NA,Time.1=NA,Time.2=NA, Temp.1=NA,Temp.2=NA, Wind.1=NA,Wind.2=NA,Observer.1=NA,Observer.2=NA)
                                             ))))))
          

          #get p(detection) estimates from top model
          det.out<-NULL
          ifelse(mod.name=="Null_P",det.out<-predict(mod1,type="det",appendData=TRUE),
                 ifelse(mod.name=="Null_NB", det.out<-predict(mod2,type="det",appendData=TRUE),
                        ifelse(mod.name=="Observer_P", det.out<-predict(mod3,type="det",appendData=TRUE),
                               ifelse(mod.name=="Observer_NB", det.out<-predict(mod4,type="det",appendData=TRUE),
                                      ifelse(mod.name=="Observer.Day_P", det.out<-predict(mod5,type="det",appendData=TRUE),
                                             ifelse(mod.name=="Observer.Day_NB", det.out<-predict(mod6,type="det",appendData=TRUE),
                                                    det.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,Year=yearName,Unit_Code=unitName,GRTS_Order=NA,Visit.1=NA,Visit.2=NA,Ord.Day.1=NA,Ord.Day.2=NA,Time.1=NA,Time.2=NA, Temp.1=NA,Temp.2=NA, Wind.1=NA,Wind.2=NA,Observer.1=NA,Observer.2=NA)
                                             ))))))
          
          
          #get model formula
          formula.out<-NULL
          ifelse(mod.name=="Null_P",formula.out<-paste(deparse(mod1@call)[1],sep=""),
                 ifelse(mod.name=="Null_NB", formula.out<-paste(deparse(mod2@call)[1],sep=""),
                        ifelse(mod.name=="Observer_P", formula.out<-paste(deparse(mod3@call)[1],sep=""),
                               ifelse(mod.name=="Observer_NB", formula.out<-paste(deparse(mod4@call)[1],sep=""),
                                      ifelse(mod.name=="Observer.Day_P", formula.out<-paste(deparse(mod5@call)[1],sep=""),
                                             ifelse(mod.name=="Observer.Day_NB", formula.out<-paste(deparse(mod6@call)[1],sep="") ,
                                                    formula.out<-NA
                                             ))))))
          
          
          
          
          
          
          #condense results
          #get model output
          abun.df<-unique(abun.out[,-5:-6])
          abun.df$Metric<-"Abundance"
          abun.df$AOU_Code<-speciesName
          
          det<-unique(det.out[,-5:-6])
          det$Metric<-"Detection"
          det$AOU_Code<-speciesName
          
          #combine output
          results.comb<-rbind(abun.df, det)
          
          #add model
          results.comb$Model.Formula<-formula.out
          },
          
          {
          #get abundance estimates from top model
          abun.out<-NULL
          ifelse(mod.name=="Null_P",abun.out<-predict(mod1,type="state",appendData=TRUE),
                 ifelse(mod.name=="Null_NB", abun.out<-predict(mod2,type="state",appendData=TRUE),
                        ifelse(mod.name=="Observer_P", abun.out<-predict(mod3,type="state",appendData=TRUE),
                               ifelse(mod.name=="Observer_NB", abun.out<-predict(mod4,type="state",appendData=TRUE),
                                      ifelse(mod.name=="Observer.Day_P", abun.out<-predict(mod5,type="state",appendData=TRUE),
                                             ifelse(mod.name=="Observer.Day_NB", abun.out<-predict(mod6,type="state",appendData=TRUE),
                                                    abun.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,y.3=NA,Year=yearName,Unit_Code=unitName,GRTS_Order=NA,Visit.1=NA,Visit.2=NA,Visit.3=NA,Ord.Day.1=NA,Ord.Day.2=NA,Ord.Day.3=NA,Time.1=NA,Time.2=NA,Time.3=NA, Temp.1=NA,Temp.2=NA,Temp.3=NA, Wind.1=NA,Wind.2=NA,Wind.3=NA,Observer.1=NA,Observer.2=NA,Observer.3=NA)
                                             ))))))
          
          
          #get p(detection) estimates from top model
          det.out<-NULL
          ifelse(mod.name=="Null_P",det.out<-predict(mod1,type="det",appendData=TRUE),
                 ifelse(mod.name=="Null_NB", det.out<-predict(mod2,type="det",appendData=TRUE),
                        ifelse(mod.name=="Observer_P", det.out<-predict(mod3,type="det",appendData=TRUE),
                               ifelse(mod.name=="Observer_NB", det.out<-predict(mod4,type="det",appendData=TRUE),
                                      ifelse(mod.name=="Observer.Day_P", det.out<-predict(mod5,type="det",appendData=TRUE),
                                             ifelse(mod.name=="Observer.Day_NB", det.out<-predict(mod6,type="det",appendData=TRUE),
                                                    det.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,y.3=NA,Year=yearName,Unit_Code=unitName,GRTS_Order=NA,Visit.1=NA,Visit.2=NA,Visit.3=NA,Ord.Day.1=NA,Ord.Day.2=NA,Ord.Day.3=NA,Time.1=NA,Time.2=NA,Time.3=NA, Temp.1=NA,Temp.2=NA,Temp.3=NA, Wind.1=NA,Wind.2=NA,Wind.3=NA,Observer.1=NA,Observer.2=NA,Observer.3=NA)
                                             ))))))
          
          
          #get model formula
          formula.out<-NULL
          ifelse(mod.name=="Null_P",formula.out<-paste(deparse(mod1@call)[1],sep=""),
                 ifelse(mod.name=="Null_NB", formula.out<-paste(deparse(mod2@call)[1],sep=""),
                        ifelse(mod.name=="Observer_P", formula.out<-paste(deparse(mod3@call)[1],sep=""),
                               ifelse(mod.name=="Observer_NB", formula.out<-paste(deparse(mod4@call)[1],sep=""),
                                      ifelse(mod.name=="Observer.Day_P", formula.out<-paste(deparse(mod5@call)[1],sep=""),
                                             ifelse(mod.name=="Observer.Day_NB", formula.out<-paste(deparse(mod6@call)[1],sep="") ,
                                                    formula.out<-NA
                                             ))))))
          
          
          
          
          
          
          #condense results
          #get model output
          abun.df<-unique(abun.out[,-5:-7])
          abun.df$Metric<-"Abundance"
          abun.df$AOU_Code<-speciesName
          
          det<-unique(det.out[,-5:-7])
          det$Metric<-"Detection"
          det$AOU_Code<-speciesName
          
          #combine output
          results.comb<-rbind(abun.df, det)
          
          #add model
          results.comb$Model.Formula<-formula.out
        })
          },error=function(cond2){
            cond2=NULL
          })
          
          #reorganize df
          results.out<-tryCatch({
            data.frame(Habitat=habitatType,Unit_Code=unitName, Year=yearName, AOU_Code=speciesName, Metric=results.comb$Metric, Predicted=results.comb$Predicted, SE=results.comb$SE, lower=results.comb$lower, upper=results.comb$upper)
          },error=function(cond2){
            cond2=data.frame(Habitat=habitatType,Unit_Code=unitName, Year=yearName, AOU_Code=speciesName, Metric=NA, Predicted=NA, SE=NA, lower=NA, upper=NA)
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