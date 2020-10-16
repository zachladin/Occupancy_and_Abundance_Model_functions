runOccu<-function(DataIn, modSel){
  require(unmarked)
  setMethod("coerce",c("NULL","data.frame"), function(from, to, strict=TRUE) as.data.frame(from))
  
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
          
          #create unmarkedFramePOccu_umf
          #make umf
          new.umf<-NULL
          new.umf<-tryCatch({
            suppressWarnings(makePOccu_umf(DataIn=year.data, SpeciesIn=speciesName))
          },error=function(cond2){
            cond2=NULL
          })
          
          
          my.umf<-as(new.umf,"data.frame")
          ifelse(length(unique(year.data$Visit))==1,
                 covs.test<-data.frame(myCov=na.exclude(my.umf$Temp)),
                 covs.test<-data.frame(myCov=na.exclude(my.umf$Temp.1)))
          
          
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
          
          #model 3 - Ord.Day (det)
          mod3<-tryCatch({occu(~Ord.Day ~1, data=new.umf)
          }, error=function(out){
            out=NULL
          })
          
          #model 4 - Observer+Ord.Day (det) 
          mod4<-tryCatch({occu(~Observer+Ord.Day ~ 1, data=new.umf)
          }, error=function(out){
            out=NULL
          })
          
          #model 5 - Observer+Ord.Day+Time (det)
          mod5<-tryCatch({occu(~Observer+Ord.Day+Time ~ 1, data=new.umf)
          }, error=function(out){
            out=NULL})

          #model 6 - Observer+Ord.Day+Time+Wind (det)
          mod6<-tryCatch({occu(~Observer+Ord.Day+Wind ~ 1, data=new.umf)
          }, error=function(out){
            out=NULL})

          ###################################################################################################
          #NOTE: Need to figure out way to add Temp, and avoid situations where All NAs in Temp columns causes R to                   crash
          #model 5 - Observer+Ord.Day+Time+Wind ifelse (check if covariates have all NAs)
         # if(length(covs.test$myCov) == 0)
         #         {mod5<-tryCatch({occu(~Observer+Ord.Day+Time+Wind ~ 1, data=new.umf)
         #         }, error=function(out){
         #           out=NULL})
         #         }else(
         #         {mod5<-tryCatch({occu(~Observer+Ord.Day+Time+Temp ~ 1, data=new.umf)
         #         }, error=function(out){
         #           out=NULL})
         #         })
         # 
         #  #model 6 - Observer+Ord.Day+Time+Temp+Wind (det)
         #  if(length(covs.test$myCov) == 0)
         #         {mod6<-tryCatch({occu(~Observer+Ord.Day+Time+Wind ~ 1, data=new.umf)
         #         }, error=function(out){
         #           out=NULL})
         #         }else(
         #         {mod6<-tryCatch({occu(~Observer+Ord.Day+Time+Temp~ 1, data=new.umf)
         #         }, error=function(out){
         #           out=NULL})
         #         })
          ###################################################################################################
          
   ####################################################################################################################
          #get detection probablilities from each model and sort to maximize detection
          
          try({
            
            ifelse(!is.null(mod1),
                   {
                     mod.1.detect.A<-NULL
                     mod.1.detect.A<-unique(suppressWarnings(predict(mod1, type="det")))
                     mod.1.detect.B<-mod.1.detect.A$Predicted[1]
                   },
                   mod.1.detect.B<-NA)
          
            ifelse(!is.null(mod2),
                   {
          mod.2.detect.A<-NULL
          mod.2.detect.A<-unique(suppressWarnings(predict(mod2, type="det")))
          mod.2.detect.A$One_minus_pred<-1-mod.2.detect.A$Predicted
          mod.2.detect.B<-1-prod(mod.2.detect.A$One_minus_pred, na.rm=TRUE)
                   },
          mod.2.detect.B<-NA)
            
            ifelse(!is.null(mod3),
                   {
          mod.3.detect.A<-NULL
          mod.3.detect.A<-unique(suppressWarnings(predict(mod3, type="det")))
          mod.3.detect.A$One_minus_pred<-1-mod.3.detect.A$Predicted
          mod.3.detect.B<-1-prod(mod.3.detect.A$One_minus_pred, na.rm=TRUE)
                   },
          mod.3.detect.B<-NA)
            
            ifelse(!is.null(mod4),
                   {
          mod.4.detect.A<-NULL
          mod.4.detect.A<-unique(suppressWarnings(predict(mod4, type="det")))
          mod.4.detect.A$One_minus_pred<-1-mod.4.detect.A$Predicted
          mod.4.detect.B<-1-prod(mod.4.detect.A$One_minus_pred, na.rm=TRUE)
            },
          mod.4.detect.B<-NA)
            
            ifelse(!is.null(mod5),
                   {
          mod.5.detect.A<-NULL
          mod.5.detect.A<-unique(suppressWarnings(predict(mod5, type="det")))
          mod.5.detect.A$One_minus_pred<-1-mod.5.detect.A$Predicted
          mod.5.detect.B<-1-prod(mod.5.detect.A$One_minus_pred, na.rm=TRUE)
            },
          mod.5.detect.B<-NA)

            ifelse(!is.null(mod6),
                   {
          mod.6.detect.A<-NULL
          mod.6.detect.A<-unique(suppressWarnings(predict(mod6, type="det")))
          mod.6.detect.A$One_minus_pred<-1-mod.6.detect.A$Predicted
          mod.6.detect.B<-1-prod(mod.6.detect.A$One_minus_pred, na.rm=TRUE)
            },
          mod.6.detect.B<-NA)
          },silent=TRUE)
          
          
          #get detection prob mod1
          mod.1.det<-tryCatch({
            mod.1.det<-data.frame(ModNames="Null",Det=mod.1.detect.B)
          },error=function(cond2){
            cond2=data.frame(ModNames="Null",Det=NA)
          })
          
          #get detection prob mod2
          mod.2.det<-tryCatch({
            mod.2.det<-data.frame(ModNames="Observer",Det=mod.2.detect.B)
          },error=function(cond2){
            cond2=data.frame(ModNames="Observer",Det=NA)
          })
          
          #get detection prob mod3
          mod.3.det<-tryCatch({
            mod.3.det<-data.frame(ModNames="Ord.Day",Det=mod.3.detect.B)
          },error=function(cond2){
            cond2=data.frame(ModNames="Ord.Day",Det=NA)
          })
          
          #get detection prob mod4
          mod.4.det<-tryCatch({
            mod.4.det<-data.frame(ModNames="Observer_Ord.Day",Det=mod.4.detect.B)
          },error=function(cond2){
            cond2=data.frame(ModNames="Observer_Ord.Day",Det=NA)
          })
          
          #get detection prob mod5
          mod.5.det<-tryCatch({
            mod.5.det<-data.frame(ModNames="Observer_Ord.Day_Time",Det=mod.5.detect.B)
          },error=function(cond2){
            cond2=data.frame(ModNames="Observer_Ord.Day_Time",Det=NA)
          })
          
          #get detection prob mod6
          mod.6.det<-tryCatch({
            mod.6.det<-data.frame(ModNames="Observer_Ord.Day_Wind",Det=mod.6.detect.B)
          },error=function(cond2){
            cond2=data.frame(ModNames="Observer_Ord.Day_Wind",Det=NA)
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
          Modnames <- c("Null","Observer","Ord.Day", "Observer_Ord.Day", "Observer_Ord.Day_Time","Observer_Ord.Day_Wind") 
          
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
          aic.summary.df.5<-as.data.frame(aic.summary.df.5)
          mod.keep<-aic.summary.df.5[1,]
          
          #make data.frame of stuff to add to model output to keep track of things
          species.data<-subset(year.data, AOU_Code==speciesName)
          keep.info<-unique(species.data[,c("Year","Unit_Code","AOU_Code")])[1,] 
          
          #add keep.info to model output
          mod.keep.2<-cbind(keep.info, mod.keep)
          
          #get model name
          
          #Let "modSel" function parameter determine whether to use maxDet or AIC to choose model predicted estimates
          mod.name<-tryCatch({
                  switch(modSel,
                           "maxDet"=mod.max.det,
                           "AIC"=as.character(mod.keep.2$Modnames),
                            mod.max.det)
          },error=function(cond2){
            cond2=NA
          })
            
            
          mod.det.df<-subset(mod.det.results.order, ModNames==mod.name)
          mod.det<-mod.det.df$Det
          
            # ifelse(modSel=="maxDet",
            #        #For maximizing detection prob
            #        mod.max.det,
            #      ifelse(modSel=="AIC",
            #               #For AIC
            #              as.character(mod.keep.2$Modnames),
            #             #Default = maximizing detection prob
            #             mod.max.det
            #      ))
          
          tryCatch({
            ifelse(habitatType=="Forest",
                 
          {       
          #get occupancy estimates from top model
          occu.out<-NULL
          ifelse(mod.name=="Null",occu.out<-suppressWarnings(predict(mod1,type="state",appendData=TRUE)),
                 ifelse(mod.name=="Observer", occu.out<-suppressWarnings(predict(mod2,type="state",appendData=TRUE)),
                        ifelse(mod.name=="Ord.Day", occu.out<-suppressWarnings(predict(mod3,type="state",appendData=TRUE)),
                               ifelse(mod.name=="Observer_Ord.Day", occu.out<-suppressWarnings(predict(mod4,type="state",appendData=TRUE)),
                                      ifelse(mod.name=="Observer_Ord.Day_Time", occu.out<-suppressWarnings(predict(mod5,type="state",appendData=TRUE)),
                                             ifelse(mod.name=="Observer_Ord.Day_Wind", occu.out<-suppressWarnings(predict(mod6,type="state",appendData=TRUE)),
                                                    occu.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,Year=yearName,Unit_Code=unitName,GRTS_Order=NA,Visit.1=NA,Visit.2=NA,Ord.Day.1=NA,Ord.Day.2=NA,Time.1=NA,Time.2=NA, Temp.1=NA,Temp.2=NA, Wind.1=NA,Wind.2=NA,Observer.1=NA,Observer.2=NA)
                                             ))))))
          

          #get p(detection) estimates from top model
          det.out<-NULL
          ifelse(mod.name=="Null",det.out<-suppressWarnings(predict(mod1,type="det",appendData=TRUE)),
                 ifelse(mod.name=="Observer", det.out<-suppressWarnings(predict(mod2,type="det",appendData=TRUE)),
                        ifelse(mod.name=="Ord.Day", det.out<-suppressWarnings(predict(mod3,type="det",appendData=TRUE)),
                               ifelse(mod.name=="Observer_Ord.Day", det.out<-suppressWarnings(predict(mod4,type="det",appendData=TRUE)),
                                      ifelse(mod.name=="Observer_Ord.Day_Time", det.out<-suppressWarnings(predict(mod5,type="det",appendData=TRUE)),
                                             ifelse(mod.name=="Observer_Ord.Day_Wind", det.out<-suppressWarnings(predict(mod6,type="det",appendData=TRUE)),
                                                    det.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,Year=yearName,Unit_Code=unitName,GRTS_Order=NA,Visit.1=NA,Visit.2=NA,Ord.Day.1=NA,Ord.Day.2=NA,Time.1=NA,Time.2=NA, Temp.1=NA,Temp.2=NA, Wind.1=NA,Wind.2=NA,Observer.1=NA,Observer.2=NA)
                                             ))))))
          
          
          #get model formula
          formula.out<-NULL
          ifelse(mod.name=="Null",formula.out<-paste(deparse(mod1@call)[1],sep=""),
                 ifelse(mod.name=="Observer", formula.out<-paste(deparse(mod2@call)[1],sep=""),
                        ifelse(mod.name=="Ord.Day", formula.out<-paste(deparse(mod3@call)[1],sep=""),
                               ifelse(mod.name=="Observer_Ord.Day", formula.out<-paste(deparse(mod4@call)[1],sep=""),
                                      ifelse(mod.name=="Observer_Ord.Day_Time", formula.out<-paste(deparse(mod5@call)[1],sep=""),
                                             ifelse(mod.name=="Observer_Ord.Day_Wind", formula.out<-paste(deparse(mod6@call)[1],sep="") ,
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
                        ifelse(mod.name=="Ord.Day", occu.out<-predict(mod3,type="state",appendData=TRUE),
                               ifelse(mod.name=="Observer_Ord.Day", occu.out<-predict(mod4,type="state",appendData=TRUE),
                                      ifelse(mod.name=="Observer_Ord.Day_Time", occu.out<-predict(mod5,type="state",appendData=TRUE),
                                             ifelse(mod.name=="Observer_Ord.Day_Wind", occu.out<-predict(mod6,type="state",appendData=TRUE),
                                                    occu.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,y.3=NA,Year=yearName,Unit_Code=unitName,GRTS_Order=NA,Visit.1=NA,Visit.2=NA,Visit.3=NA,Ord.Day.1=NA,Ord.Day.2=NA,Ord.Day.3=NA,Time.1=NA,Time.2=NA,Time.3=NA, Temp.1=NA,Temp.2=NA,Temp.3=NA, Wind.1=NA,Wind.2=NA,Wind.3=NA,Observer.1=NA,Observer.2=NA,Observer.3=NA)
                                             ))))))
          
          
          #get p(detection) estimates from top model
          det.out<-NULL
          ifelse(mod.name=="Null",det.out<-predict(mod1,type="det",appendData=TRUE),
                 ifelse(mod.name=="Observer", det.out<-predict(mod2,type="det",appendData=TRUE),
                        ifelse(mod.name=="Ord.Day", det.out<-predict(mod3,type="det",appendData=TRUE),
                               ifelse(mod.name=="Observer_Ord.Day", det.out<-predict(mod4,type="det",appendData=TRUE),
                                      ifelse(mod.name=="Observer_Ord.Day_Time", det.out<-predict(mod5,type="det",appendData=TRUE),
                                             ifelse(mod.name=="Observer_Ord.Day_Wind", det.out<-predict(mod6,type="det",appendData=TRUE),
                                                    det.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,y.3=NA,Year=yearName,Unit_Code=unitName,GRTS_Order=NA,Visit.1=NA,Visit.2=NA,Visit.3=NA,Ord.Day.1=NA,Ord.Day.2=NA,Ord.Day.3=NA,Time.1=NA,Time.2=NA,Time.3=NA, Temp.1=NA,Temp.2=NA,Temp.3=NA, Wind.1=NA,Wind.2=NA,Wind.3=NA,Observer.1=NA,Observer.2=NA,Observer.3=NA)
                                             ))))))
          
          
          #get model formula
          formula.out<-NULL
          ifelse(mod.name=="Null",formula.out<-paste(deparse(mod1@call)[1],sep=""),
                 ifelse(mod.name=="Observer", formula.out<-paste(deparse(mod2@call)[1],sep=""),
                        ifelse(mod.name=="Ord.Day", formula.out<-paste(deparse(mod3@call)[1],sep=""),
                               ifelse(mod.name=="Observer_Ord.Day", formula.out<-paste(deparse(mod4@call)[1],sep=""),
                                      ifelse(mod.name=="Observer_Ord.Day_Time", formula.out<-paste(deparse(mod5@call)[1],sep=""),
                                             ifelse(mod.name=="Observer_Ord.Day_Wind", formula.out<-paste(deparse(mod6@call)[1],sep="") ,
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
            data.frame(Habitat=habitatType,Unit_Code=unitName, Year=yearName, AOU_Code=speciesName, Metric=results.comb$Metric, Predicted=results.comb$Predicted, SE=results.comb$SE, lower=results.comb$lower, upper=results.comb$upper, OverallDet=mod.det)
          },error=function(cond2){
            cond2=data.frame(Habitat=habitatType,Unit_Code=unitName, Year=yearName, AOU_Code=speciesName, Metric=NA, Predicted=NA, SE=NA, lower=NA, upper=NA, OverallDet=NA)
          })
            
          #get unique rows
          results.out.unique<-unique(results.out)
          
          #add AICc info back in 
          
          maxDet.AIC<-
            tryCatch({
              subset(aic.summary.df.5, Modnames==mod.name)
            },error=function(cond2){
              cond2=data.frame(Modnames=NA, K=NA, AICc=NA, Delta_AICc=NA, ModelLik=NA, AICcWt=NA, LL=NA, Cum.Wt=NA)
            })
          
          ifelse(modSel=="AIC",
            results.out.unique<-
              tryCatch({
              data.frame(results.out.unique, mod.keep.2[-1:-3])
              }, error=function(cond2){
                cond2=data.frame(results.out.unique,Modnames=NA, K=NA, AICc=NA, Delta_AICc=NA, ModelLik=NA, AICcWt=NA, LL=NA, Cum.Wt=NA)
              }),
            results.out.unique<-tryCatch({
              data.frame(results.out.unique, maxDet.AIC)
            },error=function(cond2){
              cond2=data.frame(results.out.unique,Modnames=NA, K=NA, AICc=NA, Delta_AICc=NA, ModelLik=NA, AICcWt=NA, LL=NA, Cum.Wt=NA)
            })
              )
          
          
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
