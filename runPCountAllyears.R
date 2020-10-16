runPCountAllyearsNetworkWide<-function(DataIn, modSel, networkName,nGOF){
  
  require(unmarked)
  require(plyr)
  require(reshape2)
  require(AICcmodavg)
  
  setMethod("coerce",c("NULL","data.frame"), function(from, to, strict=TRUE) as.data.frame(from))
  
  data<-DataIn
  
  #convert Unit_Code to character
  data$Unit_Code<-as.character(data$Unit_Code)
  
  #change Survey_Type to Habitat
  #names(data)[names(data)=="Survey_Type"]<-"Habitat"
  
  #create habitatList
  habitatList<-unique(as.character(data$Survey_Type))
  
  habitat.save<-list()
  for(i in 1:length(habitatList)){
    
    #set habitatType
    habitatType<-as.character(habitatList[i])
    
    #subset data by habitat
    new.data<-subset(data, Survey_Type==habitatType)
    
      unitName<-as.character(networkName)
      print(unitName)
      
        #get speciesList
        speciesList<-sort(as.character(unique(new.data$AOU_Code)))

        species.save<-list()
        for(r in 1:length(speciesList)){
          
        speciesName<-speciesList[r]
        print(speciesName)
          
          #create unmarkedFramePCount_umf
          #make umf
          new.umf<-tryCatch({
            suppressWarnings(makePCount_umf(DataIn=new.data, SpeciesIn=speciesName))
          },error=function(cond2){
            cond2=NULL
          })
          
          #Fit models
          #set all models to equal NULL to clear previous results
          mod1=mod2=mod3=mod4=mod5=mod6=mod7=mod8=NULL
          
          #use convert.units=0.01 (0.01*0.01) to convert from m^2 to ha
          
          #model 1 - Null - P
          mod1<-tryCatch({pcount(~1 ~Year, data=new.umf, mixture="P")
          }, error= function(out){
            out=NULL
          })
          
          #mod1.gof
          mod1.gof<- try(Nmix.gof.test(mod1, nsim = 20, plot.hist = FALSE,report = NULL, parallel = FALSE))
          
          #get c-hat (measure of over or under dispersion of model)
          mod1.c.hat<-tryCatch({
            mod1.gof$c.hat.est
          },error=function(cond2){
            cond2=NA
          })
          
          
          #model 2 - Null - NB (det)
          mod2<-tryCatch({pcount(~1 ~Year, data=new.umf, mixture="NB")
          }, error=function(out){
            out=NULL
          })
          
          #mod2.gof
          mod2.gof<- try(Nmix.gof.test(mod2, nsim = 20, plot.hist = FALSE,report = NULL, parallel = FALSE))
          
          #get c-hat (measure of over or under dispersion of model)
          mod2.c.hat<-tryCatch({
            mod2.gof$c.hat.est
          },error=function(cond2){
            cond2=NA
          })
          
          #model 3 - Observer - P (det)
          mod3<-tryCatch({pcount(~Observer ~Year, data=new.umf, mixture="P")
          }, error=function(out){
            out=NULL
          })
          
          #mod3.gof
          mod3.gof<- try(Nmix.gof.test(mod3, nsim = 20, plot.hist = FALSE,report = NULL, parallel = FALSE))
          
          #get c-hat (measure of over or under dispersion of model)
          mod3.c.hat<-tryCatch({
            mod3.gof$c.hat.est
          },error=function(cond2){
            cond2=NA
          })
          
          
          #model 4 - Observer - NB (det) 
          mod4<-tryCatch({pcount(~Observer ~ Year, data=new.umf, mixture="NB")
          }, error=function(out){
            out=NULL
          })
          
          #mod4.gof
          mod4.gof<- try(Nmix.gof.test(mod4, nsim = 20, plot.hist = FALSE,report = NULL, parallel = FALSE))
          
          #get c-hat (measure of over or under dispersion of model)
          mod4.c.hat<-tryCatch({
            mod4.gof$c.hat.est
          },error=function(cond2){
            cond2=NA
          })
          
          
          #model 5 - Observer+Ord.Day - P
          mod5<-tryCatch({pcount(~Observer+Ord.Day ~ Year, data=new.umf, mixture="P")
          }, error=function(out){
            out=NULL
          })
          
          #mod5.gof
          mod5.gof<- try(Nmix.gof.test(mod5, nsim = 20, plot.hist = FALSE,report = NULL, parallel = FALSE))
          
          #get c-hat (measure of over or under dispersion of model)
          mod5.c.hat<-tryCatch({
            mod5.gof$c.hat.est
          },error=function(cond2){
            cond2=NA
          })
          

          #model 6 - Observer+Ord.Day - NB
          mod6<-tryCatch({pcount(~Observer+Ord.Day ~ Year, data=new.umf, mixture="NB")
          }, error=function(out){
            out=NULL
          })
          
          #mod6.gof
          mod6.gof<- try(Nmix.gof.test(mod6, nsim = 20, plot.hist = FALSE,report = NULL, parallel = FALSE))
          
          #get c-hat (measure of over or under dispersion of model)
          mod6.c.hat<-tryCatch({
            mod6.gof$c.hat.est
          },error=function(cond2){
            cond2=NA
          })
          
          # #model 7 - Null - P (No SE)
          # mod7<-tryCatch({suppressWarnings(pcount(~1 ~ Year, data=new.umf, mixture="P", se=FALSE))
          # }, error=function(out){
          #   out=NULL
          # })
          # 
          # #mod7.gof
          # mod7.gof<- try(Nmix.gof.test(mod7, nsim = 20, plot.hist = FALSE,report = NULL, parallel = FALSE))
          # 
          # #get c-hat (measure of over or under dispersion of model)
          # mod7.c.hat<-tryCatch({
          #   mod7.gof$c.hat.est
          # },error=function(cond2){
          #   cond2=NA
          # })
          # 
          # 
          # #model 8 - Null - NB (No SE)
          # mod8<-tryCatch({suppressWarnings(pcount(~1 ~ Year, data=new.umf, mixture="NB", se=FALSE))
          # }, error=function(out){
          #   out=NULL
          # })
          # 
          # #mod1.gof
          # mod8.gof<- try(Nmix.gof.test(mod8, nsim = 20, plot.hist = FALSE,report = NULL, parallel = FALSE))
          # 
          # #get c-hat (measure of over or under dispersion of model)
          # mod8.c.hat<-tryCatch({
          #   mod8.gof$c.hat.est
          # },error=function(cond2){
          #   cond2=NA
          # })
          
          
          ####################################################################################################################
          #get detection probablilities from each model and sort to maximize detection
          
          # try({
          #   
          #   ifelse(!is.null(mod1),
          #          {
          #            mod.1.detect.A<-NULL
          #            mod.1.detect.A<-unique(suppressWarnings(predict(mod1, type="det")))
          #            mod.1.detect.B<-mod.1.detect.A$Predicted[1]
          #          },
          #          mod.1.detect.B<-NA)
          #   
          #   ifelse(!is.null(mod2),
          #          {
          #            mod.2.detect.A<-NULL
          #            mod.2.detect.A<-unique(suppressWarnings(predict(mod2, type="det")))
          #            mod.2.detect.B<-mod.2.detect.A$Predicted[1]
          #          },
          #          mod.2.detect.B<-NA)
          #   
          #   ifelse(!is.null(mod3),
          #          {
          #            mod.3.detect.A<-NULL
          #            mod.3.detect.A<-unique(suppressWarnings(predict(mod3, type="det")))
          #            mod.3.detect.A$One_minus_pred<-1-mod.3.detect.A$Predicted
          #            mod.3.detect.B<-1-prod(mod.3.detect.A$One_minus_pred, na.rm=TRUE)
          #          },
          #          mod.3.detect.B<-NA)
          #   
          #   ifelse(!is.null(mod4),
          #          {
          #            mod.4.detect.A<-NULL
          #            mod.4.detect.A<-unique(suppressWarnings(predict(mod4, type="det")))
          #            mod.4.detect.A$One_minus_pred<-1-mod.4.detect.A$Predicted
          #            mod.4.detect.B<-1-prod(mod.4.detect.A$One_minus_pred, na.rm=TRUE)
          #          },
          #          mod.4.detect.B<-NA)
          #   
          #   ifelse(!is.null(mod5),
          #          {
          #            mod.5.detect.A<-NULL
          #            mod.5.detect.A<-unique(suppressWarnings(predict(mod5, type="det")))
          #            mod.5.detect.A$One_minus_pred<-1-mod.5.detect.A$Predicted
          #            mod.5.detect.B<-1-prod(mod.5.detect.A$One_minus_pred, na.rm=TRUE)
          #          },
          #          mod.5.detect.B<-NA)
          #   
          #   ifelse(!is.null(mod6),
          #          {
          #            mod.6.detect.A<-NULL
          #            mod.6.detect.A<-unique(suppressWarnings(predict(mod6, type="det")))
          #            mod.6.detect.A$One_minus_pred<-1-mod.6.detect.A$Predicted
          #            mod.6.detect.B<-1-prod(mod.6.detect.A$One_minus_pred, na.rm=TRUE)
          #          },
          #          mod.6.detect.B<-NA)
          # },silent=TRUE)
          # 
          # #get detection prob mod1
          # mod.1.det<-tryCatch({
          #   mod.1.det<-data.frame(ModNames="Null_P",Det=mod.1.detect.B)
          # },error=function(cond2){
          #   cond2=data.frame(ModNames="Null_P",Det=NA)
          # })
          # 
          # #get detection prob mod2
          # mod.2.det<-tryCatch({
          #   mod.2.det<-data.frame(ModNames="Null_NB",Det=mod.2.detect.B)
          # },error=function(cond2){
          #   cond2=data.frame(ModNames="Null_NB",Det=NA)
          # })
          # 
          # #get detection prob mod3
          # mod.3.det<-tryCatch({
          #   mod.3.det<-data.frame(ModNames="Observer_P",Det=mod.3.detect.B)
          # },error=function(cond2){
          #   cond2=data.frame(ModNames="Observer_P",Det=NA)
          # })
          # 
          # #get detection prob mod4
          # mod.4.det<-tryCatch({
          #   mod.4.det<-data.frame(ModNames="Observer_NB",Det=mod.4.detect.B)
          # },error=function(cond2){
          #   cond2=data.frame(ModNames="Observer_NB",Det=NA)
          # })
          # 
          # #get detection prob mod5
          # mod.5.det<-tryCatch({
          #   mod.5.det<-data.frame(ModNames="Observer.Day_P",Det=mod.5.detect.B)
          # },error=function(cond2){
          #   cond2=data.frame(ModNames="Observer.Day_P",Det=NA)
          # })
          # 
          # #get detection prob mod6
          # mod.6.det<-tryCatch({
          #   mod.6.det<-data.frame(ModNames="Observer.Day_NB",Det=mod.6.detect.B)
          # },error=function(cond2){
          #   cond2=data.frame(ModNames="Observer.Day_NB",Det=NA)
          # })
          # 
          # #compile detection results
          # mod.det.results<-rbind(mod.1.det, mod.2.det, mod.3.det, mod.4.det, mod.5.det, mod.6.det)
          # 
          # #remove NAs
          # mod.det.results.noNA<-na.omit(mod.det.results)
          # 
          # #order by det high to low
          # mod.det.results.order<-mod.det.results.noNA[order(mod.det.results.noNA$Det, decreasing=TRUE),]
          # 
          # #mod to keep
          # mod.max.det<-as.character(mod.det.results.order$ModNames[1])
          
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
          
          
          #if still no model works then run AIC on mod7 and mod8
          noSE.Modnames<-c("Null_P_NoSE","Null_NB_NoSE")
          
          if(is.null(aic.summary.df.4)){
            aic.summary.df.5<-tryCatch({aictab(cand.set = list(mod7,mod8), modnames = unlist(noSE.Modnames))
            },error = function(out){out=tryCatch({aictab(cand.set = list(mod7), modnames = unlist(noSE.Modnames[-2]))
            },error = function(out){out=tryCatch({aictab(cand.set = list(mod8), modnames = unlist(noSE.Modnames[-1]))
            },error = function(out){ 
              out=NULL
            })})})
          } else {
            aic.summary.df.5<-NULL
          }
          
          
          ifelse(!is.null(aic.summary.df.1), aic.summary.df.6<-aic.summary.df.1,
                 ifelse(!is.null(aic.summary.df.2), aic.summary.df.6<-aic.summary.df.2,
                        ifelse(!is.null(aic.summary.df.3), aic.summary.df.6<-aic.summary.df.3,
                               ifelse(!is.null(aic.summary.df.4), aic.summary.df.6<-aic.summary.df.4,
                                      ifelse(!is.null(aic.summary.df.5), aic.summary.df.6<-aic.summary.df.5,
                                             aic.summary.df.6<-data.frame("Modnames"=NA, "K"=NA, "AICc"=NA, "Delta_AICc"=NA,"ModelLik"=NA, 
                                                                          "AICcWt"=NA, "LL"=NA,"Cum.Wt"=NA))))))
          
          
          #select top performing model
          #take model with lowest AIC score
          mod.keep<-aic.summary.df.6[1,]
          
          #make data.frame of stuff to add to model output to keep track of things
          species.data<-subset(new.data, AOU_Code==speciesName)
          keep.info<-unique(species.data[,c("Year","Unit_Code","AOU_Code")])[1,]
          
          #add keep.info to model output
          mod.keep.2<-cbind(keep.info, mod.keep)
          
          #Let "modSel" function parameter determine whether to use maxDet or AIC to choose model predicted estimates
          mod.name<-switch(modSel,
                           "maxDet"=mod.max.det,
                           "AIC"=unique(as.character(mod.keep.2$Modnames)),
                           mod.max.det)
          

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
                                                               ifelse(mod.name=="Null_P_NoSE", abun.out<-predict(mod7,type="state",appendData=TRUE),
                                                                      ifelse(mod.name=="Null_NB_NoSE", abun.out<-predict(mod8,type="state",appendData=TRUE),
                                                                             abun.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,Year=NA,Unit_Code=unitName,GRTS_Order=NA,Visit.1=NA,Visit.2=NA,Ord.Day.1=NA,Ord.Day.2=NA,Time.1=NA,Time.2=NA, Temp.1=NA,Temp.2=NA, Wind.1=NA,Wind.2=NA,Observer.1=NA,Observer.2=NA)
                                                                      ))))))))
                     
                     
                     #get p(detection) estimates from top model
                     det.out<-NULL
                     ifelse(mod.name=="Null_P",det.out<-predict(mod1,type="det",appendData=TRUE),
                            ifelse(mod.name=="Null_NB", det.out<-predict(mod2,type="det",appendData=TRUE),
                                   ifelse(mod.name=="Observer_P", det.out<-predict(mod3,type="det",appendData=TRUE),
                                          ifelse(mod.name=="Observer_NB", det.out<-predict(mod4,type="det",appendData=TRUE),
                                                 ifelse(mod.name=="Observer.Day_P", det.out<-predict(mod5,type="det",appendData=TRUE),
                                                        ifelse(mod.name=="Observer.Day_NB", det.out<-predict(mod6,type="det",appendData=TRUE),
                                                               ifelse(mod.name=="Null_P_NoSE", det.out<-predict(mod7,type="det",appendData=TRUE),
                                                                      ifelse(mod.name=="Null_NB_NoSE", det.out<-predict(mod8,type="det",appendData=TRUE),
                                                                             det.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,Year=NA,Unit_Code=unitName,GRTS_Order=NA,Visit.1=NA,Visit.2=NA,Ord.Day.1=NA,Ord.Day.2=NA,Time.1=NA,Time.2=NA, Temp.1=NA,Temp.2=NA, Wind.1=NA,Wind.2=NA,Observer.1=NA,Observer.2=NA)
                                                                      ))))))))
                     
                     
                     #get model formula
                     formula.out<-NULL
                     ifelse(mod.name=="Null_P",formula.out<-paste(deparse(mod1@call)[1],sep=""),
                            ifelse(mod.name=="Null_NB", formula.out<-paste(deparse(mod2@call)[1],sep=""),
                                   ifelse(mod.name=="Observer_P", formula.out<-paste(deparse(mod3@call)[1],sep=""),
                                          ifelse(mod.name=="Observer_NB", formula.out<-paste(deparse(mod4@call)[1],sep=""),
                                                 ifelse(mod.name=="Observer.Day_P", formula.out<-paste(deparse(mod5@call)[1],sep=""),
                                                        ifelse(mod.name=="Observer.Day_NB", formula.out<-paste(deparse(mod6@call)[1],sep="") ,
                                                               ifelse(mod.name=="Null_P_NoSE", formula.out<-paste(deparse(mod7@call)[1],sep="") ,
                                                                      ifelse(mod.name=="Null_NB_NoSE", formula.out<-paste(deparse(mod8@call)[1],sep="") ,
                                                                             formula.out<-NA
                                                                      ))))))))
                     
                     
                     
                     
                     
                     
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
                                                               ifelse(mod.name=="Null_P_NoSE", abun.out<-predict(mod7,type="state",appendData=TRUE),
                                                                      ifelse(mod.name=="Null_NB_NoSE", abun.out<-predict(mod8,type="state",appendData=TRUE),
                                                                             abun.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,y.3=NA,Year=NA,Unit_Code=unitName,GRTS_Order=NA,Visit.1=NA,Visit.2=NA,Visit.3=NA,Ord.Day.1=NA,Ord.Day.2=NA,Ord.Day.3=NA,Time.1=NA,Time.2=NA,Time.3=NA, Temp.1=NA,Temp.2=NA,Temp.3=NA, Wind.1=NA,Wind.2=NA,Wind.3=NA,Observer.1=NA,Observer.2=NA,Observer.3=NA)
                                                                      ))))))))
                     
                     
                     #get p(detection) estimates from top model
                     det.out<-NULL
                     ifelse(mod.name=="Null_P",det.out<-predict(mod1,type="det",appendData=TRUE),
                            ifelse(mod.name=="Null_NB", det.out<-predict(mod2,type="det",appendData=TRUE),
                                   ifelse(mod.name=="Observer_P", det.out<-predict(mod3,type="det",appendData=TRUE),
                                          ifelse(mod.name=="Observer_NB", det.out<-predict(mod4,type="det",appendData=TRUE),
                                                 ifelse(mod.name=="Observer.Day_P", det.out<-predict(mod5,type="det",appendData=TRUE),
                                                        ifelse(mod.name=="Observer.Day_NB", det.out<-predict(mod6,type="det",appendData=TRUE),
                                                               ifelse(mod.name=="Null_P_NoSE", det.out<-predict(mod7,type="det",appendData=TRUE),
                                                                      ifelse(mod.name=="Null_NB_NoSE", det.out<-predict(mod8,type="det",appendData=TRUE),
                                                                             det.out<-data.frame(Predicted=NA, SE=NA, lower=NA,upper=NA,y.1=NA,y.2=NA,y.3=NA,Year=NA,Unit_Code=unitName,GRTS_Order=NA,Visit.1=NA,Visit.2=NA,Visit.3=NA,Ord.Day.1=NA,Ord.Day.2=NA,Ord.Day.3=NA,Time.1=NA,Time.2=NA,Time.3=NA, Temp.1=NA,Temp.2=NA,Temp.3=NA, Wind.1=NA,Wind.2=NA,Wind.3=NA,Observer.1=NA,Observer.2=NA,Observer.3=NA)
                                                                      ))))))))
                     
                     
                     #get model formula
                     formula.out<-NULL
                     ifelse(mod.name=="Null_P",formula.out<-paste(deparse(mod1@call)[1],sep=""),
                            ifelse(mod.name=="Null_NB", formula.out<-paste(deparse(mod2@call)[1],sep=""),
                                   ifelse(mod.name=="Observer_P", formula.out<-paste(deparse(mod3@call)[1],sep=""),
                                          ifelse(mod.name=="Observer_NB", formula.out<-paste(deparse(mod4@call)[1],sep=""),
                                                 ifelse(mod.name=="Observer.Day_P", formula.out<-paste(deparse(mod5@call)[1],sep=""),
                                                        ifelse(mod.name=="Observer.Day_NB", formula.out<-paste(deparse(mod6@call)[1],sep="") ,
                                                               ifelse(mod.name=="Null_P_NoSE", formula.out<-paste(deparse(mod7@call)[1],sep="") ,
                                                                      ifelse(mod.name=="Null_NB_NoSE", formula.out<-paste(deparse(mod8@call)[1],sep="") ,
                                                                             formula.out<-NA
                                                                      ))))))))
                     
                     
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
            data.frame(Habitat=habitatType,Unit_Code=unitName, Year=results.comb$Year, AOU_Code=speciesName, Metric=results.comb$Metric, Predicted=results.comb$Predicted, SE=results.comb$SE, lower=results.comb$lower, upper=results.comb$upper)
          },error=function(cond2){
            cond2=data.frame(Habitat=habitatType,Unit_Code=unitName, Year=NA, AOU_Code=speciesName, Metric=NA, Predicted=NA, SE=NA, lower=NA, upper=NA)
          })
          
          #get unique rows
          results.out.unique<-unique(results.out)
          
          
          #add AICc info back in 
          #results.out.unique<-data.frame(results.out.unique, unique(mod.keep.2[-1:-3]))
          
          results.out.unique.2<-NULL
          
          if(modSel=="AIC"){
                 results.out.unique.2<-
                   tryCatch({
                     data.frame(results.out.unique, unique(mod.keep.2[-1:-3]))
                   }, error=function(cond2){
                     cond2=data.frame(results.out.unique,Modnames=NA, K=NA, AICc=NA, Delta_AICc=NA, ModelLik=NA, AICcWt=NA, LL=NA, Cum.Wt=NA)
                   })}else {
                 results.out.unique.2<-tryCatch({
                   data.frame(results.out.unique, maxDet.AIC)
                 },error=function(cond2){
                   cond2=data.frame(results.out.unique,Modnames=NA, K=NA, AICc=NA, Delta_AICc=NA, ModelLik=NA, AICcWt=NA, LL=NA, Cum.Wt=NA)
                 })
                   }
          
          #add goodness of fit results
          chat.out<-NULL
          ifelse(mod.name=="Null_P",chat.out<-mod1.c.hat,
                 ifelse(mod.name=="Null_NB", chat.out<-mod2.c.hat,
                        ifelse(mod.name=="Observer_P", chat.out<-mod3.c.hat,
                               ifelse(mod.name=="Observer_NB", chat.out<-mod4.c.hat,
                                      ifelse(mod.name=="Observer.Day_P", chat.out<-mod5.c.hat,
                                             ifelse(mod.name=="Observer.Day_NB", chat.out<-mod6.c.hat ,
                                                    ifelse(mod.name=="Null_P_NoSE", chat.out<-mod7.c.hat ,
                                                           ifelse(mod.name=="Null_NB_NoSE", chat.out<-mod8.c.hat ,
                                                                  chat.out<-NA
                                                           ))))))))
          
          
          
          
          
          
          #add chat.out to data.frame
          results.out.unique.2$GOF_chat<-chat.out
          
        
          #get raw counts
          max.count<-sapply(new.umf@y, max)
          
          #make data.frame of max count per plot per year
          data.raw<-data.frame(Year=new.umf@siteCovs$Year,Plot_Name=new.umf@siteCovs$Plot_Name, MaxCount=max.count)
          #average max counts per year
          data.agg<-aggregate(MaxCount~Year+Plot_Name, data=data.raw, FUN="max")
          #summary function - mean and SE
          data.summary<-summaryFunction(dataIn=data.agg, response="MaxCount",factor="Year")
          data.summary<-data.summary[,c("Year","n","mean","SE")]
          names(data.summary)<-c("Year","n","Naive_mean","Naive_SE")
          
          #add column "Metric"
          data.summary$Metric<-"Abundance"
          
          #fix n
          samplesYear<-data.summary[,c("Year","n")]
          results.out.unique.2<-merge(results.out.unique.2, samplesYear, by="Year",all.x=TRUE)
          
          #add naive occupancy estimates
          results.out.unique.2<-merge(results.out.unique.2, data.summary[,c("Year","Naive_mean","Naive_SE","Metric")], by=c("Year","Metric"),all.x=TRUE)
          
          #fill park.output list with each iteration through park loop
          print(paste("Saving results for",habitatType, unitName, speciesName,sep=" "))
          species.save<-rbind(species.save, results.out.unique.2)

    }
    
    #fill year.output list with each iteration throught the year loop
    print(paste("Saving results for",habitatType,sep=" "))
    habitat.save<-rbind(habitat.save, species.save)
    
  }
  
  
  return(habitat.save)
}