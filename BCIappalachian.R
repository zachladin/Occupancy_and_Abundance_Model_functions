#BCI.appalachian function
BCIappalachian<-function(DataIn, unitName, yearName){
  require(dplyr)
  
  #read in Appalachian bird guild assignments
  birdGuilds.appalachian<- read.csv("./Data/BCI_Guild_Assignments/BirdGuildAssignments_Appalachian.csv", header=TRUE)
  
  #truncate birdGuilds
  birdGuilds.appalachian<-birdGuilds.appalachian[,c("AOU_Code","Migratory", "Brood.Numbers", "Primary.Habitat",
                                                    "Exotic", "Nest.Placement", "Foraging.Behavior","Trophic.Level",
                                                    "Nest.Height","Pred_Para_Desc")]

  guilds<-birdGuilds.appalachian
  unit.data<-subset(DataIn,c(Unit_Code==unitName & Year==yearName))
  
  pointList<-unique(unlist(unit.data$Plot_Name))
  
  myData<-list()
  for(i in 1:length(pointList)){
    point.data<-subset(unit.data, Plot_Name==pointList[i])
    sppList<-data.frame(unique(unlist(point.data$AOU_Code)), unique(point.data$Plot_Name))
    colnames(sppList)<-c("AOU_Code","Plot_Name")
    data.merge<-merge(sppList, guilds, by="AOU_Code")
    
    XBCI<-data.frame(Point_Name=data.merge$Plot_Name)
    XBCI<-XBCI %>% rowwise %>% 
      mutate(ForestGeneralist= sum(data.merge$Primary.Habitat=="Forest Generalist"),
             InteriorForest=sum(data.merge$Primary.Habitat=="Interior Forest Obligate"),
             ForestGroundNester=sum(data.merge$Nest.Placement=="Forest Ground Nester"),
             OpenGroundNester=sum(data.merge$Nest.Placement=="Open Ground Nester"),
             ShrubNester=sum(data.merge$Nest.Placement=="Shrub"),
             CanopyNester=sum(data.merge$Nest.Placement=="Canopy"),
             BarkProber=sum(data.merge$Foraging.Behavior=="Bark Prober" & data.merge$Trophic.Level=="Insectivore"),
             GroundGleaner=sum(data.merge$Foraging.Behavior=="Ground Gleaner" & data.merge$Trophic.Level=="Insectivore"),
             CanopyInsectivore=sum(data.merge$Foraging.Behavior=="Upper Canopy Forager" & data.merge$Trophic.Level=="Insectivore"),
             ShrubInsectivore=sum(data.merge$Foraging.Behavior=="Lower Canopy Forager" & data.merge$Trophic.Level=="Insectivore"),
             Omnivore=sum(data.merge$Trophic.Level=="Omnivore"),
             NestPredator=sum(data.merge$Pred_Para_Desc=="Yes"),
             Exotic=sum(data.merge$Exotic=="Exotic"),
             Resident=sum(data.merge$Migratory=="Resident"),
             TemperateMigrant=sum(data.merge$Migratory=="Temperate Migrant"),
             SingleBrooded=sum(data.merge$Brood.Numbers=="Single Brooded"),
             
             TotalHabitat=length(data.merge$Primary.Habitat!=""), #Totals for BCI calc
             TotalNest=length(data.merge$Nest.Placement!=""),
             TotalForage=length(data.merge$Foraging.Behavior!="" & data.merge$Trophic.Level!=""),
             TotalTrophic=length(data.merge$Trophic.Level!=""),
             TotalPredPara=length(data.merge$Pred_Para_Desc!=""),
             TotalExotic=length(data.merge$Exotic!=""),
             TotalMigratory=length(data.merge$Migratory!=""),
             TotalBrood=length(data.merge$Brood.Numbers!=""),
             Total=as.integer(length(sppList$aou.code))
      ) 
    XBCI<-XBCI %>% 
      mutate(Pro_ForestGeneralist=round(ForestGeneralist/TotalHabitat,digits=3),
             Pro_InteriorForest=round(InteriorForest/TotalHabitat,digits=3),
             Pro_ForestGroundNester=round(ForestGroundNester/TotalNest,digits=3),
             Pro_OpenGroundNester=round(OpenGroundNester/TotalNest,digits=3),
             Pro_ShrubNester=round(ShrubNester/TotalNest,digits=3),
             Pro_CanopyNester=round(CanopyNester/TotalNest,digits=3),
             Pro_BarkProber=round(BarkProber/TotalForage,digits=3),
             Pro_GroundGleaner=round(GroundGleaner/TotalForage,digits=3),
             Pro_CanopyInsectivore=round(CanopyInsectivore/TotalForage,digits=3),
             Pro_ShrubInsectivore=round(ShrubInsectivore/TotalForage,digits=3),
             Pro_Omnivore=round(Omnivore/TotalTrophic,digits=3),
             Pro_NestPredator=round(NestPredator/TotalPredPara,digits=3),
             Pro_Exotic=round(Exotic/TotalExotic, digits=3),
             Pro_Resident=round(Resident/TotalMigratory,digits=3),
             Pro_TemperateMigrant=round(TemperateMigrant/TotalMigratory,digits=3),
             Pro_SingleBrooded=round(SingleBrooded/TotalBrood,digits=3)
      )
    
    XBCI<-XBCI %>% 
      mutate(BCI_ForestGeneralist=c(4.5, 2.5)[findInterval(Pro_ForestGeneralist,vec=c(0, 0.281, 1.001))],
             BCI_InteriorForest=c(1, 1.5, 3, 4, 5)[findInterval(Pro_InteriorForest, vec=c(0,0.011,0.081,0.261,0.431,1.001))],
             BCI_ForestGroundNester=c(1, 1.5, 3, 4.5, 5)[findInterval(Pro_ForestGroundNester, vec=c(0,0.001,0.021,0.161,0.241,1.001))],
             BCI_OpenGroundNester=c(1, 2.5, 5)[findInterval(Pro_OpenGroundNester, vec=c(0,0.021,0.111,1.001))],
             BCI_ShrubNester=c(4, 1.5, 1)[findInterval(Pro_ShrubNester, vec=c(0,0.211,0.331,1.001))],
             BCI_CanopyNester=c(1.5, 2, 4.5)[findInterval(Pro_CanopyNester, vec=c(0,0.281,0.321,1.001))],
             BCI_BarkProber=c(1.5, 3, 4, 5)[findInterval(Pro_BarkProber, vec=c(0,0.061,0.111,0.171,1.001))],
             BCI_GroundGleaner=c(1.5, 2, 4.5, 5)[findInterval(Pro_GroundGleaner, vec=c(0,0.051,0.071,0.141,1.001))],
             BCI_CanopyInsectivore=c(1.5, 2, 3, 4.5, 5)[findInterval(Pro_CanopyInsectivore, vec=c(0,0.031,0.051,0.121,0.201,1.001))],
             BCI_ShrubInsectivore=c(1.5, 2.5, 5)[findInterval(Pro_ShrubInsectivore, vec=c(0,0.141,0.231,1.001))],
             BCI_Omnivore=c(5, 4, 3, 1, 2)[findInterval(Pro_Omnivore, vec=c(0,0.291,0.411,0.481,0.581,1.001))],
             BCI_NestPredator=c(5, 3.5, 2, 1)[findInterval(Pro_NestPredator, vec=c(0,0.101,0.151,0.181,1.001))],
             BCI_Exotic=c(5, 4.5, 3, 2, 1)[findInterval(Pro_Exotic, vec=c(0,0.001,0.021,0.051,0.111,1.001))],
             BCI_Resident=c(5, 3.5, 2, 1)[findInterval(Pro_Resident, vec=c(0,0.261,0.391,0.571,1.001))],
             BCI_TemperateMigrant=c(4, 2, 1)[findInterval(Pro_TemperateMigrant, vec=c(0,0.211,0.301,1.001))],
             BCI_SingleBrooded=c(1.5, 2, 3, 4, 5)[findInterval(Pro_SingleBrooded, vec=c(0,0.411,0.451,0.611,0.731,1.001))]
      )
    XBCI<-XBCI %>% 
      mutate(BCI_Structural=BCI_ForestGeneralist + BCI_InteriorForest + BCI_ForestGroundNester + BCI_OpenGroundNester + BCI_ShrubNester + 
               BCI_CanopyNester,
             BCI_Funcitonal=BCI_BarkProber + BCI_GroundGleaner + BCI_CanopyInsectivore + BCI_ShrubInsectivore + BCI_Omnivore,
             BCI_Compostional=BCI_NestPredator + BCI_Exotic + BCI_Resident + BCI_TemperateMigrant + BCI_SingleBrooded,
             BCI = BCI_Structural + BCI_Funcitonal + BCI_Compostional,
             BCI_Category=c("Low Integrity", "Medium Integrity","High Integrity","Highest Integrity")[findInterval(BCI,
                                                                                                                   vec=c(0,40.1,52.1,60.1,77.1))]
      )
    BCI.data<-data.frame(unitName,yearName,unique(XBCI))
    myData<-rbind(myData,BCI.data)
    BCI.Appalachian<-as.data.frame(myData)
  }
  return(BCI.Appalachian)
}
