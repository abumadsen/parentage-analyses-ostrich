#31/08/21
#Parentage 
#Masterbayes

pacman::p_load("reshape","tidyverse")

#Get plate summary stats
check.plate = function(platedf,plateid){
  
  #Markers with no calls due to invalid
  NodataMar <- paste(names(which(colSums(platedf == "Invalid") == "96")), collapse =";")
  NodataMar <- data.frame("InvalidMarkers" = NodataMar)
  
  #Individuals with few calls due to No Call
  LowData60Ind <- paste(platedf$id[rowSums(platedf == "No Call") >= 60], collapse =";")
  LowData60Ind <- data.frame("Individuals60nocall" = LowData60Ind)
  
  LowData80Ind <- paste(platedf$id[rowSums(platedf == "No Call") >= 80], collapse =";")
  LowData80Ind <- data.frame("Individuals80nocall" = LowData80Ind)
  
  LowData90Ind <- paste(platedf$id[rowSums(platedf == "No Call") >= 90], collapse =";")
  LowData90Ind <- data.frame("Individuals90nocall" = LowData90Ind)
  
  #Negative control
  NC <- factor(unlist(platedf[platedf$id %in% "NC",-1]), levels = unique(c("Invalid","NTC", "No Call",unlist(platedf[platedf$id %in% "NC",-1])))) #Add factor levels in case they are missing so output is standardized
  NC <- data.frame(table(NC))
  NC <- NC[grep(":", NC$NC, invert = T),] #Remove the X:X calls
  NC <- data.frame(t(NC))
  colnames(NC) <- paste("NC",unlist(NC[1,]), sep = ":")
  
  #Postive control
  PC <- factor(unlist(platedf[platedf$id %in% c("PC"),-1]), levels = unique(c("Invalid", "No Call",unlist(platedf[platedf$id %in% c("PC"),-1])))) #Add factor levels in case they are missing so output is standardized
  PC <- data.frame(table(PC))  
  PC <- PC[grep(":", PC$PC, invert = T),] #Remove the X:X calls
  PC <- data.frame(t(PC))
  colnames(PC) <- paste("PC",unlist(PC[1,]), sep = ":")
  
  NPC <- cbind(NC,PC)[2,]
  
  #Overall missing calls
  All <- data.frame("Plateid" = plateid, "NoCall" = sum(platedf == "No Call"), "Invalid" = sum(platedf == "Invalid"))
  
  #Reporting
  return(cbind(All,NodataMar,LowData60Ind,LowData80Ind,LowData90Ind,NPC))
}


#Get individual summary stats
check.individual = function(platedf,plateid){
  
  #Count nr of No Call or Invalid
  platedf$NoCall <- rowSums(platedf == "No Call")
  platedf$Invalid <- rowSums(platedf == "Invalid")

  #Reporting
  return(data.frame("Plateid" = plateid,platedf[,c("id","NoCall","Invalid")]))
}

#------------------

#Peprare plate for MasterBayes
prep.plate = function(platedf,plateid){
  
  #Remove/Fix
  platedf <- platedf[!(platedf$id == "NC" | platedf$id == "PC" | platedf$id == "Blank"),]
  platedf[platedf == "No Call" | platedf == "Invalid"] <- "NA:NA"
  
  #Split each genotype into two columns
  platedf2 = data.frame("id" = platedf$id)
  for(col in 2:ncol(platedf)){
    platedftemp <- colsplit(platedf[,col], split = ":", names = paste(colnames(platedf)[col], c("a","b"), sep = "")) #Create two alleles
    platedf2 <- cbind(platedf2,platedftemp)
  }
  
  #Make <NA> into NA = important for masterbayes
  #I also appears that genotypes need to be numeric, for the NA to be seen as a true NA, and not cause genotypes like NA/NA, which should just be a NA
  platedf2 <- data.frame(apply(platedf2, 2, as.character))
  platedf2[platedf2 == "A" & !is.na(platedf2)] <- "1"
  platedf2[platedf2 == "T" & !is.na(platedf2)] <- "2"
  platedf2[platedf2 == "C" & !is.na(platedf2)] <- "3"
  platedf2[platedf2 == "G" & !is.na(platedf2)] <- "4"
  platedf2[,colnames(platedf2) != "id"] <- data.frame(apply(platedf2[,colnames(platedf2) != "id"], 2, as.numeric))
  
  return(cbind(plateid,platedf2))
}


#------------------

#Peprare plate for sex determination
prep.plate.sex = function(platedf,plateid){
  
  #Remove/Fix
  platedfsex <- platedf[!(platedf$id == "NC" | platedf$id == "PC"),]
  #platedfsex[platedfsex == "No Call" | platedfsex == "Invalid"] <- NA
  platedfsex[platedfsex == "Invalid"] <- NA
  
  #Subset Z and W
  platedfsex <- platedfsex[,grepl("Wnonpoly_",colnames(platedfsex)) | grepl("Z_",colnames(platedfsex)) | colnames(platedfsex) == "id"]
  
  return(cbind(plateid,platedfsex))
}

