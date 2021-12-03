#31/08/21
#Parentage 
#Masterbayes
#Analyses

# 1. Parentage analyses with MasterBayes
# 2. Sex determination
# 3. Export

#Masterbayes has to be installed manually after download. Ensure file is at workdir.
#install.packages('MasterBayes_2.57.tar', lib=".", repos = NULL, type = "source")
pacman::p_load("MasterBayes","reshape","tidyverse")

#--- Notes:

#I tried running MasterBayes, letting it estimate genotypes/allele frequencies and error rates from the data
#on two plates, so leaving out the "sP=sP" and "sP<-startPed(estG=FALSE, E1=0.005, E2=0.005, A=extractA(plates.prepped.auto))".
#This however resulted in some very unconservative parentage calls, with the program being confident 
#even if almost no genotypes called, and if I mannipulated the data such that 
#If we give adult amles no gennotypes then they are also assiged as sires
#If I give two males identical genotypes then one of them is assigned (I think in alphabetical order)
#When I inspect the posterior of the potential parent pairs for offspring (e.g. model1$P$c15.344)
#All "hits" were assigned to one pair, even if I had just made two potential fathers completely identical genotype wise.
#It appears that the chains got stuck at one potential, maybe because of lack of enough data to estimate all
#these parameters.

#Estimated pedigree with the above approach:
# [,1]       [,2]      [,3]     
# [1,] "c15.338"  "a13.83"  "a13.44" 
# [2,] "c15.351"  "a13.337" "a10.599"
# [3,] "c15.64"   "a12.195" NA       
# [4,] "c15.339"  "a13.337" "a10.599"
# [5,] "c15.354"  "a10.861" "a13.212"
# [6,] "c15.367"  "a10.861" "a13.212"
# [7,] "c15.344"  "a13.358" "a10.599"
# [8,] "c15.349"  "a13.337" "a10.599"
# [9,] "c15.711"  "a13.337" "a10.599"
# [10,] "c15.699"  "a13.358" "a10.599"
# [11,] "c15.704"  "a13.358" "a10.599"

#This anti-conservativenness is fixed if I fix error-rates and the allele frequencies as also done in most
#of the tutorial:
#https://cran.r-project.org/web/packages/MasterBayes/vignettes/Tutorial.pdf
#now inspection of the potential parent pairs for offspring (e.g. model1$P$c15.344) shows how different pairs 
#get different nnumber of hits

#Later when I get more data I can revisit this, and try to let MasterBayes estimate the error rates
#From the tutorial it also appears that this requires a lot of iterations!
#But I should be cautius, and check the above problems.

#There is also the possibility to model the presence of unsampled males and females
#MyPdP<-PdataPed(formula=list(res1), data=S, USsire = T, USdam =T )
#If we think this is a possibility

###########################################
# Unsampled individuals
###########################################

# There is actually four individuals without a blood sample, two females and two males:
#   1110115: Female. There is no blood sample for this individual, but she was removed early in the season due to an injury, so it might be okay not to have her in the data set.
# 1120173: Male. We should have a sample for this individual, but the sample is not in its box.
# 1120078: Male. We should have a sample for this individual, but the sample is not in its box.
# 1150082: Female. There is no blood sample for this individual.

#Data on these individuals are available in
# Missing_adults_2015.2018.csv

###########################################
# 0. General data prep
###########################################

#Load phenotypic data
InSamples <- "geno_samples_2015.2018"
Sraw <- read.table(paste("input/", InSamples, ".csv", sep = ""), sep = ",", header = T)

#Load functions for genotypic data
source("1 dataprep functions.R")

MyPlates <- as.character(seq(1,18))
plate.info <- as.numeric()
individual.info <- as.numeric()
plates.prepped <- as.numeric()
plates.sex <- as.numeric()

for(i in MyPlates){
  P <- read.table(paste("input/Call map - Plate ", i, ".csv", sep = ""), skip = 115, sep = ",", header = T)
  P <- dplyr::rename(P, c("id" = "X.1"))
  P$X <- NULL
  
  plate.info <- rbind(plate.info,check.plate(platedf = P, plateid = i))
  individual.info <- rbind(individual.info,check.individual(platedf = P, plateid = i))
  plates.prepped <- rbind(plates.prepped,prep.plate(platedf = P, plateid = i))
  plates.sex <- rbind(plates.sex,prep.plate.sex(platedf = P, plateid = i))
}

#Remove markers with no data in any plates as this upsets MasterBayes
Nodata <- names(which(colSums(is.na(plates.prepped)) == nrow(plates.prepped))) #No markers with no data
plates.prepped <- plates.prepped[,!colnames(plates.prepped) %in% Nodata]

##----------- EXPORT info files---------##
write.table(plate.info, "output/plate.info.csv", sep = ",", row.names = F)
write.table(individual.info, "output/individual.info.csv", sep = ",", row.names = F)

###########################################
# 0.5. Project specific data prep
###########################################

#Remove individuals a07.345 (OBS! Also written a7.345) and a11.421 since they are not part of the experimental population
plates.prepped <- plates.prepped[which(plates.prepped$id != "a7.345" & plates.prepped$id != "a07.345" & plates.prepped$id != "a11.421"),]
individual.info <- individual.info[which(individual.info$id != "a7.345" & individual.info$id != "a07.345" & individual.info$id != "a11.421"),]
Sraw <- Sraw[which(Sraw$sampleID != "a7.345" & Sraw$sampleID != "a07.345" & Sraw$sampleID != "a11.421"),]

#Some samples are sometimes called a and sometimes a0, for example:
plates.prepped$id[plates.prepped$id %in% c("a09.236","a9.236")] #but they are the same sample!
#This is fixed:
Sraw$runID <- gsub("a0" ,"a",Sraw$runID)
Sraw$sampleID <- gsub("a0" ,"a",Sraw$sampleID)
plates.prepped$id <- gsub("a0" ,"a",plates.prepped$id)
individual.info$id <- gsub("a0" ,"a",individual.info$id)

#find instances where sampleID has an extra zero at the end, and not in the beginning. So, it would detect if Animal.ID = 1090058 & sampleID = a9.580,
Sraw[Sraw$sampleID %in% stringr::str_subset(Sraw$sampleID, "0$") & !(Sraw$Animal.ID %in% stringr::str_subset(Sraw$Animal.ID, "0$")) & Sraw$age_cat == "Adult",]#This line finds sampleIDs with added zeros at the end of the id
#none! -> OK

#Some have lost their ending zero and appear as two different adults, but they are one (a13.320 and a13.32)
Sraw[Sraw$Animal.ID %in% stringr::str_subset(Sraw$Animal.ID, "0$") & !(Sraw$sampleID %in% stringr::str_subset(Sraw$sampleID, "0$")),]#This line looks for samples where the sampleID has lost final zeros.
# As you probably understand from the code above; sampleID, are computed from Animal.ID.
# Animal.ID has always seven digits: the first three digits are the year of birth (eg. 2010 = 110) and the last four are the animal identification number. Since this identification number needs to have four digits, ids with less than four digits will have added zeros at the beginning (eg. 6 = 0006 and 60 = 0060). The sampleID should be: age category + year + "." + indentification number without added zeros (eg. a10.6 or a10.60).

#a10.450.4 and a10.45.18.1
Sraw[Sraw$runID %in% "a10.45.18", c("runID","sampleID")] <- c("a10.450.18","a10.450")
Sraw[Sraw$runID %in% "a10.45.1", c("runID","sampleID")] <- c("a10.450.1","a10.450")
Sraw[Sraw$runID %in% "a10.45.18.1", c("runID","sampleID")] <- c("a10.450.18.1","a10.450")
plates.prepped$id[grep("a10.45",plates.prepped$id)] <- "a10.450"
individual.info$id[grep("a10.45",individual.info$id)] <- "a10.450"

#a10.360.4 and a10.36.18.1
Sraw[Sraw$runID %in% "a10.36.18", c("runID","sampleID")] <- c("a10.360.18","a10.360")
Sraw[Sraw$runID %in% "a10.36.1", c("runID","sampleID")] <- c("a10.360.1","a10.360")
Sraw[Sraw$runID %in% "a10.36.18.1", c("runID","sampleID")] <- c("a10.360.18.1","a10.360")
plates.prepped$id[plates.prepped$id %in% "a10.36"] <- "a10.360"
individual.info$id[individual.info$id %in% "a10.36"] <- "a10.360"

#a12.6.1 and a12.600.4
Sraw[Sraw$runID %in% "a12.6.1", c("runID","sampleID")] <- c("a12.600.1","a12.600")
plates.prepped$id[grep("a12.6",plates.prepped$id)] <- "a12.600"
individual.info$id[grep("a12.6",individual.info$id)] <- "a12.600"

#13.12
Sraw[grepl("a13.12", Sraw$runID) | grepl(paste("a13.32","0", sep = ""), Sraw$runID),]
Sraw[Sraw$runID %in% "a13.12.1", c("runID","sampleID")] <- c("a13.120.1","a13.120")
Sraw[Sraw$runID %in% "a13.12.4", c("runID","sampleID")] <- c("a13.120.4","a13.120")
plates.prepped$id[plates.prepped$id %in% "a13.12"] <- "a13.120"
individual.info$id[individual.info$id %in% "a13.12"] <- "a13.120"

#a13.32 and a13.320 also appears to be the same
Sraw[Sraw$runID %in% "a13.32.1", c("runID","sampleID")] <- c("a13.320.1","a13.320")

#a13.1.1 and a13.100.3
Sraw[Sraw$runID %in% "a13.1.1", c("runID","sampleID")] <- c("a13.100.1","a13.100")
plates.prepped$id[plates.prepped$id %in% "a13.1"] <- "a13.100"
individual.info$id[individual.info$id %in% "a13.1"] <- "a13.100"

#a9.17
Sraw[grepl("a9.17", Sraw$runID) | grepl(paste("a9.17","0", sep = ""), Sraw$runID),]
Sraw[Sraw$runID %in% "a9.17.1", c("runID","sampleID")] <- c("a9.170.1","a9.170")
Sraw[Sraw$runID %in% "a9.17.4", c("runID","sampleID")] <- c("a9.017.4","a9.017") #As animalid = 1090017
#It is only this 1090017 that has a9.17 in the other datasets, this is corrected
plates.prepped$id[plates.prepped$id %in% "a9.17" & plates.prepped$plateid %in% 4] <- "a9.017"
individual.info$id[individual.info$id %in% "a9.17" & individual.info$Plateid %in% 4] <- "a9.017"

#9.38
Sraw[grepl("a9.38", Sraw$runID) | grepl(paste("a9.38","0", sep = ""), Sraw$runID),]
Sraw[Sraw$runID %in% "a9.38.1", c("runID","sampleID")] <- c("a9.380.1","a9.380")
Sraw[Sraw$runID %in% "a9.38.4", c("runID","sampleID")] <- c("a9.038.4","a9.038") #As animalid = 1090038
#It is only this 1090038 that has a9.38 in the other datasets, this is corrected
plates.prepped$id[plates.prepped$id %in% "a9.38" & plates.prepped$plateid %in% 4] <- "a9.038"
individual.info$id[individual.info$id %in% "a9.38" & individual.info$Plateid %in% 4] <- "a9.038"

#a10.55.1 and a10.550.5
Sraw[Sraw$runID %in% "a10.55.1", c("runID","sampleID")] <- c("a10.550.1","a10.550")
plates.prepped$id[plates.prepped$id %in% "a10.55"] <- "a10.550"
individual.info$id[individual.info$id %in% "a10.55"] <- "a10.550"

#Check again:
Sraw[Sraw$Animal.ID %in% stringr::str_subset(Sraw$Animal.ID, "0$") & !(Sraw$sampleID %in% stringr::str_subset(Sraw$sampleID, "0$")),]#This line looks for samples where the sampleID has lost final zeros.
#None! _> OK

#------------------  Prep autosomal geotypes ------------------ 

#Create autosomal dataset
plates.prepped.auto <- plates.prepped[,!(grepl("Wnonpoly_",colnames(plates.prepped)) | grepl("Z_",colnames(plates.prepped)))]

#plate 1 uses Tube.label and not sampleID  as ID (Julian ran this plate before creating sample IDs).
#Fix
plates.prepped.auto$id[plates.prepped.auto$plateid %in% 1] <- Sraw$sampleID[match(plates.prepped.auto$id[plates.prepped.auto$plateid %in% 1],Sraw$Tube.Label)]
individual.info$id[individual.info$Plateid %in% 1 & !individual.info$id %in% c("Blank","NC")] <- Sraw$sampleID[match(individual.info$id[individual.info$Plateid %in% 1 & !individual.info$id %in% c("Blank","NC")],Sraw$Tube.Label)]

#Some individuals (i.e. samples), where run twice due to poor genotyping
#We will use the individual.info dataset to find the ones with most genotype calls and keep these
#Before we do this we need to create a unique identifer for each of these duplicates
plates.prepped.auto$runID <- paste(plates.prepped.auto$id,plates.prepped.auto$plateid, sep = ".")
#This created placespecific ids, as most samples were rerun on another place.
#However some samples were run several times at one plate, for example
Sraw[Sraw$runID %in% "a10.315.18.1",]
Sraw[Sraw$runID %in% "a10.315.18",]
#As seen above Julian differentiated these in the phenotypic dataset using this function (addtng .X)
plates.prepped.auto$runID <- make.unique(plates.prepped.auto$runID)
#Now it is unique and matches the info in phenotypic dataset
unique(duplicated(plates.prepped.auto$runID)) #FALSE = GOOD NO DUPLICATES
#Do the same for individual info
individual.info$runID <- paste(individual.info$id,individual.info$Plateid, sep = ".")
individual.info$runID <- make.unique(individual.info$runID)
unique(duplicated(individual.info$runID)) #FALSE = GOOD NO DUPLICATES

#Find duplicates by id
iddub <- unique(individual.info$id[duplicated(individual.info$id) & !is.na(individual.info$id) & !individual.info$id %in% c("PC","NC")])
inddub <- individual.info[individual.info$id %in% iddub,]

# and check which has least genotype calls.
remove <- as.numeric()
for(i in iddub){
  temp <- inddub[inddub$id %in% i,]
  keep <- rev(temp$runID[temp$NoGeno %in% min(temp$NoGeno)])[1] #If multiple of similar missing data take the latest run
  remove <- c(remove,temp$runID[!temp$runID %in% keep])
}

#Remove them
plates.prepped.auto <- plates.prepped.auto[!plates.prepped.auto$runID %in% remove,]
individual.info <- individual.info[!individual.info$runID %in% remove,]

#Check that each individual is only present once
unique(duplicated(plates.prepped.auto$id)) # FALSE = NO DUPS = GOOD
plates.prepped.auto$id[duplicated(plates.prepped.auto$id)]

unique(duplicated(individual.info$id[!individual.info$id %in% c("NC","PC")])) # FALSE = NO DUPS = GOOD
#Check that each individual is only present once as runID
unique(duplicated(plates.prepped.auto$runID)) # FALSE = NO DUPS = GOOD
unique(duplicated(individual.info$runID[!individual.info$runID %in% c("NC","PC")])) # FALSE = NO DUPS = GOOD

#Remove plate info and id
plates.prepped.auto$plateid <- NULL
plates.prepped.auto$id <- NULL
#Make runID new id column as it can be cross-referenced with phenotypic data set
plates.prepped.auto <- dplyr::rename(plates.prepped.auto, "id" = "runID")

#------------------  Print out individuals with less than X genotypes chosen for analysis ------------------ 

MinGeno = 60
individual.info$Geno <- 94-individual.info$NoGeno
rerun <- individual.info[individual.info$Geno <= MinGeno & !individual.info$id %in% c("NC",NA,"PC"),]
rerun$age_cat <- Sraw$age_cat[match(rerun$runID, Sraw$runID)]

# write.table(rerun,paste("output/Less than", MinGeno,"genotypes.csv", sep = " "), sep= ",", row.names = F)

###########################################
# 1. Parentage analyses with MasterBayes
###########################################

#--- Check that all individuals in phnotypic dateset is repressented by a genotype sample
Present <- Sraw$Animal.ID[Sraw$runID %in% plates.prepped.auto$id & !is.na(Sraw$Animal.ID)]
length(unique(Sraw$Animal.ID[!Sraw$Animal.ID %in% Present & !is.na(Sraw$Animal.ID)])) # 0 -> GOOD!


#------------------  Prep phenotypes ------------------ 


#Make a long version, with multiple records for each adult
YC <- colsplit(Sraw$yearcamps_present, split = ";", names = seq(1,4)) #Warning is OK!

Sraw2 <- unique(rbind(data.frame(Sraw, year.camp = YC$X1),
                      data.frame(Sraw, year.camp = YC$X2),
                      data.frame(Sraw, year.camp = YC$X3),
                      data.frame(Sraw, year.camp = YC$X4)
)) #Unique removes duplicates created by colsplit when only one yearcamp

S <- Sraw2 %>%
  dplyr::rename(c("sex" = "Sex","id" = "runID","offspring" = "age_cat")) %>%
  mutate(offspring = ifelse(offspring == "Adult", 0, 1)) %>%
  mutate(year = colsplit(year.camp, split = "_", names = c("year","out"))$year) %>%
  dplyr::select("id", "offspring", "sex","year.camp","year") %>%
  filter(id %in% plates.prepped.auto$id) #only keep phenotypes of individuals that we also have genotypes on

S$sex[S$sex %in% "M"] <- "Male"
S$sex[S$sex %in% "F"] <- "Female"

#Check all individuals in genotypes has phenotypes
unique(plates.prepped.auto$id %in% S$id) #TRUE -> They have

#Check all parents have an assigned sex!
length(S$sex[S$offspring %in% "0" & is.na(S$sex)]) #0 -> They have!


###########################################
# Step 1: Estimate probability of parents outside the camp
###########################################

#Genotype data (each column is one of two alleles)
MyGdP<-GdataPed(plates.prepped.auto)

#Indivdiuals from the offspring generation are excluded as parents
res.default<-expression(varPed("offspring", restrict=0))
#Male and Female need to come from same camp (however this is detached from the other parameters, so they just need to be together at any camp at any year)
res.yearcampparents<-expression(varPed("year.camp", restrict="==", relational = "MATE"))
##Associate above criteria with phenotypic dataframe
PdP.new2<-PdataPed(formula=list(res.default,res.yearcampparents), data=S, timevar=S$year, USsire=TRUE, USdam=TRUE)

#Because our model is quite simple, we dont need to do MCMC sampling, but can do with a CERVUS approximation
X.list2<-getXlist(PdP=PdP.new2, GdP=MyGdP, E1=0.005, E2=0.005)
ped.sp<-MLE.ped(X.list2, USsire=TRUE,USdam=TRUE, nUSsire=2, nUSdam = 2, threshold=0.95)

Parentage <- data.frame("runID" = ped.sp$P[,1], "dam.runID" = ped.sp$P[,2],"sire.runID" = ped.sp$P[,3], "parentage.prob" = ped.sp$prob)
Parentage <- Parentage[grep("a",Parentage$runID, invert = T),] #Remove adults

Parentage.missing <- Parentage[is.na(Parentage$dam.runID) | is.na(Parentage$sire.runID),]
nrow(Parentage.missing) #35 missing (17 if we remove non-samples sires and dams)
Parentage.called <- Parentage[!(is.na(Parentage$dam.runID) | is.na(Parentage$sire.runID)),]
nrow(Parentage.called) #1342 assigned
Parentage.called$analysis <- "Full"

#------------- Get more info

Sout <- merge(Sraw, Parentage, by = "runID", all = T)
#Only output those ids that we surveyed
Sout <- Sout[Sout$runID %in% plates.prepped.auto$id,]

Sout$damid <- Sout$Animal.ID[match(Sout$dam.runID, Sout$runID, incomparables = NA)]
Sout$sireid <- Sout$Animal.ID[match(Sout$sire.runID, Sout$runID, incomparables = NA)]
Sout$yrcamp.damid <- Sout$yearcamps_present[match(Sout$damid, Sout$Animal.ID, incomparables = NA)]
Sout$yrcamp.sireid <- Sout$yearcamps_present[match(Sout$sireid, Sout$Animal.ID, incomparables = NA)]
Sout$Geno <- individual.info$Geno[match(Sout$runID, individual.info$runID, incomparables = NA)]
#Difference between parents yr.camp and offspring yr.camp?
Sout$CampMisMatch <- NA
for(i in 1:nrow(Sout)){
  if(Sout$age_cat[i] != "Adult" & !is.na(Sout$yrcamp.damid[i])){ #For offspring assigned a parent
    Sout$CampMisMatch[i] <- ifelse(grepl( Sout$yearcamps_present[i],Sout$yrcamp.damid[i]) | grepl( Sout$yearcamps_present[i],Sout$yrcamp.sireid[i]), "Same","DiffentCamp")
  }  
}

#-------- Inspect mismatches and missing parentage

CampMismatch <- Sout[Sout$CampMisMatch %in% c("DiffentCamp"),c("runID","yearcamps_present","yrcamp.damid","yrcamp.sireid","Lay_date","parentage.prob","Geno")]
nrow(CampMismatch)# 16 Mismatches
nrow(CampMismatch)/nrow(Sout[Sout$age_cat != "Adult",])*100 #1.06% of assignments are from different camps
#Inspect
CampMismatch #All of the assigned parents are together in the same camp (but different from offspring) that year, although this was not a requirement!

#yr.camps of offspring with no calls
NoCall <- Sout[(is.na(Sout$dam.runID) | is.na(Sout$sire.runID)) & Sout$age_cat != "Adult",c("runID","yearcamps_present","parentage.prob","Geno")]
nrow(NoCall) #35
nrow(NoCall)/nrow(Sout[Sout$age_cat != "Adult",])*100 #2.54% unasigned
table(NoCall$yearcamps_present) #No strong pattern in camps with no calls, indicating that it is not the parents that is the problem

NoCall
#Many of them have low number of genotypes
#Many of them are from camps with missing parents


###########################################
# Step 2a: Rescue unassigned chicks by restricting parents to the given camp and focusing on camps with no unsampled individuals
###########################################

offspring.missing.camp <- S$id[S$year.camp %in% c("2016_20","2016_30","2017_70","2016_70","2017_63","2018_62") & S$offspring %in% 1]
plates.prepped.auto.nomiss <- plates.prepped.auto[!plates.prepped.auto$id %in% offspring.missing.camp,]

S.nomiss <- S[!(S$year.camp %in% c("2016_20","2016_30","2017_70","2016_70","2017_63","2018_62") & S$offspring %in% 1),]

MyGdP.nomiss<-GdataPed(plates.prepped.auto.nomiss)
#Indivdiuals not from the year-camp are excluded as parents
res.yearcamp<-expression(varPed("year.camp", restrict="==", relational = "OFFSPRING"))

PdP.new.yearcamp.nomiss<-PdataPed(formula=list(res.default,res.yearcampparents,res.yearcamp), data=S.nomiss, timevar=S.nomiss$year, USsire=F, USdam=F)
X.list.yearcamp.nomiss<-getXlist(PdP=PdP.new.yearcamp.nomiss, GdP=MyGdP.nomiss, E1=0.005, E2=0.005)
ped.sp.yearcamp.nomiss<-MLE.ped(X.list.yearcamp.nomiss, USsire=F,USdam=F, threshold=0.95)

Parentage.yearcamp.nomiss <- data.frame("runID" = ped.sp.yearcamp.nomiss$P[,1], "dam.runID" = ped.sp.yearcamp.nomiss$P[,2],"sire.runID" = ped.sp.yearcamp.nomiss$P[,3], "parentage.prob" = ped.sp.yearcamp.nomiss$prob)
Parentage.yearcamp.nomiss.called <- Parentage.yearcamp.nomiss[!(is.na(Parentage.yearcamp.nomiss$dam.runID) | is.na(Parentage.yearcamp.nomiss$sire.runID)),]
nrow(Parentage.yearcamp.nomiss.called)

#How many previously unassigned did we salvage
saved1 <- Parentage.yearcamp.nomiss.called[Parentage.yearcamp.nomiss.called$runID %in% Parentage.missing$runID,]
nrow(saved1) #14
#Add to called set
saved1$analysis <- "yearcamp"
Parentage.called <- rbind(Parentage.called,saved1)
#Remove from missing set
nrow(Parentage.missing) #35 before
Parentage.missing <- Parentage.missing[!Parentage.missing$runID %in% saved1$runID,]
nrow(Parentage.missing) #now 21 missing


###########################################
# Step 2a: Run camp specific models to assign offspring to unsampled parents
###########################################

#------ Camps missing one female ------#

#----- camp 2016_20
MyGdP.camp <- GdataPed(plates.prepped.auto[plates.prepped.auto$id %in% S$id[S$year.camp %in% c("2016_20")],])
S.camp <- S[S$year.camp %in% c("2016_20"),]

#The CERVUS approximation doesnt work here, as it does not assign the unsampled dam and doesnt even give all offspring a sire, despite just one male in the group...
# PdP.camp2016_20<-PdataPed(formula=list(res.default,res.yearcampparents), data=S.camp2016_20, USsire=FALSE, USdam=TRUE)
# X.list.camp2016_20<-getXlist(PdP=PdP.camp2016_20, GdP=MyGdP.camp2016_20, E1=0.005, E2=0.005)
# ped.sp<-MLE.ped(X.list.camp2016_20, USsire=FALSE,USdam=TRUE, nUSdam=1, nUSsire=0, threshold=0.95)

# So we go with the MCMC approach
sPunsam<-startPed(estG=FALSE, E1=0.005, E2=0.005, A=extractA(plates.prepped.auto), estUSsire = FALSE, estUSdam = FALSE, USdam = 1)
MyPdP<-PdataPed(formula=list(res.default,res.yearcampparents), data=S.camp, USdam = T, USsire = F)
m.sp.unsam<-MCMCped(PdP=MyPdP, GdP=MyGdP.camp, sP=sPunsam, nitt=35000, thin=100, burnin=3000, write_postG=TRUE) 
ped.sp<-modeP(m.sp.unsam$P, threshold=0.95)
Parentage.camp <- data.frame("runID" = ped.sp$P[,1], "dam.runID" = ped.sp$P[,2],"sire.runID" = ped.sp$P[,3], "parentage.prob" = ped.sp$prob)
#Okay so 3 with NA under dam, but inspection of model object gives this information
m.sp.unsam$P$c16.111.7
m.sp.unsam$P$c16.125.7
m.sp.unsam$P$c16.158.8
#All NA's here are the unsampled dam!
Parentage.camp$dam.runID[is.na(Parentage.camp$dam.runID)] <- "unsampled.dam"

#How many previously unassigned did we salvage
saved2 <- Parentage.camp[Parentage.camp$runID %in% Parentage.missing$runID,]
nrow(saved2) #2, not 3 as one egg was already assigned to another camp...
#Add to called set
saved2$analysis <- "camp-specific"
Parentage.called <- rbind(Parentage.called,saved2)
#Remove from missing set
Parentage.missing <- Parentage.missing[!Parentage.missing$runID %in% saved2$runID,]
nrow(Parentage.missing) #now 19 missing

#----- camp 2018_62
MyGdP.camp <- GdataPed(plates.prepped.auto[plates.prepped.auto$id %in% S$id[S$year.camp %in% c("2018_62")],])
S.camp <- S[S$year.camp %in% c("2018_62"),]

sPunsam<-startPed(estG=FALSE, E1=0.005, E2=0.005, A=extractA(plates.prepped.auto), estUSsire = FALSE, estUSdam = FALSE, USdam = 1)
MyPdP<-PdataPed(formula=list(res.default,res.yearcampparents), data=S.camp, USdam = T, USsire = F)
m.sp.unsam<-MCMCped(PdP=MyPdP, GdP=MyGdP.camp, sP=sPunsam, nitt=35000, thin=100, burnin=3000, write_postG=TRUE) 
ped.sp<-modeP(m.sp.unsam$P, threshold=0.95)
Parentage.camp <- data.frame("runID" = ped.sp$P[,1], "dam.runID" = ped.sp$P[,2],"sire.runID" = ped.sp$P[,3], "parentage.prob" = ped.sp$prob)
#Okay so 5 with NA under dam, but inspection of model object gives this information
m.sp.unsam$P$c18.445.17
m.sp.unsam$P$c18.522.17
m.sp.unsam$P$c18.712.18
m.sp.unsam$P$c18.786.18
m.sp.unsam$P$c18.792.18
#All NA's here are the unsampled dam!
Parentage.camp$dam.runID[is.na(Parentage.camp$dam.runID)] <- "unsampled.dam"

#How many previously unassigned did we salvage
saved3 <- Parentage.camp[Parentage.camp$runID %in% Parentage.missing$runID,]
nrow(saved3) #4, not 5 as one egg was already assigned to another camp...
#Add to called set
saved3$analysis <- "camp-specific"
Parentage.called <- rbind(Parentage.called,saved3)
#Remove from missing set
Parentage.missing <- Parentage.missing[!Parentage.missing$runID %in% saved3$runID,]
nrow(Parentage.missing) #now 15 missing

#------ Camps missing one male ------#

#----- camp 2016_30
MyGdP.camp <- GdataPed(plates.prepped.auto[plates.prepped.auto$id %in% S$id[S$year.camp %in% c("2016_30")],])
S.camp <- S[S$year.camp %in% c("2016_30"),]

sPunsam<-startPed(estG=FALSE, E1=0.005, E2=0.005, A=extractA(plates.prepped.auto), estUSsire = FALSE, estUSdam = FALSE, USsire = 1)
MyPdP<-PdataPed(formula=list(res.default,res.yearcampparents), data=S.camp, USdam = F, USsire = T)
m.sp.unsam<-MCMCped(PdP=MyPdP, GdP=MyGdP.camp, sP=sPunsam, nitt=35000, thin=100, burnin=3000, write_postG=TRUE) 
ped.sp<-modeP(m.sp.unsam$P, threshold=0.95)
Parentage.camp <- data.frame("runID" = ped.sp$P[,1], "dam.runID" = ped.sp$P[,2],"sire.runID" = ped.sp$P[,3], "parentage.prob" = ped.sp$prob)
#All genotyped parents

#How many previously unassigned did we salvage
saved4 <- Parentage.camp[Parentage.camp$runID %in% Parentage.missing$runID,]
nrow(saved4) #0


#----- camp 2017_70
MyGdP.camp <- GdataPed(plates.prepped.auto[plates.prepped.auto$id %in% S$id[S$year.camp %in% c("2017_70")],])
S.camp <- S[S$year.camp %in% c("2017_70"),]

sPunsam<-startPed(estG=FALSE, E1=0.005, E2=0.005, A=extractA(plates.prepped.auto), estUSsire = FALSE, estUSdam = FALSE, USsire = 1)
MyPdP<-PdataPed(formula=list(res.default,res.yearcampparents), data=S.camp, USdam = F, USsire = T)
m.sp.unsam<-MCMCped(PdP=MyPdP, GdP=MyGdP.camp, sP=sPunsam, nitt=35000, thin=100, burnin=3000, write_postG=TRUE) 
ped.sp<-modeP(m.sp.unsam$P, threshold=0.95)
Parentage.camp <- data.frame("runID" = ped.sp$P[,1], "dam.runID" = ped.sp$P[,2],"sire.runID" = ped.sp$P[,3], "parentage.prob" = ped.sp$prob)
#Okay so 5 with NA under dam, but inspection of model object gives this information
m.sp.unsam$P$c16.88.8
m.sp.unsam$P$c17.132.12
m.sp.unsam$P$c17.171.12
m.sp.unsam$P$c17.340.12
m.sp.unsam$P$c17.536.13
#All NA's here are the unsampled sire!
Parentage.camp$sire.runID[is.na(Parentage.camp$sire.runID)] <- "unsampled.sire"

#How many previously unassigned did we salvage
saved5 <- Parentage.camp[Parentage.camp$runID %in% Parentage.missing$runID,]
nrow(saved5) #5
#Add to called set
saved5$analysis <- "camp-specific"
Parentage.called <- rbind(Parentage.called,saved5)
#Remove from missing set
Parentage.missing <- Parentage.missing[!Parentage.missing$runID %in% saved5$runID,]
nrow(Parentage.missing) #now 10 missing


#----- camp 2016_70
MyGdP.camp <- GdataPed(plates.prepped.auto[plates.prepped.auto$id %in% S$id[S$year.camp %in% c("2016_70")],])
S.camp <- S[S$year.camp %in% c("2016_70"),]

sPunsam<-startPed(estG=FALSE, E1=0.005, E2=0.005, A=extractA(plates.prepped.auto), estUSsire = FALSE, estUSdam = FALSE, USsire = 1)
MyPdP<-PdataPed(formula=list(res.default,res.yearcampparents), data=S.camp, USdam = F, USsire = T)
m.sp.unsam<-MCMCped(PdP=MyPdP, GdP=MyGdP.camp, sP=sPunsam, nitt=35000, thin=100, burnin=3000, write_postG=TRUE) 
ped.sp<-modeP(m.sp.unsam$P, threshold=0.95)
Parentage.camp <- data.frame("runID" = ped.sp$P[,1], "dam.runID" = ped.sp$P[,2],"sire.runID" = ped.sp$P[,3], "parentage.prob" = ped.sp$prob)
#All genotyped parents

#How many previously unassigned did we salvage
saved6 <- Parentage.camp[Parentage.camp$runID %in% Parentage.missing$runID,]
nrow(saved6) #0


#----- camp 2017_63
MyGdP.camp <- GdataPed(plates.prepped.auto[plates.prepped.auto$id %in% S$id[S$year.camp %in% c("2017_63")],])
S.camp <- S[S$year.camp %in% c("2017_63"),]

sPunsam<-startPed(estG=FALSE, E1=0.005, E2=0.005, A=extractA(plates.prepped.auto), estUSsire = FALSE, estUSdam = FALSE, USsire = 1)
MyPdP<-PdataPed(formula=list(res.default,res.yearcampparents), data=S.camp, USdam = F, USsire = T)
m.sp.unsam<-MCMCped(PdP=MyPdP, GdP=MyGdP.camp, sP=sPunsam, nitt=35000, thin=100, burnin=3000, write_postG=TRUE) 
ped.sp<-modeP(m.sp.unsam$P, threshold=0.95)
Parentage.camp <- data.frame("runID" = ped.sp$P[,1], "dam.runID" = ped.sp$P[,2],"sire.runID" = ped.sp$P[,3], "parentage.prob" = ped.sp$prob)
#All genotyped parents and one uncalled
m.sp.unsam$P$c17.158.12
#This one is sired by unsampled male, but female remains unknown
Parentage.camp$sire.runID[is.na(Parentage.camp$sire.runID)] <- "unsampled.sire"

#How many previously unassigned did we salvage
saved7 <- Parentage.camp[Parentage.camp$runID %in% Parentage.missing$runID,]
nrow(saved7) #0
#Add to called set
saved5$analysis <- "camp-specific"
Parentage.called <- rbind(Parentage.called,saved5)
#Remove from missing set
Parentage.missing <- Parentage.missing[!Parentage.missing$runID %in% saved5$runID,]
nrow(Parentage.missing) #now 10 missing


###########################################
# 2. Sex determination
###########################################

#ZW = Female (homozygote at W and at Z)
#ZZ = Male (homozygote or heterozygote at Z, with hetero being evidence, and maybe no call at W)
#We have 6 Z chromosome markers and 2 W markers



#------------------  Prep sex genotypes ------------------ 

#plate 1 uses Tube.label and not sampleID  as ID (Julian run this plate before creating sample IDs).
#Fix
plates.sex$id[plates.sex$plateid %in% 1] <- Sraw$sampleID[match(plates.sex$id[plates.sex$plateid %in% 1],Sraw$Tube.Label)]

#Some individuals (i.e. samples), where run twice due to poor genotypig,and can be differentiated by runID
plates.sex$runID <- paste(plates.sex$id,plates.sex$plateid, sep = ".")

#Using dup info from autosomes (see above)
#We will use the sex genotypes from the samples that gave the most nr of genotypes, reducing error rate
#Remove them
plates.sex <- plates.sex[!plates.sex$runID %in% remove,]

#Remove plate info and id
plates.sex$plateid <- NULL
plates.sex$id <- NULL
plates.sex <- dplyr::rename(plates.sex, "id" = "runID")


#------------------  Infer sex ------------------ 


#Also add sex of offspring, and confirmed sex of parents (can be used to check how good we are at sexing)
#Go to long format
Psex <- gather(plates.sex, marker, genotype, colnames(plates.sex)[colnames(plates.sex) != "id"], factor_key=TRUE)
Psex$Chr[grepl("Wnonpoly_",Psex$marker)] <- "W"
Psex$Chr[grepl("Z_",Psex$marker)] <- "Z"
Psex$Call[Psex$genotype %in% c("A:A","G:G","C:C","T:T") & !is.na(Psex$genotype)] <- "homo"
Psex$Call[!Psex$genotype %in% c("A:A","G:G","C:C","T:T") & !is.na(Psex$genotype)] <- "hetero"
Psex$Call[Psex$genotype %in% "No Call" & !is.na(Psex$genotype)] <- "Nocall"

#Inspect 
table(Psex$Call[Psex$Chr %in% "W"]) #Only homo
table(Psex$Call[Psex$Chr %in% "Z"]) #Both

#Summarize
Psexsum <- data.frame(table(Psex$id,Psex$Chr,Psex$Call))
colnames(Psexsum) <- c("id","chr","call","count")
Ps <- spread(Psexsum, call, count)

PsZ <- Ps[Ps$chr %in% "Z",c("id","hetero","homo","Nocall")]
colnames(PsZ) <- c("id","Z.hetero","Z.homo","Z.Nocall")
PsW <- Ps[Ps$chr %in% "W",c("id","homo","Nocall")]
colnames(PsW) <- c("id","W.homo","W.Nocall")
P <- merge(PsW,PsZ, by= "id", all = T)


summary(P) #No NAs (if nas I need to adjust conde below!)


#---- Stats to inform criteria

table(P[P$Z.hetero >= 1,"W.Nocall"]) #2 -> all individuals that have 1 hetero Z genotype, have both W as no call, so no conflict in calling these a male!
table(P[P$W.homo >= 1,"Z.hetero"]) #0 -> all individuals that have one homo W, have no hetero Z, so no conflict in calling these females

#---- Criteria
#Female: hetero at no Z and homo at >=1 W
#Male: hetero at >=1 Z and no call at both W

P$sex.ass[P$Z.hetero >= 1 & P$W.Nocall == 2] <- "M"
P$sex.ass[P$W.homo >= 1 & P$Z.hetero == 0] <- "F"

#Which are we missing
nrow(P[is.na(P$sex.ass),c("id","W.Nocall","Z.hetero")])
P[is.na(P$sex.ass),c("id","W.Nocall","Z.hetero")]
#individuals which are likely males (no W) but where no Z genotypes were heterozygotes. But we cannot be sure

#---- Transfer to phenotype data

Soutsex <- merge(Sout, P[,c("id","sex.ass")], by.x = "runID", by.y = "id", all = T)
Soutsex[,c("age_cat","Sex","sex.ass")]
Soutsex$SexMisMatch[!is.na(Soutsex$Sex) & !is.na(Soutsex$sex.ass) & Soutsex$sex.ass != Soutsex$Sex] <- "DifferentSex"
sum(Soutsex$SexMisMatch %in% "DifferentSex") #0 -> None!


###########################################
# 3. Export
###########################################


write.table(Soutsex, paste("output/",InSamples,"_parents.csv", sep = ""), sep= ",", row.names = F)


###########################################
# Old MCMC approach
###########################################


#Fixing allele frequencies (to the observed) and the error rates (maybe we can get good estimates from related parents)
sP<-startPed(estG=FALSE, E1=0.005, E2=0.005, A=extractA(plates.prepped.auto))
sPunsam<-startPed(estG=FALSE, E1=0.005, E2=0.005, A=extractA(plates.prepped.auto), estUSsire = FALSE, USsire = 2, estUSdam = FALSE, USdam = 2)

#Associate above criteria with phenotypic dataframe
MyPdP<-PdataPed(formula=list(res.default,res.yearcampparents,res.year), data=S, timevar=S$year)
MyPdPunsam<-PdataPed(formula=list(res.default,res.yearcampparents,res.year), data=S, timevar=S$year, USdam = T, USsire = T)

#Run MCMCped
#m.sp<-MCMCped(PdP=MyPdP, GdP=MyGdP, sP=sP, nitt=35000, thin=100, burnin=3000, write_postG=TRUE) 
m.sp.unsam<-MCMCped(PdP=MyPd, GdP=MyGdP, sP=sPunsam, nitt=35000, thin=100, burnin=3000, write_postG=TRUE) 
m.sp.unsam2<-MCMCped(PdP=MyPdPunsam, GdP=MyGdP, sP=sPunsam, nitt=350, thin=1, burnin=30, write_postG=TRUE) 

#Inspect posterior of potential parent pairs for an example offspring:
ped.sp$P$c17.536.13 #NB if unncertain thenn two or more pairs should have "hits"

#Extract most likely pedigree
ped.sp<-modeP(m.sp.unsam$P, threshold=0.95)

#------------- Inspect posteriors and convergence
#Genotyping error (not estimated)
#autocorr(m.yearcamp$E1)
#autocorr(m.yearcamp$E2)
#plot(m.yearcamp$E2)
#plot(m.yearcamp$E1)


#Not used as not neccesary:
#Indivdiuals not from the year are excluded as parents
#res.year<-expression(varPed("year", restrict="==", relational = c("OFFSPRING")))
