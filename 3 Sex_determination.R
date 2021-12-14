#' ---
#' title: Sex determination of ostriches by genotyping
#' author: Mads Schou
#' date: 1/08/21
#' ---

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#
# 0. summary of approach -------------------------------------------------------
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#

# Criteria:
#   Female: 
#   - 1. hetero at 0 Z and homo at 2 W
# Male: 
#   - 1. hetero at >=1 Z and no call at both W (so not called homo at any W)
# or
# - 2. hetero at >=2 Z and one homo W accepted
# 
# The last male criteria is new and made post hoc looking at individuals that we could not call sex of.
# 
# High quality:
#   - No W genotype was ever called as heterozygote 
# - Samples with more than 32 missing genotypes was not trusted (indication of low quality for all genotypes) 
# - After above filter, no sex mismatch between multiple samples of one individual
# - No mismatch with phenotypic adult sex
# 
# Data:
#   - 1428 called individuals
# - 109 called offspring with mismatch to phenotypic sex out of 1248 called offspring (=8,7%)
# - 45 offspring where we could not call genotypic sex
# - 12 offspring with no phenotypic sex where we could not call genotypic sex
# 
# Bias in our ability to call sex:
#   - Of adult samples where we cant call genotypic sex only 1 out of 9 are females
# - Of chick samples where we cant call genotypic sex only 7 out of 23 are females
# It appears we have lower power when inferring sex of males than of females. But this is not a big problem as we call sex for almost all offspring (see above).
# 
# Quick Results:
#   - sex ratio of hatched chicks: 639 F / 608 M = 1.05
# - sex ratio of unhatched chicks: 46 F / 48 M = 0.96 ( + bias towards males in the 12 non-genotyped samples = lower ratio)


#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#
# 1. General data prep -------------------------------------------------------
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#


#Load phenotypic data
InSamples <- "geno_samples_2015.2018"
Sraw <- read.table(paste("input/", InSamples, ".csv", sep = ""), sep = ",", header = T)

#Load functions for genotypic data
source("1 dataprep functions.R")

MyPlates <- as.character(seq(1,18))
plates.sex <- as.numeric()
individual.info <- as.numeric()

for(i in MyPlates){
  P <- read.table(paste("input/Call map - Plate ", i, ".csv", sep = ""), skip = 115, sep = ",", header = T)
  P <- dplyr::rename(P, c("id" = "X.1"))
  P$X <- NULL
  
  individual.info <- rbind(individual.info,check.individual(platedf = P, plateid = i))
  plates.sex <- rbind(plates.sex,prep.plate.sex(platedf = P, plateid = i))
}

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#
# 2. Project specific data prep = RENAMING  ----------------------------------
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#

#Remove individuals a07.345 (OBS! Also written a7.345) and a11.421 since they are not part of the experimental population
plates.sex <- plates.sex[which(plates.sex$id != "a7.345" & plates.sex$id != "a07.345" & plates.sex$id != "a11.421"),]
individual.info <- individual.info[which(individual.info$id != "a7.345" & individual.info$id != "a07.345" & individual.info$id != "a11.421"),]
Sraw <- Sraw[which(Sraw$sampleID != "a7.345" & Sraw$sampleID != "a07.345" & Sraw$sampleID != "a11.421"),]

#Some samples are sometimes called a and sometimes a0, for example:
plates.sex$id[plates.sex$id %in% c("a09.236","a9.236")] #but they are the same sample!
#This is fixed:
Sraw$runID <- gsub("a0" ,"a",Sraw$runID)
Sraw$sampleID <- gsub("a0" ,"a",Sraw$sampleID)
plates.sex$id <- gsub("a0" ,"a",plates.sex$id)
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
plates.sex$id[grep("a10.45",plates.sex$id)] <- "a10.450"
individual.info$id[grep("a10.45",individual.info$id)] <- "a10.450"

#a10.360.4 and a10.36.18.1
Sraw[Sraw$runID %in% "a10.36.18", c("runID","sampleID")] <- c("a10.360.18","a10.360")
Sraw[Sraw$runID %in% "a10.36.1", c("runID","sampleID")] <- c("a10.360.1","a10.360")
Sraw[Sraw$runID %in% "a10.36.18.1", c("runID","sampleID")] <- c("a10.360.18.1","a10.360")
plates.sex$id[plates.sex$id %in% "a10.36"] <- "a10.360"
individual.info$id[individual.info$id %in% "a10.36"] <- "a10.360"

#a12.6.1 and a12.600.4
Sraw[Sraw$runID %in% "a12.6.1", c("runID","sampleID")] <- c("a12.600.1","a12.600")
plates.sex$id[grep("a12.6",plates.sex$id)] <- "a12.600"
individual.info$id[grep("a12.6",individual.info$id)] <- "a12.600"

#13.12
Sraw[grepl("a13.12", Sraw$runID) | grepl(paste("a13.32","0", sep = ""), Sraw$runID),]
Sraw[Sraw$runID %in% "a13.12.1", c("runID","sampleID")] <- c("a13.120.1","a13.120")
Sraw[Sraw$runID %in% "a13.12.4", c("runID","sampleID")] <- c("a13.120.4","a13.120")
plates.sex$id[plates.sex$id %in% "a13.12"] <- "a13.120"
individual.info$id[individual.info$id %in% "a13.12"] <- "a13.120"

#a13.32 and a13.320 also appears to be the same
Sraw[Sraw$runID %in% "a13.32.1", c("runID","sampleID")] <- c("a13.320.1","a13.320")

#a13.1.1 and a13.100.3
Sraw[Sraw$runID %in% "a13.1.1", c("runID","sampleID")] <- c("a13.100.1","a13.100")
plates.sex$id[plates.sex$id %in% "a13.1"] <- "a13.100"
individual.info$id[individual.info$id %in% "a13.1"] <- "a13.100"

#a9.17
Sraw[grepl("a9.17", Sraw$runID) | grepl(paste("a9.17","0", sep = ""), Sraw$runID),]
Sraw[Sraw$runID %in% "a9.17.1", c("runID","sampleID")] <- c("a9.170.1","a9.170")
Sraw[Sraw$runID %in% "a9.17.4", c("runID","sampleID")] <- c("a9.017.4","a9.017") #As animalid = 1090017
#It is only this 1090017 that has a9.17 in the other datasets, this is corrected
plates.sex$id[plates.sex$id %in% "a9.17" & plates.sex$plateid %in% 4] <- "a9.017"
individual.info$id[individual.info$id %in% "a9.17" & individual.info$Plateid %in% 4] <- "a9.017"

#9.38
Sraw[grepl("a9.38", Sraw$runID) | grepl(paste("a9.38","0", sep = ""), Sraw$runID),]
Sraw[Sraw$runID %in% "a9.38.1", c("runID","sampleID")] <- c("a9.380.1","a9.380")
Sraw[Sraw$runID %in% "a9.38.4", c("runID","sampleID")] <- c("a9.038.4","a9.038") #As animalid = 1090038
#It is only this 1090038 that has a9.38 in the other datasets, this is corrected
plates.sex$id[plates.sex$id %in% "a9.38" & plates.sex$plateid %in% 4] <- "a9.038"
individual.info$id[individual.info$id %in% "a9.38" & individual.info$Plateid %in% 4] <- "a9.038"

#a10.55.1 and a10.550.5
Sraw[Sraw$runID %in% "a10.55.1", c("runID","sampleID")] <- c("a10.550.1","a10.550")
plates.sex$id[plates.sex$id %in% "a10.55"] <- "a10.550"
individual.info$id[individual.info$id %in% "a10.55"] <- "a10.550"

#Check again:
Sraw[Sraw$Animal.ID %in% stringr::str_subset(Sraw$Animal.ID, "0$") & !(Sraw$sampleID %in% stringr::str_subset(Sraw$sampleID, "0$")),]#This line looks for samples where the sampleID has lost final zeros.
#None! _> OK

#Remove blank
plates.sex <- plates.sex[!plates.sex$id %in% "Blank",]

#plate 1 uses Tube.label and not sampleID  as ID (Julian ran this plate before creating sample IDs).
#Fix
plates.sex$id[plates.sex$plateid %in% 1] <- Sraw$sampleID[match(plates.sex$id[plates.sex$plateid %in% 1],Sraw$Tube.Label)]
individual.info$id[individual.info$Plateid %in% 1 & !individual.info$id %in% c("Blank","NC")] <- Sraw$sampleID[match(individual.info$id[individual.info$Plateid %in% 1 & !individual.info$id %in% c("Blank","NC")],Sraw$Tube.Label)]

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#
# 3. Project specific data prep = Duplicates  ----------------------------------
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#

#Some individuals (i.e. samples), where run twice due to poor genotyping
#We will use the individual.info dataset to find the ones with most genotype calls and keep these
#Before we do this we need to create a unique identifer for each of these duplicates
plates.sex$runID <- paste(plates.sex$id,plates.sex$plateid, sep = ".")
#This created placespecific ids, as most samples were rerun on another place.
#However some samples were run several times at one plate,
#As seen above Julian differentiated these by
plates.sex$runID <- make.unique(plates.sex$runID)
#Now it is unique and matches the info in phenotypic dataset
unique(duplicated(plates.sex$runID)) #FALSE = GOOD NO DUPLICATES
#Do the same for individual info
individual.info$runID <- paste(individual.info$id,individual.info$Plateid, sep = ".")
individual.info$runID <- make.unique(individual.info$runID)
unique(duplicated(individual.info$runID)) #FALSE = GOOD NO DUPLICATES

#Remove plate info and id
plates.sex$plateid <- NULL
plates.sex$id <- NULL
plates.sex <- dplyr::rename(plates.sex, "id" = "runID")

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#
# 4. Inspect sex chr data  ----------------------------------
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#


#ZW = Female (homozygote at W and at Z)
#ZZ = Male (homozygote or heterozygote at Z, with hetero being evidence, and maybe no call at W)
#We have 6 Z chromosome markers and 2 W markers

#Go to long format
Psex <- gather(plates.sex, marker, genotype, colnames(plates.sex)[colnames(plates.sex) != "id"], factor_key=TRUE)
Psex$Chr[grepl("Wnonpoly_",Psex$marker)] <- "W"
Psex$Chr[grepl("Z_",Psex$marker)] <- "Z"
Psex$Call[Psex$genotype %in% c("A:A","G:G","C:C","T:T") & !is.na(Psex$genotype)] <- "homo"
Psex$Call[!Psex$genotype %in% c("A:A","G:G","C:C","T:T") & !is.na(Psex$genotype)] <- "hetero"
Psex$Call[Psex$genotype %in% "No Call" & !is.na(Psex$genotype)] <- "Nocall"

#Inspect 
table(Psex$Call[Psex$Chr %in% "W"]) #Only homo as expected -> very low genotyping error rate!
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

summary(P) #No NAs (if nas I need to adjust code below!)

#---- Stats to inform criteria
#Looking for males
table(P[P$Z.hetero >= 1,"W.Nocall"]) # of the (753+8+34) individuals that have 1 hetero Z genotype, 753 have both W as no call, so there is a small but noteworthy error in assigning these as males
# 0   1   2 
# 8  34 753 
table(P[P$Z.hetero >= 2,"W.Nocall"]) # 
# 0   1   2 
# 8  28 651 
# Actually not a lot better and loosing many individuals
# This may indicate that two W nocalls and just one heterozygote Z is good to call a male!

#Looking for females, we see the same picture
table(P[P$W.homo >= 1,"Z.hetero"]) #also the wast majority (832) of individuals that have one homo W, have no hetero Z
# 0   1   2   3   4   5   6 
# 832   6  15  10   6   3   2 
table(P[P$W.homo >= 2,"Z.hetero"])
# 0   2   3   4   5 
# 823   5   1   1   1
#This is much better and with low loss of individuals.
#This may indicate that two W homo and zero Z hetero is the way to call males

#---- Obtain known sex of the chicks and adults and use this to design a set of criteria for the eggs
#First get eggid
Pp <- merge(Sraw, P, by.x = "runID", by.y = "id", all.y = T)
#Load database data with sex
egg <- read.table("../../..//OstrichDatabase/Database/Intermediate files/2. Crossreferenced datasets/Eggs1998-2020_merged_cleaned.csv", sep = ";", header = TRUE)
#Transfer sex
Pp$sex.eggid <- egg$sex[match(Pp$eggid,egg$eggid)]
#also genotype calls
Pp$NoGeno <- individual.info$NoGeno[match(Pp$runID,individual.info$runID)]
#Which are we missing sex for
table(Pp[is.na(Pp$sex.eggid) & is.na(Pp$Sex),"c_hatched"]) #Mostly unhatched eggs -> makes sense
#Which non adults do we have sex for:
table(Pp[!is.na(Pp$sex.eggid) & is.na(Pp$Sex),"c_hatched"]) #Mostly hatched eggs -> makes sense
#Merge sex info
Pp$sex.eggid[Pp$sex.eggid %in% 1] <- "M"
Pp$sex.eggid[Pp$sex.eggid %in% 2] <- "F"
Pp$Sex <- as.factor(ifelse(is.na(Pp$Sex),Pp$sex.eggid, Pp$Sex))
summary(Pp)
#Missing 95 sexes

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#
# 5. Sex determination  ----------------------------------
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#

#---- Criteria 1
Pp$sex.geno <- NA
#Female: hetero at no Z and homo at >=2 W
Pp$sex.geno[Pp$W.homo >= 2 & Pp$Z.hetero == 0] <- "F"
#Male: hetero at >=1 Z and no call at both W 
Pp$sex.geno[Pp$Z.hetero >= 1 & Pp$W.Nocall == 2] <- "M"

#--- Quality filter
Pp[Pp$NoGeno > 25, c("sex.geno","W.Nocall","Z.hetero")]
#All these are called males
#It appears that a high lack of genotypes can lead to errorneous no.call of W, supporting male
#Check if any of these called as males are confirmed or rejected by another sample
LowGenoID <- Pp$sampleID[Pp$NoGeno > 25 & !is.na(Pp$sex.geno)]

for(i in LowGenoID){
  temp <- Pp[Pp$sampleID %in% i,]
  if(nrow(temp[!is.na(temp$sex.geno),]) > 1){
    print(temp)
  }
}
#4 are confirmed, and two are rejected.
#But we can also compare to phenotypic sex
Pp[!is.na(Pp$Sex) & !is.na(Pp$sex.geno) & Pp$sex.geno != Pp$Sex & Pp$NoGeno > 25,]
#Only four with a mismatch, two of these are adults with >50 missing genotypes and the other two are chicks and could therefore just have the wrong phenotypic sex
table(Pp[!is.na(Pp$Sex) & !is.na(Pp$sex.geno) & Pp$sex.geno == Pp$Sex & Pp$NoGeno > 25,"NoGeno"]) #36 are confirmed by phenotypic sex!
#In conclusion no strong evidence that some missing genotypes leads to errornoeous call of male
#However we cannot say this is not the case at very high missing genotypes, and therefore sets the filter at 66% of the genoyptes called (=32)
Pp$sex.geno[Pp$NoGeno > 32] <- NA
#Note that this number doesnt matter as the individuals we cannot call sex of, unless decreased of course

#---- Check mismatch in sex within individuals with multiple samples
for(i in unique(Pp$sampleID)){
  temp <- Pp[Pp$sampleID %in% i,]
  if(length(unique(temp$sex.geno[!is.na(temp$sex.geno)])) > 1){
    print(temp)
  }
}
#None

#--- Any bias in sexes that we could not call, so calling as NA of adutls (so we are sure on actual sex), and removing where Ngenotype are above 48
NotCalled <-Pp[!is.na(Pp$Sex) & is.na(Pp$sex.geno) & Pp$age_cat %in% "Adult" & Pp$NoGeno <= 32,]
NotCalled
#These are mostly males, this is problematic at it creates a bias in what eggs are called
#note that most are on plate 1, note also that all on plate 1 has a homozygote W, even though they are males, only one not on plate 1 is homoW. 
#So something could just be strange with this plate.

#But maybe we rescue these with accepting one W, if at least two hetero Z
Pp[Pp$Z.hetero >= 2 & Pp$W.Nocall == 1 & Pp$NoGeno <=32, c("sex.geno","Sex","NoGeno","plate","age_cat")]
#yes all with known sex of these are M
#This appears like a good idea:

#---- Criteria set 2
Pp$sex.geno <- NA
#Female: hetero at no Z and homo at >=2 W
Pp$sex.geno[Pp$W.homo >= 2 & Pp$Z.hetero == 0] <- "F"
#Male: hetero at >=1 Z and no call at both W 
Pp$sex.geno[Pp$Z.hetero >= 1 & Pp$W.Nocall == 2] <- "M"
#Male: or hetero at >=2 Z and one homo W accepted
Pp$sex.geno[Pp$Z.hetero >= 2 & Pp$W.Nocall == 1] <- "M"
#Filtering
Pp$sex.geno[Pp$NoGeno > 32] <- NA

#---- Check mismatch in sex within individuals with multiple samples
for(i in unique(Pp$sampleID)){
  temp <- Pp[Pp$sampleID %in% i,]
  if(length(unique(temp$sex.geno[!is.na(temp$sex.geno)])) > 1){
    print(temp)
  }
}
#None

#We have 100% certainty in the sex of the adults. So we can use this to check the error of our approach
#However for the chicks hear say has it that around 5% are wrongly assigned sex by the workers
KnownCalled <- unique(Pp[!is.na(Pp$sex.geno) & !is.na(Pp$Sex),c("age_cat","sex.geno","Sex","sampleID")]) #Make uniqe so individual sampled multuple times onlny counts as one
nrow(KnownCalled) #1428 called individuals
KnownMismatch <- KnownCalled[KnownCalled$sex.geno != KnownCalled$Sex,]
nrow(KnownMismatch) #109 with a mismatch
#How are these distributed across eggs and chicks?
table(KnownMismatch$age_cat) #Only chicks, so all adults that we call, we call with 100% certainty. Hence chicks that are called wrongly must be due to mistakes by workers
#What is the error rate for chicks
nrow(KnownMismatch[KnownMismatch$age_cat %in% "Chick",])/nrow(KnownCalled[KnownCalled$age_cat %in% "Chick",])*100 #8,7%
#What sexes
table(KnownMismatch$Sex) #Mosly males that are called as females, which makes sense = if they dont find the penis they say female, but it could still be there.

#--- Any bias in sexes that we could not call, so calling as NA of adutls (so we are sure on actual sex), and removing where Ngenotype are above 48
NotCalled <-Pp[!is.na(Pp$Sex) & is.na(Pp$sex.geno) & Pp$age_cat %in% "Adult" & Pp$NoGeno <= 32,]
table(NotCalled$Sex) #Still some bias
# F M 
# 1 9 
NotCalled <-Pp[!is.na(Pp$Sex) & is.na(Pp$sex.geno) & Pp$age_cat %in% "Chick" & Pp$NoGeno <= 32,]
table(NotCalled$sex.eggid) #Also a bias when looking at the chicks
# F  M 
# 7 23 

#--- But above is on a sample bias, the importance of this bias depends on how many offspring (we are 100% certain of adult sex) this applies to
#How many individuals were not called (adults)
Called <-Pp[!is.na(Pp$sex.geno),]
length(unique(Pp[!Pp$sampleID %in% Called$sampleID & !Pp$age_cat %in% "Adult","sampleID"])) #45 missing
#But we know the sex of these with 91% certainity. So only 45*0.09 = 4 of these are wrongly assigned

#How many of these dont we have sex from another source
Pp[!Pp$sampleID %in% Called$sampleID & is.na(Pp$Sex),] #7 eggs missing, because bias most of these are like males. But we are talking about 2-3 eggs in bias

#What is the sex ratio 
table(Pp$sex.geno,Pp$c_hatched)


#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#
# 5. Data check  ----------------------------------
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#

#--- Transfer to phenotype data
Soutsex <- merge(Sraw, Pp[,c("runID","sex.eggid","sex.geno")], by.x = "runID", by.y = "runID", all = T)
Soutsex[,c("age_cat","Sex","sex.ass")]
write.table(Soutsex, paste("output/",InSamples,"_sexed.csv", sep = ""), sep= ",", row.names = F)


