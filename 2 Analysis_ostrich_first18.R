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
# 1. Parentage analyses with MasterBayes
###########################################

#------------------ Load data ------------------ 


load("platesprepped.RData") #genotypes

#phenotypes
InSamples <- "geno_samples_2015.2018"
Sraw <- read.table(paste("input/", InSamples, ".csv", sep = ""), sep = ",", header = T)



#------------------  Prep autosomal geotypes ------------------ 

#plate 1 uses Tube.label and not sampleID  as ID (Julian run this plate before creating sample IDs).
#Fix
plates.prepped.auto$id[plates.prepped.auto$plateid %in% 1] <- Sraw$sampleID[match(plates.prepped.auto$id[plates.prepped.auto$plateid %in% 1],Sraw$Tube.Label)]
individual.info$id[individual.info$Plateid %in% 1] <- Sraw$sampleID[match(individual.info$id[individual.info$Plateid %in% 1],Sraw$Tube.Label)]

#id a7.345 in plate 18 is written a both a7.345 and a07.345. Change all instances to a7.345
Sraw$sampleID[Sraw$sampleID %in% "a07.345"] <- "a7.345"
plates.prepped.auto$id[plates.prepped.auto$id %in% "a07.345"] <- "a7.345"
individual.info$id[individual.info$id %in% "a07.345"] <- "a7.345"

#Some individuals (i.e. samples), where run twice due to poor genotypig,and can be differentiated by runID
plates.prepped.auto$runID <- paste(plates.prepped.auto$id,plates.prepped.auto$plateid, sep = ".")
individual.info$runID <- paste(individual.info$id,individual.info$Plateid, sep = ".")
individual.info$runID <- make.unique(individual.info$runID)
individual.info$NoGeno <- individual.info$NoCall + individual.info$Invalid

#Find duplicates by id
iddub <- individual.info$id[duplicated(individual.info$id) & !is.na(individual.info$id) & !individual.info$id %in% c("PC","NC")]
inddub <- individual.info[individual.info$id %in% iddub,]

# and check which has least genotype calls.
remove <- as.numeric()
for(i in iddub){
  temp <- inddub[inddub$id %in% i,]
  keep <- temp$runID[temp$NoGeno %in% min(temp$NoGeno)][1]
  remove <- c(remove,temp$runID[!temp$runID %in% keep])
}

#Remove them
plates.prepped.auto <- plates.prepped.auto[!plates.prepped.auto$runID %in% remove,]
individual.info <- individual.info[!individual.info$runID %in% remove,]

#Remove plate info and id
plates.prepped.auto$plateid <- NULL
plates.prepped.auto$id <- NULL
plates.prepped.auto <- dplyr::rename(plates.prepped.auto, "id" = "runID")

#------------------  Print out individuals with less than X genotypes ------------------ 

MinGeno = 60
individual.info$Geno <- 94-individual.info$NoGeno
rerun <- individual.info[individual.info$Geno <= MinGeno & !individual.info$id %in% c("NC",NA),]
write.table(rerun,paste("output/Less than", MinGeno,"genotypes.csv", sep = " "), sep= ",", row.names = F)

#------------------  Prep phenotypes ------------------ 


#Make a long version, with multiple records for each adult
YC <- colsplit(Sraw$yearcamps_present, split = ";", names = seq(1,4))

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
#Check that all individuals.plates are only present once
unique(duplicated(plates.prepped.auto$id)) # TRUE FALSE -> SHOUDL BE TRUE ERRROR!!!!!!
#ALSO NEED TO CHECK IF INDIVIDUAL IS UNIQU

#------------------  Prep model ------------------ 


#Genotype data (each column is one of two alleles)
MyGdP<-GdataPed(plates.prepped.auto)
#---- Inspect
#MyGdP$G$Auto_superscaffold21_5705b
#P2[,"Auto_superscaffold21_5705a"]
#P2[,"Auto_superscaffold21_5705b"]

#Fixing allele frequencies (to the observed) and the error rates (maybe we can get good estimates from related parents)
sP<-startPed(estG=FALSE, E1=0.005, E2=0.005, A=extractA(plates.prepped.auto))

#Indivdiuals from the offspring generation are excluded as parents
res.default<-expression(varPed("offspring", restrict=0))

#Indivdiuals not from the year-camp are excluded as parents
res.yearcamp<-expression(varPed("year.camp", restrict="==", relational = "OFFSPRING"))

#Male and Female need to come from same camp (however this is distached from the other parameters, so they just need to be together at any camp at any year)
res.yearcampparents<-expression(varPed("year.camp", restrict="==", relational = "MATE"))

#Indivdiuals not from the year are excluded as parents
res.year<-expression(varPed("year", restrict="==", relational = c("OFFSPRING")))

#Associate above criteria with phenotypic dataframe
MyPdPyearcamp<-PdataPed(formula=list(res.default,res.yearcampparents,res.year,res.yearcamp), data=S, timevar=S$year) #This uses our prior knowledge on camp origin, so doesnt allow for contamination across camps
MyPdP<-PdataPed(formula=list(res.default,res.yearcampparents,res.year), data=S, timevar=S$year)

#Run MCMCped
#model1_old<-MCMCped(PdP=MyPdP, GdP=MyGdP, nitt=15000, thin=100, burnin=3000, write_postG=TRUE) 
m.sp<-MCMCped(PdP=MyPdP, GdP=MyGdP, sP=sP, nitt=35000, thin=100, burnin=3000, write_postG=TRUE) 
#m.yearcamp<-MCMCped(PdP=MyPdPyearcamp, GdP=MyGdP, sP=sP, nitt=35000, thin=100, burnin=3000, write_postG=TRUE)

#Inspect posterior of potential parent pairs for an example offspring:
m.sp$P$c15.344 #NB if unncertain thenn two or more pairs should have "hits"

#Extract most likely pedigree
ped.sp<-modeP(m.sp$P, threshold=0.9)
#ped1.yearcamp<-modeP(m.yearcamp$P, threshold=0.8)
Parentage <- ped.sp$P
nrow(Parentage[is.na(Parentage$sire.runID),])

#------------- Inspect posteriors and convergence
#Genotyping error (not estimated)
#autocorr(model1$E1)
#autocorr(model1$E2)
#plot(model1$E2)
#plot(model1$E1)

#------------- Transfer parentage 

colnames(Parentage) <- c("runID","dam.runID","sire.runID")
Parentage <- data.frame(Parentage, "parentage.prob" = ped.sp$prob)
Sout <- merge(Sraw, Parentage, by = "runID", all = T)
#Only output those ids that we surveyed
Sout <- Sout[Sout$runID %in% plates.prepped.auto$id,]

Sout$damid <- Sout$Animal.ID[match(Sout$dam.runID, Sout$runID, incomparables = NA)]
Sout$sireid <- Sout$Animal.ID[match(Sout$sire.runID, Sout$runID, incomparables = NA)]
Sout$yrcamp.damid <- Sout$yearcamps_present[match(Sout$damid, Sout$Animal.ID, incomparables = NA)]
Sout$yrcamp.sireid <- Sout$yearcamps_present[match(Sout$sireid, Sout$Animal.ID, incomparables = NA)]
#Difference between parents yr.camp and offspring yr.camp?
Sout$CampMisMatch <- NA
for(i in 1:nrow(Sout)){
  if(Sout$age_cat[i] != "Adult" & !is.na(Sout$yrcamp.damid[i])){ #For offspring assigned a parent
    Sout$CampMisMatch[i] <- ifelse(grepl( Sout$yearcamps_present[i],Sout$yrcamp.damid[i]) | grepl( Sout$yearcamps_present[i],Sout$yrcamp.sireid[i]), "Same","DiffentCamp")
  }  
}

#-------- Inspect mismatches and missing parentage

Sout[Sout$CampMisMatch %in% c("DiffentCamp"),]
#MISMATCHES
#These are expected, the first on in camp 2015_20 because two males have almost no calls, the second the predicted dam (a12.195) also have very low genotype calls.

#yr.camps of offspring with no calls
nrow(Sout[is.na(Sout$dam.runID) & Sout$age_cat != "Adult",c("runID","yearcamps_present")]) #26
Sout[is.na(Sout$dam.runID) & Sout$age_cat != "Adult",c("runID","yearcamps_present")] 
#All the 2015_20 are due to two (a10.599 and a13.212) of three adult males being very poor quality samples

#Two are due to low quality of offspring sample
#c15.12           2015_25
#c17.840           2017_29

#Two are high quality samples but inspection of their posterior, shows that the poor quality adults (2 males and one female) somehow becomes prime candidates
#349    c15.54           2015_66 -> high quality offspring sample!
#409   c15.663           2015_66 -> high quality offspring sample!
m.sp$P$c15.54 
m.sp$P$c15.663 
#These low quality males and female are not from the same camp as the offspring, so using the fixed year.camp model gives better results:
m.yearcamp$P$c15.54 
m.yearcamp$P$c15.663 

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


