
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

