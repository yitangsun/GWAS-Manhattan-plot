library(qqman)
library(qqplotr)
library(RColorBrewer)
library(TwoSampleMR)
library(dplyr)
library(googlesheets)
library(vcfR)
library(stringr)
library(mr.raps)
'%ni%' <- Negate('%in%')
# ieugwasr::api_status()
# $`API version`
# #3.6.7
# $`Total associations`
# [1] 126335269652
# $`Total complete datasets`
# [1] 34670
# $`Total public datasets`
# [1] 34513
library(MRInstruments)
library(MVMR)
#library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
#library("SNPlocs.Hsapiens.dbSNP151.GRCh38")

Pathway_geno="/Users/yelab/GWAS/FUMA/"

#PUFA_LIST=c("PUFA","Omega3","Omega6","LA","DHA")
n=c("Omega3")

Trait1_final_infile=paste(Pathway_geno,"FUMA_",n,"/gwascatalog.txt", sep="")
Trait1_final <- read.csv(Trait1_final_infile,header=T, as.is=T, sep = "\t")

Trait_list=sort(Trait1_final[!duplicated(Trait1_final$Trait),]$Trait)


grep("acid", Trait_list, value = T)

grep("fatty", Trait_list, value = T)

final_trait=grep("polyunsaturated fatty acid", Trait_list, value = T)
final_trait1=grep("Red blood cell fatty acid levels",  Trait_list, value = T)

final_trait=append(final_trait,final_trait1)

GWAS_cata_PUFA=Trait1_final[Trait1_final$Trait %in% final_trait,]
#names(GWAS_cata_PUFA)[names(GWAS_cata_PUFA) == "IndSigSNP"] <- "rsID"
names(GWAS_cata_PUFA)[names(GWAS_cata_PUFA) == "chr"] <- "CHR"
#names(GWAS_cata_PUFA)[names(GWAS_cata_PUFA) == "GenomicLocus"] <- "NumGenomicLocus"

Trait_final_infile=paste(Pathway_geno,"FUMA_",n,"/GenomicRiskLoci.txt", sep="")
Trait_final <- read.csv(Trait_final_infile,header=T, as.is=T, sep = "\t")

Trait_final <- Trait_final %>% left_join(GWAS_cata_PUFA, by= "GenomicLocus")
Outputfile=paste(Pathway_geno,"FUMA_",n,"/",n,"_Novelty.txt", sep="")
write.table(Trait_final, file= Outputfile, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')

#PUFA_LIST=c("PUFA","Omega3","Omega6","LA","DHA")
n=c("PUFA")

Trait1_final_infile=paste(Pathway_geno,"FUMA_",n,"/gwascatalog.txt", sep="")
Trait1_final <- read.csv(Trait1_final_infile,header=T, as.is=T, sep = "\t")

Trait_list=sort(Trait1_final[!duplicated(Trait1_final$Trait),]$Trait)

grep("acid", Trait_list, value = T)

grep("fatty", Trait_list, value = T)

final_trait=grep("polyunsaturated fatty acid", Trait_list, value = T)
final_trait1=grep("Red blood cell fatty acid levels",  Trait_list, value = T)

final_trait=append(final_trait,final_trait1)

GWAS_cata_PUFA=Trait1_final[Trait1_final$Trait %in% final_trait,]
#names(GWAS_cata_PUFA)[names(GWAS_cata_PUFA) == "IndSigSNP"] <- "rsID"
names(GWAS_cata_PUFA)[names(GWAS_cata_PUFA) == "chr"] <- "CHR"
#names(GWAS_cata_PUFA)[names(GWAS_cata_PUFA) == "GenomicLocus"] <- "NumGenomicLocus"

Trait_final_infile=paste(Pathway_geno,"FUMA_",n,"/GenomicRiskLoci.txt", sep="")
Trait_final <- read.csv(Trait_final_infile,header=T, as.is=T, sep = "\t")

Trait_final <- Trait_final %>% left_join(GWAS_cata_PUFA, by= "GenomicLocus")
Outputfile=paste(Pathway_geno,"FUMA_",n,"/",n,"_Novelty.txt", sep="")
write.table(Trait_final, file= Outputfile, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')

#PUFA_LIST=c("PUFA","Omega3","Omega6","LA","DHA")
n=c("Omega6")

Trait1_final_infile=paste(Pathway_geno,"FUMA_",n,"/gwascatalog.txt", sep="")
Trait1_final <- read.csv(Trait1_final_infile,header=T, as.is=T, sep = "\t")

Trait_list=sort(Trait1_final[!duplicated(Trait1_final$Trait),]$Trait)

grep("acid", Trait_list, value = T)

grep("fatty", Trait_list, value = T)

final_trait=grep("polyunsaturated fatty acid", Trait_list, value = T)
final_trait1=grep("Red blood cell fatty acid levels",  Trait_list, value = T)

final_trait=append(final_trait,final_trait1)

GWAS_cata_PUFA=Trait1_final[Trait1_final$Trait %in% final_trait,]
#names(GWAS_cata_PUFA)[names(GWAS_cata_PUFA) == "IndSigSNP"] <- "rsID"
names(GWAS_cata_PUFA)[names(GWAS_cata_PUFA) == "chr"] <- "CHR"
#names(GWAS_cata_PUFA)[names(GWAS_cata_PUFA) == "GenomicLocus"] <- "NumGenomicLocus"

Trait_final_infile=paste(Pathway_geno,"FUMA_",n,"/GenomicRiskLoci.txt", sep="")
Trait_final <- read.csv(Trait_final_infile,header=T, as.is=T, sep = "\t")

Trait_final <- Trait_final %>% left_join(GWAS_cata_PUFA, by= "GenomicLocus")
Outputfile=paste(Pathway_geno,"FUMA_",n,"/",n,"_Novelty.txt", sep="")
write.table(Trait_final, file= Outputfile, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')

#PUFA_LIST=c("PUFA","Omega3","Omega6","LA","DHA")
n=c("LA")

Trait1_final_infile=paste(Pathway_geno,"FUMA_",n,"/gwascatalog.txt", sep="")
Trait1_final <- read.csv(Trait1_final_infile,header=T, as.is=T, sep = "\t")

Trait_list=sort(Trait1_final[!duplicated(Trait1_final$Trait),]$Trait)

grep("acid", Trait_list, value = T)

grep("fatty", Trait_list, value = T)

final_trait=grep("polyunsaturated fatty acid", Trait_list, value = T)
final_trait1=grep("Red blood cell fatty acid levels",  Trait_list, value = T)

final_trait=append(final_trait,final_trait1)

GWAS_cata_PUFA=Trait1_final[Trait1_final$Trait %in% final_trait,]
#names(GWAS_cata_PUFA)[names(GWAS_cata_PUFA) == "IndSigSNP"] <- "rsID"
names(GWAS_cata_PUFA)[names(GWAS_cata_PUFA) == "chr"] <- "CHR"
#names(GWAS_cata_PUFA)[names(GWAS_cata_PUFA) == "GenomicLocus"] <- "NumGenomicLocus"

Trait_final_infile=paste(Pathway_geno,"FUMA_",n,"/GenomicRiskLoci.txt", sep="")
Trait_final <- read.csv(Trait_final_infile,header=T, as.is=T, sep = "\t")

Trait_final <- Trait_final %>% left_join(GWAS_cata_PUFA, by= "GenomicLocus")
Outputfile=paste(Pathway_geno,"FUMA_",n,"/",n,"_Novelty.txt", sep="")
write.table(Trait_final, file= Outputfile, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')

#PUFA_LIST=c("PUFA","Omega3","Omega6","LA","DHA")
n=c("DHA")

Trait1_final_infile=paste(Pathway_geno,"FUMA_",n,"/gwascatalog.txt", sep="")
Trait1_final <- read.csv(Trait1_final_infile,header=T, as.is=T, sep = "\t")

Trait_list=sort(Trait1_final[!duplicated(Trait1_final$Trait),]$Trait)

grep("acid", Trait_list, value = T)

grep("fatty", Trait_list, value = T)

final_trait=grep("polyunsaturated fatty acid", Trait_list, value = T)
final_trait1=grep("Red blood cell fatty acid levels",  Trait_list, value = T)

final_trait=append(final_trait,final_trait1)

GWAS_cata_PUFA=Trait1_final[Trait1_final$Trait %in% final_trait,]
#names(GWAS_cata_PUFA)[names(GWAS_cata_PUFA) == "IndSigSNP"] <- "rsID"
names(GWAS_cata_PUFA)[names(GWAS_cata_PUFA) == "chr"] <- "CHR"
#names(GWAS_cata_PUFA)[names(GWAS_cata_PUFA) == "GenomicLocus"] <- "NumGenomicLocus"

Trait_final_infile=paste(Pathway_geno,"FUMA_",n,"/GenomicRiskLoci.txt", sep="")
Trait_final <- read.csv(Trait_final_infile,header=T, as.is=T, sep = "\t")

Trait_final <- Trait_final %>% left_join(GWAS_cata_PUFA, by= "GenomicLocus")
Outputfile=paste(Pathway_geno,"FUMA_",n,"/",n,"_Novelty.txt", sep="")
write.table(Trait_final, file= Outputfile, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
