
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

ao <- available_outcomes()
# data(gwas_catalog)
# head(gwas_catalog)

Pathway_SNP="/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/"
Pathway_geno="/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/Genotype/"
Pathway_out="/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/MR_result/result_GWAS_05_19/"
COVID_LIST=c("HGI_round_4_A1","HGI_round_4_A2","HGI_round_4_B1","HGI_round_4_B2","HGI_round_4_C1","HGI_round_4_C2","HGI_round_5_A2_eur","HGI_round_5_A2_eur_leave_ukbb","HGI_round_5_B1_eur","HGI_round_5_B1_eur_leave_ukbb","HGI_round_5_B2_eur","HGI_round_5_B2_eur_leave_ukbb","HGI_round_5_C2_eur","HGI_round_5_C2_eur_leave_ukbb")

NEJM_LIST=c("data_test10")

Out_Geno_filename="_All_trait_03_21.txt"
Old_Out_Geno_filename="_All_trait.csv"
Trait_filename="All_Trait_IEU_GWAS.txt"

#PUFA_LIST=c("PUFA","Omega_3","Omega_6","LA","DHA")
#PUFA 0.0616931 0.09575667 
#Omega_3 0.05280874 1.095328 
#Omega_6 0.05295491 0.1137094 
#LA 0.04952331 0.1052528 
#DHA 0.03775332 5.49533 

#######Change
for (n in c("PUFA")) {
  Trait_final_infile=paste(Pathway_out,"Sta_",n,".txt", sep="")
  Trait_final <- read.table(Trait_final_infile,header=T, as.is=T,sep = "\t")
  print(str(Trait_final))
  Trait_final$UKB_AF=sapply(strsplit(Trait_final$V10, split= ":", fixed=TRUE),"[",4)
  #Trait_final$UKB_AF=ifelse(Trait_final$Final_NS>115077,sapply(strsplit(Trait_final$V10, split= ":", fixed=TRUE),"[",4),NA)
  
  Trait_final$Final_EAF=ifelse(Trait_final$Final_NS>115077,Trait_final$UKB_AF,Trait_final$MAF)
  Trait_final$Final_EAF=as.numeric(Trait_final$Final_EAF)
  print(str(Trait_final))
  
  Trait_final$Final_R2=2*(Trait_final$Effect)*(Trait_final$Effect)*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)
  
  #######Change
  a=0.802664
  
  Trait_final$Final_R2_SD=2*(Trait_final$Effect)*(Trait_final$Effect)*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)/(a*a)
  
  b=2*(Trait_final$Effect)*(Trait_final$Effect)*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)+(Trait_final$StdErr)*(Trait_final$StdErr)*2*Trait_final$Final_NS*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)
  Trait_final$PVE=2*(Trait_final$Effect)*(Trait_final$Effect)*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)/b
  
  print(n)
  print(sum(Trait_final$Final_R2))
  print(sum(Trait_final$Final_R2_SD))
  print(sum(Trait_final$PVE))
  
  names(Trait_final)[names(Trait_final) == "SNP"] <- "rsID"
  
  Outputfile2=paste(Pathway_out,"GenomicRiskLoci_",n, ".txt", sep="")
  Trait_final2 <- read.table(Outputfile2,header=T, as.is=T,sep = "\t")
  
  Trait_final2<- Trait_final %>% right_join(Trait_final2, by= "rsID")
  print(str(Trait_final2))
  
  print(sum(Trait_final2$Final_R2))
  print(sum(Trait_final2$Final_R2_SD))
  print(sum(Trait_final2$PVE))
  
  Outputfile3=paste(Pathway_out,"R2_",n, ".txt", sep="")
  write.table(Trait_final, file= Outputfile3, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
  
  Outputfile4=paste(Pathway_out,"GenomicRiskLoci_R2_",n, ".txt", sep="")
  write.table(Trait_final2, file= Outputfile4, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
}


#######Change
for (n in c("Omega_3")) {
  Trait_final_infile=paste(Pathway_out,"Sta_",n,".txt", sep="")
  Trait_final <- read.table(Trait_final_infile,header=T, as.is=T,sep = "\t")
  print(str(Trait_final))
  Trait_final$UKB_AF=sapply(strsplit(Trait_final$V10, split= ":", fixed=TRUE),"[",4)
  #Trait_final$UKB_AF=ifelse(Trait_final$Final_NS>115077,sapply(strsplit(Trait_final$V10, split= ":", fixed=TRUE),"[",4),NA)
  
  Trait_final$Final_EAF=ifelse(Trait_final$Final_NS>115077,Trait_final$UKB_AF,Trait_final$MAF)
  Trait_final$Final_EAF=as.numeric(Trait_final$Final_EAF)
  print(str(Trait_final))
  
  Trait_final$Final_R2=2*(Trait_final$Effect)*(Trait_final$Effect)*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)
  
  #######Change
  a=0.219574
  
  Trait_final$Final_R2_SD=2*(Trait_final$Effect)*(Trait_final$Effect)*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)/(a*a)
  
  b=2*(Trait_final$Effect)*(Trait_final$Effect)*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)+(Trait_final$StdErr)*(Trait_final$StdErr)*2*Trait_final$Final_NS*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)
  Trait_final$PVE=2*(Trait_final$Effect)*(Trait_final$Effect)*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)/b
  
  print(n)
  print(sum(Trait_final$Final_R2))
  print(sum(Trait_final$Final_R2_SD))
  print(sum(Trait_final$PVE))
  
  names(Trait_final)[names(Trait_final) == "SNP"] <- "rsID"
  
  Outputfile2=paste(Pathway_out,"GenomicRiskLoci_",n, ".txt", sep="")
  Trait_final2 <- read.table(Outputfile2,header=T, as.is=T,sep = "\t")
  
  Trait_final2<- Trait_final %>% right_join(Trait_final2, by= "rsID")
  print(str(Trait_final2))
  
  print(sum(Trait_final2$Final_R2))
  print(sum(Trait_final2$Final_R2_SD))
  print(sum(Trait_final2$PVE))
  
  Outputfile3=paste(Pathway_out,"R2_",n, ".txt", sep="")
  write.table(Trait_final, file= Outputfile3, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
  
  Outputfile4=paste(Pathway_out,"GenomicRiskLoci_R2_",n, ".txt", sep="")
  write.table(Trait_final2, file= Outputfile4, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
}


#######Change
for (n in c("Omega_6")) {
  Trait_final_infile=paste(Pathway_out,"Sta_",n,".txt", sep="")
  Trait_final <- read.table(Trait_final_infile,header=T, as.is=T,sep = "\t")
  print(str(Trait_final))
  Trait_final$UKB_AF=sapply(strsplit(Trait_final$V10, split= ":", fixed=TRUE),"[",4)
  #Trait_final$UKB_AF=ifelse(Trait_final$Final_NS>115077,sapply(strsplit(Trait_final$V10, split= ":", fixed=TRUE),"[",4),NA)
  
  Trait_final$Final_EAF=ifelse(Trait_final$Final_NS>115077,Trait_final$UKB_AF,Trait_final$MAF)
  Trait_final$Final_EAF=as.numeric(Trait_final$Final_EAF)
  print(str(Trait_final))
  
  Trait_final$Final_R2=2*(Trait_final$Effect)*(Trait_final$Effect)*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)
  
  #######Change
  a=0.682425
  
  Trait_final$Final_R2_SD=2*(Trait_final$Effect)*(Trait_final$Effect)*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)/(a*a)
  
  b=2*(Trait_final$Effect)*(Trait_final$Effect)*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)+(Trait_final$StdErr)*(Trait_final$StdErr)*2*Trait_final$Final_NS*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)
  Trait_final$PVE=2*(Trait_final$Effect)*(Trait_final$Effect)*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)/b
  
  print(n)
  print(sum(Trait_final$Final_R2))
  print(sum(Trait_final$Final_R2_SD))
  print(sum(Trait_final$PVE))
  
  names(Trait_final)[names(Trait_final) == "SNP"] <- "rsID"
  
  Outputfile2=paste(Pathway_out,"GenomicRiskLoci_",n, ".txt", sep="")
  Trait_final2 <- read.table(Outputfile2,header=T, as.is=T,sep = "\t")
  
  Trait_final2<- Trait_final %>% right_join(Trait_final2, by= "rsID")
  print(str(Trait_final2))
  
  print(sum(Trait_final2$Final_R2))
  print(sum(Trait_final2$Final_R2_SD))
  print(sum(Trait_final2$PVE))
  
  Outputfile3=paste(Pathway_out,"R2_",n, ".txt", sep="")
  write.table(Trait_final, file= Outputfile3, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
  
  Outputfile4=paste(Pathway_out,"GenomicRiskLoci_R2_",n, ".txt", sep="")
  write.table(Trait_final2, file= Outputfile4, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
}

#######Change
for (n in c("LA")) {
  Trait_final_infile=paste(Pathway_out,"Sta_",n,".txt", sep="")
  Trait_final <- read.table(Trait_final_infile,header=T, as.is=T,sep = "\t")
  print(str(Trait_final))
  Trait_final$UKB_AF=sapply(strsplit(Trait_final$V10, split= ":", fixed=TRUE),"[",4)
  #Trait_final$UKB_AF=ifelse(Trait_final$Final_NS>115077,sapply(strsplit(Trait_final$V10, split= ":", fixed=TRUE),"[",4),NA)
  
  Trait_final$Final_EAF=ifelse(Trait_final$Final_NS>115077,Trait_final$UKB_AF,Trait_final$MAF)
  Trait_final$Final_EAF=as.numeric(Trait_final$Final_EAF)
  print(str(Trait_final))
  
  Trait_final$Final_R2=2*(Trait_final$Effect)*(Trait_final$Effect)*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)
  
  #######Change
  a=0.685943
  
  Trait_final$Final_R2_SD=2*(Trait_final$Effect)*(Trait_final$Effect)*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)/(a*a)
  
  b=2*(Trait_final$Effect)*(Trait_final$Effect)*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)+(Trait_final$StdErr)*(Trait_final$StdErr)*2*Trait_final$Final_NS*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)
  Trait_final$PVE=2*(Trait_final$Effect)*(Trait_final$Effect)*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)/b
  
  print(n)
  print(sum(Trait_final$Final_R2))
  print(sum(Trait_final$Final_R2_SD))
  print(sum(Trait_final$PVE))
  
  names(Trait_final)[names(Trait_final) == "SNP"] <- "rsID"
  
  Outputfile2=paste(Pathway_out,"GenomicRiskLoci_",n, ".txt", sep="")
  Trait_final2 <- read.table(Outputfile2,header=T, as.is=T,sep = "\t")
  
  Trait_final2<- Trait_final %>% right_join(Trait_final2, by= "rsID")
  print(str(Trait_final2))
  
  print(sum(Trait_final2$Final_R2))
  print(sum(Trait_final2$Final_R2_SD))
  print(sum(Trait_final2$PVE))
  
  Outputfile3=paste(Pathway_out,"R2_",n, ".txt", sep="")
  write.table(Trait_final, file= Outputfile3, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
  
  Outputfile4=paste(Pathway_out,"GenomicRiskLoci_R2_",n, ".txt", sep="")
  write.table(Trait_final2, file= Outputfile4, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
}


#######Change
for (n in c("DHA")) {
  Trait_final_infile=paste(Pathway_out,"Sta_",n,".txt", sep="")
  Trait_final <- read.table(Trait_final_infile,header=T, as.is=T,sep = "\t")
  print(str(Trait_final))
  Trait_final$UKB_AF=sapply(strsplit(Trait_final$V10, split= ":", fixed=TRUE),"[",4)
  #Trait_final$UKB_AF=ifelse(Trait_final$Final_NS>115077,sapply(strsplit(Trait_final$V10, split= ":", fixed=TRUE),"[",4),NA)
  
  Trait_final$Final_EAF=ifelse(Trait_final$Final_NS>115077,Trait_final$UKB_AF,Trait_final$MAF)
  Trait_final$Final_EAF=as.numeric(Trait_final$Final_EAF)
  print(str(Trait_final))
  
  Trait_final$Final_R2=2*(Trait_final$Effect)*(Trait_final$Effect)*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)
  
  #######Change
  a=0.0828859
  
  Trait_final$Final_R2_SD=2*(Trait_final$Effect)*(Trait_final$Effect)*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)/(a*a)
  
  b=2*(Trait_final$Effect)*(Trait_final$Effect)*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)+(Trait_final$StdErr)*(Trait_final$StdErr)*2*Trait_final$Final_NS*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)
  Trait_final$PVE=2*(Trait_final$Effect)*(Trait_final$Effect)*(Trait_final$Final_EAF)*(1-Trait_final$Final_EAF)/b
  
  print(n)
  print(sum(Trait_final$Final_R2))
  print(sum(Trait_final$Final_R2_SD))
  print(sum(Trait_final$PVE))
  
  names(Trait_final)[names(Trait_final) == "SNP"] <- "rsID"
  
  Outputfile2=paste(Pathway_out,"GenomicRiskLoci_",n, ".txt", sep="")
  Trait_final2 <- read.table(Outputfile2,header=T, as.is=T,sep = "\t")
  
  Trait_final2<- Trait_final %>% right_join(Trait_final2, by= "rsID")
  print(str(Trait_final2))
  
  print(sum(Trait_final2$Final_R2))
  print(sum(Trait_final2$Final_R2_SD))
  print(sum(Trait_final2$PVE))
  
  Outputfile3=paste(Pathway_out,"R2_",n, ".txt", sep="")
  write.table(Trait_final, file= Outputfile3, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
  
  Outputfile4=paste(Pathway_out,"GenomicRiskLoci_R2_",n, ".txt", sep="")
  write.table(Trait_final2, file= Outputfile4, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
}

PUFA_LIST=c("PUFA","Omega_3","Omega_6","LA","DHA")
#PUFA 0.0616931 0.09575667 
#Omega_3 0.05280874 1.095328 
#Omega_6 0.05295491 0.1137094 
#LA 0.04952331 0.1052528 
#DHA 0.03775332 5.49533 

#r2=2*Beta^2*EAF(1-EAF)
#PUFA 0.0616931 
#Omega_3 0.05280874 
#Omega_6 0.05295491 
#LA 0.04952331  
#DHA 0.03775332 

#r2=2*Beta^2*EAF(1-EAF)/SD
#PUFA  0.09575667 
#Omega_3  1.095328 
#Omega_6  0.1137094 
#LA  0.1052528 
#DHA  5.49533 


#r2=2*Beta^2*EAF(1-EAF)/(2*Beta^2*EAF(1-EAF)+SE^2*2*N*EAF(1-EAF))
#PUFA  0.0655968 
#Omega_3  0.0544071 
#Omega_6  0.0559526 
#LA  0.05189069 
#DHA  0.0404109 
