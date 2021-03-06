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

PUFA_LIST=c("PUFA","Omega_3","Omega_6","LA","DHA")

for (n in PUFA_LIST) {
  Trait1_final_infile=paste(Pathway_geno,"met-d-",n,".vcf", sep="")
  
  Trait1_final <- read.table(Trait1_final_infile,header=F, as.is=T,sep = "\t")
  print(table(Trait1_final$V7))
  
  
  Trait1_final$LP_TYPE=sapply(strsplit(Trait1_final$V9, split= ":", fixed=TRUE),"[",3)
  print(table(Trait1_final$LP_TYPE))
  
  Trait1_final$lp=sapply(strsplit(Trait1_final$V10, split= ":", fixed=TRUE),"[",3)
  Trait1_final$lp=as.numeric(Trait1_final$lp)
  Trait1_final$p=10^-Trait1_final$lp
  
  Trait1_final$BETA=sapply(strsplit(Trait1_final$V10, split= ":", fixed=TRUE),"[",1)
  
  Trait1_final$SEBETA=sapply(strsplit(Trait1_final$V10, split= ":", fixed=TRUE),"[",2)
  
  Trait1_final$V3=ifelse(Trait1_final$V3==".",paste(Trait1_final$V1,":",Trait1_final$V2,":",Trait1_final$V4,"/",Trait1_final$V5,sep = ""),Trait1_final$V3)
  
  Trait1_final$NS=115078
  
  
  print(str(Trait1_final))
  print(summary(Trait1_final$lp))
  print(min(Trait1_final$p))
  
  Outputfile=paste(Pathway_out,"UKB_",n, ".txt", sep="")
  write.table(Trait1_final, file= Outputfile, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
  
  # Outputfile2=paste(Pathway_out,"met-d-",n,".jpeg", sep="")
  # #pdf(file=Outputfile2, paper="a4r",width=0, height=0)
  # jpeg(file=Outputfile2,width=1200, height=700,type="cairo")
  # #manhattan(Trait1_final,chr='V1',bp="V2",snp='V3',p='p', main = n,col = c("grey", "skyblue"), ylim = c(0, 200))
  # manhattan(Trait1_final,chr='V1',bp="V2",snp='V3',p='p', main = n,col = c("grey", "skyblue"))
  # #manp <- manhattan(Trait1_final,chr='V1',bp="V2",snp='V3',p='p',annotatePval = 0.0000001, main = n,col = c("brewer.set3"))
  # 
  # #manp <-manhattan(Trait1_final,chr='V1',bp="V2",snp='V3',p='p', main = n,col = c("grey", "skyblue"))
  # #Outputfile1=paste(Pathway_out,"met-d-",n,"_ggsave.pdf", sep="")
  # #ggsave(manp, file=Outputfile1)
  # dev.off()
  # 
  # Outputfile3=paste(Pathway_out,"met-d-",n,".png", sep="")
  # png(file=Outputfile3,
  #     width=1200, height=700,type="cairo")
  # #manhattan(Trait1_final,chr='V1',bp="V2",snp='V3',p='p', main = n,col = c("grey", "skyblue"), ylim = c(0, 200))
  # manhattan(Trait1_final,chr='V1',bp="V2",snp='V3',p='p', main = n,col = c("grey", "skyblue"))
  # dev.off()
}

for (n in PUFA_LIST) {
  Trait1_final_infile=paste(Pathway_geno,"LockeAE_",n,".txt", sep="")
  
  Trait1_final <- read.table(Trait1_final_infile,header=T, as.is=T,sep = "\t")
  print(summary(Trait1_final$CALLRATE))
  
  Trait1_final$SNP=sapply(strsplit(Trait1_final$MARKER_ID, split= "_", fixed=TRUE),"[",3)
  Trait1_final$SNP1=sapply(strsplit(Trait1_final$MARKER_ID, split= "_", fixed=TRUE),"[",4)
  Trait1_final$SNP=ifelse(is.na(Trait1_final$SNP1)==F,paste(Trait1_final$SNP,Trait1_final$SNP1,sep = ":"),Trait1_final$SNP)
  Trait1_final$SNP2=sapply(strsplit(Trait1_final$MARKER_ID, split= "_", fixed=TRUE),"[",5)
  Trait1_final$SNP=ifelse(is.na(Trait1_final$SNP2)==F,paste(Trait1_final$SNP,Trait1_final$SNP2,sep = ":"),Trait1_final$SNP)
  Trait1_final$SNP3=sapply(strsplit(Trait1_final$MARKER_ID, split= "_", fixed=TRUE),"[",6)
  Trait1_final$SNP=ifelse(is.na(Trait1_final$SNP3)==F,paste(Trait1_final$SNP,Trait1_final$SNP3,sep = "/"),Trait1_final$SNP)
  
  Trait1_final$CHROM[Trait1_final$CHROM == "X"] <- 23
  print(table(Trait1_final$CHROM))
  Trait1_final$CHROM=as.numeric(Trait1_final$CHROM)
  
  Trait1_final$ALLELE=sapply(strsplit(Trait1_final$MARKER_ID, split= "_", fixed=TRUE),"[",2)
  Trait1_final$NonRefAllele=sapply(strsplit(Trait1_final$ALLELE, split= "/", fixed=TRUE),"[",2)
  Trait1_final$RefAllele=sapply(strsplit(Trait1_final$ALLELE, split= "/", fixed=TRUE),"[",1)
  
  
  print(table(sapply(strsplit(Trait1_final$MARKER_ID, split= "_", fixed=TRUE),"[",7)))
  
  print(summary(Trait1_final$END-Trait1_final$BEG))
  
  print(str(Trait1_final))
  
  print(min(Trait1_final$PVALUE))
  
  Outputfile=paste(Pathway_out,"Locke_",n, ".txt", sep="")
  write.table(Trait1_final, file= Outputfile, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
  
  # Outputfile2=paste(Pathway_out,"met-d-",n,".jpeg", sep="")
  # #pdf(file=Outputfile2, paper="a4r",width=0, height=0)
  # jpeg(file=Outputfile2,width=1200, height=700,type="cairo")
  # #manhattan(Trait1_final,chr='V1',bp="V2",snp='V3',p='p', main = n,col = c("grey", "skyblue"), ylim = c(0, 200))
  # manhattan(Trait1_final,chr='V1',bp="V2",snp='V3',p='p', main = n,col = c("grey", "skyblue"))
  # #manp <- manhattan(Trait1_final,chr='V1',bp="V2",snp='V3',p='p',annotatePval = 0.0000001, main = n,col = c("brewer.set3"))
  # 
  # #manp <-manhattan(Trait1_final,chr='V1',bp="V2",snp='V3',p='p', main = n,col = c("grey", "skyblue"))
  # #Outputfile1=paste(Pathway_out,"met-d-",n,"_ggsave.pdf", sep="")
  # #ggsave(manp, file=Outputfile1)
  # dev.off()
  # 
  # Outputfile3=paste(Pathway_out,"LockeAE_",n,".png", sep="")
  # png(file=Outputfile3,
  #     width=1200, height=700,type="cairo")
  # #manhattan(Trait1_final,chr='V1',bp="V2",snp='V3',p='p', main = n,col = c("grey", "skyblue"), ylim = c(0, 200))
  # manhattan(Trait1_final,chr='CHROM',bp="BEG",snp='SNP',p='PVALUE', main = n,col = c("grey", "skyblue"))
  # dev.off()
}
