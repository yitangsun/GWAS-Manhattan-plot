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
  Trait_final_infile=paste(Pathway_out,"META_",n,"1.txt", sep="")
  Trait_final <- read.table(Trait_final_infile,header=T, as.is=T,sep = "\t")
  
  names(Trait_final)[names(Trait_final) == "MarkerName"] <- "SNP"
  
  Outputfile2=paste(Pathway_out,"UKB_",n, ".txt", sep="")
  Trait_final2 <- read.table(Outputfile2,header=T, as.is=T,sep = "\t")
  
  names(Trait_final2)[names(Trait_final2) == "V3"] <- "SNP"
  
  Trait_final<- Trait_final %>% left_join(Trait_final2, by= "SNP")
  
  Outputfile1=paste(Pathway_out,"Locke_",n, ".txt", sep="")
  Trait_final1 <- read.table(Outputfile1,header=T, as.is=T,sep = "\t")
  
  Trait_final<- Trait_final %>% left_join(Trait_final1, by= "SNP")
  
  Trait_final$Final_chr=ifelse(Trait_final$Weight>115077,Trait_final$V1,Trait_final$CHROM)
  Trait_final$Final_pos=ifelse(Trait_final$Weight>115077,Trait_final$V2,Trait_final$BEG)
  
  Outputfile=paste(Pathway_out,"Final__",n, ".txt", sep="")
  write.table(Trait_final, file= Outputfile, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
  
  Outputfile3=paste(Pathway_out,"Final__",n,".png", sep="")
  png(file=Outputfile3,
      width=1200, height=700,type="cairo")
  #manhattan(Trait1_final,chr='V1',bp="V2",snp='V3',p='p', main = n,col = c("grey", "skyblue"), ylim = c(0, 200))
  manhattan(Trait_final,chr='Final_chr',bp="Final_pos",snp='SNP',p='P.value', main = n,col = c("grey", "skyblue"))
  dev.off()
  
  Outputfile4=paste(Pathway_out,"Clear_",n,".pdf", sep="")
  pdf(file=Outputfile4,
      width=1200, height=700)
  #manhattan(Trait1_final,chr='V1',bp="V2",snp='V3',p='p', main = n,col = c("grey", "skyblue"), ylim = c(0, 200))
  manhattan(Trait_final,chr='Final_chr',bp="Final_pos",snp='SNP',p='P.value', main = n,col = c("grey", "skyblue"))
  dev.off()
}

for (n in PUFA_LIST) {
  Outputfile=paste(Pathway_out,"Final__",n, ".txt", sep="")
  final_t <- read.table(Outputfile,header=T, as.is=T,sep = "\t")
  print(str(final_t))
  Final_1=final_t[final_t$Weight<115078,]
  
  print(str(Final_1))
  
  Final_2=final_t[final_t$P.value <5e-8,]
  
  print(str(Final_2))
  
  Outputfile1=paste(Pathway_out,"Dif_",n, ".txt", sep="")
  write.table(Final_1, file= Outputfile1, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
}


for (n in PUFA_LIST) {
  Trait_final_infile=paste(Pathway_out,"META_IVW_",n,"1.txt", sep="")
  Trait_final <- read.table(Trait_final_infile,header=T, as.is=T,sep = "\t")
  
  names(Trait_final)[names(Trait_final) == "MarkerName"] <- "SNP"
  names(Trait_final)[names(Trait_final) == "P.value"] <- "Trait_pvalue"
  
  Outputfile2=paste(Pathway_out,"UKB_",n, ".txt", sep="")
  Trait_final2 <- read.table(Outputfile2,header=T, as.is=T,sep = "\t")
  
  names(Trait_final2)[names(Trait_final2) == "V3"] <- "SNP"
  names(Trait_final2)[names(Trait_final2) == "NS"] <- "NS1"
  
  Trait_final<- Trait_final %>% left_join(Trait_final2, by= "SNP")
  
  Outputfile1=paste(Pathway_out,"Locke_",n, ".txt", sep="")
  Trait_final1 <- read.table(Outputfile1,header=T, as.is=T,sep = "\t")
  names(Trait_final1)[names(Trait_final1) == "NS"] <- "NS2"
  
  Trait_final<- Trait_final %>% left_join(Trait_final1, by= "SNP")
  
  print(str(Trait_final))
  
  # Trait_final$NS1[which(Trait_final$NS1 == "NULL")] <- 0 
  # Trait_final$NS2[which(Trait_final$NS2 == "NULL")] <- 0 
  Trait_final$NS1[is.na(Trait_final$NS1 == T)] <- 0 
  Trait_final$NS2[is.na(Trait_final$NS2 == T)] <- 0 
  #Trait_final1$NS2[which(Trait_final1$NS2 == "NULL")] <- 0 
  #Trait_final1[is.na(Trait_final1$Final_NS2 == T)]$Final_NS2 <- 0 
  Trait_final$Final_NS=Trait_final$NS1+Trait_final$NS2
  
  print(str(Trait_final))
  
  Trait_final$Final_chr=ifelse(Trait_final$Final_NS>115077,Trait_final$V1,Trait_final$CHROM)
  Trait_final$Final_pos=ifelse(Trait_final$Final_NS>115077,Trait_final$V2,Trait_final$BEG)
  
  Outputfile=paste(Pathway_out,"Sta_",n, ".txt", sep="")
  write.table(Trait_final, file= Outputfile, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
  
  Trait_final<-Trait_final%>%select(Final_chr,  Final_pos,SNP,Trait_pvalue,,Allele1,,Allele2,Effect,StdErr,Final_NS)
  
  Outputfile3=paste(Pathway_out,"Sum_sta_",n,".txt", sep="")
  write.table(Trait_final, file= Outputfile3, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
}

#Pathway_out="/Users/yelab/GWAS/result_GWAS_05_19/"

for (n in PUFA_LIST) {
  Trait_final_infile=paste(Pathway_out,"Sum_sta_",n,".txt.gz", sep="")
  Trait_final <- read.table(Trait_final_infile,header=T, as.is=T,sep = "\t")
  
  Trait_final$Final_chr[which(Trait_final$Final_chr == 23)] <- "X"
  
  Outputfile3=paste(Pathway_out,"Sum_sta_X_",n,".txt.gz", sep="")
  write.table(Trait_final, file= Outputfile3, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
}
 
                 
