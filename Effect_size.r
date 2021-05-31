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


library(genpwr)

#Pathway_out="/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/MR_result/result_GWAS_05_19/"

pw <- genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
                  N=400, Case.Rate=0.2, k=NULL,
                  MAF=seq(0.18, 0.25, 0.01), OR=c(3),Alpha=0.05,
                  True.Model=c("Dominant", "Recessive", "Additive"), 
                  Test.Model=c("Dominant", "Recessive", "Additive", "2df"))

#PUFA 0.802664 Omega_3 0.219574 Omega_6 0.682425 LA 0.685943 DHA 0.0828859
#LA 0.685943 DHA 0.0828859
#LA 0.685943 
es1 <- es.calc.linear(N=13544,power=0.8,
                     MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.685943, Alpha=5e-8,
                     True.Model='Additive', Test.Model='Additive')
es1
es2 <- es.calc.linear(N=19037,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.685943, Alpha=5e-8,
                      True.Model='Additive', Test.Model='Additive')
es2
es3 <- es.calc.linear(N=22295,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.685943, Alpha=5e-8,
                      True.Model='Additive', Test.Model='Additive')
es3
es4 <- es.calc.linear(N=19037,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.685943, Alpha=0.05,
                      True.Model='Additive', Test.Model='Additive')
es4
es5 <- es.calc.linear(N=22295,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.685943, Alpha=0.05,
                      True.Model='Additive', Test.Model='Additive')
es5
es6 <- es.calc.linear(N=108099,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.685943, Alpha=5e-8,
                      True.Model='Additive', Test.Model='Additive')
es6

library(tidyverse)
library(ggrepel)
# Set ggplot2 default theme to theme_bw()
theme_set(theme_bw())

df <- tibble::tribble(
  ~MAF, ~Group_1, ~Group_2, ~Group_3, ~Group_4,~Group_5,~Group_6,
  "0.5%",         0.37157,        0.31348,         0.28969,        0.13961,   0.12901,        0.13161,
  "1%",          0.2634,        0.22222,         0.20536,         0.09895,    0.09146,         0.09329,
  "5%",         0.12024,        0.10143,         0.09378,        0.04517,        0.04176,        0.0426,
  "10%",         0.08735,        0.07368,         0.06808,        0.03283,        0.03034,        0.03095,
  "15%",          0.07338,        0.06192,         0.05722,         0.02756,        0.02548,         0.02599,
  "20%",         0.06552,        0.05527,         0.05108,        0.02461,        0.02275,        0.02321,
  "25%",         0.06052,        0.05106,         0.04718,        0.02274,      0.02102,        0.02144,
  "30%",          0.05719,        0.04824,         0.04457,         0.02149,         0.01986,         0.02026,
  "35%",         0.05494,        0.04635,         0.04285,        0.02065,    0.01908,        0.01947,
  "40%",         0.05349,        0.04512,         0.04172,        0.0201,      0.01858,        0.01895,
  "45%",         0.05267,        0.04443,         0.04108,        0.0198,     0.0183,        0.01866,
  "50%",          0.05241,        0.04421,         0.04088,         0.0197,        0.0182,         0.01857
)
df

df_long <- df %>%
  pivot_longer(
    Group_1:Group_2:Group_3:Group_4:Group_5:Group_6,
    names_to = "Group", values_to = "MES"
  )
df_long

# df_long$MAF <- factor(df_long$MAF, levels = df_long$MAF)
# df_long$MAF <- factor(df_long$MAF, levels =df_long$MAF[order(df_long$MAF)])

df_long <- df_long %>% mutate(row = row_number())

Pathway_out="/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/MR_result/result_GWAS_05_19/"
Outputfile3=paste(Pathway_out,"LA_6_ES.png", sep="")
#png(file=Outputfile3,width=1200, height=700,type="cairo")

png(file=Outputfile3,width     = 3.25,
    height    = 3.25,
    units     = "in",
    res       = 1200,
    pointsize = 4,type="cairo")
lp <- ggplot(df_long, aes(x=reorder(MAF,row), y = MES, group = Group)) +  ylab("Minimum Effect Size")+
  geom_line(aes(color = Group), lwd=1.5) +
  geom_point() +theme(axis.text.x = element_text(colour = 'black', angle = 90, size = 10, hjust = 0.5, vjust = 0.5),axis.title.x=element_blank()) +
  theme(legend.position = c(0.8,0.79), axis.text=element_text(size=10),
        legend.title=element_text(size=0), 
        legend.text=element_text(size=10))

# Filter the last values and add onto the line plot
# Corresponds to the `virginica` species
data_ends <- df_long %>% filter(MAF == "virginica")
lp + 
  geom_text_repel(
    aes(label = MES), data = data_ends,
    fontface ="plain", color = "black", size = 3
  )
dev.off()


#PUFA 0.802664 
es1 <- es.calc.linear(N=13544,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.802664, Alpha=5e-8,
                      True.Model='Additive', Test.Model='Additive')
es1
es2 <- es.calc.linear(N=19037,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.802664, Alpha=5e-8,
                      True.Model='Additive', Test.Model='Additive')
es2
es3 <- es.calc.linear(N=22295,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.802664, Alpha=5e-8,
                      True.Model='Additive', Test.Model='Additive')
es3
es4 <- es.calc.linear(N=19037,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.802664, Alpha=0.05,
                      True.Model='Additive', Test.Model='Additive')
es4
es5 <- es.calc.linear(N=22295,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.802664, Alpha=0.05,
                      True.Model='Additive', Test.Model='Additive')
es5
es6 <- es.calc.linear(N=108099,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.802664, Alpha=5e-8,
                      True.Model='Additive', Test.Model='Additive')
es6

library(tidyverse)
library(ggrepel)
# Set ggplot2 default theme to theme_bw()
theme_set(theme_bw())

df <- tibble::tribble(
  ~MAF, ~Group_1, ~Group_2, ~Group_3, ~Group_4,~Group_5,~Group_6,
  "0.5%",         0.43479,        0.36681,         0.33898,        0.16337,   0.15097,        0.154,
  "1%",          0.30822,        0.26003,         0.2403,         0.11581,    0.10702,         0.10917,
  "5%",         0.1407,        0.11873,         0.10971,        0.05289,        0.04888,        0.04986,
  "10%",         0.1022,        0.08624,         0.07973,        0.03843,        0.03552,        0.03623,
  "15%",          0.08588,        0.07244,         0.06693,         0.03228,        0.02983,         0.03043,
  "20%",         0.07665,        0.0647,         0.0598,        0.02882,        0.02663,        0.02717,
  "25%",         0.07081,        0.05975,         0.05521,        0.0266,      0.02459,        0.02508,
  "30%",          0.06692,        0.05646,         0.05217,         0.02514,         0.02323,         0.0237,
  "35%",         0.0643,        0.05424,         0.05012,        0.02416,    0.02233,        0.02278,
  "40%",         0.0626,        0.05281,         0.0488,        0.02352,      0.02174,        0.02217,
  "45%",         0.06164,        0.052,         0.04805,        0.02316,     0.02141,        0.02184,
  "50%",          0.06133,        0.05174,         0.04781,         0.02304,        0.0213,         0.02173
)
df

df_long <- df %>%
  pivot_longer(
    Group_1:Group_2:Group_3:Group_4:Group_5:Group_6,
    names_to = "Group", values_to = "MES"
  )
df_long

# df_long$MAF <- factor(df_long$MAF, levels = df_long$MAF)
# df_long$MAF <- factor(df_long$MAF, levels =df_long$MAF[order(df_long$MAF)])

df_long <- df_long %>% mutate(row = row_number())

Pathway_out="/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/MR_result/result_GWAS_05_19/"
Outputfile3=paste(Pathway_out,"PUFA_ES.png", sep="")
#png(file=Outputfile3,width=1200, height=700,type="cairo")

png(file=Outputfile3,width     = 3.25,
    height    = 3.25,
    units     = "in",
    res       = 1200,
    pointsize = 4,type="cairo")
lp <- ggplot(df_long, aes(x=reorder(MAF,row), y = MES, group = Group)) +  ylab("Minimum Effect Size")+
  geom_line(aes(color = Group), lwd=1.5) +
  geom_point() +theme(axis.text.x = element_text(colour = 'black', angle = 90, size = 10, hjust = 0.5, vjust = 0.5),axis.title.x=element_blank()) +
  theme(legend.position = c(0.8,0.79), axis.text=element_text(size=10),
        legend.title=element_text(size=0), 
        legend.text=element_text(size=10))

# Filter the last values and add onto the line plot
# Corresponds to the `virginica` species
data_ends <- df_long %>% filter(MAF == "virginica")
lp + 
  geom_text_repel(
    aes(label = MES), data = data_ends,
    fontface ="plain", color = "black"
    #, size = 3
  )
dev.off()


#Omega_3 0.219574 
es1 <- es.calc.linear(N=13544,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.219574, Alpha=5e-8,
                      True.Model='Additive', Test.Model='Additive')
es1
es2 <- es.calc.linear(N=19037,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.219574, Alpha=5e-8,
                      True.Model='Additive', Test.Model='Additive')
es2
es3 <- es.calc.linear(N=22295,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.219574, Alpha=5e-8,
                      True.Model='Additive', Test.Model='Additive')
es3
es4 <- es.calc.linear(N=19037,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.219574, Alpha=0.05,
                      True.Model='Additive', Test.Model='Additive')
es4
es5 <- es.calc.linear(N=22295,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.219574, Alpha=0.05,
                      True.Model='Additive', Test.Model='Additive')
es5
es6 <- es.calc.linear(N=108099,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.219574, Alpha=5e-8,
                      True.Model='Additive', Test.Model='Additive')
es6

library(tidyverse)
library(ggrepel)
# Set ggplot2 default theme to theme_bw()
theme_set(theme_bw())

df <- tibble::tribble(
  ~MAF, ~Group_1, ~Group_2, ~Group_3, ~Group_4,~Group_5,~Group_6,
  "0.5%",         0.11893,        0.10032,         0.09275,        0.04468,   0.04131,        0.04214,
  "1%",          0.08431,        0.07111,         0.06576,         0.03169,    0.02928,         0.02987,
  "5%",         0.0385,        0.0325,         0.03004,        0.01447,        0.01337,        0.01364,
  "10%",         0.02794,        0.02359,         0.0218,        0.01052,        NA,        NA,
  "15%",          -0.02349,        -0.01987,         -0.01831,         NA,        NA,         NA,
  "20%",         -0.02097,        -0.01767,         -0.01635,        NA,        NA,        NA,
  "25%",         -0.01937,        -0.01634,         -0.01509,        NA,      NA,        NA,
  "30%",          -0.0183,        -0.01543,         -0.0143,         NA,         NA,         NA,
  "35%",         -0.01756,        -0.01482,         -0.01371,        NA,    NA,        NA,
  "40%",         -0.01714,        -0.01448,         -0.01335,        NA,      NA,        NA,
  "45%",         -0.01686,        -0.01425,         -0.01315,        NA,     NA,        NA,
  "50%",          -0.01678,        -0.01417,         -0.01308,         NA,        NA,         NA
)
df

df_long <- df %>%
  pivot_longer(
    Group_1:Group_2:Group_3:Group_4:Group_5:Group_6,
    names_to = "Group", values_to = "MES"
  )
df_long

# df_long$MAF <- factor(df_long$MAF, levels = df_long$MAF)
# df_long$MAF <- factor(df_long$MAF, levels =df_long$MAF[order(df_long$MAF)])

df_long <- df_long %>% mutate(row = row_number())

Pathway_out="/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/MR_result/result_GWAS_05_19/"
Outputfile3=paste(Pathway_out,"Omega_3_ES.png", sep="")
#png(file=Outputfile3,width=1200, height=700,type="cairo")

png(file=Outputfile3,width     = 3.25,
    height    = 3.25,
    units     = "in",
    res       = 1200,
    pointsize = 4,type="cairo")
lp <- ggplot(df_long, aes(x=reorder(MAF,row), y = MES, group = Group)) +  ylab("Minimum Effect Size")+
  geom_line(aes(color = Group), lwd=1.5) +
  geom_point() +theme(axis.text.x = element_text(colour = 'black', angle = 90, size = 10, hjust = 0.5, vjust = 0.5),axis.title.x=element_blank()) +
  theme(legend.position = c(0.8,0.79), axis.text=element_text(size=10),
        legend.title=element_text(size=0), 
        legend.text=element_text(size=10))

# Filter the last values and add onto the line plot
# Corresponds to the `virginica` species
data_ends <- df_long %>% filter(MAF == "virginica")
lp + 
  geom_text_repel(
    aes(label = MES), data = data_ends,
    fontface ="plain", color = "black", size = 3
  )
dev.off()


# DHA 0.0828859
es1 <- es.calc.linear(N=13544,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.0828859, Alpha=5e-8,
                      True.Model='Additive', Test.Model='Additive')
es1
es2 <- es.calc.linear(N=19037,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.0828859, Alpha=5e-8,
                      True.Model='Additive', Test.Model='Additive')
es2
es3 <- es.calc.linear(N=22295,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.0828859, Alpha=5e-8,
                      True.Model='Additive', Test.Model='Additive')
es3
es4 <- es.calc.linear(N=19037,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.0828859, Alpha=0.05,
                      True.Model='Additive', Test.Model='Additive')
es4
es5 <- es.calc.linear(N=22295,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.0828859, Alpha=0.05,
                      True.Model='Additive', Test.Model='Additive')
es5
es6 <- es.calc.linear(N=108099,power=0.8,
                      MAF=c(0.005,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), sd_y = 0.0828859, Alpha=5e-8,
                      True.Model='Additive', Test.Model='Additive')
es6

library(tidyverse)
library(ggrepel)
# Set ggplot2 default theme to theme_bw()
theme_set(theme_bw())

df <- tibble::tribble(
  ~MAF, ~Group_1, ~Group_2, ~Group_3, ~Group_4,~Group_5,~Group_6,
  "0.5%",         0.04489,        0.03789,         0.03502,        0.01687,   0.01559,        0.01591,
  "1%",          0.03185,        0.02683,         0.0248,         0.01199,    0.01103,         0.01125,
  "5%",         -0.0145,        -0.01228,         -0.01133,        NA,        NA,        NA,
  "10%",         -0.01055,        NA,         NA,        NA,      NA,        NA,
)
df

df_long <- df %>%
  pivot_longer(
    Group_1:Group_2:Group_3:Group_4:Group_5:Group_6,
    names_to = "Group", values_to = "MES"
  )
df_long

# df_long$MAF <- factor(df_long$MAF, levels = df_long$MAF)
# df_long$MAF <- factor(df_long$MAF, levels =df_long$MAF[order(df_long$MAF)])

df_long <- df_long %>% mutate(row = row_number())

Pathway_out="/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/MR_result/result_GWAS_05_19/"
Outputfile3=paste(Pathway_out,"DHA_ES.png", sep="")
#png(file=Outputfile3,width=1200, height=700,type="cairo")

png(file=Outputfile3,type="cairo")
lp <- ggplot(df_long, aes(x=reorder(MAF,row), y = MES, group = Group)) +
  geom_line(aes(color = Group)) +
  geom_point() +
  theme(legend.position = "top")

# Filter the last values and add onto the line plot
# Corresponds to the `virginica` species
data_ends <- df_long %>% filter(MAF == "virginica")
lp + 
  geom_text_repel(
    aes(label = MES), data = data_ends,
    fontface ="plain", color = "black", size = 3
  )
dev.off()



Outputfile3=paste(Pathway_out,"sd_y_1",".txt", sep="")
write.table(es, file= Outputfile3, col.names = T, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')



ss.calc.linear(
  power = 0.8,
  MAF = c(0.005,0.01,0.05,0.1,0.2,0.3,0.4),
  ES = 25.69125,
  
  sd_y = 1,
  Alpha = 0.05,
  True.Model = "All",
  Test.Model = "All"
)


pw <- power_envir.calc.linear_outcome(N=100, ES_G = 1.2, ES_E = 1.3,
                                      ES_GE = 2, Alpha = 0.05, MAF = 0.2, P_e = 0.2,
                                      sd_y = 10, True.Model = "All", Test.Model = "All")

Outputfile3=paste(Pathway_out,"Final__",n,".png", sep="")
png(file=Outputfile3,
    width=1200, height=700,type="cairo")
#manhattan(Trait1_final,chr='V1',bp="V2",snp='V3',p='p', main = n,col = c("grey", "skyblue"), ylim = c(0, 200))
manhattan(Trait_final,chr='Final_chr',bp="Final_pos",snp='SNP',p='P.value', main = n,col = c("grey", "skyblue"))
dev.off()
