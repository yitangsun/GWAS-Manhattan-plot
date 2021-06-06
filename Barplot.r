Pathway_out="/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/MR_result/result_GWAS_05_19/"
Outputfile3=paste(Pathway_out,"G_E_5e8_ES.png", sep="")
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

########
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

######## Figure 1
cat1=c("GAPDH-forward","GAPDH-forward","GAPDH-reverse","GAPDH-reverse","pp179-forward","pp179-forward","pp179-reverse","pp179-reverse","Empty","Empty")
cat2=c("Firefly","Renilla","Firefly","Renilla","Firefly","Renilla","Firefly","Renilla","Firefly","Renilla")
mean(c(8590
     ,5190
     ,9680))
mean(c(2472
     ,1522.5
     ,1777))
mean(c(561.6
     ,666.5555556
     ,841.3))
mean(c(642.2
     ,793.5555556
     ,950.5))
mean(c(2172.857143
     ,2631.428571
     ,2425.333333))
mean(c(555.2142857
     ,710
     ,618.7333333))
ratio=c(7820,0,0,1923.833,689.8185,795.4185,2409.873,627.9825,0,0)

sd(c(8590
       ,5190
       ,9680))
sd(c(2472
       ,1522.5
       ,1777))
sd(c(561.6
       ,666.5555556
       ,841.3))
sd(c(642.2
       ,793.5555556
       ,950.5))
sd(c(2172.857143
       ,2631.428571
       ,2425.333333))
sd(c(555.2142857
       ,710
       ,618.7333333))
sd=c(2341.944,0,0,491.4851,141.2937,154.1584,229.6763,77.80627,0,0)

df2=data.frame(cat1,cat2,ratio,sd)

Pathway_out="Cooperator/Review-fads -function/Barplot/"
Outputfile3=paste(Pathway_out,"Figure1.png", sep="")
#png(file=Outputfile3,width=1200, height=700,type="cairo")

png(file=Outputfile3,width     = 3,
    height    = 2.1,
    units     = "in",
    res       = 1200,
    pointsize = 4,type="cairo")

p <- ggplot(df2, aes(x=ratio, y=cat1, fill=cat2)) + xlab("Relative Luciferase Activity (RLU)")+
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(xmin=ratio-sd, xmax=ratio+sd), width=.2,
                position=position_dodge(.9))

#values=c('#E69F00','#56B4E9')
p +
  scale_fill_manual(values=c('#F8CBAD','#BDD7EE'))+
  theme_classic()+
  theme(legend.position = c(0.9,0.9), axis.text=element_text(size=8.5),
        legend.title=element_text(size=0), axis.title = element_text(size = 8.5),  
        legend.text=element_text(size=8.5),axis.title.y=element_blank())

dev.off()


######## Figure 2
cat1=c("649bp","649bp","98bp","98bp","Empty","Empty")
cat2=c("Haplotype A","Haplotype D","Haplotype A","Haplotype D","Haplotype A","Haplotype D")
mean(c(10.24153445
       ,11.91242756
       ,9.986919555
       ,12.08362863
       ,10.69320066
       ,11.48018239))
mean(c(16.07028754
       ,15.2276176
       ,14.49214366
       ,15.14016173
       ,11.86147186
       ,14.47518097))
mean(c(24.37321937
       ,25.97261465
       ,30.48351648
       ,26.11247611
       ,31.40260631
       ,27.2280874))
mean(c(48.54100946
       ,53.041988
       ,39.50251703
       ,48.798252
       ,51.05726872
       ,45.74162679))
mean(c(1.829339964
       ,1.885362015
       ,1.76038835
       ,2.587822014
       ,2.013876534
       ,2.223793887))
ratio=c(11.06632,14.54448,27.59542,47.78044,0,2.050097)

sd(c(10.24153445
       ,11.91242756
       ,9.986919555
       ,12.08362863
       ,10.69320066
       ,11.48018239))
sd(c(16.07028754
       ,15.2276176
       ,14.49214366
       ,15.14016173
       ,11.86147186
       ,14.47518097))
sd(c(24.37321937
       ,25.97261465
       ,30.48351648
       ,26.11247611
       ,31.40260631
       ,27.2280874))
sd(c(48.54100946
       ,53.041988
       ,39.50251703
       ,48.798252
       ,51.05726872
       ,45.74162679))
sd(c(1.829339964
       ,1.885362015
       ,1.76038835
       ,2.587822014
       ,2.013876534
       ,2.223793887))
sd=c(0.8839334,1.439144,2.763545,4.747982,0,0.3099784)

df2=data.frame(cat1,cat2,ratio,sd)
# df_long$MAF <- factor(df_long$MAF, levels =df_long$MAF[order(df_long$MAF)])

df2$cat1 <- factor(df2$cat1, levels = unique(df2$cat1))
df2$cat2 <- factor(df2$cat2, levels = unique(df2$cat2))

df2 <- df2[order(df2$cat1), ]

Pathway_out="Cooperator/Review-fads -function/Barplot/"
Outputfile3=paste(Pathway_out,"Figure2.png", sep="")
#png(file=Outputfile3,width=1200, height=700,type="cairo")

png(file=Outputfile3,width     = 3,
    height    = 2.1,
    units     = "in",
    res       = 1200,
    pointsize = 4,type="cairo")

p <- ggplot(df2, aes(x=ratio, y=cat1, fill=cat2)) + xlab("Relative Luciferase Activity (RLU)")+
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(xmin=ratio-sd, xmax=ratio+sd), width=.2,
                position=position_dodge(.9))

#values=c('#E69F00','#56B4E9')
p +
  scale_fill_manual(values=c('Black','Grey'))+
  theme_classic()+
  theme(legend.position = c(0.8,0.9), axis.text=element_text(size=8.5),
        legend.title=element_text(size=0), axis.title = element_text(size = 8.5),  
        legend.text=element_text(size=8.5),axis.title.y=element_blank())

dev.off()

######## Figure 3
cat1=c("GAPDH-forward","GAPDH-forward","GAPDH-reverse","GAPDH-reverse","pp179-forward","pp179-forward","pp179-reverse","pp179-reverse","Empty","Empty")
cat2=c("Firefly","Renilla","Firefly","Renilla","Firefly","Renilla","Firefly","Renilla","Firefly","Renilla")
mean(c(8590
       ,5190
       ,9680
       ,3968.115942
       ,3480
       ,4190.140845))
mean(c(0
       ,0
       ,0
       ,0
       ,0.147239264
       ,0.171428571))
mean(c(2472
       ,1522.5
       ,1777
       ,677.0072993
       ,699.3865031
       ,704.5714286))
mean(c(561.6
       ,666.5555556
       ,841.3
       ,510.7207207
       ,503.4188034
       ,651.0169492))
mean(c(642.2
       ,793.5555556
       ,950.5
       ,603.5135135
       ,569.4017094
       ,749.8305085))
mean(c(2172.857143
       ,2631.428571
       ,2425.333333
       ,1902.083333
       
       ,1700))
mean(c(555.2142857
       ,710
       ,618.7333333
       ,645.4166667
       
       ,575.1612903))
ratio=c(5849.709,0,0.05311131,1308.744,622.4353,718.1669,2166.34,620.9051,0,0)

sd(c(8590
       ,5190
       ,9680
       ,3968.115942
       ,3480
       ,4190.140845))
sd(c(0
       ,0
       ,0
       ,0
       ,0.147239264
       ,0.171428571))
sd(c(2472
       ,1522.5
       ,1777
       ,677.0072993
       ,699.3865031
       ,704.5714286))
sd(c(561.6
       ,666.5555556
       ,841.3
       ,510.7207207
       ,503.4188034
       ,651.0169492))
sd(c(642.2
       ,793.5555556
       ,950.5
       ,603.5135135
       ,569.4017094
       ,749.8305085))
sd(c(2172.857143
       ,2631.428571
       ,2425.333333
       ,1902.083333
       
       ,1700))
sd(c(555.2142857
       ,710
       ,618.7333333
       ,645.4166667
       
       ,575.1612903))
sd=c(2627.76,0,0.08263449,742.0987,127.2881,142.6269,377.734,61.13758,0,0)

df2=data.frame(cat1,cat2,ratio,sd)

Pathway_out="Cooperator/Review-fads -function/Barplot/"
Outputfile3=paste(Pathway_out,"Figure3.png", sep="")
#png(file=Outputfile3,width=1200, height=700,type="cairo")

png(file=Outputfile3,width     = 3,
    height    = 2.1,
    units     = "in",
    res       = 1200,
    pointsize = 4,type="cairo")

p <- ggplot(df2, aes(x=ratio, y=cat1, fill=cat2)) + xlab("Relative Luciferase Activity (RLU)")+
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(xmin=ratio-sd, xmax=ratio+sd), width=.2,
                position=position_dodge(.9))

#values=c('#E69F00','#56B4E9')
p +
  scale_fill_manual(values=c('#F8CBAD','#BDD7EE'))+
  theme_classic()+
  theme(legend.position = c(0.9,0.9), axis.text=element_text(size=8.5),
        legend.title=element_text(size=0), axis.title = element_text(size = 8.5),  
        legend.text=element_text(size=8.5),axis.title.y=element_blank())

dev.off()





lp <- ggplot(df_long, aes(x=reorder(MAF,row), y = MES, group = Group)) +  ylab("Minimum Effect Size")+
  geom_line(aes(color = Group), lwd=1.5) +
  geom_point() +theme(axis.text.x = element_text(colour = 'black', angle = 90, size = 10, hjust = 0.5, vjust = 0.5),axis.title.x=element_blank()) +
  theme(legend.position = c(0.8,0.79), axis.text=element_text(size=10),
        legend.title=element_text(size=0), 
        legend.text=element_text(size=10),)

# Filter the last values and add onto the line plot
# Corresponds to the `virginica` species
data_ends <- df_long %>% filter(MAF == "virginica")
lp + 
  geom_text_repel(
    aes(label = MES), data = data_ends,
    fontface ="plain", color = "black", size = 3
  )


p + labs(title="Plot of length  per dose", 
         x="Dose (mg)", y = "Length")+
  scale_fill_manual(values=c('#E69F00','#56B4E9'))+
  theme_classic()
