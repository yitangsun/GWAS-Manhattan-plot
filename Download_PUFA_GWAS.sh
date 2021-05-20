cd /scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/Genotype/

#### UKB May 19 2021
wget https://gwas.mrcieu.ac.uk/files/met-d-Omega_3/met-d-Omega_3.vcf.gz

wget https://gwas.mrcieu.ac.uk/files/met-d-Omega_6/met-d-Omega_6.vcf.gz

wget https://gwas.mrcieu.ac.uk/files/met-d-PUFA/met-d-PUFA.vcf.gz

wget https://gwas.mrcieu.ac.uk/files/met-d-LA/met-d-LA.vcf.gz

wget https://gwas.mrcieu.ac.uk/files/met-d-DHA/met-d-DHA.vcf.gz

#LockeAE May 19 2021
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008673/LockeAE_prePMID_FAw3_sex-combined.gz 

wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008673/LockeAE_prePMID_FAw6_sex-combined.gz 

wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008673/LockeAE_prePMID_PUFA_sex-combined.gz 

wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008673/LockeAE_prePMID_LA_sex-combined.gz 

wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008673/LockeAE_prePMID_DHA_sex-combined.gz

gunzip *.gz

mv LockeAE_prePMID_FAw3_sex-combined LockeAE_Omega_3.txt

mv LockeAE_prePMID_FAw6_sex-combined LockeAE_Omega_6.txt

mv LockeAE_prePMID_PUFA_sex-combined LockeAE_PUFA.txt

mv LockeAE_prePMID_LA_sex-combined LockeAE_LA.txt

mv LockeAE_prePMID_DHA_sex-combined LockeAE_DHA.txt
