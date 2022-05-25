### This script is used to perform used to generate and plot correlation scores between different clinical variables and GBS colonization phenotypes ##
### The R package 'lares' is used against the complete All.pheno.GBS dataframe consisting of both categorical and continous variables ##

#install.packages("lares")
library(lares)

## load variables and phenotypes in a dataframe
variables <- "/Users/ShwetaPipaliya/Documents/PostDoc/projects/gbs_gwas/raw_data/All.Pheno.GBS.csv"
df.GBS.variable <- read.csv(variables, header = T)

## remove rows that are not relevant variables for analysis
df.GBS.variable[ ,c('Susceptible', 'Birth_date', 'Blood_samples', 'Consent', 'Exclusion_crit', 'UTI_curr_pregn', 'All_Gr_Neg', 'S_agalactiar', 'Chr_HT', 'Gravidic_HT', 'AutoImm_thyr_dis', 'TPO', 'TG', 'TRAB', 'ANA', 'Clin_autoimm', 'Clin_autoimm_noThyr', 'Clin_inflam_dis', 'AUTO_Abs', 'UTI_curr_preg_confirmed')] <- list(NULL)

## examine colinearity between variables using the corr_cross function within the dataframe
corr_cross(df.GBS.variable , type = 1, top = 15, max_pvalue = 0.05)

## calculating correlation of different variables without any statistical method
df.GBS.variable %>% corr_var(Colonized, type = 1, max_pvalue = 0.05, ranks = TRUE)


## calculating correlation of different variables using the Kendall 
#df.GBS.variable %>% corr_var(Colonized, top = 15, method = "kendall", max_pvalue = 0.05, ranks = TRUE)

## calculating correlation of different variables without any statistical method
#df.GBS.variable %>% corr_var(Colonized, top = 15, method = "pearson", max_pvalue = 0.05)

