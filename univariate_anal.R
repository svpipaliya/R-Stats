### This script is used to perform univariate analysis to figure out which variables are significantly associated with the susceptibility phenotype

## Step 1: Run a univariate model for all your variables against your phenotype

## load variables and phenotypes in a dataframe
variables <- "/Users/ShwetaPipaliya/Documents/PostDoc/projects/gbs_gwas/raw_data/All.Pheno.GBS.csv"
df.GBS.variable <- read.csv(variables, header = T)

## remove rows that are not relevant variables for analysis
df.GBS.variable[ ,c('Birth_date', 'Blood_samples', 'Consent', 'Exclusion_crit', 'S_agalactiar')] <- list(NULL)

## Run model of variables that are kept
log.model.Age <- glm(Colonized ~ Age, data = df.GBS.variable)
summary(log.model.Age )
pv.age <- summary(log.model.Age)$coef[, "Pr(>|t|)"]
print(pv.age)

log.model.ethnicity <- glm(Colonized~ Ethnic_Group, data = df.GBS.variable)
summary(log.model.ethnicity )
pv.eth <- summary(log.model.ethnicity)$coef[, "Pr(>|t|)"]
print(pv.eth)

log.model.primiparous <- glm(Colonized ~ Primiparous, data = df.GBS.variable)
summary(log.model.primiparous )
pv.pri <- summary(log.model.primiparous)$coef[, "Pr(>|t|)"]
print(pv.pri)

log.model.ga <- glm(Colonized ~ GA_new, data = df.GBS.variable)
summary(log.model.ga)
pv.ga <- summary(log.model.ga)$coef[, "Pr(>|t|)"]
print(pv.ga)

log.model.all_gn <- glm(Colonized ~ All_Gr_Neg, data = df.GBS.variable)
summary(log.model.all_gn)
pv.uti <- summary(log.model.all_gn)$coef[, "Pr(>|t|)"]
print(pv.uti)

log.model.ecoli <- glm(Colonized ~ E.coli, data = df.GBS.variable)
summary(log.model.ecoli)
pv.ecoli <- summary(log.model.ecoli)$coef[, "Pr(>|t|)"]
print(pv.ecoli)

log.model.kleb <- glm(Colonized ~ Klebsiella, data = df.GBS.variable)
summary(log.model.kleb)
pv.kleb <- summary(log.model.kleb)$coef[, "Pr(>|t|)"]
print(pv.kleb)

log.model.infect <- glm(Colonized ~ Infections, data = df.GBS.variable)
summary(log.model.infect)
pv.infect <- summary(log.model.infect)$coef[, "Pr(>|t|)"]
print(pv.infect)
print(log.model.infect)

log.model.smoke <- glm(Colonized ~ Smoke, data = df.GBS.variable)
summary(log.model.smoke)
pv.smoke <- summary(log.model.smoke)$coef[, "Pr(>|t|)"]
print(pv.smoke)

log.model.ht <- glm(Colonized ~ All_HT, data = df.GBS.variable)
summary(log.model.ht)
pv.ht <- summary(log.model.ht)$coef[, "Pr(>|t|)"]
print(pv.ht)

log.model.dm <- glm(Colonized ~ All_DM, data = df.GBS.variable)
summary(log.model.dm)
pv.dm <- summary(log.model.dm)$coef[, "Pr(>|t|)"]
print(pv.dm)

log.model.thy <- glm(Colonized ~ All_thyr_dis, data = df.GBS.variable)
summary(log.model.thy)
pv.thy <- summary(log.model.thy)$coef[, "Pr(>|t|)"]
print(pv.thy)

log.model.obs <- glm(Colonized ~ Obesity, data = df.GBS.variable)
summary(log.model.obs)
pv.obs <- summary(log.model.obs)$coef[, "Pr(>|t|)"]
print(pv.obs)

log.model.uti <- glm(Colonized ~ UTI)

## autoimmunity variables
log.model.AutoImm_thyr <- glm(Colonized ~ AutoImm_thyr_dis, data = df.GBS.variable)
summary(log.model.AutoImm_thyr)
pv.auto_thy <- summary(log.model.AutoImm_thyr)$coef[, "Pr(>|t|)"]
print(pv.auto_thy)

log.model.TPO <- glm(Colonized ~ TPO, data = df.GBS.variable)
summary(log.model.TPO)
pv.TPO <- summary(log.model.TPO)$coef[, "Pr(>|t|)"]
print(pv.TPO)

log.model.TG <- glm(Colonized ~ TG, data = df.GBS.variable)
summary(log.model.TG)
pv.TG <- summary(log.model.TG)$coef[, "Pr(>|t|)"]
print(pv.TG)

log.model.TRAB <- glm(Colonized ~ TRAB, data = df.GBS.variable)
summary(log.model.TRAB)
pv.TRAB <- summary(log.model.TRAB)$coef[, "Pr(>|t|)"]
print(pv.TRAB)

log.model.ANA <- glm(Colonized ~ ANA, data = df.GBS.variable)
summary(log.model.ANA)
pv.ANA <- summary(log.model.ANA)$coef[, "Pr(>|t|)"]
print(pv.ANA)

log.model.Clin_autoimm_noThyr <- glm(Colonized ~ Clin_autoimm_noThyr, data = df.GBS.variable)
summary(log.model.Clin_autoimm_noThyrA)
pv.Clin_autoimm_noThyr <- summary(log.model.Clin_autoimm_noThyr)$coef[, "Pr(>|t|)"]
print(pv.Clin_autoimm_noThyr)

log.model.Clin_inflam_dis <- glm(Colonized ~ Clin_inflam_dis, data = df.GBS.variable)
summary(log.model.Clin_inflam_dis)
log.model.Clin_inflam_dis <- summary(log.model.Clin_inflam_dis)$coef[, "Pr(>|t|)"]
print(log.model.Clin_inflam_dis)

# Step 2: Select the ones with significant p-values (ex. var 2 and 3) and add them all up using “+” in a multivariate model.
log.model.final <- glm(Colonized ~ Klebsiella + 
                         Primiparous + 
                         AutoImm_thyr_dis + 
                         ANA + 
                         Clin_autoimm_noThyr +
                         Clin_inflam_dis + 
                         Obesity,
                       data = df.GBS.variable)
summary(log.model.final)
pv.final <- summary(log.model.final)$coef[, "Pr(>|t|)"]
print(pv.final)

