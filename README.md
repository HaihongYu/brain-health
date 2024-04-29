# brain-health
###Association between the frequency of adding salt to food and brain function

#Normalize and standardize

library(nortest) 
library(car)
brainfunction$PHQ4Z<-bcPower((brainfunction$PHQ4), powerTransform(brainfunction$PHQ4)$roundlam)
brainfunction$PHQ4Z<-scale(brainfunction$PHQ4Z)

#cross-sectional association

F=lm(PHQ4Z~ salt_frequency+age+sex+ethnic+education+TDI, data=brainfunction)
summary(F)
F=lm(PHQ4Z~ as.factor(salt_frequency)+age+sex+ethnic+education+TDI, data=brainfunction)
summary(F)

#longitudinal association 

library(lavaan)
model<-'

A.bl=~`salt_frequency`

A.fo=~`salt_frequency_fo`

B.bl=~`PHQ40.0`

B.fo=~`PHQ42.0`

A.fo~b1*A.bl+b2*B.bl+age+sex+ethnic+education+TDI

B.fo~b3*A.bl+b4*B.bl+age+sex+ethnic+education+TDI

A.bl~~ B.bl

A.fo~~ B.fo

'
fit<-sem(model,data=crosslagg,se="boot",bootstrap=10000,estimator="ML")

summary(fit,fit.measures=T,standardized=TRUE)

fitMeasures(fit,c("chisq","df","pvalue","gfi","cfi","rmr","srmr","rmsea"))

###Association between the frequency of adding salt to food and brain disorders

library(survival)
library(survminer)
library(ggplot2)
library(rms)

#Model 1

fit <- coxph(Surv(cox_model1$Bipolar_days,cox$Bipolar_status)~as.factor(salt_frequency)+age+sex+ethnic+education+TDI,data=cox_model1)
test.ph<-cox.zph(fit)
test.ph
summary(fit)

#Model 2

fit <- coxph(Surv(cox_model2$Dementia_days,cox_model2$Dementia_status)~as.factor(salt_frequency)+age+sex+ethnic+education+TDI
                  +BMI+smoke+alchol+cholesteral+CKD+DM+CVD+HTN,data=cox_model2)
test.ph<-cox.zph(fit)
test.ph
summary(fit)

#Model 3 

fit <-coxph(Surv(cox_diet$Sleep_days,cox_diet$Sleep_status)~as.factor(salt_frequency)+age+sex+ethnic+education+TDI
                  +BMI+smoke+alchol+cholesteral+CKD+DM+CVD+HTN+processed_meat+red_meat+fish+fruit+vegetable+energy_kcal,data=cox_diet)
test.ph<-cox.zph(fit)
test.ph
summary(fit)

##Sensitivity analysis excluding excluding excluding individuals who had major changes in their diet in last 5 years

fit <- coxph(Surv(cox_nochange$Sleep_days,cox_nochange$Sleep_status)~as.factor(salt_frequency)+age+sex+ethnic+education+TDI,data=cox_nochange)
test.ph<-cox.zph(fit)
test.ph
summary(fit)

##Subgroup analysis
#Age (<60/≥60）

fit <- coxph(Surv(cox_age$Schizophrenia_days,cox_age$Schizophrenia_status)~as.factor(salt_frequency)+sex+ethnic+education+TDI,data=cox_age)
test.ph<-cox.zph(fit)
test.ph
summary(fit)

#Sex (women/men)

fit <- coxph(Surv(cox_men$Schizophrenia_days,cox_men$Schizophrenia_status)~as.factor(salt_frequency)+age+ethnic+education+TDI,data=cox_men)
test.ph<-cox.zph(fit)
test.ph
summary(fit)

#Fruits&vegetables (< median/≥ median [5.5])

fit <- coxph(Surv(cox_diet_high$Dementia_days,cox_diet_high$Dementia_status)~as.factor(salt_frequency)+sex+age+ethnic+education+TDI
                  +BMI+smoke+alchol+cholesteral+CKD+DM+CVD+HTN+processed_meat+red_meat+fish+energy_kcal,data=cox_diet_high)
test.ph<-cox.zph(fit)
test.ph
summary(fit)

###MR analysis

library(TwoSampleMR)
exp_salt <- extract_instruments(outcomes = 'ukb-b-8121')

out_dementia <- extract_outcome_data(snps = exp_salt$SNP, outcomes = 'finn-b-F5_DEMENTIA')

out_eplepsy<-read_outcome_data(snps = exp_salt$SNP, filename = "finngen_R9_G6_EPLEPSY.gz",
                             sep = ",", snp_col = "SNP",
                             beta_col = "beta", se_col = "sebeta",
                             effect_allele_col = "ref",
                             other_allele_col = "alt",
                             pval_col = "pval")
                             
MR_dementia <- harmonise_data(exp_salt, out_dementia)

MR_dementia_result <- mr(MR_dementia)

MR_dementia_result

MR_dementia_result <-generate_odds_ratios(MR_dementia_result)

MR_dementia_result

het <- mr_heterogeneity(MR_dementia)
het

pleio <- mr_pleiotropy_test(MR_dementia)
pleio



###SEM

sem.model1<-'

         salt=~1* salt_frequency
         
         PRS=~1* Score
         
         inflammation=~1* X30710_0.0+X30730_0.0+X30000_0.0+X30140_0.0+X30720_0.0

         metabolism=~1* X23443+X23421++X23480+X23535+X23549
         
         status=~1* Dementia_status++Stroke_status+Epilepsy_status+Anxiety_status+MDD_status+Sleep_status+Bipolar_status+Schizophrenia_status
         
         status~ salt+PRS+inflammation+metabolism
         
         inflammation~ salt
         
         metabolism~ salt
         
         salt~ PRS
         
         inflammation~ PRS

         metabolism~ PRS
         
        metabolism~~ inflammation
        
        '
fit<-sem(sem.model1,data=sem5_PRS,se="boot",bootstrap=10000,estimator="ML")

summary(fit,fit.measures=T,standardized=TRUE)

fitMeasures(fit,c("chisq","df","pvalue","gfi","cfi","rmr","srmr","rmsea"))



