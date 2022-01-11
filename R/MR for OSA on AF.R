rm(list=ls())
library(TwoSampleMR)
###Reading in exposure data from a file with the column names of "SNP,beta,se,pval,effect_allele,other_allele,Phenotype,eaf,gene"
exposure_dat <-read_exposure_data(
    filename = 'Exposure data for OSA.csv',
  sep = ',',
  snp_col = 'SNP',
  beta_col = 'beta',
  se_col = 'se',
  effect_allele_col = 'effect_allele',
  phenotype_col = 'Phenotype',
  units_col = 'units',
  other_allele_col = 'other_allele',
  eaf_col = 'eaf',
  samplesize_col = 'samplesize',
  ncase_col = 'ncase',
  ncontrol_col = 'ncontrol',
  gene_col = 'gene',
  pval_col = 'pval'
)

###Extracting OSA-associated SNPs from GWASs on atrial fibrillation
outcome_dat <- extract_outcome_data(
    snps = exposure_dat$SNP,
    outcomes = c('ebi-a-GCST006414','ebi-a-GCST006061')
    )
###ebi-a-GCST006414 is the ID-number for GWAS meta-analyses of Nielsen et al.
###ebi-a-GCST006061 is the ID-number for AF GWAS from AFGEN.
###Harmonise data
dat <- harmonise_data(exposure_dat,outcome_dat, action = 1)
dat <- dat[dat$pval.outcome>5e-8,] ###To ensure no SNP directly associated with outcomes.

###Perform MR analyses with inverse variance weighted (fixed effects), 
###weighted median, MR-Egger, and MR-RAPS methods.
mr_results <- mr(dat,method_list=c("mr_ivw_fe",
                                   "mr_weighted_median",
                                   "mr_egger_regression",
                                   "mr_raps"
                                   ))
###Sensitivity analyses
het <- mr_heterogeneity(dat)###Heterogeneity statistics
ple <- mr_pleiotropy_test(dat)###Horizontal pleiotropy

###Generate odds ratios with 95% confidence intervals
OR_mr_results <- generate_odds_ratios(mr_results)

###Leave-one-out analysis
res_loo <- mr_leaveoneout(dat,method = mr_ivw_fe)
p <- mr_leaveoneout_plot(res_loo)
p


exposure_dat2 <- exposure_dat[!exposure_dat$SNP%in%"rs9937053",]
###excluding the pleiotropic SNP
dat2 <- harmonise_data(exposure_dat2,outcome_dat, action = 1)
mr_results2 <- mr(dat2,method_list= "mr_ivw_fe")
OR_mr_results2 <- generate_odds_ratios(mr_results2)


##MR-PRESSO
library(MRPRESSO)
mr_presso <- mr_presso(BetaOutcome = "beta.outcome", 
                    BetaExposure = "beta.exposure", 
                    SdOutcome = "se.outcome", 
                    SdExposure = "se.exposure", 
                    OUTLIERtest = TRUE, 
                    DISTORTIONtest = TRUE, 
                    data = dat, 
                    NbDistribution = 2500,  
                    SignifThreshold = 0.05)


rm(list=ls())
###Perform MVMR
library('MendelianRandomization')
mv_data <- read.csv(file='mvmrdat.csv',row.names = 1) 
###read in a file with the information of beta and se for OSA, AF and risk factors, respectively.

###for example, perform MVMR of OSA on AF adjusting for BMI.
MRInputObject <- mr_mvinput(bx = cbind(mv_data$beta_OSA,mv_data$beta_BMI),
                          bxse = cbind(mv_data$se_OSA,mv_data$se_BMI),
                          by = mv_data$beta_AF,
                          byse = mv_data$se_AF)
tem <- mr_mvivw(MRInputObject, model = "fixed", correl = FALSE,
               distribution = "normal", alpha = 0.05)
