rm(list=ls())
library(TwoSampleMR)
exposure_dat <-read_exposure_data(
  filename = 'Exposure data for snoring.csv',
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

outcome_dat <- extract_outcome_data(
  snps = exposure_dat$SNP,
  outcomes ='ebi-a-GCST006414')

dat <- harmonise_data(exposure_dat,outcome_dat, action = 1)
dat <- dat[dat$pval.outcome>5e-8,] 
mr_results <- mr(dat,method_list=c( "mr_ivw_fe",
                                    "mr_ivw_mre",
                                    "mr_weighted_median",
                                    "mr_egger_regression",
                                    "mr_raps"))
het <- mr_heterogeneity(dat)
ple <- mr_pleiotropy_test(dat)
OR_mr_results <- generate_odds_ratios(mr_results)

res_loo <- mr_leaveoneout(dat,method = mr_ivw_fe)
p <- mr_leaveoneout_plot(res_loo)
p

exposure_dat2 <- exposure_dat[!exposure_dat$SNP%in%c("rs2307111","rs34811474","rs8069947",
                                                     "rs2049045","rs6054427"),]
###excluding the pleiotropic SNPs
dat2 <- harmonise_data(exposure_dat2,outcome_dat, action = 1)
dat2 <- dat2[dat2$pval.outcome>5e-8,] 
mr_results2 <- mr(dat2,method_list= "mr_ivw_fe")
OR_mr_results2 <- generate_odds_ratios(mr_results2)

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
mr_presso