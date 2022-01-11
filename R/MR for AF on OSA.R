rm(list=ls())
library(TwoSampleMR)
{exposure_dat <-read_exposure_data(
  filename = 'Exposure data for AF.csv',
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
)}

###file summary_stats_finngen_R5_G6_SLEEPAPNO.txt is available for download under the phenocode G6_SLEEPAPNO at https://www.finngen.fi/en
outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "summary_stats_finngen_R5_G6_SLEEPAPNO.txt", 
  "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = 'maf',
  pval_col = "pval",
)
outcome_dat$outcome <- "obstructive sleep apnea" 

dat <- harmonise_data(exposure_dat,outcome_dat, action = 1)
dat <- dat[dat$pval.outcome>5e-8,] 
mr_results <- mr(dat,method_list=c( "mr_ivw_fe",
                                    "mr_ivw_mre",
                                    "mr_weighted_median",
                                    "mr_egger_regression",
                                    "mr_raps"))
het <- mr_heterogeneity(dat)
ple <- mr_pleiotropy_test(dat)
OR_mr_results=generate_odds_ratios(mr_results)

res_loo <- mr_leaveoneout(dat,method = mr_ivw_fe)
p <- mr_leaveoneout_plot(res_loo)
p

exposure_dat2 <- exposure_dat[!exposure_dat$SNP%in%c("rs10006327","rs12245149","rs12604076",
                                                    "rs1458038","rs1563304","rs2540949",
                                                    "rs284277","rs2885697","rs35005436",
                                                    "rs422068","rs4951258","rs56201652",
                                                    "rs60212594","rs9899183"),]
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
