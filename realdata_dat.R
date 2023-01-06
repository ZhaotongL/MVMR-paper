library(TwoSampleMR)
# TG: ebi-a-GCST002216
# LDL: ebi-a-GCST002222
# HDL: ebi-a-GCST002223
# Height: ieu-a-89 
# BMI: ieu-a-835
# DBP: ieu-b-39
# SBP: ieu-b-38
# FG: ebi-a-GCST005186


gwas_id = c(
            'ebi-a-GCST002216',
            'ebi-a-GCST002222',
            'ebi-a-GCST002223',
            'ieu-a-89',
            'ieu-a-835',
            'ieu-b-38',
            'ieu-b-39',
            'ebi-a-GCST005186')

load('/home/panwei/lin00374/cML/graph/vcffiles/ldsc_res/rho_thres_0.1.RData')
rho_8rf = rho_thres_0.1[c("TG","LDL","HDL","BMI","Height","FG","SBP","DBP"),
                  c("TG","LDL","HDL","BMI","Height","FG","SBP","DBP")]
save(rho_8rf,file='rho_mat_8rf.Rdata')

exposure_dat = mv_extract_exposures(gwas_id)
save(exposure_dat,file='8rf_TMR.Rdata')
