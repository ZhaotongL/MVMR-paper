################################################################################
#Produces simulated datasets
################################################################################
source('/home/panwei/lin00374/mvcml/GrantSIM21/simu_code.R')
library(parallel)
library(MASS)
cl = makeCluster(50)
clusterEvalQ(cl, library(MASS))

n = 50000
m = 0.3 # MAF
k = 4 # L=4
p = list('20' = 20) # m=20
sa = 0.2
sg = 0.2
g = rep(1/k, k)
bx_min = 0
bx_max = 0.22
th = list('thA' = c(0.2, 0.1, 0.3, 0.4), 'th0' = c(0, -0.1, 0.1, 0.2))
s = list('30' = 0.3, '50' = 0.5)
mua = list('S1' = 0, 'S2' = 0.1, 'S3' = 0, 'S4' = 0)
mug = list('S1' = 0, 'S2' = 0, 'S3' = 0.1, 'S4' = 0.05)
SigX = list('unc' = diag(k))
M = 500

clusterExport(cl, c('sstat', 'Data_sim_ind', 'Data_sim_irre','Data_sim_irre1','Data_sim_weak',
                    'Data_sim_ind_inside', 'Data_sim_ind_med', 'n', 'm', 'k',
                    'p', 'sa', 'sg', 'g', 'bx_min', 'bx_max', 'th', 's', 'mua',
                    'mug', 'SigX', 'M'))

################################################################################
#Primary simulations
################################################################################

################################################################################

################################################################################
#p = 20 sims
################################################################################

################################################################################

clusterSetRNGStream(cl, 20220518)
## some have alpha some irrelevant
D_S4_30_p20 = parLapply(cl, 1:M, function(j){
 Data_sim_irre1(n, m, p$'20', s = s$'30', k, mua = 0.1, mug = 0, sa = 0.1, sg = 0,
              g, theta = th$thA, bx_min, bx_max, SigX = SigX$unc)
})
D_S4_30_p20_null = parLapply(cl, 1:M, function(j){
 Data_sim_irre1(n, m, p$'20', s = s$'30', k, mua = 0.1, mug = 0, sa = 0.1, sg = 0,
              g, theta = th$th0, bx_min, bx_max, SigX = SigX$unc)
})
save(D_S4_30_p20,D_S4_30_p20_null,
   file = 'simDat_n50000/sim_p20_irre_inv_dat.Rdata')




theta = 0.2
clusterSetRNGStream(cl, 20220404)
##S1
D_S1_30_p20 = parLapply(cl, 1:M, function(j){
 Data_sim_ind(n, m, p$'20', s = s$'30', k, mua = mua$S1, mug = mug$S1, sa, sg,
              g, theta = th$thA, bx_min, bx_max, SigX = SigX$unc)
})
D_S1_50_p20 = parLapply(cl, 1:M, function(j){
 Data_sim_ind(n, m, p$'20', s = s$'50', k, mua = mua$S1, mug = mug$S1, sa, sg,
              g, theta = th$thA, bx_min, bx_max, SigX = SigX$unc)
})

##S2
D_S2_30_p20 = parLapply(cl, 1:M, function(j){
 Data_sim_ind(n, m, p$'20', s = s$'30', k, mua = mua$S2, mug = mug$S2, sa, sg,
              g, theta = th$thA, bx_min, bx_max, SigX = SigX$unc)
})
D_S2_50_p20 = parLapply(cl, 1:M, function(j){
 Data_sim_ind(n, m, p$'20', s = s$'50', k, mua = mua$S2, mug = mug$S2, sa, sg,
              g, theta = th$thA, bx_min, bx_max, SigX = SigX$unc)
})

##S3
D_S3_30_p20 = parLapply(cl, 1:M, function(j){
 Data_sim_ind(n, m, p$'20', s = s$'30', k, mua = mua$S3, mug = mug$S3, sa, sg,
              g, theta = th$thA, bx_min, bx_max, SigX = SigX$unc)
})
D_S3_50_p20 = parLapply(cl, 1:M, function(j){
 Data_sim_ind(n, m, p$'20', s = s$'50', k, mua = mua$S3, mug = mug$S3, sa, sg,
              g, theta = th$thA, bx_min, bx_max, SigX = SigX$unc)
})

save(D_S1_50_p20,D_S2_50_p20,D_S3_50_p20,D_S1_30_p20,D_S2_30_p20,D_S3_30_p20,
   file='simDat_n50000/sim_p20_dat.Rdata')

## theta = 0 ##
#S1
D_S1_30_p20_null = parLapply(cl, 1:M, function(j){
 Data_sim_ind(n, m, p$'20', s = s$'30', k, mua = mua$S1, mug = mug$S1, sa, sg,
              g, theta = th$th0, bx_min, bx_max, SigX = SigX$unc)
})
D_S1_50_p20_null = parLapply(cl, 1:M, function(j){
 Data_sim_ind(n, m, p$'20', s = s$'50', k, mua = mua$S1, mug = mug$S1, sa, sg,
              g, theta = th$th0, bx_min, bx_max, SigX = SigX$unc)
})

##S2
D_S2_30_p20_null = parLapply(cl, 1:M, function(j){
 Data_sim_ind(n, m, p$'20', s = s$'30', k, mua = mua$S2, mug = mug$S2, sa, sg,
              g, theta = th$th0, bx_min, bx_max, SigX = SigX$unc)
})
D_S2_50_p20_null = parLapply(cl, 1:M, function(j){
 Data_sim_ind(n, m, p$'20', s = s$'50', k, mua = mua$S2, mug = mug$S2, sa, sg,
              g, theta = th$th0, bx_min, bx_max, SigX = SigX$unc)
})

##S3
D_S3_30_p20_null = parLapply(cl, 1:M, function(j){
 Data_sim_ind(n, m, p$'20', s = s$'30', k, mua = mua$S3, mug = mug$S3, sa, sg,
              g, theta = th$th0, bx_min, bx_max, SigX = SigX$unc)
})
D_S3_50_p20_null = parLapply(cl, 1:M, function(j){
 Data_sim_ind(n, m, p$'20', s = s$'50', k, mua = mua$S3, mug = mug$S3, sa, sg,
              g, theta = th$th0, bx_min, bx_max, SigX = SigX$unc)
})

stopCluster(cl)

save(D_S1_50_p20_null,D_S2_50_p20_null,D_S3_50_p20_null,D_S1_30_p20_null,D_S2_30_p20_null,D_S3_30_p20_null,
   file='simDat_n50000/sim_p20_null_dat.Rdata')

