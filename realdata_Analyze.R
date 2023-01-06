library(TwoSampleMR)
library(Rcpp)
library(RcppArmadillo)
sourceCpp('/home/panwei/lin00374/mvcml/mvmr.cpp')
source('/home/panwei/lin00374/mvcml/MVMRcML.R')


args=(commandArgs(TRUE))
if(length(args)==0){
        print("No arguments supplied.")
##supply default values
}else{
        for(i in 1:length(args)){
                eval(parse(text=args[[i]]))
print(args[[i]])
                        }
}

exposure_dat = get(load('8rf_TMR.Rdata'))
rho_8rf = get(load('rho_mat_8rf.Rdata'))
id_outcome = 'ebi-a-GCST005195'
print(id_outcome)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, id_outcome)
mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)
res <- mv_multiple(mvdat)

Sig_l = Sig_inv_l = list()
se_mat = mvdat$exposure_se
for(i in 1:nrow(mvdat$exposure_beta)){
  rho_mat = diag(ncol(mvdat$exposure_beta)+1)
  rho_mat[1:ncol(mvdat$exposure_beta),1:ncol(mvdat$exposure_beta)] = rho_8rf
  Sig = rho_mat*crossprod(t(c(se_mat[i,],mvdat$outcome_se[i])))
  Sig_l[[i]] = Sig
  Sig_inv_l[[i]] = solve(Sig)
}
seed=20220428
set.seed(seed) 
mvmr_res = MVmr_cML_DP_c(mvdat$exposure_beta,as.matrix(mvdat$outcome_beta),se_mat,Sig_inv_l,
              n=133000,random_start=5,maxit=100,K_vec=0:40,num_pert=100)
mvmr_se = MVcML_SdTheta(b_exp=mvdat$exposure_beta,b_out=as.matrix(mvdat$outcome_beta),Sig_inv_l = Sig_inv_l,
                    theta = mvmr_res$BIC_theta,
                    zero_ind = setdiff(1:nrow(mvdat$exposure_beta),mvmr_res$BIC_invalid))

Res = list(mvdat = mvdat,
       rho_mat = rho_8rf,
       ivw_res = res,
       cML_res = mvmr_res,
       Sig_inv_l = Sig_inv_l,
       cML_se = mvmr_se,
       cML_pval = pnorm(-abs(mvmr_res$BIC_theta/mvmr_se))*2,
       cMLDP_pval = pnorm(-abs(mvmr_res$BIC_DP_theta/mvmr_res$BIC_DP_se))*2)
save(Res,file = paste0('8rf_',id_outcome,seed,'.Rdata'))


