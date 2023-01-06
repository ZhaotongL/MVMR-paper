source('/home/panwei/lin00374/mvcml/GrantSIM21/simu_code.R')
source('/home/panwei/lin00374/mvcml/MVMRcML.R')
library(Rcpp)
library(RcppArmadillo)
sourceCpp('/home/panwei/lin00374/mvcml/mvmr.cpp')
library(MendelianRandomization)

args=(commandArgs(TRUE))
if(length(args)==0){
        print("No arguments supplied.")
##supply default values
}else{
        for(i in 1:length(args)){
                eval(parse(text=args[[i]]))
                        }
}

array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

#load('/home/panwei/lin00374/mvcml/GrantSIM21/simDat_n50000/sim_p20_dat.Rdata')
#load('/home/panwei/lin00374/mvcml/GrantSIM21/simDat_n50000/sim_p20_null_dat.Rdata')
#load('/home/panwei/lin00374/mvcml/GrantSIM21/simDat_n50000/sim_p20_irre_inv_dat.Rdata')

k = 4
g = rep(1/k, k)
eval(parse(text=paste0('D=',scenario)))
print(scenario)
for(j in (25*(array_id-1)):(25*array_id-1)+1){
    simdat = D[[j]]
    p = length(simdat$seby)
    Sig_l = Sig_inv_l = list()
    for(i in 1:p){
      rho_mat = simdat$rho_mat
      Sig_X = rho_mat*crossprod(t(simdat$sebx[i,]))
      Sig = cbind(Sig_X,rep(0,k))
      Sig = rbind(Sig,rep(0,k+1))
      Sig[k+1,k+1] = simdat$seby[i]^2
      Sig_l[[i]] = Sig
      Sig_inv_l[[i]] = solve(Sig)
    }
    mvcmldp_res = MVmr_cML_DP_c(b_exp = simdat$bxhat,b_out = as.matrix(simdat$byhat),se_bx=simdat$sebx,Sig_inv_l,n=50000,maxit=100,random_start=10,num_pert=100)
    cMLDP_b = mvcmldp_res$BIC_DP_theta; cMLDP_se = mvcmldp_res$BIC_DP_se; cMLDP_pval = pnorm(-abs(cMLDP_b/cMLDP_se))*2
    cML_b = mvcmldp_res$BIC_theta;
    cML_se = MVcML_SdTheta(b_exp=simdat$bxhat,b_out=simdat$byhat,Sig_inv_l = Sig_inv_l,theta = cML_b,
                  zero_ind = setdiff(1:p,mvcmldp_res$BIC_invalid))
    cML_pval = pnorm(-abs(cML_b/cML_se))*2
    mrob = mr_mvinput(bx = simdat$bxhat, bxse = simdat$sebx, by = simdat$byhat, byse = simdat$seby)
    th_mrivw = mr_mvivw(mrob)
    ivw_b = th_mrivw@Estimate; ivw_se = th_mrivw@StdError; ivw_pval = th_mrivw@Pvalue
    th_egger = mr_mvegger(mrob)
    egger_b = th_egger@Estimate; egger_se = th_egger@StdError.Est; egger_pval = th_egger@Pvalue.Est
    th_robust = mvmr_robust(simdat$bxhat, simdat$byhat, simdat$seby, k.max = 1000, maxit.scale = 1000)
    robust_b = th_robust$coefficients; robust_se = th_robust$se; robust_pval = pnorm(-abs(robust_b/robust_se))*2
    th_lass = mvmr_lasso(simdat$bxhat, simdat$byhat, simdat$seby)
    lasso_b = th_lass$th_post; lasso_se = th_lass$se_post; lasso_pval = pnorm(-abs(lasso_b/lasso_se))*2
    th_qr = mvmr_median(simdat$bxhat, simdat$sebx, simdat$byhat, simdat$seby, boot=TRUE)
    median_b = th_qr$coefficients; median_se = th_qr$se; median_pval = pnorm(-abs(median_b/median_se))*2

    var_name =  c(ls(pattern='_b$'),ls(pattern='_se$'),ls(pattern='_pval$'))
    for(vn in var_name){
        t1 = paste0("write.table(t(",vn,"),'",vn,"_",array_id,".txt',quote=F,row.names=F,append=T,col.names=F)")
        eval(parse(text=t1))
    }
}
    


