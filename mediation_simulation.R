source('/home/panwei/lin00374/mvcml/GrantSIM21/simu_code.R')
source('/home/panwei/lin00374/mvcml/MVMRcML.R')
library(MendelianRandomization)
library(Rcpp)
library(RcppArmadillo)
sourceCpp('/home/panwei/lin00374/mvcml/mvmr.cpp')

generate_med_summary_gwas <- function(seed,n,m,m1,theta,k,K){
  set.seed(seed)
  gamma = sapply(1:k, function(j){runif(m, 0, 0.22)})
  bxu = 1/k
  byu = 1
  phi = alpha = rep(0,m)
  if(K>=1){
    alpha[1:K] = rnorm(K,0.1,0.2)
  }
  betax = gamma + bxu*phi
  
  betax[1:(m1+K),2] =  0.5 * betax[1:(m1+K),1]
  betax[(m1+K+1):m,2] = betax[(m1+K+1):m,2] + 0.5 * betax[(m1+K+1):m,1]

  betay = alpha + betax %*% theta + byu*phi
  
  sx = sqrt(1/n)
  sy = sqrt(1/n)
  
  rho_mat = matrix(c(1,0.5,0.5,1),nrow=2)
  betahat_x = betax + mvrnorm(n=m,mu=c(0,0),Sigma=rho_mat*sx^2)
  betahat_y = betay + rnorm(m, mean = 0, sd = sy)
  
  return(list(rho_mat = rho_mat, gamma=gamma,bx = betax, betax=betahat_x,sdx=matrix(sx,nrow=m,ncol=k),betay=as.vector(betahat_y),by=betay,sdy=rep(sy,m)))
}

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

array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# 
for(j in (50*(array_id-1)):(50*array_id-1)+1){
  simdat = generate_med_summary_gwas(n=n,m=m,m1=m1,theta=theta_vec,k=k,K=K,seed=j)
  b_exp = simdat$betax
  b_out = simdat$betay
  Sig_l = Sig_inv_l = list()
  for(i in 1:(m)){
    rho = 0.5
    rho_mat = matrix(rho,nrow=k,ncol=k)
    diag(rho_mat) = 1
    Sig_X = rho_mat*crossprod(t(simdat$sdx[i,]))
    Sig = cbind(Sig_X,rep(0,k))
    Sig = rbind(Sig,rep(0,k+1))
    Sig[k+1,k+1] = simdat$sdy[i]^2
    Sig_l[[i]] = Sig
    Sig_inv_l[[i]] = solve(Sig)
  }
  b_exp = simdat$betax
  b_out = simdat$betay
     mvcmldp_res = MVmr_cML_DP_c(b_exp = b_exp,b_out = as.matrix(b_out),se_bx = simdat$sdx,Sig_inv_l,n=n,maxit=100,random_start=10,num_pert=100)
     MVcMLDP_b = mvcmldp_res$BIC_DP_theta; MVcMLDP_se = mvcmldp_res$BIC_DP_se; MVcMLDP_pval = pnorm(-abs(MVcMLDP_b/MVcMLDP_se))*2
     MVcML_b = mvcmldp_res$BIC_theta;
     MVcML_se = MVcML_SdTheta(b_exp=b_exp,b_out=b_out,Sig_inv_l = Sig_inv_l,theta = MVcML_b,
                   zero_ind = setdiff(1:m,mvcmldp_res$BIC_invalid))
     MVcML_pval = pnorm(-abs(MVcML_b/MVcML_se))*2
  
  uvcmldp_res = MRcML::mr_cML_DP(b_exp = b_exp[,1],b_out = b_out,se_exp = simdat$sdx[,1],se_out = simdat$sdy,n=n,random_start=10)
  UVcMLDP_b = uvcmldp_res$BIC_DP_theta; UVcMLDP_se = uvcmldp_res$BIC_DP_se; UVcMLDP_pval = pnorm(-abs(UVcMLDP_b/UVcMLDP_se))*2
  UVcML_b = uvcmldp_res$BIC_theta; UVcML_se = uvcmldp_res$BIC_se; UVcML_pval = uvcmldp_res$BIC_p
  
     mv_object = MendelianRandomization::mr_mvinput(bx=b_exp,bxse=simdat$sdx,
                                                by=b_out,byse=simdat$sdy)
     mv_ivw = MendelianRandomization::mr_mvivw(mv_object)
     MVivw_b = mv_ivw@Estimate; MVivw_se = mv_ivw@StdError; MVivw_pval = mv_ivw@Pvalue
  uv_object = MendelianRandomization::mr_input(bx=b_exp[,1],bxse=simdat$sdx[,1],
                                               by=b_out,byse=simdat$sdy)
  uv_ivw = MendelianRandomization::mr_ivw(uv_object)
  UVivw_b = uv_ivw@Estimate; UVivw_se = uv_ivw@StdError; UVivw_pval = uv_ivw@Pvalue
  
  mrob = mr_mvinput(bx = b_exp, bxse = simdat$sdx, by = b_out, byse = simdat$sdy)
  th_egger = mr_mvegger(mrob)
  egger_b = th_egger@Estimate; egger_se = th_egger@StdError.Est; egger_pval = th_egger@Pvalue.Est
  th_robust = mvmr_robust(b_exp, b_out, simdat$sdy, k.max = 1000, maxit.scale = 1000)
  robust_b = th_robust$coefficients; robust_se = th_robust$se; robust_pval = pnorm(-abs(robust_b/robust_se))*2
  th_lass = tryCatch({
    mvmr_lasso(b_exp, b_out, simdat$sdy)
  },error = function(e){
    return(NULL)
  })
  if(is.null(th_lass)){lasso_b=lasso_se=lasso_pval=NA
  }else{
    lasso_b = th_lass$th_post; lasso_se = th_lass$se_post; lasso_pval = pnorm(-abs(lasso_b/lasso_se))*2}
  th_qr = mvmr_median(b_exp, simdat$sdx, b_out, simdat$sdy, boot=TRUE)
  median_b = th_qr$coefficients; median_se = th_qr$se; median_pval = pnorm(-abs(median_b/median_se))*2

  var_name =  c(ls(pattern='_b$'),ls(pattern='_se$'),ls(pattern='_pval$'))
  for(vn in var_name){
    t1 = paste0("write.table(t(",vn,"),'",vn,"_",array_id,".txt',quote=F,row.names=F,append=T,col.names=F)")
    eval(parse(text=t1))
  }
}
