source('/home/panwei/lin00374/mvcml/GrantSIM21/simu_code.R')
source('/home/panwei/lin00374/mvcml/MVMRcML.R')
library(Rcpp)
library(RcppArmadillo)
sourceCpp('/home/panwei/lin00374/mvcml/mvmr.cpp')
library(MendelianRandomization)
library(MVMR)

args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
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

pimax = 0.1
L = 45
n = 20000
delta = 0
tau_true = 0 # tau_true = 0.5  for invalid IV case


for(j in (20*(array_id-1)):(20*array_id-1)+1){
    set.seed(j)
    sim = 1
    kx = 2
    mu = c(0,0,0)
    s = matrix(c(1,0.8,-0.7,0.8,1,-0.7,-0.7,-0.7,1), ncol = 3)
    pi2 = as.matrix(c(runif((L),0,pimax)),nrow=L,ncol=n)
    if(delta==0){
        pi1 = as.matrix(c(runif((L),0,pimax)),nrow=L,ncol=n)
    }else{
        pi1 = delta*pi2 + c(0,rnorm((L-1),0,0.075))
    }
    alpha = as.matrix(c(rnorm(L/3,0,tau_true), rep(0,L-L/3)),nrow=L,ncol=1)

  
  G = matrix(nrow = n, ncol=L)
  G[,]= rbinom((n*L),2,0.5)
  V = mvrnorm(n, mu, s)
 
  x1 = G%*%pi1 + V[,2]
  x2 = G%*%pi2 + V[,3]
  x3 = G%*%alpha  #new
  y = theta1*x1 + theta2*x2 + x3 + V[,1]
  
  x = cbind(x1, x2)
    ##Second model for two sample estimation
  G2 = matrix(nrow = n, ncol=L)
  G2[,]= rbinom((n*L),2,0.5)
   V2 = mvrnorm(n, mu, s) 

  x12 = G2%*%pi1  + V2[,2]
  x22 = G2%*%pi2  + V2[,3]
  x32 = G2%*%alpha#new
  y2 = theta1*x12 + theta2*x22 +  x32 + V2[,1]
  
 pihat = stderr = matrix(nrow = L, ncol = kx)
  gammahat = segamma = NULL
 
  for(k in 1:L){
    
    gammahat[k] = summary(lm(y2~ G2[,k]))$coefficient[2,1]
    segamma[k] = summary(lm(y2~ G2[,k]))$coefficient[2,2]
    
   for(ex in 1:kx){
    
   pihat[k,ex] = summary(lm(x[,ex]~ G[,k]))$coefficient[2,1]
   stderr[k,ex] = summary(lm(x[,ex]~ G[,k]))$coefficient[2,2]
    
  }
  } 
  
  correlation = matrix(c(1, cor(x1,x2), cor(x1,x2), 1), nrow = 2, ncol = 2)
    simdat = list(bxhat=pihat,byhat=gammahat,sebx=stderr,seby=segamma,rho_mat=correlation)
    F.data <- format_mvmr(BXGs = pihat, BYG=gammahat, seBXGs=stderr, seBYG=segamma)
    Xcovmat<-phenocov_mvmr(correlation,stderr)
    sres2 <- strength_mvmr(r_input =F.data, gencov = Xcovmat)
    write.table(t(as.numeric(sres2)),paste0('Fstat_',array_id,".txt"),quote=F,row.names=F,col.names=F,append=T)


    Sig_l = Sig_inv_l = list()
    for(i in 1:L){
      rho_mat = simdat$rho_mat
      Sig_X = rho_mat*crossprod(t(simdat$sebx[i,]))
      Sig = cbind(Sig_X,rep(0,kx))
      Sig = rbind(Sig,rep(0,kx+1))
      Sig[kx+1,kx+1] = simdat$seby[i]^2
      Sig_l[[i]] = Sig
      Sig_inv_l[[i]] = solve(Sig)
    }

    mvcmldp_res = MVmr_cML_DP_c(b_exp = simdat$bxhat,b_out = as.matrix(simdat$byhat),se_bx=simdat$sebx,Sig_inv_l,n=n,maxit=100,random_start=10,num_pert=100,fix_bx=fix)
    cMLDP_b = mvcmldp_res$BIC_DP_theta; cMLDP_se = mvcmldp_res$BIC_DP_se; cMLDP_pval = pnorm(-abs(cMLDP_b/cMLDP_se))*2
    cML_b = mvcmldp_res$BIC_theta;
    cML_se = MVcML_SdTheta(b_exp=simdat$bxhat,b_out=simdat$byhat,Sig_inv_l = Sig_inv_l,theta = cML_b,
                  zero_ind = setdiff(1:length(simdat$byhat),mvcmldp_res$BIC_invalid))
    cML_pval = pnorm(-abs(cML_b/cML_se))*2
    mrob = mr_mvinput(bx = simdat$bxhat, bxse = simdat$sebx, by = simdat$byhat, byse = simdat$seby)
    th_mrivw = mr_mvivw(mrob)
    ivw_b = th_mrivw@Estimate; ivw_se = th_mrivw@StdError; ivw_pval = th_mrivw@Pvalue
    th_egger = mr_mvegger(mrob)
    egger_b = th_egger@Estimate; egger_se = th_egger@StdError.Est; egger_pval = th_egger@Pvalue.Est
    th_robust = mvmr_robust(simdat$bxhat, simdat$byhat, simdat$seby, k.max = 1000, maxit.scale = 1000)
    robust_b = th_robust$coefficients; robust_se = th_robust$se; robust_pval = pnorm(-abs(robust_b/robust_se))*2
    th_lass = tryCatch({
        mvmr_lasso(simdat$bxhat, simdat$byhat, simdat$seby)
    },error = function(e){
        return(NULL)
    })
    if(is.null(th_lass)){lasso_b=lasso_se=lasso_pval=NA
        }else{
    lasso_b = th_lass$th_post; lasso_se = th_lass$se_post; lasso_pval = pnorm(-abs(lasso_b/lasso_se))*2}
    th_qr = mvmr_median(simdat$bxhat, simdat$sebx, simdat$byhat, simdat$seby, boot=TRUE)
    median_b = th_qr$coefficients; median_se = th_qr$se; median_pval = pnorm(-abs(median_b/median_se))*2

    var_name =  c(ls(pattern='_b$'),ls(pattern='_se$'),ls(pattern='_pval$'))
    for(vn in var_name){
        t1 = paste0("write.table(t(",vn,"),'",vn,"_",array_id,".txt',quote=F,row.names=F,append=T,col.names=F)")
        eval(parse(text=t1))
    }

} 
