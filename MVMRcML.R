library(numDeriv)
l_i <- function(b_exp_i,b_out_i,Sig_inv_i,b_t_i,theta_t,r_vec_t_i){
  beta = c(b_exp_i,b_out_i) 
  bvec = c(b_t_i,sum(theta_t*b_t_i)+r_vec_t_i)
  l_i = 
    -1/2* t(beta-bvec) %*% Sig_inv_i %*% (beta-bvec)
  return(l_i)
}

fn_bi <- function(x,...){
  -l_i(b_t=x,...)
}

fn_ri <- function(x,...){
  -l_i(r_vec_t_i=x,...)
}

loglik <- function(b_exp,b_out,Sig_inv_l,b_t,theta_t,r_vec_t){
  l_i_vec = lapply(1:length(b_out),function(i) {
    l_i(b_exp[i,],b_out[i],Sig_inv_l[[i]],b_t[i,],theta_t,r_vec_t[i])})
  l = sum(unlist(l_i_vec))
  return(l)
}

pl <- function(x,b_exp_v,b_out_v,Sig_inv_v){
  k = ncol(b_exp_v)
  m_valid = length(b_out_v)
  pll = 0
  for(i in 1:m_valid){
    W = Sig_inv_v[[i]]
    b_exp_i = b_exp_v[i,]
    b_out_i = b_out_v[i]
    beta = c(b_exp_i,b_out_i)
    B = W[-(k+1),] %*% beta + c(W[(k+1),] %*% beta) * x
    A = W[-(k+1),-(k+1)] + W[-(k+1),k+1] %*% t(x) + x %*% t(W[k+1,-(k+1)]) + W[k+1,k+1] * x %*% t(x)
    bhat_xi = solve(A) %*% B
    b_i = c(bhat_xi,t(bhat_xi) %*% x)
    pll = pll  -1/2 * t(b_i) %*% W %*% b_i + t(beta) %*% W %*% b_i
  }
  return(-pll)

}


fn <- function(x,...){
  -loglik(theta_t=x,...)
}

MVcML_estimate <- function(b_exp,b_out,
                            Sig_inv_l,
                            K,initial_theta = rep(0,ncol(b_exp)),
                            initial_mu = matrix(0,nrow=nrow(b_exp),ncol=ncol(b_exp)),
                            maxit = 100)
{
  p = length(b_out)
  ### initialize
  theta = initial_theta
  theta_old = theta-1
  mu_vec = initial_mu
  ite_ind = 0
  k = length(theta)
  v_bg = rep(0,p)
  W_k.1k.1 = unlist(lapply(1:p,function(i){W_i = Sig_inv_l[[i]]; W_i[k+1,k+1]}))
  W_k.1 = Reduce(rbind,lapply(1:p,function(i){W_i = Sig_inv_l[[i]]; W_i[k+1,]}))
  pty = proc.time()
  while(sum(abs(theta-theta_old))>1e-4 & (ite_ind<maxit)){
    ite_ind = ite_ind + 1
    theta_old = theta
    if(K>0){
      A = diag(1/W_k.1k.1)
      beta_star = cbind(b_exp-mu_vec,b_out-(mu_vec %*% theta))
      B = diag(W_k.1 %*% t(beta_star))
      r_t = as.vector(A %*% B)
      v_importance = rep(NA,p)
      for(i in 1:p){
        l_r = 
          l_i(b_exp_i = b_exp[i,],b_out_i=b_out[i],Sig_inv_i = Sig_inv_l[[i]],
              theta_t=theta,r_vec_t_i = r_t[i],b_t_i = mu_vec[i,]);
        l_r0 = 
          l_i(b_exp_i = b_exp[i,],b_out_i=b_out[i],Sig_inv_i = Sig_inv_l[[i]],
              theta_t=theta,r_vec_t_i = 0,b_t_i = mu_vec[i,]);
        v_importance[i] = l_r-l_r0
      }
      nonzero_bg_ind = sort((order(v_importance,decreasing = T))[1:K])
      v_bg = rep(0,p)
      v_bg[nonzero_bg_ind] = r_t[nonzero_bg_ind]
    }else{
      v_bg = rep(0,p)
    }
    theta_cp = crossprod(t(theta))
    for(i in 1:p){
      A = 0
      B = 0
      W_i = Sig_inv_l[[i]]
      beta_star_i = c(b_exp[i,],b_out[i]) - c(rep(0,k),v_bg[i])
      #a = W_i[-(k+1),k+1] %*% t(theta)
      B =  W_i[-(k+1),] %*% (beta_star_i) + c(W_i[(k+1),] %*% (beta_star_i)) * (theta)
      A = W_i[-(k+1),-(k+1)] +  W_i[k+1,k+1]*theta_cp + W_i[-(k+1),k+1] %*% t(theta) + theta %*% t(W_i[k+1,-(k+1)]) 
      
      mu_vec[i,] = t(solve(A) %*% B)
    }
    
    
    A = 0 
    for(i in 1:p){
      W_i = Sig_inv_l[[i]]
      A = A + W_i[k+1,k+1] * crossprod(t(mu_vec[i,]))
    }

    B = 0
    for(i in 1:p){
      W_i = Sig_inv_l[[i]]
      beta_star_i = c(b_exp[i,]-mu_vec[i,],b_out[i]-v_bg[i])
      B = B + c(W_i[k+1,] %*% beta_star_i) * mu_vec[i,]
    }
    theta = as.vector(solve(A) %*% B)
  }
  return(list(theta = theta,
              b_vec = mu_vec,
              r_vec = v_bg))
}

MVcML_SdTheta <- function(b_exp,b_out,Sig_inv_l,theta,zero_ind,r_vec=NULL){
  if(!is.null(r_vec)){
      zero_ind = which(r_vec==0)
    }
  H = hessian(pl,x=theta,b_exp_v = b_exp[zero_ind,,drop=FALSE], b_out_v=b_out[zero_ind],Sig_inv_v =Sig_inv_l[zero_ind])
  se_theta = sqrt(diag(solve(H)))
  # ores = optim(theta,pl,method = 'BFGS',b_exp_v=b_exp[zero_ind,],b_out_v=b_out[zero_ind],Sig_inv_v= Sig_inv_l[zero_ind],
  #              hessian = TRUE)
  # sqrt(diag(solve(ores$hessian)))
  return(se_theta)
}


MVcML_estimate_random <- function(b_exp, b_out,
                                  Sig_inv_l,
                                  K,n,random_start = 0,
                                  maxit = 100)
{
  p = nrow(b_exp)
  k = ncol(b_exp)
  
  theta_v_RandomCandidate = NULL
  sd_v_RandomCandidate = NULL
  l_v_RandomCandidate = NULL
  invalid_RandomCandidate = NULL
  
  for(random_ind in 1:(1+random_start))
  {
    #      ptm <- proc.time()
    if(random_ind == 1)
    {
      initial_theta = rep(0,ncol(b_exp))
      initial_mu = matrix(0,nrow = nrow(b_exp),ncol=ncol(b_exp))
      #initial_mu = b_exp
    } else {
      initial_theta = runif(k,min = -0.5, max = 0.5)
#      initial_mu = rnorm(p,mean = b_exp,sd = se_exp)
      initial_mu = b_exp + rnorm(p*k,mean=0,sd=sqrt(1/n))
    }
    MLE_result =
      MVcML_estimate(b_exp,b_out,
                     Sig_inv_l,
                     K = K,initial_theta = initial_theta,
                     initial_mu = initial_mu,
                     maxit = maxit)
    
    Neg_l = -loglik(b_exp,b_out,Sig_inv_l,MLE_result$b_vec,MLE_result$theta,MLE_result$r_vec)
    

    theta_v_RandomCandidate = rbind(theta_v_RandomCandidate,MLE_result$theta)
    l_v_RandomCandidate = c(l_v_RandomCandidate,Neg_l)
    invalid_RandomCandidate = rbind(invalid_RandomCandidate,
                                    as.numeric(MLE_result$r_vec))
    
    #  print(proc.time() - ptm)
  }
  min_neg_l = which.min(l_v_RandomCandidate)
  
  theta_est = theta_v_RandomCandidate[min_neg_l,]
 # sd_est = sd_v_RandomCandidate[min_neg_l]
  l_est = l_v_RandomCandidate[min_neg_l]
  r_est = invalid_RandomCandidate[min_neg_l,]
  
  return(list(theta = theta_est,
#              se = sd_est,
              l = l_est,
              r_est = r_est
  )
  )
}

MVmr_cML <- function(b_exp,b_out,
                     Sig_inv_l,
                     K_vec = 0:(nrow(b_exp) - 2),
                     random_start = 0,
                     maxit = 100,
                     random_seed = 0,
                     n)
{
  if(random_seed)
  {
    set.seed(random_seed)
  }
  
  rand_theta = NULL
  rand_sd = NULL
  rand_l = NULL
  invalid_mat = NULL
  BIC_old = 99999
  for(K_value in K_vec)
  {
  #  print(K_value)
    rand_res = MVcML_estimate_random(b_exp = b_exp,
                                     b_out = b_out,
                                     Sig_inv_l = Sig_inv_l,
                                     K = K_value,
                                     n = n,
                                     random_start = random_start,
                                     maxit = maxit
                                     )
    rand_theta = rbind(rand_theta,rand_res$theta)
 #   rand_sd = c(rand_sd,rand_res$se)
    rand_l = c(rand_l,rand_res$l)
    invalid_mat = rbind(invalid_mat,rand_res$r_est)
  }
  
  ### get result
  theta_v = rand_theta
 # sd_v = rand_sd
  l_v = rand_l
  
  # cML-MA-BIC
  BIC_vec = log(n) * K_vec + 2 * l_v
  BIC_vec = BIC_vec - min(BIC_vec)
  weight_vec = exp(-1/2 * BIC_vec)
  weight_vec = weight_vec/sum(weight_vec)
  MA_BIC_theta = colSums(theta_v * weight_vec)
 # MA_BIC_se = sum(weight_vec * sqrt(sd_v^2 + (theta_v - MA_BIC_theta)^2),
#                  na.rm = TRUE)
#  MA_BIC_p = pnorm(-abs(MA_BIC_theta/MA_BIC_se))*2
  
  # cML-BIC
  BIC_vec = log(n) * K_vec + 2 * l_v
  BIC_vec = BIC_vec - min(BIC_vec)
  min_ind = which.min(BIC_vec)
  BIC_theta = theta_v[min_ind,]
 # BIC_se = MVcML_SdTheta(b_exp,b_out,
 #                        Sig_inv_l,
 #                        BIC_theta,
 #                        invalid_mat[min_ind,])
#  BIC_p = pnorm(-abs(BIC_theta/BIC_se))*2
  BIC_invalid = which(invalid_mat[min_ind,]!=0)
  
  
  return(list(MA_BIC_theta = MA_BIC_theta,
   #           MA_BIC_se = MA_BIC_se,
  #            MA_BIC_p = MA_BIC_p,
              BIC_theta = BIC_theta,
 #             BIC_se = BIC_se,
  #            BIC_p = BIC_p,
              BIC_invalid = BIC_invalid,
              l_vec = l_v,
              BIC_vec = log(n) * K_vec + 2 * l_v)
  )
}

MVmr_cML_DP <- function(b_exp,b_out,
                        Sig_l,
                        K_vec = 0:(length(b_exp) - 2),
                        random_start = 0,
                        random_start_pert = 0,
                        maxit = 100,
                        num_pert = 100,
                        random_seed = 0,
                        n)
{
  if(random_seed)
  {
    set.seed(random_seed)
  }
  theta_v = theta_MA_v = NULL
  Sig_inv_l = lapply(Sig_l, function(S) {solve(S)})
  cML_res = MVmr_cML(b_exp = b_exp, b_out = b_out, Sig_inv_l = Sig_inv_l,
                     random_start = random_start, random_seed = random_seed, n = n, maxit = maxit)
  p = length(b_out)
  for(pt_ind in 1:num_pert){
    epis = lapply(1:p,function(i){MASS::mvrnorm(1, mu = rep(0,ncol(b_exp)+1), Sigma = Sig_l[[i]])})
    epis = matrix(unlist(epis),ncol=ncol(b_exp)+1,byrow=T)
    b_exp_new = b_exp + epis[,1:ncol(b_exp)]
    b_out_new = b_out + epis[,ncol(b_exp)+1]
    cML_res_b = MVmr_cML(b_exp = b_exp_new, b_out = b_out_new, Sig_inv_l = Sig_inv_l,
                         random_start = random_start_pert, random_seed = random_seed, n = n, maxit = maxit)
 #   theta_MA_v = rbind(theta_MA_v,cML_res_b$MA_BIC_theta)
    theta_v = rbind(theta_v,cML_res_b$BIC_theta)
  }
  

  
  # cML-BIC-DP
  BIC_DP_theta = apply(theta_v,2,mean)
  BIC_DP_se = apply(theta_v,2,sd)
  BIC_DP_p = pnorm(-abs(BIC_DP_theta/BIC_DP_se))*2
  
  
  return(list(
              BIC_theta = cML_res$BIC_theta,
       #       BIC_se = cML_res$BIC_se,
  #            BIC_p = cML_res$BIC_p,
              BIC_invalid = cML_res$BIC_invalid,
              BIC_DP_theta = BIC_DP_theta,
              BIC_DP_se = BIC_DP_se,
              BIC_DP_p = BIC_DP_p
  ))
}







