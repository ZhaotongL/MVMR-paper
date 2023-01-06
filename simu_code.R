################################################################################
#Functions used in Sims_all_dat.R ad Sims_est_scr.R to produce the simulated
#data and estimates using the robust MVMR methods
################################################################################

Data_sim_ind = function(n, m, p, s, k, mua, mug, sa, sg, g, theta, bx_min, bx_max, SigX){
  v = sample(p, p * (1 -s))
  a = rnorm(p, mua, sa)
  a[v] = 0
  d = runif(p, 0, mug)
  d[v] = 0
  bx = sapply(1:k, function(j){runif(p, bx_min, bx_max)})
  G = sapply(1:p, function(i){rbinom(2 * n, 2, m)})
  U = G %*% d + rnorm(2 * n, 0, 1)
  X = G %*% bx + U %*% t(g) + mvrnorm(2 * n, rep(0, k), SigX)
  Y = X %*% theta + G %*% a + U + rnorm(2 * n, 0, 1)
  sampX = seq(1, n)
  sampY = seq((n+1),(2*n))
  bxhat = matrix(nrow = p, ncol = k)
  sebx = matrix(nrow = p, ncol = k)
  byhat = vector(length = p)
  seby = vector(length = p)
  for (i in 1:p){
    for (j in 1:k){
      ssx = sstat(X[sampX, j], G[sampX, i])
      bxhat[i, j] = ssx$bhat
      sebx[i, j] = ssx$se
    }
    ssy = sstat(Y[sampY], G[sampY, i])
    byhat[i] = ssy$bhat
    seby[i] = ssy$se
  }
  return(list("rho_mat" = cor(X[sampX,]), "bxhat" = bxhat, "sebx" = sebx, "byhat" = byhat, "seby" = seby))
}

Data_sim_irre1 = function(n, m, p, s, k, mua, mug, sa, sg, g, theta, bx_min, bx_max, SigX){
  v = sample(p, p * (1 -s))
  v1 = sample(p, p * (1 -s))
  a = d = rep(0,p)
  a[-v1] = rnorm(p*s, mua, sa)
  bx = sapply(1:k, function(j){runif(p, bx_min, bx_max)})
  bx[-v,] = 0
  G = sapply(1:p, function(i){rbinom(2 * n, 2, m)})
  U = G %*% d + rnorm(2 * n, 0, 1)
  X = G %*% bx + U %*% t(g) + mvrnorm(2 * n, rep(0, k), SigX)
  Y = X %*% theta + G %*% a + U + rnorm(2 * n, 0, 1)
  sampX = seq(1, n)
  sampY = seq((n+1),(2*n))
  bxhat = matrix(nrow = p, ncol = k)
  sebx = matrix(nrow = p, ncol = k)
  byhat = vector(length = p)
  seby = vector(length = p)
  for (i in 1:p){
    for (j in 1:k){
      ssx = sstat(X[sampX, j], G[sampX, i])
      bxhat[i, j] = ssx$bhat
      sebx[i, j] = ssx$se
    }
    ssy = sstat(Y[sampY], G[sampY, i])
    byhat[i] = ssy$bhat
    seby[i] = ssy$se
  }
  simdat = list("rho_mat" = cor(X[sampX,]), "bxhat" = bxhat, "sebx" = sebx, "byhat" = byhat, "seby" = seby)
  return(simdat)
}


sstat = function(Y, X, intercept = TRUE){
  n = length(Y)
  if (intercept == TRUE){
    xx = cbind(rep(1, n), X)
  }
  else {xx = X}
  mod = lm.fit(xx, Y)
  bhat= c(mod$coefficients[2])
  s = t(mod$residuals) %*% (mod$residuals) / (mod$df.residual)
  se = sqrt((c(s) * solve(t(xx) %*% xx))[2,2])
  return(list("bhat" = bhat, "se" = se))
}

Est_sim_meta = function(M, D){
  sapply(1:M, function(j){
    Dsum = D[[j]]
    mrob = mr_mvinput(bx = Dsum$bxhat, bxse = Dsum$sebx, by = Dsum$byhat, byse = Dsum$seby)
    th_mrivw = mr_mvivw(mrob)
    th_egger = mr_mvegger(mrob)
    th_robust = mvmr_robust(Dsum$bxhat, Dsum$byhat, Dsum$seby, k.max = 1000, maxit.scale = 1000)
    c(th_mrivw$Estimate, th_mrivw$StdError, th_egger$Estimate, th_egger$StdError.Est,
      th_robust$coefficients, th_robust$se)
  })
}

Est_sim_mrpresso = function(M, D){
  parLapply(cl, 1:M, function(j){
    Dsum = D[[j]]
    k = dim(Dsum$bxhat)[2]
    dat = data.frame(bx = Dsum$bxhat, sebx = Dsum$sebx, by = Dsum$byhat, seby = Dsum$seby)
    A = mr_presso(BetaOutcome = "by", BetaExposure = names(dat)[1:k], SdOutcome = "seby",
                  SdExposure = names(dat)[(k+1):(2*k)], OUTLIERtest = TRUE,
                  DISTORTIONtest = TRUE, data = dat, NbDistribution = 2000,
                  SignifThreshold = 0.05)
    return(A)
  })
}

Est_sim_med = function(M, D){
  parSapply(cl, 1:M, function(j){
    Dsum = D[[j]]
    th_qr = mvmr_qr(Dsum$bxhat, Dsum$sebx, Dsum$byhat, Dsum$seby, boot = TRUE, boot_it = 1000)
    c(th_qr$coefficients, th_qr$se)
  })
}

Est_sim_med_CI = function(M, D){
  parSapply(cl, 1:M, function(j){
    Dsum = D[[j]]
    th_qr = mvmr_qr_CI(Dsum$bxhat, Dsum$byhat, Dsum$seby, boot = TRUE, boot_it = 1000)
    c(th_qr$coefficients, th_qr$CIlow_boot, th_qr$CIupp_boot, th_qr$CIlow_rinv, th_qr$CIupp_rinv)
  })
}

Est_sim_lass = function(M, D){
  parSapply(cl, 1:M, function(j){
    Dsum = D[[j]]
    th_lass = mvmr_lass(Dsum$bxhat, Dsum$byhat, Dsum$seby)
    c(th_lass$th_post, th_lass$se_post)
  })
}

mvmr_qr_CI = function(bx, by, seby, boot = FALSE, boot_it = 1000){
  qr_mod = rq(by ~ bx - 1, weights = seby^-2)
  if (boot == TRUE){
    a = sapply(1:boot_it, function(i){
      p = length(by)
      b = sample(p, replace = TRUE)
      rq(by[b] ~ bx[b, ] - 1, weights = seby[b]^-2)$coefficients[1]
    })
    return(list("coefficients" = qr_mod$coefficients[1],
                "CIlow_boot" = sort(a)[round(boot_it/20,1)], "CIupp_boot" = sort(a)[round(19*boot_it/20,1)],
                "CIlow_rinv" = summary(qr_mod)$coefficients[1, 2], "CIupp_rinv" = summary(qr_mod)$coefficients[1, 3]))
  } else {
    return(list("coefficients" = qr_mod$coefficients[1],
                "CIlow_rinv" = summary(qr_mod)$coefficients[1, 2], "CIupp_rinv" = summary(qr_mod)$coefficients[1, 3]))
  }
}


################################################################################
#Functions to implement the MVMR-robust, MVMR-median and MVMR-lasso methods
################################################################################

library(MASS)
library(glmnet)
library(quantreg)
library(robustbase)

mvmr_med_boot = function(bx, sebx, by, seby, N){
  est = sapply(1:N, function(i){
    p = length(by)
    k = dim(bx)[2]
    Sx = lapply(1:p, function(j){diag(sebx[j, ]^2)})
    bxboot = sapply(1:p, function(j){mvrnorm(1, bx[j, ], Sx[[j]])})
    bxboot = t(bxboot)
    byboot = rnorm(p, by, seby)
    rq(byboot ~ bxboot - 1, weights = seby^-2)$coefficients
  })
  apply(est, 1, sd)
}

mvmr_median = function(bx, sebx, by, seby, boot = FALSE, boot_it = 1000){
  qr_mod = rq(by ~ bx - 1, weights = seby^-2)
  if (boot == TRUE){
    boot_se = mvmr_med_boot(bx, sebx, by, seby, boot_it)
    return(list("coefficients" = qr_mod$coefficients, "se" = boot_se))
  } else {
    return(list("coefficients" = qr_mod$coefficients))
  }
}

cv.mvmr_lasso = function(bx, by, seby){
  p = dim(bx)[1]
  k = dim(bx)[2]
  S = diag(seby^-2)
  b = S^(1/2) %*% bx
  Pb = b %*% solve(t(b) %*% b, t(b))
  xlas = (diag(p) - Pb) %*% S^(1/2)
  ylas = (diag(p) - Pb) %*% S^(1/2) %*% by
  alas = glmnet(xlas, ylas, intercept = FALSE)
  lamseq = sort(alas$lambda)
  lamlen = length(lamseq)
  rse = sapply(1:lamlen, function(j){
    av = which(alas$beta[, (lamlen - j + 1)] == 0)
    mod = lm.fit(as.matrix(S[av, av]^(1/2) %*% bx[av, ]), S[av, av]^(1/2) %*% by[av])
    c(sqrt(t(mod$residuals) %*% (mod$residuals) / (mod$df.residual)), length(av))
  })
  rse_inc = rse[1, 2:lamlen] - rse[1, 1:(lamlen-1)]
  het = which(rse[1, 2:lamlen] > 1 & rse_inc > ((qchisq(0.95, 1) / rse[2, 2:lamlen])))
  if (length(het) == 0){
    lam_pos = 1
  } else {
    lam_pos = min(het)
  }
  num_valid = rev(sapply(1:lamlen, function(j){sum(alas$beta[, j]==0)}))
  min_lam_pos = min(which(num_valid > k))
  if (lam_pos < min_lam_pos){lam_pos = min_lam_pos}
  return(list(fit = alas$beta[, (lamlen - lam_pos + 1)], lambda = lamseq[lam_pos]))
}

mvmr_lasso = function(bx, by, seby){
  p = dim(as.matrix(bx))[1]
  k = dim(as.matrix(bx))[2]
  S = diag(seby^-2)
  sn = sign(bx[, 1])
  bx_or = bx * sn
  by_or = by * sn
  cv.alas = cv.mvmr_lasso(bx_or, by_or, seby)
  a1 = cv.alas$fit
  e = by_or - a1
  thest = solve(t(bx_or) %*% S %*% bx_or, t(bx_or) %*% S %*% e)
  v = which(a1==0)
  mvmr_mod = mr_mvivw(mr_mvinput(bx = bx_or[v, ], bxse = bx_or[v, ],
                                 by = by_or[v], byse = seby[v]))
  th_post = mvmr_mod$Estimate
  se_post = mvmr_mod$StdError
  return(list(thest = thest, a = a1, lambda = cv.alas$lambda,
              th_post = th_post, se_post = se_post))
}

mvmr_robust = function(bx, by, seby, k.max = 500, maxit.scale = 500){
  robmod = lmrob(by ~ bx - 1, weights = seby^-2, k.max = k.max,
                 maxit.scale = maxit.scale)
  coefficients = summary(robmod)$coef[, 1]
  se = summary(robmod)$coef[, 2] / min(summary(robmod)$sigma, 1)
  return(list("coefficients" = coefficients, "se" = se))
}
