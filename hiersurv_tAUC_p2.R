rm(list = ls())
#-------Source data and load packages needed---------#
# source("./new_work/Simulations/data_gen_p2.R")
source("/users/hchen/new_work/Simulations/data_gen_p2.R")
library(BB)
library(glmnet)
library(BiocParallel)
library(survival)
library(nloptr)
library(timeROC)
#----------------------------------------------------#
k <- as.numeric(Sys.getenv("SGE_TASK_ID"))
#-----------projection function L1 penalty-----------#
proj <- function(r, c = 1) {
  # (c) Xiaojing Ye
  # Translated to R by Ravi Varadhan, Johns Hopkins University
  # August 8, 2012
  tr <- r
  m <- length(tr)
  u <- abs(tr)
  if (sum(u) <= c) return(tr) 
  bget <- FALSE
  sr <- sort(u, decreasing=TRUE)
  tmpsum <- 0
  for (i in 1:(m-1)){
    tmpsum <- tmpsum + sr[i]
    tmax <- (tmpsum - c) / i
    if (tmax >= sr[i+1]) {bget <- TRUE; break}
  }
  if (!bget) tmax <- (tmpsum + sr[m] - c)/m
  proj <- sign(tr) * pmax(0, u-tmax)
  return(proj)
}


#-----------fitting parameterization form 1 - 5 ---------------#
hierLasso_surv <- function(data, s, form = 1, tol = 1e-3) {
  A <- data[[3]] #---treatment----#
  X <- as.matrix(data[, -c(1, 2, 3)]) #---covariate---#
  n <- nrow(X); p <- ncol(X)
  AX <- tcrossprod(A, rep(1,p)) * X
  X_ <- cbind(A, X, AX)
  p_ <- ncol(X_)
  X2 <- cbind(X_, -X_)
  expand.data <- cbind(data, AX)
  names(expand.data)[(p + 4) : (p * 2 + 3)] <- paste0("AX", 1 : p)
  expand.data <- expand.data[order(expand.data$time, decreasing = T), ]
  #-----objective functions------#
  obj1 <- function(r) {
    ba <- r[1]
    bx <- r[2 : (p + 1)]
    d <- r[(p + 2) : (p * 2 + 1)]
    par <- c(ba, bx, bx * ba * d)
    data <- expand.data
    Pt1 <- data.matrix(data[, -c(1 : 2)]) %*% par
    Pt2 <- log(cumsum(exp(Pt1)))
    PL <- sum(data$event * (Pt1 - Pt2))
    return(-PL)
  }
  
  obj2 <- function(r) {
    ba <- r[1]
    bx <- r[2 : (p + 1)]
    d <- r[(p + 2) : (p * 2 + 1)]
    par <- c(ba, bx, bx * d)
    data <- expand.data
    Pt1 <- data.matrix(data[, -c(1 : 2)]) %*% par
    Pt2 <- log(cumsum(exp(Pt1)))
    PL <- sum(data$event * (Pt1 - Pt2))
    return(-PL)
  }
  
  obj3 <- function(w) {
    wp <- w[1:p_]
    wm <- w[(p_ + 1):(2 * p_)]
    par <- wp - wm
    data <- expand.data
    Pt1 <- data.matrix(data[, -c(1 : 2)]) %*% par
    Pt2 <- log(cumsum(exp(Pt1)))
    PL <- sum(data$event * (Pt1 - Pt2))
    return(-PL)
  }
  
  obj4 <- obj3
  
  obj5 <- function(r) {
    data <- expand.data
    Pt1 <- data.matrix(data[, -c(1 : 2)]) %*% r
    Pt2 <- log(cumsum(exp(Pt1)))
    PL <- sum(data$event * (Pt1 - Pt2))
    return(-PL)
  }
  
  #-----gradient functions------#
  GR1 <- function(r) {
    ba <- r[1]
    bx <- r[2 : ((length(r) - 1) / 2 + 1)]
    d <- r[((length(r) - 1) / 2 + 2) : length(r)]
    par <- c(ba, bx, bx * ba * d)
    data <- expand.data
    Pt1 <- data.matrix(data[, -c(1 : 2)]) %*% par
    Pt2_denom <- cumsum(exp(Pt1))
    Pt2_numer <- cumsum(exp(Pt1) * data[, -c(1 : 2)])
    Gmat1 <- c(1, rep(0, length(bx)), bx * d)
    Gmat2 <- cbind(rep(0, length(bx)), diag(length(bx)), diag(c(ba * d)))
    Gmat3 <- cbind(rep(0, length(bx)), matrix(0, length(bx), length(bx)), diag(c(ba * bx)))
    GMAT <- rbind(Gmat1, Gmat2, Gmat3)
    GR <- - GMAT %*% colSums(data$event * (data[, -c(1 : 2)] - Pt2_numer / Pt2_denom))
  }
  
  GR2 <- function(r) {
    ba <- r[1]
    bx <- r[2 : ((length(r) - 1) / 2 + 1)]
    d <- r[((length(r) - 1) / 2 + 2) : length(r)]
    par <- c(ba, bx, bx * d)
    data <- expand.data
    Pt1 <- data.matrix(data[, -c(1 : 2)]) %*% par
    Pt2_denom <- cumsum(exp(Pt1))
    Pt2_numer <- cumsum(exp(Pt1) * data[, -c(1 : 2)])
    Gmat1 <- c(1, rep(0, length(bx)), rep(0, length(d)))
    Gmat2 <- cbind(rep(0, length(bx)), diag(length(bx)), diag(d))
    Gmat3 <- cbind(rep(0, length(bx)), matrix(0, length(bx), length(bx)), diag(bx))
    GMAT <- rbind(Gmat1, Gmat2, Gmat3)
    GR <- - GMAT %*% colSums(data$event * (data[, -c(1 : 2)] - Pt2_numer / Pt2_denom))
  }
  
  GR3 <- function(w) {
    data <- cbind(data[, c(1, 2)], X2)
    names(data) <- c("time", "event", "A", paste0("X", 1 : p), paste0("AX", 1 : p), 
                     "nA", paste0("nX", 1 : p), paste0("nAX", 1 : p))
    data <- data[order(data$time, decreasing = T), ]
    Pt1 <- data.matrix(data[, -c(1 : 2)]) %*% w
    Pt2_denom <- cumsum(exp(Pt1))
    Pt2_numer <- cumsum(exp(Pt1) * data[, -c(1 : 2)])
    GR <- -colSums(data$event * (data[, -c(1 : 2)] - Pt2_numer / Pt2_denom))
    return(GR)
  }
  
  GR4 <- GR3
  
  GR5 <- function(r) {
    data <- expand.data
    Pt1 <- data.matrix(data[, -c(1 : 2)]) %*% r
    Pt2_denom <- cumsum(exp(Pt1))
    Pt2_numer <- cumsum(exp(Pt1) * data[, -c(1 : 2)])
    GR <- -colSums(data$event * (data[, -c(1 : 2)] - Pt2_numer / Pt2_denom))
    return(GR)
  }
  
  w0 <- rep(0, 2 * p_)  #----initial values for form 3 and 4------#
  
  p0 <- runif((2 * p + 1), 0, 1) #----initial values for form 1 and 2----#
  p0 <- proj(p0, c = s)
  
  #-----------------form 1-----------------------------------#
  form1 <- function() {
    ans <- spg(par = p0, fn = obj1, gr = GR1,
               project = proj, projectArgs = list(c = s), 
               control = list(maxit = 1100, maxfeval = 2000, 
                              #checkGrad = TRUE,
                              ftol = tol, gtol = tol))$par
    ans[(p + 2) : (2 * p + 1)] <- ans[(p + 2) : (2 * p + 1)] * ans[2 : (p + 1)] * ans[1]
    ans <- as.numeric(ans)
    return(ans)
  }
  
  #-----------------form 2-----------------------------------#
  form2 <- function() {
    ans <- spg(par = p0, fn = obj2, gr = GR2,
               project = proj, projectArgs = list(c = s), 
               control = list(maxit = 1100, maxfeval = 2000, 
                              #checkGrad = TRUE,
                              ftol = tol, gtol = tol))$par
    ans[(p + 2) : (2 * p + 1)] <- ans[(p + 2) : (2 * p + 1)] * ans[2 : (p + 1)]
    ans <- as.numeric(ans)
    return(ans)
  }
  
  #-----------------form 3-----------------------------------#
  form3 <- function() {
    Amat1 <- t(rep(1, 2 * p_))
    Amat2 <- cbind(rep(0, p), -diag(p), diag(p)) 
    Amat2 <- cbind(Amat2, Amat2)
    Amat3 <- cbind(rep(-1, p), matrix(0, p, p), diag(p))
    Amat3 <- cbind(Amat3, Amat3)
    Amat <- rbind(Amat1, Amat2, Amat3)
    B <- c(s, rep(0,p), rep(0,p))
    
    eval_g_ineq <- function(w) Amat %*% w - B
    eval_jac_g_ineq <- function(w) Amat
    
    w <- nloptr(x0 = w0, eval_f = obj3, eval_grad_f = GR3, 
                eval_g_ineq = eval_g_ineq, eval_jac_g_ineq = eval_jac_g_ineq,
                lb = rep(0, 2 * p_),
                opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = tol) )$solution
    
    # Amat <- -1 * Amat
    # Amat <- rbind(Amat, diag(2 * p_))
    # B <- -1 * B
    # B <- c(B, rep(0, 2 * p_))
    # w <- spg(par = w0, fn = obj3, gr = GR3,
    #          project = "projectLinear", 
    #          projectArgs = list(A = Amat, b = B, meq = 0), 
    #          control = list(maxit = 1100, maxfeval = 2000, 
    #                         checkGrad = TRUE,
    #                         ftol = tol, gtol = tol))$par
    wp <- w[1 : p_]
    wm <- w[(p_ + 1):(2 * p_)]
    wout <- as.numeric(wp - wm)
    return(wout)
  }
  
  #-----------------form 4-----------------------------------#
  form4 <- function() {
    Amat1 <- t(rep(1, 2 * p_))
    Amat2 <- cbind(rep(0, p), -diag(p), diag(p)) 
    Amat2 <- cbind(Amat2, Amat2)
    Amat3 <- t(rep(1, 2 * p_))
    Amat3[1, c(2 : (p + 1))] <- Amat3[1, (p_ + 2) : (p_ + p + 1)] <- 0
    Amat3[1, 1] <- Amat3[1, (p_ + 1)] <- -1
    Amat <- rbind(Amat1, Amat2, Amat3)
    B <- c(s, rep(0,p), 0)
    
    eval_g_ineq <- function(w) Amat %*% w - B
    eval_jac_g_ineq <- function(w) Amat
    
    w <- nloptr(x0 = w0, eval_f = obj4, eval_grad_f = GR4, 
                eval_g_ineq = eval_g_ineq, eval_jac_g_ineq = eval_jac_g_ineq,
                lb = rep(0, 2 * p_),
                opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = tol) )$solution
    
    # Amat <- -1 * Amat
    # Amat <- rbind(Amat, diag(2 * p_))
    # B <- -1 * B
    # B <- c(B, rep(0, 2 * p_))
    # w <- spg(par = w0, fn = obj4, gr = GR4,
    #          project = "projectLinear", 
    #          projectArgs = list(A = Amat, b = B, meq = 0), 
    #          control = list(maxit = 1100, maxfeval = 2000, 
    #                         checkGrad = TRUE,
    #                         ftol = tol, gtol = tol))$par
    wp <- w[1 : p_]
    wm <- w[(p_ + 1) : (2 * p_)]
    wout <- as.numeric(wp - wm)
    return(wout)
  }
  
  #-----------------form 5-----------------------------------#
  form5 <- function() {
    # fit.glmnet <- cv.glmnet(x = X_, y = Surv(data$time, data$event), 
    #                         family = "cox")
    # lambda.min <- fit.glmnet$lambda.min
    # coef.fulllasso <- as.numeric(predict(fit.glmnet, type = "coefficient", 
    #                                      s = lambda.min))
    # return(list(par.lasso = coef.fulllasso))
    ans <- spg(par = p0, fn = obj5, gr = GR5,
               project = proj, projectArgs = list(c = s), 
               control = list(maxit = 1100, maxfeval = 2000, 
                              #checkGrad = TRUE, 
                              ftol = tol, gtol = tol))$par
    ans <- as.numeric(ans)
    return(ans)
  }
  
  #---form 5 is regular lasso----#
  w_op <- switch(form, form1(), form2(), form3(), form4(), form5())
  
  return(ifelse(abs(w_op) > tol, w_op, 0))
}
#--------------------------------------------------------------#

#-----------cross validation-----------------------------------#
hierLassocv_surv <- function(data, s_seq, form = 1, k = 10, tol = 1e-3, 
                             cv.method = c("concordance", "tAUC"), quan_seq = NULL) {
  s_seq <- sort(s_seq)
  ns <- length(s_seq)  
  n <- nrow(data)
  
  samples <- split(sample(n), c(rep(seq(k), each = n %/% k), rep(1, n%%k)) )
  cv_out <- matrix(nrow = k, ncol = ns)
  colnames(cv_out) <- s_seq
  
  # CV prediction error function
  cv_pred <- function(data_cv, cv_fit) {
    A <- data_cv[[3]]
    X <- as.matrix(data_cv[,-c(1, 2, 3)])
    p <- ncol(X); n1 <- nrow(X)
    AX <- tcrossprod(A, rep(1, p)) * X
    X_ <- cbind(A, X, AX)
    pred.risk <- X_ %*% cv_fit
    if(cv.method == "concordance") {
      concord <- as.numeric(survConcordance(Surv(time, event) ~ pred.risk, data = data_cv)[[1]])
      return(concord)
    }else {
      AUC_obj <- timeROC(T = data_cv$time,
                         delta = data_cv$event, 
                         marker = pred.risk,
                         cause = 1, weighting = "marginal",
                         times = quan_seq)
      AUC_val <- mean(AUC_obj$AUC)
      return(AUC_val)
    }
    
  }
  
  for (j in seq(ns)) {
    for (i in seq(k)) {
      data_ <- data[unlist(samples[-i]),]
      data_cv <- data[samples[[i]],]
      cv_fit <- hierLasso_surv(data_, s = s_seq[j], form = form, tol = tol)
      cv_out[i,j] <- cv_pred(data_cv, cv_fit)
    }
  }
  
  cv_res <- apply(cv_out, 2, mean)
  s_cv <- as.numeric(names(which.max(cv_res)))
  par <- hierLasso_surv(data, s = s_cv, form = form, tol = tol)
  
  return (list(par = par, 
               s_cv = s_cv))
}
#--------------------------------------------------------------#

#---------------Simulations to compare formulations 1 2 3 4 5----------------------#
comp.model <- function(n = 10500, n1 = 500, p = 50, m = 25, seed = 1031, 
                       l = 5, beta.trt = 1, beta.intercept = 0.5, 
                       sigma = 1, cv.method = c("concordance", "tAUC")) {
  
  list.data <- gen.data_surv(n = n, p = p, m = m, 
                             l = l, beta.trt = beta.trt, 
                             beta.intercept = beta.intercept, 
                             sigma = sigma, seed = seed)
  obs.data <- list.data$data[1 : n1, ]
  test.data <- list.data$data[(n1 + 1) : n, ]
  truepar <- c(beta.trt, list.data$beta.main, list.data$beta.inter)
  quan_seq <- quantile(list.data$data$time, probs = seq(0.3, 0.7, 0.1))
  
  set.seed(seed)
  #-----------determine s sequence---------------#
  A <- obs.data[, 1]
  X <- obs.data[, 2 : (p + 1)]
  AX <- A * X
  obs.data_expand <- cbind(obs.data, AX)
  names(obs.data_expand)[(p + 4) : (p * 2 + 3)] <- paste0("AX", 1 : p)
  fit.surv <- coxph(Surv(time, event) ~ ., data = obs.data_expand)
  s.ub <- sum(abs(as.numeric(coef(fit.surv))))
  s_seq <- seq(as.numeric(coef(fit.surv)[1]), s.ub, length.out = 20)
  
  par.dat <- matrix(NA, nrow = (2 * p + 1), ncol = 5)
  main_terms.dat <- matrix(NA, nrow = p, ncol = 5)
  inter_terms.dat <- matrix(NA, nrow = p, ncol = 5)
  
  for(i in 1 : 5) {
    print(i)
    par.dat[, i] <- hierLassocv_surv(data = obs.data, s_seq = s_seq, form = i, k = 10, 
                                     cv.method = cv.method, quan_seq = quan_seq)$par
    main_terms.dat[, i] <- par.dat[2 : (p + 1), i]
    inter_terms.dat[, i] <- par.dat[(p + 2) : (2 * p + 1), i]
  }
  colnames(par.dat) <- paste("Form", 1 : ncol(par.dat))
  #---------------compute performance metrics-----------------#
  A <- test.data[[2]]
  X <- as.matrix(test.data[, -c(1, 2, 3)])
  p <- ncol(X)
  AX <- tcrossprod(A, rep(1,p)) * X
  #test.data_expand <- cbind(test.data, AX)
  #names(test.data_expand)[(p + 4) : (p * 2 + 3)] <- paste0("AX", 1 : p)
  X_ <- cbind(A, X, AX)
  
  if(cv.method == "concordance") {
    #-----concordance index for testset-----#
    out_val <- apply(par.dat, 2, function(a) {
      #-coxph(Surv(time, event) ~ ., data = test.data_expand, init = a, 
      #       control = list('iter.max' = 0))$loglik[2]
      pred.risk <- X_ %*% a
      concord <- as.numeric(survConcordance(Surv(time, event) ~ pred.risk, data = test.data)[[1]])
      return(concord)
      
    })
  }else {
    #-----tAUC for testset-----#
    out_val <- apply(par.dat, 2, function(a) {
      pred.risk <- X_ %*% a
      AUC_obj <- timeROC(T = test.data$time,
                         delta = test.data$event, 
                         marker = pred.risk,
                         cause = 1, weighting = "marginal",
                         times = quan_seq)
      mean(AUC_obj$AUC)
    })
  }
  
  
  #--------------Compare sensitivity and specificity for main and interaction terms--------------#
  coef.main <- list.data$beta.main
  coef.inter <- list.data$beta.inter
  
  #------for main terms-------#
  sens.main <- apply(main_terms.dat, 2, function(a){
    mean(sapply(which(coef.main != 0), function(o){
      as.numeric(a[o] != 0)       
    }))})
  
  spec.main <-  apply(main_terms.dat, 2, function(a){
    mean(sapply(which(coef.main == 0), function(o){
      as.numeric(a[o] == 0)       
    }))})
  
  count.sens.main <- apply(main_terms.dat, 2, function(a){
    sum(sapply(which(coef.main != 0), function(o){
      as.numeric(a[o] != 0)       
    }))})
  
  count.spec.main <- apply(main_terms.dat, 2, function(a){
    sum(sapply(which(coef.main == 0), function(o){
      as.numeric(a[o] == 0)       
    }))})
  
  #------for interaction terms------#
  sens.inter <- apply(inter_terms.dat, 2, function(a){
    mean(sapply(which(coef.inter != 0), function(o){
      as.numeric(a[o] != 0)       
    }))})
  
  spec.inter <- apply(inter_terms.dat, 2, function(a){
    mean(sapply(which(coef.inter == 0), function(o){
      as.numeric(a[o] == 0)       
    }))})
  
  count.sens.inter <- apply(inter_terms.dat, 2, function(a){
    sum(sapply(which(coef.inter != 0), function(o){
      as.numeric(a[o] != 0)       
    }))})
  
  count.spec.inter <- apply(inter_terms.dat, 2, function(a){
    sum(sapply(which(coef.inter == 0), function(o){
      as.numeric(a[o] == 0)       
    }))})
  
  return(list(par.dat = par.dat, out_val = out_val, truepar = truepar, 
              sens.main = sens.main, spec.main = spec.main, 
              count.sens.main = count.sens.main, count.spec.main = count.spec.main,
              sens.inter = sens.inter, spec.inter = spec.inter, 
              count.sens.inter = count.sens.inter, count.spec.inter = count.spec.inter))
  
}


comp.result <- comp.model(n = 10500, n1 = 500, p = 50, seed = k,
                          m = 25, l = 5, beta.trt = 1, beta.intercept = 0.5,
                          sigma = 1, cv.method = "tAUC")
save(comp.result, file = paste("/users/hchen/new_work/Simulations/p2_tAUC_output/", "output_tAUC", k, ".RData", sep = ""))




