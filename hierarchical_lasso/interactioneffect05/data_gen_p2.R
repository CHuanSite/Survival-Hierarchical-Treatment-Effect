projL1 <- function(r, c = 1) {
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
  #proj <- c(r[1], proj)
  return(proj)
}

#######generate fixed beta coef#########
beta.gen <- function(p = 10, m = 5, l = 3, seed = 1120){
  set.seed(seed)
  main.effect.index <- sort(sample(1 : p, m))
  inter.effect.index <- sort(sample(main.effect.index, l))
  
  beta.main <- rep(0, p)
  beta.main[main.effect.index] <- rnorm(length(main.effect.index), 0, 1)
  # beta.main[main.effect.index] <- runif(length(main.effect.index), -1, 1)
  tau <- rep(0, p)
  tau[inter.effect.index] <- runif(length(inter.effect.index), -0.5, 0.5)
  gamma <- beta.main * tau

  eta <- rnorm(p * (p - 1) / 2, 0, 1)
  t = 1
  for(i in 1 : (p - 1)) {
    for(j in (i + 1) : p){
      if(!(i %in% inter.effect.index && j %in% inter.effect.index)) {
        eta[t] = 0
      }
      t = t + 1
    }
  }
 
  return(list(main.effect.index = main.effect.index,
              inter.effect.index = inter.effect.index,
              beta = beta.main,
              gamma = gamma,
              tau = tau,
              eta = eta
	      ))
}
#######Data Generating#########
#####10 variables, 5 has main effects, 3 has interactions with treatment##########

gen.data <- function(n = 300, seed = 1032, p = 10, m = 5, 
                     l = 3, beta.trt = 1, beta.intercept = 0.5, 
                     sigma = 1) {
  ####Generate betas#####
  beta.dat <- beta.gen(p = p, m = m, l = l)
  
  ####Generate X's######
  set.seed(seed)
  vars.names <- paste0("X", 1:p)
  vars.mean <- 0 #rnorm(p, 0, 1) #sample(10:40, p)
  vars.sd <- 1#sample(5:10, p, replace = T)
  vars.dat <- NULL
  for(i in 1:p){
    vars.dat <- cbind(vars.dat, rnorm(n, vars.mean, vars.sd))
  }

  vars.dat.inter = NULL
  for(i in 1 : (p - 1)) {
    for(j in (i + 1) : p){
      vars.dat.inter <- cbind(vars.dat.inter, vars.dat[, i] * vars.dat[, j])
    }
  }

  main.effect.index <- beta.dat$main.effect.index
  #inter.effect.index <- beta.dat$inter.effect.index
  A <- sample(0:1, n, replace = T)
  beta.main <- beta.dat$beta.main
  #ksi.inter <- beta.dat$ksi.inter
  #beta.inter <- ksi.inter * beta.main * beta.trt
  beta.inter <- beta.dat$gamma
  #beta.inter <- projL1(beta.inter, beta.trt)
  inter.effect.index <- which(beta.inter != 0)
  eta.inter <- beta.dat$eta
  

  ####Generate y's
  y.mean <- beta.intercept + A * beta.trt + vars.dat %*% beta.main + (A * vars.dat) %*% beta.inter + vars.dat.inter %*% eta.inter
  y.response <- y.mean + rnorm(n, 0, sigma)
  data <- data.frame(cbind(y.response, A, vars.dat))
  names(data) <- c("y", "trt", vars.names)
  return(list(data = data, 
              main.effect.index = main.effect.index,
              inter.effect.index = inter.effect.index, 
              beta.main = beta.main, 
              #ksi.inter = ksi.inter, 
              beta.inter = beta.inter, 
              par = c(beta.intercept, beta.trt, beta.main, beta.inter)))
  #ksi.inter)))
}


######generate survival data##########
gen.data_surv <- function(n = 300, seed = 1032, p = 10, m = 5, 
                          l = 3, beta.trt = 1, beta.intercept = 0.5, 
                          sigma = 1, lambdaT = 1, lambdaC = 0.1){
  ####Generate betas#####
  beta.dat <- beta.gen(p = p, m = m, l = l)
  
  ####Generate X's######
  set.seed(seed)
  vars.names <- paste0("X", 1:p)
  vars.mean <- 0 #rnorm(p, 0, 1)#sample(10:40, p)
  vars.sd <- 1 #sample(5:10, p, replace = T)
  vars.dat <- NULL
  for(i in 1:p){
    vars.dat <- cbind(vars.dat, rnorm(n, vars.mean, vars.sd))
  }
 
  vars.dat.inter = NULL
  for(i in 1 : (p - 1)) {
    for(j in (i + 1) : p){
      vars.dat.inter <- cbind(vars.dat.inter, vars.dat[, i] * vars.dat[, j])
    }
  }


  main.effect.index <- beta.dat$main.effect.index
  #inter.effect.index <- beta.dat$inter.effect.index
  A <- sample(0:1, n, replace = T)
  beta.main <- beta.dat$beta
  #ksi.inter <- beta.dat$ksi.inter
  #beta.inter <- ksi.inter * beta.main * beta.trt
  beta.inter <- beta.dat$gamma
  #beta.inter <- projL1(beta.inter, beta.trt)
  inter.effect.index <- which(beta.inter != 0)
  eta.inter <- beta.dat$eta

  ####Generate survival data
  y.mean <- beta.intercept + A * beta.trt + vars.dat %*% beta.main + (A * vars.dat) %*% beta.inter
  #y.response <- y.mean + rnorm(n, 0, sigma)
  timeT <- rexp(n, rate = exp(y.mean) * lambdaT)
  timeC <- rexp(n, rate = lambdaC)
  time <- pmin(timeT, timeC)
  event <- timeT == time
  data <- data.frame(cbind(time, event, A, vars.dat))
  names(data) <- c("time", "event", "trt", vars.names)
  return(list(data = data, 
              main.effect.index = main.effect.index,
              inter.effect.index = inter.effect.index, 
              beta.main = beta.main, 
              #ksi.inter = ksi.inter, 
              beta.inter = beta.inter, 
              par = c(beta.intercept, beta.trt, beta.main, beta.inter)))
  #ksi.inter)))
}


