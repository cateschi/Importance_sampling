rm(list=ls())     # Clean memory
graphics.off()    # Close graphs
cat("\014")       # Clear Console



#### Packages ####

library(ucminf)



#### Functions ####


# Kalman filter and smoother
KFS <- function(y, Tmatrix, Z, H, Q, nstates, d, P10, outofsample){
  
  len <- length(y)
  x10 <- rep(0,nstates)
  Pttm1 <- lapply(seq_len(len+1), function(X) matrix(0,nstates,nstates))
  Ptt <- lapply(seq_len(len), function(X) matrix(0,nstates,nstates))
  P10 <- diag(P10,nstates,nstates)
  Pttm1[[1]] <- P10
  xtt <- matrix(0,nstates,(len))
  xttm1 <- matrix(0,nstates,(len+1))
  xttm1[,1] <- x10
  
  logl_t <- rep(NA,len)
  
  #Bulild Z:
  Z <- Z
  
  #Build T:y
  Tmatrix <- Tmatrix
  
  #initialization of loglikelihood
  logl <- 0
  
  #Start of KF recursions
  for (i in 1:len){
    epshatoutofsample <- y[i] - Z%*%xttm1[,i]
    Fmatrix <- Z%*%Pttm1[[i]]%*%t(Z) + H[i,i]
    
    Fmatrix.inv <- 1/Fmatrix
    Kg <- Tmatrix%*%Pttm1[[i]]%*%t(Z)%*%Fmatrix.inv 
    xtt[,i] <- xttm1[,i]+Pttm1[[i]]%*%t(Z)%*%Fmatrix.inv%*%epshatoutofsample 
    epshatinsample <- y[i]-Z%*%xtt[,i] 
    Ptt[[i]] <- Pttm1[[i]]-Pttm1[[i]]%*%t(Z)%*%Fmatrix.inv%*%Z%*%Pttm1[[i]] 
    
    Q <- Q     
    
    Pttm1[[i+1]] <- Tmatrix%*%Pttm1[[i]]%*%t(Tmatrix-Kg%*%Z)+Q 
    xttm1[,i+1] <- Tmatrix%*%xttm1[,i] + Kg%*%epshatoutofsample
    
    #The optimization criterion
    if (outofsample) {
      if (i <= d ){
        logl <- logl - 1/2*log(2*pi)
        if ((NaN %in% logl)==T){
          logl<- -P10[1]
        }
      } else if (i > d ){ 
        logl <- logl - 1/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatoutofsample)%*%Fmatrix.inv%*%epshatoutofsample
        if ((NaN %in% logl)==T){
          logl<- -P10[1]
        }
      }
    } else {
      if (i <= d ){
        logl <- logl - 1/2*log(2*pi)
        if ((NaN %in% logl)==T){
          logl<- -P10[1]
        }
      } else if (i > d ){ 
        logl <- logl - 1/2*log(2*pi) - 1/2*log(det(Fmatrix)) - 1/2*t(epshatinsample)%*%Fmatrix.inv%*%epshatinsample
        if ((NaN %in% logl)==T){
          logl<- -P10[1]
        }
      }
    }
    
    logl_t[i] <- logl
  } 
  
  # KF smoother
  xtm1T <- matrix(0,nstates,len) #smoothed state vector
  Ptm1T <- lapply(seq_len(len), function(X) matrix(0,nstates,nstates)) #smoothed variance of the states vector
  for (i in len:1){
    if (i==len) {
      xtm1T[,i] <- xtt[,i]
      Ptm1T[[i]] <- Ptt[[i]]
    } else{
      xtm1T[,i] <- xtt[,i] + Ptt[[i]]%*%t(Tmatrix)%*%solve(Pttm1[[i+1]])%*%(xtm1T[,i+1] - xttm1[,i+1])
      Ptm1T[[i]] <- Ptt[[i]] + Ptt[[i]]%*%t(Tmatrix)%*%solve(Pttm1[[i+1]])%*%(Ptm1T[[i+1]] - Pttm1[[i+1]])%*%t(Ptt[[i]]%*%t(Tmatrix)%*%solve(Pttm1[[i+1]]))
    }
  }
  
  return(list(logl=-logl, logl_t=logl_t, xttm1=xttm1, Pttm1=Pttm1, xtm1T=xtm1T, Ptm1T=Ptm1T))
}      


# Newton-Raphson algorithm (that uses the KFS to update g)
NR <- function(y, g, Psi, max.iter, Tmatrix, Z, Q, nstates, d, P10, outofsample){
  
  for (j in 1:max.iter){
  
    A <- - svd(diag(c(-exp(g))))$v%*%diag(1/svd(diag(c(-exp(g))))$d)%*%t(svd(diag(c(-exp(g))))$u)
    
    z <- g + A%*%(y - exp(g))
    
    Ainv <- svd(A)$v%*%diag(1/svd(A)$d)%*%t(svd(A)$u)
    
    Psinv <- svd(Psi)$v%*%diag(1/svd(Psi)$d)%*%t(svd(Psi)$u) 
    
    invMat <- svd(Psinv + Ainv)$v%*%diag(1/svd(Psinv + Ainv)$d)%*%t(svd(Psinv + Ainv)$u) 
    
    g <- t(KFS(y=g + A%*%(y - exp(g)), Tmatrix=Tmatrix, Z=Z, H=A, Q=Q, nstates=nstates, d=d, P10=P10, outofsample=outofsample)$xtm1T)      # KS step
    
    if (j > 1 && g == gplus) break
    
    gplus <- g 
  }
    
  return(list(g=g, A=A, Ainv=Ainv, z=z, Psinv=Psinv))
}      


# Importance sampling (includes simulation smoothing)
IS <- function(y, g, n.draws, Psi, Ainv, A=A, z, max.iter, Tmatrix, Z, Q, nstates, d, P10, outofsample){
  
  n <- length(g)
  
  theta.tilde <- matrix(NA,n,n.draws)
  
  w <- rep(NA,n.draws)      # importance weights
  
  mi <- rep(NA,n.draws)      # importance weights that are more numerically stable
  
  Psi <- Psi
  
  P1 <- P10
  
  for (s in 1:n.draws){
    
    theta.plus1 <- rnorm(1,0,P1)      # initialize theta.plus1 from N(a1,P1)
    
    theta.plus <- as.vector(filter(c(theta.plus1, rnorm(n-1,0,sqrt(Q))), filter=c(Tmatrix), method="recursive"))      # draw theta.plus from g(theta)=Z theta_t-1 + eta, eta ~ N(0,Q)
    
    y.plus <- theta.plus + rnorm(n=n, mean=0, sd=diag(A^(1/2)))      # use theta.plus to generate y.plus from g(y|theta.plus)
    
    theta.plus.hat <- t(KFS(y=y.plus, Tmatrix=Tmatrix, Z=Z, H=A, Q=Q, nstates=nstates, d=d, P10=P10, outofsample=outofsample)$xtm1T)
    
    theta.tilde[,s] <- g + theta.plus - theta.plus.hat      # draw from g(theta|y)
    
    g_y.theta <- rep(NA, n)      # conditional density g(y|theta)
    
    for (i in 1:n){  
      g_y.theta[i] <- exp(-1/2*log(2*pi) + 1/2*log(Ainv[i,i]) - 1/2*(z[i]-theta.tilde[i,s])*Ainv[i,i]*(z[i]-theta.tilde[i,s]))      # this is the same as dnorm(z, mean=theta.tilde[,s], sd=diag(A^(1/2)), log=F)
    }
    
    w[s] <- prod(dpois(y, lambda = exp(theta.tilde[,s]), log=F))/prod(g_y.theta)
  
    p_i <- exp(sum(log(dpois(y, lambda = exp(theta.tilde[,s]), log=F))))
      
    g_i <- exp(sum(log(g_y.theta)))
      
    mi[s] <- log(p_i) - log(g_i)
    
  }
  
  return(list(theta.tilde=theta.tilde, w=w, mi=mi))
  
}      


Everything <- function(par, y, n.draws, Z, max.iter, nstates, d, P10, outofsample, opti, initial){
  
  Q <- exp(2*par[1])    
  
  Zstack <- diag(rep(Z,n))
  
  Psi <- Zstack%*%Omega%*%t(Zstack)
  
  Tmatrix <- par[2]
  
  
  # Newton-Raphson algorithm for mode estimation:
  
  NR.results <- NR(y=y, g=rep(mean(y),n), Psi=Psi, max.iter=max.iter, Tmatrix=Tmatrix, Z=Z, Q=Q, nstates=nstates, d=d, P10=P10, outofsample=outofsample)
  
  g <- NR.results$g      # mode estimate of theta
  
  A <- NR.results$A
  
  Ainv <- NR.results$Ainv
  
  z <- NR.results$z
  
  Psinv <- NR.results$Psinv
  
  #ts.plot(cbind(theta,g), col=c("black", "red"))
  

  #### Simulation smoothing and importance sampling ####
  
  if (initial == F){
    
    set.seed(403)      # use the same random numbers in order to evaluate and maximize the log-likelihood
    
    IS_results <- IS(y=y, g=g, n.draws=n.draws, Psi=Psi, Ainv=Ainv, A=A, z=z, max.iter=max.iter, Tmatrix=Tmatrix, Z=Z, Q=Q, nstates=nstates, d=d, P10=P1, outofsample=outofsample)
    
    theta.tilde <- IS_results$theta.tilde
    
    w <- IS_results$w
    
    mi <- IS_results$mi
    
    theta.tilde <- theta.tilde[,which(is.finite(w))]
    mi <- mi[which(is.finite(w))]
    w <- w[which(is.finite(w))]
    n.draws.bis <- length(w)
    
    
    exp.x.hat <- 0
    for (s in 1:n.draws.bis){
      exp.x.hat <- exp.x.hat + exp(theta.tilde[,s])*w[s]
    }
    exp.x.hat <- exp.x.hat/sum(w)
    
    x.hat <- 0
    for (s in 1:n.draws.bis){
      x.hat <- x.hat + theta.tilde[,s]*w[s]
    }
    x.hat <- x.hat/sum(w)
    
    mbar <- mean(mi)
    x.hat.m <- 0
    for (s in 1:n.draws.bis){
      x.hat.m <- x.hat.m + theta.tilde[,s]*exp(mi[s] - mbar)
    }
    x.hat.m <- x.hat.m/sum(exp(mi - mbar))
     
  } else {
    
    g_y.theta <- rep(NA, n)      # conditional density g(y|theta)
    
    for (i in 1:n){  
      g_y.theta[i] <- exp(-1/2*log(2*pi) + 1/2*log(Ainv[i,i]) - 1/2*(z[i]-g[i])*Ainv[i,i]*(z[i]-g[i]))      # this is the same as dnorm(z, mean=theta.tilde[,s], sd=diag(A^(1/2)), log=F)
    }
    
    w <- prod(dpois(y, lambda = exp(g), log=F))/prod(g_y.theta)
    
  }
  
  
  # Evaluation of the log-likelihood
  
  logl_t <- KFS(y=z, Tmatrix=Tmatrix, Z=Z, H=A, Q=Q, nstates=nstates, d=d, P10=P1, outofsample=outofsample)$logl_t
  
  if (initial == F){
    loglp <- (sum(logl_t) + mbar - log(n.draws.bis) + log(sum(exp(mi - mbar))))/n     # logl-likelihood to be maximized
  } else {
    loglp <- (sum(logl_t) + log(w))/n      # logl-likelihood to be maximized
  }
  
  
  if (opti) {
    return(-loglp)
  }
  else {
    return(list(x.hat=x.hat, x.hat.m=x.hat.m, exp.x.hat=exp.x.hat, loglp=loglp, w=w, theta.tilde=theta.tilde, mi=mi))
  }
}
  

                  
set.seed(42)

                  
#### Generate data ####

n <- 200

Tmatrix <- 0.5

Q <- 0.2

dvec <- Tmatrix*rep(0,n)

alpha <- as.vector(filter(rnorm(n,0,Q), filter=c(Tmatrix), method="recursive") + dvec)

Z <- 1

c <- rep(0,n)

mu <- c + Z*dvec

theta <- c + Z*alpha

y <- rpois(n,lambda=exp(theta))


# Stacked notation:

Tstack <- diag(1, n, n)

for (j in 1:(n-1)){
  for (i in (j+1):n){
    Tstack[i,j] <- (Tmatrix)^(i-1)
  }
}

P1 <- Q

Omega <- Tstack%*%diag(c(P1,rep(Q,n-1)))%*%t(Tstack)

Zstack <- diag(rep(Z,n))

Psi <- Zstack%*%Omega%*%t(Zstack)
                  


#### Maximize the log-likelihood ####

n.draws <- 5      # number of draws for the importance sampling
max.iter <- 100      # maximum number of iterations for the NR algorithm
nstates <- 1      # number of state variables
d <- 0      # number of nonstationary state variables of y


# Get an initial value for the parameters by maximizing the approximate log-likelihood

objopt.init <- ucminf(par=c(log(Q)/2,Tmatrix), Everything, y=y, n.draws=n.draws, Z=Z, max.iter=max.iter, 
                 nstates=nstates, d=d, P10=P1, outofsample=T, opti=T, initial=T, hessian=2, control=list(trace=T))

init <- objopt.init$par


# Use the initial value to start the maximization of the log-likelihood

start.time <- proc.time()
objopt <- ucminf(par=init, Everything, y=y, n.draws=n.draws, Z=Z, max.iter=max.iter, 
                 nstates=nstates, d=d, P10=P1, outofsample=T, opti=T, initial=F, hessian=2, control=list(trace=T))
el.time <- proc.time() - start.time

par <- objopt$par

obj <- Everything(par=objopt$par, y=y, n.draws=n.draws, Z=Z, max.iter=max.iter, 
                  nstates=nstates, d=d, P10=P1, outofsample=T, opti=F, initial=F)

ts.plot(cbind(theta,obj$x.hat.m), col=c("black", "red"))


