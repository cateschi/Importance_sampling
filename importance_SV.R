#### Packages ####

library(ucminf)
library(matrixcalc)      



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
  Ptm1T <- lapply(seq_len(len), function(X) matrix(0,nstates,nstates)) #smoothed variance of the state vector
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
NR <- function(y, g, Psi, max.iter, Tmatrix, Z, Q, c.coef, nstates, d, P10, outofsample){
  
  for (j in 1:max.iter){
    
    n <- length(y)
    
    p.dotdot <- rep(NA, n)
    for (i in 1:n){
      p.dotdot[i] <- -1/2*y[i]^2*exp(- c.coef - g[i])
    }
    
    A <- - diag(1/p.dotdot)
    
    p.dot <- rep(NA, n)
    for (i in 1:n){
      p.dot[i] <- -1/2 + 1/2*y[i]^2*exp(- c.coef - g[i])
    }
    
    z <- g + A%*%(p.dot)
    
    Ainv <- -diag(p.dotdot)
    
    Psinv <- svd(Psi)$v%*%diag(1/svd(Psi)$d)%*%t(svd(Psi)$u) 
    
    g <- t(KFS(y=g + A%*%p.dot, Tmatrix=Tmatrix, Z=Z, H=A, Q=Q, nstates=nstates, d=d, P10=P10, outofsample=outofsample)$xtm1T)      # KS step
    
    if (j > 1 && g == gplus) break
    
    gplus <- g 
  }
  
  return(list(g=g, A=A, Ainv=Ainv, z=z, Psinv=Psinv))
}      


# Modified efficient importance sampling (weighted least squares)
MEIS <- function(y, n.draws, max.iter, Tmatrix, Z, Q, c.coef, nstates, d, P10, outofsample){
  
  n <- length(y)
  
  theta.tilde <- matrix(NA,n,n.draws)
  
  w <- rep(NA,n.draws)      # importance weights
  
  mi <- rep(NA,n.draws)      # importance weights that are more numerically stable
  
  C <- diag(1, n)      # first guess for C
  
  A <- solve(C)
  
  b <- rep(0,n)      # first guess for b
  
  z <- solve(C)%*%b      # first guess for z
  
  for (j in 1:max.iter){
    
    set.seed(281192)
  
    for (s in 1:n.draws){
      
      g <- t(KFS(y=z, Tmatrix=Tmatrix, Z=Z, H=A, Q=Q, nstates=nstates, d=d, P10=P10, outofsample=outofsample)$xtm1T)
      
      theta.plus1 <- rnorm(1,0,sqrt(P1))      # initialize theta.plus1 from N(a1,P1)
      
      theta.plus <- as.vector(filter(c(theta.plus1, rnorm(n-1,0,sqrt(Q))), filter=c(Tmatrix), method="recursive"))      # draw theta.plus from g(theta)=Z theta_t-1 + eta, eta ~ N(0,Q)
      
      y.plus <- theta.plus + rnorm(n=n, mean=0, sd=diag(A^(1/2)))      # use theta.plus to generate y.plus from g(y|theta.plus)
      
      theta.plus.hat <- t(KFS(y=y.plus, Tmatrix=Tmatrix, Z=Z, H=A, Q=Q, nstates=nstates, d=d, P10=P10, outofsample=outofsample)$xtm1T)
      
      theta.tilde[,s] <- g + theta.plus - theta.plus.hat      # draw from g(theta|y)
      
      g_y.theta <- rep(NA, n)      # conditional density g(y|theta)
      
      for (i in 1:n){  
        g_y.theta[i] <- exp(-1/2*log(2*pi) + 1/2*log(C[i,i]) - 1/2*(z[i]-theta.tilde[i,s])*C[i,i]*(z[i]-theta.tilde[i,s]))      # this is the same as dnorm(z, mean=theta.tilde[,s], sd=diag(A^(1/2)), log=F)
      }
      
      w[s] <- prod(dnorm(y, 0, sqrt(exp(c.coef + theta.tilde[,s])), log=F))/prod(g_y.theta)
      
      p_i <- exp(sum(log(dnorm(y, 0, sqrt(exp(c.coef + theta.tilde[,s])), log=F))))
      
      g_i <- exp(sum(log(g_y.theta)))
      
      mi[s] <- log(p_i) - log(g_i)
    }
    
    logp <- dnorm(y, 0, sqrt(exp(c.coef + theta.tilde)), log=T)
    
    mbar <- mean(mi)
    
    weights <- exp(mi - mbar)
    
    k <- ncol(logp)
    
    beta <- matrix(NA, nrow=n, ncol=3)
    
    
    # Weighted least squares:
    
    for (i in 1:n){
      
      y.log <- logp[i,]
      
      X <- cbind(rep(1,k), theta.tilde[i,], -1/2*theta.tilde[i,]^2)
      
      invMat <- svd(t(X)%*%diag(c(weights))%*%X)$v%*%diag(1/svd(t(X)%*%diag(c(weights))%*%X)$d)%*%t(svd(t(X)%*%diag(c(weights))%*%X)$u) 
      
      beta[i,] <- invMat%*%t(X)%*%diag(c(weights))%*%y.log
      
    }
    
    if (j > 1 && norm((beta-beta.plus)/beta.plus, type = "I") < 0.0001) break
    
    b <- beta[,2]
    
    C <- diag(beta[,3])
    
    A <- solve(C)
    
    z <- solve(C)%*%b
    
    beta.plus <- beta
    
  }
  
  g <- t(KFS(y=z, Tmatrix=Tmatrix, Z=Z, H=A, Q=Q, nstates=nstates, d=d, P10=P10, outofsample=outofsample)$xtm1T)
  
  return(list(beta=beta, g=g, theta.tilde, w=w, mi=mi, b=b, C=C, A=A, z=z, mbar=mbar))
}      

                  
# Importance sampling (includes simulation smoothing)
IS <- function(y, g, n.draws, Psi, Ainv, A=A, z, max.iter, Tmatrix, Z, Q, c.coef, nstates, d, P10, outofsample){
  
  n <- length(g)
  
  theta.tilde <- matrix(NA,n,n.draws)
  
  w <- rep(NA,n.draws)      # importance weights
  
  mi <- rep(NA,n.draws)      # importance weights that are more numerically stable
  
  Psi <- Psi
  
  P1 <- P10
  
  for (s in 1:n.draws){
    
    theta.plus1 <- rnorm(1,0,sqrt(P1))      # initialize theta.plus1 from N(a1,P1)
    
    theta.plus <- as.vector(filter(c(theta.plus1, rnorm(n-1,0,sqrt(Q))), filter=c(Tmatrix), method="recursive"))      # draw theta.plus from g(theta)=Z theta_t-1 + eta, eta ~ N(0,Q)
    
    y.plus <- theta.plus + rnorm(n=n, mean=0, sd=diag(A^(1/2)))      # use theta.plus to generate y.plus from g(y|theta.plus)
    
    theta.plus.hat <- t(KFS(y=y.plus, Tmatrix=Tmatrix, Z=Z, H=A, Q=Q, nstates=nstates, d=d, P10=P10, outofsample=outofsample)$xtm1T)
    
    theta.tilde[,s] <- g + theta.plus - theta.plus.hat      # draw from g(theta|y)
    
    g_y.theta <- rep(NA, n)      # conditional density g(y|theta)
    
    for (i in 1:n){  
      g_y.theta[i] <- exp(-1/2*log(2*pi) + 1/2*log(Ainv[i,i]) - 1/2*(z[i]-theta.tilde[i,s])*Ainv[i,i]*(z[i]-theta.tilde[i,s]))      # this is the same as dnorm(z, mean=theta.tilde[,s], sd=diag(A^(1/2)), log=F)
    }
    
    w[s] <- prod(dnorm(y, 0, sqrt(exp(c.coef + theta.tilde[,s])), log=F))/prod(g_y.theta)
    
    p_i <- exp(sum(log(dnorm(y, 0, sqrt(exp(c.coef + theta.tilde[,s])), log=F))))
    
    g_i <- exp(sum(log(g_y.theta)))
    
    mi[s] <- log(p_i) - log(g_i)
  }
  
  return(list(theta.tilde=theta.tilde, w=w, mi=mi))
  
}      


Everything <- function(par, y, n.draws, Z, max.iter, nstates, d, P10, outofsample, opti, initial, efficient){
  
  Q <- exp(2*par[1])    
  
  Zstack <- diag(rep(Z,n))
  
  Psi <- Zstack%*%Omega%*%t(Zstack)
  
  Tmatrix <- par[2]
  
  c.coef <- 1
  
  
  if (efficient == F){
  
    # Newton-Raphson algorithm for mode estimation:
    
    NR.results <- NR(y=y, g=rep(mean(y),n), Psi=Psi, max.iter=max.iter, Tmatrix=Tmatrix, Z=Z, Q=Q, c.coef=c.coef, nstates=nstates, d=d, P10=P10, outofsample=outofsample)
    
    g <- NR.results$g      # mode estimate of theta
    
    A <- NR.results$A
    
    Ainv <- NR.results$Ainv
    
    z <- NR.results$z
    
    Psinv <- NR.results$Psinv
    
    
    #ts.plot(cbind(theta,g), col=c("black", "red"))
    
    
    #### Simulation smoothing and importance sampling ####
    
    if (initial == F){
      
      set.seed(403)      # use the same random numbers in order to evaluate and maximize the log-likelihood
      
      IS_results <- IS(y=y, g=g, n.draws=n.draws, Psi=Psi, Ainv=Ainv, A=A, z=z, max.iter=max.iter, Tmatrix=Tmatrix, Z=Z, Q=Q, c.coef=c.coef, nstates=nstates, d=d, P10=P1, outofsample=outofsample)
      
      theta.tilde <- IS_results$theta.tilde
      
      w <- IS_results$w
      
      mi <- IS_results$mi
      
      theta.tilde <- theta.tilde[,which(is.finite(w))]
      mi <- mi[which(is.finite(w))]
      w <- w[which(is.finite(w))]
      n.draws.bis <- length(w)
      
      
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
      
      w <- prod(dnorm(y, 0, sqrt(exp(c.coef + g)), log=F))/prod(g_y.theta)
      
    }
    
  } else if (efficient==T){
    
    MEIS.results <- MEIS(y=y, n.draws=n.draws, max.iter=max.iter, Tmatrix=Tmatrix, Z=Z, Q=Q, c.coef=c.coef, nstates=nstates, d=d, P10=P10, outofsample=outofsample)
    
    A <- MEIS.results$A
    
    z <- MEIS.results$z
    
    g <- MEIS.results$g
    
    w <- MEIS.results$w
    
    mi <- MEIS.results$mi
    
    theta.tilde <- MEIS.results$theta.tilde
    
    mbar <- MEIS.results$mbar
    
    x.hat.m <- g
    
    n.draws.bis <- n.draws
    
  }
  
  
  # Evaluation of the log-likelihood
  
  logl_t <- KFS(y=z, Tmatrix=Tmatrix, Z=Z, H=A, Q=Q, nstates=nstates, d=d, P10=P1, outofsample=outofsample)$logl_t
  
  if (initial == F){
    loglp <- (sum(logl_t) + mbar - log(n.draws.bis) + log(sum(exp(mi - mbar))))/n       # logl-likelihood to be maximized
  } else {
    loglp <- (sum(logl_t) + log(w))/n      # logl-likelihood to be maximized
  }
  
  
  if (opti) {
    return(-loglp)
  }
  else {
    return(list(x.hat=x.hat, x.hat.m=x.hat.m, loglp=loglp, w=w, theta.tilde=theta.tilde, mi=mi))
  }
}



#### Generate data ####           
                  
set.seed(42)

n <- 200

Tmatrix <- 0.98

Q <- 0.0225

dvec <- Tmatrix*rep(0,n)

alpha <- as.vector(filter(rnorm(n,0,sqrt(Q)), filter=c(Tmatrix), method="recursive") + dvec)

Z <- 1

c <- rep(0,n)

mu <- c + Z*dvec

theta <- Z*alpha

y <- rnorm(n,0,sqrt(exp(1+theta)))

y[which(0 <= y & y <= 0.4)] <- 0.4
y[which(-0.4 <= y & y <= 0)] <- -0.4


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


n.draws <- 30      # number of draws for the importance sampling
max.iter <- 100      # maximum number of iterations for the NR algorithm
nstates <- 1      # number of state variables
d <- 0      # number of nonstationary state variables of y
efficient <- F      # if TRUE, use modified efficient importance sampling instead of SPDK


if (efficient == F){
  
  # Get an initial value for the parameters by maximizing the approximate log-likelihood
  
  objopt.init <- ucminf(par=c(log(Q)/2,Tmatrix), Everything, y=y, n.draws=n.draws, Z=Z, max.iter=max.iter, 
                        nstates=nstates, d=d, P10=P1, outofsample=T, opti=T, initial=T, efficient=efficient, hessian=2, control=list(trace=T))
  
  init <- objopt.init$par
  
} else {
  
  init <- c(log(Q)/2,Tmatrix)
  
}



# Use the initial value to start the maximization of the log-likelihood

objopt <- ucminf(par=init, Everything, y=y, n.draws=n.draws, Z=Z, max.iter=max.iter, 
                 nstates=nstates, d=d, P10=P1, outofsample=T, opti=T, initial=F, efficient=efficient, hessian=2, control=list(trace=T))

par <- objopt$par

obj <- Everything(par=objopt$par, y=y, n.draws=n.draws, Z=Z, max.iter=max.iter, 
                  nstates=nstates, d=d, P10=P1, outofsample=T, opti=F, initial=F)

ts.plot(cbind(theta,obj$x.hat.m), col=c("black", "red"))


# Variance of x.hat

n.draws.bis <- ncol(obj$theta.tilde)

x.hat.sq <- 0
for (s in 1:n.draws.bis){
  x.hat.sq <- x.hat.sq + (obj$theta.tilde[,s])^2*obj$w[s]
}
x.hat.sq <- x.hat.sq/sum(obj$w)

var <- x.hat.sq - obj$x.hat^2


ts.plot(cbind(obj$x.hat, obj$x.hat-sqrt(var), obj$x.hat+sqrt(var)), col=c("red"), lty=c(1,2,2))




