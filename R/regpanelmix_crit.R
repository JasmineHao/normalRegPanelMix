##Constraints for slsqp when p = 0
hin <- function(x){
  if(length(x)==3){
    c(x[1],x[1]*x[3] - x[2]^2,x[3])
  }else if(length(x)==6){
    c(x[1],x[1]*x[3] - x[2]^2,x[3],x[1]*x[6]-x[4]^2,x[3]*x[6]-x[5]^2,x[6]) 
  }
}

init.val.t <- function(z,q){
  if(q==0){
    t_0 <- c(max(0,z[1]),sqrt(max(0,z[1])*max(0,z[3])),max(0,z[3]))
  }else{
    s_mu <- max(0,z[1])
    s_sigma <- max(0,z[3])
    t_0 <- c(s_mu,sqrt(s_sigma*s_mu),s_sigma)
    s_beta_mu <- c()
    s_beta_sigma <- c()
    s_beta <- c()
    for (qq in 1:q){
      s_beta_qq <- max(z[(3+2*q+qq)],0)
      s_beta <- append(s_beta,s_beta_qq)
      s_beta_mu <- append(s_beta_mu,sqrt(s_beta_qq*s_mu))
      s_beta_sigma <- append(s_beta_sigma,sqrt(s_beta_qq*s_sigma))
    }
    t_0 <- append(t_0,s_beta_mu)
    t_0 <- append(t_0,s_beta_sigma)
    t_0 <- append(t_0,s_beta)
  }
  return(t_0)
}

hin_2 <- function(x){
  if(length(x)==3){
    c(x[1],x[3])}
  else if(length(x)==6){
    return(c(x[1],x[3],x[6]))
  }
}

heq <- function(x){
  if (length(x)==3){
    return(x[2]^2 - x[1]*x[3])}
  else if(length(x)==3){
    return(c(x[2]^2 - x[1]*x[3],x[4]^2-x[1]*x[6],x[5]^2-x[3]*x[6]))
  }
}

obj.func <- function(x,z,I){
  t(z - x) %*% I %*% (z - x) 
}

f.grad <- function(x,z,I){
  2*t(x - z) %*% I
}
#' @name regpanelmixCrit
#' @export
#' @title regpanelmixCrit
#' @param parallel Determines whether package \code{doParallel} is used for calculation
#' @param cl Cluster used for parallelization; if it is \code{NULL}, the system will automatically
#' @param nrep The number of replications used to compute p-values
#' @param ninits.crit The number of initial guesses to form critical values
#' @return A list with the following items:
#' \item{crit}{3 by 3 matrix of (0.1, 0.05, 0.01 critical values), jth row corresponding to k=j}
#' \item{pvals}{A vector of p-values at k = 1, 2, 3}
regpanelmixCrit <- function(y, x, parlist, z = NULL, values = NULL, parallel = TRUE,
                            cl = NULL, nrep = 1000, ninits.crit = 25, deep.check=0)
{
  # Computes the critical values of the modified EM test
  # and the estimated variance of the two-component MLE
  # Input
  #   y     : n by 1 vector of dependent variable
  #   x     : n by (q1-1) matrix of regressor not including an intercept
  #   parlist   : list including (alpha, mubeta, sigma)
  #   values (q1 by 1): values at wchich the p-values are computed
  # Output
  #   list(crit,pvals)
  #   crit = (10%, 5%, 1% critical values), pvals = p-values
  n  <- ncol(y)
  t  <- nrow(y)
  nt <- n*t
  y  <- as.vector(y)
  # n  <- length(y)
  p  <- 0

  alpha   <- parlist$alpha
  mubeta  <- parlist$mubeta
  sigma   <- parlist$sigma
  gam   <- parlist$gam
  m       <- length(alpha)

  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    y   <- y - z %*% gam
  }

  pvals <- NULL
  
  if (!is.null(x)){
    x     <- as.matrix(x)
    x1    <- cbind(1,x)
    q     <- ncol(x)
  }else{
    x1 <- rep(1,nt)
    q <- 0
  }
  # normalized data, W_p, f_p are n by m
  # W,f are nt by m.
  # f is pdf matrix.
  # W is the standardized error.
  # w is posterior matrix
  if (q !=  0){
    W  <- t(t(matrix(rep.int(y,m), ncol=m) - x1 %*% mubeta)/sigma)
  } else{
    W  <- t(t(matrix(rep.int(y,m), ncol=m) - t(matrix(rep(x1,m) * mubeta,nr=m)))/sigma)
  }
  W_T <- matrix(0,ncol = m, nrow = n) #Sum of error. Use panel version of likelihood
  f  <- t(t(exp(-W^2/2)/sqrt(2*pi))/sigma)  # pdf, nt by m
  f_T <- matrix(0,ncol=m,nrow=n)
  for (j in 1:m){
    W_T[,j] <- colSums(matrix(W[,j],nr=t))
    f_T[,j] <- exp(colSums(matrix(-(W[,j])^2/2,nr=t)))/(sqrt(2*pi)^t)/(sigma[j]^t)
  }

  f0_T <- colSums(t(f_T)*alpha)                 # data pdf, n by 1

  H <- hermite(W,sigma)

  w_m_T <- t(t(f_T)*alpha)/f0_T  # n by m matrix of w_{ji}
  if (m == 1) {
    S_alpha <- NULL
  } else {
    S_alpha <- (f_T[,1:(m-1)] - f_T[,m])/f0_T
  }
  S_mu <- matrix(0,nrow = n, ncol = m)
  S_beta  <- matrix(0, nrow=n, ncol=q*m)
  S_sigma <- matrix(0, nrow=n, ncol=m)

  for (j in 1:m){
    H_1 <- H[,j,1]
    H_2 <- H[,j,2]
    H_1_T <- colSums(matrix(H[,j,1],nr = t))
    H_2_T <- colSums(matrix(H[,j,2],nr = t))

    S_mu[,j] <- w_m_T[,j] * H_1_T
    S_sigma[,j] <-w_m_T[,j] * H_2_T
    for (qq in 1:q){
      H_1x <- colSums(matrix(x[,qq] * H[,j,1],nr=t))
      S_beta[, ((j-1)*q+qq)] <- H_1x * w_m_T[,j]
    }
  }
  # S_mu    <- w_m_T*H[,,1] #nt by m

  S_eta   <- cbind(S_alpha, S_mu, S_beta, S_sigma)


  if (!is.null(z)) {
    S_gam  <- matrix(0, nrow=n, ncol=p)

    for (j in 1:m){
      print("THIS MIGHT BE WRONG, IF YES, CHECK CASE q > 0 ")
      H_1 <- H[,j,1]
      # H_2 <- H[,j,2]
      for (pp in 1:p){
        H_1z <- rowSums(colSums(matrix(z[,pp] * H[,j,1],nr=t)))
        S_gam[,pp] =  S_gam[,pp] +  H_1z * w_m_T[,j]
      }
    }
    S_eta <- cbind(S_eta,S_gam)
  }

  n_lam <- q*(q+1)/2+ 2 * q + 3
  S_lam <- matrix(0, nrow=n, ncol=n_lam*m)
  n_eta <- q + 2
  for (j in 1:m) {
    H_1 <- matrix(H[,j,1],nr=t)
    H_2 <- matrix(H[,j,2],nr=t)
    H_3 <- matrix(H[,j,3],nr=t)
    H_4 <- matrix(H[,j,4],nr=t)

    # Need to see notes for the definition of s_lam
    s_lam_1 <- matrix(0,nr=n,ncol=1)
    s_lam_2 <- matrix(0,nr=n,ncol=1)
    s_lam_3 <- matrix(0,nr=n,ncol=1)
    
    if (q != 0){
    s_lam_4 <- matrix(0,nr=n,ncol=q) #s_beta mu
    s_lam_5 <- matrix(0,nr=n,ncol=q) #s_beta sigma
    s_lam_6 <- matrix(0,nr=n,ncol=q) #s_beta beta
    s_lam_7 <- matrix(0,nr=n,ncol=q*(q-1)/2) #s_betai betaj
    }
    
    for (nn in 1:n){
      H_1_n = H_1[,nn]
      H_2_n = H_2[,nn]
      H_3_n = H_3[,nn]
      H_4_n = H_4[,nn]
      tmp <- ((H_1_n) %*% t(H_1_n))
      s_lam_1[nn,] = w_m_T[nn,j] * (sum(H_2_n) + sum(tmp[row(tmp)!=col(tmp)])/2) 
      tmp <- ((H_1_n) %*% t(H_2_n))
      s_lam_2[nn,] = w_m_T[nn,j] * (3*sum(H_3_n)  + sum(tmp[row(tmp)!=col(tmp)]))
      tmp <- ((H_2_n) %*% t(H_2_n))
      s_lam_3[nn,] = w_m_T[nn,j] * (3*sum(H_4_n) + sum(tmp[row(tmp)!=col(tmp)])/2)
      if (q != 0){
        for (qq in 1:q){
          x_nn =  x[((nn-1)*t + 1):(nn*t),qq]
          tmp <- ((x_nn * H_1_n) %*% t( H_1_n))
          s_lam_4[nn,qq] = w_m_T[nn,j] * (2*H_2_n %*% x_nn + sum(tmp[row(tmp)!=col(tmp)])) #S_beta mu
          tmp <- ((x_nn * H_1_n) %*% t( H_2_n))
          s_lam_5[nn,qq] = w_m_T[nn,j] * (3*H_3_n %*% x_nn + 2 * sum(tmp[row(tmp)!=col(tmp)]))
          tmp <- ((x_nn * H_1_n) %*% t((x_nn * H_1_n)))
          s_lam_6[nn,qq] = w_m_T[nn,j] * ((H_2_n %*% (x_nn^2)) + sum(tmp[row(tmp)!=col(tmp)]) /2 )
          
          if (q > 1){
            qq <- 1
            for (jj in 1:(q-1)) {
              for (ii in (jj+1):q) {
                xx_i <- x[((nn-1)*t + 1):(nn*t),jj]
                xx_j <- x[((nn-1)*t + 1):(nn*t),ii]
                tmp <- (H_1_n * xx_i) %*% t((H_1_n * xx_j))
                s_lam_7[nn,qq] = w_m_T[nn,jj] * (2 * H_2_n %*% (xx_i * xx_j) + sum(tmp[row(tmp)!=col(tmp)]) )
                qq <- qq+1
              }
            }
          }
        }
      }
      
    }

    if (q == 0){
      S_lam[, ((j-1)*n_lam+1):(j*n_lam)] <- cbind(s_lam_1,s_lam_2,s_lam_3)
    }else{
      S_lam[, ((j-1)*n_lam+1):(j*n_lam)] <- cbind(s_lam_1,s_lam_2,s_lam_3,s_lam_4,s_lam_5,s_lam_6,s_lam_7)
    }
  }

  I_eta     <- t(S_eta) %*% S_eta/n
  I_lam     <- t(S_lam) %*% S_lam/n
  I_el      <- t(S_eta) %*% S_lam/n
  if (qr(I_eta)$rank == nrow(I_eta)) {
    I_eta_lam <- I_lam - t(I_el) %*% solve(I_eta,I_el)
  } else {
    stop("The critical value cannot be computed due to singularity of some matrices.
         Please try a bootstrap version, regpanelmixCritBoot and regpanelmixMEMtestBoot.")
  }


  # generate u ~ N(0,I_eta_lam)
  set.seed(123456)
  e <- eigen(I_eta_lam, symmetric=TRUE)  # eigenvalue decomposition is slower than chol but more stable
  u <- t(e$vec %*% (t(e$vec) * sqrt(e$val)) %*% matrix(rnorm(nrep*n_lam*m), nrow=n_lam*m))


  LR <- matrix(0, nrow=nrep, ncol=m)
  
  if ( (parallel) && (is.null(cl)) ) {
    # ncpus <- parallel::detectCores()
    # cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
    cl <- makeCluster(detectCores()-1)
    
  }
  
  count = 0 #Count the number of inside convex cone optimization
  for (j in 1:m) {
    I_j <- I_eta_lam[((j-1)*n_lam+1):(j*n_lam),((j-1)*n_lam+1):(j*n_lam)]
    if (qr(I_j)$rank == nrow(I_j)) {
      Z_j <- u[,((j-1)*n_lam+1):(j*n_lam)] %*% solve(I_j)    # nrep by n_lam matrix
      
      
      if(parallel){
        registerDoParallel(cl)
        t_j <- foreach(i = 1:nrep,.combine = rbind)%dopar%{
          
          t  <- slsqp(init.val.t(Z_j[i,],q),obj.func,gr = f.grad,z=Z_j[i,],I=I_j,hin = hin)$par
          cond <- hin(t)
          # print(t)
          if (!(abs(cond[3]) < 1e-5)){
            count <- count + 1
            if (deep.check == 1){
              cseq <- lapply(t , function(x){seq(0,2*x,length=5)})
              draws <- as.matrix(expand.grid(cseq))
    
              eq.opt <-  function(x){slsqp(x,obj.func,gr = f.grad,z=Z_j[i,],I=I_j,hin = hin_2,heq = heq)$par}
    
              eq.sol <- t(apply(draws,1,eq.opt))
              draw.LR <- rowSums(eq.sol %*% I_j *eq.sol)
              t <- eq.sol[which.max(draw.LR),]
            }else{
              t <- slsqp(Z_j[i,],obj.func,gr=f.grad,z=Z_j[i,],I=I_j,hin = hin_2,heq=heq)$par
            }
          }
          t
        }
        
      }
      
      else{
        t_j <- matrix(0,nr=nrow(Z_j),nc=ncol(Z_j))
        for (i in 1:nrep){
          t_0 <- init.val.t(Z_j[i,],q)
          t  <- slsqp(t_0,obj.func,gr = f.grad,z=Z_j[i,],I=I_j,hin = hin)$par

          cond <- hin(t)

          if (!(abs(cond[3]) < 1e-5)){
            count <- count + 1
            if (deep.check == 1){
              cseq <- lapply(t , function(x){seq(0,2*x,length=5)})
              draws <- as.matrix(expand.grid(cseq))

              eq.opt <-  function(x){slsqp(x,obj.func,gr = f.grad,z=Z_j[i,],I=I_j,hin = hin_2,heq = heq)$par}

              eq.sol <- t(apply(draws,1,eq.opt))
              draw.LR <- rowSums(eq.sol %*% I_j *eq.sol)
              t <- eq.sol[which.max(draw.LR),]
            }else{
              t <- slsqp(Z_j[i,],obj.func,gr=f.grad,z=Z_j[i,],I=I_j,hin = hin_2,heq=heq)$par
            }
          }

          t_j[i,] <- t
        }
        
        
      }
    } 
    else {
      stop("The critical value cannot be computed due to singularity of some matrices.
           Please try a bootstrap version, regpanelmixCritBoot and regpanelmixMEMtestBoot.")
    }
   
    
    LR[,j]  <- rowSums((t_j %*% I_j) * t_j) 
    }

  # if (parallel) { parallel::stopCluster(cl) }

  max_EM <- apply(LR, 1, max)
  max_EM_sort <- sort(max_EM)
  qc <- floor(nrep*c(0.90,0.95,0.99))
  crit <- max_EM_sort[qc]
  
  if (!is.null(values)) {
    q1 <- length(values)
    pvals <- rowMeans(t(matrix(rep.int(max_EM_sort,q1),ncol=q1)) >= values)
  }
  print(count)
  if( (parallel) && (is.null(cl)) ){
  stopCluster(cl)}
  return(list(crit = crit, pvals = pvals,test= max_EM_sort ))

  } # end function regpanelmixCrit


#' @description Computes the bootstrap critical values of the modified EM test.
#' @export
#' @title regpanelmixCritBoot
#' @name regpanelmixCritBoot
#' @param y n by 1 vector of data for y
#' @param x n by q vector of data for x
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m))
#' @param z n by p matrix of regressor associated with gamma
#' @param values 3 by 1 Vector of length 3 (k = 1, 2, 3) at which the p-values are computed
#' @param ninits The number of initial candidates to be generated
#' @param nbtsp The number of bootstrap observations; by default, it is set to be 199
#' @param parallel Determines whether package \code{doParallel} is used for calculation
#' @param cl Cluster used for parallelization (optional)
#' @return A list with the following items:
#' \item{crit}{3 by 3 matrix of (0.1, 0.05, 0.01 critical values), jth row corresponding to k=j}
#' \item{pvals}{A vector of p-values at k = 1, 2, 3}
regpanelmixCritBoot <- function (y, x, parlist, z = NULL, values = NULL, ninits = 100,
                                 nbtsp = 199, parallel = FALSE, an = 0.5, cl = NULL) {
  # if (normalregpanelmix.test.on) # initial values controlled by normalregpanelmix.test.on
  #   set.seed(normalregpanelmix.test.seed)



  n  <- ncol(y)
  t  <- nrow(y)
  nt <- n*t
  y  <- as.vector(y)
  
  
  alpha   <- parlist$alpha
  mubeta  <- parlist$mubeta
  sigma   <- parlist$sigma
  gam   <- parlist$gam
  m       <- length(alpha)
  
  
  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    y   <- y - z %*% gam
  }else{
    p <- 0
  }
  
  
  if (!is.null(x)){
    x     <- as.matrix(x)
    q     <- ncol(x)
    mu <- mubeta[1,]
    beta <- mubeta[2:(q+1),]  #CHECKED
  }else{
    q <- 0
    mu <- mubeta
    beta  <- NULL
  }
  
  an    <- anFormula(parlist = parlist, m = m, n = n, t = t, q = q)
  

  pvals <- NULL
  

  # Generate bootstrap observations
  
  ybset <- replicate(nbtsp, generateData(N = n, T = t,M=m ,alpha = alpha, mu = mu, beta = beta, sigma = sigma,p= p,q=q))
  # tmp <- lapply(seq_len(ncol(tmp)),function(i) tmp[,i])
   
  if (!is.null(z)) {
    zgam <- as.matrix(z) %*% gam
    ybset <- ybset + replicate(nbtsp, as.vector(zgam))
  }

  if (parallel) {
    if (is.null(cl)){
      cl <- makeCluster(detectCores())}
    
      registerDoParallel(cl)
      print("Parallel Bootstrap Crit")
      print(paste("an=",an))
      out <- foreach (j.btsp = 1:nbtsp) %dopar% {
      regpanelmixMEMtest (y = ybset[,j.btsp]$Y, x = ybset[,j.btsp]$X , m = m, t = t, an = an,
                          z = z, ninits = ninits, crit.method = "none") }
    on.exit(cl)
  }
  else
    {
      out <- lapply(seq_len(ncol(ybset)), 
                    function(i) regpanelmixMEMtest(y = ybset[,i]$Y,x=ybset[,i]$X,m=M,t=T,z=NULL,ninits=10,an=an,crit.method = "none"))
    # out <- apply(ybset, 3, regpanelmixMEMtest, x = x, m = m, t = t, z = z,
                 # ninits = ninits, crit.method = "none")
    }
  emstat.b <- sapply(out, "[[", "emstat")  # 3 by nbstp matrix
  emstat.b <- sort(emstat.b)
  # emstat.b <- t(apply(emstat.b, 1, sort))
  #####################################
  # FIXIT FROM HERE, NOT SURE THE MODEL IS RIGHT
  #####################################
  q <- ceiling(nbtsp*c(0.90,0.95,0.99))
  # crit <- emstat.b[, q]
  crit <- emstat.b[q]
  if (!is.null(values)) { pvals <- rowMeans(emstat.b > values) }

  return(list(crit = crit, pvals = pvals))
}  # end function regpanelmixCritBoot
