# Computes objective function used in computing LR_2
obj_zIz <- function(b,Z,I) {
  q       <- length(b)-2
  lam_mu  <- b[1]
  lam_sig <- b[2]
  lam_bet <- b[3:length(b)]
  t       <- double(2+2*q+q*(q+1)/2)
  t[1]          <- 6*lam_mu*lam_sig
  t[2:(q+1)]    <- 2*lam_bet*lam_mu
  t[q+2]        <- 12*lam_sig^2
  t[(q+3):(2+2*q)]  <- 2*lam_bet*lam_sig
  
  t[(2+2*q+1):(2+2*q+q)] <- lam_bet^2
  
  if (q>=2) {
    ii <- combn(1:q,2)
    t[(3+3*q):(2+2*q+q*(q+1)/2)] <- lam_bet[ii[1,]]*lam_bet[ii[2,]]
  }
  
  R <- chol(I)
  res <- R %*% (t-as.vector(Z))
  
  res
  
} # end function obj_zIz

# Computes Jacobian of the objective function used in computing LR_2
obj_zIz.jac <- function(b,Z,I) {
  
  q     <- length(b)-2
  lam_mu  <- b[1]
  lam_sig <- b[2]
  lam_bet <- b[3:length(b)]
  t     <- double(2+2*q+q*(q+1)/2)
  t[1]        <- 6*lam_mu*lam_sig
  t[2:(q+1)]  <- 2*lam_bet*lam_mu
  t[q+2]      <- 12*lam_sig^2
  t[(q+3):(2+2*q)]  <- 2*lam_bet*lam_sig
  
  t[(2+2*q+1):(2+2*q+q)]  <- lam_bet^2
  
  if (q >= 2) {
    ii <- combn(1:q,2)
    t[(3+3*q):(2+2*q+q*(q+1)/2)] <- lam_bet[ii[1,]]*lam_bet[ii[2,]]
  }
  
  R <- chol(I)
  
  R_sub = NULL
  if (q >= 2) {
    R_sub = cbind(lam_bet[2:q], lam_bet[1]*diag(q-1))
    if (q >= 3) {
      for (i in 2:(q-1)) {
        R_sub_2 = cbind(matrix(0,nrow=q-i,ncol=i-1),lam_bet[(i+1):q],lam_bet[i]*diag(q-i))
        R_sub = rbind(R_sub, R_sub_2)
      }
      
    }
    
  }
  
  G <- cbind(c(6*lam_sig, 2*lam_bet, double(1+q+q*(q+1)/2)),
             c(6*lam_mu, double(q), 24*lam_sig, 2*lam_bet, double(q*(q+1)/2)),
             rbind(double(q), 2*lam_mu*diag(q), double(q), 2*lam_sig*diag(q),
                   2*diag(lam_bet,nrow=length(lam_bet)), R_sub))
  
  jac <- R%*%G
  
  return(jac)
  
} # end function obj_zIz.jac

# Computes LR_2 for given Z and I, where q is dim(x)
LR_2.comp <- function(Z, I, q, ninits = 25) {
  if (normalregMix.test.on) # initial values controlled by normalregMix.test.on
    set.seed(normalregMix.test.seed)
  
  con.nls = minpack.lm::nls.lm.control(maxit = 400)
  
  n_t <- q+2
  ll <- vector("list",n_t)
  for (i in 1:(q+2)) {ll[[i]] <- c(1,-1)}
  # g <- t(as.matrix(expand.grid(ll)))
  g <- t((expand.grid(ll)))
  
  lam_sig_0 <- sqrt(abs(Z[2+q]/12))
  lam_bet_0 <- sqrt(abs(Z[(2+2*q+1):(2+2*q+q)]))
  lam_mu_0 <- sqrt(abs(Z[1]*Z[2]/12/lam_sig_0/lam_bet_0[1]))
  t0 <- cbind(c(lam_mu_0,lam_sig_0,lam_bet_0), 8*matrix(rnorm(n_t*(ninits-1)), nrow=n_t))
  
  LR_2_all <- double(ninits)
  LR_2_all.optim <- double(ninits)
  
  for (irep in 1:ninits) {
    tg <- g*t0[,irep]  # q+2 by ncol(g) matrix
    res <- apply(tg,2,obj_zIz,Z,I)
    obj_val <- colSums(res^2)
    lm.par <- tg[,which.min(colSums(res^2))]  # uses the column of tg that gives the smallest obj value
    nls.out <- minpack.lm::nls.lm(par=lm.par, fn = obj_zIz, jac = obj_zIz.jac, control = con.nls, Z=Z, I=I)
    LR_2_all[irep] <- sum((nls.out$fvec)^2)
  }
  
  LR_2 <- Z %*% I %*% Z - min(LR_2_all)
  
  LR_2
  
} # end function LR_2.comp
