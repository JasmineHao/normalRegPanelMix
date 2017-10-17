library(stargazer)
library(ggplot2)
library(reshape)
library(normalregMix)
library(foreign)
library(normalRegPanelMix)
df <- read.dta("C:/Users/haoja/Dropbox/Dropbox/workspace/R/package/DataClean/Chilean/ChileanClean.dta")

# ind.code <- ind.code[1:3] #For test purpose, only use the first three industries
# ind.code <- c(352,342,369,381,321,313,341,322,390,311,351,324,356,312)

ind.code <- c(331,382,332,384,352,342,369,381,321,313,341,322,390,311,351,324,356,312)
ind.count <- length(ind.code)
cl <- makeCluster(7)
count = 0
for (each.code in ind.code){
  t <- Sys.time()
  ind.each <- subset(df,ciiu_3d==each.code)
  ind.name <- ind.each$ciiu3d_descr[1]
  ind.each$lny <- log(ind.each$GO)
  ind.each$lnm <- log(ind.each$WI)
  ind.each$lnl <- log(ind.each$L)
  ind.each$lnk <- log(ind.each$K)
  ######################################################
  # Refine the data
  ######################################################
  ind.each <- ind.each[ind.each$L != 0 ,c("id","year","si","lny","lnm","lnl","lnk")]
  year.list <- sort(unique(ind.each$year))
  ######################################################
  # Select the data out
  ######################################################
  T.cap <- max(year.list)
  
  estimate.m.df <- matrix(0,nr=5,nc=5)
  crit.m.df <- matrix(0,nr=5,nc=5)
  result.m.df <- matrix(0,nr=5,nc=5)
  
  estimate.l.df <- matrix(0,nr=5,nc=5)
  crit.l.df <- matrix(0,nr=5,nc=5)
  result.l.df <- matrix(0,nr=5,nc=5)
  
  estimate.k.df <- matrix(0,nr=5,nc=5)
  crit.k.df <- matrix(0,nr=5,nc=5)
  result.k.df <- matrix(0,nr=5,nc=5)
  
  estimate.all.df <- matrix(0,nr=5,nc=5)
  crit.all.df <- matrix(0,nr=5,nc=5)
  result.all.df <- matrix(0,nr=5,nc=5)
  # crit.df.boot <- matrix(0,nr=5,nc=5)
  
  ######################################################
  #For cross-sectional data
  ######################################################
  ind.each.t <- ind.each[ind.each$year==T.cap,]
  ind.each.t <- ind.each.t[complete.cases(ind.each.t),]
  N <- length(ind.each.t)
  for (M in 1:5){
    ############################################
    #Estimate with regressor lnM
    ############################################
    out.h0 <- normalmixPMLE(y = ind.each.t$si, x = ind.each.t$lnm ,m=M)
    an <- normalregMix::anFormula(out.h0$parlist,M,N)
    print(paste("T=",1,"M = ",M,"an=",an))
    out.h1 <- regmixMaxPhi(y=ind.each.t$si,x = as.matrix(ind.each.t$lnm) ,parlist = out.h0$parlist,an=an)
      #Obtain critical value
    lr.crit <- regmixCrit(y=ind.each.t$si,x = as.matrix(ind.each.t$lnm),parlist = out.h0$parlist)$crit
    crit.m.df[1,M] <- paste(round(lr.crit,2),collapse = ",")
    # lr.crit.boot <- normalmixCritBoot(y=m.share.t,parlist = out.h0$parlist,parallel = FALSE,nbtsp = 200)$crit
    # crif.df.boot[1,M] <- paste(round(lr.crit.boot,2),collapse = ",")
    estimate.m.df[1,M] <- 2 * max(out.h1$penloglik - out.h0$loglik)
    result.m.df[1,M] <- (2 * max(out.h1$penloglik - out.h0$loglik) > lr.crit[2])
    
    ############################################
    #Estimate with regressor lnK
    ############################################
    out.h0 <- normalmixPMLE(y = ind.each.t$si, x = ind.each.t$lnk ,m=M)
    an <- normalregMix::anFormula(out.h0$parlist,M,N)
    print(paste("T=",1,"M = ",M,"an=",an))
    out.h1 <- regmixMaxPhi(y=ind.each.t$si,x = as.matrix( ind.each.t$lnk) ,parlist = out.h0$parlist,an=an)
    #Obtain critical value
    lr.crit <- regmixCrit(y=ind.each.t$si,x = as.matrix( ind.each.t$lnk),parlist = out.h0$parlist)$crit
    crit.k.df[1,M] <- paste(round(lr.crit,2),collapse = ",")
    # lr.crit.boot <- normalmixCritBoot(y=m.share.t,parlist = out.h0$parlist,parallel = FALSE,nbtsp = 200)$crit
    # crif.df.boot[1,M] <- paste(round(lr.crit.boot,2),collapse = ",")
    estimate.k.df[1,M] <- 2 * max(out.h1$penloglik - out.h0$loglik)
    result.k.df[1,M] <- (2 * max(out.h1$penloglik - out.h0$loglik) > lr.crit[2])
    ############################################
    #Estimate with regressor lnL
    ############################################
    out.h0 <- normalmixPMLE(y = ind.each.t$si, x = ind.each.t$lnl ,m=M)
    an <- normalregMix::anFormula(out.h0$parlist,M,N)
    print(paste("T=",1,"M = ",M,"an=",an))
    out.h1 <- regmixMaxPhi(y=ind.each.t$si,x = as.matrix( ind.each.t$lnl) ,parlist = out.h0$parlist,an=an)
    #Obtain critical value
    lr.crit <- regmixCrit(y=ind.each.t$si,x = as.matrix( ind.each.t$lnl),parlist = out.h0$parlist)$crit
    crit.l.df[1,M] <- paste(round(lr.crit,2),collapse = ",")
    # lr.crit.boot <- normalmixCritBoot(y=m.share.t,parlist = out.h0$parlist,parallel = FALSE,nbtsp = 200)$crit
    # crif.df.boot[1,M] <- paste(round(lr.crit.boot,2),collapse = ",")
    estimate.l.df[1,M] <- 2 * max(out.h1$penloglik - out.h0$loglik)
    result.l.df[1,M] <- (2 * max(out.h1$penloglik - out.h0$loglik) > lr.crit[2])
    
  }
  
  ######################################################
  #For panel data
  ######################################################
  
  
  for (T in 2:5){
    t.start <- T.cap-T+1
    #Reshape the data so that I can apply the test
    ind.each.t <- ind.each[ind.each$year>=t.start,]
    ind.each.t <- ind.each.t[complete.cases(ind.each.t),]
    ind.each.y <- cast(ind.each.t[,c("id","year","si")],id ~ year,value="si")
    id.list    <- ind.each.y[complete.cases(ind.each.y),"id"]
    #Remove the incomplete data, need balanced panel
    ind.each.t <- ind.each.t[ind.each.t$id %in% id.list,]
    ind.each.t <- ind.each.t[order(ind.each.t$id,ind.each.t$year),]
    #Reshape the Y 
    ind.each.y <- cast(ind.each.t[,c("id","year","si")],id ~ year,value="si")
    ind.each.y <- ind.each.y[,colnames(ind.each.y)!="id"]
    N <- dim(ind.each.y)[1]
    
    for (M in 1:5){
      ############################################
      #Estimate with regressor lnM
      ############################################
      data <- list(Y = t(ind.each.y), X = ind.each.t$lnm,  Z = NULL)
      out.h0 <- regpanelmixPMLE(y=data$Y,x=data$X, z = NULL,m=M,vcov.method = "none")
      # phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
      #            beta = beta, N = N, T = T, M = M, p = p, q = q)
      an <- normalRegPanelMix::anFormula(out.h0$parlist,M,N,T,q=1) 
      print(paste("T=",T,"M = ",M,"an=",an))
      
      if (is.na(an)){ 
        an <- 1.0
      }
      out.h1 <- regpanelmixMaxPhi(y=data$Y,x=data$X, z = NULL,parlist=out.h0$parlist,an=an,update.alpha = 1)
      lr.estimate <- 2 * max(out.h1$penloglik - out.h0$loglik)
      estimate.m.df[T,M] <- lr.estimate
      lr.crit <- try(regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, nbtsp = 199 ,parallel = FALSE)$crit)
      if (class(lr.crit) == "try-error"){
        lr.crit <- c(0,0,0) 
      } 
      crit.m.df[T,M] <- paste(round(lr.crit,2),collapse = ",")
      ############################################
      #Estimate with regressor lnK
      ############################################
      data <- list(Y = t(ind.each.y), X = ind.each.t$lnk,  Z = NULL)
      out.h0 <- regpanelmixPMLE(y=data$Y,x=data$X, z = NULL,m=M,vcov.method = "none")
      # phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
      #            beta = beta, N = N, T = T, M = M, p = p, q = q)
      an <- normalRegPanelMix::anFormula(out.h0$parlist,M,N,T,q=1) 
      print(paste("T=",T,"M = ",M,"an=",an))
      if (is.na(an)){ 
        an <- 1.0
      }
      out.h1 <- regpanelmixMaxPhi(y=data$Y,x=data$X, z = NULL,parlist=out.h0$parlist,an=an,update.alpha = 1)
      lr.estimate <- 2 * max(out.h1$penloglik - out.h0$loglik)
      estimate.k.df[T,M] <- lr.estimate
      lr.crit <- try(regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, nbtsp = 199 ,parallel = FALSE)$crit)
      if (class(lr.crit) == "try-error"){
        lr.crit <- c(0,0,0) 
      } 
      crit.k.df[T,M] <- paste(round(lr.crit,2),collapse = ",")
      ############################################
      #Estimate with regressor lnL
      ############################################
      data <- list(Y = t(ind.each.y), X = ind.each.t$lnm,  Z = NULL)
      out.h0 <- regpanelmixPMLE(y=data$Y,x=data$X, z = NULL,m=M,vcov.method = "none")
      # phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
      #            beta = beta, N = N, T = T, M = M, p = p, q = q)
      an <- normalRegPanelMix::anFormula(out.h0$parlist,M,N,T,q=1) 
      print(paste("T=",T,"M = ",M,"an=",an))
      if (is.na(an)){ 
        an <- 1.0
      }
      out.h1 <- regpanelmixMaxPhi(y=data$Y,x=data$X, z = NULL,parlist=out.h0$parlist,an=an,update.alpha = 1)
      lr.estimate <- 2 * max(out.h1$penloglik - out.h0$loglik)
      estimate.l.df[T,M] <- lr.estimate
      lr.crit <- try(regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, nbtsp = 199 ,parallel = FALSE)$crit)
      if (class(lr.crit) == "try-error"){
        lr.crit <- c(0,0,0) 
      } 
      crit.l.df[T,M] <- paste(round(lr.crit,2),collapse = ",")
    }
  }
  ###################################################################
  #     Output
  ###################################################################
  count = count + 1
  print(Sys.time()-t)
  print(paste(count,"/",ind.count))
  colnames(estimate.df) <- c("M=1","M=2","M=3","M=4","M=5")
  rownames(estimate.df) <- c("T=1","T=2","T=3","T=4","T=5")
  
  colnames(crit.df) <- c("M=1","M=2","M=3","M=4","M=5")
  rownames(crit.df) <- c("T=1","T=2","T=3","T=4","T=5")
  sink("C:/Users/Jasmine/Dropbox/workspace/R/package/normalRegPanelMix/experiment/DataClean/Chilean/regressorCrit.txt",append=TRUE)
  # stargazer(desc.each,type="text",title=paste("Descriptive data for Chilean Industry: ",ind.name))
  # stargazer(estimate.df,type='text',title = paste("Columbian Producer Data: Estimated LR for ",ind.name))
  # stargazer(result.df,title = ind.name)
  stargazer(crit.df,type="text",title=paste("Simulated crit for ",ind.name,each.code))
  sink()
}



# colnames(crit.df.boot) <- c("M=1","M=2","M=3","M=4","M=5")
# rownames(crit.df.boot) <- c("T=1","T=2","T=3","T=4","T=5")


# stargazer(crit.df,title=paste("estimate",ind.name))



