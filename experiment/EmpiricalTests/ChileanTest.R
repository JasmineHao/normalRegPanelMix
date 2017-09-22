library(stargazer)
library(ggplot2)
library(reshape)
library(normalregMix)
library(foreign)
library(normalRegPanelMix)
df <- read.dta("C:/Users/Jasmine/Dropbox/workspace/R/package/normalRegPanelMix-0.2/experiment/DataClean/Chilean/ChileanClean.dta")

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
  # Describe the data
  ######################################################
  
  desc.each <- ind.each[ind.each$L != 0 ,c("si","lny","lnm","lnl","lnk")]
  desc.each <- desc.each[complete.cases(desc.each),]
  
  
  
  ######################################################
  # Select the data out
  ######################################################
  m.share <- cast(ind.each,id ~ year,value="si")#Collapse the dataframe into panel form , year against firm id
  row.names(m.share) <- m.share$id 
  m.share <- m.share[,!(colnames(m.share)=="id")] 
  T.cap <- dim(m.share)[2]
  
  estimate.df <- matrix(0,nr=5,nc=5)
  crit.df <- matrix(0,nr=5,nc=5)
  result.df <- matrix(0,nr=5,nc=5)
  # crit.df.boot <- matrix(0,nr=5,nc=5)
  
  ######################################################
  #For cross-sectional data
  ######################################################
  m.share.t <- m.share[,T.cap]
  m.share.t <- m.share.t[complete.cases(m.share.t)]
  N <- length(m.share.t)
  for (M in 1:5){
    
    out.h0 <- normalmixPMLE(y = m.share.t,m=M)
    an <- normalregMix::anFormula(out.h0$parlist,M,N)
    print(paste("T=",1,"M = ",M,"an=",an))
    out.h1 <- normalmixMaxPhi(y=m.share.t,parlist = out.h0$parlist,an=an)
    #Obtain critical value
    lr.crit <- normalmixCrit(y=m.share.t,parlist = out.h0$parlist)$crit
    crit.df[1,M] <- paste(round(lr.crit,2),collapse = ",")
    # lr.crit.boot <- normalmixCritBoot(y=m.share.t,parlist = out.h0$parlist,parallel = FALSE,nbtsp = 200)$crit
    # crif.df.boot[1,M] <- paste(round(lr.crit.boot,2),collapse = ",")
    estimate.df[1,M] <- 2 * max(out.h1$penloglik - out.h0$loglik)
    result.df[1,M] <- (2 * max(out.h1$penloglik - out.h0$loglik) > lr.crit[2])
    
  }

  ######################################################
  #For panel data
  ######################################################


  for (T in 2:5){
    t.start <- T.cap-T+1
    t.seq <- seq(from=t.start,to=t.start+T-1)
    m.share.t <- m.share[,t.seq]
    data <- list(Y = t(m.share.t[complete.cases(m.share.t),]), X = NULL,  Z = NULL)
    N <- dim(data$Y)[2]
    
    for (M in 1:5){
      out.h0 <- normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
      # phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
      #            beta = beta, N = N, T = T, M = M, p = p, q = q)
      an <- normalRegPanelMix::anFormula(out.h0$parlist,M,N,T) 
      print(paste("T=",T,"M = ",M,"an=",an))
      
      if (is.na(an)){ 
        an <- 1.0
      }
      out.h1 <- normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=an,update.alpha = 1)
      lr.estimate <- 2 * max(out.h1$penloglik - out.h0$loglik)
      estimate.df[T,M] <- lr.estimate
      
      lr.crit <- try(regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, cl=cl,parallel = TRUE,nrep=1000)$crit)
      
      if (class(lr.crit) == "try-error"){
        
          lr.crit <- c(0,0,0) 
        
      } 
      crit.df[T,M] <- paste(round(lr.crit,2),collapse = ",")
      
      # lr.crit <- regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, nbtsp = 199 ,parallel = FALSE)$crit
      # crit.df.boot[T,M] <- paste(round(lr.crit,2),collapse = ",")
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
  sink("C:/Users/Jasmine/Dropbox/workspace/R/package/normalRegPanelMix-0.2/experiment/DataClean/Chilean/crit.txt",append=TRUE)
  # stargazer(desc.each,type="text",title=paste("Descriptive data for Chilean Industry: ",ind.name))
  # stargazer(estimate.df,type='text',title = paste("Columbian Producer Data: Estimated LR for ",ind.name))
  # stargazer(result.df,title = ind.name)
  stargazer(crit.df,type="text",title=paste("Simulated crit for ",ind.name,each.code))
  sink()
}



# colnames(crit.df.boot) <- c("M=1","M=2","M=3","M=4","M=5")
# rownames(crit.df.boot) <- c("T=1","T=2","T=3","T=4","T=5")


# stargazer(crit.df,title=paste("estimate",ind.name))



