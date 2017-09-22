library(stargazer)
library(ggplot2)
library(reshape)
library(normalRegPanelMix)
library(normalregMix)
library(foreign)
df <- read.dta("C:/Users/Jasmine/Dropbox/Dropbox/workspace/R/package/normalRegPanelMix-0.2/experiment/DataClean/Chilean/ChileanClean.dta")
ind.code <- unique(df$ciiu_3d)

# ind.code <- ind.code[1:3] #For test purpose, only use the first three industries


#ind.code = 2: machinery

desc.table = matrix(nc=5,nr=length(ind.code))
count = 0
for (each.code in ind.code){
  count = count + 1
  ind.each <- subset(df,ciiu_3d==each.code)
  ind.name <- ind.each$ciiu3d_descr[1]
  ind.each$lny <- log(ind.each$GO)
  ind.each$lnm <- log(ind.each$WI)
  ind.each$lnl <- log(ind.each$L)
  ind.each$lnk <- log(ind.each$K)
  
  desc.each <- ind.each[ind.each$L != 0 ,c("si","lny","lnm","lnl","lnk")]
  desc.each <- desc.each[complete.cases(desc.each),]
  
  m.share <- cast(ind.each,id ~ year,value="si")#Collapse the dataframe into panel form , year against firm id
  row.names(m.share) <- m.share$id 
  m.share <- m.share[,!(colnames(m.share)=="id")] 
  
  
  desc.table[count, ] <- c(ind.name,dim(desc.each)[1],dim(m.share)[1],round(mean(desc.each$si),2),round(sd(desc.each$si),2))
  
}
colnames(desc.table) <- c("Industry","NObs", "N","Mean","Sd")

sink("C:/Users/Jasmine/Dropbox/Dropbox/workspace/R/package/normalRegPanelMix-0.2/experiment/DataClean/Chilean/desc.table.txt")
stargazer(desc.table,type="latex",title="Descriptive statistics for Chilean producer revenue share of intermediate material")
sink()
