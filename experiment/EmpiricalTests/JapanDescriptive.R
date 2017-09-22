library(stargazer)
library(ggplot2)
library(reshape)
library(normalRegPanelMix)
library(normalregMix)
library(foreign)

ind.code <- unique(df$ciiu_3d)

ind_list <- c("food","textile", "wood","paper", "chemical", 
              "petro","plastic","ceramics","steel","othermetal",
              "metal product","machine","electronics",
              "transportation equipment","precision instrument",
              "other")

df <- read.csv(file="C:/Users/Jasmine/Dropbox/Dropbox/GNR/R/data_production_function_missing2zero.csv")
df[df==0] <- NA
df$t <- df$year
df <- df[order(df$id,df$t),]
df <- df[df$industry_2!=0,]
df <- df[!is.na(df$industry_2),]
ind.code <- unique(df$industry_2)

#ind.code = 2: machinery

desc.table = matrix(nc=5,nr=length(ind.code))
count = 0

for (each.code in ind.code){
  count = count + 1
  ind.each <- subset(df,industry_2==each.code)
  ind.name <- ind_list[each.code]
  
  desc.each <- ind.each
  desc.each <- desc.each[complete.cases(desc.each),]
  
  m.share <- cast(ind.each,id ~ year,value="lnmY_it")#Collapse the dataframe into panel form , year against firm id
  row.names(m.share) <- m.share$id 
  m.share <- m.share[,!(colnames(m.share)=="id")] 
  
  
  desc.table[count, ] <- c(ind.name,dim(desc.each)[1],dim(m.share)[1],round(mean(desc.each$lnmY_it),2),round(sd(desc.each$lnmY_it),2))
  
}

colnames(desc.table) <- c("Industry","NObs", "N","Mean","Sd")

sink("C:/Users/Jasmine/Dropbox/Dropbox/workspace/R/package/normalRegPanelMix-0.2/experiment/DataClean/Japan/desc.table.txt")
stargazer(desc.table,type="latex",title="Descriptive statistics for Japanses producer revenue share of intermediate material")
sink()
