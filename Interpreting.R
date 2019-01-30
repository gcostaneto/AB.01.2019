############################
# Interpreting results



require(plyr)
require(reshape2)


head(df.res)

# geral
model.pa <- ddply(df.res,.(pop,fold),summarise, r = cor(obs,pred))
model.pa
boxplot(r~pop,model.pa)

# por ambiente
pa.env <- ddply(df.res,.(pop,env,fold),summarise, r = cor(obs,pred))
pa.env
pa.env$pop.env <- paste(pa.env$pop,".",pa.env$env,sep="")
boxplot(r~pop.env,pa.env)
#------------------------------------------------------------------------------#
# para condicoes de manejo
df.res <- data.frame(df.res,colsplit(df.res$env,"_",c("year","site","N")))
head(df.res)

# por local
pa.site <- ddply(df.res,.(pop,site),summarise, r = cor(obs,pred))
pa.site

# por ano
pa.year <- ddply(df.res,.(pop,year),summarise, r = cor(obs,pred))
pa.year

# por ano
pa.N <- ddply(df.res,.(pop,N),summarise, r = cor(obs,pred))
pa.N

