###########################################################
# Anshu Chen
# Econ 419
# Required packages: lubridate, rugarch, reshape, RCurl, httr, gdata, zoo
###########################################################

# Creates Goyal's monthly volatility proxy

daily  <- read.csv("C:/Users/Anshu/Box Sync/junior/econ 419/paper/crsp daily.csv", 
                    na.strings="NaN", stringsAsFactors=FALSE)
daily$Date <- as.Date(as.character(daily$Date), "%Y%m%d")

daily$sq <- daily$R^2
daily$aut <- NA
for (i in 2:nrow(daily)) {
  daily$aut[i] <- 2*daily$R[i]*daily$R[i-1]
}

# Delete first "aut" of every month since we don't want last month's observations mixing with current

library(lubridate)

for (i in 2:nrow(daily)) {
  if (month(daily$Date[i])!=month(daily$Date[i-1])){
    daily$aut[i] <- NA
  }
}

bymonth <- aggregate(cbind(sq, aut)~month(Date)+year(Date), data=daily,FUN=sum)
bymonth$mvol <- bymonth$sq + bymonth$aut
bymonth$month <- bymonth$`month(Date)`
bymonth$year <- bymonth$`year(Date)`
bymonth <- bymonth[,-c(1,2)]

# Now I subset dates. Goyal used July 1962 to December 1998. 

gvol <- bymonth[bymonth$year<1999,]

# New data: July 1962 to October 2016.

nvol <- bymonth[bymonth$year>1961,]

###########################################################

# Commodity prices

comm <- read.csv("C:/Users/Anshu/Box Sync/junior/econ 419/paper/commodities.csv", 
                 na.strings="NaN", stringsAsFactors=FALSE)
comm <- comm[,c(1:4)]
comm$Date <- as.Date(comm$Date, "%m/%d/%Y")

###########################################################

# Adjust for inflation using CPI

cpi <- read.csv("C:/Users/Anshu/Box Sync/junior/econ 419/paper/cpi.csv", 
                 na.strings="NaN", stringsAsFactors=FALSE)
cpi$Date <- as.Date(cpi$Date,"%m/%d/%Y")

# CPI goes until March 2015 so that will be our benchmark. 
b <- cpi$CPI[663]
# Infl.-adj price = original(b/cpi_priceyear)
# We will have to truncate the prices at March 2015.

comm <- comm[comm$Date <= as.numeric(cpi$Date[663]),]

com <- merge(comm,cpi,by="Date")

# NOMINAL PRICES

com$Cu <- com$Cu*b/com$CPI
com$Au <- com$Au*b/com$CPI
com$Oil <- com$Oil*b/com$CPI

# ADF tests show that these are not stationary. Log. Then take 2nd difference.

com$Cu <- log(com$Cu)
com$Cu <- c(NA, NA, diff(com$Cu, differences=2))

com$Au <- log(com$Au)
com$Au <- c(NA, NA, diff(com$Au, differences=2))             

com$Oil <- log(com$Oil)
com$Oil <- c(NA, NA, diff(com$Oil, differences=2))             

com <- com[-c(1,2),]

# Lag the commodity prices for predicting next period

regs <- com
regs$Cui <- NA
regs$Aui <- NA
regs$Oili <- NA
for (i in 2:nrow(regs)) {
  regs$Cui[i] <- regs$Cu[i-1]
  regs$Aui[i] <- regs$Au[i-1]
  regs$Oili[i] <- regs$Oil[i-1]
}
regs <- regs[-1,-c(2:5)]

###########################################################

# CRSP monthly returns: used for mean model later
# I round up the monthly crsp return dates so that they are easy to handle

monthly <- read.csv("C:/Users/Anshu/Box Sync/junior/econ 419/paper/crsp monthly.csv", 
                 na.strings="NaN", stringsAsFactors=FALSE)
monthly$Date <- as.Date(as.character(monthly$Date), "%Y%m%d")
monthly$Date <- monthly$Date %m+% months(1)
day(monthly$Date) <- 1

# Since we truncated the comm prices, we have to truncate returns too. Also
# returns start later than prices.

monthly <- monthly[monthly$Date <= as.numeric(regs$Date[660]),]
regs <- regs[regs$Date >= as.numeric(monthly$Date[1]),]

# More on dates: Goyal used July 1962 to December 1998, so I'll have to truncate the tails.
# We will use Aug 1962 to December 1998

nvol <- nvol[c(1:633),]
nvol <- nvol[-1,]

rregs <- regs[regs$Date <= as.numeric(as.Date("19981201",format="%Y%m%d")),]
mmonthly <- monthly[monthly$Date <= as.numeric(as.Date("19981201",format="%Y%m%d")),]
gvol <- gvol[-1,]

###########################################################

# Now create the GARCH models.

# 1. GARCH-M

# Vanilla
library(rugarch)
specm <- ugarchspec(variance.model=list(model="sGARCH"),
                    mean.model=list(armaOrder=c(0,0),archm=TRUE,archpow=2))
fitm <- ugarchfit(data=mmonthly[,2],spec=specm)
# Commodities
specmc <- ugarchspec(variance.model=list(model="sGARCH"),
                    mean.model=list(armaOrder=c(0,0), 
                    external.regressors=as.matrix(rregs[,-1]),archm=TRUE,archpow=2))
fitmc <- ugarchfit(data=mmonthly[,2],spec=specmc)

# 2. Exponential GARCH

# Vanilla
spece <- ugarchspec(variance.model=list(model="eGARCH"),
                    mean.model=list(armaOrder=c(0,0)))
fite <- ugarchfit(data=mmonthly[,2],spec=spece)
# Commodities
specec <- ugarchspec(variance.model=list(model="eGARCH"),
                     mean.model=list(armaOrder=c(0,0), external.regressors=as.matrix(rregs[,-1])))
fitec <- ugarchfit(data=mmonthly[,2],spec=specec)

# 3. GJR-GARCH
  
# Vanilla
specg <- ugarchspec(variance.model=list(model="gjrGARCH"),
                      mean.model=list(armaOrder=c(0,0)))
fitg <- ugarchfit(data=mmonthly[,2],spec=specg)
# Commodities
specgc <- ugarchspec(variance.model=list(model="gjrGARCH"),
                     mean.model=list(armaOrder=c(0,0), external.regressors=as.matrix(rregs[,-1])))
fitgc <- ugarchfit(data=mmonthly[,2],spec=specgc)

# I did the same for the nvol longer time horizon, substituting 'monthly' and 'regs' for 'mmonthly'
# and 'rregs' 

###########################################################

# Extract and compare estimated volatilities

sfitm <- as.data.frame(sigma(fitm))
sfitmc <- as.data.frame(sigma(fitmc))
sfite <- as.data.frame(sigma(fite))
sfitec <- as.data.frame(sigma(fitec))
sfitg <- as.data.frame(sigma(fitg))
sfitgc <- as.data.frame(sigma(fitgc))

library(reshape)

# cleaning

list <- list()
list[[1]] <- sfitm
list[[2]] <- sfitmc
list[[3]] <- sfite
list[[4]] <- sfitec
list[[5]] <- sfitg
list[[6]] <- sfitgc

for (i in 1:6) {
  list[[i]] <- rename(list[[i]], c(V1="estvol"))
  list[[i]][,1] <- list[[i]][,1]^2
  rownames(list[[i]]) <- NULL
}

###########################################################
# Regress these estimated volatilities against gvol$mvol, our proxy.
# Generate White's heteroskedastic-robust std errors

library(RCurl)
library(httr)
url_robust <- "https://raw.githubusercontent.com/IsidoreBeautrelet/economictheoryblog/master/robust_summary.R"
eval(parse(text = getURL(url_robust, ssl.verifypeer = FALSE)),
     envir=.GlobalEnv)
library(gdata) 
library(zoo)


lm1 <- lm (gvol$mvol - mean(gvol$mvol) ~ list[[1]]$estvol)
lm2 <- lm (gvol$mvol - mean(gvol$mvol) ~ list[[2]]$estvol)
lm3 <- lm (gvol$mvol - mean(gvol$mvol) ~ list[[3]]$estvol)
lm4 <- lm (gvol$mvol - mean(gvol$mvol) ~ list[[4]]$estvol)
lm5 <- lm (gvol$mvol - mean(gvol$mvol) ~ list[[5]]$estvol)
lm6 <- lm (gvol$mvol - mean(gvol$mvol) ~ list[[6]]$estvol)

summary(lm1,robust = T)

lm11 <- lm (nvol$mvol - mean(nvol$mvol) ~ list[[1]]$estvol)
lm22 <- lm (nvol$mvol - mean(nvol$mvol) ~ list[[2]]$estvol)
lm33 <- lm (nvol$mvol - mean(nvol$mvol) ~ list[[3]]$estvol)
lm44 <- lm (nvol$mvol - mean(nvol$mvol) ~ list[[4]]$estvol)
lm55 <- lm (nvol$mvol - mean(nvol$mvol) ~ list[[5]]$estvol)
lm66 <- lm (nvol$mvol - mean(nvol$mvol) ~ list[[6]]$estvol)


# Now let's output these results.

lmOut <- function(res, file="test.csv", ndigit=3, writecsv=T) {
  # If summary has not been run on the model then run summary
  if (length(grep("summary", class(res)))==0) res <- summary(res, robust = T)
  co <- res$coefficients
  nvar <- nrow(co)
  ncol <- ncol(co)
  f <- res$fstatistic
  formatter <- function(x) format(round(x,ndigit),nsmall=ndigit)
  # This sets the number of rows before we start recording the coefficients
  nstats <- 4
  # G matrix stores data for output
  G <- matrix("", nrow=nvar+nstats, ncol=ncol+1)
  G[1,1] <- toString(res$call)
  # Save rownames and colnames
  G[(nstats+1):(nvar+nstats),1] <- rownames(co)
  G[nstats, 2:(ncol+1)] <- colnames(co)
  # Save Coefficients
  G[(nstats+1):(nvar+nstats), 2:(ncol+1)] <- formatter(co)
  # Save F-stat
  G[1,2] <- paste0("F(",f[2],",",f[3],")")
  G[2,2] <- formatter(f[1])
  # Save F-p value
  G[1,3] <- "Prob > P"
  G[2,3] <- formatter(1-pf(f[1],f[2],f[3]))
  # Save R2
  G[1,4] <- "R-Squared"
  G[2,4] <- formatter(res$r.squared)
  # Save Adj-R2
  G[1,5] <- "Adj-R2"
  G[2,5] <- formatter(res$adj.r.squared)
  print(G)
  if (writecsv) write.csv(G, file=file, row.names=F)
}
lmOut(lm1, file="C:/Users/Anshu/Box Sync/junior/econ 419/paper/garchmg1.csv")
lmOut(lm2, file="C:/Users/Anshu/Box Sync/junior/econ 419/paper/garchmc1.csv")
lmOut(lm3, file="C:/Users/Anshu/Box Sync/junior/econ 419/paper/garcheg1.csv")
lmOut(lm4, file="C:/Users/Anshu/Box Sync/junior/econ 419/paper/garchec1.csv")
lmOut(lm5, file="C:/Users/Anshu/Box Sync/junior/econ 419/paper/garchgg1.csv")
lmOut(lm6, file="C:/Users/Anshu/Box Sync/junior/econ 419/paper/garchgc1.csv")

###########################################################