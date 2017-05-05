# Simplified replication file for purposes of sharing and archiving. 
# Takes data previously saved to EWP_final_data_31Oct2016.RData
# and keeps necessary files in same directory as this file.

rm(list=ls())
personaldir="~/Dropbox/chad_transfer/USHMM/EWP/jaystuff2016/r/forreplication"
setwd(personaldir)

# Load required packages
library(randomForest)
library(DataCombine)
library(dplyr)
library(lubridate)

# Set seed
set.seed(20912)

# Load the compiled and transformed data
dat = load("EWP_final_data_31Oct2016.RData")

# Original MODEL FORMULAE FOR EWP STATISTICAL RISK ASSESSMENTS
# 2014-09-07
f.coup <- formula(cou.a.d.1 ~
                    imr.normed.ln + 
                    mev.regac.ln +
                    postcw +
                    cou.tries5d +
                    pol.cat.fl2 + pol.cat.fl3 + pol.cat.fl7 +
                    pol.durable.ln +
                    gdppcgrow.sr +
                    mev.civtot.ln +
                    ios.iccpr1)

f.cwar <- formula(pit.cwar.onset.1 ~
                    wdi.popsize.ln + 
                    imr.normed.ln +
                    pol.cat.fl2 + pol.cat.fl3 + pol.cat.fl7 +
                    mev.regac.ln +
                    mev.civtot.ln)

f.threat <- formula(mkl.start.1 ~ log(coup.p) + log(cwar.p))

f.pitf <- formula(pit.any.onset.1 ~
                    reg.eap + reg.eur + reg.mna + reg.sca + reg.amr +
                    imr.normed.ln +
                    pol.cat.pitf2 + pol.cat.pitf3 + pol.cat.pitf4 + pol.cat.pitf5 + pol.cat.pitf6 +
                    dis.l4pop.ln +
                    mev.regac.ln)

f.harff <- formula(mkl.start.1 ~
                     log(pitf.p) +
                     mkl.ever +
                     pit.sftpuhvl2.10.ln +
                     pol.autoc +
                     elc.eliti +
                     elc.elethc +
                     wdi.trade.ln)

f.rf <- formula(as.factor(mkl.start.1) ~
                  reg.afr + reg.eap + reg.eur + reg.mna + reg.sca + reg.amr +
                  mkl.ongoing +
                  mkl.ever +
                  countryage.ln + 
                  wdi.popsize.ln + 
                  imr.normed.ln +
                  gdppcgrow.sr +
                  wdi.trade.ln +
                  ios.iccpr1 +
                  postcw +
                  pol.cat.fl1 + pol.cat.fl2 + pol.cat.fl3 + pol.cat.fl7 +
                  pol.durable.ln +
                  dis.l4pop.ln +
                  elf.ethnicc1 + elf.ethnicc2 + elf.ethnicc3 + elf.ethnicc9 +
                  elc.eleth1 + elc.eleth2 +
                  elc.eliti +
                  cou.tries5d +
                  pit.sftpuhvl2.10.ln +
                  mev.regac.ln +
                  mev.civtot.ln)


####################################
# Model Estimation
####################################

# ATTENTION: using dat.retro.isr for dat below
dat=dat.retro.isr
rm(dat.retro.isr)

# Remove most recent 2 yrs to avoid treating in-sample ests as forecasts in cases with missing recent data
thisyear=2016
subdat = filter(dat, year < thisyear - 2)

coup <- glm(f.coup, family = binomial, data = subdat, na.action = na.exclude)
dat$coup.p <- predict(coup, newdata = dat, type = "response")
subdat$coup.p <- predict(coup, newdata = subdat, type = "response")

cwar <- glm(f.cwar, family = binomial, data = subdat, na.action = na.exclude)
dat$cwar.p <- predict(cwar, newdata = dat, type = "response")
subdat$cwar.p <- predict(cwar, newdata = subdat, type = "response")

threat <- glm(f.threat, family = binomial, data = subdat, na.action = na.exclude)
dat$threat.p <- predict(threat, newdata = dat, type = "response")

pitf <- glm(f.pitf, family = binomial, data = subset(subdat, pit.any.ongoing==0), na.action = na.exclude)
dat$pitf.p <- predict(pitf, newdata = dat, type = "response")
subdat$pitf.p <- predict(pitf, newdata = subdat, type = "response")

harff <- glm(f.harff, family = binomial, data = subdat, na.action = na.exclude)
dat$harff.p <- predict(harff, newdata = dat, type = "response")

#Original RF had 1000 trees; increase for stability when you have time.
# 20k is good for publication quality results.
rf <- randomForest(f.rf, data = subdat, na.action = na.exclude, ntree = 1000,
  mtry = 3, cutoff = c(0.2,0.8)) # Params selected thru gridded search in validation stage a couple of yrs ago
dat$rf1000.p <- predict(rf1000, type = "prob", newdata = dat, na.action = na.exclude)[,2]

###############################################
# Predictions Using Last Valid Observations
###############################################

# Create vector of labels of all variables used in all model formulae with ID vars at the front. The 'unique'
# part eliminates duplicates, and the [[3]] part restricts the call to features/independent variables. If you
# want to include targets/dependent variables from all models, too, just delete the [[3]]s.
varlist <- c(c("country", "sftgcode", "year"), unique(c(all.vars(f.coup[[3]]), all.vars(f.cwar[[3]]),
  all.vars(f.threat[[3]]), all.vars(f.pitf[[3]]), all.vars(f.harff[[3]]), all.vars(f.rf[[3]]))))

# A function to get the last non-missing value from 'dat' for each model variable for a given country code
pull <- function(i) {
  d <- subset(dat, sftgcode==i, select = varlist)
  snip <- function(x) { rev( na.omit(x) )[1] }
  xi <- lapply(d, snip)
  xi <- as.data.frame(xi)
  return(xi)
}

# Apply that function to a list of country codes and rbind the resulting vectors
latest <- as.list(unique(as.character(dat$sftgcode))) %>%
    lapply(., pull) %>%
    Reduce(function(x, y) { rbind(x, y) }, .)

# get rid of rows for Taiwan, countries that are too small, and countries that don't exist any more
latest <- latest[-which(latest$sftgcode %in% c("BHM", "BAR", "BLZ", "CAP", "TAW")),]
latest <- latest[-which(latest$year < 2015),]

# Hard code a few missing values for North Korea
# trade volume estimates $2.7 and $3.6 billion from
# http://www.nkeconwatch.com/2016/06/14/north-koreas-trade-volume-down-18-percent-in-2015/
# gdp estimate of 40.0 billion from CIA World Factbook
latest$wdi.trade.ln[latest$sftgcode == "PRK"] <- log(100 * ((2.7 + 3.6)/40.0))
# growth rate estimate (2014): http://www.nkeconwatch.com/nk-uploads/GDP_of_North_Korea_in_2014_ff.pdf
latest$gdppcgrow.sr[latest$sftgcode == "PRK"] <- sqrt(abs(1))
# population size from CIA World Factbook: 24,983,205 as mid-2015 est.
latest$wdi.popsize.ln[latest$sftgcode == "PRK"] <- log(24983205 / 1000)

# A function to get predictions from a data frame with the variables in model.formulae.r
forecast <- function(df) {

  df$coup.p <- predict(coup, newdata = df, type = "response", na.action = na.exclude)
  df$cwar.p <- predict(cwar, newdata = df, type = "response", na.action = na.exclude)
  df$pitf.p <- predict(pitf, newdata = df, type = "response", na.action = na.exclude)
  df$threat.p <- predict(threat, newdata = df, type = "response", na.action = na.exclude)
  df$harff.p <- predict(harff, newdata = df, type = "response", na.action = na.exclude)
  df$rf.p <- predict(rf, type = "prob", newdata = df, na.action = na.exclude)[,2]
  df$mean.p <- (df$harff.p + df$threat.p + df$rf.p)/3
  df$date <- as.Date(Sys.Date())
  return(df)
}

# Apply that function to latest available data to get country forecasts. NOTE: If this produced an error
# message about the number of rows not matching, it means there are still some missing values in 'latest'.
# The simplest way to spot those missing values is to call fix(latest) and scroll around the table in search
# of NAs. Every NA needs to be addressed for this function to work. In future updates, some missing values
# may best be addressed at the data ingestion stage (where the 'data.[prefix].r' scripts work, but the odd
# case can also be addressed by adding a line above to hard code the relevant input if need be.
newcast <- forecast(latest)

# Create variable identifying year to which forecasts apply (last year in original data plus 1)
newcast$forecast.year <- max(dat$year) + 1

# Move ID variables and forecasts to front of data frame
newcast <- MoveFront(newcast, c("country", "sftgcode", "forecast.year", "mean.p", "rf.p", "threat.p", "harff.p"))

# Put in descending order by mean.p
newcast <- newcast[order(-newcast$mean.p),]

# Write it out
#write.csv(newcast, file = "cjh12sep2016.csv"), row.names = FALSE)
