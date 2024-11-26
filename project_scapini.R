# ------------------------------------------------------------------------------
# Applied Macroeconometrics Project
# Author: Elia Scapini
# Date: 20/11/2024
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Questions professor
# ------------------------------------------------------------------------------
# Confidence interval
# Mean reverting, but which one
# Create a training and test data? --> optional estimate the model with the whole TS

# ------------------------------------------------------------------------------
# Connecting GitHub
# ------------------------------------------------------------------------------
library(usethis)
use_git()
# use_github()

# ------------------------------------------------------------------------------
# Load packages and functions
# ------------------------------------------------------------------------------
library(tsbox)
library(forecast)
library(quantmod)
library(xts)
library(ggplot2)   
library(seasonal)
library(CADFtest)
library(reshape2) 
library(lubridate)

# Delete all objects in the memory
rm(list=ls())

# Load user-defined commands and packages
source("UserPackages.R")

# Create the folder in the current directory
setwd("~/Personale/UNINE/Master_Applied_Economics/Sem1/Macro/Project")
mainDir <- getwd()
outDir <- makeOutDir(mainDir, "/ResultsProject")

# ------------------------------------------------------------------------------
# Import data
# ------------------------------------------------------------------------------
UR <- read.csv("Data/LRHU24TTIEM156S.csv", header = T, sep = ",")
UR <- xts(UR[, 2], order.by = as.Date(UR[, 1], format = "%d/%m/%Y"))

start_date <- index(UR)[1]
# UR <- ts_span(UR, start = start_date, end = NULL)
# Discussion: mean from 1983: 17.46%. Mean from 1999: 15.15%

# Official NBER recessions from http://www.nber.org/cycles.html
# These will be plotted along with the indicator
NBERREC <- read.table(textConnection(
  "Peak, Trough
  1857-06-01, 1858-12-01
  1860-10-01, 1861-06-01
  1865-04-01, 1867-12-01
  1869-06-01, 1870-12-01
  1873-10-01, 1879-03-01
  1882-03-01, 1885-05-01
  1887-03-01, 1888-04-01
  1890-07-01, 1891-05-01
  1893-01-01, 1894-06-01
  1895-12-01, 1897-06-01
  1899-06-01, 1900-12-01
  1902-09-01, 1904-08-01
  1907-05-01, 1908-06-01
  1910-01-01, 1912-01-01
  1913-01-01, 1914-12-01
  1918-08-01, 1919-03-01
  1920-01-01, 1921-07-01
  1923-05-01, 1924-07-01
  1926-10-01, 1927-11-01
  1929-08-01, 1933-03-01
  1937-05-01, 1938-06-01
  1945-02-01, 1945-10-01
  1948-11-01, 1949-10-01
  1953-07-01, 1954-05-01
  1957-08-01, 1958-04-01
  1960-04-01, 1961-02-01
  1969-12-01, 1970-11-01
  1973-11-01, 1975-03-01
  1980-01-01, 1980-07-01
  1981-07-01, 1982-11-01
  1990-07-01, 1991-03-01
  2001-03-01, 2001-11-01
  2007-12-01, 2009-06-01
  2020-02-01, 2020-04-01"), sep = ',',
  colClasses = c('Date', 'Date'), header = TRUE)
NBERREC <- subset(NBERREC, Peak >= as.Date(start_date))
NBERREC2 <- data.frame(Peak = as.Date("2011-11-01"), Trough = as.Date("2013-01-01"))

# ------------------------------------------------------------------------------
# Plotting and data transformation
# ------------------------------------------------------------------------------
# Plotting the TS
g <- ggplot(UR) + geom_line(aes(x = index(UR), y = UR)) + theme_minimal()
g <- g + geom_rect(data = NBERREC, aes(xmin = Peak, xmax = Trough, ymin = -Inf, ymax = +Inf), fill = 'grey', alpha = 0.5)
g <- g + geom_rect(data = NBERREC2, aes(xmin = Peak, xmax = Trough, ymin = -Inf, ymax = +Inf), fill = 'grey', alpha = 0.5)
g <- g + xlab("Years") + ylab("Youth Unemployment Rate [%]") + ggtitle("15-24 y/o Unemployment Rate", subtitle = "Monthly, seasonally adjusted")
g
ggsave(paste(outDir,"Youth_UR.pdf", sep = "/"), plot = last_plot(), width = 10, height = 8, units = c("cm"))
# Discussion: No apparent visual trend spotted 

# Check uni-variability
urootUR_drift = CADFtest(UR, max.lag.y = 10, type = "drift", criterion = "BIC")
summary(urootUR_drift)
# Discussion: the TS is not CSP
urootUR_drift = CADFtest(ts_diff(UR), max.lag.y = 10, type = "drift", criterion = "BIC")
summary(urootUR_drift)
# Discussion: once taking the first difference, TS is CSP

# Taking simple difference 
URd <- ts_ts(ts_diff(UR))

# Check the autocorrelation of the series
plotACF(URd, lag.max = 24)
# Discussion: This suggest a seasonality at quarterly frequency

# Creating the train and test data
# train <- ts_span(URd, start = NULL, end = "2011-12-01")
# test <- ts_span(URd, start = "2012-01-01", end = NULL)

# ------------------------------------------------------------------------------
# Model selection and diagnostic
# ------------------------------------------------------------------------------
maxP <- 8   # Maximum number of AR lags
maxQ <- 8   # Maximum number of MA lags

# METHOD 1: Auto Arima
model1 <- auto.arima(URd, max.p = maxP, max.q = maxQ, d = 0, ic = c("bic"), allowmean = TRUE, seasonal = FALSE, stepwise = FALSE)
summary(model1)
# The best model from the auto arima function is ARMA(4,1)

checkresiduals(model1)
# Residuals show some small significant Autocorrelation at t=5, 7 and 21, which means that there is room for improvement.

# # METHOD 2: Manual selection
# # Objects to save the criteria for every possible lag structure
# BIC = matrix(data = NA,nrow = maxP + 1, ncol = maxQ + 1)
# colnames(BIC) <- 0:maxQ
# rownames(BIC) <- 0:maxP
# 
# Loop over all possible lag orders, estimate the model, and save the criteria
# for (p in 0:maxP){
#   for (q in 0:maxQ){
# 
#     # Estimate the corresponding model
#     temp <- Arima(URd, order = c(p, 0, q), include.constant = TRUE)
# 
#     # Save the information criterion
#     BIC[p + 1, q + 1] <- temp$bic
#   }
# }
# 
# # Find the lag order with the smallest information criterion
# minCritBIC <- nameMin(BIC)
# 
# print("Lag order according to AIC and BIC (p, q)")
# 
# minCritBIC

model2 <- Arima(URd, order = c(7, 0, 0), include.constant = TRUE)
summary(model2)
checkresiduals(model2)
# Autocorrelation at t = 5 and 7 has disappear. now auto at t = 21

# ------------------------------------------------------------------------------
# Forecasting 15 months (from 10/2024, to 12/2025)
# ------------------------------------------------------------------------------
# Step in forecast using AutoArima
maxh <- 15

forecast <- forecast(model2, h = maxh)
autoplot(forecast) + theme_minimal()

# Merge the history and the forecast. Replace first observation with the first value of the UR 
# because of missing value
temp1 <- ts_bind(forecast$x, forecast$mean)
temp1_lower <- ts_bind(forecast$x, forecast$lower[, 2])
temp1_upper <- ts_bind(forecast$x, forecast$upper[, 2])

# Set first value to the first valeu of UR
temp1[1] <- UR[1]
temp1_lower[1] <- UR[1]
temp1_upper[1] <- UR[1]

# Calculate the cumulative sum of the log-differences 
temp2 <- ts_cumsum(temp1)
temp2_lower <- ts_cumsum(temp1_lower)
temp2_upper <- ts_cumsum(temp1_upper)

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------
# Split series in history and forecast for the plot
start_span <- NULL

fcst_start <- ts_summary(forecast$mean)$start
hist_end   <- ts_summary(forecast$x)$end
HistUR <- ts_span(temp2, start = start_span, hist_end)
FcstUR <- ts_span(temp2, fcst_start, NULL)
lower <- ts_span(temp2_lower, fcst_start, NULL)
upper <- ts_span(temp2_upper, fcst_start, NULL)

ts_plot(
  `History`= HistUR,
  `Forecast` = FcstUR,
  `Lower` = lower,
  `Upper` = upper,
  title = "Irish Youth Unemployment Rate and forecast",
  subtitle = "Percentage"
)
ts_save(filename = paste(outDir, "/FrcstUR.pdf", sep = ""), width = 8, height = 5, open = FALSE)

# ------------------------------------------------------------------------------
# Compute forecast error variance and simulate forecast density
# ------------------------------------------------------------------------------
# Calculate variance
sigma2 <- getForecastVariance(forecast)
sigma2

# MONTE CARLO SIMULATION : Simulate a normal distribution around each forecast point
NSim <- 1000   # Number of simulations (S)
H    <- 15   # Maximum forecast horizon
fcsth <- forecast$mean
SimFcst <- matrix(NA, nrow = H, ncol = NSim)
for (h in 1:H){
  SimFcst[h, ] <- rnorm(NSim, fcsth[h], sqrt(sigma2[h]))
}
SimFcst <- xts(SimFcst, order.by = as.Date(index(ts_xts(fcsth))))

# Check normality of Monte Carlo simulations
hist(as.numeric(SimFcst[1, ]))

plot(SimFcst)
ci95 <- t(apply(SimFcst, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE))
ci95 <- xts(ci95, order.by = as.Date(index(ts_xts(fcsth))))

all_ci <- ts_c(ci95,forecast$lower[, "95%"], forecast$upper[, "95%"])
print(all_ci)
plot(all_ci)

# ------------------------------------------------------------------

#    Calculate probability of negative GDP growth in Q1 2021
PNeg = mean(SimFcst["2025-01-01",]<0)

# Calculate probability of GDP growth larger than 2% in Q1 2021
PLarge = mean(SimFcst["2025-01-01",]>2)

# Calculate probability of GDP growth between 0 and 2% in Q1 2021
PBetween = mean(SimFcst["2025-01-01",]<=2 & SimFcst["2025-01-01",]>=0)

PNeg
PLarge
PBetween
PNeg + PLarge + PBetween

# Do a nice bar chart for over-heating and recession probabilities
PNeg = rowMeans(SimFcst<0)
PLarge = rowMeans(SimFcst>2)
PNorm = 1-PNeg-PLarge

prob.data = data.frame(index(SimFcst), PLarge, PNorm, PNeg)
colnames(prob.data) = c("Date", "Boom", "Normal", "Negative")
prob.data <- melt(prob.data,id.vars = "Date") 
ggplot(prob.data, aes(x = Date, y = value,fill=variable)) +
  geom_bar(stat='identity')+theme_minimal()+ggtitle("Probability of a boom, recession, and normal state")
ggsave(paste(outDir,"/ProbChart.pdf", sep = ""), plot = last_plot(), width = 10, height = 8, units = c("cm"))

# d) Compute probability that annual GDP growth in 2020 is negative
#    Basically, we undo the transformation we applied before and
#    then aggregate the forecast and the history to annual frequency
#    See Appendix of Ch. 5 slides for a description of the procedure.

# Step 1: Calculate the level for each simulation and aggregate to annual 
#         frequency
for (s in 1:NSim){
  
  # Add history to forecast
  temp = ts_bind(URd, SimFcst[, s])
  temp[1] <- 0

  # Aggregate to annual frequency
  temp = ts_frequency(temp, to ="year", aggregate="sum")
  # 
  # # Calculate growth rates of annual series
  # temp = ts_pc(temp)
  
  # Save result in collection of time series
  if (s == 1){
    SimFcstGrt = temp
  }else{
    SimFcstGrt = ts_c(SimFcstGrt, temp)
  }
}

# Step 2: Calculate the probability
PNeg = mean(SimFcstGrt["2021-01-01",]<0)
print(PNeg)
PNeg = mean(SimFcstGrt["2022-01-01",]<0)
print(PNeg)

# Plot some interval forecasts for annual GDP growth
ci95 = as.xts(t(apply(SimFcst, 1, quantile, probs = c(0.025, .5, 0.975),  na.rm = TRUE) ))
ts_plot(
  `Upper`  = ci95[,"97.5%"],
  `Median` = ci95[,"50%"],
  `Lower`  = ci95[,"2.5%"],
  `History`= HistUR,
  `Forecast` = FcstUR,
  title    = "Swiss GDP",
  subtitle = "Annual, growth rate"
)
ts_save(paste(outDir,"/AnnualGrowth.pdf", sep = ""), width = 10, height = 8, open = F)


# ts_plot(
#   `Youth Unemployment Rate`= UR,
#   title = "15-24 y/o Unemployment Rate",
#   subtitle = "Percentage, seasonally adjusted",
#   ylab = "Youth Unemployment rate [%]"
# )
# ts_save(filename = paste(outDir, "Youth_UR.pdf", sep = "/"), width = 8, height = 7, open = FALSE)

# # Direct forecast
# x1 <- UR
# 
# start <- ts_summary(UR)$end + months(1)   # Start date of the forecast
# end   <- ts_summary(UR)$end + months(maxh)   # End date of the forecast
# 
# # Create a time series to save the forecast in
# Forecasth <- xts(rep(NA, maxh), order.by = seq(start, end, by = "1 month"))
# 
# for (h in 1:maxh){
#   
#   Yh <- ts_lag(UR, -h) 
#   
#   # Collect everything in a data frame to estimate the model
#   Data <- data.frame(ts_c(Yh, x1))
#   DirectM <- lm(Yh ~ x1, data = Data)
#   
#   # The forecast can be extracted using the last fitted value of the series
#   # Save it in the correct date of the forecast series
#   Fitted       <- fitted(DirectM)
#   Forecasth[h] <- Fitted[length(Fitted)]
#   
# }
# 
# autoplot(Forecasth)