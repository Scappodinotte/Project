# ------------------------------------------------------------------------------
# Applied Macroeconometrics Project
# Author: Elia Scapini
# Date: 20/11/2024
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Questions professor
# ------------------------------------------------------------------------------
# 1) shrinking the data? Because before old equilibrium. What about Subprime + Euro crisis + COVID effect?
#  or use dummies to controll for the crisis
# 2) Create a training and test data? 
# 3) Build model with monthly or quarterly data?
# 4) Reporting: transform in year or quarter TS, but no yoy growth rate or log, right?


# ------------------------------------------------------------------------------
# Connecting GitHub
# ------------------------------------------------------------------------------
# library(usethis)
# use_git()
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
# outDir <- makeOutDir(mainDir, "/ResultsProject")

# ------------------------------------------------------------------------------
# Import data
# ------------------------------------------------------------------------------
UR <- read.csv("Data/LRHU24TTIEM156S.csv", header = T, sep = ",")
UR <- xts(UR[, 2], order.by = as.Date(UR[, 1], format = "%d/%m/%Y"))
URq <- ts_frequency(UR, to = "quarter", aggregate = "last")
# UR <- log(UR)
# Mean with log: 15.93

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
g <- g + xlab("Years") + ylab("Youth Unemployment Rate [%]") + ggtitle("15-24 y/o Unemployment Rate", subtitle = "Quartely, seasonally adjusted")
g
ggsave(paste(outDir,"Youth_UR.pdf", sep = "/"), plot = last_plot(), width = 10, height = 8, units = c("cm"))
# Discussion: No apparent visual trend spotted 

# ts_plot(
#   `Youth Unemployment Rate`= UR,
#   title = "15-24 y/o Unemployment Rate",
#   subtitle = "Percentage, seasonally adjusted",
#   ylab = "Youth Unemployment rate [%]"
# )
# ts_save(filename = paste(outDir, "Youth_UR.pdf", sep = "/"), width = 8, height = 7, open = FALSE)

# Check uni-variability
urootUR_drift = CADFtest(UR, max.lag.y = 10, type = "drift", criterion = "BIC")
summary(urootUR_drift)
urootUR_drift = CADFtest(ts_diff(UR), max.lag.y = 10, type = "drift", criterion = "BIC")
summary(urootUR_drift)

# Taking simple difference 
URd <- ts_diff(UR)
URd[1] <- 0
# Discussion: p-value < 0.05. So reject the null hypothesis (non-stationary around a constant mean, phi = 1)
# The TS is non CSP, after the first differenve it is. 

train <- ts_span(URd, start = NULL, end = "2011-12-01")
test <- ts_span(URd, start = "2012-01-01", end = NULL)
# ------------------------------------------------------------------------------
# Model selection and diagnostic
# ------------------------------------------------------------------------------
maxP <- 8   # Maximum number of AR lags
maxQ <- 8   # Maximum number of MA lags

model1 <- auto.arima(URd, max.p = maxP, max.q = maxQ, d = 0, ic = c("bic"), allowmean = TRUE, seasonal = FALSE, stepwise = FALSE)
summary(model1)
# The best model from the auto arima function is ARMA(4,1)

checkresiduals(model1)
# Residuals show some small significant Autocorrelation at t=5, 7 and 21, which means that there is room for improvement.

# Objects to save the criteria for every possible lag structure
BIC = matrix(data=NA,nrow=maxP+1,ncol=maxQ+1)
colnames(BIC) <- 0:maxQ
rownames(BIC) <- 0:maxP

# Loop over all possible lag orders, estimate the model, and save the criteria
# in a matrix
for (p in 0:maxP){
  for (q in 0:maxQ){

    # Estimate the corresponding model
    temp <- Arima(URd, order = c(p, 0, q), include.constant= TRUE)

    # Save the information criterion
    BIC[p+1, q+1] <- temp$bic
  }
}

# Find the lag order with the smallest information criterion (nameMin() is a
# user-defined function that is useful to find the name of the columns and rows
# of the minimum value of a matrix)
minCritBIC <- nameMin(BIC)

print("Lag order according to AIC and BIC (p, q)")
minCritBIC
model2 <- Arima(URd, order = c(7, 0, 0), include.constant= TRUE)
summary(model2)
checkresiduals(model2)
# Autocorrelation at t = 5 and 7 has disappear. now auto at t = 21

# ------------------------------------------------------------------------------
# Forecasting 15 months (from 10/2024, to 12/2025)
# ------------------------------------------------------------------------------
# Step in forecast using AutoArima
maxh <- 15

auto_forecast <- forecast(model1, h = maxh)
autoplot(auto_forecast) + theme_minimal()

# Direct forecast
x1 <- UR

start <- ts_summary(UR)$end + months(1)   # Start date of the forecast
end   <- ts_summary(UR)$end + months(maxh)   # End date of the forecast

# Create a time series to save the forecast in
Forecasth <- xts(rep(NA, maxh), order.by = seq(start, end, by = "1 month"))

for (h in 1:maxh){
  
  Yh <- ts_lag(UR, -h) 
  
  # Collect everything in a data frame to estimate the model
  Data <- data.frame(ts_c(Yh, x1))
  DirectM <- lm(Yh ~ x1, data = Data)
  
  # The forecast can be extracted using the last fitted value of the series
  # Save it in the correct date of the forecast series
  Fitted       <- fitted(DirectM)
  Forecasth[h] <- Fitted[length(Fitted)]
  
}

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------

# Merge the history and the forecast. Replace first observation with zero 
# because of missing value
temp1 <- ts_bind(auto_forecast$x, auto_forecast$mean)
# temp1[1] = 0

# Step 4: Split series in history and forecast for the plot
fcst_start <- ts_summary(auto_forecast$mean)$start
hist_end   <- ts_summary(auto_forecast$x)$end
Hist   <- ts_span(temp1, NULL, hist_end)
Fcst   <- ts_span(temp1, fcst_start, NULL)

ts_plot(
  `History`= Hist,
  `Forecast` = Fcst,
  title = "Ireland UR and forecast",
  subtitle = "Percentage"
)
ts_save(filename = paste(outDir, "/CPIFcstyoy.pdf", sep = ""), width = 8, height = 5, open = FALSE)
