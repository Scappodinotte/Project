# ------------------------------------------------------------------------------
# Applied Macroeconometrics Project
# Author: Elia Scapini
# Date: 20/11/2024
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# GitHub
# ------------------------------------------------------------------------------
library(usethis)
use_git()
use_github()


# ------------------------------------------------------------------------------
# Load packages and functions
# ------------------------------------------------------------------------------
library(tsbox)
library(forecast) # Useful package for ARMA models
library(quantmod)
library(xts)
library(ggplot2)   
library(seasonal)
library(CADFtest) # Useful package for conducting unit root tests
library(reshape2) # Useful to reshape data frames (for some specific plots)
library(lubridate)

# Delete all objects in the memory
rm(list=ls())

# Load user-defined commands and packages
source("UserPackages.R")

# Here, it executes the function and creates the folder in the current directory
setwd("~/Personale/UNINE/Master_Applied_Economics/Sem1/Macro/Project")
mainDir <- getwd()
outDir <- makeOutDir(mainDir, "/ResultsProject")

# ------------------------------------------------------------------------------
# Import data
# ------------------------------------------------------------------------------
UR <- read.csv("Data/LRHU24TTIEM156S.csv", header = T, sep = ",")
UR <- xts(UR[, 2], order.by = as.Date(UR[, 1], format = "%d/%m/%Y"))

# ------------------------------------------------------------------------------
# Plotting and data transformation
# ------------------------------------------------------------------------------
# Plotting the TS
ts_plot(
  `Youth Unemployment Rate`= UR,
  title = "Youth Unemployment Rate",
  subtitle = "Percentage, seasonally adjusted"
)
ts_save(filename = paste(outDir, "Youth_UR.pdf", sep = "/"), width = 8, height = 7, open = FALSE)
# No apparent visual trend spotted 

# Check uni-variability
urootUR_drift = CADFtest(UR, max.lag.y = 10, type = "drift", criterion = "BIC")
summary(urootUR_drift)
# p-value < 0.05. So reject the null hypothesis (non-stationary around a constant mean, phi = 1)
# The TS is CSP

# UR <- ts_span(UR, start = NULL, end = NULL)

# p <- plotACF(UR, lag.max = 24)+theme_minimal()
# p <- ggLayout(p)
# p
# ggsave(paste(outDir, "URAcf.pdf", sep = "/"), plot = last_plot(), width = 12, height = 9, units = c("cm")) 

# autoplot(decompose(ts_ts(UR), "additive"))

# ------------------------------------------------------------------------------
# Model selection and diagnostic
# ------------------------------------------------------------------------------
model1 = auto.arima(UR, max.p = 6, max.q = 6, d = 0, ic = c("bic"), allowmean = TRUE, seasonal = FALSE, stepwise = FALSE)
summary(model1)
# The best model from the auto arima function is ARMA(1,4)

checkresiduals(model1)
# Residuals show some significant Autocorrelation, which means that there is room for improvement

# maxP <- 6   # Maximum number of AR lags
# maxQ <- 6   # Maximum number of MA lags
# 
# # Objects to save the criteria for every possible lag structure
# AIC = matrix(data=NA,nrow=maxP+1,ncol=maxQ+1)
# BIC = matrix(data=NA,nrow=maxP+1,ncol=maxQ+1)
# colnames(AIC) <- 0:maxQ
# colnames(BIC) <- 0:maxQ
# rownames(AIC) <- 0:maxP
# rownames(BIC) <- 0:maxP
# 
# # Loop over all possible lag orders, estimate the model, and save the criteria
# # in a matrix
# for (p in 0:maxP){
#   for (q in 0:maxQ){
# 
#     # Estimate the corresponding model
#     temp <- Arima(UR, order = c(p, 0, q), include.constant= TRUE)
# 
#     # Save the information criterion
#     AIC[p+1, q+1] <- temp$aic
#     BIC[p+1, q+1] <- temp$bic
#   }
# }
# 
# # Find the lag order with the smallest information criterion (nameMin() is a
# # user-defined function that is useful to find the name of the columns and rows
# # of the minimum value of a matrix)
# minCritAIC <- nameMin(AIC)
# minCritBIC <- nameMin(BIC)
# 
# print("Lag order according to AIC and BIC (p, q)")
# minCritAIC
# minCritBIC

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
