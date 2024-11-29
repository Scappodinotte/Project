# ------------------------------------------------------------------------------
# Applied Macroeconometrics Project
# Author: Elia Scapini
# Date: 20/11/2024
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Questions professor
# ------------------------------------------------------------------------------


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
library(writexl)

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
# Import time series
UR <- read.csv("Data/LRHU24TTIEM156S.csv", header = T, sep = ",")
UR <- xts(UR[, 2], order.by = as.Date(UR[, 1], format = "%d/%m/%Y"))

# Fix start date
start_date <- index(UR)[1]

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
ggsave(paste(outDir,"Youth_UR.jpeg", sep = "/"), plot = last_plot(), width = 15, height = 10, units = c("cm"))
# Discussion: No apparent trend visually spotted 

# Check unit root process
urootUR_drift = CADFtest(UR, max.lag.y = 10, type = "drift", criterion = "BIC")
summary(urootUR_drift)
# Discussion: the TS is not CSP
urootUR_drift = CADFtest(ts_diff(UR), max.lag.y = 10, type = "drift", criterion = "BIC")
summary(urootUR_drift)
# Discussion: once taking the first difference, TS is CSP

# Taking simple difference 
URd <- ts_ts(ts_diff(UR))

# Check for seasonality
ggsubseriesplot(ts_ts(UR))
ggsave(paste(outDir,"Seasonality.jpeg", sep = "/"), plot = last_plot(), width = 15, height = 10, units = c("cm"))
# Discussion: no seasonality spotted, means do not change over different months

# Check the autocorrelation of the series
plotACF(URd, lag.max = 24)
# Discussion: This suggest a seasonality at quarterly frequency

# ------------------------------------------------------------------------------
# Model selection and diagnostic
# ------------------------------------------------------------------------------
maxP <- 8   # Maximum number of AR lags
maxQ <- 8   # Maximum number of MA lags

# METHOD 1: Auto Arima
model1 <- auto.arima(URd, max.p = maxP, max.q = maxQ, d = 0, ic = c("bic"), allowmean = TRUE, seasonal = FALSE, stepwise = FALSE)
summary(model1)
# The best model from the auto arima function is ARMA(4,1)

# Get R squared for model 1
cor(ts_span(URd, start = "1983-02-01"), fitted(model1))^2
# Discussion: R2 is 40%

# Check residuals model 1
checkresiduals(model1)
# Residuals show some significant Autocorrelation at t=5, 7 and 21, which means that there is room for improvement.

# # METHOD 2: Automated BIC selection
# # Objects to save the criteria for every possible lag structure
# BIC = matrix(data = NA,nrow = maxP + 1, ncol = maxQ + 1)
# colnames(BIC) <- 0:maxQ
# rownames(BIC) <- 0:maxP
# 
# #Loop over all possible lag orders, estimate the model, and save the criteria
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
# minCritBIC

# Estimate the best model from automated BIC selection 
model2 <- Arima(URd, order = c(7, 0, 0), include.constant = TRUE)
summary(model2)

# Get R2 for model 2
cor(ts_span(URd, start = "1983-02-01"), ts_span(fitted(model2), start = "1983-02-01"))^2
# The R2 is 42%

# Check residuals for model 2
checkresiduals(model2)
# Discussion: Autocorrelation at t = 5 and 7 has disappear. now auto at t = 21

# ------------------------------------------------------------------------------
# Forecasting 15 months (from 10/2024, to 12/2025)
# ------------------------------------------------------------------------------
# Iterated forecast model 2
maxh <- 15
forecast <- forecast(model2, h = maxh)

# Merge the history and the forecast
temp1 <- ts_bind(forecast$x, forecast$mean)

# Set first value to the first value of UR
temp1[1] <- UR[1]

# Calculate the cumulative sum of the difference on monthly frequency
temp1 <- ts_cumsum(temp1)

# Aggregate to yearly frequency
temp1_y <- ts_frequency(temp1, to = "year", aggregate = "mean")

# ------------------------------------------------------------------------------
# Compute forecast error variance and simulate forecast density
# ------------------------------------------------------------------------------
# Calculate variance
sigma2 <- getForecastVariance(forecast)

# MONTE CARLO SIMULATION : Simulate a normal distribution around each forecast point estimate
set.seed(5)
NSim <- 100   # Number of simulations (S)
H    <- 15   # Maximum forecast horizon

fcsth <- forecast$mean
sim_fcst <- matrix(NA, nrow = H, ncol = NSim)

for (h in 1:H){
  sim_fcst[h, ] <- rnorm(NSim, fcsth[h], sqrt(sigma2[h]))
}
sim_fcst <- xts(sim_fcst, order.by = as.Date(index(ts_xts(fcsth))))

# Calculate the level for each MONTHLY simulation and aggregate 
for (s in 1:NSim){
  
  # Add history to forecast
  temp2 <- ts_bind(URd, sim_fcst[, s])
  temp2[1] <- UR[1]
  
  # Cumulative sum
  temp2 <- ts_cumsum(temp2)
  
  # Save result in collection of time series
  if (s == 1){
    sim_fcst_cum <- temp2
  }else{
    sim_fcst_cum <- ts_c(sim_fcst_cum, temp2)
  }
}

# Calculate the level for each QUARTERLY simulation and aggregate 
for (s in 1:NSim){
  
  # Add history to forecast
  temp3 <- ts_bind(URd, sim_fcst[, s])
  temp3[1] <- UR[1]
  
  # Cumulative sum
  temp3 <- ts_cumsum(temp3)
  
  # Aggregate to quarterly frequency
  temp3 <- ts_frequency(temp3, to = "quarter", aggregate = "mean")
  
  # Save result in collection of time series
  if (s == 1){
    sim_fcst_cum_q <- temp3
  }else{
    sim_fcst_cum_q <- ts_c(sim_fcst_cum_q, temp3)
  }
}

# Calculate the level for each YEARLY simulation and aggregate 
for (s in 1:NSim){
  
  # Add history to forecast
  temp4 <- ts_bind(URd, sim_fcst[, s])
  temp4[1] <- UR[1]
  
  # Cumulative sum
  temp4 <- ts_cumsum(temp4)
  
  # Aggregate to annual frequency
  temp4 <- ts_frequency(temp4, to = "year", aggregate = "mean")
  
  # Save result in collection of time series
  if (s == 1){
    sim_fcst_cum_y <- temp4
  }else{
    sim_fcst_cum_y <- ts_c(sim_fcst_cum_y, temp4)
  }
}

# Reduce the accumulated simulations to the forecast time period
fcst_start <- ts_summary(forecast$mean)$start
hist_end   <- ts_summary(forecast$x)$end
sim_fcst_cum <- ts_span(sim_fcst_cum, start = fcst_start)
sim_fcst_cum_q <- ts_span(sim_fcst_cum_q, start = fcst_start)

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------
# Split series in history and forecast for the plot
start_span <- "2022-01-01"

HistUR <- ts_span(temp1, start = start_span, hist_end)
FcstUR <- ts_span(temp1, fcst_start, NULL)

ts_plot(
  `History`= HistUR,
  `Forecast` = FcstUR,
  title = "Irish Youth Unemployment Rate and forecast",
  subtitle = "Percentage"
)
ts_save(filename = paste(outDir, "/FrcstUR.jpeg", sep = ""), width = 10, height = 8, open = FALSE)

# Plot confidence interval forecasts
ci95 <- as.xts(t(apply(sim_fcst_cum, 1, quantile, probs = c(0.025, .5, 0.975),  na.rm = TRUE)))
ts_plot(
  `History`= HistUR,
  `Forecast` = FcstUR,
  `Upper95`  = ci95[, "97.5%"],
  `Lower95`  = ci95[, "2.5%"],
  title = "Irish Youth Unemployment Rate and forecast",
  subtitle = "Percentage"
)
ts_save(paste(outDir,"/Youth_UR_forecasted.jpeg", sep = ""), width = 10, height = 8, open = F)

# ------------------------------------------------------------------------------
# Tables Reporting
# ------------------------------------------------------------------------------
table <- data.frame(t(data.frame(Year = index(temp1_y), Youth_UR = as.matrix(temp1_y))))
write_xlsx(table, paste(outDir,"/Yearly_table.xlsx", sep = ""), col_names = F)

# Report monthly confidence interval for 2024 forecast
ci95["2024-10-01", ]
ci95["2024-11-01", ]
ci95["2024-12-01", ]

# Report quarterly confidence interval for 2025 forecast
ci95_q <- as.xts(t(apply(sim_fcst_cum_q, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE)))
ci95_q["2025-01-01", ]
ci95_q["2025-04-01", ]
ci95_q["2025-07-01", ]
ci95_q["2025-10-01", ]

# Report yearly confidence interval for 2024/25 forecast
ci95_y <- as.xts(t(apply(sim_fcst_cum_y, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE)))
ci95_24 <- ci95_y["2024-01-01", ]
ci95_25 <- ci95_y["2025-01-01", ]
ci95_24
ci95_25


# ------------------------------------------------------------------------------
# Calculate Probability
# ------------------------------------------------------------------------------
# Set historical values
Pmin <- 10.1 #min(UR[seq.Date(from = as.Date("2023-01-01"), to = as.Date("2023-12-01"), by = "month")])
Pmax <- 11.6 #max(UR[seq.Date(from = as.Date("2023-01-01"), to = as.Date("2023-12-01"), by = "month")])

# Calculate probability of UR being smaller than 10.1%, which is the 2023 min 
PNeg <- mean(sim_fcst_cum_q["2025-10-01", ] < Pmin)

# Calculate probability of UR being bigger than 11.6%, which is the 2023 max
PLarge <- mean(sim_fcst_cum_q["2025-10-01", ] > Pmax)

# Calculate probability of UR being bigger than 10.1%, and smaller than 11.6%
PBetween <- mean(sim_fcst_cum_q["2025-10-01", ] <= Pmax & sim_fcst_cum_q["2025-10-01", ] > Pmin)

PNeg
PLarge
PBetween
PNeg + PBetween

# Bar chart for exceeding 2023 values
PNeg <- rowMeans(sim_fcst_cum_q < Pmin)
PLarge <- rowMeans(sim_fcst_cum_q > Pmax)
PNorm <- 1 - PNeg - PLarge

prob.data <- data.frame(index(sim_fcst_cum_q), PLarge, PNorm, PNeg)
colnames(prob.data) = c("Date", "Higher", "Normal", "Lower")
prob.data <- melt(prob.data,id.vars = "Date") 
ggplot(prob.data, aes(x = Date, y = value, fill = variable)) +
  geom_bar(stat = 'identity') + theme_minimal() + 
  ggtitle("Probability of a higher, lower, and normal YUR compared to 2023")
ggsave(paste(outDir,"/ProbChart.jpeg", sep = ""), plot = last_plot(), width = 16, height = 11, units = c("cm"))

