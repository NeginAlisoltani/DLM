library(ggplot2)
library(forecast)
library(dlm)
library(Metrics)
library(dplyr)
library(lubridate)
library(scales)


##########################################################################

# 1_ Loading data
##########################################################################

###################### Electricity data
electricity <- read.csv("Electricity-2017-2022.csv", header = TRUE)

inf <- electricity$kW
###################### Temperature data
temp_data <- read.csv("temperature.CSV", header = TRUE)
#convert from f to c
temperature <- ((temp_data$temp) * 1.8) + 32
###################### Daylight data
daylight_data <- read.csv("daylight.CSV", header = TRUE)
daylight <- daylight_data$DaylightDuration
###################### Calendar data
calendar_data <- read.csv("calendar.CSV", header = TRUE)

vacation = matrix(as.numeric(calendar_data$vacation))
working = matrix(as.numeric(calendar_data$vacation))
lockdown = matrix(as.numeric(calendar_data$lockdown))
unlock = matrix(as.numeric(calendar_data$unlock))
crisis = matrix(as.numeric(calendar_data$crisis))

###################### Population data
population_data <- read.csv("population.csv", header = TRUE)
population <- population_data$population

##########################################################################

# 2- Modelling
##########################################################################
regressors <- cbind(war_crisis, lockdown, unlock,  temperature ,daylight)
colnames(regressors) <- c("war_crisis", "lockdown", "unlock","temperature","daylight")

fn <- function(parm, regressors) {
  
  mod <-  dlmModPoly(order = 1, dV = exp(parm[1])) + 
    dlmModSeas(frequency = 7) + 
    dlmModTrig(s = 365, q = 8) + 
    dlmModReg(regressors, addInt = FALSE, dW = exp(parm[2:9]))
  return(mod)
}

fit <- dlmMLE(inf, rep(0, 9), build = fn, hessian = TRUE, control = list(trace = TRUE, REPORT = 1), regressors = regressors)
conv <- fit$convergence  # zero for converged

loglik <- dlmLL(inf, fn(fit$par, regressors))
n.coef <- 1
r.aic <- -2 * loglik + 2 * n.coef
r.bic <- -2 * loglik + log(length(inf)) * n.coef
get_RMSE = function(original, predit){
  res = predit - original
  RSS = sum(res^2,na.rm=T)
  MSE = RSS/length(original)
  RMSE = sqrt(MSE)
  return(RMSE)
}

mod <- fn(fit$par, regressors)
obs.error.var <- V(mod)
state.error.var <- diag(W(mod))

filtered <- dlmFilter(inf, mod = mod)
smoothed <- dlmSmooth(filtered)
resids <- residuals(filtered, sd = FALSE)
mu <- dropFirst(smoothed$s[, 1])
r.rmse <- get_RMSE(inf,mu)

print(paste("log likelihood is:", loglik))
print(paste("AIC is:", r.aic))
print(paste("BIC is:", r.bic))
print(paste("RMSE is:", r.rmse))

##################################################################components
smoothed_error = dlmSvd2var(smoothed$U.S, smoothed$D.S)
smoothed$s = smoothed$s[-1,]

trend_elec = c()
err_trend_elec = c()

#yearly seasonality
yearly_elec = c()
err_yearly_elec = c()

#weekly seasonality
weekly_elec = c()
err_weekly_elec = c()

#temp
temp_elec = c()
err_temp_elec = c()

#week
week_elec = c()
err_week_elec = c()

#weekend
weekend_elec = c()
err_weekend_elec = c()

#vacation
vacation_elec = c()
err_vacation_elec = c()

#working
working_elec = c()
err_working_elec = c()

#lockdown
lockdown_elec = c()
err_lockdown_elec = c()

#unlockdown
unlock_elec = c()
err_unlock_elec = c()

#crisis
crisis_elec = c()
err_crisis_elec = c()

#daylight
daylight_elec = c()
err_daylight_elec = c()


#population
pop_elec = c()
err_pop_elec = c()


for (t in seq(1,length(inf))){
  XFF <- mod$FF
  XFF[mod$JFF != 0] <- mod$X[t, mod$JFF]
  #trend
  trend_elec = c(trend_elec,t(smoothed$s[t,1])%*%XFF[1])
  err_trend_elec = c(err_trend_elec, XFF[1] %*% smoothed_error[[t]][1,1] %*% t(XFF[1]))
  
  #weekly
  weekly_elec = c(weekly_elec,t(smoothed$s[t,2:7])%*%XFF[2:7])
  err_weekly_elec = c(err_weekly_elec, XFF[2] %*% smoothed_error[[t]][2,2] %*% t(XFF[2]))
  
  #yearly
  yearly_elec = c(yearly_elec,t(smoothed$s[t,8:23])%*%XFF[8:23])
  err_yearly_elec = c(err_yearly_elec, XFF[8] %*% smoothed_error[[t]][8,8] %*% t(XFF[8]))
  
  #exp: crisis
  crisis_elec = c(crisis_elec,t(smoothed$s[t,24])%*%XFF[24])
  err_crisis_elec = c(err_crisis_elec, XFF[24] %*% smoothed_error[[t]][24,24] %*% t(XFF[24]))
  
  #exp: lockdown
  lockdown_elec = c(lockdown_elec,t(smoothed$s[t,25])%*%XFF[25])
  err_lockdown_elec = c(err_lockdown_elec, XFF[25] %*% smoothed_error[[t]][25,25] %*% t(XFF[25]))
  
  #exp: unlockdown
  unlock_elec = c(unlock_elec,t(smoothed$s[t,26])%*%XFF[26])
  err_unlock_elec = c(err_unlock_elec, XFF[26] %*% smoothed_error[[t]][26,26] %*% t(XFF[26]))
  
  
  #exp: temp
  temp_elec = c(temp_elec,t(smoothed$s[t,27])%*%XFF[27])
  err_temp_elec = c(err_temp_elec,XFF[27] %*% smoothed_error[[t]][27,27] %*% t(XFF[27]))
  
  #exp: daylight
  daylight_elec = c(daylight_elec,t(smoothed$s[t,28])%*%XFF[28])
  err_daylight_elec = c(err_daylight_elec, XFF[28] %*% smoothed_error[[t]][28,28] %*% t(XFF[28]))
  
  #exp: vacation
  vacation_elec = c(vacation_elec,t(smoothed$s[t,29])%*%XFF[29])
  err_vacation_elec = c(err_vacation_elec, XFF[29] %*% smoothed_error[[t]][29,29] %*% t(XFF[29]))
  
  #exp: working
  working_elec = c(working_elec,t(smoothed$s[t,30])%*%XFF[30])
  err_working_elec = c(err_working_elec, XFF[30] %*% smoothed_error[[t]][30,30] %*% t(XFF[30]))
  
  #exp: pop
  pop_elec = c(pop_elec,t(smoothed$s[t,31])%*%XFF[31])
  err_pop_elec = c(err_pop_elec, XFF[31] %*% smoothed_error[[t]][31,31] %*% t(XFF[31]))
}
