# https://data.cityofnewyork.us/resource/erm2-nwe9.json

## Install the required package with:
## install.packages("RSocrata")

library("RSocrata")
library('dplyr')
library("forecast")
library('magrittr')
library("tseries")
library('DT')
# data url: https://data.cityofnewyork.us/Social-Services/311-Service-Requests-from-2010-to-Present/erm2-nwe9
# api docs: https://dev.socrata.com/docs/filtering.html

df <- read.socrata(url="https://data.cityofnewyork.us/resource/erm2-nwe9.json?$where=created_date between '2019-01-01T00:00:00.000' and '2021-04-01T00:00:00.000' AND complaint_type='Noise - Residential' AND city='NEW YORK' AND borough='MANHATTAN'")
head(df)

df$Date <- as.Date(df$created_date)
#noise_complaints <- df %>% group_by(Date) %>% tally
noise_complaints <- df %>% group_by(week = format(Date, '%Y-%U')) %>% tally %>% as.data.frame()
noise_complaints$week <- as.Date(paste0(noise_complaints$week, '-1'),'%Y-%U-%u')
noise_complaints$log_n <- log(noise_complaints$n)
noise_complaints$n_diff <- c(NA, diff(noise_complaints$n))
noise_complaints$log_n_diff <- c(NA, diff(noise_complaints$log_n))
noise_complaints$n_diff2 <- c(NA, diff(noise_complaints$n_diff))


# four plots
par(mfrow=c(2,2))
# series, diff series
plot(noise_complaints$week, noise_complaints$n, type = 'l')
plot(noise_complaints$week, noise_complaints$n_diff, type = 'l')
# log series, diff log series
plot(noise_complaints$week, noise_complaints$log_n, type = 'l')
plot(noise_complaints$week, noise_complaints$log_n_diff, type = 'l')

# looks like we don't need to log the series
# could be d=1 since it seems more mean reverting

par(mfrow=c(3,2))
# n
Acf(noise_complaints$n) # hangs
Pacf(noise_complaints$n)
# diff n
Acf(noise_complaints$n_diff) 
Pacf(noise_complaints$n_diff)
# diff2 n - overdifferenced
Acf(noise_complaints$n_diff2) 
Pacf(noise_complaints$n_diff2)

# modeling 

possible_p<-0:2
possible_q<-0:2

aicc_results <- data.frame(p=character(),
                           d=character(),
                           q=character(), 
                           constant=character(), 
                           aicc=numeric()) 

# model selection among candidates
i=0
for(p in possible_p){
        for(q in possible_q){
                # constant
                i<-i+1
                #print(paste0('Running: ARIMA(', p, ', 0, ', q, ' with constant.)'))
                aicc_results[i, 'p']<-p
                aicc_results[i, 'd']<-1
                aicc_results[i, 'q']<-q
                aicc_results[i, 'constant']<-TRUE
                aicc_results[i, 'aicc']<- Arima(noise_complaints$n_diff, c(p, 0, q), include.constant=T)$aicc
                # no constant
                i<-i+1
                #print(paste0('Running: ARIMA(', p, ', 1, ', q, ' without constant.)'))
                aicc_results[i, 'p']<-p
                aicc_results[i, 'd']<-1
                aicc_results[i, 'q']<-q
                aicc_results[i, 'constant']<-FALSE
                aicc_results[i, 'aicc']<-Arima(noise_complaints$n_diff, c(p, 0, q), include.constant=F)$aicc
                
        }
}

aicc_results$aicc <- round(aicc_results$aicc, 2)
datatable(aicc_results[order(aicc_results$aicc),], 
          options = list(pageLength = nrow(aicc_results)))

# select the arima 1, 1, 1 with AICc=1565.94

# model selected
fit_select <- Arima(noise_complaints$n_diff, c(1,0,1), include.constant=T, method="ML")
summary(fit_select)
# model for forecasting
fit.mean<-Arima(noise_complaints$n, c(1, 1, 1), include.constant=T, method="ML")
summary(fit)

resid <- residuals(fit.mean)
#tail(resid, n=10)
f <- fitted.values(fit.mean)
#tail(f, n=10)

fcasts<-forecast(fit.mean, h = 1) %>% as.data.frame()
fcasts[, c(1, 4, 5)]

plot(noise_complaints$week, resid, type="l",
     xlab="Date", ylab="Residuals Fitted")
Acf(resid)
Pacf(resid)

plot(noise_complaints$week, resid^2, type="l",
     xlab="Date", ylab="Squared Residuals Fitted")

q <- 0:10
loglik <- rep(NA, length(q))
N <- length(resid)

for (i in 1:length(q)) {
        if (q[i] == 0) {
                
                # Calculate the log likelihood for the ARCH(0) model by hand. 
                #See the handout on Estimation and Automatic Selection of ARCH models.
                
                loglik[i] <- -0.5 * N * (1 + log(2 * pi * mean(resid^2)))
                
        } else {
                fit <- garch(resid, c(0,q[i]), trace=FALSE)
                loglik[i] <- logLik(fit)
        }
}


k <- q + 1
aicc <- -2 * loglik  + 2 * k * N / (N - k - 1)

datatable(data.frame(q, loglik, aicc))


fit <- garch(resid, c(1,1), trace=FALSE)
loglik <- logLik(fit)
k <- 2
aicc <- -2 * loglik  + 2 * k * N / (N - k - 1)
datatable(data.frame(loglik, aicc))


fit.var <- garch(resid, c(1,1), trace=FALSE)
summary(fit.var)
logLik(fit.var)

f1 <- fcasts$`Point Forecast`
ht <- fit.var$fit[,1]^2
h1 <- 1.508e-07 + tail(ht, n=1) %>% as.numeric()
f1 + c(-1, 1) * 1.96 * sqrt(h1)

ht <- fit.var$fit[,1]^2
# tail(ht, n=10)
plot(noise_complaints$week, ht, type="l", col=4)

plot(noise_complaints$week, noise_complaints$n, type="l")
lines(noise_complaints$week, f + 1.96 * sqrt(ht), lty=2, col=2)
lines(noise_complaints$week, f - 1.96 * sqrt(ht), lty=2, col=2)

resid.arch <- fit.var$residuals
qqnorm(resid.arch)

sum(abs(resid.arch) > 0.5, na.rm=TRUE)
sum(abs(resid.arch) > 1.96, na.rm=TRUE)
sum(!is.na(resid.arch))
print(paste0("% miss: " , sum(abs(resid.arch) > 1.96, na.rm=TRUE)  / sum(!is.na(resid.arch))))