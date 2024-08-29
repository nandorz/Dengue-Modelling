st<-Sys.time()

# source("age stratification and seasonality decompose.R")
source("age stratification and seasonality decompose.R")
source("model.R")

### MODEL FITTING ####
# assigning timesteps
time_start <- 1
time_stop <- 360*50
deltat <- 1
tps <- seq(time_start , time_stop , by = deltat)

parms["R0A"] <- 7
parms["R0B"] <- 1

# fitted parameters
par1.st <- parms["h.dhfA"] <- 0.05
par2.st <- parms["h.dhfB"] <- 0.22

par3.st <- parms["epsilonA"] <- 7
par4.st <- parms["epsilonB"] <- 1

par5.st <- parms["epsilonA"] <- 0.18
par6.st <- parms["epsilonB"] <- 0.18

par7.st <- parms["lagA"] <- 8.505
par8.st <- parms["lagB"] <- 8.504

# fitting the model using optim()
period = 72
par.set0 <- c(par1.st, par2.st, par3.st, par4.st, par5.st, par6.st, par7.st, par8.st)
calcNLL(par.set0, month = T, day = F)

low.bound <- c(0.01, 0.01, 1, 1, 0.05, 0.05, 1, 1)
upp.bound <- c(0.3, 0.3, 10, 10, 0.35, 0.35, 12, 12)

month = T
day = F

fit <- optim(c(par.set0),
              fn = calcNLL,
              lower = low.bound,
              upper = upp.bound,
              method="L-BFGS-B",
              hessian = T)

save(fit, file = "fit.Rdata")

fit$par

parms["h.dhfA"] <- fit$par[1]
parms["h.dhfB"] <- fit$par[2]

parms["R0A"] <- fit$par[3]
parms["R0B"] <- fit$par[4]

parms["epsilonA"] <- fit$par[3]
parms["epsilonB"] <- fit$par[4]

parms["lagA"] <- fit$par[5]
parms["lagB"] <- fit$par[6]

model_output.start <- deSolve::lsoda(y = istate0, times = tps, func = model0, parms = parms)

dhf_sim.df0 <- dhf(model_output.start, period = 72, month = T, start = 360*40)

# averaging reported cases by month
avg.A <- c()
se.A <- c()

for(i in 1:12){
  avg.A[i] <- mean(c(casesA[i], casesA[i+12], casesA[i+24], casesA[i+36], casesA[i+48]))
  se.A[i] <- sd(c(casesA[i], casesA[i+12], casesA[i+24], casesA[i+36], casesA[i+48]))
}

avg.B <- c()
for(i in 1:12){
  avg.B[i] <- mean(c(casesB[i], casesB[i+12], casesB[i+24], casesB[i+36], casesB[i+48]))
}

# Create a sequence of dates from January 2018 to December 2022
start_date <- as.Date("2017-01-01")
end_date <- as.Date("2022-12-01")
all_months <- seq(from = start_date, to = end_date, by = "month")

# Convert the dates to character strings in the "yyyy-mm" format
all_months_string <- format(all_months, "%Y-%m")

par(mfrow = c(2,1))
plot(y = casesA*100000/npop.jak0.14, 
     x = all_months, 
     col = "gray30", 
     ylab = "People per 100000", xlab = "Month",
     main = "Model Fitting to DHF Incidence of Age Group A: 0-14")
lines(y = casesA*100000/npop.jak0.14, 
      x = all_months,
      col = "gray30")
lines(y = rep(avg.A, 6)*100000/npop.jak0.14, 
      x = all_months,
      type = "b", 
      col= "red")
lines(y = dhf_sim.df0$DHF.A*100000/npop.jak0.14,
      x = all_months,
      col = "blue")

plot(y = casesB*100000/npop.jak15.75, x = all_months, 
     col = "gray30", 
     ylab = "People per 100000", xlab = "Month",
     main = "Model Fitting to Incidence of Age Group B: 15-75+")
lines(y = casesB*100000/npop.jak15.75, 
      x = all_months, 
      col= "gray30")
lines(y = rep(avg.B*100000/npop.jak15.75, 6), 
      x = all_months,
      type = "b", col= "red")
lines(y = dhf_sim.df0$DHF.B*100000/npop.jak15.75, 
      x = all_months,
      col = "blue")

Sys.time() - st