################################################################################
# 1. RUNNING BASELINE SCENARIO #################################################
################################################################################
st<-Sys.time()
# source("age stratification and seasonality decompose.R")
source("age stratification and seasonality decompose.R")
source("model.R")
load("fit.RData")

# fitted parameters
parms["h.dhfA"] <- fit$par[1]
parms["h.dhfB"] <- fit$par[2]

parms["R0A"] <- fit$par[3]
parms["R0B"] <- fit$par[4]

parms["epsilonA"] <- fit$par[5]
parms["epsilonB"] <- fit$par[6]

parms["lagA"] <- fit$par[7]
parms["lagB"] <- fit$par[8]

# time steps
time_start <- 1
time_stop <- 360*20
deltat <- 1
tps <- seq(time_start , time_stop , by = deltat)

# vaccination campaign (at coverage = 0%, i.e., baseline scenario)
time_campaign <- 360*1
vac_cov <- 0
time_start_interv <- 360*2

# running the baseline model
parms["Sigma"] <- 0
parms["Pi"] <- 0

model_output0 <- deSolve::lsoda(y = istate.vaccine0, times = tps, func = model.vaccine0, parms = parms)

dhf(model_output0, period = 120, month = T, start = 360*2)

vac.df <- vacc(target_cov = 0.88, spec = 0.934, sens = 0.952, baseline = model_output0,
               screening = T, qdenga = 1, dengvaxia= 0)

################################################################################
# 2. SCENARIO ANALYSIS #########################################################
################################################################################

# fitted parameters
{
  parms["h.dhfA"] <- fit$par[1]
  parms["h.dhfB"] <- fit$par[2]
  
  parms["R0A"] <- fit$par[3]
  parms["R0B"] <- fit$par[4]
  
  parms["epsilonA"] <- fit$par[5]
  parms["epsilonB"] <- fit$par[6]
  
  parms["lagA"] <- fit$par[7]
  parms["lagB"] <- fit$par[8]
}

# timesteps
screnario.vec <- c("Scenario 1","Scenario 2","Scenario 3","Scenario 4","Scenario 5", "Baseline")
a <- length(screnario.vec)

model_output.list <- list()
dhf_output.list <- list()

for(i in 1:6){
  # time steps
  time_start <- 1
  time_stop <- 360*20
  deltat <- 1
  tps <- seq(time_start , time_stop , by = deltat)
  
  ## vaccine allocation per month during time campaign
  # assigning vaccine allocation, Sigma for Qdenga and Pi for Dengvaxia
  parms["Sigma"] <- c(0, 1, 1, 0.5, 0.5, 0)[i]
  parms["Pi"] <- c(1, 0, 0, 0.5, 0.5, 0)[i]
  
  screening <- c(T,F,T,F,T,F)[i]
  
  vac.df <- vacc(target_cov = 0.88, screening = screening, baseline = model_output0,
                 spec = 0.934, sens = 0.952,
                 qdenga = parms["Sigma"], dengvaxia= parms["Pi"])
  
  Qdenga.monthly.mat1 <- vac.df$`Qdenga with screening`
  Qdenga.monthly.mat2 <- vac.df$`Qdenga withiout screening`
  Dengvaxia.monthly.mat1 <- vac.df$`Dengvaxia with screening`
  
  # switch for qdenga with screening
  switch1 <- c(F,F,T,F,T,F)[i]       
  # switch for qdenga without screening
  switch2 <- c(F,T,F,T,F,F)[i]       
  # switch for dengvaxia with screening
  switch3 <- c(T,F,F,T,T,F)[i]       
  
  # running the screnario
  res0 <- deSolve::lsoda(y = istate.vaccine1, times = tps, func = model.vaccine1, parms = parms)
  
  model_output.list[[i]] <- res0
  
  df <- dhf(res0, 36, start = 1080)
  
  dhf_output.list[[i]] <- df
  
  print(i)
}

# saving model output
save(model_output.list, file = "model_output.list.Rdata")

# creating empty dataframe to store the calculation output
df.out.scenarioA <- data.frame(main = rep(0, 5),
                               scenario = 1:5,
                               sub_population = c(rep("A", 5)))

df.out.scenarioB <- data.frame(main = rep(0, 5),
                               scenario = 1:5,
                               sub_population = c(rep("B", 5)))

df.out.scenariototal <- data.frame(main = rep(0, 5),
                                   scenario = 1:5,
                                   sub_population = c(rep("total", 5)))

for(i in 1:5){
  # i = 1 for DENV4; 2 for DENV1; 3 for DENV2; 4 for DENV3.
  b <- colSums(dhf(model_output0, 36, start = 1080))
  c <- colSums(dhf(model_output.list[[i]], 36, start = 1080))
  
  d1 <- (b[1]-c[1])/b[1]*100
  d2 <- (b[2]-c[2])/b[2]*100
  d3 <- (sum(b) - sum(c))/sum(b)*100
  
  df.out.scenarioA[i,1] <- round(d1, 0)
  df.out.scenarioB[i,1] <- round(d2, 0)
  df.out.scenariototal[i,1] <- round(d3,0)
  
  # print(paste0("scenario: ",j, ", serotype: ",i, " = ", round(d,0)))
}

df.out.scenario <- cbind(df.out.scenarioA, df.out.scenarioB, df.out.scenariototal)

df.out.scenario

Sys.time() - st
################################################################################
# 3. SENSITIVITY ANALYSIS ######################################################
################################################################################

st<-Sys.time()

# fitted parameters
{
  parms["h.dhfA"] <- fit$par[1]
  parms["h.dhfB"] <- fit$par[2]
  
  parms["R0A"] <- fit$par[3]
  parms["R0B"] <- fit$par[4]
  
  parms["epsilonA"] <- fit$par[5]
  parms["epsilonB"] <- fit$par[6]
  
  parms["lagA"] <- fit$par[7]
  parms["lagB"] <- fit$par[8]
}

# timesteps
scenario.vec <- c("DENV 4","DENV 1","DENV 2","DENV 3")
a <- length(scenario.vec)

model_output.list.all <- list()
model_output.list2 <- list()
dhf_output.list.all <- list()
dhf_output.list2 <- list()

for(j in 1:5){
  # time steps
  time_start <- 1
  time_stop <- 360*20
  deltat <- 1
  tps <- seq(time_start , time_stop , by = deltat)
  
  ## vaccine allocation per month during time campaign
  # assigning vaccine allocation, Sigma for Qdenga and Pi for Dengvaxia
  parms["Sigma"] <- c(0, 1, 1, 0.5, 0.5)[j]
  parms["Pi"] <- c(1, 0, 0, 0.5, 0.5)[j]
  
  screening <- c(T,F,T,F,T)[j]
  
  for(i in 1:a){
    # vaccine parameters
    time_campaign <- 360*1
    vac_cov <- 0
    time_start_interv <- c(1321, 2234, 3084, 3998)[i]
    
    vac.df <- vacc(target_cov = 0.88, screening = screening, baseline = model_output0,
                   spec = 0.934, sens = 0.952,
                   qdenga = parms["Sigma"], dengvaxia= parms["Pi"])
    
    # switch for qdenga with screening
    Qdenga.monthly.mat1 <- vac.df$`Qdenga with screening`
    switch1 <- c(F,F,T,F,T)[j]       
    # switch for qdenga without screening
    Qdenga.monthly.mat2 <- vac.df$`Qdenga withiout screening`
    switch2 <- c(F,T,F,T,F)[j]       
    # switch for dengvaxia with screening
    Dengvaxia.monthly.mat1 <- vac.df$`Dengvaxia with screening`
    switch3 <- c(T,F,F,T,T)[j]     
    
    # running the screnario
    res0 <- deSolve::lsoda(y = istate.vaccine1, times = tps, func = model.vaccine1, parms = parms)
    model_output.list2[[paste0(scenario.vec[i])]] <- res0
    
    df <- dhf(model_output.list2[[i]], 36, start = time_start_interv)
    dhf_output.list2[[paste0(scenario.vec[i])]] <- df
    
    print(i)
  }
  
  model_output.list.all[[j]] <- model_output.list2
  dhf_output.list.all[[j]] <- dhf_output.list2
  print(j)
}

# saving model output
save(model_output.list.all, file = "model_output.list.all.Rdata")

# creating empty dataframe to store the calculation output
df.out.sensitivityA <- data.frame(timing_1 = rep(0, 5),
                                  timing_2 = rep(0, 5),
                                  timing_3 = rep(0, 5),
                                  timing_4 = rep(0, 5),
                                  scenario = 1:5,
                                  sub_population = c(rep("A", 5)))

df.out.sensitivityB <- data.frame(timing_1 = rep(0, 5),
                                  timing_2 = rep(0, 5),
                                  timing_3 = rep(0, 5),
                                  timing_4 = rep(0, 5),
                                  scenario = 1:5,
                                  sub_population = c(rep("B", 5)))

df.out.sensitivitytotal <- data.frame(timing_1 = rep(0, 5),
                                      timing_2 = rep(0, 5),
                                      timing_3 = rep(0, 5),
                                      timing_4 = rep(0, 5),
                                      scenario = 1:5,
                                      sub_population = c(rep("total", 5)))

# TOTAL DHF CASES
for(j in 1:5){
  for(i in 1:4){
    # i = 1 for DENV4; 2 for DENV1; 3 for DENV2; 4 for DENV3.
    b <- 1*colSums(dhf(model_output0, 36, 
                       start = c(1321, 2234, 3084, 3998)[i]))
    c <- 1*colSums(dhf(model_output.list.all[[j]][[i]], 36, 
                       start = c(1321, 2234, 3084, 3998)[i]))
    d1 <- (b[1]-c[1])/b[1]*100
    
    d2 <- (b[2]-c[2])/b[2]*100
    
    d3 <- (sum(b) - sum(c))/sum(b)*100
    
    df.out.sensitivityA[j,i] <- round(d1, 0)
    df.out.sensitivityB[j,i] <- round(d2, 0)
    df.out.sensitivitytotal[j,i] <- round(d3, 0)
  }
}

df.out.sensitivity <- cbind(df.out.sensitivityA,
                            df.out.sensitivityB,
                            df.out.sensitivitytotal)

df.out.sensitivity

Sys.time() - st

################################################################################
# 4. PLOTTING ##################################################################
################################################################################
## define function for tibble mutate ####
df.mutate <- function(x){
  y <- as_tibble(as.data.frame(x)) %>% 
    mutate(PA = SA+S1A+S2A+S3A+S4A+
             I1A+I2A+I3A+I4A+I12A+I13A+I14A+I21A+I23A+I24A+I31A+I32A+I34A+I41A+I42A+I43A+
             R1A+R2A+R3A+R4A+R12A+R13A+R14A+R21A+R23A+R24A+R31A+R32A+R34A+R41A+R42A+R43A+
             SAQ+S1AQ+S2AQ+S3AQ+S4AQ+
             I1AQ+I2AQ+I3AQ+I4AQ+I12AQ+I13AQ+I14AQ+I21AQ+I23AQ+I24AQ+I31AQ+I32AQ+I34AQ+I41AQ+I42AQ+I43AQ+
             R1AQ+R2AQ+R3AQ+R4AQ+R12AQ+R13AQ+R14AQ+R21AQ+R23AQ+R24AQ+R31AQ+R32AQ+R34AQ+R41AQ+R42AQ+R43AQ+
             SAD+S1AD+S2AD+S3AD+S4AD+
             I1AD+I2AD+I3AD+I4AD+I12AD+I13AD+I14AD+I21AD+I23AD+I24AD+I31AD+I32AD+I34AD+I41AD+I42AD+I43AD+
             R1AD+R2AD+R3AD+R4AD+R12AD+R13AD+R14AD+R21AD+R23AD+R24AD+R31AD+R32AD+R34AD+R41AD+R42AD+R43AD,
           
           PB = SB+S1B+S2B+S3B+S4B+
             I1B+I2B+I3B+I4B+I12B+I13B+I14B+I21B+I23B+I24B+I31B+I32B+I34B+I41B+I42B+I43B+
             R1B+R2B+R3B+R4B+R12B+R13B+R14B+R21B+R23B+R24B+R31B+R32B+R34B+R41B+R42B+R43B+
             SBQ+S1BQ+S2BQ+S3BQ+S4BQ+
             I1BQ+I2BQ+I3BQ+I4BQ+I12BQ+I13BQ+I14BQ+I21BQ+I23BQ+I24BQ+I31BQ+I32BQ+I34BQ+I41BQ+I42BQ+I43BQ+
             R1BQ+R2BQ+R3BQ+R4BQ+R12BQ+R13BQ+R14BQ+R21BQ+R23BQ+R24BQ+R31BQ+R32BQ+R34BQ+R41BQ+R42BQ+R43BQ+
             SBD+S1BD+S2BD+S3BD+S4BD+
             I1BD+I2BD+I3BD+I4BD+I12BD+I13BD+I14BD+I21BD+I23BD+I24BD+I31BD+I32BD+I34BD+I41BD+I42BD+I43BD+
             R1BD+R2BD+R3BD+R4BD+R12BD+R13BD+R14BD+R21BD+R23BD+R24BD+R31BD+R32BD+R34BD+R41BD+R42BD+R43BD,
           
           P = PA + PB,
           
           DENV1 = I1A+I21A+I31A+I41A + I1B+I21B+I31B+I41B + 
             I1AQ+I21AQ+I31AQ+I41AQ + I1BQ+I21BQ+I31BQ+I41BQ +
             I1AD+I21AD+I31AD+I41AD + I1BD+I21BD+I31BD+I41BD,
           
           DENV2 = I2A+I12A+I32A+I42A + I2B+I12B+I32B+I42B +
             I2AQ+I12AQ+I32AQ+I42AQ + I2BQ+I12BQ+I32BQ+I42BQ +
             I2AD+I12AD+I32AD+I42AD + I2BD+I12BD+I32BD+I42BD,
           
           DENV3 = I3A+I13A+I23A+I43A + I3B+I13B+I23B+I43B +
             I3AQ+I13AQ+I23AQ+I43AQ + I3BQ+I13BQ+I23BQ+I43BQ +
             I3AD+I13AD+I23AD+I43AD + I3BD+I13BD+I23BD+I43BD,
           
           DENV4 = I4A+I14A+I24A+I34A + I4B+I14B+I24B+I34B +
             I4AQ+I14AQ+I24AQ+I34AQ + I4BQ+I14BQ+I24BQ+I34BQ +
             I4AD+I14AD+I24AD+I34AD + I4BD+I14BD+I24BD+I34BD,
           
           SeronegativeA = SA/P,
           SeronegativeAB = (SA+SB)/P,
           
           SeropositiveA = 1 - SA/P,
           SeropositiveAB = 1 - (SA+SB)/P,
           
           DENV1_Proportion = (DENV1)/(DENV1+DENV2+DENV3+DENV4),
           DENV2_Proportion = (DENV2)/(DENV1+DENV2+DENV3+DENV4),
           DENV3_Proportion = (DENV3)/(DENV1+DENV2+DENV3+DENV4),
           DENV4_Proportion = (DENV4)/(DENV1+DENV2+DENV3+DENV4)
           
    ) %>% 
    pivot_longer(names_to = "variable", cols =! 1)
  
  return(y)
}
DENV.DYNAMICS4 <- function(x){
  df.mutate(x) %>%
    filter(variable %in% c("DENV1_Proportion", "DENV2_Proportion", "DENV3_Proportion", "DENV4_Proportion")) %>% 
    ggplot() +
    aes(x = time, y = value, colour = variable) +
    geom_rect(aes(xmin = 1321, 
                  xmax = 1321 + 360, 
                  ymin = -Inf, 
                  ymax = Inf),
              fill = "#b2e4f0", color = NA, alpha = 0.02) +
    geom_rect(aes(xmin = 1321 + 360, 
                  xmax = 1321 + 360 + 720, 
                  ymin = -Inf, 
                  ymax = Inf),
              fill = "#FFDFD3", color = NA, alpha = 0.02) +
    geom_vline(xintercept = 1321, color = "black", linetype = "dashed") +
    geom_vline(xintercept = 1321 + 360 + 720, color = "black", linetype = "dashed") +
    geom_line() +
    #    geom_point(size = 0.5) +
    scale_color_hue(direction = 1) +
    theme_minimal() + 
    theme(legend.position = "bottom") +
    labs(title = paste0("Vaccination timeline for when DENV-4 is the most dominant"),
         y = ("Proportion"))
}
DENV.DYNAMICS1 <- function(x){
  df.mutate(x) %>%
    filter(variable %in% c("DENV1_Proportion", "DENV2_Proportion", "DENV3_Proportion", "DENV4_Proportion")) %>% 
    ggplot() +
    aes(x = time, y = value, colour = variable) +
    geom_rect(aes(xmin = 2234, 
                  xmax = 2234 + 360, 
                  ymin = -Inf, 
                  ymax = Inf),
              fill = "#b2e4f0", color = NA, alpha = 0.02) +
    geom_rect(aes(xmin = 2234 + 360, 
                  xmax = 2234 + 360 + 720, 
                  ymin = -Inf, 
                  ymax = Inf),
              fill = "#FFDFD3", color = NA, alpha = 0.02) +
    geom_vline(xintercept = 2234, color = "black", linetype = "dashed") +
    geom_vline(xintercept = 2234 + 360 + 720, color = "black", linetype = "dashed") +
    geom_line() +
    #    geom_point(size = 0.5) +
    scale_color_hue(direction = 1) +
    theme_minimal() + 
    theme(legend.position = "bottom") +
    labs(title = paste0("Vaccination timeline for when DENV-1 is the most dominant"),
         y = ("Proportion"))
}
DENV.DYNAMICS2 <- function(x){
  df.mutate(x) %>%
    filter(variable %in% c("DENV1_Proportion", "DENV2_Proportion", "DENV3_Proportion", "DENV4_Proportion")) %>% 
    ggplot() +
    aes(x = time, y = value, colour = variable) +
    geom_rect(aes(xmin = 3084, 
                  xmax = 3084 + 360, 
                  ymin = -Inf, 
                  ymax = Inf),
              fill = "#b2e4f0", color = NA, alpha = 0.02) +
    geom_rect(aes(xmin = 3084 + 360, 
                  xmax = 3084 + 360 + 720, 
                  ymin = -Inf, 
                  ymax = Inf),
              fill = "#FFDFD3", color = NA, alpha = 0.02) +
    geom_vline(xintercept = 3084, color = "black", linetype = "dashed") +
    geom_vline(xintercept = 3084 + 360 + 720, color = "black", linetype = "dashed") +
    geom_line() +
    #    geom_point(size = 0.5) +
    scale_color_hue(direction = 1) +
    theme_minimal() + 
    theme(legend.position = "bottom") +
    labs(title = paste0("Vaccination timeline for when DENV-2 is the most dominant"),
         y = ("Proportion"))
}
DENV.DYNAMICS3 <- function(x){  
  df.mutate(x) %>%
    filter(variable %in% c("DENV1_Proportion", "DENV2_Proportion", "DENV3_Proportion", "DENV4_Proportion")) %>% 
    ggplot() +
    aes(x = time, y = value, colour = variable) +
    geom_rect(aes(xmin = 3998, 
                  xmax = 3998 + 360, 
                  ymin = -Inf, 
                  ymax = Inf),
              fill = "#b2e4f0", color = NA, alpha = 0.02) +
    geom_rect(aes(xmin = 3998 + 360, 
                  xmax = 3998 + 360 + 720, 
                  ymin = -Inf, 
                  ymax = Inf),
              fill = "#FFDFD3", color = NA, alpha = 0.02) +
    geom_vline(xintercept = 3998, color = "black", linetype = "dashed") +
    geom_vline(xintercept = 3998 + 360 + 720, color = "black", linetype = "dashed") +
    geom_line() +
    #    geom_point(size = 0.5) +
    scale_color_hue(direction = 1) +
    theme_minimal() + 
    theme(legend.position = "bottom") +
    labs(title = paste0("Vaccination timeline for when DENV-3 is the most dominant"),
         y = ("Proportion"))
}

baseline4 <- DENV.DYNAMICS4(model_output0)
baseline3 <- DENV.DYNAMICS3(model_output0)
baseline2 <- DENV.DYNAMICS2(model_output0)
baseline1 <- DENV.DYNAMICS1(model_output0)

scenario1.4 <- DENV.DYNAMICS4(model_output.list.all[[1]]$`DENV 4`)
scenario1.3 <- DENV.DYNAMICS3(model_output.list.all[[1]]$`DENV 3`)
scenario1.2 <- DENV.DYNAMICS2(model_output.list.all[[1]]$`DENV 2`)
scenario1.1 <- DENV.DYNAMICS1(model_output.list.all[[1]]$`DENV 1`)

scenario2.4 <- DENV.DYNAMICS4(model_output.list.all[[2]]$`DENV 4`)
scenario2.3 <- DENV.DYNAMICS3(model_output.list.all[[2]]$`DENV 3`)
scenario2.2 <- DENV.DYNAMICS2(model_output.list.all[[2]]$`DENV 2`)
scenario2.1 <- DENV.DYNAMICS1(model_output.list.all[[2]]$`DENV 1`)

scenario3.4 <- DENV.DYNAMICS4(model_output.list.all[[3]]$`DENV 4`)
scenario3.3 <- DENV.DYNAMICS3(model_output.list.all[[3]]$`DENV 3`)
scenario3.2 <- DENV.DYNAMICS2(model_output.list.all[[3]]$`DENV 2`)
scenario3.1 <- DENV.DYNAMICS1(model_output.list.all[[3]]$`DENV 1`)

scenario4.4 <- DENV.DYNAMICS4(model_output.list.all[[4]]$`DENV 4`)
scenario4.3 <- DENV.DYNAMICS3(model_output.list.all[[4]]$`DENV 3`)
scenario4.2 <- DENV.DYNAMICS2(model_output.list.all[[4]]$`DENV 2`)
scenario4.1 <- DENV.DYNAMICS1(model_output.list.all[[4]]$`DENV 1`)

scenario5.4 <- DENV.DYNAMICS4(model_output.list.all[[5]]$`DENV 4`)
scenario5.3 <- DENV.DYNAMICS3(model_output.list.all[[5]]$`DENV 3`)
scenario5.2 <- DENV.DYNAMICS2(model_output.list.all[[5]]$`DENV 2`)
scenario5.1 <- DENV.DYNAMICS1(model_output.list.all[[5]]$`DENV 1`)


z0 <- gridExtra::grid.arrange(baseline1, baseline2,
                              baseline3, baseline4,
                              ncol = 2)

z1 <- gridExtra::grid.arrange(scenario1.1, scenario1.2,
                              scenario1.3, scenario1.4,
                              ncol = 2)

z2 <- gridExtra::grid.arrange(scenario2.1, scenario2.2,
                              scenario2.3, scenario2.4,
                              ncol = 2)

z3 <- gridExtra::grid.arrange(scenario3.1, scenario3.2,
                              scenario3.3, scenario3.4,
                              ncol = 2)

z4 <- gridExtra::grid.arrange(scenario4.1, scenario4.2,
                              scenario4.3, scenario4.4,
                              ncol = 2)

z5 <- gridExtra::grid.arrange(scenario5.1, scenario5.2,
                              scenario5.3, scenario5.4,
                              ncol = 2)
