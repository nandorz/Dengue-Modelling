# loading libraries
library(pacman)
p_load(deSolve, tidyverse, gridExtra, readxl, ggplot2)

# Loading parameters ####
parameters <- read_excel("parameters.xlsx")
parameters <- as.data.frame(parameters)

parameters_value.vector <- c(parameters$value_day)  # taking the value of parameters from data frame into a vector
parameters_name.vector <- c(parameters$name)    # taking the name of parameters from data frame into a vector
names(parameters_value.vector) <- parameters_name.vector  # giving the parameters value with name

parms <- parameters_value.vector

# Initial conditions ####
# loading compartments and assigning their initial states: for non-vaccination model
compartments <- read_excel("initial conditions.xlsx", sheet = "model2")
compartments <- as.data.frame(compartments)

compartments_value.vector <- c(compartments$value)  # taking the value of parameters from data frame into a vector
compartments_name.vector <- c(compartments$name)    # taking the name of parameters from data frame into a vector
names(compartments_value.vector) <- compartments_name.vector  # giving the parameters value with name

istate <- compartments_value.vector
istate0 <- c(istate, 
             CumIncDHFA = 0, CumIncDHFB = 0,
             CumIncDF = 0, CumIncVCD = 0, CumIncAsymptomatic = 0, 
             CumIncHospDHFA = 0, CumIncHospDHFB = 0, deathA = 0, ddeathB = 0, 
             lambda1A = 0, lambda2A = 0, lambda3A = 0, lambda4A = 0,
             lambda12A = 0, lambda21A = 0, lambda31A = 0, lambda41A = 0)

# loading compartments and assigning their initial states: for vaccination model
compartments <- read_excel("initial conditions.xlsx",
                           sheet = "double vaccination")
compartments <- as.data.frame(compartments)

compartments_value.vector <- c(compartments$value)  # taking the value of parameters from data frame into a vector
compartments_name.vector <- c(compartments$name)    # taking the name of parameters from data frame into a vector
names(compartments_value.vector) <- compartments_name.vector  # giving the parameters value with name

# for baseline
istate <- compartments_value.vector
istate.vaccine0 <- c(istate, 
                     CumIncVCDA = 0, CumIncVCDB = 0, 
                     CumIncHospDHFA = 0, CumIncHospDHFB = 0, deathA = 0, deathB = 0,
                     PA_unvacc = 0, PA_vacc = 0, PA = 0, PB = 0, P=0,
                     lambda1AQ = 0, lambda2AQ = 0, lambda3AQ = 0, lambda4AQ = 0, 
                     lambda1AD = 0, lambda2AD = 0, lambda3AD = 0, lambda4AD = 0)

# for intervention scenario
istate <- compartments_value.vector
istate.vaccine1<- c(istate, 
                    CumIncVCDA = 0, CumIncVCDB = 0, 
                    CumIncHospDHFA = 0, CumIncHospDHFB = 0, deathA = 0, deathB = 0,
                    PA_unvacc = 0, PA_vacc = 0, PA = 0, PB = 0, P=0,
                    Cumvaccsn = 0, Cumvaccsp = 0,
                    lambda1AQ = 0, lambda2AQ = 0, lambda3AQ = 0, lambda4AQ = 0, 
                    lambda1AD = 0, lambda2AD = 0, lambda3AD = 0, lambda4AD = 0)

# 1. Dengue Vaccination Model ####
# model construction
model.vaccine0 <- function(t, state, parms) {
  with(as.list(c(state, parms)),{
    PA <- SA+S1A+S2A+S3A+S4A+I1A+I2A+I3A+I4A+I12A+I13A+I14A+I21A+I23A+I24A+I31A+I32A+I34A+I41A+I42A+I43A+R1A+R2A+R3A+R4A+R12A+R13A+R14A+R21A+R23A+R24A+R31A+R32A+R34A+R41A+R42A+R43A+SAQ+S1AQ+S2AQ+S3AQ+S4AQ+I1AQ+I2AQ+I3AQ+I4AQ+I12AQ+I13AQ+I14AQ+I21AQ+I23AQ+I24AQ+I31AQ+I32AQ+I34AQ+I41AQ+I42AQ+I43AQ+R1AQ+R2AQ+R3AQ+R4AQ+R12AQ+R13AQ+R14AQ+R21AQ+R23AQ+R24AQ+R31AQ+R32AQ+R34AQ+R41AQ+R42AQ+R43AQ+SAD+S1AD+S2AD+S3AD+S4AD+I1AD+I2AD+I3AD+I4AD+I12AD+I13AD+I14AD+I21AD+I23AD+I24AD+I31AD+I32AD+I34AD+I41AD+I42AD+I43AD+R1AD+R2AD+R3AD+R4AD+R12AD+R13AD+R14AD+R21AD+R23AD+R24AD+R31AD+R32AD+R34AD+R41AD+R42AD+R43AD
    PB <- SB+S1B+S2B+S3B+S4B+I1B+I2B+I3B+I4B+I12B+I13B+I14B+I21B+I23B+I24B+I31B+I32B+I34B+I41B+I42B+I43B+R1B+R2B+R3B+R4B+R12B+R13B+R14B+R21B+R23B+R24B+R31B+R32B+R34B+R41B+R42B+R43B+SBQ+S1BQ+S2BQ+S3BQ+S4BQ+I1BQ+I2BQ+I3BQ+I4BQ+I12BQ+I13BQ+I14BQ+I21BQ+I23BQ+I24BQ+I31BQ+I32BQ+I34BQ+I41BQ+I42BQ+I43BQ+R1BQ+R2BQ+R3BQ+R4BQ+R12BQ+R13BQ+R14BQ+R21BQ+R23BQ+R24BQ+R31BQ+R32BQ+R34BQ+R41BQ+R42BQ+R43BQ+SBD+S1BD+S2BD+S3BD+S4BD+I1BD+I2BD+I3BD+I4BD+I12BD+I13BD+I14BD+I21BD+I23BD+I24BD+I31BD+I32BD+I34BD+I41BD+I42BD+I43BD+R1BD+R2BD+R3BD+R4BD+R12BD+R13BD+R14BD+R21BD+R23BD+R24BD+R31BD+R32BD+R34BD+R41BD+R42BD+R43BD
    P <- PA + PB
    
    # unvaccinated population in age group A
    PA_unvacc <- SA+S1A+S2A+S3A+S4A+
      I1A+I2A+I3A+I4A+I12A+I13A+I14A+I21A+I23A+I24A+I31A+I32A+I34A+I41A+I42A+I43A+
      R1A+R2A+R3A+R4A+R12A+R13A+R14A+R21A+R23A+R24A+R31A+R32A+R34A+R41A+R42A+R43A
    
    PA_vacc <- SAQ+S1AQ+S2AQ+S3AQ+S4AQ+
      I1AQ+I2AQ+I3AQ+I4AQ+I12AQ+I13AQ+I14AQ+I21AQ+I23AQ+I24AQ+I31AQ+I32AQ+I34AQ+I41AQ+I42AQ+I43AQ+
      R1AQ+R2AQ+R3AQ+R4AQ+R12AQ+R13AQ+R14AQ+R21AQ+R23AQ+R24AQ+R31AQ+R32AQ+R34AQ+R41AQ+R42AQ+R43AQ+
      SAD+S1AD+S2AD+S3AD+S4AD+
      I1AD+I2AD+I3AD+I4AD+I12AD+I13AD+I14AD+I21AD+I23AD+I24AD+I31AD+I32AD+I34AD+I41AD+I42AD+I43AD+
      R1AD+R2AD+R3AD+R4AD+R12AD+R13AD+R14AD+R21AD+R23AD+R24AD+R31AD+R32AD+R34AD+R41AD+R42AD+R43AD
    
    ### SEASONALITY
    # accounting seasonality forcing, variability within a year, and inter-annual variability
    seasA <- 1 - epsilonA*cos(2*pi*t/(360) + lagA*30)
    seasB <- 1 - epsilonB*cos(2*pi*t/(360) + lagB*30)
    
    betaA <- R0A*(tau + mu_A)
    betaB <- R0B*(tau + mu_B)
    
    ### FORCE OF INFECTION ESTIMATION
    # infection importation for each group
    ma <- m*PA
    mb <- m*PB
    
    # FOI for age group A: 0-14 years old
    lambda1A <- betaA*seasA*(I1A+I21A+I31A+I41A + I1B+I21B+I31B+I41B + 
                               I1AQ+I21AQ+I31AQ+I41AQ + I1BQ+I21BQ+I31BQ+I41BQ + 
                               I1AD+I21AD+I31AD+I41AD + I1BD+I21BD+I31BD+I41BD + ma)/P
    lambda2A <- betaA*seasA*(I2A+I12A+I32A+I42A + I2B+I12B+I32B+I42B + 
                               I2AQ+I12AQ+I32AQ+I42AQ + I2BQ+I12BQ+I32BQ+I42BQ + 
                               I2AD+I12AD+I32AD+I42AD + I2BD+I12BD+I32BD+I42BD + ma)/P
    lambda3A <- betaA*seasA*(I3A+I13A+I23A+I43A + I3B+I13B+I23B+I43B + 
                               I3AQ+I13AQ+I23AQ+I43AQ + I3BQ+I13BQ+I23BQ+I43BQ + 
                               I3AD+I13AD+I23AD+I43AD + I3BD+I13BD+I23BD+I43BD + ma)/P
    lambda4A <- betaA*seasA*(I4A+I14A+I24A+I34A + I4B+I14B+I24B+I34B + 
                               I4AQ+I14AQ+I24AQ+I34AQ + I4BQ+I14BQ+I24BQ+I34BQ + 
                               I4AD+I14AD+I24AD+I34AD + I4BD+I14BD+I24BD+I34BD + ma)/P
    lambda12A <- lambda2A
    lambda13A <- lambda3A
    lambda14A <- lambda4A
    lambda21A <- lambda1A
    lambda23A <- lambda3A
    lambda24A <- lambda4A
    lambda31A <- lambda1A
    lambda32A <- lambda2A
    lambda34A <- lambda4A
    lambda41A <- lambda1A
    lambda42A <- lambda2A
    lambda43A <- lambda3A
    
    # FOI for age group B: 15+ years old
    lambda1B <- betaB*seasB*(I1A+I21A+I31A+I41A + I1B+I21B+I31B+I41B + 
                               I1AQ+I21AQ+I31AQ+I41AQ + I1BQ+I21BQ+I31BQ+I41BQ + 
                               I1AD+I21AD+I31AD+I41AD + I1BD+I21BD+I31BD+I41BD + mb)/P
    lambda2B <- betaB*seasB*(I2A+I12A+I32A+I42A + I2B+I12B+I32B+I42B + 
                               I2AQ+I12AQ+I32AQ+I42AQ + I2BQ+I12BQ+I32BQ+I42BQ + 
                               I2AD+I12AD+I32AD+I42AD + I2BD+I12BD+I32BD+I42BD + mb)/P
    lambda3B <- betaB*seasB*(I3A+I13A+I23A+I43A + I3B+I13B+I23B+I43B + 
                               I3AQ+I13AQ+I23AQ+I43AQ + I3BQ+I13BQ+I23BQ+I43BQ + 
                               I3AD+I13AD+I23AD+I43AD + I3BD+I13BD+I23BD+I43BD + mb)/P
    lambda4B <- betaB*seasB*(I4A+I14A+I24A+I34A + I4B+I14B+I24B+I34B + 
                               I4AQ+I14AQ+I24AQ+I34AQ + I4BQ+I14BQ+I24BQ+I34BQ + 
                               I4AD+I14AD+I24AD+I34AD + I4BD+I14BD+I24BD+I34BD + mb)/P
    lambda12B <- lambda2B
    lambda13B <- lambda3B
    lambda14B <- lambda4B
    lambda21B <- lambda1B
    lambda23B <- lambda3B
    lambda24B <- lambda4B
    lambda31B <- lambda1B
    lambda32B <- lambda2B
    lambda34B <- lambda4B
    lambda41B <- lambda1B
    lambda42B <- lambda2B
    lambda43B <- lambda3B
    
    # FOI for age group A: 0-14 years old -- QDENGA
    lambda1AQ <- (1-v_eff_qdenga_sn1)*lambda1A
    lambda2AQ <- (1-v_eff_qdenga_sn2)*lambda2A
    lambda3AQ <- (1-v_eff_qdenga_sn3)*lambda3A
    lambda4AQ <- (1-v_eff_qdenga_sn4)*lambda4A
    lambda12AQ <- (1-v_eff_qdenga_sp2)*lambda12A
    lambda13AQ <- (1-v_eff_qdenga_sp3)*lambda13A
    lambda14AQ <- (1-v_eff_qdenga_sp4)*lambda14A
    lambda21AQ <- (1-v_eff_qdenga_sp1)*lambda21A
    lambda23AQ <- (1-v_eff_qdenga_sp3)*lambda23A
    lambda24AQ <- (1-v_eff_qdenga_sp4)*lambda24A
    lambda31AQ <- (1-v_eff_qdenga_sp1)*lambda31A
    lambda32AQ <- (1-v_eff_qdenga_sp2)*lambda32A
    lambda34AQ <- (1-v_eff_qdenga_sp4)*lambda34A
    lambda41AQ <- (1-v_eff_qdenga_sp1)*lambda41A
    lambda42AQ <- (1-v_eff_qdenga_sp2)*lambda42A
    lambda43AQ <- (1-v_eff_qdenga_sp3)*lambda43A
    
    # FOI for age group B: 15+ years old -- QDENGA
    lambda1BQ <- lambda1B
    lambda2BQ <- lambda2B
    lambda3BQ <- lambda3B
    lambda4BQ <- lambda4B
    lambda12BQ <- lambda2B
    lambda13BQ <- lambda3B
    lambda14BQ <- lambda4B
    lambda21BQ <- lambda1B
    lambda23BQ <- lambda3B
    lambda24BQ <- lambda4B
    lambda31BQ <- lambda1B
    lambda32BQ <- lambda2B
    lambda34BQ <- lambda4B
    lambda41BQ <- lambda1B
    lambda42BQ <- lambda2B
    lambda43BQ <- lambda3B
    
    # FOI for age group A: 0-14 years old -- DENGVAXIA
    lambda1AD <- (1-v_eff_dengvaxia_sn1)*lambda1A
    lambda2AD <- (1-v_eff_dengvaxia_sn2)*lambda2A
    lambda3AD <- (1-v_eff_dengvaxia_sn3)*lambda3A
    lambda4AD <- (1-v_eff_dengvaxia_sn4)*lambda4A
    lambda12AD <- (1-v_eff_dengvaxia_sp2)*lambda12A
    lambda13AD <- (1-v_eff_dengvaxia_sp3)*lambda13A
    lambda14AD <- (1-v_eff_dengvaxia_sp4)*lambda14A
    lambda21AD <- (1-v_eff_dengvaxia_sp1)*lambda21A
    lambda23AD <- (1-v_eff_dengvaxia_sp3)*lambda23A
    lambda24AD <- (1-v_eff_dengvaxia_sp4)*lambda24A
    lambda31AD <- (1-v_eff_dengvaxia_sp1)*lambda31A
    lambda32AD <- (1-v_eff_dengvaxia_sp2)*lambda32A
    lambda34AD <- (1-v_eff_dengvaxia_sp4)*lambda34A
    lambda41AD <- (1-v_eff_dengvaxia_sp1)*lambda41A
    lambda42AD <- (1-v_eff_dengvaxia_sp2)*lambda42A
    lambda43AD <- (1-v_eff_dengvaxia_sp3)*lambda43A
    
    # FOI for age group B: 15+ years old -- DENGVAXIA
    lambda1BD <- lambda1B
    lambda2BD <- lambda2B
    lambda3BD <- lambda3B
    lambda4BD <- lambda4B
    lambda12BD <- lambda2B
    lambda13BD <- lambda3B
    lambda14BD <- lambda4B
    lambda21BD <- lambda1B
    lambda23BD <- lambda3B
    lambda24BD <- lambda4B
    lambda31BD <- lambda1B
    lambda32BD <- lambda2B
    lambda34BD <- lambda4B
    lambda41BD <- lambda1B
    lambda42BD <- lambda2B
    lambda43BD <- lambda3B
    
    birth <- mu_A*PA + mu_B*PB
    
    screening <- 0
    fp <- (1-spec)/(1 - spec + sens) # false positivity rate for Qdenga screening
    tp <- sens/(1 - spec + sens) # true positivity rate for Qdenga screening
    
    ### DENGUE MODEL
    # the age group A: 0-14 years old
    # extended SIR model for the primary infection
    # Totally Naive individuals, i.e., no prior infection with any serotypes
    dSA <- birth + omicron*(SAQ+SAD) - (lambda1A + lambda2A + lambda3A + lambda4A + mu_A + nu)*SA - 
      (screening*fp*Sigma + screening*fp*Pi)*SA # when screening is applied, they are vaccinated only when tested positive
    
    # Infected compartments during primary infection
    dI1A <- lambda1A*SA + omicron*(I1AQ+I1AD) - (tau + mu_A + nu + screening*(Sigma*tp + Pi*tp))*I1A
    dI2A <- lambda2A*SA + omicron*(I2AQ+I2AD) - (tau + mu_A + nu + screening*(Sigma*tp + Pi*tp))*I2A
    dI3A <- lambda3A*SA + omicron*(I3AQ+I3AD) - (tau + mu_A + nu + screening*(Sigma*tp + Pi*tp))*I3A
    dI4A <- lambda4A*SA + omicron*(I4AQ+I4AD) - (tau + mu_A + nu + screening*(Sigma*tp + Pi*tp))*I4A
    
    # Recovered compartments after primary infection
    dR1A <- tau*I1A + omicron*(R1AQ+R1AD) - (omega + mu_A + nu + lambda12A + lambda13A + lambda14A + screening*(Sigma*tp + Pi*tp))*R1A
    dR2A <- tau*I2A + omicron*(R2AQ+R2AD) - (omega + mu_A + nu + lambda21A + lambda23A + lambda24A + screening*(Sigma*tp + Pi*tp))*R2A
    dR3A <- tau*I3A + omicron*(R3AQ+R3AD) - (omega + mu_A + nu + lambda31A + lambda32A + lambda34A + screening*(Sigma*tp + Pi*tp))*R3A
    dR4A <- tau*I4A + omicron*(R4AQ+R4AD) - (omega + mu_A + nu + lambda41A + lambda42A + lambda43A + screening*(Sigma*tp + Pi*tp))*R4A
    
    # extended SIR model for the secondary infection
    # Populations who are now susceptible with the other serotypes after previous infection with DENV    
    dS1A <- omega*R1A + omicron*(S1AQ+S1AD) - (mu_A + nu + lambda12A + lambda13A + lambda14A + screening*(Sigma*tp + Pi*tp))*S1A
    dS2A <- omega*R2A + omicron*(S2AQ+S2AD) - (mu_A + nu + lambda21A + lambda23A + lambda24A + screening*(Sigma*tp + Pi*tp))*S2A
    dS3A <- omega*R3A + omicron*(S3AQ+S3AD) - (mu_A + nu + lambda31A + lambda32A + lambda34A + screening*(Sigma*tp + Pi*tp))*S3A
    dS4A <- omega*R4A + omicron*(S4AQ+S4AD) - (mu_A + nu + lambda41A + lambda42A + lambda43A + screening*(Sigma*tp + Pi*tp))*S4A
    
    # Previously infected with DENV-1 and now are infected with the remaining DENV serotypes
    dI12A <- lambda12A*S1A + omicron*(I12AQ+I12AD) - (tau + mu_A + nu + screening*(Sigma*tp + Pi*tp))*I12A
    dI13A <- lambda13A*S1A + omicron*(I13AQ+I13AD) - (tau + mu_A + nu + screening*(Sigma*tp + Pi*tp))*I13A
    dI14A <- lambda14A*S1A + omicron*(I14AQ+I14AD) - (tau + mu_A + nu + screening*(Sigma*tp + Pi*tp))*I14A
    
    # Previously infected with DENV-2 and now are infected with the remaining DENV serotypes
    dI21A <- lambda21A*S2A + omicron*(I21AQ+I21AD) - (tau + mu_A + nu + screening*(Sigma*tp + Pi*tp))*I21A
    dI23A <- lambda23A*S2A + omicron*(I23AQ+I23AD) - (tau + mu_A + nu + screening*(Sigma*tp + Pi*tp))*I23A
    dI24A <- lambda24A*S2A + omicron*(I24AQ+I24AD) - (tau + mu_A + nu + screening*(Sigma*tp + Pi*tp))*I24A
    
    # Previously infected with DENV-3 and now are infected with the remaining DENV serotypes
    dI31A <- lambda31A*S3A + omicron*(I31AQ+I31AD) - (tau + mu_A + nu + screening*(Sigma*tp + Pi*tp))*I31A
    dI32A <- lambda32A*S3A + omicron*(I32AQ+I32AD) - (tau + mu_A + nu + screening*(Sigma*tp + Pi*tp))*I32A
    dI34A <- lambda34A*S3A + omicron*(I34AQ+I34AD) - (tau + mu_A + nu + screening*(Sigma*tp + Pi*tp))*I34A
    
    # Previously infected with DENV-4 and now are infected with the remaining DENV serotypes
    dI41A <- lambda41A*S4A + omicron*(I41AQ+I41AD) - (tau + mu_A + nu + screening*(Sigma*tp + Pi*tp))*I41A
    dI42A <- lambda42A*S4A + omicron*(I42AQ+I42AD) - (tau + mu_A + nu + screening*(Sigma*tp + Pi*tp))*I42A
    dI43A <- lambda43A*S4A + omicron*(I43AQ+I43AD) - (tau + mu_A + nu + screening*(Sigma*tp + Pi*tp))*I43A
    
    # Previously infected with DENV-1 and now are recovered from the remaining DENV serotypes
    dR12A <- lambda12A*R1A + tau*I12A + omicron*(R12AQ+R12AD) - (mu_A + nu + screening*(Sigma*tp + Pi*tp))*R12A
    dR13A <- lambda13A*R1A + tau*I13A + omicron*(R13AQ+R13AD) - (mu_A + nu + screening*(Sigma*tp + Pi*tp))*R13A
    dR14A <- lambda14A*R1A + tau*I14A + omicron*(R14AQ+R14AD) - (mu_A + nu + screening*(Sigma*tp + Pi*tp))*R14A
    
    # Previously infected with DENV-2 and now are recovered from the remaining DENV serotypes
    dR21A <- lambda21A*R2A + tau*I21A + omicron*(R21AQ+R21AD) - (mu_A + nu + screening*(Sigma*tp + Pi*tp))*R21A
    dR23A <- lambda23A*R2A + tau*I23A + omicron*(R23AQ+R23AD) - (mu_A + nu + screening*(Sigma*tp + Pi*tp))*R23A
    dR24A <- lambda24A*R2A + tau*I24A + omicron*(R24AQ+R24AD) - (mu_A + nu + screening*(Sigma*tp + Pi*tp))*R24A
    
    # Previously infected with DENV-3 and now are recovered from the remaining DENV serotypes
    dR31A <- lambda31A*R3A + tau*I31A + omicron*(R31AQ+R31AD) - (mu_A + nu + screening*(Sigma*tp + Pi*tp))*R31A
    dR32A <- lambda32A*R3A + tau*I32A + omicron*(R32AQ+R32AD) - (mu_A + nu + screening*(Sigma*tp + Pi*tp))*R32A
    dR34A <- lambda34A*R3A + tau*I34A + omicron*(R34AQ+R34AD) - (mu_A + nu + screening*(Sigma*tp + Pi*tp))*R34A
    
    # Previously infected with DENV-4 and now are recovered from the remaining DENV serotypes
    dR41A <- lambda41A*R4A + tau*I41A + omicron*(R41AQ+R41AD) - (mu_A + nu + screening*(Sigma*tp + Pi*tp))*R41A
    dR42A <- lambda42A*R4A + tau*I42A + omicron*(R42AQ+R42AD) - (mu_A + nu + screening*(Sigma*tp + Pi*tp))*R42A
    dR43A <- lambda43A*R4A + tau*I43A + omicron*(R43AQ+R43AD) - (mu_A + nu + screening*(Sigma*tp + Pi*tp))*R43A
    
    # the age group B: 15+ years old
    # extended SIR model for the primary infection
    
    # Totally Naive individuals, i.e., never been infected with any serotypes
    dSB <- nu*SA + omicron*(SBQ+SBD) - (lambda1B + lambda2B + lambda3B + lambda4B + mu_B)*SB
    
    # Infected compartments during primary infection
    dI1B <- nu*I1A + lambda1B*SB + omicron*(I1BQ+I1BD) - (tau + mu_B)*I1B
    dI2B <- nu*I2A + lambda2B*SB + omicron*(I2BQ+I2BD) - (tau + mu_B)*I2B
    dI3B <- nu*I3A + lambda3B*SB + omicron*(I3BQ+I3BD) - (tau + mu_B)*I3B
    dI4B <- nu*I4A + lambda4B*SB + omicron*(I4BQ+I4BD) - (tau + mu_B)*I4B
    
    # Recovered compartments after primary infection
    dR1B <- nu*R1A + tau*I1B + omicron*(R1BQ+R1BD) - (omega + mu_B + lambda12B + lambda13B + lambda14B)*R1B
    dR2B <- nu*R2A + tau*I2B + omicron*(R2BQ+R2BD) - (omega + mu_B + lambda21B + lambda23B + lambda24B)*R2B
    dR3B <- nu*R3A + tau*I3B + omicron*(R3BQ+R3BD) - (omega + mu_B + lambda31B + lambda32B + lambda34B)*R3B
    dR4B <- nu*R4A + tau*I4B + omicron*(R4BQ+R4BD) - (omega + mu_B + lambda41B + lambda42B + lambda43B)*R4B
    
    # extended SIR model for the secondary infection
    # Populations who are now susceptible with the other serotypes after previous infection with DENV
    dS1B <- nu*S1A + omega*R1B + omicron*(S1BQ+S1BD) - (lambda12B + lambda13B + lambda14B + mu_B)*S1B
    dS2B <- nu*S2A + omega*R2B + omicron*(S2BQ+S2BD) - (lambda21B + lambda23B + lambda24B + mu_B)*S2B
    dS3B <- nu*S3A + omega*R3B + omicron*(S3BQ+S3BD) - (lambda31B + lambda32B + lambda34B + mu_B)*S3B
    dS4B <- nu*S4A + omega*R4B + omicron*(S4BQ+S4BD) - (lambda41B + lambda42B + lambda43B + mu_B)*S4B
    
    # Previously infected with DENV-1 and now are infected with the remaining DENV serotypes
    dI12B <- nu*I12A + lambda12B*S1B + omicron*(I12BQ+I12BD) - (tau + mu_B)*I12B
    dI13B <- nu*I13A + lambda13B*S1B + omicron*(I13BQ+I13BD) - (tau + mu_B)*I13B
    dI14B <- nu*I14A + lambda14B*S1B + omicron*(I14BQ+I14BD) - (tau + mu_B)*I14B
    
    # Previously infected with DENV-2 and now are infected with the remaining DENV serotypes
    dI21B <- nu*I21A + lambda21B*S2B + omicron*(I21BQ+I21BD) - (tau + mu_B)*I21B
    dI23B <- nu*I23A + lambda23B*S2B + omicron*(I23BQ+I23BD) - (tau + mu_B)*I23B
    dI24B <- nu*I24A + lambda24B*S2B + omicron*(I24BQ+I24BD) - (tau + mu_B)*I24B
    
    # Previously infected with DENV-3 and now are infected with the remaining DENV serotypes
    dI31B <- nu*I31A + lambda31B*S3B + omicron*(I31BQ+I31BD) - (tau + mu_B)*I31B
    dI32B <- nu*I32A + lambda32B*S3B + omicron*(I32BQ+I32BD) - (tau + mu_B)*I32B
    dI34B <- nu*I34A + lambda34B*S3B + omicron*(I34BQ+I34BD) - (tau + mu_B)*I34B
    
    # Previously infected with DENV-4 and now are infected with the remaining DENV serotypes
    dI41B <- nu*I41A + lambda41B*S4B + omicron*(I41BQ+I41BD) - (tau + mu_B)*I41B
    dI42B <- nu*I42A + lambda42B*S4B + omicron*(I42BQ+I42BD) - (tau + mu_B)*I42B
    dI43B <- nu*I43A + lambda43B*S4B + omicron*(I43BQ+I43BD) - (tau + mu_B)*I43B
    
    # Previously infected with DENV-1 and now are recovered from the remaining DENV serotypes
    dR12B <- nu*R12A + tau*I12B + lambda12B*R1B + omicron*(R12BQ+R12BD) - (mu_B)*R12B 
    dR13B <- nu*R13A + tau*I13B + lambda13B*R1B + omicron*(R13BQ+R13BD) - (mu_B)*R13B 
    dR14B <- nu*R14A + tau*I14B + lambda14B*R1B + omicron*(R14BQ+R14BD) - (mu_B)*R14B
    
    # Previously infected with DENV-2 and now are recovered from the remaining DENV serotypes
    dR21B <- nu*R21A + tau*I21B + lambda21B*R2B + omicron*(R21BQ+R21BD) - (mu_B)*R21B
    dR23B <- nu*R23A + tau*I23B + lambda23B*R2B + omicron*(R23BQ+R23BD) - (mu_B)*R23B 
    dR24B <- nu*R24A + tau*I24B + lambda24B*R2B + omicron*(R24BQ+R24BD) - (mu_B)*R24B
    
    # Previously infected with DENV-3 and now are recovered from the remaining DENV serotypes
    dR31B <- nu*R31A + tau*I31B + lambda31B*R3B + omicron*(R31BQ+R31BD) - (mu_B)*R31B
    dR32B <- nu*R32A + tau*I32B + lambda32B*R3B + omicron*(R32BQ+R32BD) - (mu_B)*R32B
    dR34B <- nu*R34A + tau*I34B + lambda34B*R3B + omicron*(R34BQ+R34BD) - (mu_B)*R34B
    
    # Previously infected with DENV-4 and now are recovered from the remaining DENV serotypes
    dR41B <- nu*R41A + tau*I41B + lambda41B*R4B + omicron*(R41BQ+R41BD) - (mu_B)*R41B 
    dR42B <- nu*R42A + tau*I42B + lambda42B*R4B + omicron*(R42BQ+R42BD) - (mu_B)*R42B
    dR43B <- nu*R43A + tau*I43B + lambda43B*R4B + omicron*(R43BQ+R43BD) - (mu_B)*R43B
    
    # Qdenga vaccination  Model for Age group A (0 to 14 years old)
    dSAQ <- fp*screening*Sigma*SA - (lambda1AQ + lambda2AQ + lambda3AQ + lambda4AQ + mu_A + nu + omicron)*SAQ
    
    dI1AQ <- lambda1AQ*SAQ + tp*screening*Sigma*I1A - (tau + mu_A + nu + omicron)*I1AQ
    dI2AQ <- lambda2AQ*SAQ + tp*screening*Sigma*I2A - (tau + mu_A + nu + omicron)*I2AQ
    dI3AQ <- lambda3AQ*SAQ + tp*screening*Sigma*I3A - (tau + mu_A + nu + omicron)*I3AQ
    dI4AQ <- lambda4AQ*SAQ + tp*screening*Sigma*I4A - (tau + mu_A + nu + omicron)*I4AQ
    
    dR1AQ <- tau*I1AQ + tp*screening*Sigma*R1A - (omega + mu_A + nu + omicron + lambda12AQ + lambda13AQ + lambda14AQ)*R1AQ
    dR2AQ <- tau*I2AQ + tp*screening*Sigma*R2A - (omega + mu_A + nu + omicron + lambda21AQ + lambda23AQ + lambda24AQ)*R2AQ
    dR3AQ <- tau*I3AQ + tp*screening*Sigma*R3A - (omega + mu_A + nu + omicron + lambda31AQ + lambda32AQ + lambda34AQ)*R3AQ
    dR4AQ <- tau*I4AQ + tp*screening*Sigma*R4A - (omega + mu_A + nu + omicron + lambda41AQ + lambda42AQ + lambda43AQ)*R4AQ
    
    dS1AQ <- omega*R1AQ + tp*screening*Sigma*S1A - (mu_A + nu + omicron + lambda12AQ + lambda13AQ + lambda14AQ)*S1AQ
    dS2AQ <- omega*R2AQ + tp*screening*Sigma*S2A - (mu_A + nu + omicron + lambda21AQ + lambda23AQ + lambda24AQ)*S2AQ
    dS3AQ <- omega*R3AQ + tp*screening*Sigma*S3A - (mu_A + nu + omicron + lambda31AQ + lambda32AQ + lambda34AQ)*S3AQ
    dS4AQ <- omega*R4AQ + tp*screening*Sigma*S4A - (mu_A + nu + omicron + lambda41AQ + lambda42AQ + lambda43AQ)*S4AQ
    
    dI12AQ <- lambda12AQ*S1AQ + tp*screening*Sigma*I12A - (tau + mu_A + nu + omicron)*I12AQ
    dI13AQ <- lambda13AQ*S1AQ + tp*screening*Sigma*I13A - (tau + mu_A + nu + omicron)*I13AQ
    dI14AQ <- lambda14AQ*S1AQ + tp*screening*Sigma*I14A - (tau + mu_A + nu + omicron)*I14AQ
    
    dI21AQ <- lambda21AQ*S2AQ + tp*screening*Sigma*I21A - (tau + mu_A + nu + omicron)*I21AQ
    dI23AQ <- lambda23AQ*S2AQ + tp*screening*Sigma*I23A - (tau + mu_A + nu + omicron)*I23AQ
    dI24AQ <- lambda24AQ*S2AQ + tp*screening*Sigma*I24A - (tau + mu_A + nu + omicron)*I24AQ
    
    dI31AQ <- lambda31AQ*S3AQ + tp*screening*Sigma*I31A - (tau + mu_A + nu + omicron)*I31AQ
    dI32AQ <- lambda32AQ*S3AQ + tp*screening*Sigma*I32A - (tau + mu_A + nu + omicron)*I32AQ
    dI34AQ <- lambda34AQ*S3AQ + tp*screening*Sigma*I34A - (tau + mu_A + nu + omicron)*I34AQ
    
    dI41AQ <- lambda41AQ*S4AQ + tp*screening*Sigma*I41A - (tau + mu_A + nu + omicron)*I41AQ
    dI42AQ <- lambda42AQ*S4AQ + tp*screening*Sigma*I42A - (tau + mu_A + nu + omicron)*I42AQ
    dI43AQ <- lambda43AQ*S4AQ + tp*screening*Sigma*I43A - (tau + mu_A + nu + omicron)*I43AQ
    
    dR12AQ <- lambda12AQ*R1AQ + tau*I12AQ + tp*screening*Sigma*R12A - (mu_A + nu + omicron)*R12AQ
    dR13AQ <- lambda13AQ*R1AQ + tau*I13AQ + tp*screening*Sigma*R13A - (mu_A + nu + omicron)*R13AQ
    dR14AQ <- lambda14AQ*R1AQ + tau*I14AQ + tp*screening*Sigma*R14A - (mu_A + nu + omicron)*R14AQ
    
    dR21AQ <- lambda21AQ*R2AQ + tau*I21AQ + tp*screening*Sigma*R21A - (mu_A + nu + omicron)*R21AQ
    dR23AQ <- lambda23AQ*R2AQ + tau*I23AQ + tp*screening*Sigma*R23A - (mu_A + nu + omicron)*R23AQ
    dR24AQ <- lambda24AQ*R2AQ + tau*I24AQ + tp*screening*Sigma*R24A - (mu_A + nu + omicron)*R24AQ
    
    dR31AQ <- lambda31AQ*R3AQ + tau*I31AQ + tp*screening*Sigma*R31A - (mu_A + nu + omicron)*R31AQ
    dR32AQ <- lambda32AQ*R3AQ + tau*I32AQ + tp*screening*Sigma*R32A - (mu_A + nu + omicron)*R32AQ
    dR34AQ <- lambda34AQ*R3AQ + tau*I34AQ + tp*screening*Sigma*R34A - (mu_A + nu + omicron)*R34AQ
    
    dR41AQ <- lambda41AQ*R4AQ + tau*I41AQ + tp*screening*Sigma*R41A - (mu_A + nu + omicron)*R41AQ
    dR42AQ <- lambda42AQ*R4AQ + tau*I42AQ + tp*screening*Sigma*R42A - (mu_A + nu + omicron)*R42AQ
    dR43AQ <- lambda43AQ*R4AQ + tau*I43AQ + tp*screening*Sigma*R43A - (mu_A + nu + omicron)*R43AQ
    
    # Qdenga screeningination Model for Age group B (15+ years old)  
    dSBQ <- nu*SAQ - (lambda1BQ + lambda2BQ + lambda3BQ + lambda4BQ + mu_B + omicron)*SBQ
    
    dI1BQ <- lambda1BQ*SBQ + 0*I1B + nu*I1AQ - (tau + mu_B + omicron)*I1BQ
    dI2BQ <- lambda2BQ*SBQ + 0*I2B + nu*I2AQ - (tau + mu_B + omicron)*I2BQ
    dI3BQ <- lambda3BQ*SBQ + 0*I3B + nu*I3AQ - (tau + mu_B + omicron)*I3BQ
    dI4BQ <- lambda4BQ*SBQ + 0*I4B + nu*I4AQ - (tau + mu_B + omicron)*I4BQ
    
    dR1BQ <- tau*I1BQ + 0*R1B + nu*R1AQ - (omega + mu_B + omicron + lambda12BQ + lambda13BQ + lambda14BQ)*R1BQ
    dR2BQ <- tau*I2BQ + 0*R2B + nu*R2AQ - (omega + mu_B + omicron + lambda21BQ + lambda23BQ + lambda24BQ)*R2BQ
    dR3BQ <- tau*I3BQ + 0*R3B + nu*R3AQ - (omega + mu_B + omicron + lambda31BQ + lambda32BQ + lambda34BQ)*R3BQ
    dR4BQ <- tau*I4BQ + 0*R4B + nu*R4AQ - (omega + mu_B + omicron + lambda41BQ + lambda42BQ + lambda43BQ)*R4BQ
    
    dS1BQ <- omega*R1BQ + 0*S1B + nu*S1AQ - (mu_B + omicron + lambda12BQ + lambda13BQ + lambda14BQ)*S1BQ
    dS2BQ <- omega*R2BQ + 0*S2B + nu*S2AQ - (mu_B + omicron + lambda21BQ + lambda23BQ + lambda24BQ)*S2BQ
    dS3BQ <- omega*R3BQ + 0*S3B + nu*S3AQ - (mu_B + omicron + lambda31BQ + lambda32BQ + lambda34BQ)*S3BQ
    dS4BQ <- omega*R4BQ + 0*S4B + nu*S4AQ - (mu_B + omicron + lambda41BQ + lambda42BQ + lambda43BQ)*S4BQ
    
    dI12BQ <- lambda12BQ*S1BQ + 0*I12B + nu*I12AQ - (tau + mu_B + omicron)*I12BQ
    dI13BQ <- lambda13BQ*S1BQ + 0*I13B + nu*I13AQ - (tau + mu_B + omicron)*I13BQ
    dI14BQ <- lambda14BQ*S1BQ + 0*I14B + nu*I14AQ - (tau + mu_B + omicron)*I14BQ
    
    dI21BQ <- lambda21BQ*S2BQ + 0*I21B + nu*I21AQ - (tau + mu_B + omicron)*I21BQ
    dI23BQ <- lambda23BQ*S2BQ + 0*I23B + nu*I23AQ - (tau + mu_B + omicron)*I23BQ
    dI24BQ <- lambda24BQ*S2BQ + 0*I24B + nu*I24AQ - (tau + mu_B + omicron)*I24BQ
    
    dI31BQ <- lambda31BQ*S3BQ + 0*I31B + nu*I31AQ - (tau + mu_B + omicron)*I31BQ
    dI32BQ <- lambda32BQ*S3BQ + 0*I32B + nu*I32AQ - (tau + mu_B + omicron)*I32BQ
    dI34BQ <- lambda34BQ*S3BQ + 0*I34B + nu*I34AQ - (tau + mu_B + omicron)*I34BQ
    
    dI41BQ <- lambda41BQ*S4BQ + 0*I41B + nu*I41AQ - (tau + mu_B + omicron)*I41BQ
    dI42BQ <- lambda42BQ*S4BQ + 0*I42B + nu*I42AQ - (tau + mu_B + omicron)*I42BQ
    dI43BQ <- lambda43BQ*S4BQ + 0*I43B + nu*I43AQ - (tau + mu_B + omicron)*I43BQ
    
    dR12BQ <- lambda12BQ*R1BQ + tau*I12BQ + 0*R12B + nu*R12AQ - (mu_B + omicron)*R12BQ
    dR13BQ <- lambda13BQ*R1BQ + tau*I13BQ + 0*R13B + nu*R13AQ - (mu_B + omicron)*R13BQ
    dR14BQ <- lambda14BQ*R1BQ + tau*I14BQ + 0*R14B + nu*R14AQ - (mu_B + omicron)*R14BQ
    
    dR21BQ <- lambda21BQ*R2BQ + tau*I21BQ + 0*R21B + nu*R21AQ - (mu_B + omicron)*R21BQ
    dR23BQ <- lambda23BQ*R2BQ + tau*I23BQ + 0*R23B + nu*R23AQ - (mu_B + omicron)*R23BQ
    dR24BQ <- lambda24BQ*R2BQ + tau*I24BQ + 0*R24B + nu*R24AQ - (mu_B + omicron)*R24BQ
    
    dR31BQ <- lambda31BQ*R3BQ + tau*I31BQ + 0*R31B + nu*R31AQ - (mu_B + omicron)*R31BQ
    dR32BQ <- lambda32BQ*R3BQ + tau*I32BQ + 0*R32B + nu*R32AQ - (mu_B + omicron)*R32BQ
    dR34BQ <- lambda34BQ*R3BQ + tau*I34BQ + 0*R34B + nu*R34AQ - (mu_B + omicron)*R34BQ
    
    dR41BQ <- lambda41BQ*R4BQ + tau*I41BQ + 0*R41B + nu*R41AQ - (mu_B + omicron)*R41BQ
    dR42BQ <- lambda42BQ*R4BQ + tau*I42BQ + 0*R42B + nu*R42AQ - (mu_B + omicron)*R42BQ
    dR43BQ <- lambda43BQ*R4BQ + tau*I43BQ + 0*R43B + nu*R43AQ - (mu_B + omicron)*R43BQ
    
    # Dengvaxia vaccination Model for Age group A (0 - 14 years old)    
    dSAD <- fp*screening*Pi*SA - (lambda1AD + lambda2AD + lambda3AD + lambda4AD + mu_A + nu + omicron)*SAD
    
    dI1AD <- lambda1AD*SAD + tp*screening*Pi*I1A - (tau + mu_A + nu + omicron)*I1AD
    dI2AD <- lambda2AD*SAD + tp*screening*Pi*I2A - (tau + mu_A + nu + omicron)*I2AD
    dI3AD <- lambda3AD*SAD + tp*screening*Pi*I3A - (tau + mu_A + nu + omicron)*I3AD
    dI4AD <- lambda4AD*SAD + tp*screening*Pi*I4A - (tau + mu_A + nu + omicron)*I4AD
    
    dR1AD <- tau*I1AD + tp*screening*Pi*R1A - (omega + mu_A + nu + omicron + lambda12AD + lambda13AD + lambda14AD)*R1AD
    dR2AD <- tau*I2AD + tp*screening*Pi*R2A - (omega + mu_A + nu + omicron + lambda21AD + lambda23AD + lambda24AD)*R2AD
    dR3AD <- tau*I3AD + tp*screening*Pi*R3A - (omega + mu_A + nu + omicron + lambda31AD + lambda32AD + lambda34AD)*R3AD
    dR4AD <- tau*I4AD + tp*screening*Pi*R4A - (omega + mu_A + nu + omicron + lambda41AD + lambda42AD + lambda43AD)*R4AD
    
    dS1AD <- omega*R1AD + tp*screening*Pi*S1A - (mu_A + nu + omicron + lambda12AD + lambda13AD + lambda14AD)*S1AD
    dS2AD <- omega*R2AD + tp*screening*Pi*S2A - (mu_A + nu + omicron + lambda21AD + lambda23AD + lambda24AD)*S2AD
    dS3AD <- omega*R3AD + tp*screening*Pi*S3A - (mu_A + nu + omicron + lambda31AD + lambda32AD + lambda34AD)*S3AD
    dS4AD <- omega*R4AD + tp*screening*Pi*S4A - (mu_A + nu + omicron + lambda41AD + lambda42AD + lambda43AD)*S4AD
    
    dI12AD <- lambda12AD*S1AD + tp*screening*Pi*I12A - (tau + mu_A + nu + omicron)*I12AD
    dI13AD <- lambda13AD*S1AD + tp*screening*Pi*I13A - (tau + mu_A + nu + omicron)*I13AD
    dI14AD <- lambda14AD*S1AD + tp*screening*Pi*I14A - (tau + mu_A + nu + omicron)*I14AD
    
    dI21AD <- lambda21AD*S2AD + tp*screening*Pi*I21A - (tau + mu_A + nu + omicron)*I21AD
    dI23AD <- lambda23AD*S2AD + tp*screening*Pi*I23A - (tau + mu_A + nu + omicron)*I23AD
    dI24AD <- lambda24AD*S2AD + tp*screening*Pi*I24A - (tau + mu_A + nu + omicron)*I24AD
    
    dI31AD <- lambda31AD*S3AD + tp*screening*Pi*I31A - (tau + mu_A + nu + omicron)*I31AD
    dI32AD <- lambda32AD*S3AD + tp*screening*Pi*I32A - (tau + mu_A + nu + omicron)*I32AD
    dI34AD <- lambda34AD*S3AD + tp*screening*Pi*I34A - (tau + mu_A + nu + omicron)*I34AD
    
    dI41AD <- lambda41AD*S4AD + tp*screening*Pi*I41A - (tau + mu_A + nu + omicron)*I41AD
    dI42AD <- lambda42AD*S4AD + tp*screening*Pi*I42A - (tau + mu_A + nu + omicron)*I42AD
    dI43AD <- lambda43AD*S4AD + tp*screening*Pi*I43A - (tau + mu_A + nu + omicron)*I43AD
    
    dR12AD <- lambda12AD*R1AD + tau*I12AD + tp*screening*Pi*R12A - (mu_A + nu + omicron)*R12AD
    dR13AD <- lambda13AD*R1AD + tau*I13AD + tp*screening*Pi*R13A - (mu_A + nu + omicron)*R13AD
    dR14AD <- lambda14AD*R1AD + tau*I14AD + tp*screening*Pi*R14A - (mu_A + nu + omicron)*R14AD
    
    dR21AD <- lambda21AD*R2AD + tau*I21AD + tp*screening*Pi*R21A - (mu_A + nu + omicron)*R21AD
    dR23AD <- lambda23AD*R2AD + tau*I23AD + tp*screening*Pi*R23A - (mu_A + nu + omicron)*R23AD
    dR24AD <- lambda24AD*R2AD + tau*I24AD + tp*screening*Pi*R24A - (mu_A + nu + omicron)*R24AD
    
    dR31AD <- lambda31AD*R3AD + tau*I31AD + tp*screening*Pi*R31A - (mu_A + nu + omicron)*R31AD
    dR32AD <- lambda32AD*R3AD + tau*I32AD + tp*screening*Pi*R32A - (mu_A + nu + omicron)*R32AD
    dR34AD <- lambda34AD*R3AD + tau*I34AD + tp*screening*Pi*R34A - (mu_A + nu + omicron)*R34AD
    
    dR41AD <- lambda41AD*R4AD + tau*I41AD + tp*screening*Pi*R41A - (mu_A + nu + omicron)*R41AD
    dR42AD <- lambda42AD*R4AD + tau*I42AD + tp*screening*Pi*R42A - (mu_A + nu + omicron)*R42AD
    dR43AD <- lambda43AD*R4AD + tau*I43AD + tp*screening*Pi*R43A - (mu_A + nu + omicron)*R43AD
    
    # Dengvaxia Vaccination Model for Age group B (15+ years old)    
    dSBD <- nu*SAD - (lambda1BD + lambda2BD + lambda3BD + lambda4BD + mu_B + omicron)*SBD
    
    dI1BD <- lambda1BD*SBD + 0*I1B + nu*I1AD - (tau + mu_B + omicron)*I1BD
    dI2BD <- lambda2BD*SBD + 0*I2B + nu*I2AD - (tau + mu_B + omicron)*I2BD
    dI3BD <- lambda3BD*SBD + 0*I3B + nu*I3AD - (tau + mu_B + omicron)*I3BD
    dI4BD <- lambda4BD*SBD + 0*I4B + nu*I4AD - (tau + mu_B + omicron)*I4BD
    
    dR1BD <- tau*I1BD + 0*R1B + nu*R1AD - (omega + mu_B + omicron + lambda12BD + lambda13BD + lambda14BD)*R1BD
    dR2BD <- tau*I2BD + 0*R2B + nu*R2AD - (omega + mu_B + omicron + lambda21BD + lambda23BD + lambda24BD)*R2BD
    dR3BD <- tau*I3BD + 0*R3B + nu*R3AD - (omega + mu_B + omicron + lambda31BD + lambda32BD + lambda34BD)*R3BD
    dR4BD <- tau*I4BD + 0*R4B + nu*R4AD - (omega + mu_B + omicron + lambda41BD + lambda42BD + lambda43BD)*R4BD
    
    dS1BD <- omega*R1BD + 0*S1B + nu*S1AD - (mu_B + omicron + lambda12BD + lambda13BD + lambda14BD)*S1BD
    dS2BD <- omega*R2BD + 0*S2B + nu*S2AD - (mu_B + omicron + lambda21BD + lambda23BD + lambda24BD)*S2BD
    dS3BD <- omega*R3BD + 0*S3B + nu*S3AD - (mu_B + omicron + lambda31BD + lambda32BD + lambda34BD)*S3BD
    dS4BD <- omega*R4BD + 0*S4B + nu*S4AD - (mu_B + omicron + lambda41BD + lambda42BD + lambda43BD)*S4BD
    
    dI12BD <- lambda12BD*S1BD + 0*I12B + nu*I12AD - (tau + mu_B + omicron)*I12BD
    dI13BD <- lambda13BD*S1BD + 0*I13B + nu*I13AD - (tau + mu_B + omicron)*I13BD
    dI14BD <- lambda14BD*S1BD + 0*I14B + nu*I14AD - (tau + mu_B + omicron)*I14BD
    
    dI21BD <- lambda21BD*S2BD + 0*I21B + nu*I21AD - (tau + mu_B + omicron)*I21BD
    dI23BD <- lambda23BD*S2BD + 0*I23B + nu*I23AD - (tau + mu_B + omicron)*I23BD
    dI24BD <- lambda24BD*S2BD + 0*I24B + nu*I24AD - (tau + mu_B + omicron)*I24BD
    
    dI31BD <- lambda31BD*S3BD + 0*I31B + nu*I31AD - (tau + mu_B + omicron)*I31BD
    dI32BD <- lambda32BD*S3BD + 0*I32B + nu*I32AD - (tau + mu_B + omicron)*I32BD
    dI34BD <- lambda34BD*S3BD + 0*I34B + nu*I34AD - (tau + mu_B + omicron)*I34BD
    
    dI41BD <- lambda41BD*S4BD + 0*I41B + nu*I41AD - (tau + mu_B + omicron)*I41BD
    dI42BD <- lambda42BD*S4BD + 0*I42B + nu*I42AD - (tau + mu_B + omicron)*I42BD
    dI43BD <- lambda43BD*S4BD + 0*I43B + nu*I43AD - (tau + mu_B + omicron)*I43BD
    
    dR12BD <- lambda12BD*R1BD + tau*I12BD + 0*R12B + nu*R12AD - (mu_B + omicron)*R12BD
    dR13BD <- lambda13BD*R1BD + tau*I13BD + 0*R13B + nu*R13AD - (mu_B + omicron)*R13BD
    dR14BD <- lambda14BD*R1BD + tau*I14BD + 0*R14B + nu*R14AD - (mu_B + omicron)*R14BD
    
    dR21BD <- lambda21BD*R2BD + tau*I21BD + 0*R21B + nu*R21AD - (mu_B + omicron)*R21BD
    dR23BD <- lambda23BD*R2BD + tau*I23BD + 0*R23B + nu*R23AD - (mu_B + omicron)*R23BD
    dR24BD <- lambda24BD*R2BD + tau*I24BD + 0*R24B + nu*R24AD - (mu_B + omicron)*R24BD
    
    dR31BD <- lambda31BD*R3BD + tau*I31BD + 0*R31B + nu*R31AD - (mu_B + omicron)*R31BD
    dR32BD <- lambda32BD*R3BD + tau*I32BD + 0*R32B + nu*R32AD - (mu_B + omicron)*R32BD
    dR34BD <- lambda34BD*R3BD + tau*I34BD + 0*R34B + nu*R34AD - (mu_B + omicron)*R34BD
    
    dR41BD <- lambda41BD*R4BD + tau*I41BD + 0*R41B + nu*R41AD - (mu_B + omicron)*R41BD
    dR42BD <- lambda42BD*R4BD + tau*I42BD + 0*R42B + nu*R42AD - (mu_B + omicron)*R42BD
    dR43BD <- lambda43BD*R4BD + tau*I43BD + 0*R43B + nu*R43AD - (mu_B + omicron)*R43BD
    
    ### OUTCOMES
    ## VCD Incindence
    # total VCD incidence = incidence of 1st VCD from both age-group A and B + incidence of 2nd VCD from both age-group A and B
    dCumIncVCDA <- tau*pi_i*(I1A+I2A+I3A+I4A) + tau*pi_ij*(I12A+I13A+I14A + I21A+I23A+I24A + I31A+I32A+I34A + I41A+I42A+I43A) +
      # accounting for vaccination, it is assumed that first infection after vaccination to cause ADE
      tau*pi_ij*(I1AQ+I2AQ+I3AQ+I4AQ + I12AQ+I13AQ+I14AQ + I21AQ+I23AQ+I24AQ + I31AQ+I32AQ+I34AQ + I41AQ+I42AQ+I43AQ) +
      # for qdenga, it is assumed that first infection after vaccination does not cause ADE (assumption 2)
      tau*pi_i*(I1AD+I2AD+I3AD+I4AD) + tau*pi_ij*(I12AD+I13AD+I14AD + I21AD+I23AD+I24AD + I31AD+I32AD+I34AD + I41AD+I42AD+I43AD)
    
    
    dCumIncVCDB <- tau*pi_i*(I1B+I2B+I3B+I4B) + tau*pi_ij*(I12B+I13B+I14B + I21B+I23B+I24B + I31B+I32B+I34B + I41B+I42B+I43B) +
      # accounting for vaccination
      # for dengvaxia, it is assumed that first infection after vaccination to cause ADE
      tau*pi_ij*(I1BQ+I2BQ+I3BQ+I4BQ + I12BQ+I13BQ+I14BQ + I21BQ+I23BQ+I24BQ + I31BQ+I32BQ+I34BQ + I41BQ+I42BQ+I43BQ) +
      # for qdenga, it is assumed that first infection after vaccination does not cause ADE (assumption 2)
      tau*pi_i*(I1BD+I2BD+I3BD+I4BD) + tau*pi_ij*(I12BD+I13BD+I14BD + I21BD+I23BD+I24BD + I31BD+I32BD+I34BD + I41BD+I42BD+I43BD)
    
    ## Hospitalised DHF
    # The proportion of hospitalised primary DHF = probability of being hospitalised when having DHF, h.dhf, x proportion of primary DHF from infection, pi_i
    h.dhf_iA <- h.dhfA*pi_i
    h.dhf_iB <- h.dhfB*pi_i
    
    # The proportion of hospitalised secondary DHF = probability of being hospitalised when having DHF, h.dhf, x proportion of secondary DHF from infection, pi_ij
    h.dhf_ijA <- h.dhfA*pi_ij
    h.dhf_ijB <- h.dhfB*pi_ij
    
    # total hospitalised DHF incidence = incidence of primary hospitalised DHF from both age-group A and B + incidence of secondary hospitalised DHF from both age-group A and B
    dCumIncHospDHFA <- tau*h.dhf_iA*(I1A+I2A+I3A+I4A) + tau*h.dhf_ijA*(I12A+I13A+I14A + I21A+I23A+I24A + I31A+I32A+I34A + I41A+I42A+I43A) +
      # accounting for vaccination, it is assumed that first infection after vaccination to cause ADE
      tau*h.dhf_ijA*(I1AD+I2AD+I3AD+I4AD + I12AD+I13AD+I14AD + I21AD+I23AD+I24AD + I31AD+I32AD+I34AD + I41AD+I42AD+I43AD) +
      # for qdenga, it is assumed that first infection after vaccination does not cause ADE (assumption 2)
      tau*h.dhf_iA*(I1AQ+I2AQ+I3AQ+I4AQ) + tau*h.dhf_ijA*(I12AQ+I13AQ+I14AQ + I21AQ+I23AQ+I24AQ + I31AQ+I32AQ+I34AQ + I41AQ+I42AQ+I43AQ)
    
    dCumIncHospDHFB <- tau*h.dhf_iB*(I1B+I2B+I3B+I4B) + tau*h.dhf_ijB*(I12B+I13B+I14B + I21B+I23B+I24B + I31B+I32B+I34B + I41B+I42B+I43B) +
      # accounting for vaccination
      # for dengvaxia, it is assumed that first infection after vaccination to cause ADE
      tau*h.dhf_ijB*(I1BD+I2BD+I3BD+I4BD + I12BD+I13BD+I14BD + I21BD+I23BD+I24BD + I31BD+I32BD+I34BD + I41BD+I42BD+I43BD) +
      # for qdenga, it is assumed that first infection after vaccination does not cause ADE (assumption 2)
      tau*h.dhf_iB*(I1BQ+I2BQ+I3BQ+I4BQ) + tau*h.dhf_ijB*(I12BQ+I13BQ+I14BQ + I21BQ+I23BQ+I24BQ + I31BQ+I32BQ+I34BQ + I41BQ+I42BQ+I43BQ)
    
    ## Deaths
    # Deaths from cases are assumed to be as a result of DHF infection with the rate of mu_dhf
    # Meanwhile, DHF is resulted from the proportion of infection at time t
    ddeathA <- mu_dhf*h.dhf_iA*(I1A+I2A+I3A+I4A + 
                                  I1AQ+I2AQ+I3AQ+I4AQ + I1AD+I2AD+I3AD+I4AD) + 
      mu_dhf*h.dhf_ijA*(I12A+I13A+I14A + I21A+I23A+I24A + I31A+I32A+I34A + I41A+I42A+I43A +
                          I12AQ+I13AQ+I14AQ + I21AQ+I23AQ+I24AQ + I31AQ+I32AQ+I34AQ + I41AQ+I42AQ+I43AQ +
                          I12AD+I13AD+I14AD + I21AD+I23AD+I24AD + I31AD+I32AD+I34AD + I41AD+I42AD+I43AD)    
    
    ddeathB <- mu_dhf*h.dhf_iB*(I1B+I2B+I3B+I4B + 
                                  I1BQ+I2BQ+I3BQ+I4BQ + I1BD+I2BD+I3BD+I4BD) + 
      mu_dhf*h.dhf_ijB*(I12B+I13B+I14B + I21B+I23B+I24B + I31B+I32B+I34B + I41B+I42B+I43B +
                          I12BQ+I13BQ+I14BQ + I21BQ+I23BQ+I24BQ + I31BQ+I32BQ+I34BQ + I41BQ+I42BQ+I43BQ +
                          I12BD+I13BD+I14BD + I21BD+I23BD+I24BD + I31BD+I32BD+I34BD + I41BD+I42BD+I43BD)    
    
    ### RETURNING SIMULATED OUTPUT FROM THE MODEL
    list(c(dSA, dS1A, dS2A, dS3A, dS4A, 
           dI1A, dI2A, dI3A, dI4A, dI12A, dI13A, dI14A, dI21A, dI23A, dI24A, dI31A, dI32A, dI34A, dI41A, dI42A, dI43A, 
           dR1A, dR2A, dR3A, dR4A, dR12A, dR13A, dR14A, dR21A, dR23A, dR24A, dR31A, dR32A, dR34A, dR41A, dR42A, dR43A, 
           dSB, dS1B, dS2B, dS3B, dS4B, 
           dI1B, dI2B, dI3B, dI4B, dI12B, dI13B, dI14B, dI21B, dI23B, dI24B, dI31B, dI32B, dI34B, dI41B, dI42B, dI43B, 
           dR1B, dR2B, dR3B, dR4B, dR12B, dR13B, dR14B, dR21B, dR23B, dR24B, dR31B, dR32B, dR34B, dR41B, dR42B, dR43B, 
           dSAQ, dS1AQ, dS2AQ, dS3AQ, dS4AQ, 
           dI1AQ, dI2AQ, dI3AQ, dI4AQ, dI12AQ, dI13AQ, dI14AQ, dI21AQ, dI23AQ, dI24AQ, dI31AQ, dI32AQ, dI34AQ, dI41AQ, dI42AQ, dI43AQ, 
           dR1AQ, dR2AQ, dR3AQ, dR4AQ, dR12AQ, dR13AQ, dR14AQ, dR21AQ, dR23AQ, dR24AQ, dR31AQ, dR32AQ, dR34AQ, dR41AQ, dR42AQ, dR43AQ, 
           dSBQ, dS1BQ, dS2BQ, dS3BQ, dS4BQ, 
           dI1BQ, dI2BQ, dI3BQ, dI4BQ, dI12BQ, dI13BQ, dI14BQ, dI21BQ, dI23BQ, dI24BQ, dI31BQ, dI32BQ, dI34BQ, dI41BQ, dI42BQ, dI43BQ, 
           dR1BQ, dR2BQ, dR3BQ, dR4BQ, dR12BQ, dR13BQ, dR14BQ, dR21BQ, dR23BQ, dR24BQ, dR31BQ, dR32BQ, dR34BQ, dR41BQ, dR42BQ, dR43BQ, 
           dSAD, dS1AD, dS2AD, dS3AD, dS4AD, 
           dI1AD, dI2AD, dI3AD, dI4AD, dI12AD, dI13AD, dI14AD, dI21AD, dI23AD, dI24AD, dI31AD, dI32AD, dI34AD, dI41AD, dI42AD, dI43AD, 
           dR1AD, dR2AD, dR3AD, dR4AD, dR12AD, dR13AD, dR14AD, dR21AD, dR23AD, dR24AD, dR31AD, dR32AD, dR34AD, dR41AD, dR42AD, dR43AD, 
           dSBD, dS1BD, dS2BD, dS3BD, dS4BD, 
           dI1BD, dI2BD, dI3BD, dI4BD, dI12BD, dI13BD, dI14BD, dI21BD, dI23BD, dI24BD, dI31BD, dI32BD, dI34BD, dI41BD, dI42BD, dI43BD, 
           dR1BD, dR2BD, dR3BD, dR4BD, dR12BD, dR13BD, dR14BD, dR21BD, dR23BD, dR24BD, dR31BD, dR32BD, dR34BD, dR41BD, dR42BD, dR43BD,
           dCumIncVCDA, dCumIncVCDB, dCumIncHospDHFA, dCumIncHospDHFB, ddeathA, ddeathB, PA_unvacc, PA_vacc, PA, PB, P,
           lambda1AQ, lambda2AQ, lambda3AQ, lambda4AQ, lambda1AD, lambda2AD, lambda3AD, lambda4AD))})
}

{
  # loading compartments and assigning their initial states
  compartments <- read_excel("initial conditions.xlsx",
                             sheet = "double vaccination")
  compartments <- as.data.frame(compartments)
  
  compartments_value.vector <- c(compartments$value)  # taking the value of parameters from data frame into a vector
  compartments_name.vector <- c(compartments$name)    # taking the name of parameters from data frame into a vector
  names(compartments_value.vector) <- compartments_name.vector  # giving the parameters value with name
  
  istate <- compartments_value.vector
  istate <- c(istate, 
              CumIncVCDA = 0, CumIncVCDB = 0, 
              CumIncHospDHFA = 0, CumIncHospDHFB = 0, deathA = 0, deathB = 0,
              lambda1A = 0, lambda2A = 0, lambda3A = 0, lambda4A = 0,
              lambda1AQ = 0, lambda1BQ = 0, lambda1AD = 0, lambda1BD = 0,
              vacc = 0, screening = 0, CumscreennumQ = 0, 
              CumscreennumD = 0, CumvaccnumQ = 0, CumvaccnumD = 0, 
              PA_unvacc = 0, PA_vacc = 0, PA = 0, PB = 0, P=0,
              vaccine.num = 0)
}

# 2. Dengue Vaccination Model ####
# model construction
model.vaccine1 <- function(t, state, parms) {
  with(as.list(c(state, parms)),{
    PA <- SA+S1A+S2A+S3A+S4A+I1A+I2A+I3A+I4A+I12A+I13A+I14A+I21A+I23A+I24A+I31A+I32A+I34A+I41A+I42A+I43A+R1A+R2A+R3A+R4A+R12A+R13A+R14A+R21A+R23A+R24A+R31A+R32A+R34A+R41A+R42A+R43A+SAQ+S1AQ+S2AQ+S3AQ+S4AQ+I1AQ+I2AQ+I3AQ+I4AQ+I12AQ+I13AQ+I14AQ+I21AQ+I23AQ+I24AQ+I31AQ+I32AQ+I34AQ+I41AQ+I42AQ+I43AQ+R1AQ+R2AQ+R3AQ+R4AQ+R12AQ+R13AQ+R14AQ+R21AQ+R23AQ+R24AQ+R31AQ+R32AQ+R34AQ+R41AQ+R42AQ+R43AQ+SAD+S1AD+S2AD+S3AD+S4AD+I1AD+I2AD+I3AD+I4AD+I12AD+I13AD+I14AD+I21AD+I23AD+I24AD+I31AD+I32AD+I34AD+I41AD+I42AD+I43AD+R1AD+R2AD+R3AD+R4AD+R12AD+R13AD+R14AD+R21AD+R23AD+R24AD+R31AD+R32AD+R34AD+R41AD+R42AD+R43AD
    PB <- SB+S1B+S2B+S3B+S4B+I1B+I2B+I3B+I4B+I12B+I13B+I14B+I21B+I23B+I24B+I31B+I32B+I34B+I41B+I42B+I43B+R1B+R2B+R3B+R4B+R12B+R13B+R14B+R21B+R23B+R24B+R31B+R32B+R34B+R41B+R42B+R43B+SBQ+S1BQ+S2BQ+S3BQ+S4BQ+I1BQ+I2BQ+I3BQ+I4BQ+I12BQ+I13BQ+I14BQ+I21BQ+I23BQ+I24BQ+I31BQ+I32BQ+I34BQ+I41BQ+I42BQ+I43BQ+R1BQ+R2BQ+R3BQ+R4BQ+R12BQ+R13BQ+R14BQ+R21BQ+R23BQ+R24BQ+R31BQ+R32BQ+R34BQ+R41BQ+R42BQ+R43BQ+SBD+S1BD+S2BD+S3BD+S4BD+I1BD+I2BD+I3BD+I4BD+I12BD+I13BD+I14BD+I21BD+I23BD+I24BD+I31BD+I32BD+I34BD+I41BD+I42BD+I43BD+R1BD+R2BD+R3BD+R4BD+R12BD+R13BD+R14BD+R21BD+R23BD+R24BD+R31BD+R32BD+R34BD+R41BD+R42BD+R43BD
    P <- PA + PB
    
    # unvaccinated population in age group A
    PA_unvacc <- SA+S1A+S2A+S3A+S4A+
      I1A+I2A+I3A+I4A+I12A+I13A+I14A+I21A+I23A+I24A+I31A+I32A+I34A+I41A+I42A+I43A+
      R1A+R2A+R3A+R4A+R12A+R13A+R14A+R21A+R23A+R24A+R31A+R32A+R34A+R41A+R42A+R43A
    
    PA_vacc <- SAQ+S1AQ+S2AQ+S3AQ+S4AQ+
      I1AQ+I2AQ+I3AQ+I4AQ+I12AQ+I13AQ+I14AQ+I21AQ+I23AQ+I24AQ+I31AQ+I32AQ+I34AQ+I41AQ+I42AQ+I43AQ+
      R1AQ+R2AQ+R3AQ+R4AQ+R12AQ+R13AQ+R14AQ+R21AQ+R23AQ+R24AQ+R31AQ+R32AQ+R34AQ+R41AQ+R42AQ+R43AQ+
      SAD+S1AD+S2AD+S3AD+S4AD+
      I1AD+I2AD+I3AD+I4AD+I12AD+I13AD+I14AD+I21AD+I23AD+I24AD+I31AD+I32AD+I34AD+I41AD+I42AD+I43AD+
      R1AD+R2AD+R3AD+R4AD+R12AD+R13AD+R14AD+R21AD+R23AD+R24AD+R31AD+R32AD+R34AD+R41AD+R42AD+R43AD
    
    ### SEASONALITY
    # accounting seasonality forcing, variability within a year, and inter-annual variability
    seasA <- 1 - epsilonA*cos(2*pi*t/(360) + lagA*30)
    seasB <- 1 - epsilonB*cos(2*pi*t/(360) + lagB*30)
    
    betaA <- R0A*(tau + mu_A)
    betaB <- R0B*(tau + mu_B)
    
    ### FORCE OF INFECTION ESTIMATION
    # infection importation for each group
    ma <- m*PA
    mb <- m*PB
    
    # FOI for age group A: 0-14 years old
    lambda1A <- betaA*seasA*(I1A+I21A+I31A+I41A + I1B+I21B+I31B+I41B + 
                               I1AQ+I21AQ+I31AQ+I41AQ + I1BQ+I21BQ+I31BQ+I41BQ + 
                               I1AD+I21AD+I31AD+I41AD + I1BD+I21BD+I31BD+I41BD + ma)/P
    lambda2A <- betaA*seasA*(I2A+I12A+I32A+I42A + I2B+I12B+I32B+I42B + 
                               I2AQ+I12AQ+I32AQ+I42AQ + I2BQ+I12BQ+I32BQ+I42BQ + 
                               I2AD+I12AD+I32AD+I42AD + I2BD+I12BD+I32BD+I42BD + ma)/P
    lambda3A <- betaA*seasA*(I3A+I13A+I23A+I43A + I3B+I13B+I23B+I43B + 
                               I3AQ+I13AQ+I23AQ+I43AQ + I3BQ+I13BQ+I23BQ+I43BQ + 
                               I3AD+I13AD+I23AD+I43AD + I3BD+I13BD+I23BD+I43BD + ma)/P
    lambda4A <- betaA*seasA*(I4A+I14A+I24A+I34A + I4B+I14B+I24B+I34B + 
                               I4AQ+I14AQ+I24AQ+I34AQ + I4BQ+I14BQ+I24BQ+I34BQ + 
                               I4AD+I14AD+I24AD+I34AD + I4BD+I14BD+I24BD+I34BD + ma)/P
    lambda12A <- lambda2A
    lambda13A <- lambda3A
    lambda14A <- lambda4A
    lambda21A <- lambda1A
    lambda23A <- lambda3A
    lambda24A <- lambda4A
    lambda31A <- lambda1A
    lambda32A <- lambda2A
    lambda34A <- lambda4A
    lambda41A <- lambda1A
    lambda42A <- lambda2A
    lambda43A <- lambda3A
    
    # FOI for age group B: 15+ years old
    lambda1B <- betaB*seasB*(I1A+I21A+I31A+I41A + I1B+I21B+I31B+I41B + 
                               I1AQ+I21AQ+I31AQ+I41AQ + I1BQ+I21BQ+I31BQ+I41BQ + 
                               I1AD+I21AD+I31AD+I41AD + I1BD+I21BD+I31BD+I41BD + mb)/P
    lambda2B <- betaB*seasB*(I2A+I12A+I32A+I42A + I2B+I12B+I32B+I42B + 
                               I2AQ+I12AQ+I32AQ+I42AQ + I2BQ+I12BQ+I32BQ+I42BQ + 
                               I2AD+I12AD+I32AD+I42AD + I2BD+I12BD+I32BD+I42BD + mb)/P
    lambda3B <- betaB*seasB*(I3A+I13A+I23A+I43A + I3B+I13B+I23B+I43B + 
                               I3AQ+I13AQ+I23AQ+I43AQ + I3BQ+I13BQ+I23BQ+I43BQ + 
                               I3AD+I13AD+I23AD+I43AD + I3BD+I13BD+I23BD+I43BD + mb)/P
    lambda4B <- betaB*seasB*(I4A+I14A+I24A+I34A + I4B+I14B+I24B+I34B + 
                               I4AQ+I14AQ+I24AQ+I34AQ + I4BQ+I14BQ+I24BQ+I34BQ + 
                               I4AD+I14AD+I24AD+I34AD + I4BD+I14BD+I24BD+I34BD + mb)/P
    lambda12B <- lambda2B
    lambda13B <- lambda3B
    lambda14B <- lambda4B
    lambda21B <- lambda1B
    lambda23B <- lambda3B
    lambda24B <- lambda4B
    lambda31B <- lambda1B
    lambda32B <- lambda2B
    lambda34B <- lambda4B
    lambda41B <- lambda1B
    lambda42B <- lambda2B
    lambda43B <- lambda3B
    
    # FOI for age group A: 0-14 years old -- QDENGA
    lambda1AQ <- (1-v_eff_qdenga_sn1)*lambda1A
    lambda2AQ <- (1-v_eff_qdenga_sn2)*lambda2A
    lambda3AQ <- (1-v_eff_qdenga_sn3)*lambda3A
    lambda4AQ <- (1-v_eff_qdenga_sn4)*lambda4A
    lambda12AQ <- (1-v_eff_qdenga_sp2)*lambda12A
    lambda13AQ <- (1-v_eff_qdenga_sp3)*lambda13A
    lambda14AQ <- (1-v_eff_qdenga_sp4)*lambda14A
    lambda21AQ <- (1-v_eff_qdenga_sp1)*lambda21A
    lambda23AQ <- (1-v_eff_qdenga_sp3)*lambda23A
    lambda24AQ <- (1-v_eff_qdenga_sp4)*lambda24A
    lambda31AQ <- (1-v_eff_qdenga_sp1)*lambda31A
    lambda32AQ <- (1-v_eff_qdenga_sp2)*lambda32A
    lambda34AQ <- (1-v_eff_qdenga_sp4)*lambda34A
    lambda41AQ <- (1-v_eff_qdenga_sp1)*lambda41A
    lambda42AQ <- (1-v_eff_qdenga_sp2)*lambda42A
    lambda43AQ <- (1-v_eff_qdenga_sp3)*lambda43A
    
    # FOI for age group B: 15+ years old -- QDENGA
    lambda1BQ <- lambda1B
    lambda2BQ <- lambda2B
    lambda3BQ <- lambda3B
    lambda4BQ <- lambda4B
    lambda12BQ <- lambda2B
    lambda13BQ <- lambda3B
    lambda14BQ <- lambda4B
    lambda21BQ <- lambda1B
    lambda23BQ <- lambda3B
    lambda24BQ <- lambda4B
    lambda31BQ <- lambda1B
    lambda32BQ <- lambda2B
    lambda34BQ <- lambda4B
    lambda41BQ <- lambda1B
    lambda42BQ <- lambda2B
    lambda43BQ <- lambda3B
    
    # FOI for age group A: 0-14 years old -- DENGVAXIA
    lambda1AD <- (1-v_eff_dengvaxia_sn1)*lambda1A
    lambda2AD <- (1-v_eff_dengvaxia_sn2)*lambda2A
    lambda3AD <- (1-v_eff_dengvaxia_sn3)*lambda3A
    lambda4AD <- (1-v_eff_dengvaxia_sn4)*lambda4A
    lambda12AD <- (1-v_eff_dengvaxia_sp2)*lambda12A
    lambda13AD <- (1-v_eff_dengvaxia_sp3)*lambda13A
    lambda14AD <- (1-v_eff_dengvaxia_sp4)*lambda14A
    lambda21AD <- (1-v_eff_dengvaxia_sp1)*lambda21A
    lambda23AD <- (1-v_eff_dengvaxia_sp3)*lambda23A
    lambda24AD <- (1-v_eff_dengvaxia_sp4)*lambda24A
    lambda31AD <- (1-v_eff_dengvaxia_sp1)*lambda31A
    lambda32AD <- (1-v_eff_dengvaxia_sp2)*lambda32A
    lambda34AD <- (1-v_eff_dengvaxia_sp4)*lambda34A
    lambda41AD <- (1-v_eff_dengvaxia_sp1)*lambda41A
    lambda42AD <- (1-v_eff_dengvaxia_sp2)*lambda42A
    lambda43AD <- (1-v_eff_dengvaxia_sp3)*lambda43A
    
    # FOI for age group B: 15+ years old -- DENGVAXIA
    lambda1BD <- lambda1B
    lambda2BD <- lambda2B
    lambda3BD <- lambda3B
    lambda4BD <- lambda4B
    lambda12BD <- lambda2B
    lambda13BD <- lambda3B
    lambda14BD <- lambda4B
    lambda21BD <- lambda1B
    lambda23BD <- lambda3B
    lambda24BD <- lambda4B
    lambda31BD <- lambda1B
    lambda32BD <- lambda2B
    lambda34BD <- lambda4B
    lambda41BD <- lambda1B
    lambda42BD <- lambda2B
    lambda43BD <- lambda3B
    
    birth <- mu_A*PA + mu_B*PB
    
    ### DENGUE MODEL
    # the age group A: 0-14 years old
    # extended SIR model for the primary infection
    # Totally Naive individuals, i.e., no prior infection with any serotypes
    dSA <- birth + omicron*(SAQ+SAD) - (lambda1A + lambda2A + lambda3A + lambda4A + mu_A + nu)*SA -
      Qdenga.monthly.mat1[t,"SA"]*switch1 -
      Qdenga.monthly.mat2[t,"SA"]*switch2 -
      Dengvaxia.monthly.mat1[t,"SA"]*switch3
    
    # Infected compartments during primary infection
    dI1A <- lambda1A*SA + omicron*(I1AQ+I1AD) - (tau + mu_A + nu)*I1A - 
      Qdenga.monthly.mat1[t,"I1A"]*switch1 - 
      Qdenga.monthly.mat2[t,"I1A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"I1A"]*switch3
    
    dI2A <- lambda2A*SA + omicron*(I2AQ+I2AD) - (tau + mu_A + nu)*I2A - 
      Qdenga.monthly.mat1[t,"I2A"]*switch1 - 
      Qdenga.monthly.mat2[t,"I2A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"I2A"]*switch3
    
    dI3A <- lambda3A*SA + omicron*(I3AQ+I3AD) - (tau + mu_A + nu)*I3A - 
      Qdenga.monthly.mat1[t,"I3A"]*switch1 - 
      Qdenga.monthly.mat2[t,"I3A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"I3A"]*switch3
    
    dI4A <- lambda4A*SA + omicron*(I4AQ+I4AD) - (tau + mu_A + nu)*I4A - 
      Qdenga.monthly.mat1[t,"I4A"]*switch1 - 
      Qdenga.monthly.mat2[t,"I4A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"I4A"]*switch3
    
    
    # Recovered compartments after primary infection
    dR1A <- tau*I1A + omicron*(R1AQ+R1AD) - (omega + mu_A + nu + lambda12A + lambda13A + lambda14A)*R1A - 
      Qdenga.monthly.mat1[t,"R1A"]*switch1 - 
      Qdenga.monthly.mat2[t,"R1A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"R1A"]*switch3
    
    dR2A <- tau*I2A + omicron*(R2AQ+R2AD) - (omega + mu_A + nu + lambda21A + lambda23A + lambda24A)*R2A - 
      Qdenga.monthly.mat1[t,"R2A"]*switch1 - 
      Qdenga.monthly.mat2[t,"R2A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"R2A"]*switch3
    
    dR3A <- tau*I3A + omicron*(R3AQ+R3AD) - (omega + mu_A + nu + lambda31A + lambda32A + lambda34A)*R3A - 
      Qdenga.monthly.mat1[t,"R3A"]*switch1 - 
      Qdenga.monthly.mat2[t,"R3A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"R3A"]*switch3
    
    dR4A <- tau*I4A + omicron*(R4AQ+R4AD) - (omega + mu_A + nu + lambda41A + lambda42A + lambda43A)*R4A - 
      Qdenga.monthly.mat1[t,"R4A"]*switch1 - 
      Qdenga.monthly.mat2[t,"R4A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"R4A"]*switch3
    
    
    # extended SIR model for the secondary infection
    # Populations who are now susceptible with the other serotypes after previous infection with DENV    
    dS1A <- omega*R1A + omicron*(S1AQ+S1AD) - (mu_A + nu + lambda12A + lambda13A + lambda14A)*S1A - 
      Qdenga.monthly.mat1[t,"S1A"]*switch1 - 
      Qdenga.monthly.mat2[t,"S1A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"S1A"]*switch3
    
    dS2A <- omega*R2A + omicron*(S2AQ+S2AD) - (mu_A + nu + lambda21A + lambda23A + lambda24A)*S2A - 
      Qdenga.monthly.mat1[t,"S2A"]*switch1 - 
      Qdenga.monthly.mat2[t,"S2A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"S2A"]*switch3
    
    dS3A <- omega*R3A + omicron*(S3AQ+S3AD) - (mu_A + nu + lambda31A + lambda32A + lambda34A)*S3A - 
      Qdenga.monthly.mat1[t,"S3A"]*switch1 - 
      Qdenga.monthly.mat2[t,"S3A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"S3A"]*switch3
    
    dS4A <- omega*R4A + omicron*(S4AQ+S4AD) - (mu_A + nu + lambda41A + lambda42A + lambda43A)*S4A - 
      Qdenga.monthly.mat1[t,"S4A"]*switch1 - 
      Qdenga.monthly.mat2[t,"S4A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"S4A"]*switch3
    
    
    # Previously infected with DENV-1 and now are infected with the remaining DENV serotypes
    dI12A <- lambda12A*S1A + omicron*(I12AQ+I12AD) - (tau + mu_A + nu)*I12A - 
      Qdenga.monthly.mat1[t,"I12A"]*switch1 - 
      Qdenga.monthly.mat2[t,"I12A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"I12A"]*switch3
    
    dI13A <- lambda13A*S1A + omicron*(I13AQ+I13AD) - (tau + mu_A + nu)*I13A - 
      Qdenga.monthly.mat1[t,"I13A"]*switch1 - 
      Qdenga.monthly.mat2[t,"I13A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"I13A"]*switch3
    
    dI14A <- lambda14A*S1A + omicron*(I14AQ+I14AD) - (tau + mu_A + nu)*I14A - 
      Qdenga.monthly.mat1[t,"I14A"]*switch1 - 
      Qdenga.monthly.mat2[t,"I14A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"I14A"]*switch3
    
    
    # Previously infected with DENV-2 and now are infected with the remaining DENV serotypes
    dI21A <- lambda21A*S2A + omicron*(I21AQ+I21AD) - (tau + mu_A + nu)*I21A - 
      Qdenga.monthly.mat1[t,"I21A"]*switch1 - 
      Qdenga.monthly.mat2[t,"I21A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"I21A"]*switch3
    
    dI23A <- lambda23A*S2A + omicron*(I23AQ+I23AD) - (tau + mu_A + nu)*I23A - 
      Qdenga.monthly.mat1[t,"I23A"]*switch1 - 
      Qdenga.monthly.mat2[t,"I23A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"I23A"]*switch3
    
    dI24A <- lambda24A*S2A + omicron*(I24AQ+I24AD) - (tau + mu_A + nu)*I24A - 
      Qdenga.monthly.mat1[t,"I24A"]*switch1 - 
      Qdenga.monthly.mat2[t,"I24A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"I24A"]*switch3
    
    
    # Previously infected with DENV-3 and now are infected with the remaining DENV serotypes
    dI31A <- lambda31A*S3A + omicron*(I31AQ+I31AD) - (tau + mu_A + nu)*I31A - 
      Qdenga.monthly.mat1[t,"I31A"]*switch1 - 
      Qdenga.monthly.mat2[t,"I31A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"I31A"]*switch3
    
    dI32A <- lambda32A*S3A + omicron*(I32AQ+I32AD) - (tau + mu_A + nu)*I32A - 
      Qdenga.monthly.mat1[t,"I32A"]*switch1 - 
      Qdenga.monthly.mat2[t,"I32A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"I32A"]*switch3
    
    dI34A <- lambda34A*S3A + omicron*(I34AQ+I34AD) - (tau + mu_A + nu)*I34A - 
      Qdenga.monthly.mat1[t,"I34A"]*switch1 - 
      Qdenga.monthly.mat2[t,"I34A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"I34A"]*switch3
    
    
    # Previously infected with DENV-4 and now are infected with the remaining DENV serotypes
    dI41A <- lambda41A*S4A + omicron*(I41AQ+I41AD) - (tau + mu_A + nu)*I41A - 
      Qdenga.monthly.mat1[t,"I41A"]*switch1 - 
      Qdenga.monthly.mat2[t,"I41A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"I41A"]*switch3
    
    dI42A <- lambda42A*S4A + omicron*(I42AQ+I42AD) - (tau + mu_A + nu)*I42A - 
      Qdenga.monthly.mat1[t,"I42A"]*switch1 - 
      Qdenga.monthly.mat2[t,"I42A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"I42A"]*switch3
    
    dI43A <- lambda43A*S4A + omicron*(I43AQ+I43AD) - (tau + mu_A + nu)*I43A - 
      Qdenga.monthly.mat1[t,"I43A"]*switch1 - 
      Qdenga.monthly.mat2[t,"I43A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"I43A"]*switch3
    
    
    # Previously infected with DENV-1 and now are recovered from the remaining DENV serotypes
    dR12A <- lambda12A*R1A + tau*I12A + omicron*(R12AQ+R12AD) - (mu_A + nu)*R12A - 
      Qdenga.monthly.mat1[t,"R12A"]*switch1 - 
      Qdenga.monthly.mat2[t,"R12A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"R12A"]*switch3
    
    dR13A <- lambda13A*R1A + tau*I13A + omicron*(R13AQ+R13AD) - (mu_A + nu)*R13A - 
      Qdenga.monthly.mat1[t,"R13A"]*switch1 - 
      Qdenga.monthly.mat2[t,"R13A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"R13A"]*switch3
    
    dR14A <- lambda14A*R1A + tau*I14A + omicron*(R14AQ+R14AD) - (mu_A + nu)*R14A - 
      Qdenga.monthly.mat1[t,"R14A"]*switch1 - 
      Qdenga.monthly.mat2[t,"R14A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"R14A"]*switch3
    
    
    # Previously infected with DENV-2 and now are recovered from the remaining DENV serotypes
    dR21A <- lambda21A*R2A + tau*I21A + omicron*(R21AQ+R21AD) - (mu_A + nu)*R21A - 
      Qdenga.monthly.mat1[t,"I21A"]*switch1 - 
      Qdenga.monthly.mat2[t,"I21A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"I21A"]*switch3
    
    dR23A <- lambda23A*R2A + tau*I23A + omicron*(R23AQ+R23AD) - (mu_A + nu)*R23A - 
      Qdenga.monthly.mat1[t,"I23A"]*switch1 - 
      Qdenga.monthly.mat2[t,"I23A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"I23A"]*switch3
    
    dR24A <- lambda24A*R2A + tau*I24A + omicron*(R24AQ+R24AD) - (mu_A + nu)*R24A - 
      Qdenga.monthly.mat1[t,"I24A"]*switch1 - 
      Qdenga.monthly.mat2[t,"I24A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"I24A"]*switch3
    
    
    # Previously infected with DENV-3 and now are recovered from the remaining DENV serotypes
    dR31A <- lambda31A*R3A + tau*I31A + omicron*(R31AQ+R31AD) - (mu_A + nu)*R31A - 
      Qdenga.monthly.mat1[t,"R31A"]*switch1 - 
      Qdenga.monthly.mat2[t,"R31A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"R31A"]*switch3
    
    dR32A <- lambda32A*R3A + tau*I32A + omicron*(R32AQ+R32AD) - (mu_A + nu)*R32A - 
      Qdenga.monthly.mat1[t,"R32A"]*switch1 - 
      Qdenga.monthly.mat2[t,"R32A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"R32A"]*switch3
    
    dR34A <- lambda34A*R3A + tau*I34A + omicron*(R34AQ+R34AD) - (mu_A + nu)*R34A - 
      Qdenga.monthly.mat1[t,"R34A"]*switch1 - 
      Qdenga.monthly.mat2[t,"R34A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"R34A"]*switch3
    
    
    # Previously infected with DENV-4 and now are recovered from the remaining DENV serotypes
    dR41A <- lambda41A*R4A + tau*I41A + omicron*(R41AQ+R41AD) - (mu_A + nu)*R41A - 
      Qdenga.monthly.mat1[t,"R41A"]*switch1 - 
      Qdenga.monthly.mat2[t,"R41A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"R41A"]*switch3
    
    dR42A <- lambda42A*R4A + tau*I42A + omicron*(R42AQ+R42AD) - (mu_A + nu)*R42A - 
      Qdenga.monthly.mat1[t,"R42A"]*switch1 - 
      Qdenga.monthly.mat2[t,"R42A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"R42A"]*switch3
    
    dR43A <- lambda43A*R4A + tau*I43A + omicron*(R43AQ+R43AD) - (mu_A + nu)*R43A - 
      Qdenga.monthly.mat1[t,"R43A"]*switch1 - 
      Qdenga.monthly.mat2[t,"R43A"]*switch2 - 
      Dengvaxia.monthly.mat1[t,"R43A"]*switch3
    
    
    # the age group B: 15+ years old
    # extended SIR model for the primary infection
    
    # Totally Naive individuals, i.e., never been infected with any serotypes
    dSB <- nu*SA + omicron*(SBQ+SBD) - (lambda1B + lambda2B + lambda3B + lambda4B + mu_B)*SB
    
    # Infected compartments during primary infection
    dI1B <- nu*I1A + lambda1B*SB + omicron*(I1BQ+I1BD) - (tau + mu_B)*I1B
    dI2B <- nu*I2A + lambda2B*SB + omicron*(I2BQ+I2BD) - (tau + mu_B)*I2B
    dI3B <- nu*I3A + lambda3B*SB + omicron*(I3BQ+I3BD) - (tau + mu_B)*I3B
    dI4B <- nu*I4A + lambda4B*SB + omicron*(I4BQ+I4BD) - (tau + mu_B)*I4B
    
    # Recovered compartments after primary infection
    dR1B <- nu*R1A + tau*I1B + omicron*(R1BQ+R1BD) - (omega + mu_B + lambda12B + lambda13B + lambda14B)*R1B
    dR2B <- nu*R2A + tau*I2B + omicron*(R2BQ+R2BD) - (omega + mu_B + lambda21B + lambda23B + lambda24B)*R2B
    dR3B <- nu*R3A + tau*I3B + omicron*(R3BQ+R3BD) - (omega + mu_B + lambda31B + lambda32B + lambda34B)*R3B
    dR4B <- nu*R4A + tau*I4B + omicron*(R4BQ+R4BD) - (omega + mu_B + lambda41B + lambda42B + lambda43B)*R4B
    
    # extended SIR model for the secondary infection
    # Populations who are now susceptible with the other serotypes after previous infection with DENV
    dS1B <- nu*S1A + omega*R1B + omicron*(S1BQ+S1BD) - (lambda12B + lambda13B + lambda14B + mu_B)*S1B
    dS2B <- nu*S2A + omega*R2B + omicron*(S2BQ+S2BD) - (lambda21B + lambda23B + lambda24B + mu_B)*S2B
    dS3B <- nu*S3A + omega*R3B + omicron*(S3BQ+S3BD) - (lambda31B + lambda32B + lambda34B + mu_B)*S3B
    dS4B <- nu*S4A + omega*R4B + omicron*(S4BQ+S4BD) - (lambda41B + lambda42B + lambda43B + mu_B)*S4B
    
    # Previously infected with DENV-1 and now are infected with the remaining DENV serotypes
    dI12B <- nu*I12A + lambda12B*S1B + omicron*(I12BQ+I12BD) - (tau + mu_B)*I12B
    dI13B <- nu*I13A + lambda13B*S1B + omicron*(I13BQ+I13BD) - (tau + mu_B)*I13B
    dI14B <- nu*I14A + lambda14B*S1B + omicron*(I14BQ+I14BD) - (tau + mu_B)*I14B
    
    # Previously infected with DENV-2 and now are infected with the remaining DENV serotypes
    dI21B <- nu*I21A + lambda21B*S2B + omicron*(I21BQ+I21BD) - (tau + mu_B)*I21B
    dI23B <- nu*I23A + lambda23B*S2B + omicron*(I23BQ+I23BD) - (tau + mu_B)*I23B
    dI24B <- nu*I24A + lambda24B*S2B + omicron*(I24BQ+I24BD) - (tau + mu_B)*I24B
    
    # Previously infected with DENV-3 and now are infected with the remaining DENV serotypes
    dI31B <- nu*I31A + lambda31B*S3B + omicron*(I31BQ+I31BD) - (tau + mu_B)*I31B
    dI32B <- nu*I32A + lambda32B*S3B + omicron*(I32BQ+I32BD) - (tau + mu_B)*I32B
    dI34B <- nu*I34A + lambda34B*S3B + omicron*(I34BQ+I34BD) - (tau + mu_B)*I34B
    
    # Previously infected with DENV-4 and now are infected with the remaining DENV serotypes
    dI41B <- nu*I41A + lambda41B*S4B + omicron*(I41BQ+I41BD) - (tau + mu_B)*I41B
    dI42B <- nu*I42A + lambda42B*S4B + omicron*(I42BQ+I42BD) - (tau + mu_B)*I42B
    dI43B <- nu*I43A + lambda43B*S4B + omicron*(I43BQ+I43BD) - (tau + mu_B)*I43B
    
    # Previously infected with DENV-1 and now are recovered from the remaining DENV serotypes
    dR12B <- nu*R12A + tau*I12B + lambda12B*R1B + omicron*(R12BQ+R12BD) - (mu_B)*R12B 
    dR13B <- nu*R13A + tau*I13B + lambda13B*R1B + omicron*(R13BQ+R13BD) - (mu_B)*R13B 
    dR14B <- nu*R14A + tau*I14B + lambda14B*R1B + omicron*(R14BQ+R14BD) - (mu_B)*R14B
    
    # Previously infected with DENV-2 and now are recovered from the remaining DENV serotypes
    dR21B <- nu*R21A + tau*I21B + lambda21B*R2B + omicron*(R21BQ+R21BD) - (mu_B)*R21B
    dR23B <- nu*R23A + tau*I23B + lambda23B*R2B + omicron*(R23BQ+R23BD) - (mu_B)*R23B 
    dR24B <- nu*R24A + tau*I24B + lambda24B*R2B + omicron*(R24BQ+R24BD) - (mu_B)*R24B
    
    # Previously infected with DENV-3 and now are recovered from the remaining DENV serotypes
    dR31B <- nu*R31A + tau*I31B + lambda31B*R3B + omicron*(R31BQ+R31BD) - (mu_B)*R31B
    dR32B <- nu*R32A + tau*I32B + lambda32B*R3B + omicron*(R32BQ+R32BD) - (mu_B)*R32B
    dR34B <- nu*R34A + tau*I34B + lambda34B*R3B + omicron*(R34BQ+R34BD) - (mu_B)*R34B
    
    # Previously infected with DENV-4 and now are recovered from the remaining DENV serotypes
    dR41B <- nu*R41A + tau*I41B + lambda41B*R4B + omicron*(R41BQ+R41BD) - (mu_B)*R41B 
    dR42B <- nu*R42A + tau*I42B + lambda42B*R4B + omicron*(R42BQ+R42BD) - (mu_B)*R42B
    dR43B <- nu*R43A + tau*I43B + lambda43B*R4B + omicron*(R43BQ+R43BD) - (mu_B)*R43B
    
    # Qdenga vaccination  Model for Age group A (0 to 14 years old)
    dSAQ <- - (lambda1AQ + lambda2AQ + lambda3AQ + lambda4AQ + mu_A + nu + omicron)*SAQ +
      Qdenga.monthly.mat1[t,"SA"]*switch1 + 
      Qdenga.monthly.mat2[t,"SA"]*switch2 
    
    dI1AQ <- lambda1AQ*SAQ - (tau + mu_A + nu + omicron)*I1AQ +
      Qdenga.monthly.mat1[t,"I1A"]*switch1 + 
      Qdenga.monthly.mat2[t,"I1A"]*switch2 
    
    dI2AQ <- lambda2AQ*SAQ - (tau + mu_A + nu + omicron)*I2AQ +
      Qdenga.monthly.mat1[t,"I2A"]*switch1 + 
      Qdenga.monthly.mat2[t,"I2A"]*switch2 
    
    dI3AQ <- lambda3AQ*SAQ - (tau + mu_A + nu + omicron)*I3AQ +
      Qdenga.monthly.mat1[t,"I3A"]*switch1 + 
      Qdenga.monthly.mat2[t,"I3A"]*switch2 
    
    dI4AQ <- lambda4AQ*SAQ - (tau + mu_A + nu + omicron)*I4AQ +
      Qdenga.monthly.mat1[t,"I4A"]*switch1 + 
      Qdenga.monthly.mat2[t,"I4A"]*switch2 
    
    
    dR1AQ <- tau*I1AQ - (omega + mu_A + nu + omicron + lambda12AQ + lambda13AQ + lambda14AQ)*R1AQ +
      Qdenga.monthly.mat1[t,"R1A"]*switch1 + 
      Qdenga.monthly.mat2[t,"R1A"]*switch2 
    
    dR2AQ <- tau*I2AQ - (omega + mu_A + nu + omicron + lambda21AQ + lambda23AQ + lambda24AQ)*R2AQ +
      Qdenga.monthly.mat1[t,"R2A"]*switch1 + 
      Qdenga.monthly.mat2[t,"R2A"]*switch2 
    
    dR3AQ <- tau*I3AQ - (omega + mu_A + nu + omicron + lambda31AQ + lambda32AQ + lambda34AQ)*R3AQ +
      Qdenga.monthly.mat1[t,"R3A"]*switch1 + 
      Qdenga.monthly.mat2[t,"R3A"]*switch2 
    
    dR4AQ <- tau*I4AQ - (omega + mu_A + nu + omicron + lambda41AQ + lambda42AQ + lambda43AQ)*R4AQ +
      Qdenga.monthly.mat1[t,"R4A"]*switch1 + 
      Qdenga.monthly.mat2[t,"R4A"]*switch2 
    
    
    dS1AQ <- omega*R1AQ - (mu_A + nu + omicron + lambda12AQ + lambda13AQ + lambda14AQ)*S1AQ +
      Qdenga.monthly.mat1[t,"S1A"]*switch1 + 
      Qdenga.monthly.mat2[t,"S1A"]*switch2 
    
    dS2AQ <- omega*R2AQ - (mu_A + nu + omicron + lambda21AQ + lambda23AQ + lambda24AQ)*S2AQ +
      Qdenga.monthly.mat1[t,"S2A"]*switch1 + 
      Qdenga.monthly.mat2[t,"S2A"]*switch2 
    
    dS3AQ <- omega*R3AQ - (mu_A + nu + omicron + lambda31AQ + lambda32AQ + lambda34AQ)*S3AQ +
      Qdenga.monthly.mat1[t,"S3A"]*switch1 + 
      Qdenga.monthly.mat2[t,"S3A"]*switch2 
    
    dS4AQ <- omega*R4AQ - (mu_A + nu + omicron + lambda41AQ + lambda42AQ + lambda43AQ)*S4AQ +
      Qdenga.monthly.mat1[t,"S4A"]*switch1 + 
      Qdenga.monthly.mat2[t,"S4A"]*switch2 
    
    
    dI12AQ <- lambda12AQ*S1AQ - (tau + mu_A + nu + omicron)*I12AQ +
      Qdenga.monthly.mat1[t,"I12A"]*switch1 + 
      Qdenga.monthly.mat2[t,"I12A"]*switch2 
    
    dI13AQ <- lambda13AQ*S1AQ - (tau + mu_A + nu + omicron)*I13AQ +
      Qdenga.monthly.mat1[t,"I13A"]*switch1 + 
      Qdenga.monthly.mat2[t,"I13A"]*switch2 
    
    dI14AQ <- lambda14AQ*S1AQ - (tau + mu_A + nu + omicron)*I14AQ +
      Qdenga.monthly.mat1[t,"I14A"]*switch1 + 
      Qdenga.monthly.mat2[t,"I14A"]*switch2 
    
    
    dI21AQ <- lambda21AQ*S2AQ - (tau + mu_A + nu + omicron)*I21AQ +
      Qdenga.monthly.mat1[t,"I21A"]*switch1 + 
      Qdenga.monthly.mat2[t,"I21A"]*switch2 
    
    dI23AQ <- lambda23AQ*S2AQ - (tau + mu_A + nu + omicron)*I23AQ +
      Qdenga.monthly.mat1[t,"I23A"]*switch1 + 
      Qdenga.monthly.mat2[t,"I23A"]*switch2 
    
    dI24AQ <- lambda24AQ*S2AQ - (tau + mu_A + nu + omicron)*I24AQ +
      Qdenga.monthly.mat1[t,"I24A"]*switch1 + 
      Qdenga.monthly.mat2[t,"I24A"]*switch2 
    
    
    dI31AQ <- lambda31AQ*S3AQ - (tau + mu_A + nu + omicron)*I31AQ +
      Qdenga.monthly.mat1[t,"I31A"]*switch1 + 
      Qdenga.monthly.mat2[t,"I31A"]*switch2 
    
    dI32AQ <- lambda32AQ*S3AQ - (tau + mu_A + nu + omicron)*I32AQ +
      Qdenga.monthly.mat1[t,"I32A"]*switch1 + 
      Qdenga.monthly.mat2[t,"I32A"]*switch2 
    
    dI34AQ <- lambda34AQ*S3AQ - (tau + mu_A + nu + omicron)*I34AQ +
      Qdenga.monthly.mat1[t,"I34A"]*switch1 + 
      Qdenga.monthly.mat2[t,"I34A"]*switch2 
    
    
    dI41AQ <- lambda41AQ*S4AQ - (tau + mu_A + nu + omicron)*I41AQ +
      Qdenga.monthly.mat1[t,"I41A"]*switch1 + 
      Qdenga.monthly.mat2[t,"I41A"]*switch2 
    
    dI42AQ <- lambda42AQ*S4AQ - (tau + mu_A + nu + omicron)*I42AQ +
      Qdenga.monthly.mat1[t,"I42A"]*switch1 + 
      Qdenga.monthly.mat2[t,"I42A"]*switch2 
    
    dI43AQ <- lambda43AQ*S4AQ - (tau + mu_A + nu + omicron)*I43AQ +
      Qdenga.monthly.mat1[t,"I43A"]*switch1 + 
      Qdenga.monthly.mat2[t,"I43A"]*switch2 
    
    
    dR12AQ <- lambda12AQ*R1AQ + tau*I12AQ - (mu_A + nu + omicron)*R12AQ +
      Qdenga.monthly.mat1[t,"R12A"]*switch1 + 
      Qdenga.monthly.mat2[t,"R12A"]*switch2 
    
    dR13AQ <- lambda13AQ*R1AQ + tau*I13AQ - (mu_A + nu + omicron)*R13AQ +
      Qdenga.monthly.mat1[t,"R13A"]*switch1 + 
      Qdenga.monthly.mat2[t,"R13A"]*switch2 
    
    dR14AQ <- lambda14AQ*R1AQ +tau*I14AQ - (mu_A + nu + omicron)*R14AQ +
      Qdenga.monthly.mat1[t,"R14A"]*switch1 + 
      Qdenga.monthly.mat2[t,"R14A"]*switch2 
    
    
    dR21AQ <- lambda21AQ*R2AQ + tau*I21AQ - (mu_A + nu + omicron)*R21AQ +
      Qdenga.monthly.mat1[t,"R21A"]*switch1 + 
      Qdenga.monthly.mat2[t,"R21A"]*switch2 
    
    dR23AQ <- lambda23AQ*R2AQ + tau*I23AQ - (mu_A + nu + omicron)*R23AQ +
      Qdenga.monthly.mat1[t,"R23A"]*switch1 + 
      Qdenga.monthly.mat2[t,"R23A"]*switch2 
    
    dR24AQ <- lambda24AQ*R2AQ + tau*I24AQ - (mu_A + nu + omicron)*R24AQ +
      Qdenga.monthly.mat1[t,"R24A"]*switch1 + 
      Qdenga.monthly.mat2[t,"R24A"]*switch2 
    
    
    dR31AQ <- lambda31AQ*R3AQ + tau*I31AQ - (mu_A + nu + omicron)*R31AQ +
      Qdenga.monthly.mat1[t,"R31A"]*switch1 + 
      Qdenga.monthly.mat2[t,"R31A"]*switch2 
    
    dR32AQ <- lambda32AQ*R3AQ + tau*I32AQ - (mu_A + nu + omicron)*R32AQ +
      Qdenga.monthly.mat1[t,"R32A"]*switch1 + 
      Qdenga.monthly.mat2[t,"R32A"]*switch2 
    
    dR34AQ <- lambda34AQ*R3AQ + tau*I34AQ - (mu_A + nu + omicron)*R34AQ +
      Qdenga.monthly.mat1[t,"R34A"]*switch1 + 
      Qdenga.monthly.mat2[t,"R34A"]*switch2 
    
    
    dR41AQ <- lambda41AQ*R4AQ + tau*I41AQ - (mu_A + nu + omicron)*R41AQ +
      Qdenga.monthly.mat1[t,"R41A"]*switch1 + 
      Qdenga.monthly.mat2[t,"R41A"]*switch2 
    
    dR42AQ <- lambda42AQ*R4AQ + tau*I42AQ - (mu_A + nu + omicron)*R42AQ +
      Qdenga.monthly.mat1[t,"R42A"]*switch1 + 
      Qdenga.monthly.mat2[t,"R42A"]*switch2 
    
    dR43AQ <- lambda43AQ*R4AQ + tau*I43AQ - (mu_A + nu + omicron)*R43AQ +
      Qdenga.monthly.mat1[t,"R43A"]*switch1 + 
      Qdenga.monthly.mat2[t,"R43A"]*switch2 
    
    
    # Qdenga screeningination Model for Age group B (15+ years old)  
    dSBQ <- nu*SAQ - (lambda1BQ + lambda2BQ + lambda3BQ + lambda4BQ + mu_B + omicron)*SBQ
    
    dI1BQ <- lambda1BQ*SBQ + nu*I1AQ - (tau + mu_B + omicron)*I1BQ
    dI2BQ <- lambda2BQ*SBQ + nu*I2AQ - (tau + mu_B + omicron)*I2BQ
    dI3BQ <- lambda3BQ*SBQ + nu*I3AQ - (tau + mu_B + omicron)*I3BQ
    dI4BQ <- lambda4BQ*SBQ + nu*I4AQ - (tau + mu_B + omicron)*I4BQ
    
    dR1BQ <- tau*I1BQ + nu*R1AQ - (omega + mu_B + omicron + lambda12BQ + lambda13BQ + lambda14BQ)*R1BQ
    dR2BQ <- tau*I2BQ + nu*R2AQ - (omega + mu_B + omicron + lambda21BQ + lambda23BQ + lambda24BQ)*R2BQ
    dR3BQ <- tau*I3BQ + nu*R3AQ - (omega + mu_B + omicron + lambda31BQ + lambda32BQ + lambda34BQ)*R3BQ
    dR4BQ <- tau*I4BQ + nu*R4AQ - (omega + mu_B + omicron + lambda41BQ + lambda42BQ + lambda43BQ)*R4BQ
    
    dS1BQ <- omega*R1BQ + nu*S1AQ - (mu_B + omicron + lambda12BQ + lambda13BQ + lambda14BQ)*S1BQ
    dS2BQ <- omega*R2BQ + nu*S2AQ - (mu_B + omicron + lambda21BQ + lambda23BQ + lambda24BQ)*S2BQ
    dS3BQ <- omega*R3BQ + nu*S3AQ - (mu_B + omicron + lambda31BQ + lambda32BQ + lambda34BQ)*S3BQ
    dS4BQ <- omega*R4BQ + nu*S4AQ - (mu_B + omicron + lambda41BQ + lambda42BQ + lambda43BQ)*S4BQ
    
    dI12BQ <- lambda12BQ*S1BQ + nu*I12AQ - (tau + mu_B + omicron)*I12BQ
    dI13BQ <- lambda13BQ*S1BQ + nu*I13AQ - (tau + mu_B + omicron)*I13BQ
    dI14BQ <- lambda14BQ*S1BQ + nu*I14AQ - (tau + mu_B + omicron)*I14BQ
    
    dI21BQ <- lambda21BQ*S2BQ + nu*I21AQ - (tau + mu_B + omicron)*I21BQ
    dI23BQ <- lambda23BQ*S2BQ + nu*I23AQ - (tau + mu_B + omicron)*I23BQ
    dI24BQ <- lambda24BQ*S2BQ + nu*I24AQ - (tau + mu_B + omicron)*I24BQ
    
    dI31BQ <- lambda31BQ*S3BQ + nu*I31AQ - (tau + mu_B + omicron)*I31BQ
    dI32BQ <- lambda32BQ*S3BQ + nu*I32AQ - (tau + mu_B + omicron)*I32BQ
    dI34BQ <- lambda34BQ*S3BQ + nu*I34AQ - (tau + mu_B + omicron)*I34BQ
    
    dI41BQ <- lambda41BQ*S4BQ + nu*I41AQ - (tau + mu_B + omicron)*I41BQ
    dI42BQ <- lambda42BQ*S4BQ + nu*I42AQ - (tau + mu_B + omicron)*I42BQ
    dI43BQ <- lambda43BQ*S4BQ + nu*I43AQ - (tau + mu_B + omicron)*I43BQ
    
    dR12BQ <- lambda12BQ*R1BQ + tau*I12BQ + nu*R12AQ - (mu_B + omicron)*R12BQ
    dR13BQ <- lambda13BQ*R1BQ + tau*I13BQ + nu*R13AQ - (mu_B + omicron)*R13BQ
    dR14BQ <- lambda14BQ*R1BQ + tau*I14BQ + nu*R14AQ - (mu_B + omicron)*R14BQ
    
    dR21BQ <- lambda21BQ*R2BQ + tau*I21BQ + nu*R21AQ - (mu_B + omicron)*R21BQ
    dR23BQ <- lambda23BQ*R2BQ + tau*I23BQ + nu*R23AQ - (mu_B + omicron)*R23BQ
    dR24BQ <- lambda24BQ*R2BQ + tau*I24BQ + nu*R24AQ - (mu_B + omicron)*R24BQ
    
    dR31BQ <- lambda31BQ*R3BQ + tau*I31BQ + nu*R31AQ - (mu_B + omicron)*R31BQ
    dR32BQ <- lambda32BQ*R3BQ + tau*I32BQ + nu*R32AQ - (mu_B + omicron)*R32BQ
    dR34BQ <- lambda34BQ*R3BQ + tau*I34BQ + nu*R34AQ - (mu_B + omicron)*R34BQ
    
    dR41BQ <- lambda41BQ*R4BQ + tau*I41BQ + nu*R41AQ - (mu_B + omicron)*R41BQ
    dR42BQ <- lambda42BQ*R4BQ + tau*I42BQ + nu*R42AQ - (mu_B + omicron)*R42BQ
    dR43BQ <- lambda43BQ*R4BQ + tau*I43BQ + nu*R43AQ - (mu_B + omicron)*R43BQ
    
    # Dengvaxia vaccination Model for Age group A (0 - 14 years old)    
    dSAD <- - (lambda1AD + lambda2AD + lambda3AD + lambda4AD + mu_A + nu + omicron)*SAD + 
      Dengvaxia.monthly.mat1[t,"SA"]*switch3 
    
    dI1AD <- lambda1AD*SAD - (tau + mu_A + nu + omicron)*I1AD + 
      Dengvaxia.monthly.mat1[t,"I1A"]*switch3 
    
    dI2AD <- lambda2AD*SAD - (tau + mu_A + nu + omicron)*I2AD + 
      Dengvaxia.monthly.mat1[t,"I2A"]*switch3 
    
    dI3AD <- lambda3AD*SAD - (tau + mu_A + nu + omicron)*I3AD + 
      Dengvaxia.monthly.mat1[t,"I3A"]*switch3 
    
    dI4AD <- lambda4AD*SAD - (tau + mu_A + nu + omicron)*I4AD + 
      Dengvaxia.monthly.mat1[t,"I4A"]*switch3 
    
    
    dR1AD <- tau*I1AD - (omega + mu_A + nu + omicron + lambda12AD + lambda13AD + lambda14AD)*R1AD + 
      Dengvaxia.monthly.mat1[t,"R1A"]*switch3 
    
    dR2AD <- tau*I2AD - (omega + mu_A + nu + omicron + lambda21AD + lambda23AD + lambda24AD)*R2AD + 
      Dengvaxia.monthly.mat1[t,"R2A"]*switch3 
    
    dR3AD <- tau*I3AD - (omega + mu_A + nu + omicron + lambda31AD + lambda32AD + lambda34AD)*R3AD + 
      Dengvaxia.monthly.mat1[t,"R3A"]*switch3 
    
    dR4AD <- tau*I4AD - (omega + mu_A + nu + omicron + lambda41AD + lambda42AD + lambda43AD)*R4AD + 
      Dengvaxia.monthly.mat1[t,"R4A"]*switch3 
    
    dS1AD <- omega*R1AD - (mu_A + nu + omicron + lambda12AD + lambda13AD + lambda14AD)*S1AD + 
      Dengvaxia.monthly.mat1[t,"S1A"]*switch3 
    
    dS2AD <- omega*R2AD - (mu_A + nu + omicron + lambda21AD + lambda23AD + lambda24AD)*S2AD + 
      Dengvaxia.monthly.mat1[t,"S2A"]*switch3 
    
    dS3AD <- omega*R3AD - (mu_A + nu + omicron + lambda31AD + lambda32AD + lambda34AD)*S3AD + 
      Dengvaxia.monthly.mat1[t,"S3A"]*switch3 
    
    dS4AD <- omega*R4AD - (mu_A + nu + omicron + lambda41AD + lambda42AD + lambda43AD)*S4AD + 
      Dengvaxia.monthly.mat1[t,"S4A"]*switch3 
    
    
    dI12AD <- lambda12AD*S1AD - (tau + mu_A + nu + omicron)*I12AD + 
      Dengvaxia.monthly.mat1[t,"I12A"]*switch3 
    
    dI13AD <- lambda13AD*S1AD - (tau + mu_A + nu + omicron)*I13AD + 
      Dengvaxia.monthly.mat1[t,"I13A"]*switch3 
    
    dI14AD <- lambda14AD*S1AD - (tau + mu_A + nu + omicron)*I14AD + 
      Dengvaxia.monthly.mat1[t,"I14A"]*switch3 
    
    
    dI21AD <- lambda21AD*S2AD - (tau + mu_A + nu + omicron)*I21AD + 
      Dengvaxia.monthly.mat1[t,"I21A"]*switch3 
    
    dI23AD <- lambda23AD*S2AD - (tau + mu_A + nu + omicron)*I23AD + 
      Dengvaxia.monthly.mat1[t,"I23A"]*switch3 
    
    dI24AD <- lambda24AD*S2AD - (tau + mu_A + nu + omicron)*I24AD + 
      Dengvaxia.monthly.mat1[t,"I24A"]*switch3 
    
    
    dI31AD <- lambda31AD*S3AD - (tau + mu_A + nu + omicron)*I31AD + 
      Dengvaxia.monthly.mat1[t,"I31A"]*switch3 
    
    dI32AD <- lambda32AD*S3AD - (tau + mu_A + nu + omicron)*I32AD + 
      Dengvaxia.monthly.mat1[t,"I32A"]*switch3 
    
    dI34AD <- lambda34AD*S3AD - (tau + mu_A + nu + omicron)*I34AD + 
      Dengvaxia.monthly.mat1[t,"I34A"]*switch3 
    
    
    dI41AD <- lambda41AD*S4AD - (tau + mu_A + nu + omicron)*I41AD + 
      Dengvaxia.monthly.mat1[t,"I41A"]*switch3 
    
    dI42AD <- lambda42AD*S4AD - (tau + mu_A + nu + omicron)*I42AD + 
      Dengvaxia.monthly.mat1[t,"I42A"]*switch3 
    
    dI43AD <- lambda43AD*S4AD - (tau + mu_A + nu + omicron)*I43AD + 
      Dengvaxia.monthly.mat1[t,"I43A"]*switch3
    
    dR12AD <- lambda12AD*R1AD + tau*I12AD - (mu_A + nu + omicron)*R12AD + 
      Dengvaxia.monthly.mat1[t,"R12A"]*switch3 
    
    dR13AD <- lambda13AD*R1AD + tau*I13AD - (mu_A + nu + omicron)*R13AD + 
      Dengvaxia.monthly.mat1[t,"R13A"]*switch3 
    
    dR14AD <- lambda14AD*R1AD + tau*I14AD - (mu_A + nu + omicron)*R14AD + 
      Dengvaxia.monthly.mat1[t,"R14A"]*switch3 
    
    
    dR21AD <- lambda21AD*R2AD + tau*I21AD - (mu_A + nu + omicron)*R21AD + 
      Dengvaxia.monthly.mat1[t,"R21A"]*switch3 
    
    dR23AD <- lambda23AD*R2AD + tau*I23AD - (mu_A + nu + omicron)*R23AD + 
      Dengvaxia.monthly.mat1[t,"R23A"]*switch3 
    
    dR24AD <- lambda24AD*R2AD + tau*I24AD - (mu_A + nu + omicron)*R24AD + 
      Dengvaxia.monthly.mat1[t,"R24A"]*switch3 
    
    
    dR31AD <- lambda31AD*R3AD + tau*I31AD - (mu_A + nu + omicron)*R31AD + 
      Dengvaxia.monthly.mat1[t,"R31A"]*switch3 
    
    dR32AD <- lambda32AD*R3AD + tau*I32AD - (mu_A + nu + omicron)*R32AD + 
      Dengvaxia.monthly.mat1[t,"R32A"]*switch3 
    
    dR34AD <- lambda34AD*R3AD + tau*I34AD - (mu_A + nu + omicron)*R34AD + 
      Dengvaxia.monthly.mat1[t,"R34A"]*switch3 
    
    
    dR41AD <- lambda41AD*R4AD + tau*I41AD - (mu_A + nu + omicron)*R41AD + 
      Dengvaxia.monthly.mat1[t,"R41A"]*switch3 
    
    dR42AD <- lambda42AD*R4AD + tau*I42AD - (mu_A + nu + omicron)*R42AD + 
      Dengvaxia.monthly.mat1[t,"R42A"]*switch3 
    
    dR43AD <- lambda43AD*R4AD + tau*I43AD - (mu_A + nu + omicron)*R43AD + 
      Dengvaxia.monthly.mat1[t,"R43A"]*switch3 
    
    
    # Dengvaxia Vaccination Model for Age group B (15+ years old)    
    dSBD <- nu*SAD - (lambda1BD + lambda2BD + lambda3BD + lambda4BD + mu_B + omicron)*SBD
    
    dI1BD <- lambda1BD*SBD + nu*I1AD - (tau + mu_B + omicron)*I1BD
    dI2BD <- lambda2BD*SBD + nu*I2AD - (tau + mu_B + omicron)*I2BD
    dI3BD <- lambda3BD*SBD + nu*I3AD - (tau + mu_B + omicron)*I3BD
    dI4BD <- lambda4BD*SBD + nu*I4AD - (tau + mu_B + omicron)*I4BD
    
    dR1BD <- tau*I1BD + nu*R1AD - (omega + mu_B + omicron + lambda12BD + lambda13BD + lambda14BD)*R1BD
    dR2BD <- tau*I2BD + nu*R2AD - (omega + mu_B + omicron + lambda21BD + lambda23BD + lambda24BD)*R2BD
    dR3BD <- tau*I3BD + nu*R3AD - (omega + mu_B + omicron + lambda31BD + lambda32BD + lambda34BD)*R3BD
    dR4BD <- tau*I4BD + nu*R4AD - (omega + mu_B + omicron + lambda41BD + lambda42BD + lambda43BD)*R4BD
    
    dS1BD <- omega*R1BD + nu*S1AD - (mu_B + omicron + lambda12BD + lambda13BD + lambda14BD)*S1BD
    dS2BD <- omega*R2BD + nu*S2AD - (mu_B + omicron + lambda21BD + lambda23BD + lambda24BD)*S2BD
    dS3BD <- omega*R3BD + nu*S3AD - (mu_B + omicron + lambda31BD + lambda32BD + lambda34BD)*S3BD
    dS4BD <- omega*R4BD + nu*S4AD - (mu_B + omicron + lambda41BD + lambda42BD + lambda43BD)*S4BD
    
    dI12BD <- lambda12BD*S1BD + nu*I12AD - (tau + mu_B + omicron)*I12BD
    dI13BD <- lambda13BD*S1BD + nu*I13AD - (tau + mu_B + omicron)*I13BD
    dI14BD <- lambda14BD*S1BD + nu*I14AD - (tau + mu_B + omicron)*I14BD
    
    dI21BD <- lambda21BD*S2BD + nu*I21AD - (tau + mu_B + omicron)*I21BD
    dI23BD <- lambda23BD*S2BD + nu*I23AD - (tau + mu_B + omicron)*I23BD
    dI24BD <- lambda24BD*S2BD + nu*I24AD - (tau + mu_B + omicron)*I24BD
    
    dI31BD <- lambda31BD*S3BD + nu*I31AD - (tau + mu_B + omicron)*I31BD
    dI32BD <- lambda32BD*S3BD + nu*I32AD - (tau + mu_B + omicron)*I32BD
    dI34BD <- lambda34BD*S3BD + nu*I34AD - (tau + mu_B + omicron)*I34BD
    
    dI41BD <- lambda41BD*S4BD + nu*I41AD - (tau + mu_B + omicron)*I41BD
    dI42BD <- lambda42BD*S4BD + nu*I42AD - (tau + mu_B + omicron)*I42BD
    dI43BD <- lambda43BD*S4BD + nu*I43AD - (tau + mu_B + omicron)*I43BD
    
    dR12BD <- lambda12BD*R1BD + tau*I12BD + nu*R12AD - (mu_B + omicron)*R12BD
    dR13BD <- lambda13BD*R1BD + tau*I13BD + nu*R13AD - (mu_B + omicron)*R13BD
    dR14BD <- lambda14BD*R1BD + tau*I14BD + nu*R14AD - (mu_B + omicron)*R14BD
    
    dR21BD <- lambda21BD*R2BD + tau*I21BD + nu*R21AD - (mu_B + omicron)*R21BD
    dR23BD <- lambda23BD*R2BD + tau*I23BD + nu*R23AD - (mu_B + omicron)*R23BD
    dR24BD <- lambda24BD*R2BD + tau*I24BD + nu*R24AD - (mu_B + omicron)*R24BD
    
    dR31BD <- lambda31BD*R3BD + tau*I31BD + nu*R31AD - (mu_B + omicron)*R31BD
    dR32BD <- lambda32BD*R3BD + tau*I32BD + nu*R32AD - (mu_B + omicron)*R32BD
    dR34BD <- lambda34BD*R3BD + tau*I34BD + nu*R34AD - (mu_B + omicron)*R34BD
    
    dR41BD <- lambda41BD*R4BD + tau*I41BD + nu*R41AD - (mu_B + omicron)*R41BD
    dR42BD <- lambda42BD*R4BD + tau*I42BD + nu*R42AD - (mu_B + omicron)*R42BD
    dR43BD <- lambda43BD*R4BD + tau*I43BD + nu*R43AD - (mu_B + omicron)*R43BD

    ### INTERVENTION
    dCumvaccsn <- Qdenga.monthly.mat1[t, 2]*switch1 +
      Qdenga.monthly.mat2[t, 2]*switch2 +
      Dengvaxia.monthly.mat1[t, 2]*switch3
    
    dCumvaccsp <- sum(Qdenga.monthly.mat1[t, 3:38]*switch1) +
      sum(Qdenga.monthly.mat2[t, 3:38]*switch2) +
      sum(Dengvaxia.monthly.mat1[t, 3:38]*switch3)
    
    ### OUTCOMES
    ## VCD Incindence
    # total VCD incidence = incidence of 1st VCD from both age-group A and B + incidence of 2nd VCD from both age-group A and B
    dCumIncVCDA <- tau*pi_i*(I1A+I2A+I3A+I4A) + tau*pi_ij*(I12A+I13A+I14A + I21A+I23A+I24A + I31A+I32A+I34A + I41A+I42A+I43A) +
      # accounting for vaccination, it is assumed that first infection after vaccination to cause ADE
      tau*pi_ij*(I1AQ+I2AQ+I3AQ+I4AQ + I12AQ+I13AQ+I14AQ + I21AQ+I23AQ+I24AQ + I31AQ+I32AQ+I34AQ + I41AQ+I42AQ+I43AQ) +
      # for qdenga, it is assumed that first infection after vaccination does not cause ADE (assumption 2)
      tau*pi_i*(I1AD+I2AD+I3AD+I4AD) + tau*pi_ij*(I12AD+I13AD+I14AD + I21AD+I23AD+I24AD + I31AD+I32AD+I34AD + I41AD+I42AD+I43AD)
    
    
    dCumIncVCDB <- tau*pi_i*(I1B+I2B+I3B+I4B) + tau*pi_ij*(I12B+I13B+I14B + I21B+I23B+I24B + I31B+I32B+I34B + I41B+I42B+I43B) +
      # accounting for vaccination
      # for dengvaxia, it is assumed that first infection after vaccination to cause ADE
      tau*pi_ij*(I1BQ+I2BQ+I3BQ+I4BQ + I12BQ+I13BQ+I14BQ + I21BQ+I23BQ+I24BQ + I31BQ+I32BQ+I34BQ + I41BQ+I42BQ+I43BQ) +
      # for qdenga, it is assumed that first infection after vaccination does not cause ADE (assumption 2)
      tau*pi_i*(I1BD+I2BD+I3BD+I4BD) + tau*pi_ij*(I12BD+I13BD+I14BD + I21BD+I23BD+I24BD + I31BD+I32BD+I34BD + I41BD+I42BD+I43BD)
    
    ## Hospitalised DHF
    # The proportion of hospitalised primary DHF = probability of being hospitalised when having DHF, h.dhf, x proportion of primary DHF from infection, pi_i
    h.dhf_iA <- h.dhfA*pi_i
    h.dhf_iB <- h.dhfB*pi_i
    
    # The proportion of hospitalised secondary DHF = probability of being hospitalised when having DHF, h.dhf, x proportion of secondary DHF from infection, pi_ij
    h.dhf_ijA <- h.dhfA*pi_ij
    h.dhf_ijB <- h.dhfB*pi_ij
    
    # total hospitalised DHF incidence = incidence of primary hospitalised DHF from both age-group A and B + incidence of secondary hospitalised DHF from both age-group A and B
    dCumIncHospDHFA <- tau*h.dhf_iA*(I1A+I2A+I3A+I4A) + tau*h.dhf_ijA*(I12A+I13A+I14A + I21A+I23A+I24A + I31A+I32A+I34A + I41A+I42A+I43A) +
      # accounting for vaccination, it is assumed that first infection after vaccination to cause ADE
      tau*h.dhf_iA*(I1AD+I2AD+I3AD+I4AD) + tau*h.dhf_ijA*(I12AD+I13AD+I14AD + I21AD+I23AD+I24AD + I31AD+I32AD+I34AD + I41AD+I42AD+I43AD) +
      # for qdenga, it is assumed that first infection after vaccination does not cause ADE (assumption 2)
      tau*h.dhf_iA*(I1AQ+I2AQ+I3AQ+I4AQ) + tau*h.dhf_ijA*(I12AQ+I13AQ+I14AQ + I21AQ+I23AQ+I24AQ + I31AQ+I32AQ+I34AQ + I41AQ+I42AQ+I43AQ)
    
    dCumIncHospDHFB <- tau*h.dhf_iB*(I1B+I2B+I3B+I4B) + tau*h.dhf_ijB*(I12B+I13B+I14B + I21B+I23B+I24B + I31B+I32B+I34B + I41B+I42B+I43B) +
      # accounting for vaccination
      # for dengvaxia, it is assumed that first infection after vaccination to cause ADE
      tau*h.dhf_iB*(I1BD+I2BD+I3BD+I4BD) + tau*h.dhf_ijB*(I12BD+I13BD+I14BD + I21BD+I23BD+I24BD + I31BD+I32BD+I34BD + I41BD+I42BD+I43BD) +
      # for qdenga, it is assumed that first infection after vaccination does not cause ADE (assumption 2)
      tau*h.dhf_iB*(I1BQ+I2BQ+I3BQ+I4BQ) + tau*h.dhf_ijB*(I12BQ+I13BQ+I14BQ + I21BQ+I23BQ+I24BQ + I31BQ+I32BQ+I34BQ + I41BQ+I42BQ+I43BQ)
    
    ## Deaths
    # Deaths from cases are assumed to be as a result of DHF infection with the rate of mu_dhf
    # Meanwhile, DHF is resulted from the proportion of infection at time t
    ddeathA <- mu_dhf*h.dhf_iA*(I1A+I2A+I3A+I4A + 
                                 I1AQ+I2AQ+I3AQ+I4AQ + I1AD+I2AD+I3AD+I4AD) + 
      mu_dhf*h.dhf_ijB*(I12A+I13A+I14A + I21A+I23A+I24A + I31A+I32A+I34A + I41A+I42A+I43A +
                         I12AQ+I13AQ+I14AQ + I21AQ+I23AQ+I24AQ + I31AQ+I32AQ+I34AQ + I41AQ+I42AQ+I43AQ +
                         I12AD+I13AD+I14AD + I21AD+I23AD+I24AD + I31AD+I32AD+I34AD + I41AD+I42AD+I43AD)    
    
    ddeathB <- mu_dhf*h.dhf_iA*(I1B+I2B+I3B+I4B + 
                                 I1BQ+I2BQ+I3BQ+I4BQ + I1BD+I2BD+I3BD+I4BD) + 
      mu_dhf*h.dhf_ijB*(I12B+I13B+I14B + I21B+I23B+I24B + I31B+I32B+I34B + I41B+I42B+I43B +
                         I12BQ+I13BQ+I14BQ + I21BQ+I23BQ+I24BQ + I31BQ+I32BQ+I34BQ + I41BQ+I42BQ+I43BQ +
                         I12BD+I13BD+I14BD + I21BD+I23BD+I24BD + I31BD+I32BD+I34BD + I41BD+I42BD+I43BD)    
    
    ### RETURNING SIMULATED OUTPUT FROM THE MODEL
    list(c(dSA, dS1A, dS2A, dS3A, dS4A, 
           dI1A, dI2A, dI3A, dI4A, dI12A, dI13A, dI14A, dI21A, dI23A, dI24A, dI31A, dI32A, dI34A, dI41A, dI42A, dI43A, 
           dR1A, dR2A, dR3A, dR4A, dR12A, dR13A, dR14A, dR21A, dR23A, dR24A, dR31A, dR32A, dR34A, dR41A, dR42A, dR43A, 
           dSB, dS1B, dS2B, dS3B, dS4B, 
           dI1B, dI2B, dI3B, dI4B, dI12B, dI13B, dI14B, dI21B, dI23B, dI24B, dI31B, dI32B, dI34B, dI41B, dI42B, dI43B, 
           dR1B, dR2B, dR3B, dR4B, dR12B, dR13B, dR14B, dR21B, dR23B, dR24B, dR31B, dR32B, dR34B, dR41B, dR42B, dR43B, 
           dSAQ, dS1AQ, dS2AQ, dS3AQ, dS4AQ, 
           dI1AQ, dI2AQ, dI3AQ, dI4AQ, dI12AQ, dI13AQ, dI14AQ, dI21AQ, dI23AQ, dI24AQ, dI31AQ, dI32AQ, dI34AQ, dI41AQ, dI42AQ, dI43AQ, 
           dR1AQ, dR2AQ, dR3AQ, dR4AQ, dR12AQ, dR13AQ, dR14AQ, dR21AQ, dR23AQ, dR24AQ, dR31AQ, dR32AQ, dR34AQ, dR41AQ, dR42AQ, dR43AQ, 
           dSBQ, dS1BQ, dS2BQ, dS3BQ, dS4BQ, 
           dI1BQ, dI2BQ, dI3BQ, dI4BQ, dI12BQ, dI13BQ, dI14BQ, dI21BQ, dI23BQ, dI24BQ, dI31BQ, dI32BQ, dI34BQ, dI41BQ, dI42BQ, dI43BQ, 
           dR1BQ, dR2BQ, dR3BQ, dR4BQ, dR12BQ, dR13BQ, dR14BQ, dR21BQ, dR23BQ, dR24BQ, dR31BQ, dR32BQ, dR34BQ, dR41BQ, dR42BQ, dR43BQ, 
           dSAD, dS1AD, dS2AD, dS3AD, dS4AD, 
           dI1AD, dI2AD, dI3AD, dI4AD, dI12AD, dI13AD, dI14AD, dI21AD, dI23AD, dI24AD, dI31AD, dI32AD, dI34AD, dI41AD, dI42AD, dI43AD, 
           dR1AD, dR2AD, dR3AD, dR4AD, dR12AD, dR13AD, dR14AD, dR21AD, dR23AD, dR24AD, dR31AD, dR32AD, dR34AD, dR41AD, dR42AD, dR43AD, 
           dSBD, dS1BD, dS2BD, dS3BD, dS4BD, 
           dI1BD, dI2BD, dI3BD, dI4BD, dI12BD, dI13BD, dI14BD, dI21BD, dI23BD, dI24BD, dI31BD, dI32BD, dI34BD, dI41BD, dI42BD, dI43BD, 
           dR1BD, dR2BD, dR3BD, dR4BD, dR12BD, dR13BD, dR14BD, dR21BD, dR23BD, dR24BD, dR31BD, dR32BD, dR34BD, dR41BD, dR42BD, dR43BD,
           dCumIncVCDA, dCumIncVCDB, dCumIncHospDHFA, dCumIncHospDHFB, ddeathA, ddeathB, PA_unvacc, PA_vacc, PA, PB, P, dCumvaccsn, dCumvaccsp,
           lambda1AQ, lambda2AQ, lambda3AQ, lambda4AQ, lambda1AD, lambda2AD, lambda3AD, lambda4AD))})
}

# 3. Dengue Non-Vaccination Model ####
# model construction
model0 <- function(t, state, parms) {
  with(as.list(c(state, parms)),{
    ### TOTAL POPULATION 
    # is the sum of initial population
    PA <- SA+S1A+S2A+S3A+S4A+I1A+I2A+I3A+I4A+I12A+I13A+I14A+I21A+I23A+I24A+I31A+I32A+I34A+I41A+I42A+I43A+R1A+R2A+R3A+R4A+R12A+R13A+R14A+R21A+R23A+R34A+R31A+R32A+R34A+R41A+R42A+R43A
    PB <- SB+S1B+S2B+S3B+S4B+I1B+I2B+I3B+I4B+I12B+I13B+I14B+I21B+I23B+I24B+I31B+I32B+I34B+I41B+I42B+I43B+R1B+R2B+R3B+R4B+R12B+R13B+R14B+R21B+R23B+R34B+R31B+R32B+R34B+R41B+R42B+R43B
    P <- PA + PB
    
    seasA <- 1 - epsilonA*cos(2*pi*t/(360) + lagA*30)
    seasB <- 1 - epsilonB*cos(2*pi*t/(360) + lagB*30)
    
    # R0A <- chiA
    # R0B <- chiB
    
    betaA <- R0A*(tau + mu_A)
    betaB <- R0B*(tau + mu_B)

    # cases imported
    ma <- m*PA
    mb <- m*PB
    
    # FOI for age group A: 0-14 years old
    lambda1A <- betaA*seasA*(I1A+I21A+I31A+I41A + I1B+I21B+I31B+I41B + ma)/P
    lambda2A <- betaA*seasA*(I2A+I12A+I32A+I42A + I2B+I12B+I32B+I42B + ma)/P
    lambda3A <- betaA*seasA*(I3A+I13A+I23A+I43A + I3B+I13B+I23B+I43B + ma)/P
    lambda4A <- betaA*seasA*(I4A+I14A+I24A+I34A + I4B+I14B+I24B+I34B + ma)/P
    lambda12A <- lambda2A
    lambda13A <- lambda3A
    lambda14A <- lambda4A
    lambda21A <- lambda1A
    lambda23A <- lambda3A
    lambda24A <- lambda4A
    lambda31A <- lambda1A
    lambda32A <- lambda2A
    lambda34A <- lambda4A
    lambda41A <- lambda1A
    lambda42A <- lambda2A
    lambda43A <- lambda3A
    
    # FOI for age group B: 15+ years old
    lambda1B <- betaB*seasB*(I1A+I21A+I31A+I41A + I1B+I21B+I31B+I41B + mb)/P
    lambda2B <- betaB*seasB*(I2A+I12A+I32A+I42A + I2B+I12B+I32B+I42B + mb)/P
    lambda3B <- betaB*seasB*(I3A+I13A+I23A+I43A + I3B+I13B+I23B+I43B + mb)/P
    lambda4B <- betaB*seasB*(I4A+I14A+I24A+I34A + I4B+I14B+I24B+I34B + mb)/P
    lambda12B <- lambda2B
    lambda13B <- lambda3B
    lambda14B <- lambda4B
    lambda21B <- lambda1B
    lambda23B <- lambda3B
    lambda24B <- lambda4B
    lambda31B <- lambda1B
    lambda32B <- lambda2B
    lambda34B <- lambda4B
    lambda41B <- lambda1B
    lambda42B <- lambda2B
    lambda43B <- lambda3B
    
    birth <- mu_A*PA + mu_B*PB
    
    ### DENGUE MODEL
    # the age group A: 0-14 years old
    # extended SIR model for the primary infection
    # Totally Naive individuals, i.e., never been infected with any serotypes
    dSA <- birth - (lambda1A + lambda2A + lambda3A + lambda4A + mu_A + nu)*SA
    
    # Infected compartments during primary infection
    dI1A <- lambda1A*SA - (tau + mu_A + nu)*I1A
    dI2A <- lambda2A*SA - (tau + mu_A + nu)*I2A
    dI3A <- lambda3A*SA - (tau + mu_A + nu)*I3A
    dI4A <- lambda4A*SA - (tau + mu_A + nu)*I4A
    
    # Recovered compartments after primary infection
    dR1A <- tau*I1A - (omega + mu_A + nu + lambda12A + lambda13A + lambda14A)*R1A
    dR2A <- tau*I2A - (omega + mu_A + nu + lambda21A + lambda23A + lambda24A)*R2A
    dR3A <- tau*I3A - (omega + mu_A + nu + lambda31A + lambda32A + lambda34A)*R3A
    dR4A <- tau*I4A - (omega + mu_A + nu + lambda41A + lambda42A + lambda43A)*R4A
    
    # extended SIR model for the secondary infection
    # Populations who are now susceptible with the other serotypes after previous infection with DENV    
    dS1A <- omega*R1A - (mu_A + nu + lambda12A + lambda13A + lambda14A)*S1A
    dS2A <- omega*R2A - (mu_A + nu + lambda21A + lambda23A + lambda24A)*S2A
    dS3A <- omega*R3A - (mu_A + nu + lambda31A + lambda32A + lambda34A)*S3A
    dS4A <- omega*R4A - (mu_A + nu + lambda41A + lambda42A + lambda43A)*S4A
    
    # Previously infected with DENV-1 and now are infected with the remaining DENV serotypes
    dI12A <- lambda12A*S1A - (tau + mu_A + nu)*I12A
    dI13A <- lambda13A*S1A - (tau + mu_A + nu)*I13A
    dI14A <- lambda14A*S1A - (tau + mu_A + nu)*I14A
    
    # Previously infected with DENV-2 and now are infected with the remaining DENV serotypes
    dI21A <- lambda21A*S2A - (tau + mu_A + nu)*I21A
    dI23A <- lambda23A*S2A - (tau + mu_A + nu)*I23A
    dI24A <- lambda24A*S2A - (tau + mu_A + nu)*I24A
    
    # Previously infected with DENV-3 and now are infected with the remaining DENV serotypes
    dI31A <- lambda31A*S3A - (tau + mu_A + nu)*I31A
    dI32A <- lambda32A*S3A - (tau + mu_A + nu)*I32A
    dI34A <- lambda34A*S3A - (tau + mu_A + nu)*I34A
    
    # Previously infected with DENV-4 and now are infected with the remaining DENV serotypes
    dI41A <- lambda41A*S4A - (tau + mu_A + nu)*I41A
    dI42A <- lambda42A*S4A - (tau + mu_A + nu)*I42A
    dI43A <- lambda43A*S4A - (tau + mu_A + nu)*I43A
    
    # Previously infected with DENV-1 and now are recovered from the remaining DENV serotypes
    dR12A <- lambda12A*R1A + tau*I12A - (mu_A + nu)*R12A
    dR13A <- lambda13A*R1A + tau*I13A - (mu_A + nu)*R13A
    dR14A <- lambda14A*R1A + tau*I14A - (mu_A + nu)*R14A
    
    # Previously infected with DENV-2 and now are recovered from the remaining DENV serotypes
    dR21A <- lambda21A*R2A + tau*I21A - (mu_A + nu)*R21A
    dR23A <- lambda23A*R2A + tau*I23A - (mu_A + nu)*R23A
    dR24A <- lambda24A*R2A + tau*I24A - (mu_A + nu)*R24A
    
    # Previously infected with DENV-3 and now are recovered from the remaining DENV serotypes
    dR31A <- lambda31A*R3A + tau*I31A - (mu_A + nu)*R31A
    dR32A <- lambda32A*R3A + tau*I32A - (mu_A + nu)*R32A
    dR34A <- lambda34A*R3A + tau*I34A - (mu_A + nu)*R34A
    
    # Previously infected with DENV-4 and now are recovered from the remaining DENV serotypes
    dR41A <- lambda41A*R4A + tau*I41A - (mu_A + nu)*R41A
    dR42A <- lambda42A*R4A + tau*I42A - (mu_A + nu)*R42A
    dR43A <- lambda43A*R4A + tau*I43A - (mu_A + nu)*R43A
    
    # the age group B: 15+ years old
    # extended SIR model for the primary infection
    
    # Totally Naive individuals, i.e., never been infected with any serotypes
    dSB <- nu*SA - (lambda1B + lambda2B + lambda3B + lambda4B + mu_B)*SB
    
    # Infected compartments during primary infection
    dI1B <- nu*I1A + lambda1B*SB - (tau + mu_B)*I1B
    dI2B <- nu*I2A + lambda2B*SB - (tau + mu_B)*I2B
    dI3B <- nu*I3A + lambda3B*SB - (tau + mu_B)*I3B
    dI4B <- nu*I4A + lambda4B*SB - (tau + mu_B)*I4B
    
    # Recovered compartments after primary infection
    dR1B <- nu*R1A + tau*I1B - (omega + mu_B + lambda12B + lambda13B + lambda14B)*R1B
    dR2B <- nu*R2A + tau*I2B - (omega + mu_B + lambda21B + lambda23B + lambda24B)*R2B
    dR3B <- nu*R3A + tau*I3B - (omega + mu_B + lambda31B + lambda32B + lambda34B)*R3B
    dR4B <- nu*R4A + tau*I4B - (omega + mu_B + lambda41B + lambda42B + lambda43B)*R4B
    
    # extended SIR model for the secondary infection
    # Populations who are now susceptible with the other serotypes after previous infection with DENV
    dS1B <- nu*S1A + omega*R1B - (lambda12B + lambda13B + lambda14B + mu_B)*S1B
    dS2B <- nu*S2A + omega*R2B - (lambda21B + lambda23B + lambda24B + mu_B)*S2B
    dS3B <- nu*S3A + omega*R3B - (lambda31B + lambda32B + lambda34B + mu_B)*S3B
    dS4B <- nu*S4A + omega*R4B - (lambda41B + lambda42B + lambda43B + mu_B)*S4B
    
    # Previously infected with DENV-1 and now are infected with the remaining DENV serotypes
    dI12B <- nu*I12A + lambda12B*S1B - (tau + mu_B)*I12B
    dI13B <- nu*I13A + lambda13B*S1B - (tau + mu_B)*I13B
    dI14B <- nu*I14A + lambda14B*S1B - (tau + mu_B)*I14B
    
    # Previously infected with DENV-2 and now are infected with the remaining DENV serotypes
    dI21B <- nu*I21A + lambda21B*S2B - (tau + mu_B)*I21B
    dI23B <- nu*I23A + lambda23B*S2B - (tau + mu_B)*I23B
    dI24B <- nu*I24A + lambda24B*S2B - (tau + mu_B)*I24B
    
    # Previously infected with DENV-3 and now are infected with the remaining DENV serotypes
    dI31B <- nu*I31A + lambda31B*S3B - (tau + mu_B)*I31B
    dI32B <- nu*I32A + lambda32B*S3B - (tau + mu_B)*I32B
    dI34B <- nu*I34A + lambda34B*S3B - (tau + mu_B)*I34B
    
    # Previously infected with DENV-4 and now are infected with the remaining DENV serotypes
    dI41B <- nu*I41A + lambda41B*S4B - (tau + mu_B)*I41B
    dI42B <- nu*I42A + lambda42B*S4B - (tau + mu_B)*I42B
    dI43B <- nu*I43A + lambda43B*S4B - (tau + mu_B)*I43B
    
    # Previously infected with DENV-1 and now are recovered from the remaining DENV serotypes
    dR12B <- nu*R12A + tau*I12B + lambda12B*R1B - mu_B*R12B 
    dR13B <- nu*R13A + tau*I13B + lambda13B*R1B - mu_B*R13B 
    dR14B <- nu*R14A + tau*I14B + lambda14B*R1B - mu_B*R14B
    
    # Previously infected with DENV-2 and now are recovered from the remaining DENV serotypes
    dR21B <- nu*R21A + tau*I21B + lambda21B*R2B - mu_B*R21B
    dR23B <- nu*R23A + tau*I23B + lambda23B*R2B - mu_B*R23B 
    dR24B <- nu*R24A + tau*I24B + lambda24B*R2B - mu_B*R24B
    
    # Previously infected with DENV-3 and now are recovered from the remaining DENV serotypes
    dR31B <- nu*R31A + tau*I31B + lambda31B*R3B - mu_B*R31B
    dR32B <- nu*R32A + tau*I32B + lambda32B*R3B - mu_B*R32B
    dR34B <- nu*R34A + tau*I34B + lambda34B*R3B - mu_B*R34B
    
    # Previously infected with DENV-4 and now are recovered from the remaining DENV serotypes
    dR41B <- nu*R41A + tau*I41B + lambda41B*R4B - mu_B*R41B 
    dR42B <- nu*R42A + tau*I42B + lambda42B*R4B - mu_B*R42B
    dR43B <- nu*R43A + tau*I43B + lambda43B*R4B - mu_B*R43B
    
    ### OUTCOMES
    ## DHF Incindence
    # total DHF incidence = incidence of 1st DHF from both age-group A and B + incidence of 2nd DHF from both age-group A and B
    dCumIncDHFA <- SA*theta_i*(lambda1A+lambda2A+lambda3A+lambda4A) +
      S1A*theta_ij*(lambda12A+lambda13A+lambda14A) +
      S2A*theta_ij*(lambda21A+lambda23A+lambda24A) +
      S3A*theta_ij*(lambda31A+lambda32A+lambda34A) +
      S4A*theta_ij*(lambda41A+lambda42A+lambda43A)
    
    dCumIncDHFB <- SB*theta_i*(lambda1B+lambda2B+lambda3B+lambda4B) +
      S1B*theta_ij*(lambda12B+lambda13B+lambda14B) +
      S2B*theta_ij*(lambda21B+lambda23B+lambda24B) +
      S3B*theta_ij*(lambda31B+lambda32B+lambda34B) +
      S4B*theta_ij*(lambda41B+lambda42B+lambda43B)
    
    # DF Incindence
    # total DF incidence = incidence of 1st DF from both age-group A and B + incidence of 2nd DF from both age-group A and B
    dCumIncDF <- (SA*phi_i)*(lambda1A+lambda2A+lambda3A+lambda4A) + (SB*phi_i)*(lambda1B+lambda2B+lambda3B+lambda4B) +
      (S1A*phi_ij)*(lambda12A+lambda13A+lambda14A) + (S1B*phi_ij)*(lambda12B+lambda13B+lambda14B)
    (S2A*phi_ij)*(lambda21A+lambda23A+lambda24A) + (S2B*phi_ij)*(lambda21B+lambda23B+lambda24B)
    (S3A*phi_ij)*(lambda31A+lambda32A+lambda34A) + (S3B*phi_ij)*(lambda31B+lambda32B+lambda34B)
    (S4A*phi_ij)*(lambda41A+lambda42A+lambda43A) + (S4B*phi_ij)*(lambda41B+lambda42B+lambda43B)
    
    ## VCD Incindence
    # total VCD incidence = incidence of 1st VCD from both age-group A and B + incidence of 2nd VCD from both age-group A and B
    dCumIncVCD <- (SA*pi_i)*(lambda1A+lambda2A+lambda3A+lambda4A) + (SB*pi_i)*(lambda1B+lambda2B+lambda3B+lambda4B) +
      (S1A*pi_ij)*(lambda12A+lambda13A+lambda14A) + (S1B*pi_ij)*(lambda12B+lambda13B+lambda14B)
    (S2A*pi_ij)*(lambda21A+lambda23A+lambda24A) + (S2B*pi_ij)*(lambda21B+lambda23B+lambda24B)
    (S3A*pi_ij)*(lambda31A+lambda32A+lambda34A) + (S3B*pi_ij)*(lambda31B+lambda32B+lambda34B)
    (S4A*pi_ij)*(lambda41A+lambda42A+lambda43A) + (S4B*pi_ij)*(lambda41B+lambda42B+lambda43B)
    
    ## Asymptomatic Dengue Incindence
    # total Asymptomatic Dengue incidence = incidence of 1st Asymptomatic Dengue from both age-group A and B + incidence of 2nd Asymptomatic Dengue from both age-group A and B
    dCumIncAsymptomatic <- (SA*rho_i)*(lambda1A+lambda2A+lambda3A+lambda4A) + (SB*rho_i)*(lambda1B+lambda2B+lambda3B+lambda4B) +
      (S1A*rho_ij)*(lambda12A+lambda13A+lambda14A) + (S1B*rho_ij)*(lambda12B+lambda13B+lambda14B)
    (S2A*rho_ij)*(lambda21A+lambda23A+lambda24A) + (S2B*rho_ij)*(lambda21B+lambda23B+lambda24B)
    (S3A*rho_ij)*(lambda31A+lambda32A+lambda34A) + (S3B*rho_ij)*(lambda31B+lambda32B+lambda34B)
    (S4A*rho_ij)*(lambda41A+lambda42A+lambda43A) + (S4B*rho_ij)*(lambda41B+lambda42B+lambda43B)
    
    ## Hospitalised DHF
    ## Hospitalised DHF
    # The proportion of hospitalised primary DHF = probability of being hospitalised when having DHF, h.dhf, x proportion of primary DHF from infection, pi_i
    h.dhf_iA <- h.dhfA*pi_i
    h.dhf_iB <- h.dhfB*pi_i
    
    # The proportion of hospitalised secondary DHF = probability of being hospitalised when having DHF, h.dhf, x proportion of secondary DHF from infection, pi_ij
    h.dhf_ijA <- h.dhfA*pi_ij
    h.dhf_ijB <- h.dhfB*pi_ij
    
    # total hospitalised DHF incidence = incidence of primary hospitalised DHF from both age-group A and B + incidence of secondary hospitalised DHF from both age-group A and B
    dCumIncHospDHFA <- tau*h.dhf_iA*(I1A+I2A+I3A+I4A) + tau*h.dhf_ijA*(I12A+I13A+I14A + I21A+I23A+I24A + I31A+I32A+I34A + I41A+I42A+I43A)
    
    dCumIncHospDHFB <- tau*h.dhf_iB*(I1B+I2B+I3B+I4B) + tau*h.dhf_ijB*(I12B+I13B+I14B + I21B+I23B+I24B + I31B+I32B+I34B + I41B+I42B+I43B)
    
    ## Deaths
    # Deaths from cases are assumed to be as a result of DHF infection with the rate of mu_dhf
    # Meanwhile, DHF is resulted from the proportion of infection at time t
    ddeathA <- mu_dhf*h.dhf_iA*(I1A+I2A+I3A+I4A) + 
      mu_dhf*h.dhf_ijA*(I12A+I13A+I14A + I21A+I23A+I24A + I31A+I32A+I34A + I41A+I42A+I43A)    
    
    ddeathB <- mu_dhf*h.dhf_iB*(I1B+I2B+I3B+I4B) + 
      mu_dhf*h.dhf_ijB*(I12B+I13B+I14B + I21B+I23B+I24B + I31B+I32B+I34B + I41B+I42B+I43B)    
    
    ### RETURNING SIMULATED OUTPUT FROM THE MODEL
    list(c(dSA, dS1A, dS2A, dS3A, dS4A, 
           dI1A, dI2A, dI3A, dI4A, dI12A, dI13A, dI14A, dI21A, dI23A, dI24A, dI31A, dI32A, dI34A, dI41A, dI42A, dI43A, 
           dR1A, dR2A, dR3A, dR4A, dR12A, dR13A, dR14A, dR21A, dR23A, dR24A, dR31A, dR32A, dR34A, dR41A, dR42A, dR43A, 
           dSB, dS1B, dS2B, dS3B, dS4B, 
           dI1B, dI2B, dI3B, dI4B, dI12B, dI13B, dI14B, dI21B, dI23B, dI24B, dI31B, dI32B, dI34B, dI41B, dI42B, dI43B, 
           dR1B, dR2B, dR3B, dR4B, dR12B, dR13B, dR14B, dR21B, dR23B, dR24B, dR31B, dR32B, dR34B, dR41B, dR42B, dR43B,
           dCumIncDHFA, dCumIncDHFB, dCumIncDF, dCumIncVCD, dCumIncAsymptomatic, dCumIncHospDHFA, dCumIncHospDHFB, ddeathA, ddeathB,
           lambda1A, lambda2A, lambda3A, lambda4A, lambda12A, lambda21A, lambda31A, lambda41A))})
}

# 4. DHF function ####
dhf <- function(dat, period, start, month = T, day = F){
  m <- dat
  
  ## construct the new infection
  if(day == T){  
    num.timepoints <- dim(m)[1]
    
    index.day <- (start - 1):(start + period)
    len <- length(index.day)
    
    index.day1 <- index.month[-c(len)]
    index.day2 <- index.month[-c(1)]
    
    new.infectionsA <- round(c(m[index.day2,"CumIncHospDHFA"] - m[index.day1,"CumIncHospDHFA"]))    # taking the difference of the cumulative incidence returns to the incidence
    
    new.infectionsB <- round(c(m[index.day2,"CumIncHospDHFB"] - m[index.day1,"CumIncHospDHFB"]))    # taking the difference of the cumulative incidence returns to the incidence
  }
  
  # extracting the number of row for the model output
  if(month == T){ 
    index.month <- which((tps %% 30 == 1) & (tps >= (start - 30)) & tps < (period*30 + (start)))
    len <- length(index.month)
    
    index.month1 <- index.month[-c(len)]
    index.month2 <- index.month[-c(1)]
    
    # taking the difference of the cumulative incidence returns to the incidence
    new.infectionsA <- round(c(m[index.month2,"CumIncHospDHFA"] - m[index.month1,"CumIncHospDHFA"]))    # taking the difference of the cumulative incidence returns to the incidence
    
    new.infectionsB <- round(c(m[index.month2,"CumIncHospDHFB"] - m[index.month1,"CumIncHospDHFB"]))    # taking the difference of the cumulative incidence returns to the incidence
  }
  
  df <- data.frame("DHF A" = new.infectionsA, "DHF B" = new.infectionsB)
  return(df)
}

# 5. NLL ####
calcNLL <- function(params, month = T, day = F){
  # returns NLL assuming possion distribution
  # with initial conditions
  parms["h.dhfA"] <- params[1]
  parms["h.dhfB"] <- params[2]
  # constructing the model 
  parms["R0A"] <- params[3]
  parms["R0B"] <- params[4]
  parms["epsilonA"] <- params[5]
  parms["epsilonB"] <- params[6]
  parms["lagA"] <- params[7]
  parms["lagB"] <- params[8]
  
  model_output.df <- (deSolve::lsoda(y = istate0, times = tps, func = model0, parms = parms))
  
  ## construct the new infection: lambda
  # extracting the number of row for the model output
  num.timepoints <- dim(model_output.df)[1]
  
  # taking the difference of the cumulative incidence returns to the incidence
  m <- model_output.df
  
  df <- dhf(m, period = period,  start = 360*40)
  
  lambdaA <- df$DHF.A
  lambdaB <- df$DHF.B
  xa <- rep(avg.A, period/12)
  xb <- rep(avg.B, period/12)
  
  ## returning the output NLL value given the inputs
  # the output fitted is at equilibrium, thus the first 45 years are discarded
  # then, monthly cases in the last 5 years are then fitted to the data
  NLLA <- -sum(dpois(round(xa), lambdaA, log = TRUE))
  NLLB <- -sum(dpois(round(xb), lambdaB, log = TRUE))
  NLL <- NLLA + NLLB
  
  # par(mfrow = c(1,2))
  # plot(y = casesA*100000/npop.jak0.14, x = 1:60, col = "gray", ylab = "People per 100000", xlab = "month",
  #      main = "age 0-14")
  # lines(casesA*100000/npop.jak0.14, col= "gray")
  # lines(xa*100000/npop.jak0.14, type = "b", col= "red")
  # lines(df$DHF.A*100000/npop.jak0.14, col = "blue")
  # 
  # 
  # plot(y = casesB*100000/npop.jak15.75, x = 1:60, col = "gray", ylab = "People per 100000", xlab = "month",
  #      main = "age 15-75+")
  # lines(casesB*100000/npop.jak15.75, col= "gray")
  # lines(xb*100000/npop.jak15.75, type = "b", col= "red")
  # lines(df$DHF.B*100000/npop.jak15.75, col = "blue")
   
  print(c(params[1], params[2], params[3], params[4],
        params[5], params[6], params[7], params[8]))
  print(NLL)
  return(NLL)
}

# 6. ALLOCATING VACCINES AND SCREENING ####
### MODELLING VACCINATION AS A RATE
vacc <- function(target_cov, spec, sens, baseline = model_output0,
                 screening = T, dengvaxia = 0.5, qdenga = 0.5){
  
  model_output <- model_output0[tps,]
  
  timeline.vector <- (1:time_stop >= time_start_interv & 1:time_stop < (time_start_interv + time_campaign))
  
  SA0 <- (model_output[1,"SA"])            # SA at starting condition, taking the average of 12-months of vaccination period
  PA0 <- (sum(model_output[1,2:38]))   # PA at starting condition, taking the average of 12-months of vaccination period
  SPA <- 1 - (SA0/PA0)                   # SPA = 1 - SA/P at starting condition
  vac_cov <- target_cov*SPA
  vaccine.num <- vac_cov*(PA0)
  
  # vaccine allocation per month given the time campaign
  vac.alloc <- vaccine.num/time_campaign*timeline.vector
  vacc.rate <- vac.alloc/PA0
  
  ### MODELLING SCREENING AS A RATE
  # as the vaccination coverage depends on the screening coverage and the diagnostic
  # performance (e.g. false positivity rate and true positivity rate), hence, to
  # calculate how much is the screening coverage to achieve such number of vaccination
  # to be delivered in the population can be calculated as follows:
  # target_coverage_new = vaccine_coverage/SPA_new
  
  # diagnostic performance
  # false positive for those who are seronegative, yet tested as seropositive and  receive vaccine
  fp <- (1-spec)/(1 - spec + sens) # false positivity rate for Qdenga screening
  # true positive for those who are seropositive, and testes as seropositive and receive vaccine
  tp <- sens/(1 - spec + sens) # true positivity rate for Qdenga screening
  
  SPA_new <- mean(fp*model_output[timeline.vector,"SA"] + tp*rowSums(model_output[timeline.vector,3:38]))
  
  target_coverage_new <- target_cov*SPA*PA0/(SPA_new)
  
  # now the new target coverage become the new target coverage for screening
  screen.coverage0 <- target_coverage_new
  
  # as the screening coverage cannot be more than 1, hence:
  screen.coverage <- min(screen.coverage0,1)
  screen.coverage
  
  # calculating the number of screening needed to be performed
  screen.num <- screen.coverage*PA0
  
  screen.per.month <- screen.num/time_campaign*timeline.vector
  screen.rate <- screen.per.month/PA0
  
  ## screen rate for SA
  parms["Sigma"]
  parms["Pi"]
  # vaccine allocated to seronegative individuals
  snQ.vec <- screen.rate*fp*model_output[,"SA"]*parms["Sigma"]*timeline.vector # qdenga
  snD.vec <- screen.rate*fp*model_output[,"SA"]*parms["Pi"]*timeline.vector # dengvaxia
  # vaccine allocated to seropositive individuals
  spQ.vec <- screen.rate*tp*rowSums(model_output[,3:38])*parms["Sigma"]*timeline.vector # qdenga
  spD.vec <- screen.rate*tp*rowSums(model_output[,3:38])*parms["Pi"]*timeline.vector # dengvaxia
  
  sum(snQ.vec, snD.vec, spQ.vec, spD.vec) # estimated allocated vaccine at start = 1792146
  
  # at this point, it has ben demonstrated that to achieve a vaccination coverage
  # of vac_cov (68.0%) , the screening coverage of screen.coverage (92.2%) is needed
  
  # a. screening scenario
  # creating a matrix of vaccine allocation to corresponding to compartments at each timestep
  # QDENGA
  Qdenga.monthly.mat1 <- matrix(0, nrow = time_stop, ncol =38)
  for(i in 1:time_stop){
    # for seronegative as there is only 1 compartment representing it, "SA
    Qdenga.monthly.mat1[i,2] <- snQ.vec[i]
    
    for(j in 3:38){
      # this is for the rest of the compartments (the seropositive compartments)
      Qdenga.monthly.mat1[i,j] <- model_output[i,j]*spQ.vec[i]*timeline.vector[i]/sum(model_output[i,3:38])
    }
  }
  
  colnames(Qdenga.monthly.mat1) <- dimnames(model_output)[[2]][1:38]
  rownames(Qdenga.monthly.mat1) <- 1:time_stop
  
  #DENGVAXIA
  Dengvaxia.monthly.mat1 <- matrix(0, nrow = time_stop, ncol =38)
  for(i in 1:time_stop){
    # for seronegative as there is only 1 compartment representing it, "SA
    
    Dengvaxia.monthly.mat1[i,2] <- snD.vec[i]
    for(j in 3:38){
      # this is for the rest of the compartments (the seropositive compartments)
      Dengvaxia.monthly.mat1[i,j] <- model_output[i,j]*spD.vec[i]*timeline.vector[i]/sum(model_output[i,3:38])
    }
  }
  
  colnames(Dengvaxia.monthly.mat1) <- dimnames(model_output)[[2]][1:38]
  rownames(Dengvaxia.monthly.mat1) <- 1:time_stop
  
  # b. without screening scenario only in some particular scenarios involving Qdenga vaccines
  # creating a matrix of vaccine allocation to corresponding to compartments at each timestep
  # QDENGA
  # first summing up total number of vaccines
  Qdenga.monthly2 <- (snQ.vec+spQ.vec)
  
  # then distributing it according to the population size in each compartment
  Qdenga.monthly.mat2 <- matrix(0, nrow = time_stop, ncol =38)
  for(i in 1:time_stop){
    for(j in 2:38){
      Qdenga.monthly.mat2[i,j] <- model_output[i,j]*Qdenga.monthly2[i]*timeline.vector[i]/sum(model_output[i,2:38])
    }
  }
  
  colnames(Qdenga.monthly.mat2) <- dimnames(model_output)[[2]][1:38]
  rownames(Qdenga.monthly.mat2) <- 1:time_stop
  
  arr <- list("Qdenga with screening"=Qdenga.monthly.mat1, 
              "Qdenga withiout screening"=Qdenga.monthly.mat2,
              "Dengvaxia with screening"=Dengvaxia.monthly.mat1)
  print(screen.num)
  return(arr)
}
