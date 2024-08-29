library(readxl)
library(tidyverse)
library(ggplot2)
library(dplyr)

### 1. AGE STRATIFICATION ####
# loading monthly incidence
data.age <- as.data.frame(read_excel("data.xlsx", sheet = "age-stratified incidence"))

data.pop <- as.data.frame(read_excel("data.xlsx", sheet = "age structure 2022"))
rownames(data.pop) <- data.pop$Subdistrict
data.pop <- data.pop[,-1]

npop.jak <- data.pop["Grand Total", "Grand Total"]
npop.jak0.4 <- data.pop["Grand Total", "0-4"]
npop.jak5.9 <- data.pop["Grand Total", "5-9"]
npop.jak10.14 <- data.pop["Grand Total", "10-14"]
npop.jak15.19 <- data.pop["Grand Total", "15-19"]
npop.jak20.24 <- data.pop["Grand Total", "20-24"]
npop.jak25.29 <- data.pop["Grand Total", "25-29"]
npop.jak30.34 <- data.pop["Grand Total", "30-34"]
npop.jak35.39 <- data.pop["Grand Total", "35-39"]
npop.jak40.44 <- data.pop["Grand Total", "40-44"]
npop.jak45.49 <- data.pop["Grand Total", "45-49"]
npop.jak50.54 <- data.pop["Grand Total", "50-54"]
npop.jak55.59 <- data.pop["Grand Total", "55-59"]
npop.jak60.64 <- data.pop["Grand Total", "60-64"]
npop.jak65.69 <- data.pop["Grand Total", "65-69"]
npop.jak70.74 <- data.pop["Grand Total", "70-74"]
npop.jak75 <- data.pop["Grand Total", "75+"]

npop.jak0.14 <- sum(npop.jak0.4, npop.jak5.9, npop.jak10.14)
npop.jak15.75 <- sum(npop.jak15.19, npop.jak20.24, npop.jak25.29,
                     npop.jak30.34, npop.jak35.39, npop.jak40.44,
                     npop.jak45.49, npop.jak50.54, npop.jak55.59,
                     npop.jak60.64, npop.jak65.69, npop.jak70.74,
                     npop.jak75)

## DATA BY AGE GROUP AND CITY ####
# extracting data based on the age group
a0 <- data.age[1:6,-c(1,2)]
a1.4 <- data.age[7:12,-c(1,2)]
a5.9 <- data.age[13:18,-c(1,2)]
a10.14 <- data.age[19:24,-c(1,2)]
a15.19 <- data.age[25:30,-c(1,2)]
a20.44 <- data.age[31:36,-c(1,2)]
a45.54 <- data.age[37:42,-c(1,2)]
a55.64 <- data.age[43:48,-c(1,2)]
a65.74 <- data.age[49:54,-c(1,2)]
a75  <- data.age[55:60,-c(1,2)]

## MONTHLY INCIDENCE ####
# monthly data of age group A: 0-14 years old
cities <- c("Central Jakarta", "North Jakarta", "West Jakarta", "South Jakarta", "East Jakarta", "Thousands Islands")

a0.4 <- a0 + a1.4
monthly0.14 <- a0.4 + a5.9

# taking the incidence per 100k population
rownames(monthly0.14) <- cities

# monthly data of age group B: 0-14 years old
monthly15.75 <- a15.19 + a20.44 + a45.54 + a55.64 + a65.74 + a75
rownames(monthly15.75) <- cities

# transnsforming monthly DHF incidence data
monthly0.14_tf <- data.frame()
monthly15.75_tf <- data.frame()

start_date <- as.Date("2017-01-01")
end_date <- as.Date("2023-06-01")
monthly_dates <- seq(from = start_date, to = end_date, by = "month")            # Create a sequence of dates from January 2017 to December 2022
cities <- c("Central Jakarta", "North Jakarta", "West Jakarta", "South Jakarta", "East Jakarta", "Thousands Islands")      # Create a sequence of cities

for(i in 1:ncol(monthly0.14)){
  monthly0.14_tf[i, "month"] <- monthly_dates[i]
  monthly0.14_tf[i, "age_group"] <- "0-14"                      # assigning age group
  for(j in cities){
    monthly0.14_tf[i,j] <- monthly0.14[j, i]                    # transforming data [j,i] to [i,j] for each city
  }
  monthly0.14_tf[i,"Jakarta"] <- sum(monthly0.14_tf[i,cities])  # assigning total cases for Jakarta
}
monthly0.14_tf

for(i in 1:ncol(monthly15.75)){
  monthly15.75_tf[i, "month"] <- monthly_dates[i]
  monthly15.75_tf[i, "age_group"] <- "15-75+"                      # assigning age group
  for(j in cities){
    monthly15.75_tf[i,j] <- monthly15.75[j, i]                    # transforming data [j,i] to [i,j] for each city
  }
  monthly15.75_tf[i,"Jakarta"] <- sum(monthly15.75_tf[i,cities])  # assigning total cases for Jakarta
}
monthly15.75_tf

# combining plot of both age groups
monthly0.75_tf <- rbind(monthly0.14_tf, monthly15.75_tf)

## CUMULATIVE INCIDENCE ####
# calculating the cumulative incidence of DHF data for each city
length_data <- ncol(a0)

for(i in 1:6){
  a0.4[i,"total"] <- sum(a1.4[i,1:length_data])
  a5.9[i,"total"] <- sum(a5.9[i,1:length_data])
  a10.14[i,"total"] <- sum(a10.14[i,1:length_data])
  a15.19[i,"total"] <- sum(a15.19[i,1:length_data])
  a20.44[i,"total"] <- sum(a20.44[i,1:length_data])
  a45.54[i,"total"] <- sum(a45.54[i,1:length_data])
  a55.64[i,"total"] <- sum(a55.64[i,1:length_data])
  a65.74[i,"total"] <- sum(a65.74[i,1:length_data])
  a75[i,"total"] <- sum(a75[i,1:length_data])
}

age.group <- c("0-4", "5-9", "10-14", "15-19", 
               "20-24", "25-29", "30-34", "35-39", "40-44",
               "45-49", "50-54",
               "55-64", "60-64",
               "65-69", "70-74",
               "75+")

a.tot <- rbind(a0.4$total, a5.9$total, a10.14$total, a15.19$total, 
               a20.44$total/5, a20.44$total/5, a20.44$total/5, a20.44$total/5, a20.44$total/5, # because the age range fo this age group is 25 years, so each 5-year age group (i.e. 20-24, 25-29, etc.) is interpoalted by taking its 5-year average value 
               a45.54$total/2, a45.54$total/2, # the same approach is used for this age group
               a55.64$total/2, a55.64$total/2, # the same approach is used for this age group
               a65.74$total/2, a65.74$total/2, # the same approach is used for this age group
               rep(a75$total))
colnames(a.tot) <- cities 
rownames(a.tot) <- age.group

a.tot

# creating histogram plot for cumulative incidence by age group
yyy <- as_tibble(as.data.frame(a.tot)) %>% 
  mutate(province = `Central Jakarta` + `North Jakarta` + `West Jakarta` + `South Jakarta` +`East Jakarta` + `Thousands Islands`,
         age_group = as_factor(age.group)) %>%
  pivot_longer(names_to = "cities", cols = c(1:6))

# histogram plot for provincial cumulative incidence
p1 <- ggplot(yyy) +
  aes(x = age_group, fill = province, weight = province) +
  geom_bar() +
  scale_fill_gradient(low = "gray", high = "gray") +
  labs(
    x = "Age Group",
    y = "Cumulative Incidence",
    title = "Cumulative Incidence from January 2017 to December 2022 by Age Group: Jakarta"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# histogram plot for cumulative incidence by cities
p2<- ggplot(yyy) +
  aes(x = age_group, y = value, fill = as_factor(cities)) +
  geom_col() +
  scale_fill_brewer(palette = "Pastel2", direction = 1) +
  labs(x = "Age Group", y = "Cumulative Incidence", fill = "Cities",
       title = "Cumulative Incidence from January 2017 to December 2022 by Age Group: Stratified by Cities") +
  theme_minimal() +
  facet_wrap(vars(as_factor(cities)), scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

gridExtra::grid.arrange(p1,p2)

# plotting monthly incidence
cases.jakarta.plot <- pivot_longer(monthly0.75_tf, col = !c(1, 2)) %>%
  filter((name %in% "Jakarta")) %>%
  ggplot() +
  aes(x = month, y = value, colour = age_group) +
  geom_line() +
  scale_color_brewer(palette = "Set1", direction = 1) +
  labs(x = "Year", y = "Monthly Incidence", color = "Age Group",
       title = "Monthly Incidence from January 2017 to December 2022 by Age Group: Jakarta") +
  theme_minimal() +
  facet_wrap(vars(as_factor(name)), scales = "free_y")

cases.cities.plot <- pivot_longer(monthly0.75_tf, col = !c(1, 2)) %>%
  filter(!(name %in% "Jakarta")) %>%
  ggplot() +
  aes(x = month, y = value, colour = age_group) +
  geom_line() +
  scale_color_brewer(palette = "Set1", direction = 1) +
  labs(x = "Year", y = "Monthly Incidence", color = "Age Group",
       title = "Monthly Incidence from January 2017 to December 2022 by Age Group: Stratified by Cities") +
  theme_minimal() +
  facet_wrap(vars(as_factor(name)), scales = "free_y")

gridExtra::grid.arrange(cases.jakarta.plot, cases.cities.plot)

### 2. DATA SEASONALITY ####
# recalling monthly data for age group A
monthly0.14
# data for jakarta:
N1A <- monthly0.14 %>% colSums() %>% as.numeric()
N1A

# recalling monthly data for age group A
monthly15.75
# data for jakarta:
N1B <- monthly15.75 %>% colSums() %>% as.numeric()
N1B

# total provincial incidence for both age group
N1 <- N1A + N1B

## decompose seasonality ####
# 1. Total cases 
tsN1 <- ts(N1, frequency = 12, start = 2017) # freq 12 => Monthly data. 
decomposedRes1 <- decompose(tsN1, type="multiplicative") # use type = "additive" for additive components
decomposedRes1$seasonal
decomposedRes1$trend
plot(decomposedRes1)

v1 <- decomposedRes1$seasonal
v1
par(mfrow=c(1,1))
plot(v1, type="l", main = "seasonality x trend")
# since the first 6 months of 2017 trend cannot be extracted, 
# only the seasonality and trend from Jan 2018 to Dec 2023 will be taken
seas.vec0 <- v1[c(13:(72))]

# Assuming the cycle occurs every 5 years, and the run-time for the model is 15 years
# The seasonality vector is replicated a few times
seas.vec <- rep(as.numeric(seas.vec0), 4)

par(mfrow = c(1,2))
plot(seas.vec0, type = "l", main = "5-year seasonality cycle")
plot(seas.vec, type = "l", main = "50-year seasonality cycle ")

# 2. Age group A
tsN1A <- ts(N1A, frequency = 12, start = 2017) # freq 12 => Monthly data. 
decomposedRes1A <- decompose(tsN1A, type="multiplicative") # use type = "additive" for additive components
decomposedRes1A$type <- "multiplicative of group A"
decomposedRes1A$seasonal
decomposedRes1A$trend
plot(decomposedRes1A)

# since the first 6 months of 2017 trend cannot be extracted, 
# only the seasonality and trend from Jan 2018 to Dec 2023 will be taken
seas.vec0A <- decomposedRes1A$seasonal[13:72]*(0.8 + 0.2*decomposedRes1A$trend[13:72]/sum(decomposedRes1A$trend[13:72])*60)

# Assuming the cycle occurs every 5 years, and the run-time for the model is 15 years
# The seasonality vector is replicated a few times
seas.vecA <- rep(as.numeric(seas.vec0A), 4)

# 3. Age group B
tsN1B <- ts(N1B, frequency = 12, start = 2017) # freq 12 => Monthly data. 
decomposedRes1B <- decompose(tsN1B, type="multiplicative") # use type = "additive" for additive components
decomposedRes1B$type <- "multiplicative of group B"
decomposedRes1B$seasonal
decomposedRes1B$trend
plot(decomposedRes1B)

# since the first 6 months of 2017 trend cannot be extracted, 
# only the seasonality and trend from Jan 2018 to Dec 2023 will be taken
seas.vec0B <- decomposedRes1B$seasonal[13:72]*(0.2 + 0.2*decomposedRes1B$trend[13:72]/sum(decomposedRes1B$trend[13:72])*60)

# Assuming the cycle occurs every 5 years, and the run-time for the model is 15 years
# The seasonality vector is replicated a few times
seas.vecB <- rep(as.numeric(seas.vec0B), 4)

## AVERAGE INCIDENCE BY MONTH ####
# assigning data vector of provincial monthly incidence from 2018 to 2022
casesA <- N1A[1:72]
casesB <- N1B[1:72]

avg.A <- c()
for(i in 1:12){
  avg.A[i] <- mean(c(casesA[i], casesA[i+12], casesA[i+24], casesA[i+36], casesA[i+48], casesA[i+60]))
}
avg.A <- rep(avg.A)

avg.B <- c()
for(i in 1:12){
  avg.B[i] <- mean(c(casesB[i], casesB[i+12], casesB[i+24], casesB[i+36], casesB[i+48], casesB[i+60]))
}
avg.B <- rep(avg.B)

## AVERAGE INCIDENCE BY MONTH PER 100,000 PEOPLE ####
avg.A.per.100000 <- avg.A*100000/npop.jak0.14
avg.B.per.100000 <- avg.B*100000/npop.jak15.75