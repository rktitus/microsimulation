################################################################################
################# MALARIA MICROSIMULATION MODEL  ##### 2021 ####################
################################################################################
#
# This is a microsimulation model for healthcare cost estimation. 
#
# It estimates cost of treating malaria in Kenya.
#
# It was developed as part of the analysis for chapter two of a PhD thesis on
# healthcare costing in Kenya.
#
# It was developed and ran on the following:
# - R version 4.0.3 &
# - RStudio version 1.4.1106
#
# Although extensive use of packages (except for graphical representation of 
# results) have been avoided wherever possible, it cannot be guaranteed that
# it will run on earlier versions. Data.tables have been used to speed up the
# modeling/simulations.
#
# This model has 4 main sections:
# - Settings which prepares the R environment for modeling
# - The inputs which sets the input data and variables for the model
#   - This part can be modified at will by the user
# - The model section which
#   - estimates the resistance to ACT (AL)
#   - uses the supplied data inputs to calculate other model inputs and/or
#   - the microsimulation model which runs the actual microsimulation
# - The output section which organizes the output of the microsimulation model
#   in a clear and presentable format.

#############################  SETTINGS  #######################################

rm(list = ls())                   # remove any variables in R's memory

options(scipen = 999, digits=5)   # Disabling the scientific notation in R

# Installing/Loading required packages (can omit rstudioapi and set wd manually)
if (!require(data.table)) install.packages('data.table') ; library(data.table)
if (!require(rstudioapi)) install.packages('rstudioapi') ; library(rstudioapi)
if (!require(tidyverse)) install.packages('tidyverse') ; library(tidyverse)
if (!require(reshape2)) install.packages('reshape2') ; library(reshape2)
if (!require(pryr)) install.packages('pryr') ; library(pryr)
if (!require(rms)) install.packages('rms') ; library(rms)
if (!require(deSolve)) install.packages('deSolve') ; library(deSolve)
#if (!require(fitdistrplus)) install.packages('fitdistrplus') ; library(fitdistrplus)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))

# set.seed(1)       # Make the model output reproducible

#Creating new environment (for functions)

my_env <- new.env(parent = globalenv())

#Loading functions to R

path_dir <- list.files(pattern = "[.]R$", path = "./functions/", full.names = T)
sapply(path_dir, FUN=source)

# Increase Memory Limit (3 times)

memory.limit(memory.limit()*3)

##################################### INPUTS ###################################

########################### Model Data #########################################

# Mortality, incidence, fertility, population proportions, sex ratios and
# real (inflation/interest) rates

load("./Data/Out/Out.Data")

######################### Input Parameters#####################################

start_year <- y_0 <- 2010 -1              # -1 to allow for pregnancy
n_days <- 364                             # Number of days in a year

#Check if this is a valid year based on the input data
ifelse(start_year %in% unique(mort$year),"VALID start_year",
stop("Invalid start_year: Try a different start_year"))

n_i <- n_0 <- 10^5*1                      # Number of individuals
n_t <- 2                                  # FOllow-up period (in years)

# Check if the followup period is valid based on the input data
ifelse((start_year+n_t) %in% unique(mort$year),"Valid n_t",
       stop("Invalid n_t: Try a different n_t"))

step_size = 14                          # Modelling block/step size (in days)

n_b <- (n_t+1)/(step_size/n_days)     # Total blocks of follow-up period
Ages <- unique(mort$age)                # Ages to be considered

# Get the age distribution for each gender
age_distr <- setkey(pop_prop, age, year)[.(Ages)][year == start_year,]

# Proportionately sample from each gender the ages to include based on the
# respective age distribution.
male_ages <- sample(Ages, replace = T, 
                    round(n_i*(sex_ratio[year == start_year & sex == "M", prop])), 
                    prob = age_distr[year == start_year & sex == "M", prop])
female_ages <- sample(Ages, replace = T, 
                    round(n_i*(sex_ratio[year == start_year & sex == "F", prop])), 
                    prob = age_distr[year == start_year & sex == "F", prop])

ages_com <- rbind(data.table(age = male_ages, sex = "M"), 
                   data.table(age = female_ages, sex = "F"))

# Draw sample to use in the microsimulation
ages_comb <- ages_fin <- cbind(ID = 1:n_i, ages_com[sample(nrow(ages_com)), ])

##################### age specific probabilities  #############################

# Under 5
P1_U_f <- 0.346; P1_U_m <- 0.375; P4_U_f <- 0.722; P4_U_m <- 0.716
P5_U <- 0.383; P35_U = 0.029; P37_U = 0.5; P38_U = 0.0132; P39_U = 0.192
P49_U <- 0.2

# Over 5

P1_O <- 0.31; P4_O <- 0.64; P5_O <- 0.68; P35_5_14 <- rbeta(1, 69.799, 11398.63)
P35_15_49 <- rbeta(1, 1, 12585.95); P35_50plus <- rbeta(1, 1, 4744.006)
P37_O <- 0.25; P38_O <- 0.005; P39_O <- 0.1; P49_O <- runif(1, 0.1, 0.25)

####################### non age specific probabilities #######################

P2 <- 0.711; P8 <- 0.929; P9 <- 0.912; P10 <- 0.545; P11 <- 0.1; P18 <- 0.367; 
P19 <- 0.25; P20 <- 0.129; P21 <- 0.37; P22 <- 0.325; P23 <- 0.864; P24 <- 0.442; 
P25 <- 0.697; P26 <- 0.2; P27 <- 0.9065; P28 <- 0.872; P29 <- 0.953; P30 <- 0.973; 
P31 <- 0.9204; P32 <- 0.9733; P33 <- 0.667; P34 <- 0.97; P36 <- 0.48; 
P40 <- runif(1, 0.4, 0.6); P41 <- 0.5; P43 <- 0.3; P44 <- 5; P45 <- 0.17; 
P46 <- runif(1, 0.5, 0.75); P50 <- 0.2; P51 <- runif(1, 0.991,1); P53 <- 0.1

####################### Derived Probability Inputs ##########################

P6 <- P23*(1-((1-P24)*(1-P25))); P14 <- P8 * (1-P33); P15 <- P9*(1-P33)
P16 <- P10 * (1-P33)*(1-P8); P17 <- P10 * (1-P33)*(1-P9)

#################### Resistance Model #######################################

# Parameters input

params <- list(
  lambda = 0.025,    #Population growth rate pa
  delta = 1/65,  #Overall annual mortality rate
  gamma = 0.1,   #Rate at which immunes become  susceptible again pa
  f = mean(c(P4_U_f, P4_U_m, P4_O)), #Fraction of infected population that receive treatment
  theta_s = 0.7, #Recovery rate of susceptible individuals to join immune pa
  theta_r = 2.5, #Recovery rate of individuals with resistant strain to join immune pa
  kappa = 12, #Excess rate of recovery of treated individuals pa
  m= 4,      #Number of female mosquitoes per human (per night)
  alpha = 537/n_days,   #Biting rate (per night)
  beta_1 = 0.8,     #Infectiousness of humans to mosquitoes
  beta_2 = 0.8,     #Susceptibility of mosquitoes to humans
  mu = 36/n_days,      #Mortality rate of mosquitoes (Daily mortality of mosquitoes)
  tau = 8.75,     #Incubation period of parasites in the mosquito
  step_size, #Modelling block/step size (in days)
  y = 0.01, #initial malaria prevalence
  r = 10^{-5},  #initial drug resistance
  z = 0,   #proportion initially immune
  start_year = 2010,  # The year antimalaria (here ACT/AL) was first used
  years = 5,  #Number of burn years. May want to choose this value by running 
              #resistance_0() with different values
  years_fit = n_t+2 #Number of years to simulate spread of resistance. Add
                    # 2 to allow for final year + initial 19 months
)

# Estimating Resistance

Rt <- resistance(parameters = params)$out2[,"Y_r"]/
      (resistance(parameters = params)$out2[,"Y_r"] + 
      resistance(parameters = params)$out2[,"Y_s"])

######################## Additional probabilities ########################

pmp   <- 0.1      # Increased susceptibility to malaria due to pregnancy
qp    <- 0.1      # Increased mortality due to pregnancy
qNS   <- 0.1      # Increased mortality due to NS
r_imm <- (1-exp(-params$theta_s*step_size/n_days))   # Prob of immunity gain

######################## Cost Inputs #####################################

# Table 4.4
Costs <- data.table(Cost = c("C1", "C2", "C3", "C4", "C8", "C9", "C20", "C21"),
                B.E = c(0, 0.74, 3.9, 0.72, 0.56, 0.99, 3.84, 14.15),
                Distr = c(NA, rlnorm(1, -0.357, 0.335), rlnorm(1, 1.286, 0.338),
                          rlnorm(1, -0.351, 0.213), rlnorm(1, 4.025, 0.027)/100,
                          rlnorm(1, 4.595, 0.034)/100, rlnorm(1, 1.504, 0.3),
                          rlnorm(1, 2.387, 0.293)),
                Base = c(NA, 1995, 2002, 2002, 2021, 2021, 2004, 2004))

####################### The Microsimulation Model #######################

# In this microsimulation model, we're interested in 5 costs:
#     - The cost of false positive malaria
#     - The costs for hospital and health center visits
#     - The cost of delayed malaria treatment
#     - The cost following drug failure and
#     - The overall (total) cost for treating malaria in Kenya.
#
# We'd therefore want to track the "health states visited by an individual. 
# To do this, for each individual at each simulation, we store:
#     - The cost of accessing healthcare - Cost
#     - Healthcare access path history (where they've visited) - Trace
#     - The health state after each simulation - State_Matrix
#     - Those recovering with NS - NS

Cost <- array(0,dim = c(n_i, n_b),dimnames = list(paste("ID",   1:n_i, sep =" "),
                                         paste("cycle", 1:n_b, 
                                               sep =" ")))
Trace <- State_Matrix <- array(NA_character_,
                            dim = c(n_i, n_b),
                            dimnames = list(paste("ID",   1:n_i, sep =" "),
                                            paste("cycle", 1:n_b, 
                                                  sep =" ")))

Cost %>% data.frame() %>% setDT(.)
Trace %>% data.frame() %>% setDT(.)

NS <- 0

##### Initial filters for incorporating pregnancy in the model #########
setkey(ages_comb, age, sex)
ferty <- cbind(fert, sex = "F")

#Filtering only those who can get pregnant - start off as healthy
can_preg <- merge(ages_comb,ferty[year == start_year])[,
            rate:=1-exp(log(1-rate)*step_size/n_days)][,.(ID, age, sex,rate)]
preg <- data.table(ID = can_preg$ID, tp = 0)  # tp - time in pregancy
preg_st <- data.table(ID = can_preg$ID, state="H")
#can_preg_1 <- unique(can_preg$ID)    # Use to track/add those newly 
                                      # Entering the pregnancy age range
change_year <- 0                      # Track when year changes
preg_55 <- 0                          # Track those pregant at age 55
n_c_babies <- c(0)  # number of babies in each cycle
baby_ind <- n_i # Index for adding new babies to the sample
die <- ages_comb[age > 85,]     # track those dying, including at age 85.

# Track the states over time
state_track <- data.table(ages_comb[,c(1:3)], current = "Healthy")

inci <- inci[age %chin% Ages,]

################################## THE MODEL ##########################

start_time <- Sys.time()

for (i in 1:n_b)
{
  #Update the eligibility criteria - remove those who die(d) and above 85.
  # First store those dying, including at 85
  die <- rbind(die, ages_comb[ID %chin% 
                              state_track[current == "Die",]$ID,])
  ages_comb <- ages_comb[ID %chin% state_track[current != "Die",]$ID,]
  preg_st <- preg_st[ID %chin% state_track[current != "Die",]$ID,]
  preg <- preg[ID %chin% state_track[current != "Die",]$ID,]

  ########################## Evaluation of Pregnancy ##################
  
  # Allowing for pregnancy
  can_preg <- merge(ages_comb,ferty[year == start_year])[,
              rate:=1-exp(log(1-rate)*step_size/n_days)][,
              .(ID, age, sex,rate)]
  
  # Add new entries to the pregnancy state as Healthy (and not pregnant)
  # when the year changes
  if(change_year>0){
    # Adding those aged 10
    new_preg <- data.table(ID = can_preg[age == 10,]$ID, tp = 0)
    new_preg_st <- data.table(ID = can_preg[age == 10,]$ID, state="H")
    preg <- rbind(new_preg,
                  preg[ID %chin% ages_comb[ID %in% preg$ID & age < 55,]$ID, ], 
                  preg[ID %chin% preg_55,])
    preg_st <- rbind(new_preg_st,
                     preg_st[ID %chin% ages_comb[ID %in% preg$ID & 
                                            age < 55,]$ID, ],
                     preg_st[ID %chin% preg_55,])
    change_year <- 0
  }
  
  # Track time during pregnancy (increments of 2 wks)
  preg[ID %in% preg_st[state=="P",]$ID,]$tp <- 
    preg[ID %in% preg_st[state=="P",]$ID,]$tp + 1
  
  # Transition matrix from Healthy (not pregnant) to pregnant only
  # for those not yet pregnant
  preg_prob <- can_preg[,`:=`(P=rate, H=1-rate)][
                ID %chin% preg[preg$tp==0,]$ID,.(P,H)]
  
  # Sampling the next state (Healthy/pregnant)
  new_st <- my_env$samplev(preg_prob, 1)
  preg_st[ID %chin% preg[preg$tp==0,]$ID, state := new_st] 

  #Adding the babies to the sample data
  n_babies <- length(preg[preg$tp>19,]$ID)
  if (n_babies > 0) {
    new_babies <- data.table(ID = c(baby_ind+1):c(baby_ind + n_babies),
                             age = rep(0, n_babies),
                             sex = sample(c("F", "M"), prob = c(0.5, 0.5),
                                          replace = T, size = n_babies))
    ages_comb <- rbind(ages_comb, new_babies)
    
    # Adding the babies to the state and cost tables
    n_c_babies <- c(n_c_babies, n_babies)
    Cost <- rbind(Cost, matrix(0, nrow=n_babies, ncol=n_b))
    Trace <- rbind(Trace, matrix(NA_character_, nrow=n_babies, ncol=n_b))
    State_Matrix <- rbind(State_Matrix, matrix(NA_character_, nrow=n_babies, 
                                               ncol=n_b))
    state_track <- rbind(state_track, data.table(new_babies[,c(1:3)], 
                                                 current = "Healthy"))
    
    # Updating the sample size
    n_i <- nrow(Cost)
    baby_ind <- baby_ind + n_babies
  }
  
  # Allowing for delivery on the 9th month = 39 wks - step size = 2wks.
  # Return those deliverying to healthy state
  preg_st[ID %in% subset(preg, tp>19)$ID,]$state <- "H"
  preg <- subset(preg, tp<=19)    # Remove those deliverying

  #Drop those giving birth at age 55
  preg[-c(ID %in% preg[c(tp == 0 & ID %in% preg_55),]$ID),]
  
  ######################## Evaluating the Decision Treee #######################
  
  # Updating state history based on the new sample
  state_track <- state_track[state_track$ID %in% ages_comb$ID,]
  
  # Defining filters to be used with age/gender specific probabilities
  under_5_female <- ages_comb[ages_comb$age<=5 & ages_comb$sex=="F",]
  under_5_male <- ages_comb[ages_comb$age<=5 & ages_comb$sex=="M",]
  under_5 <- rbind(under_5_female, under_5_male)
  over_5 <- ages_comb[ages_comb$age>5,]
  age_5_14 <- ages_comb[ages_comb$age>5 & ages_comb$age<15,]
  age_15_49 <- ages_comb[ages_comb$age>=15 & ages_comb$age<50,]
  age_50plus <- ages_comb[ages_comb$age>=50,]
  
  ############# Estimating age-dependent probabilities ###############
  P1 <- rbind(cbind(under_5_female, P1 = P1_U_f),cbind(under_5_male, P1 = P1_U_m),
              cbind(over_5, P1 = P1_O))
  
  P3 <- merge(ages_comb, inci[year == start_year])[,P3 := (P2*pm)/P1$P1][,
                                                          .(ID, age, sex, P3)]
  
  P4 <- rbind(cbind(under_5_female, P4 = P4_U_f), cbind(under_5_male, P4 = P4_U_m),
              cbind(over_5, P4 = P4_O))
  
  P5 <- rbind(cbind(rbind(under_5_female, under_5_male), P5 = P5_U),
              cbind(over_5, P5 = P5_O))
  
  P7 <- merge(ages_comb, mort[year == start_year])[,P7 := 1-(1-q)^((q-qm)/q)][,
                                                  .(ID, age, sex, P7)]
  
  q <- merge(ages_comb, mort[year == start_year])[,.(ID, age, sex, q)]
  
  P35 <- rbind(cbind(under_5, P35 = P35_U),
               cbind(age_5_14, P35 = P35_5_14),
               cbind(age_15_49, P35 = P35_15_49),
               cbind(age_50plus, P35 = P35_50plus))
  
  P37 <- rbind(cbind(under_5, P37 = P37_U),cbind(over_5, P37 = P37_O))
  
  P38 <- rbind(cbind(under_5, P38 = P38_U),cbind(over_5, P38 = P38_O))
  
  P39 <- rbind(cbind(under_5, P39 = P39_U),cbind(over_5, P39 = P39_O))
  
  P42 <- cbind(P37[,P37. := (1-P40)*P37*P41],P39[,"P39"])[, 
                                P42:=P39 + P37.][, -c(4:6)]
  
  P47 <- rbind(cbind(under_5, P47 = 0.579),cbind(over_5, P47 = 0.421))
  
  P48 <- rbind(cbind(under_5, P48 = runif(1,0.1,0.25)),cbind(over_5, P48 = 0.2))
  P49 <- rbind(cbind(under_5, P49 = P49_U), cbind(over_5, P49 = P49_O))
  
  P52 <- rbind(cbind(under_5, P52 = 0.0132),cbind(over_5, P52 = 0.005))
  P53 <- P18 * (1- (P23*(P25*P27*(1-P24) + P28*P24*(1-P25) + 0.5*P25*P24*
                      (1-(1-P27)*(1-P28))))) + (1-P10*(1-P8))* (P23*(P25*P27*
                    (1-P24) + P28*P24*(1-P25) + 0.5*P25*P24*(1-(1-P27)*(1-P28))))
  P54 <- P10*(1-P8)*(1- (P23*(P25*P27*(1-P24) + P28*P24*(1-P25) + 0.5*P25*P24*
                     (1-P27)*(1-P28))))
  P55 <- P18*(1-P34) + P34*(1-P10*(1-P8))
  P56 <- P34 * P10 *(1-P8)
  P57 <- cbind(P47[,.(ID, age, sex)], P57 = 1-((1-Rt[i+1])*((1-P47$P47)*P48$P48 + 
                                                         P47$P47)))
  P58 <- cbind(P49[,.(ID, age, sex)], P58 = 1-(P51*(P49$P49 + P50*(1-P49$P49))))
  P59 <- P18 * (1- (P23*(P25*P27*P31*(1-P24) + P28*P24*(1-P25) + 0.5*P25*P24*
                           (1-(1-P27)*(1-P28))))) + (P23*(P25*P27*P31*(1-P24) + 
                            P28*P24*(1-P25) + 0.5*P25*P24*
                            (1-(1-P27)*(1-P28))))*(1-P10*(1-P9))
  P60 <- (P23*(P25*P27*P31*(1-P24) + P28*P24*(1-P25) + 0.5*P25*P24*
                 (1-(1-P27)*(1-P28)))) * P10 * (1-P9)
  P61 <- P18*(1-P34) + P34*(1-P10*(1-P9))
  P62 <- P34*P10*(1-P9)
    
  ################### Drug Costs ###################
  
  # ACT Costs: Artemether + Lumefatrine (Table 4.1)
  AL <- rbind(cbind(ages_comb[ages_comb$age<=3], price = 0.5810),
              cbind(ages_comb[ages_comb$age>3 & ages_comb$age<=7], price = 1.1619),
              cbind(ages_comb[ages_comb$age>7 ], price = 0.9212))
  
  # Quinine COsts: Quinine + Dextrose cost (Table 4.1)
  QN <- rbind(cbind(ages_comb[ages_comb$age<=3], price = runif(1,0.0777, 0.3177)),
              cbind(ages_comb[ages_comb$age>3 & ages_comb$age<=7], 
                    price = runif(1,.3404, 0.5445)),
              cbind(ages_comb[ages_comb$age>7 & ages_comb$age<=11], 
                    price = runif(1,.6372, 0.7642)),
              cbind(ages_comb[ages_comb$age>11 & ages_comb$age<=15], 
                    price = runif(1,.7762, 0.8844)),
              cbind(ages_comb[ages_comb$age>15 & ages_comb$age<=19], 
                    price = runif(1,.8964, 1.0646)),
              cbind(ages_comb[ages_comb$age>19 ], price = 1.0766))
  
  # Ceftriaxone (Table 4.1)
  Cef <- 0.792
  
  # DHA-PPQ Costs (Table 4.1)
  DHA <- rbind(cbind(ages_comb[ages_comb$age<=3], price = 0.9282),
               cbind(ages_comb[ages_comb$age>3 & ages_comb$age<=5], 
                     price = 1.8565),
               cbind(ages_comb[ages_comb$age>5 & ages_comb$age<=11], 
                     price = 0.8116),
               cbind(ages_comb[ages_comb$age>11 & ages_comb$age<=16], 
                     price = 1.6233),
               cbind(ages_comb[ages_comb$age>16 ], price = 2.0174))
  
  ###################### Estimating real costs ################
  
 # if (i == 1 |i%%(n_days/step_size) == 0) {
    # Only evaluate these for first cycle, and new year(s)
  C1 <- Costs[Cost == "C1",B.E]
  C2 <- Costs[Cost == "C2",B.E] * real[year == start_year, real]/
    real[year == (Costs[Cost == "C2",Base]-1), real]
  C3 <- Costs[Cost == "C3",B.E] * real[year == start_year, real]/
    real[year == (Costs[Cost == "C3",Base]-1), real]
  C4 <- Costs[Cost == "C4",B.E] * real[year == start_year, real]/
    real[year == (Costs[Cost == "C4",Base]-1), real]
  C8 <- Costs[Cost == "C8",B.E] * real[year == start_year, real]/
    real[year == (Costs[Cost == "C8",Base]-1), real]
  C9 <- Costs[Cost == "C9",B.E] * real[year == start_year, real]/
    real[year == (Costs[Cost == "C9",Base]-1), real]
  C20 <- Costs[Cost == "C20",B.E] * real[year == start_year, real]/
    real[year == (Costs[Cost == "C20",Base]-1), real]
  C21 <- Costs[Cost == "C21",B.E] * real[year == start_year, real]/
    real[year == (Costs[Cost == "C21",Base]-1), real]
  
  C10 <- AL[,.price := price*real[year == start_year, real]/
              real[year == 2019, real]][,.(ID, age, sex, .price)]
  C17 <- QN[,.price := price*real[year == start_year, real]/
              real[year == 2019, real]][,.(ID, age, sex, .price)]
  C23 <- Cef*real[year == start_year, real]/real[year == 2019, real]
  C26 <- DHA[,.price := price*real[year == start_year, real]/
              real[year == 2019, real]][,.(ID, age, sex, .price)]
  
  ####################### Derived Cost Inputs ##################################
  
  C5 <- C2 + C3*((1-P21)*(1-P22)); C6 <- C2 + C4*((1-P21)*(1-P22))
  C7 <- C8*P24*(1-P25) + (C9*P25*(1-P24)) + (P24*P25*0.5*(C8 + C9))
  C11 <- C10[, C11 := .price * (1+P26)*P8 + (.price*(1+P26)*P8*
        (P25*((1-P29) + (P29*(1-P32))) + (P24*(1-P30))))][,.(ID, age, sex, C11)]
  C12 <- C10[, C12 := .price * P18 * (1 + P26) * (1 + P19)][,.(ID, age, sex, C12)]
  C13 <- C10[, C13 := P18 * .price * (1+ P26) * 
               (1+P19) * (1+P11)][,.(ID, age, sex, C13)]
  C14 <- C17[, C14 := .price * (1+P26)][,.(ID, age, sex, C14)]
  C15 <- C10[, C15 := .price * (1+P26) * (1+ ((1+P19) * (1-P8)))][,
                                                          .(ID, age, sex, C15)]
  C16 <- C10[, C16 := .price * (1+P25)*P9 + (.price*(1+P25)*P9*
        (P25*((1-P29) + (P29*(1-P32))) + (P24*(1-P30))))][,.(ID, age, sex, C16)]
  C18 <- C10[, C18 := .price * (1+P26) * (1+ ((1+P19) * (1-P9)))][,
                                                          .(ID, age, sex, C18)]
  C19 <- C10[, C19:= .price * P18*(1+P26)*(1+P19)*(1+P11)][,.(ID, age, sex, C19)]
  C22 <- C10[, C22:= .price*P9*(1+P26) + C6*(1-P40) + C20 + (P44*C21*(1-P45)*(1-P46)) +
               P43*(C23*(1+P26))+C17[, C22:= .price*P10*(1+P26)*(1-P40)]$C22][,
                                                          .(ID, age, sex, C22)]
  C24 <- C10[, C24:= .price*P8*(1+P26) + C20 + (P44*C21*(1-P45)*(1-P46)) +
               P43*(C23*(1+P26))][, .(ID, age, sex, C24)]
  C25 <- C26[, C25:= C5 + C8 + .price*P47$P47*(1+P26) + C10[,
          x:= .price*(1+P26)*(1-P47$P47)]$x][, .(ID, age, sex, C25)]
  C27 <- C10[, C27:= .price*(1+P26)*(1 + P18*(1 + P19))][, .(ID, age, sex, C27)]
  
 # }

  
  #################### Malaria Microsimulation #################################
  
  ############################### Entry Point ##################################
  
  state_track <- my_env$choose_alt3(states = c("Fever", "Die", "Healthy"), 
                                     prob = list(1-exp(log(1-P1$P1)*step_size/n_days), 
                                                 1-exp(log(1-q$q)*step_size/n_days)), 
                                     track = state_track, 
                                     current_node = "Healthy",
                                     pmp = pmp,qp = qp, qNS = qNS)
  
  # Update for those already sick
  state_track <- my_env$choose_alt3(states = c("Fever", "Die", "Sick"), 
                                     prob = list(1-exp(log(1-P1$P1)*step_size/n_days), 
                                                 1-exp(log(1-q$q)*step_size/n_days)), 
                                     track = state_track, 
                                     current_node = "Sick",
                                     pmp = pmp,qp = qp, qNS = qNS)
  
  # Evaluating the Fever node
  state_track <- my_env$choose_alt(states = c("No Malaria", "Malaria"), 
                                   prob = P3$P3, 
                                   track = state_track, 
                                   current_node = "Fever")
  Trace[state_track[current == "Malaria",ID],i] <- "Malaria"
   
  ##################### Malaria Negative Febrile Cases #########################
  
   my_env$no_Mal()
  
  ##################### Malaria Positive Febrile Cases #########################
  
  ########### Malaria-Positive Non-Care Seeking cases ##########################
  
  state_track <- my_env$choose_alt(states = c("No Care", "Care"), 
                                   prob = P4$P4, 
                                   track = state_track, 
                                   current_node = "Malaria")
  
  Trace[state_track[current == "No Care",ID],i] <- "MNC"
  
  state_track <- my_env$choose_alt3(states = c( "Severe Malaria", "Die", "Recover"), 
                                    prob = list(1-exp(-P35$P35),
                                                1-exp(log(1-q$q)*step_size/n_days)), 
                                    track = state_track, 
                                    current_node = "No Care", pmp = pmp,
                                    qp = qp, qNS = qNS)
  
  State_Matrix[state_track[current == "Die",ID],i] <- "D"
  State_Matrix[state_track[current == "Recover",ID],i] <- "UM"
  
  # All recovered enter the Sick state
  state_track[state_track$current=="Recover","current"] <- "Sick"
  
  my_env$severe_malaria(trace = c("MNCSMNC", "MNCSMH", "MNCSMHC"))
  
  ############### Malaria-Positive Care Seeking cases ##########################
  
  state_track <- my_env$choose_alt(states = c("Hospital", "Health Center"), 
                                    prob = P5$P5, 
                                    track = state_track, 
                                    current_node = "Care")
  
  #####################  Hospital Evaluation ###################################
  
  my_env$Mal_H()
  
  #####################  Health Center Evaluation ##############################
  
  out <- my_env$Mal_HC()
  
  # state_track <- my_env$choose_alt(states = c("No Diagnosis", "Diagnosis"), 
  #                                   prob = P6, 
  #                                   track = state_track, 
  #                                   current_node = "Hospital")
  # 
  # Trace[state_track[current == "Diagnosis",ID],i] <- "M_HD"
  # Trace[state_track[current == "No Diagnosis",ID],i] <- "M_HND"
  # 
  # state_track <- my_env$choose_alt3(states = c("Quinine", "ACT", "No Treatment"), 
  #                                    prob = list(P16, P8*(1-P33)), 
  #                                    track = state_track, 
  #                                    current_node = "Diagnosis")
  # 
  # out <- my_env$failure() 
  
  Cost <- out[[1]]; state_track <- out[[2]]
  State_Matrix <- out[[3]]; Trace <- out[[4]]
  
  # Moving to the next year
  if (i%%(n_days/step_size) == 0) {
    start_year = start_year + 1           # Moving to the next year
    change_year <- change_year + 1        # Track the year change
    ages_comb$age = ages_comb$age + 1     # Increasing age by 1 every year
    can_preg$age <- can_preg$age + 1
    die$age <- die$age + 1
    state_track$age <- state_track$age + 1
    state_track[age>85,]$current <- "Die" # No one lives beyond 85.
    
    # Updating eligibility of getting pregnant + keeping those already pregnant
    # can_preg_1 <- can_preg[ID %chin% state_track[current != "Die",]$ID,][
    #            age <= max(ferty[year == start_year]$age) |
    #              ID %in% preg[preg$tp>0,]$ID,]$ID
    
    # Those aged 55 now and are pregnant
    preg_55 <- can_preg[ID %chin% state_track[current != "Die",]$ID,][
      age ==55, ][ID %in% preg[preg$tp>0,]$ID,]$ID
  } 
  cat('\r', paste(round(i/n_b * 100, 0) , "% done", ": current year =",
                  start_year, sep = " "))
  
}


end_time <- Sys.time()
run_time <- end_time - start_time
run_time



# Reformatting State_Matrix for those who die
for (j in 1:(n_b-1))
{
  State_Matrix[State_Matrix[,j] %like% "D", j+1] = "D"
}

# Combine all subjects (died + alive) used in the microsimulation

age_filter <- rbind(ages_comb, die)
age_filter <- age_filter[ order(ID), ]

# Checking missing cost values

colSums(is.na(Cost)) %>% as.data.frame() %>% unique()

# Creating yearly cost for each of the simulation years

n_w <- n_days/step_size

yearly_cost <- Cost %>% 
  cbind(age_filter[,c("ID", "age")], rowSums(.[,1:n_w])) %>%
  filter(age>=(n_t+1)) %>% select((n_b+1):(n_b+3)) %>%
  mutate(age = age-n_t-1)

for (i in 1:n_t)
{ 
  yearly_cost <- Cost %>%
    cbind(age_filter[,c("ID", "age")], rowSums(.[,((n_w*i+1):(n_w*(i+1)))])) %>%
    filter(age>=(n_t+1)) %>% select((n_b+1):(n_b+3)) %>%
    mutate(age = age-n_t+i-1) %>% select(3) %>% bind_cols(yearly_cost,.) 
}

colnames(yearly_cost) <- c("ID", "age", as.character(y_0: (y_0+n_t)))

# Creating age-wise cost in each of the simulation years

age_specific_cost <- list()
age_fil <- age_filter[ID %chin% 1:(n_0 + sum(n_c_babies[1:(length(n_c_babies)
                                            %%(n_days/step_size))])),]
age_specific_cost[[1]] <-  Cost %>% 
  cbind(age_filter[,c("ID", "age")], rowSums(.[,1:n_w])) %>%
  select((n_b+1):(n_b+3)) %>% .[!(State_Matrix[,n_w] %like% "D"),] %>%
  .[ID %chin% age_fil$ID] %>% mutate(age = age-n_t-1) %>% .[age %chin% Ages]
for (i in 1:n_t)
{
  age_fil <- age_filter[ID %chin% 1:(n_0 + 
              sum(n_c_babies[1:(length(n_c_babies)%%
                                  (n_days/step_size) + (n_w*i))])),]
  age_specific_cost [[i+1]] <- Cost %>%
    cbind(age_filter[,c("ID", "age")], rowSums(.[,((n_w*i+1):(n_w*(i+1)))])) %>%
    select((n_b+1):(n_b+3)) %>% .[!(State_Matrix[,(n_w*i)] %like% "D"),] %>% 
    .[ID %chin% age_fil$ID] %>% mutate(age = age-n_t+i-1) %>% 
    .[age %chin% Ages]
}

# age-wise cost for patients experiencing treatment failure

Trt_F %<d-% matrix(0, nrow = n_i, ncol = n_b)
for (i in 1:n_b)
{
  Trt_F[grep("F", Trace[,i]),i] <- as.numeric(Cost[grep("F", Trace[,i]),i]) 
}

age_specific_Trt_F_cost <- list()
age_fil <- age_filter[ID %chin% 1:(n_0 + sum(n_c_babies[1:(length(n_c_babies)
                                                           %%(n_days/step_size))])),]
age_specific_Trt_F_cost[[1]] <-  Trt_F %>% 
  cbind(age_filter[,c("ID", "age")], rowSums(.[,1:n_w])) %>%
  select((n_b+1):(n_b+3)) %>% .[!(State_Matrix[,n_w] %like% "D"),] %>%
  .[ID %chin% age_fil$ID] %>% mutate(age = age-n_t-1) %>% .[age %chin% Ages]
for (i in 1:n_t)
{
  age_fil <- age_filter[ID %chin% 1:(n_0 + 
                                       sum(n_c_babies[1:(length(n_c_babies)%%
                                                           (n_days/step_size) + (n_w*i))])),]
  age_specific_Trt_F_cost [[i+1]] <- Trt_F %>%
    cbind(age_filter[,c("ID", "age")], rowSums(.[,((n_w*i+1):(n_w*(i+1)))])) %>%
    select((n_b+1):(n_b+3)) %>% .[!(State_Matrix[,(n_w*i)] %like% "D"),] %>% 
    .[ID %chin% age_fil$ID] %>% mutate(age = age-n_t+i-1) %>% 
    .[age %chin% Ages]
}

# age-wise cost for patients experiencing severe malaria

Trt_F %<d-% matrix(0, nrow = n_i, ncol = n_b)
for (i in 1:n_b)
{
  Trt_F[grep("SM", Trace[,i]),i] <- as.numeric(Cost[grep("SM", Trace[,i]),i]) 
}

age_specific_SM_cost <- list()
age_fil <- age_filter[ID %chin% 1:(n_0 + sum(n_c_babies[1:(length(n_c_babies)
                                                           %%(n_days/step_size))])),]
age_specific_SM_cost[[1]] <-  Trt_F %>% 
  cbind(age_filter[,c("ID", "age")], rowSums(.[,1:n_w])) %>%
  select((n_b+1):(n_b+3)) %>% .[!(State_Matrix[,n_w] %like% "D"),] %>%
  .[ID %chin% age_fil$ID] %>% mutate(age = age-n_t-1) %>% .[age %chin% Ages]
for (i in 1:n_t)
{
  age_fil <- age_filter[ID %chin% 1:(n_0 + 
                                       sum(n_c_babies[1:(length(n_c_babies)%%
                                                           (n_days/step_size) + (n_w*i))])),]
  age_specific_SM_cost [[i+1]] <- Trt_F %>%
    cbind(age_filter[,c("ID", "age")], rowSums(.[,((n_w*i+1):(n_w*(i+1)))])) %>%
    select((n_b+1):(n_b+3)) %>% .[!(State_Matrix[,(n_w*i)] %like% "D"),] %>% 
    .[ID %chin% age_fil$ID] %>% mutate(age = age-n_t+i-1) %>% 
    .[age %chin% Ages]
}

# age-wise cost for patients experiencing treatment failure & Severe malaria

Trt_F %<d-% matrix(0, nrow = n_i, ncol = n_b)
for (i in 1:n_b)
{
  Trt_F[grep("F.*SM", Trace[,i]),i] <- as.numeric(Cost[grep("F.*SM", Trace[,i]),i]) 
}

age_specific_F_SM_cost <- list()
age_fil <- age_filter[ID %chin% 1:(n_0 + sum(n_c_babies[1:(length(n_c_babies)
                                                           %%(n_days/step_size))])),]
age_specific_F_SM_cost[[1]] <-  Trt_F %>% 
  cbind(age_filter[,c("ID", "age")], rowSums(.[,1:n_w])) %>%
  select((n_b+1):(n_b+3)) %>% .[!(State_Matrix[,n_w] %like% "D"),] %>%
  .[ID %chin% age_fil$ID] %>% mutate(age = age-n_t-1) %>% .[age %chin% Ages]
for (i in 1:n_t)
{
  age_fil <- age_filter[ID %chin% 1:(n_0 + 
                                       sum(n_c_babies[1:(length(n_c_babies)%%
                                                           (n_days/step_size) + (n_w*i))])),]
  age_specific_F_SM_cost [[i+1]] <- Trt_F %>%
    cbind(age_filter[,c("ID", "age")], rowSums(.[,((n_w*i+1):(n_w*(i+1)))])) %>%
    select((n_b+1):(n_b+3)) %>% .[!(State_Matrix[,(n_w*i)] %like% "D"),] %>% 
    .[ID %chin% age_fil$ID] %>% mutate(age = age-n_t+i-1) %>% 
    .[age %chin% Ages]
}


#################### Visualizations ############################################

###################### Yearly Costs ############################################

# Overall Mean Cost
yearly_cost %>% gather (.[,as.character((y_0+1):(start_year-1))], key = "year", value = "cost") %>%
  group_by(age, year)  %>% summarise_at(vars(cost), mean) %>%
  ggplot(.) +
  geom_line(aes(age, cost, colour = year))+
  theme_classic() + guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous("Average Cost (USD)")+
  labs (colour = "") +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        legend.text = element_text( face = "bold")) 



# Mean Cost for those seeking formal care
yearly_cost %>% gather (.[,as.character((y_0+1):(start_year-1))], key = "year", value = "cost") %>% 
  filter(cost>0) %>%
  group_by(age, year)  %>% summarise_at(vars(cost), mean) %>%
  ggplot(.) +
  geom_line(aes(age, cost, colour = year))+
  theme_classic() + guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous("Average Cost (USD)")+
  labs (colour = "") +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        legend.text = element_text( face = "bold")) 

# Total Cost for those seeking formal care
yearly_cost %>% gather (.[,as.character((y_0+1):(start_year-1))], key = "year", value = "cost") %>% 
  filter(cost>0) %>%
  group_by(age, year)  %>% summarise_at(vars(cost), sum) %>%
  ggplot(.) +
  geom_line(aes(age, cost, colour = year))+
  theme_classic() + guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous("Total Cost (USD)")+
  labs (colour = "") +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        legend.text = element_text( face = "bold")) 

###################### Age-Specific Costs ######################################

# Overal Mean Cost

avg <- function(x){group_by(x, age) %>% 
    summarise_at(vars(V3), mean) %>% select(-age)}

dplyr::bind_cols(lapply(age_specific_cost, avg)) %>% `colnames<-` (y_0: (y_0+n_t)) %>%
  bind_cols(age = Ages,.) %>% 
  as.data.table() %>% select(-as.character(y_0)) %>%
  pivot_longer (.[,as.character((y_0+1):(start_year-1))], 
                names_to = "year", values_to = "cost") %>%
  ggplot(.) +
  geom_line(aes(age, cost, colour = year))+
  theme_classic() + guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous("Average Cost (USD)")+
  labs (colour = "") +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        legend.text = element_text( face = "bold")) 

# Mean Cost for those seeking formal care

avg1 <- function(x){group_by(x, age) %>% filter(V3>0) %>%
    summarise_at(vars(V3), mean) %>% select(-age)}

dplyr::bind_cols(lapply(age_specific_cost, avg1)) %>% `colnames<-` (y_0: (y_0+n_t)) %>%
  bind_cols(age = Ages,.) %>% 
  as.data.table() %>% select(-as.character(y_0)) %>%
  pivot_longer (.[,as.character((y_0+1):(start_year-1))], 
                names_to = "year", values_to = "cost") %>% filter(cost>0) %>%
  ggplot(.) +
  geom_line(aes(age, cost, colour = year))+
  theme_classic() + guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous("Average Cost (USD)")+
  labs (colour = "") +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        legend.text = element_text( face = "bold")) 

# Total Cost for those seeking formal care

tot <- function(x){group_by(x, age) %>% filter(V3>0) %>%
    summarise_at(vars(V3), sum) %>% select(-age)}

dplyr::bind_cols(lapply(age_specific_cost, tot)) %>% `colnames<-` (y_0: (y_0+n_t)) %>%
  bind_cols(age = Ages,.) %>% 
  as.data.table() %>% select(-as.character(y_0)) %>%
  pivot_longer (.[,as.character((y_0+1):(start_year-1))], 
                names_to = "year", values_to = "cost") %>% filter(cost>0) %>%
  ggplot(.) +
  geom_line(aes(age, cost, colour = year))+
  theme_classic() + guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous("Total Cost (USD)")+
  labs (colour = "") +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        legend.text = element_text( face = "bold")) 

###################### Treatment Failure Age-Specific Costs ####################

# Mean Cost for those seeking formal care
mat %<d-% matrix(0, nrow=length(Ages), ncol = (n_t+1))
for (i in 1:(n_t+1)){
  cc <- group_by(age_specific_Trt_F_cost[[i]],age) %>% filter(V3>0) %>%
    summarise_at(vars(V3), mean) 
  mat[cc$age+1,i] <- cc$V3 
}

mat %>% `colnames<-` (y_0: (y_0+n_t)) %>% cbind(age = Ages, .) %>%
  as.data.table() %>% select(-as.character(y_0)) %>%
  pivot_longer (.[,as.character((y_0+1):(start_year-1))], 
                names_to = "year", values_to = "cost") %>% filter(cost>0) %>%
  ggplot(.) +
  geom_line(aes(age, cost, colour = year))+
  theme_classic() + guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous("Average Cost (USD) - Treatment Failure")+
  labs (colour = "") +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        legend.text = element_text( face = "bold")) 

# Total Cost for those seeking formal care
mat %<d-% matrix(0, nrow=length(Ages), ncol = (n_t+1))
for (i in 1:(n_t+1)){
  cc <- group_by(age_specific_Trt_F_cost[[i]],age) %>% filter(V3>0) %>%
    summarise_at(vars(V3), sum) 
  mat[cc$age+1,i] <- cc$V3 
}

mat %>% `colnames<-` (y_0: (y_0+n_t)) %>% cbind(age = Ages, .) %>%
  as.data.table() %>% select(-as.character(y_0)) %>%
  pivot_longer (.[,as.character((y_0+1):(start_year-1))], 
                names_to = "year", values_to = "cost") %>% filter(cost>0) %>%
  ggplot(.) +
  geom_line(aes(age, cost, colour = year))+
  theme_classic() + guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous("Total Cost (USD) - Treatment Failure")+
  labs (colour = "") +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        legend.text = element_text( face = "bold")) 

###################### Severe Malaria Age-Specific Costs ####################

# Mean Cost for those seeking formal care
mat %<d-% matrix(0, nrow=length(Ages), ncol = (n_t+1))
for (i in 1:(n_t+1)){
  cc <- group_by(age_specific_SM_cost[[i]],age) %>% filter(V3>0) %>%
    summarise_at(vars(V3), mean) 
  mat[cc$age+1,i] <- cc$V3 
}

mat %>% `colnames<-` (y_0: (y_0+n_t)) %>% cbind(age = Ages, .) %>%
  as.data.table() %>% select(-as.character(y_0)) %>%
  pivot_longer (.[,as.character((y_0+1):(start_year-1))], 
                names_to = "year", values_to = "cost") %>% filter(cost>0) %>%
  ggplot(.) +
  geom_line(aes(age, cost, colour = year))+
  theme_classic() + guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous("Average Cost (USD) - Severe Malaria")+
  labs (colour = "") +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        legend.text = element_text( face = "bold")) 

# Total Cost for those seeking formal care
mat %<d-% matrix(0, nrow=length(Ages), ncol = (n_t+1))
for (i in 1:(n_t+1)){
  cc <- group_by(age_specific_SM_cost[[i]],age) %>% filter(V3>0) %>%
    summarise_at(vars(V3), sum) 
  mat[cc$age+1,i] <- cc$V3 
}

mat %>% `colnames<-` (y_0: (y_0+n_t)) %>% cbind(age = Ages, .) %>%
  as.data.table() %>% select(-as.character(y_0)) %>%
  pivot_longer (.[,as.character((y_0+1):(start_year-1))], 
                names_to = "year", values_to = "cost") %>% filter(cost>0) %>%
  ggplot(.) +
  geom_line(aes(age, cost, colour = year))+
  theme_classic() + guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous("Total Cost (USD) - Severe Malaria")+
  labs (colour = "") +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        legend.text = element_text( face = "bold")) 

############## Treatment Failure & Severe Malaria Age-Specific Costs ###########

# Mean Cost for those seeking formal care
mat %<d-% matrix(0, nrow=length(Ages), ncol = (n_t+1))
for (i in 1:(n_t+1)){
  cc <- group_by(age_specific_F_SM_cost[[i]],age) %>% filter(V3>0) %>%
    summarise_at(vars(V3), mean) 
  mat[cc$age+1,i] <- cc$V3 
}

mat %>% `colnames<-` (y_0: (y_0+n_t)) %>% cbind(age = Ages, .) %>%
  as.data.table() %>% select(-as.character(y_0)) %>%
  pivot_longer (.[,as.character((y_0+1):(start_year-1))], 
                names_to = "year", values_to = "cost") %>% filter(cost>0) %>%
  ggplot(.) +
  geom_line(aes(age, cost, colour = year))+
  theme_classic() + guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous("Average Cost (USD) - Treatment Failure & Severe")+
  labs (colour = "") +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        legend.text = element_text( face = "bold")) 

# Total Cost for those seeking formal care
mat %<d-% matrix(0, nrow=length(Ages), ncol = (n_t+1))
for (i in 1:(n_t+1)){
  cc <- group_by(age_specific_F_SM_cost[[i]],age) %>% filter(V3>0) %>%
    summarise_at(vars(V3), sum) 
  mat[cc$age+1,i] <- cc$V3 
}

mat %>% `colnames<-` (y_0: (y_0+n_t)) %>% cbind(age = Ages, .) %>%
  as.data.table() %>% select(-as.character(y_0)) %>%
  pivot_longer (.[,as.character((y_0+1):(start_year-1))], 
                names_to = "year", values_to = "cost") %>% filter(cost>0) %>%
  ggplot(.) +
  geom_line(aes(age, cost, colour = year))+
  theme_classic() + guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous("Total Cost (USD) - Treatment Failure & Severe Malaria")+
  labs (colour = "") +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        legend.text = element_text( face = "bold")) 

######################## Care Seeking Proportions ##############################

# Plot health state trace
m_TR_count <- t(apply(State_Matrix[,-c(1:n_w)], 2, function(x) table(factor(x, levels = 
                                c("H", "D", "UM","SM"), ordered = TRUE))))

ggplot(as.data.frame(m_TR_count)) +
  geom_line(aes(1:(n_b-n_w), H, colour = "Malaria Negative")) +
  geom_line(aes(1:(n_b-n_w), D, colour = "Dead")) +
  geom_line(aes(1:(n_b-n_w), UM, colour = "Uncomplicated Malaria")) +
  geom_line(aes(1:(n_b-n_w), SM, colour = "Severe Malaria"))+
  theme_classic() +
  scale_x_continuous("Cycle")+
  scale_y_continuous("Number of individuals")+
  labs (colour = "") +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        legend.text = element_text(face = "bold")) 

# Plot health state trace
m_TR_count <- t(apply(Trace[,-c(1:n_w)], 2, function(x) table(factor(x, 
                            levels = c(unique(Trace[,1])[-1]), ordered = TRUE))))

# Proportion of Malaria Negative Fever

data.table(NM = rowSums(m_TR_count[,c("NNC", "NHD", "NHND", "NHCD", "NHCND")]),
      M = rowSums(m_TR_count) - 
        rowSums(m_TR_count[,c("NNC", "NHD", 
                              "NHND", "NHCD", "NHCND")]))[,Tot := NM + M][,`:=`(NM = NM/Tot,
                                                                          M = M/Tot)] %>%
  ggplot() + geom_line(aes(1:(n_b-n_w), NM, colour = "NM")) +
  theme_classic() + scale_color_manual(values=c("#00CC00"))  + 
  scale_x_continuous("Cycle")+
  scale_y_continuous("Proportion with Malaria Negative Fever", labels = scales::percent)+
  labs (colour = "") +
  theme(legend.position = "",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        legend.text = element_text(face = "bold")) 

# Proportion of CareSeeking by Malaria Negative Fever Patients

data.table(C = m_TR_count[,"NNC"],
             NC = rowSums(m_TR_count[,c("NHD", "NHND", "NHCD", 
               "NHCND")]))[,Tot := C + NC][,`:=`(C = C/Tot, NC = NC/Tot)] %>%
  ggplot() + geom_line(aes(1:(n_b-n_w), C, colour = "C")) +
  theme_classic() + scale_color_manual(values=c("#00CC00"))  + 
  scale_x_continuous("Cycle")+
  scale_y_continuous("Proportion with Malaria Negative Fever Seeking Formal Care", labels = scales::percent)+
  labs (colour = "") +
  theme(legend.position = "",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        legend.text = element_text(face = "bold")) 

# Treatment Failures leading to severe malaria

data.table(F = m_TR_count[, "C_F_UM_SM_HC"]) %>%
  ggplot() + geom_line(aes(1:(n_b-n_w), F, colour = "F")) + theme_classic() +
  scale_x_continuous("Cycle") + scale_color_manual(values=c("#00CC00"))  +
  scale_y_continuous("Number with Severe Malaria after Treatment Failure")+
  labs (colour = "") +
  theme(legend.position = "",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        legend.text = element_text(face = "bold")) 

# Copy # Total Cost for those seeking formal care
mat %<d-% matrix(0, nrow=length(Ages), ncol = (n_t+1))
for (i in 1:(n_t+1)){
  cc <- group_by(age_specific_F_SM_cost[[i]],age) %>% filter(V3>0) %>%
    summarise_at(vars(V3), sum) 
  mat[cc$age+1,i] <- cc$V3 
}

mat %>% `colnames<-` (y_0: (y_0+n_t)) %>% cbind(age = Ages, .) %>%
  as.data.table() %>% select(-as.character(y_0)) %>%
  pivot_longer (.[,as.character((y_0+1):(start_year-1))], 
                names_to = "year", values_to = "cost") %>% filter(cost>0) %>%
  ggplot(.) +
  geom_line(aes(age, cost, colour = year))+
  theme_classic() + guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous("Total Cost (USD) - Treatment Failure & Severe Malaria")+
  labs (colour = "") +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        legend.text = element_text( face = "bold")) 

# Summary of Costs
# Get unique entries of Trace

my.names <- Trace %>% as.matrix() %>% as.vector() %>% unique()
sum.cost <- c(mean(unlist(sapply((n_days/step_size+1):n_b, function(x) 
  as.numeric(Cost[grep(my.names[-1][1], Trace[,x]),x])))),
  sum(unlist(sapply((n_w+1):length(Trace[1,]), function(x) 
    as.numeric(Cost[grep(my.names[-1][1], Trace[,x]),x]))))/n_t)
  for (i in 2: length(my.names[-1]))
  {
    sum.cost <- rbind(sum.cost, 
                      c(mean(unlist(sapply((n_days/step_size+1):n_b, function(x) 
                        as.numeric(Cost[grep(my.names[-1][i], Trace[,x]),x])))),
                        sum(unlist(sapply((n_days/step_size+1):n_b, function(x) 
                          as.numeric(Cost[grep(my.names[-1][i], Trace[,x]),x]))))/n_t))
  }

colnames(sum.cost) <- c("Mean Cost", "Total Cost")
rownames(sum.cost) <- my.names[-1]

sum.cost

#END
