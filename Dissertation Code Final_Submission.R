
# Clear the current R environment of any existing variables to avoid naming conflicts
rm(list = ls())

# Install and load tidyverse if not already installed
if (!require(tidyverse)) install.packages("tidyverse")
library(tidyverse)

# Install and load deSolve if not already installed
if (!require(deSolve)) install.packages("deSolve")
library(deSolve)

# Set the working directory to the specified path where the data and outputs will be managed
setwd("Insert Working Directory Here")

######################################################################
#LOADING DATA 
# Load Lao PDR Age-Stratified Census Data and Projections for (1995, 2005, 2015, 2020, 2025, 2030, 2035, 2040, 2045)
#https://laosis.lsb.gov.la/board/BoardList.do?bbs_bbsid=B404
popstruc <- as.data.frame(read.csv("Laos Population Structure Age Breakdown.csv", header = TRUE)) #creating 6 month age group to more actually capture vertical transmission
A <- length(popstruc[, 1]) #number of age groups within the model 

# Load Lao PDR Age-Stratified Census Data and Projections for (1995, 2005, 2015, 2020, 2025, 2030, 2035, 2040, 2045)
popstruc_1985 <- as.data.frame(read.csv("Laos Population Structure Age Breakdown - 1985.csv", header = TRUE)) #for initalizing the model 

#Load Lao PDR HBsAG Prevelance Data (2004, 2006, 2014, 2019, 2021)
HBsAG_prev_data <- read.csv("HBsAG Prevalence.csv", header = TRUE)

######################################################################
#ASSEMBLING AGEING MATRIX 
######################################################################

# Assemble Ageing Matrix
age_groups <- as.numeric(popstruc[,1]) #lower-limit of age group boundaries
age_group_durations <- diff(age_groups) * 365.25 #convert age group duration from years to days
ageing_matrix <- matrix(0, nrow = A, ncol = A) #create matrix of zeroes with number of rows and columns equal to number of age groups

# Populate Ageing Matrix with rates of ageing (first 5 columns of age matrix are for males ageging)
for (i in 1:A-1) {
  ageing_matrix[i,i] <- -1/age_group_durations[i] #rate of leaving an age group
  ageing_matrix[i + 1,i] <- 1/age_group_durations[i] #rate of entering the subsequent age group
}

# Assign names to rows and columns
age_group_names <- c( "0to6mnths", "7mnthsto4yrs", #creating 6 month age group to more actually capture vertical transmission
                      "5to9yrs", "10to14yrs", 
                      "15to19yrs", "20to24yrs", 
                      "25to29yrs", "30to34yrs", 
                      "35to39yrs", "40to44yrs", 
                      "45to49yrs", "50to54yrs", 
                      "55to59yrs", "60to64yrs", 
                      "65to69yrs", "70plusyrs")

rownames(ageing_matrix) <- age_group_names
colnames(ageing_matrix) <- age_group_names

######################################################################
#ASSEMBLING CONTACTS MATRICES
######################################################################

#Home Contact Matrix 
contacts_home <- matrix(c( # Prem et. al 2021, projected home contact matrix for Lao PDR adjusted for age-groups susceptible to infection at home (0to6mnths, 7mnthstoyrs )
  0.58759, 1.07285, 0.80480, 0.52906, 0.68695, 0.82761, 0.74483, 0.59712, 0.42232, 0.29616, 0.34158, 0.32728, 0.25837, 0.19890, 0.15114, 0.13466, #focus on contacts between children and older adults within the home
  0.58759, 1.07285, 0.80480, 0.52906, 0.68695, 0.82761, 0.74483, 0.59712, 0.42232, 0.29616, 0.34158, 0.32728, 0.25837, 0.19890, 0.15114, 0.13466,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000
), nrow = 16, byrow = TRUE)

rownames(contacts_home) <- age_group_names
colnames(contacts_home) <- age_group_names

#School Contact Matrix 
contacts_school <- matrix(c(# Prem et. al 2021, projected home contact matrix for Lao PDR adjusted for age-groups susceptible to infection at school (5to9yrs, 10to14yrs)
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 3.48906, 0.19761, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.75217, 4.43910, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000
), nrow = 16, byrow = TRUE)

rownames(contacts_school) <- age_group_names
colnames(contacts_school) <- age_group_names

#Sexual Contact Matrix 
contacts_sex <- matrix(c( # Calvert et. al 2023, average PAR of men in Cambodia (8.13) and women in Cambodia (3.76) to get 5.945 rounded to 6 adjusted for age-groups susceptible to infection from sexual contact (15to49yrs)
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
), nrow = 16, byrow = TRUE)
contacts_sex <- (0.6*contacts_sex) + ((0.4*contacts_sex)/(A - 1)) # Chantavilay et. al 2016, formula for sexual contact matrix for HPV in Lao PDR
rownames(contacts_sex) <- age_group_names
colnames(contacts_sex) <- age_group_names

###########################################################################
# DEFINE THE INDICES FOR EACH VARIABLE (each index represents age group within a health state)
###########################################################################

# Health State Compartment Indicies within the model 

#Indices for Male Population 
S_index_m <- 1 : A #susceptible male population
Acu_index_m <- (A + 1) : (2 * A) #acute male infection
Chron1_index_m <- (2 * A + 1) : (3 * A) #immune tolerant chronic male infection
Chron2_index_m <- (3 * A + 1) : (4 * A) #immune clearance chronic male infection
Chron3_index_m <- (4 * A + 1) : (5 * A) #inactive carrier chronic male infection 
Chron4_index_m <- (5 * A + 1) : (6 * A) #immune reactivation chronic male infection 
CompC_index_m <- (6 * A + 1) : (7 * A) #compensated cirrhosis in males
DeCompC_index_m <- (7 * A + 1) : (8 * A) #decompensated cirrhosis in males
HCC_index_m <- (8 * A + 1) : (9 * A) #hepatocellular carcinoma in males 
R_index_m <- (9 * A + 1) : (10 * A) #recovered males
Vac_index_m <- (10 * A + 1) : (11 * A) #vaccination in males 

#Indices for Female Population 
S_index_f <- (11 * A + 1) : (12 * A) #susceptible female population
Acu_index_f <- (12 * A + 1) : (13 * A)   #acute female infection 
Chron1_index_f <- (13 * A + 1) : (14 * A) #immune tolerant chronic female infection 
Chron2_index_f <- (14 * A + 1) : (15 * A)  #immune clearance chronic female infection
Chron3_index_f <- (15 * A + 1) : (16 * A) #inactive carrier chronic female infection 
Chron4_index_f <- (16 * A + 1) : (17 * A) #immune reactivation chronic female infection 
CompC_index_f <- (17 * A + 1) : (18 * A) #compensated cirrhosis in females
DeCompC_index_f <- (18 * A + 1) : (19 * A) #decompensated cirrhosis in females 
HCC_index_f <- (19 * A + 1) : (20 * A) #hepatocellular carcinoma in females
R_index_f <- (20 * A + 1) : (21 * A) #recovered females
Vac_index_f <- (21 * A + 1) : (22 * A) #vaccination in females 

#Counters for tracking various dynamics within the model
CMort_index <- (22 * A + 1) : (23 * A) #cumulative disease-induced mortality 
CInc_H_index <-  (23 * A + 1) : (24 * A) #cumulative incidence of all new cases (all new cases start as acute)
CInc_V_index <-  (24 * A + 1) : (25 * A) #cumulative incidence of all new cases (all new cases start as acute)
CChron_index <- (25 * A + 1) : (26 * A) #cumulative incidence of chronic cases (all new case chronic cases start as immune tolerant)
CCompC_index <- (26 * A + 1) : (27 * A) # cumulative incidence of compensated cirrhosis 
CDeCompC_index <- (27 * A + 1) : (28 * A) #cumulative incidence of decompensated cirrhosis 
CHCC_index <- (28 * A + 1) : (29 * A) # cumulative incidence of hepatocellular carcinoma
CVac_index <- (29 * A + 1) : (30 * A) #cumulative vaccinated population

###########################################################################
# MODEL INITIAL CONDITIONS 
##########################################################################

# Initial conditions for males
initS_m <- popstruc_1985[,3] - c(1000, # "0to6mnths"
                                 1000, # "7mnthsto4yrs"
                                 0.03 * popstruc_1985[3,3], # "5to9yrs" 
                                 0.08 * popstruc_1985[4,3], # "10to14yrs" 
                                 0.11 * popstruc_1985[5,3], # "15to19yrs"
                                 0.07 * popstruc_1985[6,3] + 300, # "20to24yrs"
                                 0.05 * popstruc_1985[7,3] + 300, # "25to29yrs"
                                 0.05 * popstruc_1985[8,3] + 900, # "30to34yrs"
                                 0.05 * popstruc_1985[9,3] + 900, #  "35to39yrs"
                                 0.05 * popstruc_1985[10,3] + 1800, # "40to44yrs" 
                                 0.05 * popstruc_1985[11,3] + 1800, # "45to49yrs" 
                                 0.05 * popstruc_1985[12,3] + 1800, # "50to54yrs" 
                                 0.05 * popstruc_1985[13,3] + 900, #  "55to59yrs"  
                                 0.05 * popstruc_1985[14,3] + 600, # "60to64yrs"
                                 0.05 * popstruc_1985[15,3] + 300, # "65to69yrs" 
                                 0.05 * popstruc_1985[16,3] + 300) # 70plusyrs"
initAcu_m <- rep(0, A)
initAcu_m[1:2] <- 1000
initChron1_m <- c(0, 
                  0, 
                  0.00 * popstruc_1985[3,3], 
                  0.00 * popstruc_1985[4,3], 
                  0.02 * popstruc_1985[5,3], 
                  0.02 * popstruc_1985[6,3], 
                  0.00 * popstruc_1985[7,3], 
                  0.00 * popstruc_1985[8,3], 
                  0.00 * popstruc_1985[9,3], 
                  0.00 * popstruc_1985[10,3],  
                  0.00 * popstruc_1985[11,3],  
                  0.00 * popstruc_1985[12,3], 
                  0.00 * popstruc_1985[13,3],  
                  0.00 * popstruc_1985[14,3],  
                  0.00 * popstruc_1985[15,3],  
                  0.00 * popstruc_1985[16,3])
initChron2_m <- c(0, 0, 
                  0.03 * popstruc_1985[3,3], 
                  0.03 * popstruc_1985[4,3], 
                  0.04 * popstruc_1985[5,3], 
                  0.04 * popstruc_1985[6,3], 
                  0.02 * popstruc_1985[7,3], 
                  0.02 * popstruc_1985[8,3], 
                  0.02 * popstruc_1985[9,3], 
                  0.02 * popstruc_1985[10,3],  
                  0.02 * popstruc_1985[11,3],  
                  0.02 * popstruc_1985[12,3], 
                  0.02 * popstruc_1985[13,3],  
                  0.02 * popstruc_1985[14,3],  
                  0.02 * popstruc_1985[15,3],  
                  0.02 * popstruc_1985[16,3])
initChron3_m <- c(0, 0, 
                  0.00 * popstruc_1985[3,3], 
                  0.00 * popstruc_1985[4,3], 
                  0.02 * popstruc_1985[5,3], 
                  0.02 * popstruc_1985[6,3], 
                  0.02 * popstruc_1985[7,3], 
                  0.02 * popstruc_1985[8,3], 
                  0.02 * popstruc_1985[9,3], 
                  0.02 * popstruc_1985[10,3],  
                  0.02 * popstruc_1985[11,3],  
                  0.02 * popstruc_1985[12,3], 
                  0.02 * popstruc_1985[13,3],  
                  0.02 * popstruc_1985[14,3],  
                  0.02 * popstruc_1985[15,3],  
                  0.02 * popstruc_1985[16,3])
initChron4_m <- c(0, 0, 
                  0.00 * popstruc_1985[3,3], 
                  0.00 * popstruc_1985[4,3], 
                  0.00 * popstruc_1985[5,3], 
                  0.07 * popstruc_1985[6,3], 
                  0.01 * popstruc_1985[7,3], 
                  0.01 * popstruc_1985[8,3], 
                  0.01 * popstruc_1985[9,3], 
                  0.01 * popstruc_1985[10,3],  
                  0.01 * popstruc_1985[11,3],  
                  0.01 * popstruc_1985[12,3], 
                  0.01 * popstruc_1985[13,3],  
                  0.01 * popstruc_1985[14,3],  
                  0.01 * popstruc_1985[15,3],  
                  0.01 * popstruc_1985[16,3])
initCompC_m <- c(0, 0, 0, 0, 0, 100, 100, 300, 300, 600, 600, 600, 300, 200, 100, 100)
initDeCompC_m <- c(0, 0, 0, 0, 0, 100, 100, 300, 300, 600, 600, 600, 300, 200, 100, 100)
initHCC_m <- c(0, 0, 0, 0, 0, 100, 100, 300, 300, 600, 600, 600, 300, 200, 100, 100)
initR_m <- rep(0, A)
initVac_m <- rep(0, A)

# Initial conditions for females
initS_f <- popstruc_1985[,4] - c(1000, 1000, 
                                 0.03 * popstruc_1985[3,4], 
                                 0.075 * popstruc_1985[4,4], 
                                 0.11 * popstruc_1985[5,4], 
                                 0.07 * popstruc_1985[6,4] + 100, 
                                 0.05 * popstruc_1985[7,4] + 100, 
                                 0.05 * popstruc_1985[8,4] + 300, 
                                 0.05 * popstruc_1985[9,4] + 300, 
                                 0.05 * popstruc_1985[10,4] + 600,  
                                 0.05 * popstruc_1985[11,4] + 600,  
                                 0.05 * popstruc_1985[12,4] + 600, 
                                 0.05 * popstruc_1985[13,4] + 300,  
                                 0.05 * popstruc_1985[14,4] + 200,  
                                 0.05 * popstruc_1985[15,4] + 100,  
                                 0.05 * popstruc_1985[16,4] + 100)
initAcu_f <- rep(0, A)
initAcu_f[1:2] <- 1000
initChron1_f <- c(0, 0, 
                  0.00 * popstruc_1985[3,4], 
                  0.00 * popstruc_1985[4,4], 
                  0.02 * popstruc_1985[5,4], 
                  0.02 * popstruc_1985[6,4], 
                  0.00 * popstruc_1985[7,4], 
                  0.00 * popstruc_1985[8,4], 
                  0.00 * popstruc_1985[9,4], 
                  0.00 * popstruc_1985[10,4],  
                  0.00 * popstruc_1985[11,4],  
                  0.00 * popstruc_1985[12,4], 
                  0.00 * popstruc_1985[13,4],  
                  0.00 * popstruc_1985[14,4],  
                  0.00 * popstruc_1985[15,4],  
                  0.00 * popstruc_1985[16,4])
initChron2_f <- c(0, 0, 
                  0.03 * popstruc_1985[3,4], 
                  0.03 * popstruc_1985[4,4], 
                  0.04 * popstruc_1985[5,4], 
                  0.04 * popstruc_1985[6,4], 
                  0.02 * popstruc_1985[7,4], 
                  0.02 * popstruc_1985[8,4], 
                  0.02 * popstruc_1985[9,4], 
                  0.02 * popstruc_1985[10,4],  
                  0.02 * popstruc_1985[11,4],  
                  0.02 * popstruc_1985[12,4], 
                  0.02 * popstruc_1985[13,4],  
                  0.02 * popstruc_1985[14,4],  
                  0.02 * popstruc_1985[15,4],  
                  0.02 * popstruc_1985[16,4])
initChron3_f <- c(0, 0, 
                  0.00 * popstruc_1985[3,4], 
                  0.00 * popstruc_1985[4,4], 
                  0.02 * popstruc_1985[5,4], 
                  0.02 * popstruc_1985[6,4], 
                  0.02 * popstruc_1985[7,4], 
                  0.02 * popstruc_1985[8,4], 
                  0.02 * popstruc_1985[9,4], 
                  0.02 * popstruc_1985[10,4],  
                  0.02 * popstruc_1985[11,4],  
                  0.02 * popstruc_1985[12,4], 
                  0.02 * popstruc_1985[13,4],  
                  0.02 * popstruc_1985[14,4],  
                  0.02 * popstruc_1985[15,4],  
                  0.02 * popstruc_1985[16,4])
initChron4_f <- c(0, 0, 
                  0.00 * popstruc_1985[3,4], 
                  0.00 * popstruc_1985[4,4], 
                  0.00 * popstruc_1985[5,4], 
                  0.07 * popstruc_1985[6,4], 
                  0.01 * popstruc_1985[7,4], 
                  0.01 * popstruc_1985[8,4], 
                  0.01 * popstruc_1985[9,4], 
                  0.01 * popstruc_1985[10,4],  
                  0.01 * popstruc_1985[11,4],  
                  0.01 * popstruc_1985[12,4], 
                  0.01 * popstruc_1985[13,4],  
                  0.01 * popstruc_1985[14,4],  
                  0.01 * popstruc_1985[15,4],  
                  0.01 * popstruc_1985[16,4])
initCompC_f <- c(0, 0, 0, 0, 0, 33, 33, 100, 100, 200, 200, 200, 100, 67, 33, 33)
initDeCompC_f <- c(0, 0, 0, 0, 0, 33, 33, 100, 100, 200, 200, 200, 100, 67, 33, 33)
initHCC_f <- c(0, 0, 0, 0, 0, 33, 33, 100, 100, 200, 200, 200, 100, 67, 33, 33)
initR_f <- rep(0, A)
initVac_f <- rep(0, A)

#Counters
initCMort <- rep(0, A)
initCInc_H <- rep(0, A)
initCInc_V <- rep(0, A)
initCChron <- rep(0, A)
initCCompC <- rep(0, A)
initCDeCompC <- rep(0, A)
initCHCC <- rep(0, A)
initCVac <- rep(0, A)

init_laopdr_HBV_model_bc1 <- c(
  initS_m, initAcu_m, initChron1_m, initChron2_m, initChron3_m, 
  initChron4_m, initCompC_m, initDeCompC_m, initHCC_m, initR_m, initVac_m, 
  initS_f, initAcu_f, initChron1_f, initChron2_f, initChron3_f, 
  initChron4_f, initCompC_f, initDeCompC_f, initHCC_f, initR_f, initVac_f,
  initCMort, initCInc_H, initCInc_V, initCChron, initCCompC, 
  initCDeCompC, initCHCC, initCVac
)


###########################################################################
# MODEL PARAMETERS 
##########################################################################

# Model Parameters List (each parameter is a vector of length A (16) corresponding to the number of age groups in the model)
parms <- list(
  muDCC = rep(4.45E-04, A),                                                                     # HBV-induced mortality (decompensated cirrhosis)
  muHCC = c(rep(0, 2), rep(1.03E-03, (A-2))),                                                   # HBV-induced mortality (hepatocellular carcinoma)
  gamma1 = rep(1/180, A),                                                                       # duration of acute infection in an individual aged 0 - 4 years of age
  gamma2 = rep(1/180, A),                                                                       # duration of acute infection in an individual 5 years of age or older
  rho_c2 = c(rep((0.0019/365.25), 4), rep((0.004/365.25), 7), rep((0.006/365.25), 5)),          # rate of HBsAG seroclearance from Chron 2
  rho_c3 = c(rep((0.0056/365.25), 4), rep((0.0116/365.25), 7), rep((0.0183/365.25), 5)),        # rate of HBsAG seroclearance from Chron 3
  rho_c4 = c(rep((0.0056/365.25), 4), rep((0.0116/365.25), 7), rep((0.0183/365.25), 5)),                           # rate of HBsAG seroclearance from Chron 4
  r12 = c(rep((0.07/365.25), 2), rep((0.018/365.25), 2), rep((0.039/365.25), 7), rep((0.06/365.25), 5)),           # rate of progression from immune tolerant (Chron 1) to immune clearance (Chron 2)
  r23 = c(rep((0.05/365.25), 2), rep((0.09/365.25), 2), rep((0.15/365.35), 7), rep((0.11/365.25), 5)),             # rate of progression from immune immune clearance (Chron 2) to inactive carrier (Chron 3)
  r23_treated = c(rep((0.025/365.25), 2), rep((0.045/365.25), 2), rep((0.075/365.35), 7), rep((0.055/365.25), 5)), # rate of progression from immune immune clearance (Chron 2) to inactive carrier (Chron 3) on treatment
  r32 = c(rep((0.005/365.25), 2), rep((0.009/365.25), 2), rep((0.015/365.25), 7), rep((0.011/365.25), 5)),         # rate of de-progression from inactive carrier (Chron 3) back to immune clearance (Chron 2)
  r34 = c(rep((0.089/365.25), 2), rep((0.0089/365.25), 2), rep((0.0145/365.25), 7), rep((0.0151/365.25), 5)),      # rate of progression from inactive carrier (Chron 3) to immune reactivation (Chron 4)
  r34_treated = c(rep((0.045/365.25), 2), rep((0.0045/365.25), 2), rep((0.00725/365.25), 7), rep((0.00755/365.25), 5)), # rate of progression from inactive carrier (Chron 3) to immune reactivation (Chron 4) on treatment 
  epsilon1 = c(rep(0, 2), (0.0001/365.25), (0.0001/365.25), rep((0.0003/365.25), 7), rep((0.001/365.25), 5)),           # rate of progression from chronic 1 directly to hepatocellular carcinoma (HCC)
  epsilon1_treated = c(rep(0, 2), (0.00005/365.25), (0.00005/365.25), rep((0.00015/365.25), 7), rep((0.0005/365.25), 5)), # rate of progression from chronic 1 directly to hepatocellular carcinoma (HCC) on treatment  
  epsilon2_m = c(rep(0, 2), (0.0005/365.25), (0.0005/365.25), (0.0063/365.25), (0.0063/365.25), (0.0063/365.25), (0.0063/365.25), (0.0063/365.25), (0.0063/365.25), (0.0063/365.25), rep((0.0107/365.25), 4), (0.018/365.25)), # rate of progression from chronic 2 directly to hepatocellular carcinoma (HCC) in males
  epsilon2_m_treated = c(rep(0, 2), (0.00025/365.25), (0.00025/365.25), (0.00315/365.25), (0.00315/365.25), (0.00315/365.25), (0.00315/365.25), (0.00315/365.25), (0.00315/365.25), (0.00315/365.25), rep((0.00535/365.25), 4), (0.009/365.25)), # rate of progression from chronic 2 directly to hepatocellular carcinoma (HCC) in males on treatment 
  epsilon2_f = c(rep(0, 2), (0.0005/365.25), (0.0005/365.25), (0.0021/365.25), (0.0021/365.25), (0.0021/365.25), (0.0021/365.25), (0.0021/365.25), (0.0021/365.25), (0.0021/365.25), rep((0.0036/365.25), 4), (0.0061/365.25)),    # rate of progression from chronic 2 directly to hepatocellular carcinoma (HCC) in females
  epsilon2_f_treated = c(rep(0, 2), (0.00025/365.25), (0.00025/365.25), (0.00105/365.25), (0.00105/365.25), (0.00105/365.25), (0.00105/365.25), (0.00105/365.25), (0.00105/365.25), (0.00105/365.25), rep((0.0018/365.25), 4), (0.00305/365.25)),  # rate of progression from chronic 2 directly to hepatocellular carcinoma (HCC) in females on treatment
  epsilon3_m = c(rep(0, 2), (0.0001/365.25), (0.0001/365.25), (0.00075/365.25), (0.00075/365.25), (0.00075/365.25), (0.00075/365.25), (0.00075/365.25), (0.00075/365.25), (0.00075/365.25), rep((0.00128/365.25), 4), (0.0022/365.25)),    # rate of progression from chronic 3 directly to hepatocellular carcinoma (HCC) in males
  epsilon3_m_treated = c(rep(0, 2), (0.00005/365.25), (0.00005/365.25), (0.000375/365.25), (0.000375/365.25), (0.000375/365.25), (0.000375/365.25), (0.000375/365.25), (0.000375/365.25), (0.000375/365.25), rep((0.00064/365.25), 4), (0.0011/365.25)), # rate of progression from chronic 3 directly to hepatocellular carcinoma (HCC) in males on treatment
  epsilon3_f = c(rep(0, 2), (0.0005/365.25), (0.0005/365.25), (0.00025/365.25), (0.00025/365.25), (0.00025/365.25), (0.00025/365.25), (0.00025/365.25), (0.00025/365.25), (0.00025/365.25), rep((0.00043/365.25), 4), (0.00073/365.25)),    # rate of progression from chronic 3 directly to hepatocellular carcinoma (HCC) in females
  epsilon3_f_treated = c(rep(0, 2), (0.00025/365.25), (0.00025/365.25), (0.000125/365.25), (0.000125/365.25), (0.000125/365.25), (0.000125/365.25), (0.000125/365.25), (0.000125/365.25), (0.000125/365.25), rep((0.000215/365.25), 4), (0.000365/365.25)), # rate of progression from chronic 3 directly to hepatocellular carcinoma (HCC) in females on treatment 
  epsilon4_m = c(rep(0, 2), (0.0005/365.25), (0.0005/365.25), (0.0063/365.25), (0.0063/365.25), (0.0063/365.25), (0.0063/365.25), (0.0063/365.25), (0.0063/365.25), (0.0063/365.25), rep((0.0107/365.25), 4), (0.018/365.25)),    # rate of progression from chronic 4 directly to hepatocellular carcinoma (HCC) in males
  epsilon4_m_treated = c(rep(0, 2), (0.00025/365.25), (0.00025/365.25), (0.00315/365.25), (0.00315/365.25), (0.00315/365.25), (0.00315/365.25), (0.00315/365.25), (0.00315/365.25), (0.00315/365.25), rep((0.00535/365.25), 4), (0.009/365.25)),  # rate of progression from chronic 4 directly to hepatocellular carcinoma (HCC) in males on treatment 
  epsilon4_f = c(rep(0, 2), (0.0005/365.25), (0.0005/365.25), (0.0021/365.25), (0.0021/365.25), (0.0021/365.25), (0.0021/365.25), (0.0021/365.25), (0.0021/365.25), (0.0021/365.25), rep((0.0036/365.25), 4), (0.0061/365.21)), # rate of progression from chronic 4 directly to hepatocellular carcinoma (HCC) in females
  epsilon4_f_treated = c(rep(0, 2), (0.00025/365.25), (0.00025/365.25), (0.00105/365.25), (0.00105/365.25), (0.00105/365.25), (0.00105/365.25), (0.00105/365.25), (0.00105/365.25), (0.00105/365.25), rep((0.0018/365.25), 4), (0.00305/365.21)), # rate of progression from chronic 4 directly to hepatocellular carcinoma (HCC) in females on treatment
  delta2_m = c((0.0065/365.25), (0.0065/365.25), (0.0065/365.25), (0.0065/365.25), (0.0124/365.25), (0.0124/365.25), (0.0124/365.25), (0.0124/365.25), (0.0124/365.25), (0.0124/365.25), (0.0124/365.25), rep((0.0364/365.25), 5)), # rate of progression from immunce clearance chronic infection to Compensated Cirrhosis (CC) in males
  delta2_m_treated  = c((0.00325/365.25), (0.00325/365.25), (0.00325/365.25), (0.00325/365.25), (0.0062/365.25), (0.0062/365.25), (0.0062/365.25), (0.0062/365.25), (0.0062/365.25), (0.0062/365.25), (0.0062/365.25), rep((0.0182/365.25), 5)),  # rate of progression from immunce clearance chronic infection to Compensated Cirrhosis (CC) in males on treatment   
  delta2_f = c((0.0023/365.25), (0.0023/365.25), (0.0023/365.25), (0.0023/365.25), (0.0064/365.25), (0.0064/365.25), (0.0064/365.25), (0.0064/365.25), (0.0064/365.25), (0.0064/365.25), (0.0064/365.25), rep((0.0215/365.25), 5)), # rate of progression from immunce clearance chronic infection to Compensated Cirrhosis (CC) in females
  delta2_f_treated = c((0.00115/365.25), (0.00115/365.25), (0.00115/365.25), (0.0023/365.25), (0.0032/365.25), (0.0032/365.25), (0.0032/365.25), (0.0032/365.25), (0.0032/365.25), (0.0032/365.25), (0.0032/365.25), rep((0.01075/365.25), 5)), # rate of progression from immunce clearance chronic infection to Compensated Cirrhosis (CC) in females on treatment 
  delta4_m = c((0.0147/365.25), (0.0147/365.25), (0.0147/365.25), (0.0147/365.25), (0.0278/365.25), (0.0278/365.25), (0.0278/365.25), (0.0278/365.25), (0.0278/365.25), (0.0278/365.25), (0.0278/365.25), rep((0.082/365.25), 5)), # rate of progression from immune reactivation chronic infection to Compensated Cirrhosis (CC) in males
  delta4_m_treated = c((0.0147/365.25), (0.0147/365.25), (0.0147/365.25), (0.0147/365.25), (0.0278/365.25), (0.0278/365.25), (0.0278/365.25), (0.0278/365.25), (0.0278/365.25), (0.0278/365.25), (0.0278/365.25), rep((0.082/365.25), 5)), # rate of progression from immune reactivation chronic infection to Compensated Cirrhosis (CC) in males on treatment 
  delta4_f = c((0.053/365.25),(0.0530/365.25), (0.0530/365.25), (0.0530/365.25), (0.0143/365.25), (0.0143/365.25), (0.0143/365.25), (0.0143/365.25), (0.0143/365.25), (0.0143/365.25), (0.0143/365.25),  rep((0.0483/365.25), 5)),  # rate of progression from immune reactivation chronic infection to Compensated Cirrhosis (CC) in females
  delta4_f_treated = c((0.0265/365.25), (0.0265/365.25), (0.0265/365.25), (0.0265/365.25), (0.00715/365.25), (0.00715/365.25), (0.00715/365.25), (0.00715/365.25), (0.00715/365.25), (0.00715/365.25), (0.00715/365.25),  rep((0.02415/365.25), 5)),  # rate of progression from immune reactivation chronic infection to Compensated Cirrhosis (CC) in females on treatment
  alpha = rep((0.039/365.25), A),          # rate of progression from CC to Decompensated Cirrhosis (DCC)                                                          
  alpha_treated = rep((0.0195/365.25), A), # rate of progression from CC to Decompensated Cirrhosis (DCC) on treatment 
  theta_m = c(rep(0, 2), rep((0.0007/365.25), 2), rep((0.0271/365.25), 7), rep((0.0435/365.25), 4), (0.074/365.25)), # rate of progression from CC to HCC in males
  theta_f = c(rep(0, 2), rep((0.0007/365.25), 2), rep((0.0123/365.25), 7), rep((0.0356/365.25), 4), (0.0605/365.24)), # rate of progression from CC to HCC in females
  kappa = c(rep(0, 2), rep((0.04/365.25), (A-2))),         # rate of progression from DCC to HC
  kappa_treated = c(rep(0, 2), rep((0.02/365.25), (A-2))), # rate of progression from DCC to HC o treatmetn 
  p_chron1 = c(0.885, rep(0,(A-1))),                       # probability of moving from an acute infetion to Chronic Infection (Stage 1) (only possible if acquired during vertical transmission at birth)
  p_chron2 = c(0, 0.38, 0.17, rep(0.05, (A-3))),           # probability of moving from an acute infetion to Chronic Infection (Stage 2)
  p_MTCT = 0.75,                                          # probability of mother-to-child transmission
  prop_m = 0.5,                                           # proportion of births that are male
  prop_f = 0.5,                                           # proportion of births that are female
  p_infect = 0.0004,                                      # proportion of infections that result in successful transmission 
  zetaAcu = 0.1,                                          # reduced transmissibility of acute infections 
  vac_eff = 0.71,                                         # vaccination efficacy (no immunoglobin)
  vac_cov_2000 = 0.159, #vaccination coverage 2000 (assumed)
  vac_cov_2006 = 0.318, #vaccination coverage 2006 (Multiple Indicator Cluster Survey 2006)
  vac_cov_2011 = 0.515, #vaccination coverage 2011 (Lao Social Indicator Survey 2011)
  vac_cov_2017 = 0.608, #vaccination coverage 2017 (Lao Social Indicator Survey 2017)
  vac_cov_2023 = 0.613, #vaccination coverage 2011 (Lao Social Indicator Survey 2023)
  vac_cov_2029 = 0.613, # in base case we are assuming status quo but these values are altered in different scenarios 
  vac_cov_2035 = 0.613,
  vac_cov_2041 = 0.613,
  vac_cov_2047 = 0.613,
  vac_cov_2053 = 0.613,
  vac_cov_2059 = 0.613,
  vac_cov_2065 = 0.613,
  vac_cov_2071 = 0.613,
  vac_cov_2077 = 0.613,
  vac_cov_2083 = 0.613,
  vac_cov_2089 = 0.613,
  vac_cov_2095 = 0.613,
  prop_treat_2000 = 0.01, # in base case we are assuming status quo but these values are altered in different scenarios 
  prop_treat_2006 = 0.01, 
  prop_treat_2011 = 0.01,
  prop_treat_2017 = 0.01,
  prop_treat_2023 = 0.01,
  prop_treat_2029 = 0.01, 
  prop_treat_2035 = 0.01, 
  prop_treat_2041 = 0.01,
  prop_treat_2047 = 0.01,
  prop_treat_2053 = 0.01,
  prop_treat_2059 = 0.01,
  prop_treat_2065 = 0.01,
  prop_treat_2071 = 0.01,
  prop_treat_2077 = 0.01,
  prop_treat_2083 = 0.01,
  prop_treat_2089 = 0.01,
  prop_treat_2095 = 0.01
)

######################################################################
#BASE CASE CODE FOR DIFFERENTIAL EQUATIONS
######################################################################

# Differential Equations 
laopdr_HBV_model_bc1 <- function(t, Y, parms) {
  with(as.list(c(parms, Y)), {
    
    # Extract compartments from Y
    S_m <- Y[S_index_m]
    Acu_m  <- Y[Acu_index_m ]
    Chron1_m  <- Y[Chron1_index_m]
    Chron2_m  <- Y[Chron2_index_m]
    Chron3_m <- Y[Chron3_index_m]
    Chron4_m  <- Y[Chron4_index_m]
    CompC_m  <- Y[CompC_index_m]
    DeCompC_m <- Y[DeCompC_index_m]
    HCC_m <- Y[HCC_index_m]
    R_m <- Y[R_index_m]
    Vac_m <- Y[Vac_index_m]
    S_f <- Y[S_index_f]
    Acu_f <- Y[Acu_index_f]
    Chron1_f <- Y[Chron1_index_f]
    Chron2_f <- Y[Chron2_index_f]
    Chron3_f <- Y[Chron3_index_f]
    Chron4_f <- Y[Chron4_index_f]
    CompC_f <- Y[CompC_index_f]
    DeCompC_f <- Y[DeCompC_index_f]
    HCC_f <- Y[HCC_index_f]
    R_f <- Y[R_index_f]
    Vac_f <- Y[Vac_index_f]
    CMort <- Y[CMort_index]
    CInc_H <- Y[CInc_H_index]
    CInc_V <- Y[CInc_V_index]
    CChron <- Y[CChron_index]
    CCompC <- Y[CCompC_index]
    CDeCompC <- Y[CDeCompC_index]
    CHCC <- Y[CHCC_index]
    CVac <- Y[CVac_index]
    
    # Total population is the sum of all individuals in all compartments   
    P_m <- S_m + Acu_m + Chron1_m + Chron2_m + Chron3_m + Chron4_m + CompC_m + DeCompC_m + HCC_m + R_m + Vac_m #total ppoulation of males across all age groups and health states 
    P_f <- S_f + Acu_f + Chron1_f + Chron2_f + Chron3_f + Chron4_f + CompC_f + DeCompC_f + HCC_f + R_f + Vac_f #total population of females across all age groups and health states 
    P_total <- P_m + P_f #total population 
    
    vaccination <- 0 #initializing vaccines 
    treatment <- 0 #initializing treatment 
    
    #Time-Varying Birth and Age-Specific Death Rates 
    if (t < (365.25*15)) { #birth rate, age-specific death rates, vaccination rates and proportion treated between 1985 and 2000
      mu_b <- c(((1148/1000)/365.25), rep(0,A-1)) #people are only born into the first age group
      mu_d <- c(((75/1000)/365.25), ((97/1000)/365.25), ((21/1000)/365.25), ((21/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((30/1000)/365.25), ((30/1000)/365.25), ((30/1000)/365.25), ((30/1000)/365.25), ((150/1000)/365.25)) #Each index of the vector represents a different age group
    } else if ((365.25*15) && t < (365.25*21)) {  #birth rate, age-specific death rates, vaccination rates and proportion treated between 2000 and 2006
      mu_b <- c(((1090/1000)/365.25), rep(0,A-1))
      mu_d <- c(((75/1000)/365.25), ((97/1000)/365.25), ((21/1000)/365.25), ((21/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((30/1000)/365.25), ((30/1000)/365.25), ((30/1000)/365.25), ((30/1000)/365.25), ((150/1000)/365.25))
      vaccination <- (vac_cov_2000/365.25) * vac_eff
      treatment <- prop_treat_2000
    } else if ((365.25*21) && t < (365.25*26)) {  #birth rate, age-specific death rates, vaccination rates and proportion treated between 2006 and 2011
      mu_b <- c(((1065/1000)/365.25), rep(0,A-1))
      mu_d <- c(((75/1000)/365.25), ((79/1000)/365.25), ((16/1000)/365.25), ((16/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((30/1000)/365.25), ((30/1000)/365.25), ((30/1000)/365.25), ((30/1000)/365.25), ((150/1000)/365.25))
      vaccination <- (vac_cov_2006/365.25) * vac_eff
      treatment <- prop_treat_2006
    } else if ((365.25*26) && t < (365.25*32)) {  #birth rate, age-specific death rates, vaccination rates and proportion treated between 2011 and 2017
      mu_b <- c(((1050/1000)/365.25), rep(0,A-1))
      mu_d <- c(((68/1000)/365.25), ((79/1000)/365.25), ((11/1000)/365.25), ((11/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((150/1000)/365.25))
      vaccination <- (vac_cov_2011/365.35) * vac_eff
      treatment <- prop_treat_2011
    } else if ((365.25*32) && t < (365.25*38)) {  #birth rate, age-specific death rates, vaccination rates and proportion treated between 2017 and 2023
      mu_b <- c(((1045/1000)/365.25), rep(0,A-1))
      mu_d <- c(((40/1000)/365.25), ((46/1000)/365.25), ((6/1000)/365.25),((6/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((150/1000)/365.25))
      vaccination <- (vac_cov_2017/365.25) * vac_eff
      treatment <- prop_treat_2017
    } else if ((365.25*38) && t < (365.25*44)) {  #birth rate, age-specific death rates, vaccination rates and proportion treated between 2023 and 2029
      mu_b <- c(((1020/1000)/365.25), rep(0,A-1))
      mu_d <- c(((25/1000)/365.25), ((28/1000)/365.25), ((3/1000)/365.25), ((3/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((150/1000)/365.25))
      vaccination <- (vac_cov_2023/365.25) * vac_eff
      treatment <- prop_treat_2023
    } else if ((365.25*44) && t < (365.25*50)) {  #birth rate, age-specific death rates, vaccination rates and proportion treated between 2029 and 2035
      mu_b <- c(((1020/1000)/365.25), rep(0,A-1))
      mu_d <- c(((25/1000)/365.25), ((28/1000)/365.25), ((3/1000)/365.25), ((3/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((150/1000)/365.25))
      vaccination <- (vac_cov_2029/365.25) * vac_eff
      treatment <- prop_treat_2029
    } else if ((365.25*50) && t < (365.25*56)) {  #birth rate, age-specific death rates, vaccination rates and proportion treated between 2035 and 2041
      mu_b <- c(((1020/1000)/365.25), rep(0,A-1))
      mu_d <- c(((25/1000)/365.25), ((28/1000)/365.25), ((3/1000)/365.25), ((3/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((150/1000)/365.25))
      vaccination <- (vac_cov_2035/365.25) * vac_eff
      treatment <- prop_treat_2035 
    } else if ((365.25*56) && t < (365.25*62)) {  #birth rate, age-specific death rates, vaccination rates and proportion treated between 2041 and 2047
      mu_b <- c(((1020/1000)/365.25), rep(0,A-1))
      mu_d <- c(((25/1000)/365.25), ((28/1000)/365.25), ((3/1000)/365.25), ((3/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((150/1000)/365.25))
      vaccination <- (vac_cov_2041/365.25) * vac_eff
      treatment <- prop_treat_2041
    } else if ((365.25*62) && t < (365.25*68)) {  #birth rate, age-specific death rates, vaccination rates and proportion treated between 2047 and 2053
      mu_b <- c(((1020/1000)/365.25), rep(0,A-1))
      mu_d <- c(((25/1000)/365.25), ((28/1000)/365.25), ((3/1000)/365.25), ((3/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((150/1000)/365.25))
      vaccination <- (vac_cov_2047/365.25) * vac_eff
      treatment <- prop_treat_2047
    } else if ((365.25*38) && t < (365.25*44)) {  #birth rate, age-specific death rates, vaccination rates and proportion treated between 2053 and 2059
      mu_b <- c(((1020/1000)/365.25), rep(0,A-1))
      mu_d <- c(((25/1000)/365.25), ((28/1000)/365.25), ((3/1000)/365.25), ((3/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((150/1000)/365.25))
      vaccination <- (vac_cov_2053/365.25) * vac_eff
      treatment <- prop_treat_2053
    } else if ((365.25*38) && t < (365.25*44)) {  #birth rate, age-specific death rates, vaccination rates and proportion treated between 2059 and 2065
      mu_b <- c(((1020/1000)/365.25), rep(0,A-1))
      mu_d <- c(((25/1000)/365.25), ((28/1000)/365.25), ((3/1000)/365.25), ((3/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((150/1000)/365.25))
      vaccination <- (vac_cov_2059/365.25) * vac_eff
      treatment <- prop_treat_2059
    } else if ((365.25*38) && t < (365.25*44)) {  #birth rate, age-specific death rates, vaccination rates and proportion treated between 2065 and 2071
      mu_b <- c(((1020/1000)/365.25), rep(0,A-1))
      mu_d <- c(((25/1000)/365.25), ((28/1000)/365.25), ((3/1000)/365.25), ((3/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((150/1000)/365.25))
      vaccination <- (vac_cov_2065/365.25) * vac_eff
      treatment <- prop_treat_2065
    } else if ((365.25*38) && t < (365.25*44)) {  #birth rate, age-specific death rates, vaccination rates and proportion treated between 2071 and 2077
      mu_b <- c(((1020/1000)/365.25), rep(0,A-1))
      mu_d <- c(((25/1000)/365.25), ((28/1000)/365.25), ((3/1000)/365.25), ((3/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((150/1000)/365.25))
      vaccination <- (vac_cov_2071/365.25) * vac_eff
      treatment <- prop_treat_2071
    } else if ((365.25*38) && t < (365.25*44)) {  #birth rate, age-specific death rates, vaccination rates and proportion treated between 2077 and 2083
      mu_b <- c(((1020/1000)/365.25), rep(0,A-1))
      mu_d <- c(((25/1000)/365.25), ((28/1000)/365.25), ((3/1000)/365.25), ((3/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((150/1000)/365.25))
      vaccination <- (vac_cov_2077/365.25) * vac_eff
      treatment <- prop_treat_2077
    } else if ((365.25*38) && t < (365.25*44)) {  #birth rate, age-specific death rates, vaccination rates and proportion treated between 2083 and 2089
      mu_b <- c(((1020/1000)/365.25), rep(0,A-1))
      mu_d <- c(((25/1000)/365.25), ((28/1000)/365.25), ((3/1000)/365.25), ((3/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((150/1000)/365.25))
      vaccination <- (vac_cov_2083/365.25) * vac_eff
      treatment <- prop_treat_2083
    } else if ((365.25*38) && t < (365.25*44)) {  #birth rate, age-specific death rates, vaccination rates and proportion treated between 2089 and 2095
      mu_b <- c(((1020/1000)/365.25), rep(0,A-1))
      mu_d <- c(((25/1000)/365.25), ((28/1000)/365.25), ((3/1000)/365.25), ((3/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((150/1000)/365.25))
      vaccination <- (vac_cov_2089/365.25) * vac_eff
      treatment <- prop_treat_2089
    } else {  #birth rate, age-specific death rates, vaccination rates and proportion treated between 2089 and 2095
      mu_b <- c(((1020/1000)/365.25), rep(0,A-1))
      mu_d <- c(((25/1000)/365.25), ((28/1000)/365.25), ((3/1000)/365.25), ((3/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((2.7/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((22.6/1000)/365.25), ((150/1000)/365.25))
      vaccination <- (vac_cov_2095/365.25) * vac_eff
      treatment <- prop_treat_2095
    }
    
    #Dynamically calculated parameters 
    Infectious <- zetaAcu * (Acu_m + Acu_f) + Chron1_m + Chron1_f  + Chron2_m + Chron2_f + Chron3_m + Chron3_f + Chron4_m + Chron4_f
    lambda <- (p_infect * (contacts_home + contacts_school + contacts_sex)) %*% (Infectious/P_total)
    HBsAG_prev_15to49f <- (sum(Chron1_f[5:10]) + sum(Chron2_f[5:10]) + sum(Chron3_f[5:10]) + sum(Chron4_f[5:10])) / sum(P_f[5:10]) #prevalence of chronically infected women in the age group 15 - 49 (ages when women are likely to become pregnant)
    Susceptible_birth_prop <- ((1 - HBsAG_prev_15to49f) + (HBsAG_prev_15to49f * (1 - p_MTCT))) 
    Acute_birth_prop <- (HBsAG_prev_15to49f * p_MTCT)
    vaccination_vector <- c(vaccination, rep(0, (A -1))) #only people in the first age group are being vaccinated 
  
    # Differential Equations for Males
    dS_m <- Susceptible_birth_prop * prop_m * mu_b * P_total - (lambda + vaccination_vector + mu_d) * S_m + ageing_matrix %*% S_m
    dAcu_m <- Acute_birth_prop * prop_m * mu_b * P_total + lambda * S_m - (gamma1 + gamma2 + vaccination_vector + mu_d) * Acu_m + ageing_matrix %*% Acu_m 
    dChron1_m <- p_chron1 * gamma1 * Acu_m - (r12 + ((1 - treatment) * epsilon1) + (treatment * epsilon1_treated) + mu_d) * Chron1_m + ageing_matrix %*% Chron1_m  
    dChron2_m <- p_chron2 * gamma2 * Acu_m + r12 * Chron1_m + r32 * Chron3_m - (((1 - treatment) * r23) + (treatment * r23_treated) + ((1 - treatment) * delta2_m) + (treatment * delta2_m_treated) + ((1 - treatment) * epsilon2_m) + (treatment * epsilon2_m_treated) + rho_c2 + mu_d) * Chron2_m + ageing_matrix %*% Chron2_m  
    dChron3_m <- ((1 - treatment) * r23) * Chron2_m + (treatment * r23_treated) * Chron2_m - (((1 - treatment) * r34) + (treatment * r34_treated) + r32 + ((1 - treatment) * epsilon3_m) + (treatment * epsilon3_m_treated) + rho_c3 + mu_d) * Chron3_m + ageing_matrix %*% Chron3_m 
    dChron4_m <- ((1 - treatment) * r34) * Chron3_m + (treatment * r34_treated) * Chron3_m - (((1 - treatment) * delta4_m) + (treatment * delta4_m_treated) + ((1 - treatment) * epsilon4_m) + (treatment * epsilon4_m_treated) + rho_c4 + mu_d) * Chron4_m + ageing_matrix %*% Chron4_m 
    dCompC_m <- ((1 - treatment) * delta2_m) * Chron2_m + (treatment * delta2_m_treated) * Chron2_m + ((1 - treatment) * delta4_m) * Chron4_m + (treatment * delta4_m_treated) * Chron4_m - (((1 - treatment) * alpha) + (treatment * alpha_treated) + theta_m + mu_d) * CompC_m + ageing_matrix %*% CompC_m
    dDeCompC_m <- ((1 - treatment) * alpha) * CompC_m + (treatment * alpha_treated) * CompC_m - (((1 - treatment) * kappa) + (treatment * kappa_treated) + mu_d + muDCC) * DeCompC_m + ageing_matrix %*% DeCompC_m  
    dHCC_m <- ((1 - treatment) * epsilon1) * Chron1_m + (treatment * epsilon1_treated) * Chron1_m + ((1 - treatment) * epsilon2_m) * Chron2_m + (treatment * epsilon2_m_treated) * Chron2_m + ((1 - treatment) * epsilon3_m) * Chron3_m + (treatment * epsilon3_m_treated) * Chron3_m + ((1 - treatment) * epsilon4_m) * Chron4_m + (treatment * epsilon4_m_treated) * Chron4_m + theta_m * CompC_m + ((1 - treatment) * kappa) * DeCompC_m + (treatment * kappa_treated) * DeCompC_m - (mu_d + muHCC) * HCC_m + ageing_matrix %*% HCC_m 
    dR_m <- ((1 - p_chron1) * gamma1 + (1 - p_chron2) * gamma2) * Acu_m + rho_c2 * Chron2_m + rho_c3 * Chron3_m + rho_c4 * Chron4_m - mu_d * R_m + ageing_matrix %*% R_m  
    dVac_m <- vaccination_vector * S_m + vaccination_vector * Acu_m - mu_d * Vac_m + ageing_matrix %*% Vac_m
    
    # Differential Equations for Females
    dS_f <- Susceptible_birth_prop * prop_f * mu_b * P_total - (lambda + vaccination_vector + mu_d) * S_f + ageing_matrix %*% S_f
    dAcu_f <- Acute_birth_prop * prop_f * mu_b * P_total + lambda * S_f - (gamma1 + gamma2 + vaccination_vector + mu_d) * Acu_f + ageing_matrix %*% Acu_f 
    dChron1_f <- p_chron1 * gamma1 * Acu_f - (r12 + ((1 - treatment) * epsilon1) + (treatment * epsilon1_treated) + mu_d) * Chron1_f + ageing_matrix %*% Chron1_f  
    dChron2_f <- p_chron2 * gamma2 * Acu_f + r12 * Chron1_f + r32 * Chron3_f - (((1 - treatment) * r23) + (treatment * r23_treated) + ((1 - treatment) * delta2_f) + (treatment * delta2_f_treated) + ((1 - treatment) * epsilon2_f) + (treatment * epsilon2_f_treated) + rho_c2 + mu_d) * Chron2_f + ageing_matrix %*% Chron2_f  
    dChron3_f <- ((1 - treatment) * r23) * Chron2_f + (treatment * r23_treated) * Chron2_f - (((1 - treatment) * r34) + (treatment * r34_treated) + r32 + ((1 - treatment) * epsilon3_f) + (treatment * epsilon3_f_treated) + rho_c3 + mu_d) * Chron3_f + ageing_matrix %*% Chron3_f 
    dChron4_f <- ((1 - treatment) * r34) * Chron3_f + (treatment * r34_treated) * Chron3_f - (((1 - treatment) * delta4_f) + (treatment * delta4_f_treated) + ((1 - treatment) * epsilon4_f) + (treatment * epsilon4_f_treated) + rho_c4 + mu_d) * Chron4_f + ageing_matrix %*% Chron4_f 
    dCompC_f <- ((1 - treatment) * delta2_f) * Chron2_f + (treatment * delta2_f_treated) * Chron2_f + ((1 - treatment) * delta4_f) * Chron4_f + (treatment * delta4_f_treated) * Chron4_f - (((1 - treatment) * alpha) + (treatment * alpha_treated) + theta_f + mu_d) * CompC_f + ageing_matrix %*% CompC_f
    dDeCompC_f <- ((1 - treatment) * alpha) * CompC_f + (treatment * alpha_treated) * CompC_f - (((1 - treatment) * kappa) + (treatment * kappa_treated) + mu_d + muDCC) * DeCompC_f + ageing_matrix %*% DeCompC_f  
    dHCC_f <- ((1 - treatment) * epsilon1) * Chron1_f + (treatment * epsilon1_treated) * Chron1_f + ((1 - treatment) * epsilon2_f) * Chron2_f + (treatment * epsilon2_f_treated) * Chron2_f + ((1 - treatment) * epsilon3_f) * Chron3_f + (treatment * epsilon3_f_treated) * Chron3_f + ((1 - treatment) * epsilon4_f) * Chron4_f + (treatment * epsilon4_f_treated) * Chron4_f + theta_f * CompC_f + ((1 - treatment) * kappa) * DeCompC_f + (treatment * kappa_treated) * DeCompC_f - (mu_d + muHCC) * HCC_f + ageing_matrix %*% HCC_f 
    dR_f <- ((1 - p_chron1) * gamma1 + (1 - p_chron2) * gamma2) * Acu_f + rho_c2 * Chron2_f + rho_c3 * Chron3_f + rho_c4 * Chron4_f - mu_d * R_f + ageing_matrix %*% R_f  
    dVac_f <- vaccination_vector * S_f + vaccination_vector * Acu_f - mu_d * Vac_f + ageing_matrix %*% Vac_f
    
    
    # Counters for tracking other dynamics within the model (summing the number of males and females in all cumulative counters)
    dCMort <- muDCC * (DeCompC_m + DeCompC_f) + muHCC * (HCC_m + HCC_f) 
    dCInc_H <- lambda * (S_m + S_f)
    dCInc_V <- HBsAG_prev_15to49f * p_MTCT * mu_b * P_total
    dCChron <- gamma1 * (Acu_m + Acu_f) + gamma2 * (Acu_m + Acu_f)
    dCCompC <- delta2_m * Chron2_m + delta2_f * Chron2_f + delta4_m * Chron4_m + delta4_f * Chron4_f
    dCDeCompC <- (alpha + alpha_treated) * (CompC_m + CompC_f)
    dCHCC <- (epsilon1 + epsilon1_treated) * (Chron1_m + Chron1_f) + ((epsilon2_m + epsilon2_m_treated) * Chron2_m) + ((epsilon2_f + epsilon2_f_treated)  * Chron2_f) + ((epsilon3_m + epsilon3_m_treated) * Chron3_m) + ((epsilon3_f + epsilon3_f_treated)  * Chron3_f) + ((epsilon4_m + epsilon4_m_treated)  * Chron4_m) + ((epsilon4_f + epsilon4_f_treated)  * Chron4_f) + (theta_m * CompC_m) + (theta_f * CompC_f) + (kappa + kappa_treated) * (DeCompC_m + DeCompC_f)
    dCVac <- vaccination_vector * (S_m + S_f) + vaccination_vector * (Acu_m + Acu_f)
    
    # Outputs an ordered list of derivatives of the state variables
    list(c(dS_m, dAcu_m, dChron1_m, dChron2_m, dChron3_m,
           dChron4_m, dCompC_m, dDeCompC_m, dHCC_m, dR_m, dVac_m,  
           dS_f, dAcu_f, dChron1_f, dChron2_f, dChron3_f, 
           dChron4_f, dCompC_f, dDeCompC_f, dHCC_f, dR_f, dVac_f,
           dCMort, dCInc_H, dCInc_V, dCChron, dCCompC, dCDeCompC, dCHCC, dCVac))
  })
}

######################################################################
# MODEL TIMES 
######################################################################
times_laopdr_HBV_model_bc1 <- seq(0, (115 * 365.25), 1)

######################################################################
# ORDINARY DIFFERENTIAL EQUATION SOLVER
######################################################################
laopdr_HBV_output_bc1 <- ode(times = times_laopdr_HBV_model_bc1, 
                             y = init_laopdr_HBV_model_bc1, 
                             func = laopdr_HBV_model_bc1, 
                             parms = parms)

df_laopdr_HBV_output_bc1 <- as.data.frame(laopdr_HBV_output_bc1) #store in data frame
#View(df_laopdr_HBV_output_bc1)

# Rename columns in data frame
health_states <- c(
  "S_m", "Acu_m", "Chron1_m", "Chron2_m", "Chron3_m", 
  "Chron4_m", "CompC_m", "DeCcompC_m", "HCC_m", "R_m", "Vac_m",
  "S_f", "Acu_f", "Chron1_f", "Chron2_f", "Chron3_f",
  "Chron4_f", "CompC_f", "DeCompC_f", "HCC_f", "R_f", "Vac_f",
  "CMort", "CInc_H", "CInc_V", "CChron", "CCompC", "CDecompC", "CHCC", "CVac"
)

# column_names <- c("time", unlist(sapply(health_states, function(hs) paste(age_group_names, hs, sep = "_"))))
column_names <- c("time", unlist(sapply(health_states, function(hs) paste(hs, age_group_names, sep = "_")))) # Generate the column names
colnames(df_laopdr_HBV_output_bc1) <- unlist(column_names) # Rename columns in df_combined
columns_to_sum <- colnames(df_laopdr_HBV_output_bc1)[2:length(df_laopdr_HBV_output_bc1)] #function that sums specified columns

#View(df_laopdr_HBV_output_bc1)


######################################################################
#COMPUTING ANNUAL METRICS
######################################################################

# Calculate Total Population Population by Age Group and by Health State
df_laopdr_HBV_output_bc1 <- df_laopdr_HBV_output_bc1 %>%
  mutate(
    S_total = rowSums(select(., columns_to_sum[c((1:A), ((11 * A + 1) : (12 * A)))])),# summing the number of males and females in each health state
    Acu_total= rowSums(select(., columns_to_sum[c(((A + 1) : (2 * A) ), ((12 * A + 1) : (13 * A)))])),
    Chron1_total = rowSums(select(., columns_to_sum[c(((2 * A + 1) : (3 * A) ), ((13 * A + 1) : (14 * A)))])),
    Chron2_total = rowSums(select(., columns_to_sum[c(((3 * A + 1) : (4 * A) ), ((14 * A + 1) : (15 * A)))])),
    Chron3_total = rowSums(select(., columns_to_sum[c(((4 * A + 1) : (5 * A) ), ((15 * A + 1) : (16 * A)))])),
    Chron4_total = rowSums(select(., columns_to_sum[c(((5 * A + 1) : (6 * A) ), ((16 * A + 1) : (17 * A)))])),
    CompC_total = rowSums(select(., columns_to_sum[c(((6 * A + 1) : (7 * A) ), ((17 * A + 1) : (18 * A)))])),
    DecompC_total = rowSums(select(., columns_to_sum[c(((7 * A + 1) : (8 * A) ), ((18 * A + 1) : (19 * A)))])),
    HCC_total= rowSums(select(., columns_to_sum[c(((8 * A + 1) : (9 * A) ), ((19 * A + 1) : (20 * A)))])),
    R_total = rowSums(select(., columns_to_sum[c((( 9 * A + 1) : (10 * A) ), ((20 * A + 1) : (21 * A)))])),
    Vac_total = rowSums(select(., columns_to_sum[c((( 10 * A + 1) : (11 * A) ), ((21 * A + 1) : (22 * A)))])),
    CMort_total = rowSums(select(., columns_to_sum[(22 * A + 1) : (23 * A)])), # males and females in cumulative counters have already been summed
    CInc_H_total = rowSums(select(., columns_to_sum[(23 * A + 1) : (24 * A)])),
    CInc_V_total = rowSums(select(., columns_to_sum[(24 * A + 1) : (25 * A)])),
    CInc_total = rowSums(select(., columns_to_sum[(23 * A + 1) : (25 * A)])), # keeping track of horizontal and vertical transmission total to track percentage coming from each
    CChron_total = rowSums(select(., columns_to_sum[(25 * A + 1) : (26 * A)])),
    CCompC_total = rowSums(select(., columns_to_sum[(26 * A + 1) : (27 * A)])),
    CDecompC_total = rowSums(select(., columns_to_sum[(27 * A + 1) : (28 * A)])),
    CHCC_total = rowSums(select(., columns_to_sum[(28 * A + 1) : (29 * A)])),
    CVac_total = rowSums(select(., columns_to_sum[(29 * A + 1) : (30 * A)])),
    years = 1985 + floor(time / 365.25),
    Prev_total = (rowSums(select(., columns_to_sum[c(((2 * A + 1) : (3 * A) ), ((13 * A + 1) : (14 * A)), ((3 * A + 1) : (4 * A) ), ((14 * A + 1) : (15 * A)), ((4 * A + 1) : (5 * A) ), ((15 * A + 1) : (16 * A)), ((5 * A + 1) : (6 * A)), ((16 * A + 1) : (17 * A)))]))) / (rowSums(select(., all_of(columns_to_sum[1 : (22 * A)])))),
    Prev_0to4yrs = rowSums(select(., columns_to_sum[c((2 * A + 1), (2 * A + 2), (3 * A + 1), (3 * A + 2), (4 * A + 1), (4 * A + 2), (5 * A + 1), (5 * A + 2), (6 * A + 1), (6 * A + 2), (13 * A + 1), (13 * A + 2), (14 * A + 1), (14 * A + 2), (15 * A + 1), (15 * A + 2), (16 * A + 1), (16 * A + 2), (17 * A + 1), (17 * A) + 2)]))  / #combine 0 - 4 years (0 - 6 months + 7 months - 4 years)
      (rowSums(select(., columns_to_sum[1])) + rowSums(select(., columns_to_sum[2])) + # Combine 0 - 4 years (0 - 6 months + 7 months - 4 years)
         rowSums(select(., columns_to_sum[A + 1])) + rowSums(select(., columns_to_sum[A + 2])) +
         rowSums(select(., columns_to_sum[2 * A + 1])) + rowSums(select(., columns_to_sum[2 * A + 2])) +
         rowSums(select(., columns_to_sum[3 * A + 1])) + rowSums(select(., columns_to_sum[3 * A + 2])) +
         rowSums(select(., columns_to_sum[4 * A + 1])) + rowSums(select(., columns_to_sum[4 * A + 2])) +
         rowSums(select(., columns_to_sum[5 * A + 1])) + rowSums(select(., columns_to_sum[5 * A + 2])) +
         rowSums(select(., columns_to_sum[6 * A + 1])) + rowSums(select(., columns_to_sum[6 * A + 2])) +
         rowSums(select(., columns_to_sum[7 * A + 1])) + rowSums(select(., columns_to_sum[7 * A + 2])) +
         rowSums(select(., columns_to_sum[8 * A + 1])) + rowSums(select(., columns_to_sum[8 * A + 2])) +
         rowSums(select(., columns_to_sum[9 * A + 1])) + rowSums(select(., columns_to_sum[9 * A + 2])) +
         rowSums(select(., columns_to_sum[10 * A + 1])) + rowSums(select(., columns_to_sum[10 * A + 2])) +
         rowSums(select(., columns_to_sum[11 * A + 1])) + rowSums(select(., columns_to_sum[11 * A + 2])) +
         rowSums(select(., columns_to_sum[12 * A + 1])) + rowSums(select(., columns_to_sum[12 * A + 2])) +
         rowSums(select(., columns_to_sum[13 * A + 1])) + rowSums(select(., columns_to_sum[13 * A + 2])) +
         rowSums(select(., columns_to_sum[14 * A + 1])) + rowSums(select(., columns_to_sum[14 * A + 2])) +
         rowSums(select(., columns_to_sum[15 * A + 1])) + rowSums(select(., columns_to_sum[15 * A + 2])) +
         rowSums(select(., columns_to_sum[16 * A + 1])) + rowSums(select(., columns_to_sum[16 * A + 2])) +
         rowSums(select(., columns_to_sum[17 * A + 1])) + rowSums(select(., columns_to_sum[17 * A + 2])) +
         rowSums(select(., columns_to_sum[18 * A + 1])) + rowSums(select(., columns_to_sum[18 * A + 2])) +
         rowSums(select(., columns_to_sum[19 * A + 1])) + rowSums(select(., columns_to_sum[19 * A + 2])) +
         rowSums(select(., columns_to_sum[20 * A + 1])) + rowSums(select(., columns_to_sum[20 * A + 2])) +
         rowSums(select(., columns_to_sum[21 * A + 1])) + rowSums(select(., columns_to_sum[21 * A + 2])) +
         rowSums(select(., columns_to_sum[22 * A + 1])) + rowSums(select(., columns_to_sum[22 * A + 2]))),
    Prev_5to9yrs = rowSums(select(., columns_to_sum[c((2 * A + 3), (3 * A + 3), (4 * A + 3), (5 * A + 3), (6 * A + 3), (13 * A + 3), (14 * A + 3), (15 * A + 3), (16 * A + 3), (17 * A + 3))]))  /
      (rowSums(select(., columns_to_sum[3])) + rowSums(select(., columns_to_sum[A + 3])) +
         rowSums(select(., columns_to_sum[2 * A + 3])) + rowSums(select(., columns_to_sum[3 * A + 3])) +
         rowSums(select(., columns_to_sum[4 * A + 3])) + rowSums(select(., columns_to_sum[5 * A + 3])) +
         rowSums(select(., columns_to_sum[6 * A + 3])) + rowSums(select(., columns_to_sum[7 * A + 3])) +
         rowSums(select(., columns_to_sum[8 * A + 3])) + rowSums(select(., columns_to_sum[9 * A + 3])) +
         rowSums(select(., columns_to_sum[10 * A + 3])) + rowSums(select(., columns_to_sum[11 * A + 3])) +
         rowSums(select(., columns_to_sum[12 * A + 3])) + rowSums(select(., columns_to_sum[13 * A + 3])) +
         rowSums(select(., columns_to_sum[14 * A + 3])) + rowSums(select(., columns_to_sum[15 * A + 3])) +
         rowSums(select(., columns_to_sum[16 * A + 3])) + rowSums(select(., columns_to_sum[17 * A + 3])) +
         rowSums(select(., columns_to_sum[18 * A + 3])) + rowSums(select(., columns_to_sum[19 * A + 3])) +
         rowSums(select(., columns_to_sum[20 * A + 3])) + rowSums(select(., columns_to_sum[21 * A + 3])) +
         rowSums(select(., columns_to_sum[22 * A + 3]))),
    Prev_10to14yrs = rowSums(select(., columns_to_sum[c((2 * A + 4), (3 * A + 4), (4 * A + 4), (5 * A + 4), (6 * A + 4), (13 * A + 4), (14 * A + 4), (15 * A + 4), (16 * A + 4), (17 * A + 4))]))  /
      (rowSums(select(., columns_to_sum[4])) + rowSums(select(., columns_to_sum[A + 4])) +
         rowSums(select(., columns_to_sum[2 * A + 4])) + rowSums(select(., columns_to_sum[3 * A + 4])) +
         rowSums(select(., columns_to_sum[4 * A + 4])) + rowSums(select(., columns_to_sum[5 * A + 4])) +
         rowSums(select(., columns_to_sum[6 * A + 4])) + rowSums(select(., columns_to_sum[7 * A + 4])) +
         rowSums(select(., columns_to_sum[8 * A + 4])) + rowSums(select(., columns_to_sum[9 * A + 4])) +
         rowSums(select(., columns_to_sum[10 * A + 4])) + rowSums(select(., columns_to_sum[11 * A + 4])) +
         rowSums(select(., columns_to_sum[12 * A + 4])) + rowSums(select(., columns_to_sum[13 * A + 4])) +
         rowSums(select(., columns_to_sum[14 * A + 4])) + rowSums(select(., columns_to_sum[15 * A + 4])) +
         rowSums(select(., columns_to_sum[16 * A + 4])) + rowSums(select(., columns_to_sum[17 * A + 4])) +
         rowSums(select(., columns_to_sum[18 * A + 4])) + rowSums(select(., columns_to_sum[19 * A + 4])) +
         rowSums(select(., columns_to_sum[20 * A + 4])) + rowSums(select(., columns_to_sum[21 * A + 4])) +
         rowSums(select(., columns_to_sum[22 * A + 4]))),
    Prev_15to19yrs = rowSums(select(., columns_to_sum[c((2 * A + 5), (3 * A + 5), (4 * A + 5), (5 * A + 5), (6 * A + 5), (13 * A + 5), (14 * A + 5), (15 * A + 5), (16 * A + 5), (17 * A + 5))]))  /
      (rowSums(select(., columns_to_sum[5])) + rowSums(select(., columns_to_sum[A + 5])) +
         rowSums(select(., columns_to_sum[2 * A + 5])) + rowSums(select(., columns_to_sum[3 * A + 5])) +
         rowSums(select(., columns_to_sum[4 * A + 5])) + rowSums(select(., columns_to_sum[5 * A + 5])) +
         rowSums(select(., columns_to_sum[6 * A + 5])) + rowSums(select(., columns_to_sum[7 * A + 5])) +
         rowSums(select(., columns_to_sum[8 * A + 5])) + rowSums(select(., columns_to_sum[9 * A + 5])) +
         rowSums(select(., columns_to_sum[10 * A + 5])) + rowSums(select(., columns_to_sum[11 * A + 5])) +
         rowSums(select(., columns_to_sum[12 * A + 5])) + rowSums(select(., columns_to_sum[13 * A + 5])) +
         rowSums(select(., columns_to_sum[14 * A + 5])) + rowSums(select(., columns_to_sum[15 * A + 5])) +
         rowSums(select(., columns_to_sum[16 * A + 5])) + rowSums(select(., columns_to_sum[17 * A + 5])) +
         rowSums(select(., columns_to_sum[18 * A + 5])) + rowSums(select(., columns_to_sum[19 * A + 5])) +
         rowSums(select(., columns_to_sum[20 * A + 5])) + rowSums(select(., columns_to_sum[21 * A + 5])) +
         rowSums(select(., columns_to_sum[22 * A + 5]))),
    Prev_20to24yrs = rowSums(select(., columns_to_sum[c((2 * A + 6), (3 * A + 6), (4 * A + 6), (5 * A + 6), (6 * A + 6), (13 * A + 6), (14 * A + 6), (15 * A + 6), (16 * A + 6), (17 * A + 6))]))  /
      (rowSums(select(., columns_to_sum[6])) + rowSums(select(., columns_to_sum[A + 6])) +
         rowSums(select(., columns_to_sum[2 * A + 6])) + rowSums(select(., columns_to_sum[3 * A + 6])) +
         rowSums(select(., columns_to_sum[4 * A + 6])) + rowSums(select(., columns_to_sum[5 * A + 6])) +
         rowSums(select(., columns_to_sum[6 * A + 6])) + rowSums(select(., columns_to_sum[7 * A + 6])) +
         rowSums(select(., columns_to_sum[8 * A + 6])) + rowSums(select(., columns_to_sum[9 * A + 6])) +
         rowSums(select(., columns_to_sum[10 * A + 6])) + rowSums(select(., columns_to_sum[11 * A + 6])) +
         rowSums(select(., columns_to_sum[12 * A + 6])) + rowSums(select(., columns_to_sum[13 * A + 6])) +
         rowSums(select(., columns_to_sum[14 * A + 6])) + rowSums(select(., columns_to_sum[15 * A + 6])) +
         rowSums(select(., columns_to_sum[16 * A + 6])) + rowSums(select(., columns_to_sum[17 * A + 6])) +
         rowSums(select(., columns_to_sum[18 * A + 6])) + rowSums(select(., columns_to_sum[19 * A + 6])) +
         rowSums(select(., columns_to_sum[20 * A + 6])) + rowSums(select(., columns_to_sum[21 * A + 6])) +
         rowSums(select(., columns_to_sum[22 * A + 6]))),
    Prev_25to29yrs = rowSums(select(., columns_to_sum[c((2 * A + 7), (3 * A + 7), (4 * A + 7), (5 * A + 7), (6 * A + 7), (13 * A + 7), (14 * A + 7), (15 * A + 7), (16 * A + 7), (17 * A + 7))]))  /
      (rowSums(select(., columns_to_sum[7])) + rowSums(select(., columns_to_sum[A + 7])) +
         rowSums(select(., columns_to_sum[2 * A + 7])) + rowSums(select(., columns_to_sum[3 * A + 7])) +
         rowSums(select(., columns_to_sum[4 * A + 7])) + rowSums(select(., columns_to_sum[5 * A + 7])) +
         rowSums(select(., columns_to_sum[6 * A + 7])) + rowSums(select(., columns_to_sum[7 * A + 7])) +
         rowSums(select(., columns_to_sum[8 * A + 7])) + rowSums(select(., columns_to_sum[9 * A + 7])) +
         rowSums(select(., columns_to_sum[10 * A + 7])) + rowSums(select(., columns_to_sum[11 * A + 7])) +
         rowSums(select(., columns_to_sum[12 * A + 7])) + rowSums(select(., columns_to_sum[13 * A + 7])) +
         rowSums(select(., columns_to_sum[14 * A + 7])) + rowSums(select(., columns_to_sum[15 * A + 7])) +
         rowSums(select(., columns_to_sum[16 * A + 7])) + rowSums(select(., columns_to_sum[17 * A + 7])) +
         rowSums(select(., columns_to_sum[18 * A + 7])) + rowSums(select(., columns_to_sum[19 * A + 7])) +
         rowSums(select(., columns_to_sum[20 * A + 7])) + rowSums(select(., columns_to_sum[21 * A + 7])) +
         rowSums(select(., columns_to_sum[22 * A + 7]))),
    Prev_30to34yrs = rowSums(select(., columns_to_sum[c((2 * A + 8), (3 * A + 8), (4 * A + 8), (5 * A + 8), (6 * A + 8), (13 * A + 8), (14 * A + 8), (15 * A + 8), (16 * A + 8), (17 * A + 8))]))  /
      (rowSums(select(., columns_to_sum[8])) + rowSums(select(., columns_to_sum[A + 8])) +
         rowSums(select(., columns_to_sum[2 * A + 8])) + rowSums(select(., columns_to_sum[3 * A + 8])) +
         rowSums(select(., columns_to_sum[4 * A + 8])) + rowSums(select(., columns_to_sum[5 * A + 8])) +
         rowSums(select(., columns_to_sum[6 * A + 8])) + rowSums(select(., columns_to_sum[7 * A + 8])) +
         rowSums(select(., columns_to_sum[8 * A + 8])) + rowSums(select(., columns_to_sum[9 * A + 8])) +
         rowSums(select(., columns_to_sum[10 * A + 8])) + rowSums(select(., columns_to_sum[11 * A + 8])) +
         rowSums(select(., columns_to_sum[12 * A + 8])) + rowSums(select(., columns_to_sum[13 * A + 8])) +
         rowSums(select(., columns_to_sum[14 * A + 8])) + rowSums(select(., columns_to_sum[15 * A + 8])) +
         rowSums(select(., columns_to_sum[16 * A + 8])) + rowSums(select(., columns_to_sum[17 * A + 8])) +
         rowSums(select(., columns_to_sum[18 * A + 8])) + rowSums(select(., columns_to_sum[19 * A + 8])) +
         rowSums(select(., columns_to_sum[20 * A + 8])) + rowSums(select(., columns_to_sum[21 * A + 8])) +
         rowSums(select(., columns_to_sum[22 * A + 8]))),
    Prev_35to39yrs = rowSums(select(., columns_to_sum[c((2 * A + 9), (3 * A + 9), (4 * A + 9), (5 * A + 9), (6 * A + 9), (13 * A + 9), (14 * A + 9), (15 * A + 9), (16 * A + 9), (17 * A + 9))]))  /
      (rowSums(select(., columns_to_sum[9])) + rowSums(select(., columns_to_sum[A + 9])) +
         rowSums(select(., columns_to_sum[2 * A + 9])) + rowSums(select(., columns_to_sum[3 * A + 9])) +
         rowSums(select(., columns_to_sum[4 * A + 9])) + rowSums(select(., columns_to_sum[5 * A + 9])) +
         rowSums(select(., columns_to_sum[6 * A + 9])) + rowSums(select(., columns_to_sum[7 * A + 9])) +
         rowSums(select(., columns_to_sum[8 * A + 9])) + rowSums(select(., columns_to_sum[9 * A + 9])) +
         rowSums(select(., columns_to_sum[10 * A + 9])) + rowSums(select(., columns_to_sum[11 * A + 9])) +
         rowSums(select(., columns_to_sum[12 * A + 9])) + rowSums(select(., columns_to_sum[13 * A + 9])) +
         rowSums(select(., columns_to_sum[14 * A + 9])) + rowSums(select(., columns_to_sum[15 * A + 9])) +
         rowSums(select(., columns_to_sum[16 * A + 9])) + rowSums(select(., columns_to_sum[17 * A + 9])) +
         rowSums(select(., columns_to_sum[18 * A + 9])) + rowSums(select(., columns_to_sum[19 * A + 9])) +
         rowSums(select(., columns_to_sum[20 * A + 9])) + rowSums(select(., columns_to_sum[21 * A + 9])) +
         rowSums(select(., columns_to_sum[22 * A + 9]))),
    Prev_40to44yrs = rowSums(select(., columns_to_sum[c((2 * A + 10), (3 * A + 10), (4 * A + 10), (5 * A + 10), (6 * A + 10), (13 * A + 10), (14 * A + 10), (15 * A + 10), (16 * A + 10), (17 * A + 10))]))  /
      (rowSums(select(., columns_to_sum[10])) + rowSums(select(., columns_to_sum[A + 10])) +
         rowSums(select(., columns_to_sum[2 * A + 10])) + rowSums(select(., columns_to_sum[3 * A + 10])) +
         rowSums(select(., columns_to_sum[4 * A + 10])) + rowSums(select(., columns_to_sum[5 * A + 10])) +
         rowSums(select(., columns_to_sum[6 * A + 10])) + rowSums(select(., columns_to_sum[7 * A + 10])) +
         rowSums(select(., columns_to_sum[8 * A + 10])) + rowSums(select(., columns_to_sum[9 * A + 10])) +
         rowSums(select(., columns_to_sum[10 * A + 10])) + rowSums(select(., columns_to_sum[11 * A + 10])) +
         rowSums(select(., columns_to_sum[12 * A + 10])) + rowSums(select(., columns_to_sum[13 * A + 10])) +
         rowSums(select(., columns_to_sum[14 * A + 10])) + rowSums(select(., columns_to_sum[15 * A + 10])) +
         rowSums(select(., columns_to_sum[16 * A + 10])) + rowSums(select(., columns_to_sum[17 * A + 10])) +
         rowSums(select(., columns_to_sum[18 * A + 10])) + rowSums(select(., columns_to_sum[19 * A + 10])) +
         rowSums(select(., columns_to_sum[20 * A + 10])) + rowSums(select(., columns_to_sum[21 * A + 10])) +
         rowSums(select(., columns_to_sum[22 * A + 10]))),
    Prev_45to49yrs = rowSums(select(., columns_to_sum[c((2 * A + 11), (3 * A + 11), (4 * A + 11), (5 * A + 11), (6 * A + 11), (13 * A + 11), (14 * A + 11), (15 * A + 11), (16 * A + 11), (17 * A + 11))]))  /
      (rowSums(select(., columns_to_sum[11])) + rowSums(select(., columns_to_sum[A + 11])) +
         rowSums(select(., columns_to_sum[2 * A + 11])) + rowSums(select(., columns_to_sum[3 * A + 11])) +
         rowSums(select(., columns_to_sum[4 * A + 11])) + rowSums(select(., columns_to_sum[5 * A + 11])) +
         rowSums(select(., columns_to_sum[6 * A + 11])) + rowSums(select(., columns_to_sum[7 * A + 11])) +
         rowSums(select(., columns_to_sum[8 * A + 11])) + rowSums(select(., columns_to_sum[9 * A + 11])) +
         rowSums(select(., columns_to_sum[10 * A + 11])) + rowSums(select(., columns_to_sum[11 * A + 11])) +
         rowSums(select(., columns_to_sum[12 * A + 11])) + rowSums(select(., columns_to_sum[13 * A + 11])) +
         rowSums(select(., columns_to_sum[14 * A + 11])) + rowSums(select(., columns_to_sum[15 * A + 11])) +
         rowSums(select(., columns_to_sum[16 * A + 11])) + rowSums(select(., columns_to_sum[17 * A + 11])) +
         rowSums(select(., columns_to_sum[18 * A + 11])) + rowSums(select(., columns_to_sum[19 * A + 11])) +
         rowSums(select(., columns_to_sum[20 * A + 11])) + rowSums(select(., columns_to_sum[21 * A + 11])) +
         rowSums(select(., columns_to_sum[22 * A + 11]))),
    Prev_50to54yrs = rowSums(select(., columns_to_sum[c((2 * A + 12), (3 * A + 12), (4 * A + 12), (5 * A + 12), (6 * A + 12), (13 * A + 12), (14 * A + 12), (15 * A + 12), (16 * A + 12), (17 * A + 12))]))  /
      (rowSums(select(., columns_to_sum[12])) + rowSums(select(., columns_to_sum[A + 12])) +
         rowSums(select(., columns_to_sum[2 * A + 12])) + rowSums(select(., columns_to_sum[3 * A + 12])) +
         rowSums(select(., columns_to_sum[4 * A + 12])) + rowSums(select(., columns_to_sum[5 * A + 12])) +
         rowSums(select(., columns_to_sum[6 * A + 12])) + rowSums(select(., columns_to_sum[7 * A + 12])) +
         rowSums(select(., columns_to_sum[8 * A + 12])) + rowSums(select(., columns_to_sum[9 * A + 12])) +
         rowSums(select(., columns_to_sum[10 * A + 12])) + rowSums(select(., columns_to_sum[11 * A + 12])) +
         rowSums(select(., columns_to_sum[12 * A + 12])) + rowSums(select(., columns_to_sum[13 * A + 12])) +
         rowSums(select(., columns_to_sum[14 * A + 12])) + rowSums(select(., columns_to_sum[15 * A + 12])) +
         rowSums(select(., columns_to_sum[16 * A + 12])) + rowSums(select(., columns_to_sum[17 * A + 12])) +
         rowSums(select(., columns_to_sum[18 * A + 12])) + rowSums(select(., columns_to_sum[19 * A + 12])) +
         rowSums(select(., columns_to_sum[20 * A + 12])) + rowSums(select(., columns_to_sum[21 * A + 12])) +
         rowSums(select(., columns_to_sum[22 * A + 12]))),
    Prev_55to59yrs = rowSums(select(., columns_to_sum[c((2 * A + 13), (3 * A + 13), (4 * A + 13), (5 * A + 13), (6 * A + 13), (13 * A + 13), (14 * A + 13), (15 * A + 13), (16 * A + 13), (17 * A + 13))]))  /
      (rowSums(select(., columns_to_sum[13])) + rowSums(select(., columns_to_sum[A + 13])) +
         rowSums(select(., columns_to_sum[2 * A + 13])) + rowSums(select(., columns_to_sum[3 * A + 13])) +
         rowSums(select(., columns_to_sum[4 * A + 13])) + rowSums(select(., columns_to_sum[5 * A + 13])) +
         rowSums(select(., columns_to_sum[6 * A + 13])) + rowSums(select(., columns_to_sum[7 * A + 13])) +
         rowSums(select(., columns_to_sum[8 * A + 13])) + rowSums(select(., columns_to_sum[9 * A + 13])) +
         rowSums(select(., columns_to_sum[10 * A + 13])) + rowSums(select(., columns_to_sum[11 * A + 13])) +
         rowSums(select(., columns_to_sum[12 * A + 13])) + rowSums(select(., columns_to_sum[13 * A + 13])) +
         rowSums(select(., columns_to_sum[14 * A + 13])) + rowSums(select(., columns_to_sum[15 * A + 13])) +
         rowSums(select(., columns_to_sum[16 * A + 13])) + rowSums(select(., columns_to_sum[17 * A + 13])) +
         rowSums(select(., columns_to_sum[18 * A + 13])) + rowSums(select(., columns_to_sum[19 * A + 13])) +
         rowSums(select(., columns_to_sum[20 * A + 13])) + rowSums(select(., columns_to_sum[21 * A + 13])) +
         rowSums(select(., columns_to_sum[22 * A + 13]))),
    Prev_60to64yrs = rowSums(select(., columns_to_sum[c((2 * A + 14), (3 * A + 14), (4 * A + 14), (5 * A + 14), (6 * A + 14), (13 * A + 14), (14 * A + 14), (15 * A + 14), (16 * A + 14), (17 * A + 14))])) /
      (rowSums(select(., columns_to_sum[14])) + rowSums(select(., columns_to_sum[A + 14])) +
         rowSums(select(., columns_to_sum[2 * A + 14])) + rowSums(select(., columns_to_sum[3 * A + 14])) +
         rowSums(select(., columns_to_sum[4 * A + 14])) + rowSums(select(., columns_to_sum[5 * A + 14])) +
         rowSums(select(., columns_to_sum[6 * A + 14])) + rowSums(select(., columns_to_sum[7 * A + 14])) +
         rowSums(select(., columns_to_sum[8 * A + 14])) + rowSums(select(., columns_to_sum[9 * A + 14])) +
         rowSums(select(., columns_to_sum[10 * A + 14])) + rowSums(select(., columns_to_sum[11 * A + 14])) +
         rowSums(select(., columns_to_sum[12 * A + 14])) + rowSums(select(., columns_to_sum[13 * A + 14])) +
         rowSums(select(., columns_to_sum[14 * A + 14])) + rowSums(select(., columns_to_sum[15 * A + 14])) +
         rowSums(select(., columns_to_sum[16 *A + 14])) + rowSums(select(., columns_to_sum[17 * A + 14])) +
         rowSums(select(., columns_to_sum[18 * A + 14])) + rowSums(select(., columns_to_sum[19 * A + 14])) +
         rowSums(select(., columns_to_sum[20 * A + 14])) + rowSums(select(., columns_to_sum[21 * A + 14])) +
         rowSums(select(., columns_to_sum[22 * A + 14]))),
    Prev_65to69yrs = rowSums(select(., columns_to_sum[c((2 * A + 15), (3 * A + 15), (4 * A + 15), (5 * A + 15), (6 * A + 15), (13 * A + 15), (14 * A + 15), (15 * A + 15), (16 * A + 15), (17 * A + 15))]))  /
      (rowSums(select(., columns_to_sum[15])) + rowSums(select(., columns_to_sum[A + 15])) +
         rowSums(select(., columns_to_sum[2 * A + 15])) + rowSums(select(., columns_to_sum[3 * A + 15])) +
         rowSums(select(., columns_to_sum[4 * A + 15])) + rowSums(select(., columns_to_sum[5 * A + 15])) +
         rowSums(select(., columns_to_sum[6 * A + 15])) + rowSums(select(., columns_to_sum[7 * A + 15])) +
         rowSums(select(., columns_to_sum[8 * A + 15])) + rowSums(select(., columns_to_sum[9 * A + 15])) +
         rowSums(select(., columns_to_sum[10 * A + 15])) + rowSums(select(., columns_to_sum[11 * A + 15])) +
         rowSums(select(., columns_to_sum[12 * A + 15])) + rowSums(select(., columns_to_sum[13 * A + 15])) +
         rowSums(select(., columns_to_sum[14 * A + 15])) + rowSums(select(., columns_to_sum[15 * A + 15])) +
         rowSums(select(., columns_to_sum[16 * A + 15])) + rowSums(select(., columns_to_sum[17 * A + 15])) +
         rowSums(select(., columns_to_sum[18 * A + 15])) + rowSums(select(., columns_to_sum[19 * A + 15])) +
         rowSums(select(., columns_to_sum[20 * A + 15])) + rowSums(select(., columns_to_sum[21 * A + 15])) +
         rowSums(select(., columns_to_sum[22 * A + 15]))),
    Prev_70plusyrs = rowSums(select(., columns_to_sum[c((2 * A + 16), (3 * A + 16), (4 * A + 16), (5 * A + 16), (6 * A + 16), (13 * A + 16), (14 * A + 16), (15 * A + 16), (16 * A + 16), (17 * A + 16))])) /
      (rowSums(select(., columns_to_sum[16])) + rowSums(select(., columns_to_sum[A + 16])) +
         rowSums(select(., columns_to_sum[2 * A + 16])) + rowSums(select(., columns_to_sum[3 * A + 16])) +
         rowSums(select(., columns_to_sum[4 * A + 16])) + rowSums(select(., columns_to_sum[5 * A + 16])) +
         rowSums(select(., columns_to_sum[6 * A + 16])) + rowSums(select(., columns_to_sum[7 * A + 16])) +
         rowSums(select(., columns_to_sum[8 * A + 16])) + rowSums(select(., columns_to_sum[9 * A + 16])) +
         rowSums(select(., columns_to_sum[10 * A + 16])) + rowSums(select(., columns_to_sum[11 * A + 16])) +
         rowSums(select(., columns_to_sum[12 * A + 16])) + rowSums(select(., columns_to_sum[13 * A + 16])) +
         rowSums(select(., columns_to_sum[14 * A + 16])) + rowSums(select(., columns_to_sum[15 * A + 16])) +
         rowSums(select(., columns_to_sum[16 * A + 16])) + rowSums(select(., columns_to_sum[17 * A + 16])) +
         rowSums(select(., columns_to_sum[18 * A + 16])) + rowSums(select(., columns_to_sum[19 * A + 16])) +
         rowSums(select(., columns_to_sum[20 * A + 16])) + rowSums(select(., columns_to_sum[21 * A + 16])) +
         rowSums(select(., columns_to_sum[22 * A + 16]))),
    P_Total = rowSums(select(., all_of(columns_to_sum[1 : (22 * A)]))),
    P_0to4yrs = (rowSums(select(., columns_to_sum[1])) + rowSums(select(., columns_to_sum[2])) + # Combine 0 - 4 years (0 - 6 months + 7 months - 4 years)
                   rowSums(select(., columns_to_sum[A + 1])) + rowSums(select(., columns_to_sum[A + 2])) +
                   rowSums(select(., columns_to_sum[2 * A + 1])) + rowSums(select(., columns_to_sum[2 * A + 2])) +
                   rowSums(select(., columns_to_sum[3 * A + 1])) + rowSums(select(., columns_to_sum[3 * A + 2])) +
                   rowSums(select(., columns_to_sum[4 * A + 1])) + rowSums(select(., columns_to_sum[4 * A + 2])) +
                   rowSums(select(., columns_to_sum[5 * A + 1])) + rowSums(select(., columns_to_sum[5 * A + 2])) +
                   rowSums(select(., columns_to_sum[6 * A + 1])) + rowSums(select(., columns_to_sum[6 * A + 2])) +
                   rowSums(select(., columns_to_sum[7 * A + 1])) + rowSums(select(., columns_to_sum[7 * A + 2])) +
                   rowSums(select(., columns_to_sum[8 * A + 1])) + rowSums(select(., columns_to_sum[8 * A + 2])) +
                   rowSums(select(., columns_to_sum[9 * A + 1])) + rowSums(select(., columns_to_sum[9 * A + 2])) +
                   rowSums(select(., columns_to_sum[10 * A + 1])) + rowSums(select(., columns_to_sum[10 * A + 2])) +
                   rowSums(select(., columns_to_sum[11 * A + 1])) + rowSums(select(., columns_to_sum[11 * A + 2])) +
                   rowSums(select(., columns_to_sum[12 * A + 1])) + rowSums(select(., columns_to_sum[12 * A + 2])) +
                   rowSums(select(., columns_to_sum[13 * A + 1])) + rowSums(select(., columns_to_sum[13 * A + 2])) +
                   rowSums(select(., columns_to_sum[14 * A + 1])) + rowSums(select(., columns_to_sum[14 * A + 2])) +
                   rowSums(select(., columns_to_sum[15 * A + 1])) + rowSums(select(., columns_to_sum[15 * A + 2])) +
                   rowSums(select(., columns_to_sum[16 * A + 1])) + rowSums(select(., columns_to_sum[16 * A + 2])) +
                   rowSums(select(., columns_to_sum[17 * A + 1])) + rowSums(select(., columns_to_sum[17 * A + 2])) +
                   rowSums(select(., columns_to_sum[18 * A + 1])) + rowSums(select(., columns_to_sum[18 * A + 2])) +
                   rowSums(select(., columns_to_sum[19 * A + 1])) + rowSums(select(., columns_to_sum[19 * A + 2])) +
                   rowSums(select(., columns_to_sum[20 * A + 1])) + rowSums(select(., columns_to_sum[20 * A + 2])) +
                   rowSums(select(., columns_to_sum[21 * A + 1])) + rowSums(select(., columns_to_sum[21 * A + 2])) +
                   rowSums(select(., columns_to_sum[22 * A + 1])) + rowSums(select(., columns_to_sum[22 * A + 2]))),
    P_5to9yrs = (rowSums(select(., columns_to_sum[3])) + rowSums(select(., columns_to_sum[A + 3])) +
                   rowSums(select(., columns_to_sum[2 * A + 3])) + rowSums(select(., columns_to_sum[3 * A + 3])) +
                   rowSums(select(., columns_to_sum[4 * A + 3])) + rowSums(select(., columns_to_sum[5 * A + 3])) +
                   rowSums(select(., columns_to_sum[6 * A + 3])) + rowSums(select(., columns_to_sum[7 * A + 3])) +
                   rowSums(select(., columns_to_sum[8 * A + 3])) + rowSums(select(., columns_to_sum[9 * A + 3])) +
                   rowSums(select(., columns_to_sum[10 * A + 3])) + rowSums(select(., columns_to_sum[11 * A + 3])) +
                   rowSums(select(., columns_to_sum[12 * A + 3])) + rowSums(select(., columns_to_sum[13 * A + 3])) +
                   rowSums(select(., columns_to_sum[14 * A + 3])) + rowSums(select(., columns_to_sum[15 * A + 3])) +
                   rowSums(select(., columns_to_sum[16 * A + 3])) + rowSums(select(., columns_to_sum[17 * A + 3])) +
                   rowSums(select(., columns_to_sum[18 * A + 3])) + rowSums(select(., columns_to_sum[19 * A + 3])) +
                   rowSums(select(., columns_to_sum[20 * A + 3])) + rowSums(select(., columns_to_sum[21 * A + 3])) +
                   rowSums(select(., columns_to_sum[22 * A + 3]))),
    P_10to14yrs = (rowSums(select(., columns_to_sum[4])) + rowSums(select(., columns_to_sum[A + 4])) +
                     rowSums(select(., columns_to_sum[2 * A + 4])) + rowSums(select(., columns_to_sum[3 * A + 4])) +
                     rowSums(select(., columns_to_sum[4 * A + 4])) + rowSums(select(., columns_to_sum[5 * A + 4])) +
                     rowSums(select(., columns_to_sum[6 * A + 4])) + rowSums(select(., columns_to_sum[7 * A + 4])) +
                     rowSums(select(., columns_to_sum[8 * A + 4])) + rowSums(select(., columns_to_sum[9 * A + 4])) +
                     rowSums(select(., columns_to_sum[10 * A + 4])) + rowSums(select(., columns_to_sum[11 * A + 4])) +
                     rowSums(select(., columns_to_sum[12 * A + 4])) + rowSums(select(., columns_to_sum[13 * A + 4])) +
                     rowSums(select(., columns_to_sum[14 * A + 4])) + rowSums(select(., columns_to_sum[15 * A + 4])) +
                     rowSums(select(., columns_to_sum[16 * A + 4])) + rowSums(select(., columns_to_sum[17 * A + 4])) +
                     rowSums(select(., columns_to_sum[18 * A + 4])) + rowSums(select(., columns_to_sum[19 * A + 4])) +
                     rowSums(select(., columns_to_sum[20 * A + 4])) + rowSums(select(., columns_to_sum[21 * A + 4])) +
                     rowSums(select(., columns_to_sum[22 * A + 4]))),
    P_15to19yrs = (rowSums(select(., columns_to_sum[5])) + rowSums(select(., columns_to_sum[A + 5])) +
                     rowSums(select(., columns_to_sum[2 * A + 5])) + rowSums(select(., columns_to_sum[3 * A + 5])) +
                     rowSums(select(., columns_to_sum[4 * A + 5])) + rowSums(select(., columns_to_sum[5 * A + 5])) +
                     rowSums(select(., columns_to_sum[6 * A + 5])) + rowSums(select(., columns_to_sum[7 * A + 5])) +
                     rowSums(select(., columns_to_sum[8 * A + 5])) + rowSums(select(., columns_to_sum[9 * A + 5])) +
                     rowSums(select(., columns_to_sum[10 * A + 5])) + rowSums(select(., columns_to_sum[11 * A + 5])) +
                     rowSums(select(., columns_to_sum[12 * A + 5])) + rowSums(select(., columns_to_sum[13 * A + 5])) +
                     rowSums(select(., columns_to_sum[14 * A + 5])) + rowSums(select(., columns_to_sum[15 * A + 5])) +
                     rowSums(select(., columns_to_sum[16 * A + 5])) + rowSums(select(., columns_to_sum[17 * A + 5])) +
                     rowSums(select(., columns_to_sum[18 * A + 5])) + rowSums(select(., columns_to_sum[19 * A + 5])) +
                     rowSums(select(., columns_to_sum[20 * A + 5])) + rowSums(select(., columns_to_sum[21 * A + 5])) +
                     rowSums(select(., columns_to_sum[22 * A + 5]))),
    P_20to24yrs = (rowSums(select(., columns_to_sum[6])) + rowSums(select(., columns_to_sum[A + 6])) +
                     rowSums(select(., columns_to_sum[2 * A + 6])) + rowSums(select(., columns_to_sum[3 * A + 6])) +
                     rowSums(select(., columns_to_sum[4 * A + 6])) + rowSums(select(., columns_to_sum[5 * A + 6])) +
                     rowSums(select(., columns_to_sum[6 * A + 6])) + rowSums(select(., columns_to_sum[7 * A + 6])) +
                     rowSums(select(., columns_to_sum[8 * A + 6])) + rowSums(select(., columns_to_sum[9 * A + 6])) +
                     rowSums(select(., columns_to_sum[10 * A + 6])) + rowSums(select(., columns_to_sum[11 * A + 6])) +
                     rowSums(select(., columns_to_sum[12 * A + 6])) + rowSums(select(., columns_to_sum[13 * A + 6])) +
                     rowSums(select(., columns_to_sum[14 * A + 6])) + rowSums(select(., columns_to_sum[15 * A + 6])) +
                     rowSums(select(., columns_to_sum[16 * A + 6])) + rowSums(select(., columns_to_sum[17 * A + 6])) +
                     rowSums(select(., columns_to_sum[18 * A + 6])) + rowSums(select(., columns_to_sum[19 * A + 6])) +
                     rowSums(select(., columns_to_sum[20 * A + 6])) + rowSums(select(., columns_to_sum[21 * A + 6])) +
                     rowSums(select(., columns_to_sum[22 * A + 6]))),
    P_25to29yrs = (rowSums(select(., columns_to_sum[7])) + rowSums(select(., columns_to_sum[A + 7])) +
                     rowSums(select(., columns_to_sum[2 * A + 7])) + rowSums(select(., columns_to_sum[3 * A + 7])) +
                     rowSums(select(., columns_to_sum[4 * A + 7])) + rowSums(select(., columns_to_sum[5 * A + 7])) +
                     rowSums(select(., columns_to_sum[6 * A + 7])) + rowSums(select(., columns_to_sum[7 * A + 7])) +
                     rowSums(select(., columns_to_sum[8 * A + 7])) + rowSums(select(., columns_to_sum[9 * A + 7])) +
                     rowSums(select(., columns_to_sum[10 * A + 7])) + rowSums(select(., columns_to_sum[11 * A + 7])) +
                     rowSums(select(., columns_to_sum[12 * A + 7])) + rowSums(select(., columns_to_sum[13 * A + 7])) +
                     rowSums(select(., columns_to_sum[14 * A + 7])) + rowSums(select(., columns_to_sum[15 * A + 7])) +
                     rowSums(select(., columns_to_sum[16 * A + 7])) + rowSums(select(., columns_to_sum[17 * A + 7])) +
                     rowSums(select(., columns_to_sum[18 * A + 7])) + rowSums(select(., columns_to_sum[19 * A + 7])) +
                     rowSums(select(., columns_to_sum[20 * A + 7])) + rowSums(select(., columns_to_sum[21 * A + 7])) +
                     rowSums(select(., columns_to_sum[22 * A + 7]))),
    P_30to34yrs = (rowSums(select(., columns_to_sum[8])) + rowSums(select(., columns_to_sum[A + 8])) +
                     rowSums(select(., columns_to_sum[2 * A + 8])) + rowSums(select(., columns_to_sum[3 * A + 8])) +
                     rowSums(select(., columns_to_sum[4 * A + 8])) + rowSums(select(., columns_to_sum[5 * A + 8])) +
                     rowSums(select(., columns_to_sum[6 * A + 8])) + rowSums(select(., columns_to_sum[7 * A + 8])) +
                     rowSums(select(., columns_to_sum[8 * A + 8])) + rowSums(select(., columns_to_sum[9 * A + 8])) +
                     rowSums(select(., columns_to_sum[10 * A + 8])) + rowSums(select(., columns_to_sum[11 * A + 8])) +
                     rowSums(select(., columns_to_sum[12 * A + 8])) + rowSums(select(., columns_to_sum[13 * A + 8])) +
                     rowSums(select(., columns_to_sum[14 * A + 8])) + rowSums(select(., columns_to_sum[15 * A + 8])) +
                     rowSums(select(., columns_to_sum[16 * A + 8])) + rowSums(select(., columns_to_sum[17 * A + 8])) +
                     rowSums(select(., columns_to_sum[18 * A + 8])) + rowSums(select(., columns_to_sum[19 * A + 8])) +
                     rowSums(select(., columns_to_sum[20 * A + 8])) + rowSums(select(., columns_to_sum[21 * A + 8])) +
                     rowSums(select(., columns_to_sum[22 * A + 8]))),
    P_35to39yrs = (rowSums(select(., columns_to_sum[9])) + rowSums(select(., columns_to_sum[A + 9])) +
                     rowSums(select(., columns_to_sum[2 * A + 9])) + rowSums(select(., columns_to_sum[3 * A + 9])) +
                     rowSums(select(., columns_to_sum[4 * A + 9])) + rowSums(select(., columns_to_sum[5 * A + 9])) +
                     rowSums(select(., columns_to_sum[6 * A + 9])) + rowSums(select(., columns_to_sum[7 * A + 9])) +
                     rowSums(select(., columns_to_sum[8 * A + 9])) + rowSums(select(., columns_to_sum[9 * A + 9])) +
                     rowSums(select(., columns_to_sum[10 * A + 9])) + rowSums(select(., columns_to_sum[11 * A + 9])) +
                     rowSums(select(., columns_to_sum[12 * A + 9])) + rowSums(select(., columns_to_sum[13 * A + 9])) +
                     rowSums(select(., columns_to_sum[14 * A + 9])) + rowSums(select(., columns_to_sum[15 * A + 9])) +
                     rowSums(select(., columns_to_sum[16 * A + 9])) + rowSums(select(., columns_to_sum[17 * A + 9])) +
                     rowSums(select(., columns_to_sum[18 * A + 9])) + rowSums(select(., columns_to_sum[19 * A + 9])) +
                     rowSums(select(., columns_to_sum[20 * A + 9])) + rowSums(select(., columns_to_sum[21 * A + 9])) +
                     rowSums(select(., columns_to_sum[22 * A + 9]))),
    P_40to44yrs = (rowSums(select(., columns_to_sum[10])) + rowSums(select(., columns_to_sum[A + 10])) +
                     rowSums(select(., columns_to_sum[2 * A + 10])) + rowSums(select(., columns_to_sum[3 * A + 10])) +
                     rowSums(select(., columns_to_sum[4 * A + 10])) + rowSums(select(., columns_to_sum[5 * A + 10])) +
                     rowSums(select(., columns_to_sum[6 * A + 10])) + rowSums(select(., columns_to_sum[7 * A + 10])) +
                     rowSums(select(., columns_to_sum[8 * A + 10])) + rowSums(select(., columns_to_sum[9 * A + 10])) +
                     rowSums(select(., columns_to_sum[10 * A + 10])) + rowSums(select(., columns_to_sum[11 * A + 10])) +
                     rowSums(select(., columns_to_sum[12 * A + 10])) + rowSums(select(., columns_to_sum[13 * A + 10])) +
                     rowSums(select(., columns_to_sum[14 * A + 10])) + rowSums(select(., columns_to_sum[15 * A + 10])) +
                     rowSums(select(., columns_to_sum[16 * A + 10])) + rowSums(select(., columns_to_sum[17 * A + 10])) +
                     rowSums(select(., columns_to_sum[18 * A + 10])) + rowSums(select(., columns_to_sum[19 * A + 10])) +
                     rowSums(select(., columns_to_sum[20 * A + 10])) + rowSums(select(., columns_to_sum[21 * A + 10])) +
                     rowSums(select(., columns_to_sum[22 * A + 10]))),
    P_45to49yrs = (rowSums(select(., columns_to_sum[11])) + rowSums(select(., columns_to_sum[A + 11])) +
                     rowSums(select(., columns_to_sum[2 * A + 11])) + rowSums(select(., columns_to_sum[3 * A + 11])) +
                     rowSums(select(., columns_to_sum[4 * A + 11])) + rowSums(select(., columns_to_sum[5 * A + 11])) +
                     rowSums(select(., columns_to_sum[6 * A + 11])) + rowSums(select(., columns_to_sum[7 * A + 11])) +
                     rowSums(select(., columns_to_sum[8 * A + 11])) + rowSums(select(., columns_to_sum[9 * A + 11])) +
                     rowSums(select(., columns_to_sum[10 * A + 11])) + rowSums(select(., columns_to_sum[11 * A + 11])) +
                     rowSums(select(., columns_to_sum[12 * A + 11])) + rowSums(select(., columns_to_sum[13 * A + 11])) +
                     rowSums(select(., columns_to_sum[14 * A + 11])) + rowSums(select(., columns_to_sum[15 * A + 11])) +
                     rowSums(select(., columns_to_sum[16 * A + 11])) + rowSums(select(., columns_to_sum[17 * A + 11])) +
                     rowSums(select(., columns_to_sum[18 * A + 11])) + rowSums(select(., columns_to_sum[19 * A + 11])) +
                     rowSums(select(., columns_to_sum[20 * A + 11])) + rowSums(select(., columns_to_sum[21 * A + 11])) +
                     rowSums(select(., columns_to_sum[22 * A + 11]))),
    P_50to54yrs = (rowSums(select(., columns_to_sum[12])) + rowSums(select(., columns_to_sum[A + 12])) +
                     rowSums(select(., columns_to_sum[2 * A + 12])) + rowSums(select(., columns_to_sum[3 * A + 12])) +
                     rowSums(select(., columns_to_sum[4 * A + 12])) + rowSums(select(., columns_to_sum[5 * A + 12])) +
                     rowSums(select(., columns_to_sum[6 * A + 12])) + rowSums(select(., columns_to_sum[7 * A + 12])) +
                     rowSums(select(., columns_to_sum[8 * A + 12])) + rowSums(select(., columns_to_sum[9 * A + 12])) +
                     rowSums(select(., columns_to_sum[10 * A + 12])) + rowSums(select(., columns_to_sum[11 * A + 12])) +
                     rowSums(select(., columns_to_sum[12 * A + 12])) + rowSums(select(., columns_to_sum[13 * A + 12])) +
                     rowSums(select(., columns_to_sum[14 * A + 12])) + rowSums(select(., columns_to_sum[15 * A + 12])) +
                     rowSums(select(., columns_to_sum[16 * A + 12])) + rowSums(select(., columns_to_sum[17 * A + 12])) +
                     rowSums(select(., columns_to_sum[18 * A + 12])) + rowSums(select(., columns_to_sum[19 * A + 12])) +
                     rowSums(select(., columns_to_sum[20 * A + 12])) + rowSums(select(., columns_to_sum[21 * A + 12])) +
                     rowSums(select(., columns_to_sum[22 * A + 12]))),
    P_55to59yrs = (rowSums(select(., columns_to_sum[13])) + rowSums(select(., columns_to_sum[A + 13])) +
                     rowSums(select(., columns_to_sum[2 * A + 13])) + rowSums(select(., columns_to_sum[3 * A + 13])) +
                     rowSums(select(., columns_to_sum[4 * A + 13])) + rowSums(select(., columns_to_sum[5 * A + 13])) +
                     rowSums(select(., columns_to_sum[6 * A + 13])) + rowSums(select(., columns_to_sum[7 * A + 13])) +
                     rowSums(select(., columns_to_sum[8 * A + 13])) + rowSums(select(., columns_to_sum[9 * A + 13])) +
                     rowSums(select(., columns_to_sum[10 * A + 13])) + rowSums(select(., columns_to_sum[11 * A + 13])) +
                     rowSums(select(., columns_to_sum[12 * A + 13])) + rowSums(select(., columns_to_sum[13 * A + 13])) +
                     rowSums(select(., columns_to_sum[14 * A + 13])) + rowSums(select(., columns_to_sum[15 * A + 13])) +
                     rowSums(select(., columns_to_sum[16 * A + 13])) + rowSums(select(., columns_to_sum[17 * A + 13])) +
                     rowSums(select(., columns_to_sum[18 * A + 13])) + rowSums(select(., columns_to_sum[19 * A + 13])) +
                     rowSums(select(., columns_to_sum[20 * A + 13])) + rowSums(select(., columns_to_sum[21 * A + 13])) +
                     rowSums(select(., columns_to_sum[22 * A + 13]))),
    P_60to64yrs = (rowSums(select(., columns_to_sum[14])) + rowSums(select(., columns_to_sum[A + 14])) +
                     rowSums(select(., columns_to_sum[2 * A + 14])) + rowSums(select(., columns_to_sum[3 * A + 14])) +
                     rowSums(select(., columns_to_sum[4 * A + 14])) + rowSums(select(., columns_to_sum[5 * A + 14])) +
                     rowSums(select(., columns_to_sum[6 * A + 14])) + rowSums(select(., columns_to_sum[7 * A + 14])) +
                     rowSums(select(., columns_to_sum[8 * A + 14])) + rowSums(select(., columns_to_sum[9 * A + 14])) +
                     rowSums(select(., columns_to_sum[10 * A + 14])) + rowSums(select(., columns_to_sum[11 * A + 14])) +
                     rowSums(select(., columns_to_sum[12 * A + 14])) + rowSums(select(., columns_to_sum[13 * A + 14])) +
                     rowSums(select(., columns_to_sum[14 * A + 14])) + rowSums(select(., columns_to_sum[15 * A + 14])) +
                     rowSums(select(., columns_to_sum[16 *A + 14])) + rowSums(select(., columns_to_sum[17 * A + 14])) +
                     rowSums(select(., columns_to_sum[18 * A + 14])) + rowSums(select(., columns_to_sum[19 * A + 14])) +
                     rowSums(select(., columns_to_sum[20 * A + 14])) + rowSums(select(., columns_to_sum[21 * A + 14])) +
                     rowSums(select(., columns_to_sum[22 * A + 14]))),
    P_65to69yrs = (rowSums(select(., columns_to_sum[15])) + rowSums(select(., columns_to_sum[A + 15])) +
                     rowSums(select(., columns_to_sum[2 * A + 15])) + rowSums(select(., columns_to_sum[3 * A + 15])) +
                     rowSums(select(., columns_to_sum[4 * A + 15])) + rowSums(select(., columns_to_sum[5 * A + 15])) +
                     rowSums(select(., columns_to_sum[6 * A + 15])) + rowSums(select(., columns_to_sum[7 * A + 15])) +
                     rowSums(select(., columns_to_sum[8 * A + 15])) + rowSums(select(., columns_to_sum[9 * A + 15])) +
                     rowSums(select(., columns_to_sum[10 * A + 15])) + rowSums(select(., columns_to_sum[11 * A + 15])) +
                     rowSums(select(., columns_to_sum[12 * A + 15])) + rowSums(select(., columns_to_sum[13 * A + 15])) +
                     rowSums(select(., columns_to_sum[14 * A + 15])) + rowSums(select(., columns_to_sum[15 * A + 15])) +
                     rowSums(select(., columns_to_sum[16 * A + 15])) + rowSums(select(., columns_to_sum[17 * A + 15])) +
                     rowSums(select(., columns_to_sum[18 * A + 15])) + rowSums(select(., columns_to_sum[19 * A + 15])) +
                     rowSums(select(., columns_to_sum[20 * A + 15])) + rowSums(select(., columns_to_sum[21 * A + 15])) +
                     rowSums(select(., columns_to_sum[22 * A + 15]))),
    P_70plusyrs = (rowSums(select(., columns_to_sum[16])) + rowSums(select(., columns_to_sum[A + 16])) +
                     rowSums(select(., columns_to_sum[2 * A + 16])) + rowSums(select(., columns_to_sum[3 * A + 16])) +
                     rowSums(select(., columns_to_sum[4 * A + 16])) + rowSums(select(., columns_to_sum[5 * A + 16])) +
                     rowSums(select(., columns_to_sum[6 * A + 16])) + rowSums(select(., columns_to_sum[7 * A + 16])) +
                     rowSums(select(., columns_to_sum[8 * A + 16])) + rowSums(select(., columns_to_sum[9 * A + 16])) +
                     rowSums(select(., columns_to_sum[10 * A + 16])) + rowSums(select(., columns_to_sum[11 * A + 16])) +
                     rowSums(select(., columns_to_sum[12 * A + 16])) + rowSums(select(., columns_to_sum[13 * A + 16])) +
                     rowSums(select(., columns_to_sum[14 * A + 16])) + rowSums(select(., columns_to_sum[15 * A + 16])) +
                     rowSums(select(., columns_to_sum[16 * A + 16])) + rowSums(select(., columns_to_sum[17 * A + 16])) +
                     rowSums(select(., columns_to_sum[18 * A + 16])) + rowSums(select(., columns_to_sum[19 * A + 16])) +
                     rowSums(select(., columns_to_sum[20 * A + 16])) + rowSums(select(., columns_to_sum[21 * A + 16])) +
                     rowSums(select(., columns_to_sum[22 * A + 16]))),
  )

#View(df_laopdr_HBV_output_bc1)

######################################################################
# FITTING - PLOTTING POPULATION DYNAMICS
######################################################################

# Creating vectors of year data points to plot population demographics against
popyearsvector_data <- c(1995, 2005, 2015) # The years between 1995 - 2015 are official Lao Statisitcs Bureau Census Data collected once every 10 years
popyearsvector_projections <- c(2020, 2025, 2030, 2035, 2040, 2045) # The years beyond 2020 are official Lao Statistics Bureau Population Projections performed in July 2018

#Population Data
popvector_0to4yrs_data  <- c(sum(popstruc[1:2,2]), sum(popstruc[1:2,5]), sum(popstruc[1:2,8]))
popvector_5to9yrs_data  <- c(popstruc[3,2], popstruc[3,5], popstruc[3,8])
popvector_10to14yrs_data  <- c(popstruc[4,2], popstruc[4,5], popstruc[4,8])
popvector_15to19yrs_data  <- c(popstruc[5,2], popstruc[5,5], popstruc[5,8])
popvector_20to24yrs_data  <- c(popstruc[6,2], popstruc[6,5], popstruc[6,8])
popvector_25to29yrs_data  <- c(popstruc[7,2], popstruc[7,5], popstruc[7,8])
popvector_30to34yrs_data  <- c(popstruc[8,2], popstruc[8,5], popstruc[8,8])
popvector_35to39yrs_data  <- c(popstruc[9,2], popstruc[9,5], popstruc[9,8])
popvector_40to44yrs_data  <- c(popstruc[10,2], popstruc[10,5], popstruc[10,8])
popvector_45to49yrs_data  <- c(popstruc[11,2], popstruc[11,5], popstruc[11,8])
popvector_50to54yrs_data  <- c(popstruc[12,2], popstruc[12,5], popstruc[12,8])
popvector_55to59yrs_data  <- c(popstruc[13,2], popstruc[13,5], popstruc[13,8])
popvector_60to64yrs_data  <- c(popstruc[14,2], popstruc[14,5], popstruc[14,8])
popvector_65to69yrs_data  <- c(popstruc[15,2], popstruc[15,5], popstruc[15,8])
popvector_70plusyrs_data  <- c(popstruc[16,2], popstruc[16,5], popstruc[16,8])

# Low Population Projections
popvector_0to4yrs_low <- c(sum(popstruc[1:2,11]), sum(popstruc[1:2,20]), sum(popstruc[1:2,29]), sum(popstruc[1:2,38]), sum(popstruc[1:2,47]), sum(popstruc[1:2,56]))
popvector_5to9yrs_low <- c(popstruc[3,11], popstruc[3,20], popstruc[3,29], popstruc[3,38], popstruc[3,47], popstruc[3,56])
popvector_10to14yrs_low <- c(popstruc[4,11], popstruc[4,20], popstruc[4,29], popstruc[4,38], popstruc[4,47], popstruc[4,56])
popvector_15to19yrs_low <- c(popstruc[5,11], popstruc[5,20], popstruc[5,29], popstruc[5,38], popstruc[5,47], popstruc[5,56])
popvector_20to24yrs_low <- c(popstruc[6,11], popstruc[6,20], popstruc[6,29], popstruc[6,38], popstruc[6,47], popstruc[6,56])
popvector_25to29yrs_low <- c(popstruc[7,11], popstruc[7,20], popstruc[7,29], popstruc[7,38], popstruc[7,47], popstruc[7,56])
popvector_30to34yrs_low <- c(popstruc[8,11], popstruc[8,20], popstruc[8,29], popstruc[8,38], popstruc[8,47], popstruc[8,56])
popvector_35to39yrs_low <- c(popstruc[9,11], popstruc[9,20], popstruc[9,29], popstruc[9,38], popstruc[9,47], popstruc[9,56])
popvector_40to44yrs_low <- c(popstruc[10,11], popstruc[10,20], popstruc[10,29], popstruc[10,38], popstruc[10,47], popstruc[10,56])
popvector_45to49yrs_low <- c(popstruc[11,11], popstruc[11,20], popstruc[11,29], popstruc[11,38], popstruc[11,47], popstruc[11,56])
popvector_50to54yrs_low <- c(popstruc[12,11], popstruc[12,20], popstruc[12,29], popstruc[12,38], popstruc[12,47], popstruc[12,56])
popvector_55to59yrs_low <- c(popstruc[13,11], popstruc[13,20], popstruc[13,29], popstruc[13,38], popstruc[13,47], popstruc[13,56])
popvector_60to64yrs_low <- c(popstruc[14,11], popstruc[14,20], popstruc[14,29], popstruc[14,38], popstruc[14,47], popstruc[14,56])
popvector_65to69yrs_low <- c(popstruc[15,11], popstruc[15,20], popstruc[15,29], popstruc[15,38], popstruc[15,47], popstruc[15,56])
popvector_70plusyrs_low <- c(popstruc[16,11], popstruc[16,20], popstruc[16,29], popstruc[16,38], popstruc[16,47], popstruc[16,56])

# Medium Population Projections
popvector_0to4yrs_medium <- c(sum(popstruc[1:2,14]), sum(popstruc[1:2,23]), sum(popstruc[1:2,32]), sum(popstruc[1:2,41]), sum(popstruc[1:2,50]), sum(popstruc[1:2,59]))
popvector_5to9yrs_medium <- c(popstruc[3,14], popstruc[3,23], popstruc[3,32], popstruc[3,41], popstruc[3,50], popstruc[3,59])
popvector_10to14yrs_medium <- c(popstruc[4,14], popstruc[4,23], popstruc[4,32], popstruc[4,41], popstruc[4,50], popstruc[4,59])
popvector_15to19yrs_medium <- c(popstruc[5,14], popstruc[5,23], popstruc[5,32], popstruc[5,41], popstruc[5,50], popstruc[5,59])
popvector_20to24yrs_medium <- c(popstruc[6,14], popstruc[6,23], popstruc[6,32], popstruc[6,41], popstruc[6,50], popstruc[6,59])
popvector_25to29yrs_medium <- c(popstruc[7,14], popstruc[7,23], popstruc[7,32], popstruc[7,41], popstruc[7,50], popstruc[7,59])
popvector_30to34yrs_medium <- c(popstruc[8,14], popstruc[8,23], popstruc[8,32], popstruc[8,41], popstruc[8,50], popstruc[8,59])
popvector_35to39yrs_medium <- c(popstruc[9,14], popstruc[9,23], popstruc[9,32], popstruc[9,41], popstruc[9,50], popstruc[9,59])
popvector_40to44yrs_medium <- c(popstruc[10,14], popstruc[10,23], popstruc[10,32], popstruc[10,41], popstruc[10,50], popstruc[10,59])
popvector_45to49yrs_medium <- c(popstruc[11,14], popstruc[11,23], popstruc[11,32], popstruc[11,41], popstruc[11,50], popstruc[11,59])
popvector_50to54yrs_medium <- c(popstruc[12,14], popstruc[12,23], popstruc[12,32], popstruc[12,41], popstruc[12,50], popstruc[12,59])
popvector_55to59yrs_medium <- c(popstruc[13,14], popstruc[13,23], popstruc[13,32], popstruc[13,41], popstruc[13,50], popstruc[13,59])
popvector_60to64yrs_medium <- c(popstruc[14,14], popstruc[14,23], popstruc[14,32], popstruc[14,41], popstruc[14,50], popstruc[14,59])
popvector_65to69yrs_medium <- c(popstruc[15,14], popstruc[15,23], popstruc[15,32], popstruc[15,41], popstruc[15,50], popstruc[15,59])
popvector_70plusyrs_medium <- c(popstruc[16,14], popstruc[16,23], popstruc[16,32], popstruc[16,41], popstruc[16,50], popstruc[16,59])

# High Population Projections
popvector_0to4yrs_high <- c(sum(popstruc[1:2,17]), sum(popstruc[1:2,26]), sum(popstruc[1:2,35]), sum(popstruc[1:2,44]), sum(popstruc[1:2,53]), sum(popstruc[1:2,62]))
popvector_5to9yrs_high <- c(popstruc[3,17], popstruc[3,26], popstruc[3,35], popstruc[3,44], popstruc[3,53], popstruc[3,62])
popvector_10to14yrs_high <- c(popstruc[4,17], popstruc[4,26], popstruc[4,35], popstruc[4,44], popstruc[4,53], popstruc[4,62])
popvector_15to19yrs_high <- c(popstruc[5,17], popstruc[5,26], popstruc[5,35], popstruc[5,44], popstruc[5,53], popstruc[5,62])
popvector_20to24yrs_high <- c(popstruc[6,17], popstruc[6,26], popstruc[6,35], popstruc[6,44], popstruc[6,53], popstruc[6,62])
popvector_25to29yrs_high <- c(popstruc[7,17], popstruc[7,26], popstruc[7,35], popstruc[7,44], popstruc[7,53], popstruc[7,62])
popvector_30to34yrs_high <- c(popstruc[8,17], popstruc[8,26], popstruc[8,35], popstruc[8,44], popstruc[8,53], popstruc[8,62])
popvector_35to39yrs_high <- c(popstruc[9,17], popstruc[9,26], popstruc[9,35], popstruc[9,44], popstruc[9,53], popstruc[9,62])
popvector_40to44yrs_high <- c(popstruc[10,17], popstruc[10,26], popstruc[10,35], popstruc[10,44], popstruc[10,53], popstruc[10,62])
popvector_45to49yrs_high <- c(popstruc[11,17], popstruc[11,26], popstruc[11,35], popstruc[11,44], popstruc[11,53], popstruc[11,62])
popvector_50to54yrs_high <- c(popstruc[12,17], popstruc[12,26], popstruc[12,35], popstruc[12,44], popstruc[12,53], popstruc[12,62])
popvector_55to59yrs_high <- c(popstruc[13,17], popstruc[13,26], popstruc[13,35], popstruc[13,44], popstruc[13,53], popstruc[13,62])
popvector_60to64yrs_high <- c(popstruc[14,17], popstruc[14,26], popstruc[14,35], popstruc[14,44], popstruc[14,53], popstruc[14,62])
popvector_65to69yrs_high <- c(popstruc[15,17], popstruc[15,26], popstruc[15,35], popstruc[15,44], popstruc[15,53], popstruc[15,62])
popvector_70plusyrs_high <- c(popstruc[16,17], popstruc[16,26], popstruc[16,35], popstruc[16,44], popstruc[16,53], popstruc[16,62])

#Storing population data in a data frame 
pop_data <- data.frame(
  Year = rep(popyearsvector_data, 15),
  Population = c(popvector_0to4yrs_data, popvector_5to9yrs_data, popvector_10to14yrs_data, popvector_15to19yrs_data,
                 popvector_20to24yrs_data, popvector_25to29yrs_data, popvector_30to34yrs_data, popvector_35to39yrs_data,
                 popvector_40to44yrs_data, popvector_45to49yrs_data, popvector_50to54yrs_data, popvector_55to59yrs_data,
                 popvector_60to64yrs_data, popvector_65to69yrs_data, popvector_70plusyrs_data),
  Age_Group = rep(c("0to4yrs", "5to9yrs", "10to14yrs", "15to19yrs", "20to24yrs", "25to29yrs", "30to34yrs", "35to39yrs",
                    "40to44yrs", "45to49yrs", "50to54yrs", "55to59yrs", "60to64yrs", "65to69yrs", "70plusyrs"), each = length(popyearsvector_data)),
  Projection = 'Data'
)

#Storing low population projections in a data frame 
pop_projections_low <- data.frame(
  Year = rep(popyearsvector_projections, 15),
  Population = c(popvector_0to4yrs_low, popvector_5to9yrs_low, popvector_10to14yrs_low, popvector_15to19yrs_low,
                 popvector_20to24yrs_low, popvector_25to29yrs_low, popvector_30to34yrs_low, popvector_35to39yrs_low,
                 popvector_40to44yrs_low, popvector_45to49yrs_low, popvector_50to54yrs_low, popvector_55to59yrs_low,
                 popvector_60to64yrs_low, popvector_65to69yrs_low, popvector_70plusyrs_low),
  Age_Group = rep(c("0to4yrs", "5to9yrs", "10to14yrs", "15to19yrs", "20to24yrs", "25to29yrs", "30to34yrs", "35to39yrs",
                    "40to44yrs", "45to49yrs", "50to54yrs", "55to59yrs", "60to64yrs", "65to69yrs", "70plusyrs"), each = length(popyearsvector_projections)),
  Projection = 'Low'
)

#Storing medium population projections in a data frame 
pop_projections_medium <- data.frame(
  Year = rep(popyearsvector_projections, 15),
  Population = c(popvector_0to4yrs_medium, popvector_5to9yrs_medium, popvector_10to14yrs_medium, popvector_15to19yrs_medium,
                 popvector_20to24yrs_medium, popvector_25to29yrs_medium, popvector_30to34yrs_medium, popvector_35to39yrs_medium,
                 popvector_40to44yrs_medium, popvector_45to49yrs_medium, popvector_50to54yrs_medium, popvector_55to59yrs_medium,
                 popvector_60to64yrs_medium, popvector_65to69yrs_medium, popvector_70plusyrs_medium),
  Age_Group = rep(c("0to4yrs", "5to9yrs", "10to14yrs", "15to19yrs", "20to24yrs", "25to29yrs", "30to34yrs", "35to39yrs",
                    "40to44yrs", "45to49yrs", "50to54yrs", "55to59yrs", "60to64yrs", "65to69yrs", "70plusyrs"), each = length(popyearsvector_projections)),
  Projection = 'Medium'
)

#Storing high population projections in a data frame 
pop_projections_high <- data.frame(
  Year = rep(popyearsvector_projections, 15),
  Population = c(popvector_0to4yrs_high, popvector_5to9yrs_high, popvector_10to14yrs_high, popvector_15to19yrs_high,
                 popvector_20to24yrs_high, popvector_25to29yrs_high, popvector_30to34yrs_high, popvector_35to39yrs_high,
                 popvector_40to44yrs_high, popvector_45to49yrs_high, popvector_50to54yrs_high, popvector_55to59yrs_high,
                 popvector_60to64yrs_high, popvector_65to69yrs_high, popvector_70plusyrs_high),
  Age_Group = rep(c("0to4yrs", "5to9yrs", "10to14yrs", "15to19yrs", "20to24yrs", "25to29yrs", "30to34yrs", "35to39yrs",
                    "40to44yrs", "45to49yrs", "50to54yrs", "55to59yrs", "60to64yrs", "65to69yrs", "70plusyrs"), each = length(popyearsvector_projections)),
  Projection = 'High'
)


#####################################################
#PLOT OF POPULATION BY AGE GROUP WITH THE MODEL
#####################################################
annual_summary_popagegroups <- df_laopdr_HBV_output_bc1 %>%
  group_by(years) %>%
  filter(years >= 1985 & years <= 2100) %>%
  summarise(
    AnnualPopulation_0to4yrs = mean(P_0to4yrs, na.rm = TRUE),
    AnnualPopulation_5to9yrs = mean(P_5to9yrs, na.rm = TRUE),
    AnnualPopulation_10to14yrs = mean(P_10to14yrs, na.rm = TRUE),
    AnnualPopulation_15to19yrs = mean(P_15to19yrs, na.rm = TRUE),
    AnnualPopulation_20to24yrs = mean(P_20to24yrs, na.rm = TRUE),
    AnnualPopulation_25to29yrs = mean(P_25to29yrs, na.rm = TRUE),
    AnnualPopulation_30to34yrs = mean(P_30to34yrs, na.rm = TRUE),
    AnnualPopulation_35to39yrs = mean(P_35to39yrs, na.rm = TRUE),
    AnnualPopulation_40to44yrs = mean(P_40to44yrs, na.rm = TRUE),
    AnnualPopulation_45to49yrs = mean(P_45to49yrs, na.rm = TRUE),
    AnnualPopulation_50to54yrs = mean(P_50to54yrs, na.rm = TRUE),
    AnnualPopulation_55to59yrs = mean(P_55to59yrs, na.rm = TRUE),
    AnnualPopulation_60to64yrs = mean(P_60to64yrs, na.rm = TRUE),
    AnnualPopulation_65to69yrs = mean(P_65to69yrs, na.rm = TRUE),
    AnnualPopulation_70plusyrs = mean(P_70plusyrs, na.rm = TRUE),
  )

# Reshape the annual summary data to a long format for plotting
annual_summary_popagegroups_long <- annual_summary_popagegroups %>%
  pivot_longer(cols = starts_with("AnnualPopulation"),
               names_to = "Age_Group",
               values_to = "Population") %>%
  mutate(Age_Group = factor(gsub("AnnualPopulation_", "", Age_Group),
                            levels = c("0to4yrs", "5to9yrs", "10to14yrs", "15to19yrs",
                                       "20to24yrs", "25to29yrs", "30to34yrs", "35to39yrs",
                                       "40to44yrs", "45to49yrs", "50to54yrs", "55to59yrs",
                                       "60to64yrs", "65to69yrs", "70plusyrs")))


#Assigning levels so age groups stay in the same order
pop_data$Age_Group <- factor(pop_data$Age_Group,
                             levels = c("0to4yrs", "5to9yrs", "10to14yrs", "15to19yrs",
                                        "20to24yrs", "25to29yrs", "30to34yrs", "35to39yrs",
                                        "40to44yrs", "45to49yrs", "50to54yrs", "55to59yrs",
                                        "60to64yrs", "65to69yrs", "70plusyrs"))
pop_projections_low$Age_Group <- factor(pop_projections_low$Age_Group,
                                        levels = c("0to4yrs", "5to9yrs", "10to14yrs", "15to19yrs",
                                                   "20to24yrs", "25to29yrs", "30to34yrs", "35to39yrs",
                                                   "40to44yrs", "45to49yrs", "50to54yrs", "55to59yrs",
                                                   "60to64yrs", "65to69yrs", "70plusyrs"))
pop_projections_medium$Age_Group <- factor(pop_projections_medium$Age_Group,
                                           levels = c("0to4yrs", "5to9yrs", "10to14yrs", "15to19yrs",
                                                      "20to24yrs", "25to29yrs", "30to34yrs", "35to39yrs",
                                                      "40to44yrs", "45to49yrs", "50to54yrs", "55to59yrs",
                                                      "60to64yrs", "65to69yrs", "70plusyrs"))
pop_projections_high$Age_Group <- factor(pop_projections_high$Age_Group,
                                         levels = c("0to4yrs", "5to9yrs", "10to14yrs", "15to19yrs",
                                                    "20to24yrs", "25to29yrs", "30to34yrs", "35to39yrs",
                                                    "40to44yrs", "45to49yrs", "50to54yrs", "55to59yrs",
                                                    "60to64yrs", "65to69yrs", "70plusyrs"))

# Plotting
ggplot() +
  geom_line(data = annual_summary_popagegroups_long, aes(x = years, y = Population, color = "Total Population"), size = 1) +
  geom_point(data = pop_data, aes(x = Year, y = Population, color = "Historical Data"), size = 2) +
  geom_point(data = pop_projections_low, aes(x = Year, y = Population, color = "Low Projection"), size = 2, shape = 18) +
  geom_point(data = pop_projections_medium, aes(x = Year, y = Population, color = "Medium Projection"), size = 2, shape = 18) +
  geom_point(data = pop_projections_high, aes(x = Year, y = Population, color = "High Projection"), size = 2, shape = 18) +
  facet_wrap(~ Age_Group, scales = "fixed") +
  labs(title = "Population Projections Over Time by Age Group",
       x = "Year",
       y = "Population") +
  scale_x_continuous(limits = c(1985, NA)) +
  scale_color_manual(values = c("Total Population" = "gray",
                                "Historical Data" = "black",
                                "Low Projection" = "skyblue",
                                "Medium Projection" = "blue",
                                "High Projection" = "navyblue")) +
  theme_minimal() +
  theme(legend.position = "bottom")

####################################################
#PLOT OF TOTAL POPULATION WITHIN THE MODEL
#####################################################

annual_summary_totalpop <- df_laopdr_HBV_output_bc1 %>%
  group_by(years) %>%
  filter(years >= 1985 & years <= 2100) %>%
  summarise(
    AnnualPopulation_Total = mean(P_Total, na.rm = TRUE),
  )

# Create the pop_data_total dataframe
pop_data_total <- data.frame(
  Year = c(1995, 2005, 2015),
  Population = c(sum(popstruc[, 2]), sum(popstruc[, 5]), sum(popstruc[, 8]))
)

# Create the pop_projection_total_low dataframe
pop_projection_total_low <- data.frame(
  Year = c(2020, 2025, 2030, 2035, 2040, 2045),
  Population = c(sum(popstruc[, 11]), sum(popstruc[, 20]), sum(popstruc[, 29]), sum(popstruc[, 38]), sum(popstruc[, 47]), sum(popstruc[, 56]))
)

# Create the pop_projection_total_medium dataframe
pop_projection_total_medium <- data.frame(
  Year = c(2020, 2025, 2030, 2035, 2040, 2045),
  Population = c(sum(popstruc[, 14]), sum(popstruc[, 23]), sum(popstruc[, 32]), sum(popstruc[, 41]), sum(popstruc[, 50]), sum(popstruc[, 59]))
)

# Create the pop_projection_total_high dataframe
pop_projection_total_high <- data.frame(
  Year = c(2020, 2025, 2030, 2035, 2040, 2045),
  Population = c(sum(popstruc[, 17]), sum(popstruc[, 26]), sum(popstruc[, 35]), sum(popstruc[, 44]), sum(popstruc[, 53]), sum(popstruc[, 62]))
)

# Plotting
ggplot() +
  geom_line(data = annual_summary_totalpop, aes(x = years, y = AnnualPopulation_Total, color = "Total Population"), size = 1) +
  geom_point(data = pop_data_total, aes(x = Year, y = Population, color = "Historical Data"), size = 3) +
  geom_point(data = pop_projection_total_low, aes(x = Year, y = Population, color = "Low Projection"), size = 4, shape = 18) +
  geom_point(data = pop_projection_total_medium, aes(x = Year, y = Population, color = "Medium Projection"), size = 4, shape = 18) +
  geom_point(data = pop_projection_total_high, aes(x = Year, y = Population, color = "High Projection"), size = 4, shape = 18) +
  labs(title = "Total Population Projections Over Time",
       x = "Year",
       y = "Population") +
  scale_color_manual(values = c("Historical Data" = "black",
                                "Low Projection" = "skyblue",
                                "Medium Projection" = "blue",
                                "High Projection" = "navyblue")) +
  theme_minimal() +
  theme(legend.position = "bottom")

 
# ######################################################################
# #PLOTTING PREVALENCE
# ######################################################################

#Computing Annual Incidence, Deaths, Clinical Cases, Severe Cases, Hospitalizations, Treatments, Deaths and Vaccines
annual_summary_prevalence_agegroups <- df_laopdr_HBV_output_bc1 %>%
  group_by(years) %>%
  filter(years >= 1985 & years <= 2045) %>%
  summarise(
    AnnualPrevalence_0to4yrs = mean(Prev_0to4yrs, na.rm = TRUE),
    AnnualPrevalence_5to9yrs = mean(Prev_5to9yrs, na.rm = TRUE),
    AnnualPrevalence_10to14yrs = mean(Prev_10to14yrs, na.rm = TRUE),
    AnnualPrevalence_15to19yrs = mean(Prev_15to19yrs, na.rm = TRUE),
    AnnualPrevalence_20to24yrs = mean(Prev_20to24yrs, na.rm = TRUE),
    AnnualPrevalence_25to29yrs = mean(Prev_25to29yrs, na.rm = TRUE),
    AnnualPrevalence_30to34yrs = mean(Prev_30to34yrs, na.rm = TRUE),
    AnnualPrevalence_35to39yrs = mean(Prev_35to39yrs, na.rm = TRUE),
    AnnualPrevalence_40to44yrs = mean(Prev_40to44yrs, na.rm = TRUE),
    AnnualPrevalence_45to49yrs = mean(Prev_45to49yrs, na.rm = TRUE),
    AnnualPrevalence_50to54yrs = mean(Prev_50to54yrs, na.rm = TRUE),
    AnnualPrevalence_55to59yrs = mean(Prev_55to59yrs, na.rm = TRUE),
    AnnualPrevalence_60to64yrs = mean(Prev_60to64yrs, na.rm = TRUE),
    AnnualPrevalence_65to69yrs = mean(Prev_65to69yrs, na.rm = TRUE),
    AnnualPrevalence_70plusyrs = mean(Prev_70plusyrs, na.rm = TRUE),
  )

# Reshape the annual summary data to a long format for plotting
annual_summary_prevalence_agegroups_long <- annual_summary_prevalence_agegroups %>%
  pivot_longer(cols = starts_with("AnnualPrevalence"),
               names_to = "Age_Group",
               values_to = "Prevalence") %>%
  mutate(Age_Group = factor(gsub("AnnualPrevalence_", "", Age_Group),
                            levels = c("0to4yrs", "5to9yrs", "10to14yrs", "15to19yrs",
                                       "20to24yrs", "25to29yrs", "30to34yrs", "35to39yrs",
                                       "40to44yrs", "45to49yrs", "50to54yrs", "55to59yrs",
                                       "60to64yrs", "65to69yrs", "70plusyrs")))

# Creating vectors of year data points to plot prevalence data
prevyearsvector_blooddonors_NBC <- c(2004, 2006) # blood donors at National Blood Centre in Vientiane (Jutavijittum et al., 2007 & Jutavijittum et al., 2014)
prevyearsvector_blooddonors_8provinces <- c(2014) # blood donors at locations in Lao PDR (Nouanthong et al., 2021)
prevyearsvector_clustersampling <- c(2019) # cluster sampling at random villages across Lao PDR (Miyano et al., 2022)
prevyearsvector_crosssectional <- c(2021) # cross-sectional seroprevalence survey (Sitbounlang et al., 2022)

# Prevalence Data from Blood Donors in Vientiane
prevvector_15to19yrs_blooddonors_NBC <- c(HBsAG_prev_data[4,2], HBsAG_prev_data[4,3])
prevvector_20to24yrs_blooddonors_NBC <- c(HBsAG_prev_data[5,2], HBsAG_prev_data[5,3])
prevvector_25to29yrs_blooddonors_NBC <- c(HBsAG_prev_data[6,2], HBsAG_prev_data[6,3])
prevvector_30to34yrs_blooddonors_NBC <- c(HBsAG_prev_data[7,2], HBsAG_prev_data[7,3])
prevvector_35to39yrs_blooddonors_NBC <- c(HBsAG_prev_data[8,2], HBsAG_prev_data[8,3])
prevvector_40to44yrs_blooddonors_NBC <- c(HBsAG_prev_data[9,2], HBsAG_prev_data[9,3])
prevvector_45to49yrs_blooddonors_NBC <- c(HBsAG_prev_data[10,2], HBsAG_prev_data[10,3])
prevvector_50to54yrs_blooddonors_NBC <- c(HBsAG_prev_data[11,2], HBsAG_prev_data[11,3])


# Prevalence Data from Blood Donors Across 8 Provinces in Lao PDR
prevvector_15to19yrs_blooddonors_8provinces <- c(HBsAG_prev_data[4,4])
prevvector_20to24yrs_blooddonors_8provinces <- c(HBsAG_prev_data[5,4])
prevvector_25to29yrs_blooddonors_8provinces <- c(HBsAG_prev_data[6,4])
prevvector_30to34yrs_blooddonors_8provinces <- c(HBsAG_prev_data[7,4])
prevvector_35to39yrs_blooddonors_8provinces <- c(HBsAG_prev_data[8,4])
prevvector_40to44yrs_blooddonors_8provinces <- c(HBsAG_prev_data[9,4])
prevvector_45to49yrs_blooddonors_8provinces <- c(HBsAG_prev_data[10,4])
prevvector_50to54yrs_blooddonors_8provinces <- c(HBsAG_prev_data[11,4])

# Prevalence Data from National Cross-Sectional Survey
prevvector_15to19yrs_clustersampling <- c(HBsAG_prev_data[4,5])
prevvector_20to24yrs_clustersampling <- c(HBsAG_prev_data[5,5])
prevvector_25to29yrs_clustersampling <- c(HBsAG_prev_data[6,5])
prevvector_30to34yrs_clustersampling <- c(HBsAG_prev_data[7,5])
prevvector_35to39yrs_clustersampling <- c(HBsAG_prev_data[8,5])
prevvector_40to44yrs_clustersampling <- c(HBsAG_prev_data[9,5])
prevvector_45to49yrs_clustersampling <- c(HBsAG_prev_data[10,5])
prevvector_50to54yrs_clustersampling <- c(HBsAG_prev_data[11,5])

# Prevalence Data from Cluster Sampling from National Cross-Sectional Survey
prevvector_15to19yrs_crosssectional <- c(HBsAG_prev_data[4,6])
prevvector_20to24yrs_crosssectional <- c(HBsAG_prev_data[5,6])
prevvector_25to29yrs_crosssectional <- c(HBsAG_prev_data[6,6])
prevvector_30to34yrs_crosssectional <- c(HBsAG_prev_data[7,6])
prevvector_35to39yrs_crosssectional <- c(HBsAG_prev_data[8,6])
prevvector_40to44yrs_crosssectional <- c(HBsAG_prev_data[9,6])
prevvector_45to49yrs_crosssectional <- c(HBsAG_prev_data[10,6])
prevvector_50to54yrs_crosssectional <- c(HBsAG_prev_data[11,6])

# Data frame creation
prev_data_blooddonors_NBC <- data.frame(
  Year = rep(prevyearsvector_blooddonors_NBC, 8),
  Prevalence = c(
    prevvector_15to19yrs_blooddonors_NBC,
    prevvector_20to24yrs_blooddonors_NBC,
    prevvector_25to29yrs_blooddonors_NBC,
    prevvector_30to34yrs_blooddonors_NBC,
    prevvector_35to39yrs_blooddonors_NBC,
    prevvector_40to44yrs_blooddonors_NBC,
    prevvector_45to49yrs_blooddonors_NBC,
    prevvector_50to54yrs_blooddonors_NBC
  ),
  Age_Group = rep(c("15to19yrs", "20to24yrs", "25to29yrs", "30to34yrs", "35to39yrs",
                    "40to44yrs", "45to49yrs", "50to54yrs"), each = length(prevyearsvector_blooddonors_NBC))

)

prev_data_blooddonors_8provinces <- data.frame(
  Year = rep(prevyearsvector_blooddonors_8provinces, 8),
  Prevalence = c(#prevvector_0to4yrs_blooddonors_8provinces,
    prevvector_15to19yrs_blooddonors_8provinces,
    prevvector_20to24yrs_blooddonors_8provinces,
    prevvector_25to29yrs_blooddonors_8provinces,
    prevvector_30to34yrs_blooddonors_8provinces,
    prevvector_35to39yrs_blooddonors_8provinces,
    prevvector_40to44yrs_blooddonors_8provinces,
    prevvector_45to49yrs_blooddonors_8provinces,
    prevvector_50to54yrs_blooddonors_8provinces
  ),
  
  Age_Group = rep(c("15to19yrs", "20to24yrs", "25to29yrs", "30to34yrs", "35to39yrs",
                    "40to44yrs", "45to49yrs", "50to54yrs"), each = length(prevyearsvector_blooddonors_8provinces))
)

prev_data_clustersampling <- data.frame(
  Year = rep(prevyearsvector_clustersampling, 8),
  Prevalence = c(
    prevvector_15to19yrs_clustersampling,
    prevvector_20to24yrs_clustersampling,
    prevvector_25to29yrs_clustersampling,
    prevvector_30to34yrs_clustersampling,
    prevvector_35to39yrs_clustersampling,
    prevvector_40to44yrs_clustersampling,
    prevvector_45to49yrs_clustersampling,
    prevvector_50to54yrs_clustersampling

  ),
  Age_Group = rep(c("15to19yrs", "20to24yrs", "25to29yrs", "30to34yrs", "35to39yrs",
                    "40to44yrs", "45to49yrs", "50to54yrs"), each = length(prevyearsvector_clustersampling))

)

prev_data_crosssectional <- data.frame(
  Year = rep(prevyearsvector_crosssectional, 8),
  Prevalence = c(
    prevvector_15to19yrs_crosssectional,
    prevvector_20to24yrs_crosssectional,
    prevvector_25to29yrs_crosssectional,
    prevvector_30to34yrs_crosssectional,
    prevvector_35to39yrs_crosssectional,
    prevvector_40to44yrs_crosssectional,
    prevvector_45to49yrs_crosssectional,
    prevvector_50to54yrs_crosssectional
  ),
  Age_Group = rep(c("15to19yrs", "20to24yrs", "25to29yrs", "30to34yrs", "35to39yrs",
                    "40to44yrs", "45to49yrs", "50to54yrs"), each = length(prevyearsvector_crosssectional))
)

#Assigning levels so age groups stay in the same order
prev_data_blooddonors_NBC$Age_Group <- factor(prev_data_blooddonors_NBC$Age_Group,
                                              levels = c("15to19yrs",
                                                         "20to24yrs", "25to29yrs", "30to34yrs", "35to39yrs",
                                                         "40to44yrs", "45to49yrs", "50to54yrs"))
prev_data_blooddonors_8provinces$Age_Group <- factor(prev_data_blooddonors_8provinces$Age_Group,
                                                     levels = c("15to19yrs",
                                                                "20to24yrs", "25to29yrs", "30to34yrs", "35to39yrs",
                                                                "40to44yrs", "45to49yrs", "50to54yrs"))
prev_data_clustersampling$Age_Group <- factor(prev_data_clustersampling$Age_Group,
                                              levels = c("15to19yrs",
                                                         "20to24yrs", "25to29yrs", "30to34yrs", "35to39yrs",
                                                         "40to44yrs", "45to49yrs", "50to54yrs"
                                              ))
prev_data_crosssectional$Age_Group <- factor(prev_data_crosssectional$Age_Group,
                                             levels = c("15to19yrs",
                                                        "20to24yrs", "25to29yrs", "30to34yrs", "35to39yrs",
                                                        "40to44yrs", "45to49yrs", "50to54yrs"
                                             ))

# Define the specific age groups you want to plot
selected_age_groups <- c("15to19yrs","20to24yrs", "25to29yrs", "30to34yrs", "35to39yrs", "40to44yrs", "45to49yrs", "50to54yrs")

# Filter the data for the selected age groups
filtered_annual_summary <- annual_summary_prevalence_agegroups_long %>%
  filter(Age_Group %in% selected_age_groups)

# Plotting Age-Specific Prevalence for selected age groups
ggplot() +
  geom_line(data = filtered_annual_summary, aes(x = years, y = Prevalence, color = "Modelled Prevalence"), size = 1) +
  geom_point(data = prev_data_blooddonors_NBC, aes(x = Year, y =  Prevalence, color = "Blood Donors - NBC"), size = 2) +
  geom_point(data = prev_data_blooddonors_8provinces, aes(x = Year, y =  Prevalence, color = "Blood Donors - 8 Provinces"), size = 3, shape = 18) +
  geom_point(data = prev_data_clustersampling, aes(x = Year, y =  Prevalence, color = "Cluster Sampling"), size = 3, shape = 18) +
  geom_point(data = prev_data_crosssectional, aes(x = Year, y =  Prevalence, color = "Cross Sectional"), size = 3, shape = 18) +
  facet_wrap(~ Age_Group, scales = "fixed") +
  labs(title = "Prevalence Over Time by Age Group",
       x = "Year",
       y = "Population") +
  scale_x_continuous(limits = c(2000, 2025)) +
  scale_color_manual(values = c("Modelled Prevalence" = "gray",
                                "Blood Donors - NBC" = "black",
                                "Blood Donors - 8 Provinces" = "skyblue",
                                "Cluster Sampling" = "blue",
                                "Cross Sectional" = "navyblue")) +
  theme_minimal() +
  theme(legend.position = "bottom")



#########################################################
#RESULTS - SCENARIO ANALYSIS
#########################################################

# Creat lists of changes to vaccination coverage and proportion of population with access to treatment over time in the 5 different scenarios 
parms_variations <- list(
  list(vac_cov = c(0.613, 0.613, 0.613, 0.613, 0.613, 0.613, 0.613, 0.613, 0.613, 0.613, 0.613, 0.613),
       prop_treat = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)),  # Scenario 1
  list(vac_cov = c(0.663, 0.713, 0.763, 0.813, 0.863, 0.913, 0.963, 0.993, 0.993, 0.993, 0.993, 0.993),
       prop_treat = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)),  # Scenario 2
  list(vac_cov = c(0.613, 0.613, 0.613, 0.613, 0.613, 0.613, 0.613, 0.613, 0.613, 0.613, 0.613, 0.613),
       prop_treat = c(0.06, 0.11, 0.16, 0.21, 0.26, 0.31, 0.36, 0.41, 0.46, 0.51, 0.56, 0.61)),  # Scenario 3
  list(vac_cov = c(0.613, 0.638, 0.663, 0.688, 0.713, 0.763, 0.788, 0.813, 0.838, 0.863, 0.880, 0.913),
       prop_treat = c(0.035, 0.060, 0.085, 0.110, 0.135, 0.160, 0.185, 0.210, 0.235, 0.260, 0.285, 0.310)), # Scenario 4
  list(vac_cov = c(0.583, 0.553, 0.523, 0.493, 0.463, 0.433, 0.403, 0.373, 0.343, 0.313, 0.383, 0.253),
       prop_treat = c(0.06, 0.11, 0.16, 0.21, 0.26, 0.31, 0.36, 0.41, 0.46, 0.51, 0.56, 0.61))   # Scenario 5
)

# Define a function to run the model with modified parameters
run_model <- function(vac_cov, prop_treat) {

  # Modify the parameters
  parms$vac_cov_2029 <- vac_cov[1]
  parms$vac_cov_2035 <- vac_cov[2]
  parms$vac_cov_2041 <- vac_cov[3]
  parms$vac_cov_2047 <- vac_cov[4]
  parms$vac_cov_2053 <- vac_cov[5]
  parms$vac_cov_2059 <- vac_cov[6]
  parms$vac_cov_2065 <- vac_cov[7]
  parms$vac_cov_2071 <- vac_cov[8]
  parms$vac_cov_2077 <- vac_cov[9]
  parms$vac_cov_2083 <- vac_cov[10]
  parms$vac_cov_2089 <- vac_cov[11]
  parms$vac_cov_2095 <- vac_cov[12]
  parms$prop_treat_2029 <- prop_treat[1]
  parms$prop_treat_2035 <- prop_treat[2]
  parms$prop_treat_2041 <- prop_treat[3]
  parms$prop_treat_2047 <- prop_treat[4]
  parms$prop_treat_2053 <- prop_treat[5]
  parms$prop_treat_2059 <- prop_treat[6]
  parms$prop_treat_2065 <- prop_treat[7]
  parms$prop_treat_2071 <- prop_treat[8]
  parms$prop_treat_2077 <- prop_treat[9]
  parms$prop_treat_2083 <- prop_treat[10]
  parms$prop_treat_2089 <- prop_treat[11]
  parms$prop_treat_2095 <- prop_treat[12]

  # Run the differential equation solver (e.g., ode function from deSolve package)
  laopdr_HBV_output_scenarios <- ode(times = times_laopdr_HBV_model_bc1,
                                     y = init_laopdr_HBV_model_bc1,
                                     func = laopdr_HBV_model_bc1,
                                     parms = parms)

  return(laopdr_HBV_output_scenarios)
}

# Loop over the parameter variations and run the model
laopdr_HBV_output_scenarios_results <- lapply(parms_variations, function(p) {
                                        run_model(p$vac_cov, p$prop_treat)

})

#Store results of Scenarios in data frames 
scenario1dataframe <- as.data.frame(laopdr_HBV_output_scenarios_results[[1]])
scenario2dataframe <- as.data.frame(laopdr_HBV_output_scenarios_results[[2]])
scenario3dataframe <- as.data.frame(laopdr_HBV_output_scenarios_results[[3]])
scenario4dataframe <- as.data.frame(laopdr_HBV_output_scenarios_results[[4]])
scenario5dataframe <- as.data.frame(laopdr_HBV_output_scenarios_results[[5]])

#Assign names to columns using function defined earlier in the code 
column_names <- c("time", unlist(sapply(health_states, function(hs) paste(hs, age_group_names, sep = "_")))) # Generate the column names
colnames(scenario1dataframe) <- unlist(column_names) # Rename columns in df_combined

column_names <- c("time", unlist(sapply(health_states, function(hs) paste(hs, age_group_names, sep = "_")))) # Generate the column names
colnames(scenario2dataframe) <- unlist(column_names) # Rename columns in df_combined

column_names <- c("time", unlist(sapply(health_states, function(hs) paste(hs, age_group_names, sep = "_")))) # Generate the column names
colnames(scenario3dataframe) <- unlist(column_names) # Rename columns in df_combined

column_names <- c("time", unlist(sapply(health_states, function(hs) paste(hs, age_group_names, sep = "_")))) # Generate the column names
colnames(scenario4dataframe) <- unlist(column_names) # Rename columns in df_combined

column_names <- c("time", unlist(sapply(health_states, function(hs) paste(hs, age_group_names, sep = "_")))) # Generate the column names
colnames(scenario5dataframe) <- unlist(column_names) # Rename columns in df_combined

#Define functuons to sum columsn based on column names 
columns_to_sum_s1 <- colnames(scenario1dataframe)[2:length(scenario1dataframe)]
columns_to_sum_s2 <- colnames(scenario2dataframe)[2:length(scenario2dataframe)]
columns_to_sum_s3 <- colnames(scenario3dataframe)[2:length(scenario3dataframe)]
columns_to_sum_s4 <- colnames(scenario4dataframe)[2:length(scenario4dataframe)]
columns_to_sum_s5 <- colnames(scenario5dataframe)[2:length(scenario5dataframe)]

#Attached summed results to the end of existing data frames 
scenario1_results <- scenario1dataframe %>%
  mutate(
    CMort_total = rowSums(select(., columns_to_sum_s1[(22 * A + 1) : (23 * A)])),
    CInc_total = rowSums(select(., columns_to_sum_s1[(23 * A + 1) : (25 * A)]))
  )

scenario2_results <- scenario2dataframe %>%
  mutate(
    CMort_total = rowSums(select(., columns_to_sum_s2[(22 * A + 1) : (23 * A)])),
    CInc_total = rowSums(select(., columns_to_sum_s2[(23 * A + 1) : (25 * A)]))
  )

scenario3_results <- scenario3dataframe %>%
  mutate(
    CMort_total = rowSums(select(., columns_to_sum_s3[(22 * A + 1) : (23 * A)])),
    CInc_total = rowSums(select(., columns_to_sum_s3[(23 * A + 1) : (25 * A)]))
  )

scenario4_results <- scenario4dataframe %>%
  mutate(
    CMort_total = rowSums(select(., columns_to_sum_s4[(22 * A + 1) : (23 * A)])),
    CInc_total = rowSums(select(., columns_to_sum_s4[(23 * A + 1) : (25 * A)]))
  )

scenario5_results <- scenario5dataframe %>%
  mutate(
    CMort_total = rowSums(select(., columns_to_sum_s5[(22 * A + 1) : (23 * A)])),
    CInc_total = rowSums(select(., columns_to_sum_s5[(23 * A + 1) : (25 * A)]))
  )

# Convert df_laopdr_HBV_output_bc1 to tibble and calculate differences at each year
scenario1_results_final <- as_tibble(scenario1_results) %>%
  mutate(
    years = 1985 + floor(time/365.25),
    Mort = c(0, diff(CMort_total)), 
    Inc = c(0, diff(CInc_total))
  )

scenario2_results_final <- as_tibble(scenario2_results) %>%
  mutate(
    years = 1985 + floor(time/365.25),
    Mort = c(0, diff(CMort_total)),
    Inc = c(0, diff(CInc_total))
  )

scenario3_results_final <- as_tibble(scenario3_results) %>%
  mutate(
    years = 1985 + floor(time/365.25),
    Mort = c(0, diff(CMort_total)),
    Inc = c(0, diff(CInc_total))
  )

scenario4_results_final <- as_tibble(scenario4_results) %>%
  mutate(
    years = 1985 + floor(time/365.25),
    Mort = c(0, diff(CMort_total)),
    Inc = c(0, diff(CInc_total))
  )

scenario5_results_final <- as_tibble(scenario5_results) %>%
  mutate(
    years = 1985 + floor(time/365.25),
    Mort = c(0, diff(CMort_total)),
    Inc = c(0, diff(CInc_total))
  )

#Collect results for modelled period (2025 - 2100)
annual_summary_scenario1 <- scenario1_results_final %>%
  group_by(years) %>%
  filter(years >= 2025 & years <= 2100) %>%
  summarise(
    AnnualHepBDeaths = sum(Mort, na.rm = TRUE),
    AnnualIncidence = sum(Inc, na.rm = TRUE),
  )

annual_summary_scenario2 <- scenario2_results_final %>%
  group_by(years) %>%
  filter(years >= 2025 & years <= 2100) %>%
  summarise(
    AnnualHepBDeaths = sum(Mort, na.rm = TRUE),
    AnnualIncidence = sum(Inc, na.rm = TRUE)
  )

annual_summary_scenario3 <- scenario3_results_final %>%
  group_by(years) %>%
  filter(years >= 2025 & years <= 2100) %>%
  summarise(
    AnnualHepBDeaths = sum(Mort, na.rm = TRUE),
    AnnualIncidence = sum(Inc, na.rm = TRUE)
  )

annual_summary_scenario4 <- scenario4_results_final %>%
  group_by(years) %>%
  filter(years >= 2025 & years <= 2100) %>%
  summarise(
    AnnualHepBDeaths = sum(Mort, na.rm = TRUE),
    AnnualIncidence = sum(Inc, na.rm = TRUE)
  )

annual_summary_scenario5 <- scenario5_results_final %>%
  group_by(years) %>%
  filter(years >= 2025 & years <= 2100) %>%
  summarise(
    AnnualHepBDeaths = sum(Mort, na.rm = TRUE),
    AnnualIncidence = sum(Inc, na.rm = TRUE)
  )

# Plotting Annual Mortality Over Time
ggplot() +
  geom_line(data = annual_summary_scenario1, aes(x = years, y = AnnualHepBDeaths, color = "Scenario 1"), size = 1) +
  geom_line(data = annual_summary_scenario2, aes(x = years, y = AnnualHepBDeaths, color = "Scenario 2"), size = 1) +
  geom_line(data = annual_summary_scenario3, aes(x = years, y = AnnualHepBDeaths, color = "Scenario 3"), size = 1) +
  geom_line(data = annual_summary_scenario4, aes(x = years, y = AnnualHepBDeaths, color = "Scenario 4"), size = 1) +
  geom_line(data = annual_summary_scenario5, aes(x = years, y =  AnnualHepBDeaths, color = "Scenario 5"), size = 1) +
  labs(title = "HBV Mortality Over Time",
       x = "Year",
       y = "Annual HBV Mortality") +
  scale_color_manual(values = c("Scenario 1" = "black",
                                "Scenario 2" = "skyblue",
                                "Scenario 3" = "purple",
                                "Scenario 4" = "yellow",
                                "Scenario 5" = "green"
                                )) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Plotting Annual Incidence Over Time
ggplot() +
  geom_line(data = annual_summary_scenario1, aes(x = years, y =  AnnualIncidence, color = "Scenario 1"), size = 1) +
  geom_line(data = annual_summary_scenario2, aes(x = years, y =  AnnualIncidence, color = "Scenario 2"), size = 1) +
  geom_line(data = annual_summary_scenario3, aes(x = years, y =  AnnualIncidence, color = "Scenario 3"), size = 1) +
  geom_line(data = annual_summary_scenario4, aes(x = years, y =  AnnualIncidence, color = "Scenario 4"), size = 1) +
  geom_line(data = annual_summary_scenario5, aes(x = years, y =  AnnualIncidence, color = "Scenario 5"), size = 1) +
  labs(title = "HBV Incidence Over Time",
       x = "Year",
       y = "Annual HBV Mortality") +
  scale_color_manual(values = c("Scenario 1" = "black",
                                "Scenario 2" = "skyblue",
                                "Scenario 3" = "purple",
                                "Scenario 4" = "yellow",
                                "Scenario 5" = "green"
                                )) +
  theme_minimal() +
  theme(legend.position = "bottom")





















