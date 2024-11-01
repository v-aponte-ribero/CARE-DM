# Load packages and files --------------------------------------------------

#Load packages
if (!require(pacman)) install.packages('pacman')
library(pacman)

p_load(tidyverse,mice)

#Load data
source("YOUR FILEPATH")

#Set seed
set.seed(2023)


# Define parameters -------------------------------------------------------

seed <- 202302
maxit <- 20
m <- 10

#Variables used
vars <- c("ID",
          "Study",
          "Country", #Included as auxiliary variable
          "Region_mod",
          "Region_high",
          "cAge",
          "Gender",
          "Black", #Systematically missing from SHARE, but included as auxiliary variable for data from other studies
          "DiabetesDurCat",
          "cDiabetesAge", #Included to allow calculation of SCORE2-Diabetes
          "AntiDiab_adj",
          "Insulin", #Systematically missing from SHARE, but included as auxiliary variable for data from other studies
          "cHbA1c",
          "clneGFR",
          "Smoker",
          "PrevSmoker", #Included as auxiliary variable
          "CurDrink",
          "cBMI",
          "TreatHyperten_adj",
          "cSBP", #Systematically missing from SHARE, but imputed using other Studies to allow calculation of SCORE2-Diabetes
          "cTOTC",
          "cHDLC",
          "Choltreat",
          "CVD_event","TTE_CVDDeath","na_CVD",
          "Other_death","TTE_Death","na_OtherDeath"
)

#Build predictor matrix
predmat <- matrix(1,nrow=length(vars),ncol=length(vars),dimnames=list(vars,vars))
diag(predmat) <- 0

predmat["ID",] <- 0 
predmat[,"ID"] <- 0
predmat["Study",] <- 0
predmat["cAge",] <- 0
predmat["Gender",] <- 0

predmat["PrevSmoker",] <- 0
predmat[,"PrevSmoker"] <- 0

predmat["cDiabetesAge","DiabetesDurCat"] <- 0
predmat["DiabetesDurCat","cDiabetesAge"] <- 0

predmat["AntiDiab_adj","Insulin"] <- 0

predmat["CVD_event",] <- 0
predmat["Other_death",] <- 0
predmat["TTE_CVDDeath",] <- 0
predmat[,"TTE_CVDDeath"] <- 0
predmat["TTE_Death",] <- 0
predmat[,"TTE_Death"] <- 0
predmat["na_CVD",] <- 0
predmat["na_OtherDeath",] <- 0

#Defined imputation methods for each variable
meth <- c(
  ID="",
  Study="",
  Country="",
  Region_mod="",
  Region_high="",
  cAge="",
  Gender="",
  Black="pmm",
  DiabetesDurCat="pmm",
  cDiabetesAge="pmm",
  AntiDiab_adj="pmm",
  Insulin="pmm",
  cHbA1c="pmm",
  clneGFR="pmm",
  Smoker="pmm",
  PrevSmoker="pmm",
  CurDrink="pmm",
  cBMI="pmm",
  TreatHyperten_adj="pmm",
  cSBP="pmm",
  cTOTC="pmm",
  cHDLC="pmm",
  Choltreat="pmm",
  TTE_CVDDeath="",
  CVD_event="",
  na_CVD="",
  Other_death="",
  TTE_Death="",
  na_OtherDeath=""
)

# Prepare dataset
dat <- full_dat %>% select(any_of(vars))

## Add nelson aalen estimator
dat$na_CVD <- nelsonaalen(data = dat,
                          timevar = "TTE_CVDDeath",
                          statusvar = "CVD_event")

dat$na_OtherDeath <- nelsonaalen(data = dat,
                                 timevar = "TTE_Death",
                                 statusvar = "Other_death")

# Update predictor matrix, methods, and visit sequence according to degree of missingness
predmat1 <- predmat
meth1 <- meth

missing_table <- function(dat){
  cbind("# NA" = sort(colSums(is.na(dat))),
        "% NA" = round(sort(colMeans(is.na(dat))) * 100, 2))
}

mt <- as.data.frame(missing_table(dat))
mt <- rownames_to_column(mt)
mtNA <- mt[mt[2]==0,]
visseq <- mt[,1]

for(n in 1:length(mtNA$rowname)){
  predmat1[mtNA$rowname[n],] <- 0
  meth1[mtNA$rowname[n]] <- ""
}


# Impute ------------------------------------------------------------------

imp_tot <- mice(dat,
                predictorMatrix = predmat1,
                visitSequence = visseq,
                method = meth1,
                maxit = maxit, m = m,
                seed = seed)

imp_tot$loggedEvents
plot(imp_tot, layout = c(6, 6))
densityplot(imp_tot)
bwplot(imp_tot)

