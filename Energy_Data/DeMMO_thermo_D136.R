install.packages("CHNOSZ")


#install developer update:
#> install.packages("/Users/Caitlin/Downloads/CHNOSZ_1.1.3-63.tgz", repos = NULL)
library(CHNOSZ)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(plyr)
library(grid)

#create thermo database
data(thermo)

#add pyrolusite to database
pyrolusite <- mod.obigt("pyrolusite", G=-111100, H=-124283, S=12.61, V=17.3, formula="MnO2", state="cr", a1.a=12.55, a2.b=9.761, a3.c=-2.105)

#add ferrihydrite to database
ferrihydrite <- mod.obigt("ferrihydrite", G=-111200, H=-127800, S=16.7, V=20.88, formula="FeOOH", state="cr", a1.a=8.70, a2.b=36.71, a3.c=-1.0146)

#add manganite to database 
manganite <- mod.obigt("manganite", G=-133300, formula="MnOOH", state="cr")

#store database in a dataframe for reference 
thermo_db <- data.frame(thermo$obigt)

#set temperature units to Celcius
T.units("K")
E.units("J")

#import DeMMO mineral reactions
reactions <- read.csv("DeMMO_CHNOSZ.csv", header = TRUE)

#import DeMMO fluid activities 
activities <- read.csv("DeMMO_SpecE8.csv", header=TRUE)
colnames(activities) <- c("Site","Ca+2", "acetate","CH4","Fe+2",	"H+",	"H2",	"HCO3-","HS-","Mn+2","NH4+","NO2-","NO3-","SO4-2")

DeMMO_thermo<- cbind(Site=c("DeMMO1", "DeMMO2", "DeMMO3", "DeMMO4", "DeMMO5", "DeMMO6"))
                     
#create the function DeMMO_logK to calulate logK values for every reaction for six DeMMO sites 
DeMMO_logK <- function(data){
  rxn.number <- data[1]
  react.prod <- as.vector(unlist(c(data[4:12])))
  react.prod <-react.prod[react.prod != ""] 
  react.coeffs <- as.vector(na.omit(as.numeric(data[14:22])))
  react.states <- as.vector(unlist(c(data[23:31])))
  react.states <- react.states[react.states != ""]
  DeMMO_T <- c(283.45, 285.55, 289.35, 295.65, 304.85, 294.65)
  rxn <- subcrt(react.prod, react.coeffs, react.states, T=DeMMO_T)
  DeMMO_thermo[rxn.number] <- rxn$out$logK
}

#store the output of the function applied to the reactions dataframe in new dataframe DeMMO_thermo
DeMMO_thermo <- apply(reactions, 1, DeMMO_logK)

#calculate logK at 25C for manganite rxns
manganite_reactions <- reactions[c(14,22,31),]

DeMMO_manganite.rxns_logK <- function(data){
  rxn.number <- data[1]
  react.prod <- as.vector(unlist(c(data[4:12])))
  react.prod <-react.prod[react.prod != ""] 
  react.coeffs <- as.vector(na.omit(as.numeric(data[14:22])))
  react.states <- as.vector(unlist(c(data[23:31])))
  react.states <- react.states[react.states != ""]
  #DeMMO_T <- c(283.45, 285.55, 289.35, 295.65, 304.85, 294.65)
  rxn <- subcrt(react.prod, react.coeffs, react.states, T=c(rep(298.15, 6)))
  DeMMO_thermo[rxn.number] <- rxn$out$logK
}

DeMMO_manganite.logK <- apply(manganite_reactions, 1,DeMMO_manganite.rxns_logK)

#add manganite values to DeMMO_thermo
DeMMO_thermo[,14] <- DeMMO_manganite.logK[,1]
DeMMO_thermo[,22] <- DeMMO_manganite.logK[,2]
DeMMO_thermo[,31] <- DeMMO_manganite.logK[,3]

#format the DeMMO_thermo dataframe
colnames(DeMMO_thermo) <- c(1:33)
colnames(DeMMO_thermo) <- paste("rxn", colnames(DeMMO_thermo), "logK", sep = ".")
DeMMO_thermo<- cbind(Site=c("DeMMO1", "DeMMO2", "DeMMO3", "DeMMO4", "DeMMO5", "DeMMO6"), DeMMO_thermo)


DeMMO1_logQ <- function(data){
  rxn.reactants <- as.vector(unlist(c(data[4:7])))
  rxn.reactants<- gsub("[\\+]", "[+]", rxn.reactants)
  rxn.reactants <- gsub("[\\-]", "[-]", rxn.reactants)
  rxn.reactants <-rxn.reactants[rxn.reactants != ""]
  rxn.react.act <- activities[,grepl(paste(rxn.reactants, collapse = "|"), colnames(activities))]
  rxn.react.act.len <- length(rxn.react.act)
  if (rxn.react.act.len > 1) {
    rxn.react.act.index <- match(grep(paste(rxn.reactants, collapse = "|"), names(activities), value=TRUE), as.vector(unlist(c(data[4:7]))))
    rxn.react.coeff <- abs(as.vector(na.omit(as.numeric(data[14:17]))))[rxn.react.act.index]
    react.vec.len <- length(as.vector(grep(paste(rxn.reactants, collapse = "|"), names(activities), value=TRUE)))
    if (react.vec.len > 1) {
      log.act.reactants <- rxn.react.act[1,]*rxn.react.coeff
    } else{
      log.act.reactants <- rxn.react.act[1]*rxn.react.coeff
    }
  }else{
    log.act.reactants <- 1
  }


  rxn.products <- as.vector(unlist(c(data[8:12])))
  rxn.products<- gsub("[\\+]", "[+]", rxn.products)
  rxn.products <- gsub("[\\-]", "[-]", rxn.products)
  rxn.products <-rxn.products[rxn.products != ""] 
  rxn.prod.act <- activities[,grepl(paste(rxn.products, collapse = "|"), colnames(activities))]
  rxn.prod.act.index <- match(grep(paste(rxn.products, collapse = "|"), names(activities), value=TRUE), as.vector(unlist(c(data[8:12]))))
  rxn.prod.coeff <- abs(as.vector(na.omit(as.numeric(data[18:22]))))[rxn.prod.act.index]
  prod.vec.len <- length(as.vector(grep(paste(rxn.products, collapse = "|"), names(activities), value=TRUE)))
  if (prod.vec.len > 1) {
    log.act.products <- rxn.prod.act[1,]*rxn.prod.coeff
  } else{
    log.act.products <- rxn.prod.act[1]*rxn.prod.coeff
  }
  rxn.number <- data[1]
  sum(log.act.products) + sum(-(log.act.reactants))
}



DeMMO1_logQ <- apply(reactions, 1, DeMMO1_logQ)
DeMMO_thermo_logQ<- as.data.frame(DeMMO1_logQ)

#calculate logQ for all rxns at DeMMO2
DeMMO2_logQ <- function(data){
  rxn.reactants <- as.vector(unlist(c(data[4:7])))
  rxn.reactants<- gsub("[\\+]", "[+]", rxn.reactants)
  rxn.reactants <- gsub("[\\-]", "[-]", rxn.reactants)
  rxn.reactants <-rxn.reactants[rxn.reactants != ""]
  rxn.react.act <- activities[,grepl(paste(rxn.reactants, collapse = "|"), colnames(activities))]
  rxn.react.act.len <- length(rxn.react.act)
  if (rxn.react.act.len > 1) {
    rxn.react.act.index <- match(grep(paste(rxn.reactants, collapse = "|"), names(activities), value=TRUE), as.vector(unlist(c(data[4:7]))))
    rxn.react.coeff <- abs(as.vector(na.omit(as.numeric(data[14:17]))))[rxn.react.act.index]
    react.vec.len <- length(as.vector(grep(paste(rxn.reactants, collapse = "|"), names(activities), value=TRUE)))
    if (react.vec.len > 1) {
      log.act.reactants <- rxn.react.act[2,]*rxn.react.coeff
    } else{
      log.act.reactants <- rxn.react.act[2]*rxn.react.coeff
    }
  }else{
    log.act.reactants <- 1
  }
  
  
  rxn.products <- as.vector(unlist(c(data[8:12])))
  rxn.products<- gsub("[\\+]", "[+]", rxn.products)
  rxn.products <- gsub("[\\-]", "[-]", rxn.products)
  rxn.products <-rxn.products[rxn.products != ""] 
  rxn.prod.act <- activities[,grepl(paste(rxn.products, collapse = "|"), colnames(activities))]
  rxn.prod.act.index <- match(grep(paste(rxn.products, collapse = "|"), names(activities), value=TRUE), as.vector(unlist(c(data[8:12]))))
  rxn.prod.coeff <- abs(as.vector(na.omit(as.numeric(data[18:22]))))[rxn.prod.act.index]
  prod.vec.len <- length(as.vector(grep(paste(rxn.products, collapse = "|"), names(activities), value=TRUE)))
  if (prod.vec.len > 1) {
    log.act.products <- rxn.prod.act[2,]*rxn.prod.coeff
  } else{
    log.act.products <- rxn.prod.act[2]*rxn.prod.coeff
  }
  rxn.number <- data[1]
  sum(log.act.products) + sum(-(log.act.reactants))
}



DeMMO2_logQ <- apply(reactions, 1, DeMMO2_logQ)
DeMMO_thermo_logQ$DeMMO2_logQ <- DeMMO2_logQ

DeMMO2_logQ <- function(data){
  rxn.reactants <- as.vector(unlist(c(data[4:7])))
  rxn.reactants<- gsub("[\\+]", "[+]", rxn.reactants)
  rxn.reactants <- gsub("[\\-]", "[-]", rxn.reactants)
  rxn.reactants <-rxn.reactants[rxn.reactants != ""]
  rxn.react.act <- activities[,grepl(paste(rxn.reactants, collapse = "|"), colnames(activities))]
  rxn.react.act.len <- length(rxn.react.act)
  if (rxn.react.act.len > 1) {
    rxn.react.act.index <- match(grep(paste(rxn.reactants, collapse = "|"), names(activities), value=TRUE), as.vector(unlist(c(data[4:7]))))
    rxn.react.coeff <- abs(as.vector(na.omit(as.numeric(data[14:17]))))[rxn.react.act.index]
    react.vec.len <- length(as.vector(grep(paste(rxn.reactants, collapse = "|"), names(activities), value=TRUE)))
    if (react.vec.len > 1) {
      log.act.reactants <- rxn.react.act[2,]*rxn.react.coeff
    } else{
      log.act.reactants <- rxn.react.act[2]*rxn.react.coeff
    }
  }else{
    log.act.reactants <- 1
  }
  
  
  rxn.products <- as.vector(unlist(c(data[8:12])))
  rxn.products<- gsub("[\\+]", "[+]", rxn.products)
  rxn.products <- gsub("[\\-]", "[-]", rxn.products)
  rxn.products <-rxn.products[rxn.products != ""] 
  rxn.prod.act <- activities[,grepl(paste(rxn.products, collapse = "|"), colnames(activities))]
  rxn.prod.act.index <- match(grep(paste(rxn.products, collapse = "|"), names(activities), value=TRUE), as.vector(unlist(c(data[8:12]))))
  rxn.prod.coeff <- abs(as.vector(na.omit(as.numeric(data[18:22]))))[rxn.prod.act.index]
  prod.vec.len <- length(as.vector(grep(paste(rxn.products, collapse = "|"), names(activities), value=TRUE)))
  if (prod.vec.len > 1) {
    log.act.products <- rxn.prod.act[2,]*rxn.prod.coeff
  } else{
    log.act.products <- rxn.prod.act[2]*rxn.prod.coeff
  }
  rxn.number <- data[1]
  sum(log.act.products) + sum(-(log.act.reactants))
}



DeMMO2_logQ <- apply(reactions, 1, DeMMO2_logQ)
DeMMO_thermo_logQ$DeMMO2_logQ <- DeMMO2_logQ


#calculate logQ for all reactions at DeMMO3
DeMMO3_logQ <- function(data){
  rxn.reactants <- as.vector(unlist(c(data[4:7])))
  rxn.reactants<- gsub("[\\+]", "[+]", rxn.reactants)
  rxn.reactants <- gsub("[\\-]", "[-]", rxn.reactants)
  rxn.reactants <-rxn.reactants[rxn.reactants != ""]
  rxn.react.act <- activities[,grepl(paste(rxn.reactants, collapse = "|"), colnames(activities))]
  rxn.react.act.len <- length(rxn.react.act)
  if (rxn.react.act.len > 1) {
    rxn.react.act.index <- match(grep(paste(rxn.reactants, collapse = "|"), names(activities), value=TRUE), as.vector(unlist(c(data[4:7]))))
    rxn.react.coeff <- abs(as.vector(na.omit(as.numeric(data[14:17]))))[rxn.react.act.index]
    react.vec.len <- length(as.vector(grep(paste(rxn.reactants, collapse = "|"), names(activities), value=TRUE)))
    if (react.vec.len > 1) {
      log.act.reactants <- rxn.react.act[3,]*rxn.react.coeff
    } else{
      log.act.reactants <- rxn.react.act[3]*rxn.react.coeff
    }
  }else{
    log.act.reactants <- 1
  }
  
  
  rxn.products <- as.vector(unlist(c(data[8:12])))
  rxn.products<- gsub("[\\+]", "[+]", rxn.products)
  rxn.products <- gsub("[\\-]", "[-]", rxn.products)
  rxn.products <-rxn.products[rxn.products != ""] 
  rxn.prod.act <- activities[,grepl(paste(rxn.products, collapse = "|"), colnames(activities))]
  rxn.prod.act.index <- match(grep(paste(rxn.products, collapse = "|"), names(activities), value=TRUE), as.vector(unlist(c(data[8:12]))))
  rxn.prod.coeff <- abs(as.vector(na.omit(as.numeric(data[18:22]))))[rxn.prod.act.index]
  prod.vec.len <- length(as.vector(grep(paste(rxn.products, collapse = "|"), names(activities), value=TRUE)))
  if (prod.vec.len > 1) {
    log.act.products <- rxn.prod.act[3,]*rxn.prod.coeff
  } else{
    log.act.products <- rxn.prod.act[3]*rxn.prod.coeff
  }
  rxn.number <- data[1]
  sum(log.act.products) + sum(-(log.act.reactants))
}



DeMMO3_logQ <- apply(reactions, 1, DeMMO3_logQ)
DeMMO_thermo_logQ$DeMMO3_logQ <- DeMMO3_logQ


#calculate logQ for all reactions at DeMMO4
DeMMO4_logQ <- function(data){
  rxn.reactants <- as.vector(unlist(c(data[4:7])))
  rxn.reactants<- gsub("[\\+]", "[+]", rxn.reactants)
  rxn.reactants <- gsub("[\\-]", "[-]", rxn.reactants)
  rxn.reactants <-rxn.reactants[rxn.reactants != ""]
  rxn.react.act <- activities[,grepl(paste(rxn.reactants, collapse = "|"), colnames(activities))]
  rxn.react.act.len <- length(rxn.react.act)
  if (rxn.react.act.len > 1) {
    rxn.react.act.index <- match(grep(paste(rxn.reactants, collapse = "|"), names(activities), value=TRUE), as.vector(unlist(c(data[4:7]))))
    rxn.react.coeff <- abs(as.vector(na.omit(as.numeric(data[14:17]))))[rxn.react.act.index]
    react.vec.len <- length(as.vector(grep(paste(rxn.reactants, collapse = "|"), names(activities), value=TRUE)))
    if (react.vec.len > 1) {
      log.act.reactants <- rxn.react.act[4,]*rxn.react.coeff
    } else{
      log.act.reactants <- rxn.react.act[4]*rxn.react.coeff
    }
  }else{
    log.act.reactants <- 1
  }
  
  
  rxn.products <- as.vector(unlist(c(data[8:12])))
  rxn.products<- gsub("[\\+]", "[+]", rxn.products)
  rxn.products <- gsub("[\\-]", "[-]", rxn.products)
  rxn.products <-rxn.products[rxn.products != ""] 
  rxn.prod.act <- activities[,grepl(paste(rxn.products, collapse = "|"), colnames(activities))]
  rxn.prod.act.index <- match(grep(paste(rxn.products, collapse = "|"), names(activities), value=TRUE), as.vector(unlist(c(data[8:12]))))
  rxn.prod.coeff <- abs(as.vector(na.omit(as.numeric(data[18:22]))))[rxn.prod.act.index]
  prod.vec.len <- length(as.vector(grep(paste(rxn.products, collapse = "|"), names(activities), value=TRUE)))
  if (prod.vec.len > 1) {
    log.act.products <- rxn.prod.act[4,]*rxn.prod.coeff
  } else{
    log.act.products <- rxn.prod.act[4]*rxn.prod.coeff
  }
  rxn.number <- data[1]
  sum(log.act.products) + sum(-(log.act.reactants))
}



DeMMO4_logQ <- apply(reactions, 1, DeMMO4_logQ)
DeMMO_thermo_logQ$DeMMO4_logQ <- DeMMO4_logQ


#calculate logQ for all reactions at DeMMO5
DeMMO5_logQ <- function(data){
  rxn.reactants <- as.vector(unlist(c(data[4:7])))
  rxn.reactants<- gsub("[\\+]", "[+]", rxn.reactants)
  rxn.reactants <- gsub("[\\-]", "[-]", rxn.reactants)
  rxn.reactants <-rxn.reactants[rxn.reactants != ""]
  rxn.react.act <- activities[,grepl(paste(rxn.reactants, collapse = "|"), colnames(activities))]
  rxn.react.act.len <- length(rxn.react.act)
  if (rxn.react.act.len > 1) {
    rxn.react.act.index <- match(grep(paste(rxn.reactants, collapse = "|"), names(activities), value=TRUE), as.vector(unlist(c(data[4:7]))))
    rxn.react.coeff <- abs(as.vector(na.omit(as.numeric(data[14:17]))))[rxn.react.act.index]
    react.vec.len <- length(as.vector(grep(paste(rxn.reactants, collapse = "|"), names(activities), value=TRUE)))
    if (react.vec.len > 1) {
      log.act.reactants <- rxn.react.act[5,]*rxn.react.coeff
    } else{
      log.act.reactants <- rxn.react.act[5]*rxn.react.coeff
    }
  }else{
    log.act.reactants <- 1
  }
  
  
  rxn.products <- as.vector(unlist(c(data[8:12])))
  rxn.products<- gsub("[\\+]", "[+]", rxn.products)
  rxn.products <- gsub("[\\-]", "[-]", rxn.products)
  rxn.products <-rxn.products[rxn.products != ""] 
  rxn.prod.act <- activities[,grepl(paste(rxn.products, collapse = "|"), colnames(activities))]
  rxn.prod.act.index <- match(grep(paste(rxn.products, collapse = "|"), names(activities), value=TRUE), as.vector(unlist(c(data[8:12]))))
  rxn.prod.coeff <- abs(as.vector(na.omit(as.numeric(data[18:22]))))[rxn.prod.act.index]
  prod.vec.len <- length(as.vector(grep(paste(rxn.products, collapse = "|"), names(activities), value=TRUE)))
  if (prod.vec.len > 1) {
    log.act.products <- rxn.prod.act[5,]*rxn.prod.coeff
  } else{
    log.act.products <- rxn.prod.act[5]*rxn.prod.coeff
  }
  rxn.number <- data[1]
  sum(log.act.products) + sum(-(log.act.reactants))
}



DeMMO5_logQ <- apply(reactions, 1, DeMMO5_logQ)
DeMMO_thermo_logQ$DeMMO5_logQ <- DeMMO5_logQ


#calculate logQ for all reactions at DeMMO6
DeMMO6_logQ <- function(data){
  rxn.reactants <- as.vector(unlist(c(data[4:7])))
  rxn.reactants<- gsub("[\\+]", "[+]", rxn.reactants)
  rxn.reactants <- gsub("[\\-]", "[-]", rxn.reactants)
  rxn.reactants <-rxn.reactants[rxn.reactants != ""]
  rxn.react.act <- activities[,grepl(paste(rxn.reactants, collapse = "|"), colnames(activities))]
  rxn.react.act.len <- length(rxn.react.act)
  if (rxn.react.act.len > 1) {
    rxn.react.act.index <- match(grep(paste(rxn.reactants, collapse = "|"), names(activities), value=TRUE), as.vector(unlist(c(data[4:7]))))
    rxn.react.coeff <- abs(as.vector(na.omit(as.numeric(data[14:17]))))[rxn.react.act.index]
    react.vec.len <- length(as.vector(grep(paste(rxn.reactants, collapse = "|"), names(activities), value=TRUE)))
    if (react.vec.len > 1) {
      log.act.reactants <- rxn.react.act[6,]*rxn.react.coeff
    } else{
      log.act.reactants <- rxn.react.act[6]*rxn.react.coeff
    }
  }else{
    log.act.reactants <- 1
  }
  
  
  rxn.products <- as.vector(unlist(c(data[8:12])))
  rxn.products<- gsub("[\\+]", "[+]", rxn.products)
  rxn.products <- gsub("[\\-]", "[-]", rxn.products)
  rxn.products <-rxn.products[rxn.products != ""] 
  rxn.prod.act <- activities[,grepl(paste(rxn.products, collapse = "|"), colnames(activities))]
  rxn.prod.act.index <- match(grep(paste(rxn.products, collapse = "|"), names(activities), value=TRUE), as.vector(unlist(c(data[8:12]))))
  rxn.prod.coeff <- abs(as.vector(na.omit(as.numeric(data[18:22]))))[rxn.prod.act.index]
  prod.vec.len <- length(as.vector(grep(paste(rxn.products, collapse = "|"), names(activities), value=TRUE)))
  if (prod.vec.len > 1) {
    log.act.products <- rxn.prod.act[6,]*rxn.prod.coeff
  } else{
    log.act.products <- rxn.prod.act[6]*rxn.prod.coeff
  }
  rxn.number <- data[1]
  sum(log.act.products) + sum(-(log.act.reactants))
}



DeMMO6_logQ <- apply(reactions, 1, DeMMO6_logQ)
DeMMO_thermo_logQ$DeMMO6_logQ <- DeMMO6_logQ
DeMMO_thermo <- t(DeMMO_thermo)
DeMMO_thermo <- DeMMO_thermo[-1,]
colnames(DeMMO_thermo) <- c(1:6)
colnames(DeMMO_thermo) <- paste("DeMMO", colnames(DeMMO_thermo), "logK", sep = ".")



DeMMO_thermo <- cbind(DeMMO_thermo, DeMMO_thermo_logQ)
row.names(DeMMO_thermo) <- c(1:33)
DeMMO_thermo <- as.data.frame(DeMMO_thermo)
DeMMO_thermo$e.transfer <- reactions$e.transfer







DeMMO1_DeltaG <- function(data){
  DeMMO1_DeltaG <- -2.303*8.314*283.45*(as.numeric(as.vector(data[1]))-as.numeric(as.vector(data[7])))
}

DeMMO_thermo$DeMMO1_DeltaG <- apply(DeMMO_thermo, 1, DeMMO1_DeltaG)

DeMMO2_DeltaG <- function(data){
DeMMO2_DeltaG <- -2.303*8.314*285.55*(as.numeric(as.vector(data[2]))-as.numeric(as.vector(data[8])))
}
DeMMO_thermo$DeMMO2_DeltaG <- apply(DeMMO_thermo, 1, DeMMO2_DeltaG)

DeMMO3_DeltaG <- function(data){
DeMMO3_DeltaG <- -2.303*8.314*289.35*(as.numeric(as.vector(data[3]))-as.numeric(as.vector(data[9])))
}
DeMMO_thermo$DeMMO3_DeltaG <- apply(DeMMO_thermo, 1, DeMMO3_DeltaG)

DeMMO4_DeltaG <- function(data){
DeMMO4_DeltaG <- -2.303*8.314*295.65*(as.numeric(as.vector(data[4]))-as.numeric(as.vector(data[10])))
}
DeMMO_thermo$DeMMO4_DeltaG <- apply(DeMMO_thermo, 1, DeMMO4_DeltaG)

DeMMO5_DeltaG <- function(data){
DeMMO5_DeltaG <- -2.303*8.314*304.85*(as.numeric(as.vector(data[5]))-as.numeric(as.vector(data[11])))
}
DeMMO_thermo$DeMMO5_DeltaG <- apply(DeMMO_thermo, 1, DeMMO5_DeltaG)

DeMMO6_DeltaG <- function(data){
DeMMO6_DeltaG <- -2.303*8.314*294.65*(as.numeric(as.vector(data[6]))-as.numeric(as.vector(data[12])))
}
DeMMO_thermo$DeMMO6_DeltaG <- apply(DeMMO_thermo, 1, DeMMO6_DeltaG)


#normalize to the number of e- transferred in units of kJ/mole e- transferred 
DeMMO1_DeltaG_norm_fun <- function(data){
  as.numeric(data[14])/(as.numeric(as.vector(data[13]))*1000)
}

DeMMO_thermo$DeMMO1_DeltaG_norm <- apply(DeMMO_thermo, 1, DeMMO1_DeltaG_norm_fun)

DeMMO2_DeltaG_norm_fun <- function(data){
  as.numeric(data[15])/(as.numeric(as.vector(data[13]))*1000)
}

DeMMO_thermo$DeMMO2_DeltaG_norm <- apply(DeMMO_thermo, 1, DeMMO2_DeltaG_norm_fun)

DeMMO3_DeltaG_norm_fun <- function(data){
  as.numeric(data[16])/(as.numeric(as.vector(data[13]))*1000)
}

DeMMO_thermo$DeMMO3_DeltaG_norm <- apply(DeMMO_thermo, 1, DeMMO3_DeltaG_norm_fun)

DeMMO4_DeltaG_norm_fun <- function(data){
  as.numeric(data[17])/(as.numeric(as.vector(data[13]))*1000)
}

DeMMO_thermo$DeMMO4_DeltaG_norm <- apply(DeMMO_thermo, 1, DeMMO4_DeltaG_norm_fun)

DeMMO5_DeltaG_norm_fun <- function(data){
  as.numeric(data[18])/(as.numeric(as.vector(data[13]))*1000)
}

DeMMO_thermo$DeMMO5_DeltaG_norm <- apply(DeMMO_thermo, 1, DeMMO5_DeltaG_norm_fun)

DeMMO6_DeltaG_norm_fun <- function(data){
  as.numeric(data[19])/(as.numeric(as.vector(data[13]))*1000)
}

DeMMO_thermo$DeMMO6_DeltaG_norm <- apply(DeMMO_thermo, 1, DeMMO6_DeltaG_norm_fun)


###calculate energy density by normalizing to e-donor (ED) or e-acceptor (EA)

reactions <- cbind(DeMMO_thermo[,14:19], reactions)
DeMMO_E_dens <- data.frame("rxn.number"=c(1:33))
DeMMO1_ED_dens_fun <- function(data){
  e.donor <- as.vector(data[9])
  e.donor<- gsub("[\\+]", "[+]", e.donor)
  e.donor <- gsub("[\\-]", "[-]", e.donor)
  e.donor <-e.donor[e.donor != ""]
  rxn.number <- reactions[7]
  e.donor.act <- as.vector(unlist(c(activities[,grepl(paste(e.donor, collapse = "|"), colnames(activities))])))
  e.donor.act.len <- length(e.donor.act[1])
  if (e.donor.act.len > 0){
    print(log10(-(10^e.donor.act[1])*(as.numeric(data[1])/abs(as.numeric(data[21])))))
  }
}

DeMMO_E_dens$DeMMO1_EDdens <- apply(reactions, 1, DeMMO1_ED_dens_fun)
  
DeMMO2_ED_dens_fun <- function(data){
  e.donor <- as.vector(data[9])
  e.donor<- gsub("[\\+]", "[+]", e.donor)
  e.donor <- gsub("[\\-]", "[-]", e.donor)
  e.donor <-e.donor[e.donor != ""]
  rxn.number <- reactions[7]
  e.donor.act <- as.vector(unlist(c(activities[,grepl(paste(e.donor, collapse = "|"), colnames(activities))])))
  e.donor.act.len <- length(e.donor.act[2])
  if (e.donor.act.len > 0){
    print(log10(-(10^e.donor.act[2])*(as.numeric(data[2])/abs(as.numeric(data[21])))))
  }
}

DeMMO_E_dens$DeMMO2_EDdens <- apply(reactions, 1, DeMMO2_ED_dens_fun)


DeMMO3_ED_dens_fun <- function(data){
  e.donor <- as.vector(data[9])
  e.donor<- gsub("[\\+]", "[+]", e.donor)
  e.donor <- gsub("[\\-]", "[-]", e.donor)
  e.donor <-e.donor[e.donor != ""]
  rxn.number <- reactions[7]
  e.donor.act <- as.vector(unlist(c(activities[,grepl(paste(e.donor, collapse = "|"), colnames(activities))])))
  e.donor.act.len <- length(e.donor.act[3])
  if (e.donor.act.len > 0){
    print(log10(-(10^e.donor.act[3])*(as.numeric(data[3])/abs(as.numeric(data[21])))))
  }
}

DeMMO_E_dens$DeMMO3_EDdens <- apply(reactions, 1, DeMMO3_ED_dens_fun)

DeMMO4_ED_dens_fun <- function(data){
  e.donor <- as.vector(data[9])
  e.donor<- gsub("[\\+]", "[+]", e.donor)
  e.donor <- gsub("[\\-]", "[-]", e.donor)
  e.donor <-e.donor[e.donor != ""]
  rxn.number <- reactions[7]
  e.donor.act <- as.vector(unlist(c(activities[,grepl(paste(e.donor, collapse = "|"), colnames(activities))])))
  e.donor.act.len <- length(e.donor.act[4])
  if (e.donor.act.len > 0){
    print(log10(-(10^e.donor.act[4])*(as.numeric(data[4])/abs(as.numeric(data[21])))))
  }
}

DeMMO_E_dens$DeMMO4_EDdens <- apply(reactions, 1, DeMMO4_ED_dens_fun)

DeMMO5_ED_dens_fun <- function(data){
  e.donor <- as.vector(data[9])
  e.donor<- gsub("[\\+]", "[+]", e.donor)
  e.donor <- gsub("[\\-]", "[-]", e.donor)
  e.donor <-e.donor[e.donor != ""]
  rxn.number <- reactions[7]
  e.donor.act <- as.vector(unlist(c(activities[,grepl(paste(e.donor, collapse = "|"), colnames(activities))])))
  e.donor.act.len <- length(e.donor.act[5])
  if (e.donor.act.len > 0){
    print(log10(-(10^e.donor.act[5])*(as.numeric(data[5])/abs(as.numeric(data[21])))))
  }
}

DeMMO_E_dens$DeMMO5_EDdens <- apply(reactions, 1, DeMMO5_ED_dens_fun)

DeMMO6_ED_dens_fun <- function(data){
  e.donor <- as.vector(data[9])
  e.donor<- gsub("[\\+]", "[+]", e.donor)
  e.donor <- gsub("[\\-]", "[-]", e.donor)
  e.donor <-e.donor[e.donor != ""]
  rxn.number <- reactions[7]
  e.donor.act <- as.vector(unlist(c(activities[,grepl(paste(e.donor, collapse = "|"), colnames(activities))])))
  e.donor.act.len <- length(e.donor.act[6])
  if (e.donor.act.len > 0){
    print(log10(-(10^e.donor.act[6])*(as.numeric(data[6])/abs(as.numeric(data[21])))))
  }
}

DeMMO_E_dens$DeMMO6_EDdens <- apply(reactions, 1, DeMMO6_ED_dens_fun)

DeMMO1_EA_dens_fun <- function(data){
  e.acceptor <- as.vector(data[8])
  e.acceptor<- gsub("[\\+]", "[+]", e.acceptor)
  e.acceptor <- gsub("[\\-]", "[-]", e.acceptor)
  e.acceptor <-e.acceptor[e.acceptor != ""]
  rxn.number <- data[7]
  e.acceptor.act <- as.vector(unlist(c(activities[,grepl(paste(e.acceptor, collapse = "|"), colnames(activities))])))
  e.acceptor.act.len <- length(e.acceptor.act[1])
  if (e.acceptor.act.len > 0){
    print(log10(-(10^e.acceptor.act[1])*(as.numeric(data[1])/abs(as.numeric(data[21])))))
  }
}

DeMMO_E_dens$DeMMO1_EAdens <- apply(reactions, 1, DeMMO1_EA_dens_fun)


DeMMO2_EA_dens_fun <- function(data){
  e.acceptor <- as.vector(data[8])
  e.acceptor<- gsub("[\\+]", "[+]", e.acceptor)
  e.acceptor <- gsub("[\\-]", "[-]", e.acceptor)
  e.acceptor <-e.acceptor[e.acceptor != ""]
  rxn.number <- data[7]
  e.acceptor.act <- as.vector(unlist(c(activities[,grepl(paste(e.acceptor, collapse = "|"), colnames(activities))])))
  e.acceptor.act.len <- length(e.acceptor.act[2])
  if (e.acceptor.act.len > 0){
    print(log10(-(10^e.acceptor.act[2])*(as.numeric(data[2])/abs(as.numeric(data[21])))))
  }
}

DeMMO_E_dens$DeMMO2_EAdens <- apply(reactions, 1, DeMMO2_EA_dens_fun)

DeMMO3_EA_dens_fun <- function(data){
  e.acceptor <- as.vector(data[8])
  e.acceptor<- gsub("[\\+]", "[+]", e.acceptor)
  e.acceptor <- gsub("[\\-]", "[-]", e.acceptor)
  e.acceptor <-e.acceptor[e.acceptor != ""]
  rxn.number <- data[7]
  e.acceptor.act <- as.vector(unlist(c(activities[,grepl(paste(e.acceptor, collapse = "|"), colnames(activities))])))
  e.acceptor.act.len <- length(e.acceptor.act[3])
  if (e.acceptor.act.len > 0){
    print(log10(-(10^e.acceptor.act[3])*(as.numeric(data[3])/abs(as.numeric(data[21])))))
  }
}

DeMMO_E_dens$DeMMO3_EAdens <- apply(reactions, 1, DeMMO3_EA_dens_fun)

DeMMO4_EA_dens_fun <- function(data){
  e.acceptor <- as.vector(data[8])
  e.acceptor<- gsub("[\\+]", "[+]", e.acceptor)
  e.acceptor <- gsub("[\\-]", "[-]", e.acceptor)
  e.acceptor <-e.acceptor[e.acceptor != ""]
  rxn.number <- data[7]
  e.acceptor.act <- as.vector(unlist(c(activities[,grepl(paste(e.acceptor, collapse = "|"), colnames(activities))])))
  e.acceptor.act.len <- length(e.acceptor.act[4])
  if (e.acceptor.act.len > 0){
    print(log10(-(10^e.acceptor.act[4])*(as.numeric(data[4])/abs(as.numeric(data[21])))))
  }
}

DeMMO_E_dens$DeMMO4_EAdens <- apply(reactions, 1, DeMMO4_EA_dens_fun)

DeMMO5_EA_dens_fun <- function(data){
  e.acceptor <- as.vector(data[8])
  e.acceptor<- gsub("[\\+]", "[+]", e.acceptor)
  e.acceptor <- gsub("[\\-]", "[-]", e.acceptor)
  e.acceptor <-e.acceptor[e.acceptor != ""]
  rxn.number <- data[7]
  e.acceptor.act <- as.vector(unlist(c(activities[,grepl(paste(e.acceptor, collapse = "|"), colnames(activities))])))
  e.acceptor.act.len <- length(e.acceptor.act[5])
  if (e.acceptor.act.len > 0){
    print(log10(-(10^e.acceptor.act[5])*(as.numeric(data[5])/abs(as.numeric(data[21])))))
  }
}

DeMMO_E_dens$DeMMO5_EAdens <- apply(reactions, 1, DeMMO5_EA_dens_fun)

DeMMO6_EA_dens_fun <- function(data){
  e.acceptor <- as.vector(data[8])
  e.acceptor<- gsub("[\\+]", "[+]", e.acceptor)
  e.acceptor <- gsub("[\\-]", "[-]", e.acceptor)
  e.acceptor <-e.acceptor[e.acceptor != ""]
  rxn.number <- data[7]
  e.acceptor.act <- as.vector(unlist(c(activities[,grepl(paste(e.acceptor, collapse = "|"), colnames(activities))])))
  e.acceptor.act.len <- length(e.acceptor.act[6])
  if (e.acceptor.act.len > 0){
    print(log10(-(10^e.acceptor.act[6])*(as.numeric(data[6])/abs(as.numeric(data[21])))))
  }
}

DeMMO_E_dens$DeMMO6_EAdens <- apply(reactions, 1, DeMMO6_EA_dens_fun)



#format data for plotting 
DeMMO_thermo_data <- gather(DeMMO_thermo, Site, DeltaG_norm, DeMMO1_DeltaG_norm:DeMMO6_DeltaG_norm, factor_key=TRUE)
DeMMO_thermo_data <- DeMMO_thermo_data[20:21]
DeMMO_thermo_data$Site <- str_replace(DeMMO_thermo_data$Site, "_DeltaG_norm", " ")
DeMMO_thermo_data$rxn.number <- c(1:33)
DeMMO_thermo_data$mineral <- reactions$reactant.a

DeMMO_thermo_data <- DeMMO_thermo_data %>%
  filter(Site=="DeMMO1 " | Site=="DeMMO3 " | Site=="DeMMO6 ")

deltaG_plot <- ggplot(DeMMO_thermo_data, aes(DeltaG_norm, reorder(rxn.number, -DeltaG_norm), shape=Site, group=rxn.number)) +
  theme_gray() +
  geom_line(aes(color=mineral), size=2.5, alpha=0.6) +
  geom_point() + 
  scale_shape_manual(values = c(15,16,17)) + 
  scale_x_reverse() +
  xlab("ΔGr (kJ/mol e-)") + 
  ylab("Reaction #") +
  geom_vline(xintercept = 0, linetype="dotted", color = "black") +
  theme(legend.position = c(.1, .84), legend.text=element_text(size=6), legend.title = element_text(size=8, face="bold")) +
  theme(legend.key.size =  unit(0.1, "in"))

#DeMMO_E_dens <- DeMMO_E_dens %>% select(-DeMMO2_EDdens, -DeMMO4_EDdens, -DeMMO5_EDdens, -DeMMO2_EAdens, -DeMMO4_EAdens, -DeMMO5_EAdens)

Edens_plot_data <- gather(DeMMO_E_dens, Site, E_dens, DeMMO1_EDdens:DeMMO6_EAdens, factor_key=TRUE)
Edens_plot_data$Site <- str_replace(Edens_plot_data$Site, "_EDdens", " ")
Edens_plot_data$Site <- str_replace(Edens_plot_data$Site, "_EAdens", " ")
Edens_plot_data$Type <- c(rep("ED", 198), rep("EA", 198))
E.donors <- as.vector(reactions$e.donor)
E.acceptors <- as.vector(reactions$e.acceptor)
Edens_plot_data$E.donor.acceptor <- c(rep(E.donors, 6), rep(E.acceptors, 6))
Edens_plot_data$DeltaG_norm <- DeMMO_thermo_data$DeltaG_norm
Edens_plot_data$E.donor.acceptor <- str_replace(Edens_plot_data$E.donor.acceptor, c("sulfur|pyrite|siderite|hematite|gypsum|magnetite|pyrolusite|Mn[+]2"), " ")

#convert E_dens column from list to numeric 
Edens_plot_data$E_dens <- as.numeric(as.character(Edens_plot_data$E_dens))
#Edens_plot_data <- na.omit(Edens_plot_data)

Edens_plot_data$E.donor.acceptor <- factor(Edens_plot_data$E.donor.acceptor,  levels = c("SO4-2", "NO3-", "CH4", "NH4+", "HS-", "HCO3-","H2", "Fe+2", "acetate"))

Edens_plot_data <- Edens_plot_data %>%
  filter(Site=="DeMMO1 " | Site=="DeMMO3 " | Site=="DeMMO6 ")

Edens_plot <- ggplot(Edens_plot_data, aes(E_dens, reorder(rxn.number, -DeltaG_norm), shape=Site, group=rxn.number)) +
  theme_gray() +
  geom_line(aes(color=E.donor.acceptor), size=2.5, alpha=0.6) +
  geom_point(show.legend = FALSE) + 
  scale_shape_manual(values = c(15,16,17)) + 
  coord_cartesian(xlim = c(-5, 2)) +
  xlab("ΔGr (log J/kg H2O)") + 
  geom_hline(yintercept = c(17,18,21,31), color = "black", linetype="dotted", size=1.5, alpha=0.3) +
  theme(panel.grid.minor = element_line(colour="white", size=0.5)) +
  theme(axis.title.y = element_blank()) +
  scale_y_discrete(breaks = seq(0, 33, 1)) +
  theme(legend.position = c(.15, .89), legend.text=element_text(size=6), legend.title = element_text(size=8, face="bold")) +
  theme(legend.key.size =  unit(0.005, "in"))


grid.newpage()
grid.draw(cbind(ggplotGrob(deltaG_plot), ggplotGrob(Edens_plot), size = "last"))
