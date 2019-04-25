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
library(dplyr)
library(grid)


#add pyrolusite to database
pyrolusite <- mod.obigt("pyrolusite", G=-111100, H=-124283, S=12.61, V=17.3, formula="MnO2", state="cr", a1.a=12.55, a2.b=9.761, a3.c=-2.105)

#add ferrihydrite to database
ferrihydrite <- mod.obigt("ferrihydrite", G=-111200, H=-127800, S=16.7, V=20.88, formula="FeOOH", state="cr", a1.a=8.70, a2.b=36.71, a3.c=-1.0146)

#add manganite to database 
manganite <- mod.obigt("manganite", G=-133300, formula="MnOOH", state="cr")

#store database in a dataframe for reference 
thermo_db <- data.frame(thermo()$obigt)

#set temperature units to Celcius
T.units("K")
E.units("J")

#import DeMMO mineral reactions
reactions <- read.csv("reactions_aq_gas.csv", header = TRUE)

#import DeMMO fluid activities 
activities <- read.csv("DeMMO_SpecE8_aqueousGas.csv", header=TRUE)
colnames(activities) <- c("Site","Ca+2", "acetate","methane","Fe+2",	"H+",	"H2",	"HCO3-","HS-","Mn+2","NH4+","NO2-","NO3-","SO4-2", "CO")

#initialize logK dataframe
DeMMO_logK <- data.frame()   

#store average fluid temperatures for DeMMO1-DeMMO6
DeMMO_T <- c(283.45, 285.55, 289.35, 295.65, 304.85, 294.65)

#create function to calulate logK values for every reaction for six DeMMO sites 
DeMMO_logK_fun <- function(data){
  rxn.number <- data[1]
  react.prod <- as.vector(unlist(c(data[4:12])))
  react.prod <-react.prod[react.prod != ""] 
  react.coeffs <- as.vector(na.omit(as.numeric(data[14:22])))
  react.states <- as.vector(unlist(c(data[23:31])))
  react.states <- react.states[react.states != ""]
  for(temperature in DeMMO_T){
  if(data[9] != "manganite"){
    rxn <- subcrt(react.prod, react.coeffs, react.states, T=temperature)
    DeMMO_logK[rxn.number, match(temperature,DeMMO_T)] <<- rxn$out$logK
  }else{
    #calculate logK at 25C for manganite rxns
    rxn <- subcrt(react.prod, react.coeffs, react.states, T=298.15)
    DeMMO_logK[rxn.number, match(temperature,DeMMO_T)] <<- rxn$out$logK
    colnames(DeMMO_logK)[match(temperature,DeMMO_T)] <<- paste0("DeMMO", match(temperature,DeMMO_T), ".logK")
  }
  }
}

#store the output of the function applied to the reactions dataframe in new dataframe DeMMO_thermo
apply(reactions, 1, DeMMO_logK_fun)

#initialize 
DeMMO_logQ <- data.frame()

DeMMO_logQ_fun <- function(data){
  for(site.number in 1:6){
  rxn.reactants <- as.vector(unlist(c(data[4:7])))
  rxn.reactants<- gsub("[\\+]", "[+]", rxn.reactants)
  rxn.reactants <- gsub("[\\-]", "[-]", rxn.reactants)
  rxn.reactants <-rxn.reactants[rxn.reactants != ""]
  rxn.react.act <- activities[,grepl(paste(rxn.reactants, collapse = "|"), colnames(activities))]
  rxn.react.act.len <- length(rxn.react.act)
  if (rxn.react.act.len > 1) {
    rxn.react.act.index <- match(grep(paste(rxn.reactants, collapse = "\\b|"), names(activities), value=TRUE), as.vector(unlist(c(data[4:7]))))
    rxn.react.coeff <- abs(as.vector(na.omit(as.numeric(data[14:17]))))[rxn.react.act.index]
    react.vec.len <- length(as.vector(grep(paste(rxn.reactants, collapse = "|"), names(activities), value=TRUE)))
    if (react.vec.len > 1) {
      log.act.reactants <- rxn.react.act[site.number,]*rxn.react.coeff
    } else{
      log.act.reactants <- rxn.react.act[site.number]*rxn.react.coeff
    }
  }else{
    log.act.reactants <- 1
  }
    
  rxn.products <- as.vector(unlist(c(data[8:12])))
  rxn.products<- gsub("[\\+]", "[+]", rxn.products)
  rxn.products <- gsub("[\\-]", "[-]", rxn.products)
  rxn.products <-rxn.products[rxn.products != ""] 
  rxn.prod.act <- activities[,grepl(paste(rxn.products, collapse = "|"), colnames(activities))]
  rxn.prod.act.index <- match(grep(paste(rxn.products, collapse = "\\b|"), names(activities), value=TRUE), as.vector(unlist(c(data[8:12]))))
  rxn.prod.coeff <- abs(as.vector(na.omit(as.numeric(data[18:22]))))[rxn.prod.act.index]
  prod.vec.len <- length(as.vector(grep(paste(rxn.products, collapse = "|"), names(activities), value=TRUE)))
  if (prod.vec.len > 1) {
    log.act.products <- rxn.prod.act[site.number,]*rxn.prod.coeff
  } else{
    log.act.products <- rxn.prod.act[site.number]*rxn.prod.coeff
  }
  DeMMO_logQ[data[1], site.number] <<- sum(log.act.products) + sum(-(log.act.reactants))
  colnames(DeMMO_logQ)[site.number] <<-paste0("DeMMO", site.number, ".logQ")
  }
  }

apply(reactions, 1, DeMMO_logQ_fun)

#combine logK and logQ values into one dataframe to make it easier to iterate over rows to calculate deltaG
DeMMO_thermo <- cbind(DeMMO_logK, DeMMO_logQ)
DeMMO_thermo$e.transfer <- reactions$e.transfer
DeMMO_thermo$rxn.number <- reactions$rxn.number


#initialize DeltaG dataframe
DeMMO_DeltaG <- data.frame()

DeMMO_DeltaG_fun <- function(data){
  for(site.number in 1:6){
    deltaGrxn <- -2.303*8.314*283.45*(as.numeric(as.vector(data[site.number]))-as.numeric(as.vector(data[site.number+6])))
    DeMMO_DeltaG[data[14],site.number] <<- deltaGrxn
    colnames(DeMMO_DeltaG)[site.number] <<- paste0("DeMMO", site.number, ".DeltaG")
  }
  for(site.number in 1:6){
  DeMMO_DeltaG[data[14],site.number+6] <<- DeMMO_DeltaG[data[14],site.number]/(as.numeric(as.vector(data[13]))*1000)
  colnames(DeMMO_DeltaG)[site.number+6] <<- paste0("DeMMO", site.number, ".DeltaG.norm")
  }
}

#calculate deltaG for each reaction
apply(DeMMO_thermo, 1, DeMMO_DeltaG_fun)

#add all data to DeMMO_thermo dataframe
DeMMO_thermo <-cbind(reactions, DeMMO_DeltaG, DeMMO_thermo[,-c(13:14)])

###calculate energy density by normalizing to e-donor (ED) or e-acceptor (EA)
DeMMO_E_dens <- data.frame()
DeMMO_ED_dens_fun <- function(data){
  for(site.number in 1:6){
  e.donor <- as.vector(data[3])
  e.donor<- gsub("[\\+]", "[+]", e.donor)
  e.donor <- gsub("[\\-]", "[-]", e.donor)
  e.donor <-e.donor[e.donor != ""]
  e.donor.act <- as.vector(unlist(c(activities[,grepl(paste(e.donor, collapse = "\\b|"), colnames(activities))])))
  e.donor.act.len <- length(e.donor.act[site.number])
  if (e.donor.act.len > 0){
    DeMMO_E_dens[data[1], site.number] <<- log10(-(10^e.donor.act[site.number])*(as.numeric(data[ncol(reactions)+site.number])/abs(as.numeric(data[13]))))
  }else{
    DeMMO_E_dens[data[1], site.number] <<- "NA"
  }
  colnames(DeMMO_E_dens)[site.number] <<- paste0("DeMMO", site.number, ".ED_dens")
  }
  for(site.number in 1:6){
  e.acceptor <- as.vector(data[2])
  e.acceptor<- gsub("[\\+]", "[+]", e.acceptor)
  e.acceptor <- gsub("[\\-]", "[-]", e.acceptor)
  e.acceptor <-e.acceptor[e.acceptor != ""]
  e.acceptor.act <- as.vector(unlist(c(activities[,grepl(paste(e.acceptor, collapse = "|"), colnames(activities))])))
  e.acceptor.act.len <- length(e.acceptor.act[site.number])
  if (e.acceptor.act.len > 0){
    DeMMO_E_dens[data[1], site.number+6] <<- log10(-(10^e.acceptor.act[site.number])*(as.numeric(data[ncol(reactions)+site.number])/abs(as.numeric(data[13]))))
  }else{
    DeMMO_E_dens[data[1], site.number+6] <<-NA
  }
  colnames(DeMMO_E_dens)[site.number+6] <<- paste0("DeMMO", site.number, ".EA_dens")
}
}

#calculate energy density for electron donors and acceptors at all sites 
apply(DeMMO_thermo, 1, DeMMO_ED_dens_fun)



#format data for plotting 
DeMMO_thermo_data <- gather(DeMMO_thermo, Site, DeltaG_norm, DeMMO1.DeltaG.norm:DeMMO6.DeltaG.norm, factor_key=TRUE)
DeMMO_thermo_data <- DeMMO_thermo_data[50:51]
DeMMO_thermo_data$Site <- str_replace(DeMMO_thermo_data$Site, ".DeltaG.norm", "")
DeMMO_thermo_data$rxn.number <- c(1:nrow(reactions))
DeMMO_thermo_data$mineral <- reactions$reactant.a

DeMMO_thermo_data <- DeMMO_thermo_data %>%
  filter(Site=="DeMMO1" | Site=="DeMMO3" | Site=="DeMMO6")

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

#gather Edens data into long dataframe for plotting
DeMMO_Edens_data <- gather(DeMMO_E_dens, Site, E_dens, DeMMO1.ED_dens:DeMMO6.EA_dens, factor_key=TRUE, na.rm = FALSE)
#format Site names
DeMMO_Edens_data$Site <- str_replace(DeMMO_Edens_data$Site, ".ED_dens", "")
DeMMO_Edens_data$Site <- str_replace(DeMMO_Edens_data$Site, ".EA_dens", "")
#filter data for only DeMMO1,3,and 6 for this study
DeMMO_Edens_data <- DeMMO_Edens_data %>% 
  filter(Site=="DeMMO1"|Site=="DeMMO3"|Site=="DeMMO6")
DeMMO_Edens_data$Type <- c(rep("ED", nrow(DeMMO_Edens_data)/2), rep("EA", nrow(DeMMO_Edens_data)/2))

#create lists of e- donors and acceptors 
E.donors <- as.vector(reactions$e.donor)
E.acceptors <- as.vector(reactions$e.acceptor)

#add e- donors and acceptors to Edens plotting dataframe 
DeMMO_Edens_data$E.donor.acceptor <- c(rep(E.donors, 3), rep(E.acceptors, 3))

#add deltaGnorm values to Edens plotting dataframe 
DeMMO_Edens_data$DeltaG_norm <- DeMMO_thermo_data$DeltaG_norm

#add reaction numbers to Edens data for plotting 
DeMMO_Edens_data <- cbind(rxn.number = DeMMO_thermo$rxn.number, DeMMO_Edens_data)

sulfur.rxns <- as.numeric(unlist(unique(DeMMO_Edens_data %>% 
  filter(E.donor.acceptor == "sulfur") %>% 
  select(rxn.number))))

ordered.data <- levels(reorder(DeMMO_Edens_data$rxn.number, -DeMMO_Edens_data$DeltaG_norm))

#convert E_dens column from list to numeric 
DeMMO_Edens_data$E_dens <- as.numeric(as.character(DeMMO_Edens_data$E_dens))

#set aqueous e- donors/acceptors to factors for plot legend
DeMMO_Edens_data$E.donor.acceptor <- factor(DeMMO_Edens_data$E.donor.acceptor,  levels = c("SO4-2", "NO3-", "methane", "NH4+", "HS-", "HCO3-","H2", "Fe+2", "acetate", "Mn+2", "CO"))

Edens_plot <- ggplot(DeMMO_Edens_data, aes(E_dens, reorder(rxn.number, -DeltaG_norm), shape=Site, group=rxn.number)) +
  theme_gray() +
  geom_line(aes(color=E.donor.acceptor), size=2.5, alpha=0.6) +
  geom_point(show.legend = FALSE) + 
  scale_shape_manual(values = c(15,16,17)) + 
  coord_cartesian(xlim = c(-6.5, 3)) +
  xlab("ΔGr (log J/kg H2O)") + 
  geom_hline(yintercept = match(sulfur.rxns, ordered.data), color = "black", linetype="dotted", size=1.5, alpha=0.3) +
  theme(panel.grid.minor = element_line(colour="white", size=0.5)) +
  theme(axis.title.y = element_blank()) +
  scale_y_discrete(breaks = seq(0, nrow(reactions), 1)) +
  theme(legend.position = c(.15, .2), legend.text=element_text(size=6), legend.title = element_text(size=8, face="bold")) +
  theme(legend.key.size =  unit(0.005, "in"))


grid.newpage()
grid.draw(cbind(ggplotGrob(deltaG_plot), ggplotGrob(Edens_plot), size = "last"))





