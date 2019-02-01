install.packages("CHNOSZ")
#install developer update:
#> install.packages("/Users/Caitlin/Downloads/CHNOSZ_1.1.3-63.tgz", repos = NULL)
library(CHNOSZ)
library(data.table)

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
T.units("C")

#store array of DeMMO 1-6 fluid temperatures in Celcius in DeMOMO_T
DeMMO_T <- c(10.3, 12.4, 16.2, 22.5, 31.7, 21.5)

#pyrolusite reactions
pyrolusite_H2 <- subcrt(c("pyrolusite", "H2", "H+", "Mn+2", "H2O"), c(-1, -1, -2, 1, 2), c("cr", "gas","aq","aq","liq"), T=DeMMO_T)

#magnetite reactions 
magnetite_H2 <- subcrt(c("magnetite", "H2", "H+","Fe+2","H2O"), c(-1, -1, -6, 3, 4), c("cr","gas","aq", "aq", "liq"), T=DeMMO_T)

magnetite_HS <- subcrt(c("magnetite", "HS-", "H+", "Fe+2", "sulfur", "H2O"), c(-1, -1, -7, 3, 1, 4), c("cr","aq", "aq", "aq", "cr", "liq"), T=DeMMO_T)

#check that rxn coefficients are set properly by running magnetite_H2$reaction$coeff
magnetite_H2$reaction$coeff
###output shows coefficients do not match my manual input above:
###> magnetite_H2$reaction$coeff
###[1] -1 -1 -1 -6  3

#gypsum reactions 
gypsum_H2 <- subcrt(c("gypsum","H2", "H+", "HS-", "Ca+2", "H2O"), c(-1, -4, -1, 1, 1, 6), c("cr", "gas","aq", "aq", "aq", "liq"), T=DeMMO_T)

###output:
###Error in subcrt(c("gypsum", "H2", "H+", "HS-", "Ca+2", "H2O"), c(-1, -4,  : 
                                                                   ###-9.3 are not valid properties
                                                                 ###try rho, logK, G, H, S, V, Cp, kT, or E-11.4 are not valid properties
                                                                 ###try rho, logK, G, H, S, V, Cp, kT, or E-15.2 are not valid properties
                                                                 ###try rho, logK, G, H, S, V, Cp, kT, or E-21.5 are not valid properties
                                                                 ###try rho, logK, G, H, S, V, Cp, kT, or E-30.7 are not valid properties
                                                                 ###try rho, logK, G, H, S, V, Cp, kT, or E-20.5 are not valid properties
                                                                 ###try rho, logK, G, H, S, V, Cp, kT, or E


#pyrite reactions

#siderite reactions 
siderite_NO3 <- subcrt(c("siderite","NO3-","H2O", "goethite", "NO2-", "HCO3-","H+"), c(-2, -1, -3,2,1,2,2), c("cr","aq","liq","cr","aq","aq","aq"), T=DeMMO_T)

#hematite reactions 
hematite_H2 <- subcrt(c("hematite", "H2", "H+", "Fe+2", "H2O"), c(-1, -1, -4, 2, 3), c("cr", "gas", "aq", "aq", "liq"), T=DeMMO_T)

#check that rxn coefficients are set properly by running hematite_H2$reaction$coeff
hematite_H2$reaction$coeff
###output shows coefficients do not match my manual input above:
###> hematite_H2$reaction$coeff
###[1] -1 -1 -1 -1 -4

#output dataframe
DeMMO_thermo <- data.frame(Site=c("DeMMO1", "DeMMO2", "DeMMO3", "DeMMO4", "DeMMO5", "DeMMO6"))
DeMMO_thermo$pyrolusite_H2_logK <- pyrolusite_H2$out$logK
DeMMO_thermo$magnetite_H2_logK <- magnetite_H2$out$logK
DeMMO_thermo$magnetite_HS_logK <- magnetite_HS$out$logK
DeMMO_thermo$siderite_NO3_logK <- siderite_NO3$out$logK
DeMMO_thermo$hematite_H2_logK <- hematite_H2$out$logK


##with developer update, logK output is NA for rxns with magnetite and hematite 
