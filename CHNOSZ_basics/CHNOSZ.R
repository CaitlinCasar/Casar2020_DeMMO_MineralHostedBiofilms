install.packages("CHNOSZ")
library(CHNOSZ)
library(data.table)

data(thermo)
db<- data.frame(thermo$obigt)
write.csv(db, file = "CHNOSZ.csv")

#add a mineral to the database
icaitlin <- mod.obigt("caitlin", formula = "MnO22",H=77, S=22, G=20970)
info(icaitlin)

#set units, default units of temperature, pressure, and energy are °C, bar, and calories
T.units("K")
P.units("MPa")
E.units("J")

#calculate Gibbs energy in J/mol for aqueous methane at 298.15 K and 0.1 MPa
subcrt("methane", T = 298.15, P = 0.1)$out$methane$G

#can also use convert function to display Gibbs energy in J, this modifies db standard energy unit
convert(info(info("methane"))$G, "J")

#restore database defaults 
data(thermo)

#use subcrt function to calculate standard molal properties of adenine at 100C
T.units("C")
subcrt("adenine", T = 100)

#define basis species with basis function to automatically balance reactions
basis(c("H+", "H2", "H2O"), c("aq","gas","liq"))


#store balanced reaction with magnetite in variable 'magnetite_reduction'
magnetite_reduction <- subcrt(c("magnetite", "H2", "H+","Fe+2","H2O"), c(-1, -1, -6, 3, 4), c("cr","gas","aq", "aq", "liq"),T=c(294.65)$out$logK

gypsum_reduction <- subcrt(c("gypsum", "H2","H+", "HS-", "Ca+2", "H2O"), c(-1, -4, -1, 1,1, 6))

siderite_oxidation <- subcrt(c("siderite","NO3-","H2O", "goethite", "NO2-", "HCO3-","H+"), c(-2, -1, -3,2,1,2,2), c("cr","aq","liq","cr","aq","aq","aq"), T=294.65)$out$logK
                             

gypsum_acetate <- subcrt(c("gypsum", "H+", "acetate", "HCO3-", "HS-", "Ca+2", "H2O"), c(-9,-9,-8,16,9,9,22))

#display the balanced reaction for magnetite reduction
#this displays the wrong reaction for some reason...
plot(0, 0, type = "n", axes = FALSE, ann=FALSE, xlim=c(0, 5), ylim=c(5.2, -0.2))
text(0, 0, "magnetite reduction", adj = 0)
text(5, 1, describe.reaction(magnetite_reduction$reaction), adj = 1)
text(0, 1, "siderite oxidation", adj = 0)
text(5, 2, describe.reaction(siderite_oxidation$reaction), adj = 1)
text(0, 2, "gypsum reduction", adj = 0)
text(5, 3, describe.reaction(gypsum_reduction$reaction), adj = 1)


#define the basis species for auto balancing reactions with acetate and methane
basis(c("CO2", "H2", "H2O", "H+"))

#store balanced reactions for acetoclastic methanogenesis and acetate oxidation
acetoclastic <- subcrt(c("acetate", "methane"), c(-1, 1))
acetate_oxidation <- subcrt("acetate", -1)

#display the balanced reactions for acetoclastic methanogenesis and acetate oxidation
plot(0, 0, type = "n", axes = FALSE, ann=FALSE, xlim=c(0, 5), ylim=c(5.2, -0.2))
text(0, 0, "acetoclastic methanogenesis", adj = 0)
text(5, 1, describe.reaction(acetoclastic$reaction), adj = 1)
text(0, 2, "acetate oxidation", adj = 0)
text(5, 3, describe.reaction(acetate_oxidation$reaction), adj = 1)


#activities of basis species can be modified with basis(), and those of the other species using the logact argument in subcrt()
#calculate logarithms of activity coefficients (loggam) for three aqueous species at two ionic strengths and temperatures:
subcrt(c("Fe+2", "Na+", "K+"),
       T = c(25, 100), IS = c(0, 0.25), property = "G")$out


#calcuate ionic strength I = 1/2 Σmizi^2 where m is concentration in moles/L and z is charge
DeMMO_chemistry <- read.csv("DeMMO_Chemistry_Dec15toSep18.csv")


#add pyrolusite to database
pyrolusite <- mod.obigt("pyrolusite", G=-111100, H=-124283, S=12.61, V=17.3, formula="MnO2", state="cr", a1.a=12.55, a2.b=9.761, a3.c=-2.105)

#add ferrihydrite to database
ferrihydrite <- mod.obigt("ferrihydrite", G=-111200, H=-127800, S=16.7, V=20.88, formula="FeOOH", state="cr", a1.a=8.70, a2.b=36.71, a3.c=-1.0146)

#add manganite to database 
manganite <- mod.obigt("manganite", G=-133300, formula="MnOOH", state="cr")


#calculate logK for magnetite reduction at DeMMO 1-6 (change T for each site)
magnetite_reduction_logK <- subcrt(c("magnetite", "H2", "H+","Fe+2","H2O"), c(-1, -1, -6, 3, 4), c("cr","gas","aq", "aq", "liq"),T=294.65)$out$logK



#testing output for logK to compare with Doug's SUPCRT output
methane_ox_logK <- subcrt(c("methane", "oxygen", "carbon dioxide","H2O"), c(-1, -2, 1, 2), c("gas", "gas", "gas", "liq"), T=408)

