library(ggplot2)
library(tidyr)
library(stringr)


#read in csv files (these were generated with separate code, I removed the other code to simplify here)
reactions <- read.csv("reactions.csv", header = TRUE)
DeMMO_thermo <- read.csv("DeMMO_thermo.csv", header=TRUE)


#format data for plotting (this is to give an example of how to reshape data from wide to long format for plotting)
DeMMO_thermo_data <- gather(DeMMO_thermo, Site, DeltaG_norm, DeMMO1_DeltaG_norm:DeMMO6_DeltaG_norm, factor_key=TRUE)
DeMMO_thermo_data <- DeMMO_thermo_data[21:22]
DeMMO_thermo_data$Site <- str_replace(DeMMO_thermo_data$Site, "_DeltaG_norm", " ")
DeMMO_thermo_data$rxn.number <- c(1:33)
DeMMO_thermo_data$mineral <- reactions$reactant.a

#store plot in variable deltaG_plot, to execute type deltaG_plot in the console 
##reorder function orders the DeltaG values so that the most exergonic are on top, alpha value in geom_line sets transparency of lines 
deltaG_plot <- ggplot(DeMMO_thermo_data, aes(DeltaG_norm, reorder(rxn.number, -DeltaG_norm), shape=Site, group=rxn.number)) +
  geom_line(aes(color=mineral), size=2.5, alpha=0.6) +
  geom_point() + 
  scale_shape_manual(values = c(17, 25,18, 19, 1, 15)) + 
  scale_x_reverse() +
  xlab("ΔGr (kJ/mol e-)") + 
  ylab("Reaction #") +
  geom_vline(xintercept = 0, linetype="dotted", color = "black")