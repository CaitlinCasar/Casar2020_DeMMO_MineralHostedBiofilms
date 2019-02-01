library(ggplot2)
library(scales) 
library(cowplot)
library(gridExtra)
install.packages("extrafont")
data <- read.csv("SiteGeochem.csv")

temp <- data[1:3,]
redox<- data[7:21,]
orp<-data[4:6,]
doc_dic<- data[22:27,]
cell_dens<-data[28:30,]
gas<-data[31:36,]

redox_data <- read.csv("redox_geochem.csv")




temp_plot<- ggplot(temp, aes(x=measurement, y=depth_position, color=variable)) + 
  geom_path(linetype = "dotted", color="black") +
  geom_point(size=5) + 
  scale_x_continuous(name="Temp oC", limits=c(10, 25)) + 
  scale_y_reverse() +
  theme(axis.title.y = element_blank()) +
  theme(legend.position="none") +
  ggtitle("Temperature")

orp_plot<- ggplot(orp, aes(x=measurement, y=depth_position, color=variable)) + 
  geom_path(linetype = "dotted", color="black") +
  geom_point(size=5) + 
  scale_x_continuous(name="ORP", limits=c(-250, -25)) +
  scale_y_reverse() +
  theme(axis.title.y = element_blank()) +
  theme(legend.position="none") +
  ggtitle("Redox Potential")

redox_plot<-ggplot(redox, aes(x=measurement, y=depth_position, group=variable, color=variable)) + 
  geom_path(linetype = "dotted", color="black") +
  geom_point(size=5) +
  scale_x_log10(name="mg/L | S2- ug/L",breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_reverse() +
  theme(axis.title.y = element_blank()) +
  theme(legend.position="none") +
  ggtitle("Redox Species")

doc_dic_plot<-ggplot(doc_dic, aes(x=measurement, y=depth_position, group=variable, color=variable)) + 
  geom_path(linetype = "dotted", color="black") +
  geom_point(size=5) +
  scale_x_log10(name="mg/L | mM",breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_reverse() +
  theme(axis.title.y = element_blank()) +
  theme(legend.position="none") +
  ggtitle("DOC | DIC")


cell_dens_plot<- ggplot(cell_dens, aes(x=measurement, y=depth_position, color=variable)) + 
  geom_path(linetype = "dotted", color="black") +
  geom_point(size=5) + 
  scale_x_log10(name="Cells/mL", breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_reverse() +
  theme(axis.title.y = element_blank()) +
  theme(legend.position="none") +
  ggtitle("Cell Density")

gases_plot<- ggplot(gas, aes(x=measurement, y=depth_position, color=variable)) + 
  geom_path(linetype = "dotted", color="black") +
  geom_point(size=5) + 
  scale_x_log10(name="nM", breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_reverse() +
  theme(axis.title.y = element_blank()) +
  theme(legend.position="none") +
  ggtitle("H2 | CH4")

redox_plot_facet<-ggplot(redox_data, aes(x=measurement, y=depth_position, group=variable, color=variable)) + 
  geom_path(linetype = "dotted", color="black") +
  geom_point(size=5) +
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x))
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_reverse() +
  theme(axis.title.y = element_blank()) +
  ggtitle("Redox Species")
redox_plot_facet
redox_plot_facet + facet_wrap( ~ facet_order, ncol = 8)


dissolved_species<-redox_data[1:18,]
gasesous_species<-redox_data[19:24,]

dissolved_plot<-ggplot(dissolved_species, aes(x=measurement, y=depth_position, group=variable, color=variable)) + 
  geom_path(linetype = "dotted", color="black") +
  geom_point(size=5) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_reverse() +
  theme(axis.title.y = element_blank()) +
  ggtitle("Dissolved Species") +
  theme(legend.position="none") 

gaseous_plot<-ggplot(gasesous_species, aes(x=measurement, y=depth_position, group=variable, color=variable)) + 
  geom_path(linetype = "dotted", color="black") +
  geom_point(size=5) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_reverse() +
  theme(axis.title.y = element_blank()) +
  ggtitle("Gaseous Species") +
  theme(legend.position="none") 

plot_grid(gaseous_plot, dissolved_plot, labels = "AUTO", align = 'h')

grid.arrange(
  temp_plot,
  orp_plot,
  redox_plot,
  doc_dic_plot,
  cell_dens_plot,
  gases_plot,
  nrow = 1)

grid.arrange(
  orp_plot,
  redox_plot,
  doc_dic_plot,
  gases_plot,
  nrow = 1)

#p + facet_wrap( ~ facet_order, ncol = 6, scales = "free")

#p + facet_grid(. ~ variable)


             