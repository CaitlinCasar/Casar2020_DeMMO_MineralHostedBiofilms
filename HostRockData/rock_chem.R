library(tidyverse)
library(ggplot2)

data <- read.csv("Demmo_Rocks.csv", header=TRUE)

data <- gather(data, Element, ppm, Ti:Ca, factor_key = TRUE)
data <- data %>% filter(Element %in% c("S", "Fe", "Mn"))
data$Rock <- gsub("Demmo3", "Poorman", data$Rock) 
data$Rock <- gsub("Demmo1", "Homestake", data$Rock) 
data$Rock <- gsub("Demmo2", "Ellison", data$Rock)


rock_chem.plot <- ggplot(data, aes(Element, ppm, color=Rock)) +
  geom_point(size=3) +
  scale_y_continuous(trans = "log10") + 
  theme_gray()
  #labs(x=expression('Cell density ('~10^4~'cells/'~mm^2~')')) +
  #scale_x_continuous(
    #breaks = c(0, 20000, 40000, 60000),
    #label = c(0, 2, 4, 6)) +
  #theme(axis.title.y="element_blank()") +
#legend.position = "none") 

rock_chem.plot +  facet_grid(cols = vars(Element))