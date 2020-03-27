#load dependencies 
pacman::p_load(tidyverse, readr)

geochem_data <- read_csv("../orig_data/site_geochem.csv")

avg_geochem_data <- geochem_data %>% 
  mutate(sulfide = sulfide*.001) %>% #convert sulfide from ug/L to mg/L
  gather(solute, value, methane:DOC) %>%
  group_by(Site, solute) %>%
  summarise(mean = mean(value, na.rm = T),
            st_dev = sd(value, na.rm = T)) %>%
  ungroup() %>%
  mutate(solute = factor(solute, levels=c("hydrogen", "CO", "methane", "sulfide", "ammonia", "nitrate", "DOC", "ferrous_iron", "sulfate")))
  

geochem_plot <- ggplot(avg_geochem_data, aes(solute, mean, color=solute)) +
  geom_point(size=5) + 
  scale_y_log10() + 
  #ggplot2::geom_errorbar(ggplot2::aes(ymin=mean-std_dev, ymax=mean+std_dev), width=0.2) + 
  ggplot2::coord_flip() +
  ggplot2::geom_line(ggplot2::aes(group=solute), color="gray", linetype="dotted") +
  facet_wrap(~Site) +
  theme(legend.position = "none", 
        axis.title.y=ggplot2::element_blank()) +
  ylab("Concentration (mg/L)") +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) 
