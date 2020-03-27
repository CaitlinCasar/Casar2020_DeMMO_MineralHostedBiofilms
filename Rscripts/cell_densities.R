#load dependencies 
pacman::p_load(tidyverse, readr, plyr, plotly, lubridate, gridExtra, scales)

cell_densities <- read_csv("../orig_data/cellCounts_sqcm.csv")

mineral_color_dict <- c("#f8766d", "#d28f00", "#93aa00", "#00ba38", "#00c19e", "#00b9e3", "#5e9aff", "#db6efb", "#fd70c8")
names(mineral_color_dict) <- c("fluid", "calcite", "muscovite", "siderite", "pyrite", "hematite", "magnetite", "pyrolusite", "sand")

shape_dict <- c(16, 8, 15)
names(shape_dict) <- c("fluid", "control", "mineral")


cell_density_plot <- cell_densities %>% dplyr::filter(!Cell.density.sq.cm %in% c(NA, 0)) %>% 
  #dplyr::filter(!Experiment == "control" && !Mineral == "sand") %>%
  dplyr::group_by(Site, Experiment, Mineral) %>% dplyr::summarise(density=mean(Cell.density.sq.cm)) %>%
  dplyr::mutate(Mineral = factor(Mineral, levels = names(mineral_color_dict))) %>% 
  ggplot2::ggplot(aes(density, Mineral)) +
  geom_point(aes(color=Mineral, shape=Experiment), size=5) + 
  scale_shape_manual(values=shape_dict) +
  labs(x=expression('Cell density (cells/'~cm^2~')')) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme(axis.title.y=element_blank(),
        legend.position = "none")

cell.density.plot +  facet_grid(cols = vars(Site))
