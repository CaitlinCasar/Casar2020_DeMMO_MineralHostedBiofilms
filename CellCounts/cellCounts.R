library(ggplot2)
library(plotly)
library(lubridate)
library(scales)


#import cell count data
cell.density.data <- read.csv("cellCounts_sqcm.csv", header=TRUE)
cell.density.data$Site <- gsub('DeMMO ', 'D', cell.density.data$Site)
cell.density.data$Site <- gsub('Demmo ', 'D', cell.density.data$Site)
cell.density.data$Mineral <- tolower(cell.density.data$Mineral)
cell.density.data$Mineral <- gsub('control sand', 'sand', cell.density.data$Mineral)
cell.density.data$Exp.Type <- paste(cell.density.data$Mineral, cell.density.data$Experiment, sep='.')
cell.density.data$Exp.Type <- tolower(cell.density.data$Exp.Type)
cell.density.data$Date <-format(as.Date(cell.density.data$Date.Deployed), "%Y-%m")
cell.density.data$Site.Date <- paste(cell.density.data$Site, cell.density.data$Date, sep='.')

cell.density.data <- cell.density.data[order(cell.density.data$Site.Date),] 
#cell.density.data$Site.Date <- factor(cell.density.data$Site.Date, levels=c("D3.2017-02", "D6.2017-02", "D1.2017-08","D3.2017-08", "D1.2017-11", "D3.2017-11"))

#reorder experiments for plot 
exp.order <- rev(c("sand.control","pyrolusite.mineral","pyrolusite.control","magnetite.mineral","magnetite.control",
                   "hematite.mineral","hematite.control","pyrite.mineral","pyrite.control",
                   "siderite.mineral","siderite.control","muscovite.mineral","muscovite.control","calcite.mineral","calcite.control", "fluid"))

#remove controls from data for plot
#exp.order <- Filter(function(x) !any(grepl("control", x)), exp.order)

cell.density.data <- cell.density.data %>% dplyr::filter(!Cell.density.sq.cm %in% c(NA, 0)) %>%
  #dplyr::filter(!Type %in% c("Internal Control", "internal Control")) %>%
  dplyr::group_by(Site, Exp.Type) %>%
  dplyr::summarise(density = mean(Cell.density.sq.cm))

cell.density.data$Exp.Type <- gsub(".fluid", "", cell.density.data$Exp.Type)
cell.density.data$site.experiment <- paste0(cell.density.data$Site, '.', cell.density.data$Exp.Type)
cell.density.data$site.experiment <-gsub("sand", "inert.control", cell.density.data$site.experiment)
cell.density.data$Exp.Type <- factor(cell.density.data$Exp.Type, levels = exp.order)

rarefied_alpha <- read.csv("rarefied_alpha_data.csv", header = TRUE, row.names = 1)

cell.count.alpha.div <- merge(rarefied_alpha_data, cell.density.data, by = "site.experiment")
cell.count.alpha.div$experiment <- plyr::mapvalues(cell.count.alpha.div$site.experiment, from=metadata$site.experiment, to=metadata$experiment)
  
diversity_vs_duration <- ggplot2::ggplot(cell.count.alpha.div, ggplot2::aes(density, duration))


cell.density.plot <- ggplot2::ggplot(cell.density.data, aes(density, Exp.Type)) +
  #geom_line(aes(color=Exp.Type), size=2.5, alpha=0.6) +
  #geom_point(aes(shape=Site, label=Date.Deployed)) + 
  geom_point(aes(color=Exp.Type), size=3) + 
  #scale_shape_manual(values=c(15,16,17)) +
  #labs(x=expression('Cell density ('~10^6~'cells/'~cm^2~')')) +
  labs(x=expression('Cell density (cells/'~cm^2~')')) +
  #scale_x_continuous(
  #breaks = c(0, 2000000, 4000000, 6000000),
  #label = c(0, 2, 4, 6)) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  #theme_gray() +
  theme(axis.title.y=element_blank(),
  legend.position = "none")

cell.density.plot +  facet_grid(cols = vars(Site))

density_vs_diversity <- ggplot(cell.count.alpha.div, aes(avg.abundance, density, color=Site.x, shape=experiment)) +
  geom_point() +
  scale_x_log10() +
  geom_smooth(method='lm',formula=y~x)

cell.density.plot <- ggplot(cell.density.data, aes(Cell.density.sq.cm, Exp.Type), shape=Site, group=Mineral) +
  geom_line(aes(color=Mineral), size=2.5, alpha=0.6) +
  #geom_point(aes(shape=Site, label=Date.Deployed)) + 
  geom_point(aes(label=Date.Deployed)) + 
  #scale_shape_manual(values=c(15,16,17)) +
  #labs(x=expression('Cell density ('~10^6~'cells/'~cm^2~')')) +
  labs(x=expression('Cell density (cells/'~cm^2~')')) +
  #scale_x_continuous(
      #breaks = c(0, 2000000, 4000000, 6000000),
      #label = c(0, 2, 4, 6)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_hline(yintercept=c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5, 15.5), colour="#3B3838", linetype="dotted") +
  #theme_gray() +
  theme(axis.title.y=element_blank())
        #legend.position = "none") 

cell.density.plot +  facet_grid(cols = vars(Site))
ggplotly(cell.density.plot +  facet_grid(cols = vars(Site)))


density.vs.date <- ggplot(cell.density.data, aes(Site.Date, Cell.density.sq.mm), group=Mineral) +
  geom_boxplot(color="gray") +
  geom_point(aes(y=Cell.density.sq.mm, color=Exp.Type)) + 
  coord_flip()


density.vs.site <- ggplot(cell.density.data, aes(Site, Cell.density.sq.mm)) +
  geom_boxplot(color="gray") +
  geom_point(aes(y=Cell.density.sq.mm, color=Exp.Type)) + 
  coord_flip()

density.vs.substrate <- ggplot(cell.density.data, aes(Exp.Type, Cell.density.sq.mm)) +
  geom_boxplot(color="gray") +
  geom_point(aes(y=Cell.density.sq.mm, color=Site)) +
  coord_flip()

density.vs.date.scatter <- ggplot(cell.density.data, aes(Date, Cell.density.sq.mm), group=Mineral) +
  geom_point(aes(color=Exp.Type, shape=Site)) + 
  geom_smooth(method='lm')


#cell densitites vs. thermo models
DeMMO_thermo_data <- read.csv("DeMMO_thermo_data.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

pyrolusite_densities <- dplyr::filter(cell.density.data, Mineral == "pyrolusite" )[,9]
names(pyrolusite_densities) <- dplyr::filter(cell.density.data, Mineral == "pyrolusite" )[,17]
pyrolusite_thermo <- dplyr::filter(DeMMO_thermo_data, mineral == "pyrolusite")[,2]
names(pyrolusite_thermo) <- rep(dplyr::filter(reactions, e.acceptor == "pyrolusite")[,3], 3)

pyrolusite_data <- expand.grid(density = pyrolusite_densities, delta_G = pyrolusite_thermo)
pyrolusite_data$site.exp <- plyr::mapvalues(pyrolusite_data$density, from=pyrolusite_densities, to=names(pyrolusite_densities))
pyrolusite_data$e.donor <- plyr::mapvalues(pyrolusite_data$delta_G, from=pyrolusite_thermo, to=names(pyrolusite_thermo))
pyrolusite_data <- pyrolusite_data %>% tidyr::separate(site.exp, c("site", "experiment"))

pyrolusite_plot <- ggplot(pyrolusite_data, aes(delta_G, density)) +
  geom_abline(intercept = 0, slope = 500) +
  geom_point(aes(shape=site, color=e.donor)) +
  scale_x_reverse() 

