library(ggplot2)
library(plotly)
library(lubridate)


#import cell count data
cell.density.data <- read.csv("cellCounts.csv", header=TRUE)
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

cell.density.plot <- ggplot(cell.density.data, aes(Cell.density.sq.mm, reorder(Exp.Type, -Cell.density.sq.mm)), shape=Site, group=Mineral) +
  geom_line(aes(color=Mineral), size=2.5, alpha=0.6) +
  geom_point(aes(shape=Site, label=Date.Deployed)) + 
  scale_shape_manual(values=c(15,16,17)) +
  xlab("Cell density (10^4 cells/mm^2)") +
  scale_x_continuous(
      breaks = c(0, 20000, 40000, 60000),
      label = c(0, 2, 4, 6)) +
  theme_gray()
  #theme(axis.title.y=element_blank(), legend.position = "none") 

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

