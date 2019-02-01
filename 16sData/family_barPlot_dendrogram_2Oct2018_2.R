library("ggplot2")
library(tidyr)
library(ecodist)
#install.packages("dendextend")
library(dendextend)
#install.packages("cowplot")
library(cowplot)

#code source: http://www.cookbook-r.com/Manipulating_data/Converting_data_between_wide_and_long_format/
#code source: https://stackoverflow.com/questions/44646488/stacked-bar-plot-with-hierarchical-clustering-dendrogram

#read in data 
#family_bar_data<-read.csv("taxa_family_bar_Fe.csv")
family_bar_data<-read.csv("taxa_family_all_D3MayRemoved_collapsed.csv")2
#family <- read.csv("taxa_dendro_Fe.csv", row.names = 1)
family <- read.csv("taxa_family_all_D3MayRemoved_collapsed.csv", row.names = 1)

#generate dendrogram using bray curtis dissimilarity metric
family_dendrogram <- as.dendrogram(hclust(bcdist(family[,6:850])))

#generate long-format table of data
family_bar_long <- gather(family_bar_data, Family, Abundance, Unassigned:Less_Abundant_Taxa, factor_key=TRUE)

#add a column called Sample to store labeal, this will ensure the bar plot and dendrogram align on labels
family_bar_long$Sample <- factor(family_bar_long$Sample, levels = labels(family_dendrogram))

#create plots 
dendro_plot <- ggplot(family_dendrogram, horiz = T)

abundance_plot <- ggplot(family_bar_long, aes(fill=Family, y=Abundance, x=Sample)) + 
  geom_bar(stat='identity', position='fill') + 
  coord_flip() #+
  #theme(legend.position="none")

plot_grid(dendro_plot, abundance_plot, align = "h")

legend <- get_legend(abundance_plot + theme(legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom"))

plot(legend)
