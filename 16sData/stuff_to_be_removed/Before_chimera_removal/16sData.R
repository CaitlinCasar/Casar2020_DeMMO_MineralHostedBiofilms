library(tidyverse)
library(ggplot2)
library(vegan)
library(grid)
library(lubridate)
library(qdap)
library(ecodist)
library(dendextend)
library(cowplot)
library(RColorBrewer)
library(randomcoloR)
library(plotly)

#import rarefied family-level OTU table from Qiime 
raw.data <- read.csv("DeMMO_Dec2015_to_Aug2017_d3581_L5.csv", header=TRUE, row.names = 1)

#clean up sample ID's and remove unwanted columns 
colnames(raw.data) = gsub("X", "", colnames(raw.data))
colnames(raw.data) = gsub("NU.", "", colnames(raw.data))
colnames(raw.data) = gsub("10dash1", "DeMMO1", colnames(raw.data))
colnames(raw.data) = gsub("24790", "DeMMO3", colnames(raw.data))
colnames(raw.data) = gsub("DuselB", "DeMMO6", colnames(raw.data))
colnames(raw.data) = gsub(".D1.", ".", colnames(raw.data))
colnames(raw.data) = gsub(".Steri.", ".fluid.", colnames(raw.data), ignore.case=TRUE)
raw.data <- raw.data[, -grep("DuselD", colnames(raw.data))]
colnames(raw.data)[112] <- "24.D6.sand.051017"
raw.data <- raw.data[, -grep("bottom", colnames(raw.data))]
raw.data <- raw.data[, -grep("AfterPacker", colnames(raw.data))]
colnames(raw.data) = gsub("AfterDrilling.", "", colnames(raw.data))
colnames(raw.data) <- gsub("[.]A[.]top[.]|[.]SC1[.]top[.]|[.]T1[.]top[.]|[.]T7[.]top[.]", ".pyrolusite.", colnames(raw.data))
colnames(raw.data) <- gsub("[.]T2[.]top[.]|[.]T8[.]top[.]|[.]SC2[.]top[.]", ".siderite.", colnames(raw.data))
colnames(raw.data) <- gsub("[.]C[.]top[.]|[.]SC3[.]top[.]|[.]T3[.]top[.]|[.]T9[.]top[.]", ".pyrite.", colnames(raw.data))
colnames(raw.data) <- gsub("[.]D[.]top[.]|[.]SC4[.]top[.]|[.]T10[.]top[.]|[.]T4[.]top[.]", ".hematite.", colnames(raw.data))
colnames(raw.data) <- gsub("[.]E[.]top[.]|[.]SC5[.]top[.]|[.]T11[.]top[.]|[.]T5[.]top[.]", ".magnetite.", colnames(raw.data))
colnames(raw.data) <- gsub("[.]F[.]top[.]|[.]SC10[.]top[.]|[.]T12[.]top[.]|[.]T6[.]top[.]|[.]12[.]top[.]|[.]6[.]top[.]", ".sand.", colnames(raw.data))
colnames(raw.data) <- gsub("[.]1[.]top.|[.]7[.]top[.]", ".calcite.", colnames(raw.data))
colnames(raw.data) <- gsub("[.]2[.]top.|[.]8[.]top[.]", ".gypsum.", colnames(raw.data))
colnames(raw.data) <- gsub("[.]3[.]top[.]|[.]9[.]top[.]", ".muscovite.", colnames(raw.data))
colnames(raw.data) <- gsub("[.]SC7[.]top[.]", ".3mmPyrex.", colnames(raw.data))
colnames(raw.data) <- gsub("[.]SC8[.]top[.]", ".5mmPyrex.", colnames(raw.data))
colnames(raw.data) <- gsub("[.]SC7[.]top[.]|[.]SC9[.]top[.]", ".wool.", colnames(raw.data))
raw.data <- raw.data[, -grep("[.]B[.]top[.]|[.]4[.]top[.]|[.]5[.]top[.]|[.]10[.]top[.]|[.]11[.]top[.]", colnames(raw.data))]
colnames(raw.data) <- gsub("[.]2[.]|[.]1[.]", ".", colnames(raw.data))
colnames(raw.data) <- gsub("1[.]0", "1.fluid.0", colnames(raw.data))
colnames(raw.data) <- gsub("3[.]0", "3.fluid.0", colnames(raw.data))
colnames(raw.data) <- gsub("6[.]0", "6.fluid.0", colnames(raw.data))
colnames(raw.data) <- gsub("DeMMO", "D", colnames(raw.data))
row.names(raw.data) <- gsub("D_0__Archaea;|D_0__Bacteria;|D_1__Proteobacteria;|D_1__|D_2__|D_3__|D_4__|D_5__", "", row.names(raw.data))
row.names(raw.data) <- gsub(";", ".",row.names(raw.data))
row.names(raw.data)[844] <- "Bacteria.Other"
row.names(raw.data)[845] <- "Unassigned"
raw.data = raw.data[ rowSums(raw.data)!=0, ] 

#convert absolute abundanes to relative abundances - note that colsums = 100 (not 1)
data.rel.abundance <- raw.data %>%
  mutate_at(vars(1:91), funs(as.numeric(paste0(100*./sum(.)))))
rownames(data.rel.abundance) <- rownames(raw.data)

#transform data and separate name into metadata columns
data <- as.data.frame(t(data.rel.abundance))
data <- cbind(Sample.ID = rownames(data), data)
data <- data %>% separate(Sample.ID, c("Sample.Number", "Site", "Sample.Type", "Sample.Date"))
data$Sample.Date <- paste0(lubridate::month(mdy(data$Sample.Date), label = TRUE),".", year(mdy(data$Sample.Date)))

#add sample type to data frame 
categorize.samples <- function(x){
  minerals <- c("pyrolusite","pyrite","hematite","magnetite", "siderite","calcite","gypsum","muscovite")
  inert.controls <- c("sand", "wool", "5mmPyrex", "3mmPyrex")
  cont.controls <- c("DitchFluid","SlideBlank")
  fluid <- "fluid"
  if (x[3] %in% minerals) {
    print("mineral")
  } else if (x[3] %in% inert.controls){
    print("inert.control")
  } else if (x[3] %in% cont.controls){
    print("cont.control")
  } else {
    print("fluid")
  }
}

#add categories to data frame
data <- cbind(experiment = apply(data, 1,categorize.samples), data)

#transform data to binary presence/absence
data.pres.abs <- + (data[6:ncol(data)] > 0)
data.pres.abs <- cbind(data[1:5], data.pres.abs)

#square root transform data
data.sqrt <- sqrt(data[6:ncol(data)])
data.sqrt <- cbind(data[1:5], data.sqrt)

#1/x transform data
data.1.x <- 1/(data[6:ncol(data)])
data.1.x <- do.call(data.frame, lapply(data.1.x, function(x) {
  replace(x, is.infinite(x), 0)
})
)
data.1.x <- cbind(data[1:5], data.1.x)

###NMDS Plot

NMDS_fun <- function(z){
  #generate NMDS coordinates, exclude first 5 columns of metadata
  NMDS <- metaMDS(z[,6:ncol(data)],k=2)

  #save results in data.frame
  NMDS.frame = data.frame(MDS1 = NMDS$points[,1], MDS2 = NMDS$points[,2])

  #combine NMDS coordinates with metadata
  merged_NMDS <- merge(NMDS.frame, data[,1:5], by="row.names", all=TRUE)

  #format sample dates 
  #merged_NMDS$Sample.Date <- paste0(lubridate::month(mdy(merged_NMDS$Sample.Date), label = TRUE),".", year(mdy(merged_NMDS$Sample.Date)))


  #fit vectors to family
  family.vectors <-envfit(NMDS$points, data[,6:ncol(data)], perm=1000)
  family.vectors.df<-as.data.frame(family.vectors$vectors$arrows*sqrt(family.vectors$vectors$r))
  family.vectors.df$pval <- family.vectors$vectors$pvals

  #sort vectors by smallest p value
  sig_vectors <- subset(family.vectors.df, pval <= 0.05, select=c(MDS1, MDS2, pval))
  sig_vectors <- sig_vectors[order(sig_vectors$pval),] 
  sig_vectors$family<-rownames(sig_vectors)
  
  #add phylum column 
  sig_vectors$phylum <- gsub( "[.].*$", "", sig_vectors$family)

#create new data frame with unique phyla with lowest pvals - these will be labels on NMDS plot
  sig.vectors.phylum <- sig_vectors %>% 
    group_by(phylum) %>% 
    slice(which.min(pval))

  #find row number of Desulfobulbaceae
  highlight.taxa <- c(which(grepl("Desulfobulbaceae", sig_vectors$family)), which(grepl("Thermodesulfovibrionia", sig_vectors$family)))
  highlight.pyrolusite <- which(grepl("pyrolusite", merged_NMDS$Sample.Type))

  site.palette <- gray.colors(length(unique(merged_NMDS$Site)), start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL)
  
  #Now, plot them like a badass
    ggplot(merged_NMDS, aes(x=MDS2, y=MDS1)) + 
    geom_point(data=merged_NMDS[which(grepl("D1", merged_NMDS$Site)),],inherit.aes = FALSE,aes(x=MDS2, y=MDS1),color=site.palette[1], size=2,alpha=0.8) +
    geom_point(data=merged_NMDS[which(grepl("D3", merged_NMDS$Site)),],inherit.aes = FALSE,aes(x=MDS2, y=MDS1),color=site.palette[2], size=2,alpha=0.8) +
    geom_point(data=merged_NMDS[which(grepl("D6", merged_NMDS$Site)),],inherit.aes = FALSE,aes(x=MDS2, y=MDS1),color=site.palette[3], size=2, alpha=0.8) +
    geom_point(data=merged_NMDS[which(grepl("800", merged_NMDS$Site)),],inherit.aes = FALSE,aes(x=MDS2, y=MDS1),color=site.palette[4], size=2, alpha=0.8) +
    geom_point(data=merged_NMDS[which(grepl("4800", merged_NMDS$Site)),],inherit.aes = FALSE,aes(x=MDS2, y=MDS1),color=site.palette[5], size=2, alpha=0.8) +
    geom_text(aes(label=paste0(Site, ".", Sample.Type),hjust = 1, vjust = 1),  size=2, color="black") +
    geom_segment(data=sig_vectors,inherit.aes = FALSE, aes(x=0,xend=MDS2,y=0,yend=MDS1, color=phylum, label=family), alpha=0.3)+
    geom_text(data=sig_vectors[highlight.taxa,],inherit.aes = FALSE,aes(x=MDS2, y=MDS1,label=family),size=4)+
    geom_text(data=sig.vectors.phylum,inherit.aes = FALSE,aes(x=MDS2, y=MDS1,label=phylum),size=2, color="black", alpha=0.3)+
    geom_segment(data=sig_vectors[highlight.taxa,],inherit.aes = FALSE, aes(x=0,xend=MDS2,y=0,yend=MDS1), color="black", size=1, linetype = "dotted")+
    geom_point(data=merged_NMDS[highlight.pyrolusite,],inherit.aes = FALSE,aes(x=MDS2, y=MDS1),color="black", size=2, stroke=2) +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) + 
    #theme(legend.position="bottom", legend.box = "horizontal") +
    stat_ellipse(data=merged_NMDS[which(grepl("D1", merged_NMDS$Site)),],inherit.aes = FALSE,aes(x=MDS2, y=MDS1), color=site.palette[1]) +
    stat_ellipse(data=merged_NMDS[which(grepl("D3", merged_NMDS$Site)),],inherit.aes = FALSE,aes(x=MDS2, y=MDS1), color=site.palette[2]) +
    stat_ellipse(data=merged_NMDS[which(grepl("D6", merged_NMDS$Site)),],inherit.aes = FALSE,aes(x=MDS2, y=MDS1), color=site.palette[3]) +
    stat_ellipse(data=merged_NMDS[which(grepl("800", merged_NMDS$Site)),],inherit.aes = FALSE,aes(x=MDS2, y=MDS1), color=site.palette[4]) +
    stat_ellipse(data=merged_NMDS[which(grepl("4800", merged_NMDS$Site)),],inherit.aes = FALSE,aes(x=MDS2, y=MDS1), color=site.palette[5]) +
      theme_grey() 
    }

#generate NMDS plots for data + transformed data 
NMDS_plot.data <- NMDS_fun(data)
NMDS_plot.pres.abs <- NMDS_fun(data.pres.abs)
NMDS_plot.data.sqrt <- NMDS_fun(data.sqrt)
NMDS_plot.data.1.x <- NMDS_fun(data.1.x)

#create interactive html plots 
ggplotly(NMDS_plot.data.1.x  + theme(legend.position="none"))

###Dendrogram and bar plots

#generate dendrogram using bray curtis dissimilarity metric
family.dendrogram <- as.dendrogram(hclust(bcdist(data[,6:ncol(data)])))

#generate dendrogram
dendro.plot <- ggplot(family.dendrogram, horiz = T)

#select rows where min abundance is greater than 15%
family.bar <-  data[6:ncol(data)] %>%
   select_if(function(col) max(col) > 15)
family.bar$Less.Abundant.Taxa <- 100-rowSums(family.bar)

family.bar <- cbind(rownames(data), family.bar)
colnames(family.bar)[1] <- "Sample"

#convert family.bar to long format   
family.bar <- gather(family.bar, Family, Abundance, colnames(family.bar[2]):colnames(family.bar[ncol(family.bar)]), factor_key=TRUE)

#add a column called Sample to store labels, this will ensure the bar plot and dendrogram align on labels
family.bar$Sample <- factor(family.bar$Sample, levels = labels(family.dendrogram))

#arrange bars in by size 
#family.bar$Sample <- reorder(family.bar$Sample, family.bar$Abundance)
#family.bar$Sample <- factor(family.bar$Sample, levels=rev(levels(family.bar$Sample)))

family.bar.plot.15percent <- ggplot(family.bar, aes(fill=Family, y=Abundance, x=Sample)) + 
  geom_col() + 
  geom_bar(stat='identity', position='fill') +
  #scale_fill_viridis_d() +
  coord_flip()

#filter data for families greater than 1% community - low abundance taxa >15% will be plotted in grayscale
family.bar.low.ab <-  data[6:ncol(data)] %>%
  select_if(function(col) max(col) > 1)
family.bar.low.ab$Less.Abundant.Taxa <- 100-rowSums(family.bar.low.ab)

family.bar.low.ab <- cbind(rownames(data), family.bar.low.ab)
colnames(family.bar.low.ab )[1] <- "Sample"

#convert family.bar.low.ab to long format   
family.bar.low.ab <- gather(family.bar.low.ab, Family, Abundance, colnames(family.bar.low.ab[2]):colnames(family.bar.low.ab[ncol(family.bar.low.ab)]), factor_key=TRUE)

#add a column called Sample to store labels, this will ensure the bar plot and dendrogram align on labels
family.bar.low.ab$Sample <- factor(family.bar.low.ab$Sample, levels = labels(family.dendrogram))

#arrange bars by size 
family.bar.low.ab$Sample <- reorder(family.bar.low.ab$Sample, family.bar.low.ab$Abundance)
family.bar.low.ab$Sample <- factor(family.bar.low.ab$Sample, levels=rev(levels(family.bar.low.ab$Sample)))

#generate color palette for family >15% community
n <- length(unique(family.bar$Family))
palette <- distinctColorPalette(n)
pie(rep(1, n), col=palette)

#generate grayscale palette for families < 15% community
gray.palette <- gray.colors(length(unique(family.bar.low.ab$Family))-n, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL)

#combine color and grayscale palettes 
mixed.palette <- c(gray.palette,palette)
#pie(rep(1, length(unique(family.bar.low.ab$Family))), col=mixed.palette)

#determine order of family abundances 
fam.abundances = aggregate(family.bar.low.ab$Abundance, list(family.bar.low.ab$Family), max)
fam.abundances <- fam.abundances[order(fam.abundances[2]),] 
fam.abundances$Group.1 <- as.character(fam.abundances$Group.1)

#assign colors to families if > 15% of community, else grayscale
assigned.mixed.palette <- setNames(as.list(as.character(mixed.palette)), fam.abundances$Group.1)

#generate bar plot using mixed color palette
family.bar.plot <- ggplot(family.bar.low.ab, aes(fill=Family, y=Abundance, x=Sample)) + 
  geom_col() + 
  geom_bar(stat='identity', position='fill') +
  #scale_fill_viridis_d() +
  coord_flip() #+
  #scale_x_discrete(limits = rev(levels(family.bar.plot$Sample))) +
  #scale_fill_manual(values=assigned.mixed.palette) +
  #theme(legend.position="none")

#plot dendrogram and bar plot side by side 
dendro.bar.plot <- plot_grid(dendro.plot, family.bar.plot, align = "h")


###Area Plot

#select rows where min abundance is greater than 1%
family.area <-  data[6:ncol(data)] %>%
  select_if(function(col) max(col) > 1)
family.area$Less.Abundant.Taxa <- 100-rowSums(family.area)

family.area <- cbind(data[1:5], family.area)

family.area <- cbind(rownames(data), family.area)
colnames(family.area)[1] <- "Sample"

#convert family.area to long format   
family.area <- gather(family.area, Family, Abundance, colnames(family.area[7]):colnames(family.area[ncol(family.area)]), factor_key=TRUE)


library(randomcoloR)
n <- length(unique(family.area$Family))
palette <- distinctColorPalette(n)
#pie(rep(1, n), col=palette)

#organize data frame by family and site.exp
family.area$site.exp <- paste0(family.area$Site, ".", family.area$experiment)
family.area <- family.area[
  with(family.area, order(Family, desc(site.exp))),
  ]


#add column of unique numeric sample ID's for area chart (cannot use strings)
family.area <- cbind(c(1:91), family.area)


colnames(family.area)[1] <- "Sample.ID"


Sample.labels <- with(family.area,setNames(as.character(unique(Sample)),as.character(unique(Sample.ID))))

#set colors for each experiment type label on plot 
experiment.color.fun <- function(q){
  if(q[3]=="fluid"){
    print("red")
  } else if (q[3] == "inert.control"){
    print("blue")
  } else if(q[3] =="mineral"){
    print("orange")
  } else {
    print("purple")
  }
}
family.area <- cbind(exp.label.color = apply(family.area, 1, experiment.color.fun),family.area)
exp.label.color <- with(family.area,setNames(as.character(exp.label.color),as.character(Sample.ID)))



#generate family area plot 
family.area.plot <- ggplot(family.area[order(family.area$Family),], aes(x=Sample.ID, y=Abundance, fill=Family)) +
  scale_fill_manual(values=palette) +
  geom_area() +
  theme(legend.position="bottom", legend.box = "horizontal") +
  scale_x_discrete(limits=1:length(unique(family.area$Sample.ID)), breaks=as.character(unique(family.area$Sample.ID)),labels=Sample.labels) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color=exp.label.color, size=3)) +
  geom_vline(xintercept=seq(1, length(unique(family.area$Sample.ID)), by=1), color="white", linetype="dotted")


###alternative color palette methods

#set numer of colors needed for color palette to represent families, store in col_vector
n <- ncol(family.bar)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#visualize colors stored in col_vector
pie(rep(1,n), col=sample(col_vector, n))

#color palettes for color blindness 
library(viridis)
n=100
viridis.palette <- viridis::viridis_pal(option = "D")(n)
pie(rep(1, n), viridis.palette)




