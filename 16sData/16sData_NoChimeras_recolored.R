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
library(directlabels)

#import rarefied family-level OTU table from Qiime 
raw.data <- read.csv("DeMMO136_Dec2015toApril2018_noChimera_otuTable_withTaxa_d10000_categorized_L5.csv", header=TRUE, row.names = 1)
raw.data.unrarefied <- read.csv("DeMMO136_Dec2015toApril2018_noChimera_otuTable_withTaxa_categorized_L5.csv", header=TRUE, row.names = 1)

#identify taxa removed after rarefaction
rarefied.taxa <- as.data.frame(cbind(taxonomy = rownames(raw.data.unrarefied)[(!rownames(raw.data.unrarefied) %in% rownames(raw.data))]))

#clean up sample ID's and remove unwanted columns 
colnames(raw.data) = gsub("X", "", colnames(raw.data))
colnames(raw.data) = gsub("NU.", "", colnames(raw.data))
colnames(raw.data) = gsub("10dash1", "DeMMO1", colnames(raw.data))
colnames(raw.data) = gsub("24790", "DeMMO3", colnames(raw.data))
colnames(raw.data) = gsub("DuselB", "DeMMO6", colnames(raw.data))
colnames(raw.data) = gsub(".D1.", ".", colnames(raw.data))
colnames(raw.data) = gsub(".Steri.", ".fluid.", colnames(raw.data), ignore.case=TRUE)
raw.data <- raw.data[, -grep("DuselD", colnames(raw.data))]
colnames(raw.data) = gsub("24.DeMMO6.T6.bottom.051017", "24.D6.sand.051017", colnames(raw.data))
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

#format taxa names 
row.names(raw.data) <- gsub("D_0__Archaea;Other;Other;Other;Other", "Archaea.Other", row.names(raw.data))
row.names(raw.data) <- gsub("D_0__Bacteria;Other;Other;Other;Other", "Bacteria.Other", row.names(raw.data))
row.names(raw.data) <- gsub("D_0__Bacteria;D_1__Proteobacteria;Other;Other;Other", "Proteobacteria.Other", row.names(raw.data))
row.names(raw.data) <- gsub("D_0__Archaea;|D_0__Bacteria;|D_1__Proteobacteria;|D_1__|D_2__|D_3__|D_4__|D_5__", "", row.names(raw.data))
row.names(raw.data) <- gsub("Unknown;Other;Other;Other;Other", "Unassigned",row.names(raw.data))
row.names(raw.data) <- gsub(";", ".",row.names(raw.data))
row.names(raw.data) <- gsub("Gammaproteobacteria.Betaproteobacteriales", "Betaproteobacteria.Betaproteobacteriales",row.names(raw.data))

raw.data = raw.data[ rowSums(raw.data)!=0, ] 

#convert absolute abundanes to relative abundances - note that colsums = 100 (not 1)
data.rel.abundance <- raw.data %>%
  mutate_at(vars(1:length(colnames(raw.data))), funs(as.numeric(paste0(100*./sum(.)))))
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

data$Sample.Number <- paste0(data$Site, '.', data$experiment)
colnames(data)[2] <- 'Site.experiment'


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

#define shapes palette for nmds plot
shapes <- c(0, 15, 15, 1, 19, 19, 2, 17, 17, 5, 5, 5, 5)
names(shapes) <- c("D1.fluid", "D1.inert.control", "D1.mineral", "D3.fluid", "D3.inert.control", "D3.mineral", "D6.fluid", "D6.inert.control", "D6.mineral","D3.cont.control", "4800.cont.control", "800.cont.control", "4100L.fluid")

#import color ID's for phyla
phylum.colors.data <- as.data.frame(read.csv("phylumcolors_all_LM.csv", header=TRUE))
phylum.colors.data$hex.color <- as.character(phylum.colors.data$hex.color)

#generate phylum color key
phylum.colors <-c(rep(c('#1c1c1c', '#cccccc', '#4d4d4d'),3), rep('#afafae',4), paste0('#', phylum.colors.data$hex.color))
names(phylum.colors) <- c(names(shapes), as.character(unique(data$name)), as.character(unique(data$site)),as.character(phylum.colors.data$full.name))

#double check that all phyla have color assignments in the dataset
data_info <- cbind(phylum = gsub( "[.].*$", "", rownames(raw.data)))
rownames(data_info) <- rownames(raw.data)
data_info <- unique(data_info[,1])
#print out names of any phyla missing int the color key
data_info[(!data_info %in% phylum.colors.data$full.name)]

#check that all hex colors are valid
test <- paste0('#', as.character(phylum.colors.data$hex.color))
#if this throws an error, a hex code is invalid and needs to be corrected
col2rgb(test)




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
  
  
  #find row number of Desulfobulbaceae
  highlight.taxa <- c(which(grepl("Desulfobulbaceae", sig_vectors$family)), which(grepl("Thermodesulfovibrionia", sig_vectors$family)))
  highlight.taxa.df <- sig_vectors[highlight.taxa,]
  highlight.pyrolusite <- which(grepl("pyrolusite", merged_NMDS$Sample.Type))
  
  #filter sig_vectors for only taxa with lowest pvals
  #sig_vectors <- sig_vectors %>% filter(pval <= min(sig_vectors$pval))
  
  #subset sig-vectors for pvals < .001
  sig_vectors <- subset(sig_vectors, pval < 0.006)
 
  #Now, plot them like a badass
    plot <- ggplot(merged_NMDS, aes(x=MDS2, y=MDS1)) + 
    geom_point(aes(shape=Site.experiment, color=Site.experiment),size=2, alpha=0.8) +
    #geom_text(aes(label=paste0(Site, ".", Sample.Type),hjust = 1, vjust = 1),  size=2, color="black") +
    geom_segment(data=sig_vectors,inherit.aes = FALSE, aes(x=0,xend=MDS2,y=0,yend=MDS1, color=phylum, label=family), alpha=0.3)+
    geom_text(data=highlight.taxa.df,inherit.aes = FALSE,aes(x=MDS2, y=MDS1,label=family),size=4)+
    #geom_text(data=sig.vectors.phylum,inherit.aes = FALSE,aes(x=MDS2, y=MDS1,label=phylum),size=2, color="black", alpha=0.3)+
    geom_segment(data=highlight.taxa.df,inherit.aes = FALSE, aes(x=0,xend=MDS2,y=0,yend=MDS1), color="black", size=1, linetype = "dotted")+
    geom_point(data=merged_NMDS[highlight.pyrolusite,],inherit.aes = FALSE,aes(x=MDS2, y=MDS1),color="black", size=2, stroke=2) +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) + 
    #theme(legend.position="bottom", legend.box = "horizontal") +
    theme(legend.position = "none") +
    stat_ellipse(data=merged_NMDS, aes(color=Site.experiment)) +
    scale_shape_manual(values=shapes) + scale_color_manual(values=phylum.colors) +
    theme_grey()  
}

#generate NMDS plots for data + transformed data 
NMDS_plot.data <- NMDS_fun(data)
NMDS_plot.pres.abs <- NMDS_fun(data.pres.abs)
NMDS_plot.data.sqrt <- NMDS_fun(data.sqrt)
NMDS_plot.data.1.x <- NMDS_fun(data.1.x)

###Dendrogram and bar plots

#generate dendrogram using bray curtis dissimilarity metric
subset.data <- data[,6:ncol(data)] %>%
  rownames_to_column('Sample')  %>%
  filter(rownames(data)=="12.D1.fluid.041818" |rownames(data)=="26.D1.pyrolusite.041818" | rownames(data)=="22.D1.pyrite.041818"
         | rownames(data)=="23.D1.hematite.041818" | rownames(data)=="24.D1.magnetite.041818" | rownames(data)=="27.D1.siderite.041818"
         | rownames(data)=="34.D1.sand.041818" | rownames(data)=="45.D1.calcite.041818" | rownames(data)=="46.D1.gypsum.041818"
         | rownames(data)=="47.D1.muscovite.041818" | rownames(data)=="14.D3.fluid.041718" | rownames(data)=="51.D3.pyrolusite.041718"
         | rownames(data)=="10.D3.siderite.083017" | rownames(data)=="53.D3.pyrite.041718" | rownames(data)=="54.D3.hematite.041718"
         | rownames(data)=="55.D3.magnetite.041718" | rownames(data)=="39.D3.calcite.041718" | rownames(data)=="40.D3.gypsum.041718"
         | rownames(data)=="41.D3.muscovite.041718" | rownames(data)=="56.D3.sand.041718"  | rownames(data)=="12.D6.fluid.051017" | rownames(data)=="13.D6.pyrolusite.051017" 
         | rownames(data)=="15.D6.siderite.051017" | rownames(data)=="17.D6.pyrite.051017" | rownames(data)=="19.D6.hematite.051017"
         | rownames(data)=="21.D6.magnetite.051017" | rownames(data)=="24.D6.sand.051017" ) %>%
  column_to_rownames('Sample')

rownames(subset.data) <- gsub(".*[.]D", 'D', rownames(subset.data))
rownames(subset.data) <- gsub("[.]0.*", '', rownames(subset.data))

family.dendrogram <- as.dendrogram(hclust(bcdist(subset.data)))

#generate dendrogram
dendro.plot <- ggplot(family.dendrogram, horiz = T) + theme_gray()

#select rows where min abundance is greater than 10%
family.bar <-  subset.data %>%
   select_if(function(col) max(col) > 10)
family.bar$Less.Abundant.Taxa <- 100-rowSums(family.bar)

family.bar <- cbind(rownames(family.bar), family.bar)
colnames(family.bar)[1] <- "Sample"

#convert family.bar to long format   
family.bar <- gather(family.bar, Family, Abundance, colnames(family.bar[2]):colnames(family.bar[ncol(family.bar)]), factor_key=TRUE)

#add a column called Sample to store labels, this will ensure the bar plot and dendrogram align on labels
family.bar$Sample <- factor(family.bar$Sample, levels = labels(family.dendrogram))

#arrange bars in by size 
#family.bar$Sample <- reorder(family.bar$Sample, family.bar$Abundance)
#family.bar$Sample <- factor(family.bar$Sample, levels=rev(levels(family.bar$Sample)))

#generate color palette for family >10% community
n <- length(unique(family.bar$Family))
palette <- distinctColorPalette(n)
pie(rep(1, n), col=palette)

#import color ID's for families
family.colors.data <- as.data.frame(read.csv("familycolors.csv", header=TRUE))

#generate phylum color key
family.colors <-paste0('#', family.colors.data$color)
names(family.colors) <- as.character(family.colors.data$Family)

unique(family.bar$Family)[(!unique(family.bar$Family) %in% family.colors.data$Family)]

#generate family bar plot 
family.bar.plot.10percent <- ggplot(family.bar, aes(fill=Family, y=Abundance, x=Sample)) + 
  theme_gray() +
  geom_bar(stat='identity', position='fill') +
  #scale_fill_viridis_d() +
  scale_fill_manual(values=family.colors) +
  coord_flip() +
  theme(text = element_text(size=5), legend.position = "none")

#plot dendrogram and bar plot side by side 
dendro.bar.plot <- plot_grid(dendro.plot, family.bar.plot.10percent, align = "h")

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


###correlogram 
library(ggcorrplot)
library(GGally)
library(reshape2)
data$sample.ID <- rownames(data)
biofilm.data <-  data[3:ncol(data)] %>%
  filter(data$Sample.Type == "sand" | data$Sample.Type == "pyrolusite" |
              data$Sample.Type == "pyrite" | data$Sample.Type == "hematite" |
              data$Sample.Type == "magnetite" | data$Sample.Type == "siderite" |
              data$Sample.Type == "calcite" | data$Sample.Type == "muscovite") %>%
  select_if(function(col) max(col) >= 0.5) %>%
  filter(Sample.Date != "Sep.2016")

biofilm.data$Sample.Date <- c(rep("Nov.2017",4), rep("Aug.2017", 9), rep("Feb.2017", 6), rep("Aug.2017", 8), rep("Nov.2017", 5), rep("Feb.2017", 6))
biofilm.data$ID <- paste0(biofilm.data$Site, '.', biofilm.data$Sample.Type, '.', biofilm.data$Sample.Date)

cell.dens.data  <- read.csv("cellCounts.csv", header=TRUE)
cell.dens.data$Site <- gsub('DeMMO ', 'D', cell.dens.data$Site)
cell.dens.data$Site <- gsub('Demmo ', 'D', cell.dens.data$Site)
cell.dens.data$Mineral <- tolower(cell.dens.data$Mineral)
cell.dens.data$Mineral <- gsub('control sand', 'sand', cell.density.data$Mineral)

cell.dens.data$Type <- tolower(cell.dens.data$Type)
cell.dens.data$Type <- gsub('internal control', 'int.ctrl', cell.dens.data$Type)
cell.dens.data$Type <- gsub('control rep', 'rep', cell.dens.data$Type)
cell.dens.data$Type <- gsub('control', 'exp', cell.dens.data$Type)
cell.dens.data$Date <- paste0(lubridate::month(as.Date(cell.dens.data$Date.Deployed), label = TRUE),".", year(as.Date(cell.dens.data$Date.Deployed)))
cell.dens.data <- na.omit(cell.dens.data)
cell.dens.data$ID <- paste0(cell.dens.data$Site, '.', cell.dens.data$Mineral, '.', cell.dens.data$Date)

cell.dens.data <- dcast(cell.dens.data, ID ~ Type, value.var='Cell.density.sq.mm',fun.aggregate = toString)

biofilm.cell.dens.data <- merge(biofilm.data, cell.dens.data, by="ID", all=TRUE)
biofilm.cell.dens.data <- na.omit(biofilm.cell.dens.data)
drops <- c("sample.ID","ID")
biofilm.cell.dens.data<- biofilm.cell.dens.data[ , !(names(biofilm.cell.dens.data) %in% drops)]


biofilm.cell.dens.data <- biofilm.cell.dens.data[4:ncol(biofilm.cell.dens.data)] %>%
  select_if(function(col) max(col) > 0) 

biofilm.cell.dens.data$exp<-as.numeric(biofilm.cell.dens.data$exp)
biofilm.cell.dens.data$int.ctrl<-as.numeric(biofilm.cell.dens.data$int.ctrl)
biofilm.cell.dens.data$rep<-as.numeric(biofilm.cell.dens.data$rep)


sd <- as.data.frame(apply(biofilm.cell.dens.data[4:104], 2, sd))

ggcorr(biofilm.cell.dens.data)

corr <- round(cor(biofilm.cell.dens.data[4:104]), 1)
p.mat <- cor_pmat(biofilm.cell.dens.data[4:104])

DeMMO.correlogram <- 
  ggcorrplot(corr, hc.order = TRUE, type = "lower",
             outline.col = "white",
             ggtheme = ggplot2::theme_gray,
             colors = c("#6D9EC1", "white", "#E46726"))


### scatterplot of taxa variance
library(Rfast)
taxa.var <- as.data.frame(colVars(as.matrix(data[6:310])))
taxa.var$taxa <- colnames(data[6:310])
colnames(taxa.var)[1] <- "var"
ggplot(taxa.var, aes(x=reorder(taxa, var), y=var)) +
  theme_gray() +
  geom_bar(stat='identity') +
  expand_limits(x = 0, y = 0) +
  #geom_text(aes(label=taxa), vjust=0) +
  coord_flip()
  #geom_point() +

# rarefaction curves
rarefaction.data <- read.csv("alpha_rare/collated_observed_otus/observed_otus.csv", header=TRUE)

colnames(rarefaction.data) = gsub("X", "", colnames(rarefaction.data))
colnames(rarefaction.data) = gsub("NU.", "", colnames(rarefaction.data))
colnames(rarefaction.data) = gsub("10dash1", "DeMMO1", colnames(rarefaction.data))
colnames(rarefaction.data) = gsub("24790", "DeMMO3", colnames(rarefaction.data))
colnames(rarefaction.data) = gsub("DuselB", "DeMMO6", colnames(rarefaction.data))
colnames(rarefaction.data) = gsub(".D1.", ".", colnames(rarefaction.data))
colnames(rarefaction.data) = gsub(".Steri.", ".fluid.", colnames(rarefaction.data), ignore.case=TRUE)
rarefaction.data <- rarefaction.data[, -grep("DuselD", colnames(rarefaction.data))]
colnames(rarefaction.data) = gsub("24.DeMMO6.T6.bottom.051017", "24.D6.sand.051017", colnames(rarefaction.data))
rarefaction.data <- rarefaction.data[, -grep("bottom", colnames(rarefaction.data))]
rarefaction.data <- rarefaction.data[, -grep("AfterPacker", colnames(rarefaction.data))]
colnames(rarefaction.data) = gsub("AfterDrilling.", "", colnames(rarefaction.data))
colnames(rarefaction.data) <- gsub("[.]A[.]top[.]|[.]SC1[.]top[.]|[.]T1[.]top[.]|[.]T7[.]top[.]", ".pyrolusite.", colnames(rarefaction.data))
colnames(rarefaction.data) <- gsub("[.]T2[.]top[.]|[.]T8[.]top[.]|[.]SC2[.]top[.]", ".siderite.", colnames(rarefaction.data))
colnames(rarefaction.data) <- gsub("[.]C[.]top[.]|[.]SC3[.]top[.]|[.]T3[.]top[.]|[.]T9[.]top[.]", ".pyrite.", colnames(rarefaction.data))
colnames(rarefaction.data) <- gsub("[.]D[.]top[.]|[.]SC4[.]top[.]|[.]T10[.]top[.]|[.]T4[.]top[.]", ".hematite.", colnames(rarefaction.data))
colnames(rarefaction.data) <- gsub("[.]E[.]top[.]|[.]SC5[.]top[.]|[.]T11[.]top[.]|[.]T5[.]top[.]", ".magnetite.", colnames(rarefaction.data))
colnames(rarefaction.data) <- gsub("[.]F[.]top[.]|[.]SC10[.]top[.]|[.]T12[.]top[.]|[.]T6[.]top[.]|[.]12[.]top[.]|[.]6[.]top[.]", ".sand.", colnames(rarefaction.data))
colnames(rarefaction.data) <- gsub("[.]1[.]top.|[.]7[.]top[.]", ".calcite.", colnames(rarefaction.data))
colnames(rarefaction.data) <- gsub("[.]2[.]top.|[.]8[.]top[.]", ".gypsum.", colnames(rarefaction.data))
colnames(rarefaction.data) <- gsub("[.]3[.]top[.]|[.]9[.]top[.]", ".muscovite.", colnames(rarefaction.data))
colnames(rarefaction.data) <- gsub("[.]SC7[.]top[.]", ".3mmPyrex.", colnames(rarefaction.data))
colnames(rarefaction.data) <- gsub("[.]SC8[.]top[.]", ".5mmPyrex.", colnames(rarefaction.data))
colnames(rarefaction.data) <- gsub("[.]SC7[.]top[.]|[.]SC9[.]top[.]", ".wool.", colnames(rarefaction.data))
rarefaction.data <- rarefaction.data[, -grep("[.]B[.]top[.]|[.]4[.]top[.]|[.]5[.]top[.]|[.]10[.]top[.]|[.]11[.]top[.]", colnames(rarefaction.data))]
colnames(rarefaction.data) <- gsub("[.]2[.]|[.]1[.]", ".", colnames(rarefaction.data))
colnames(rarefaction.data) <- gsub("1[.]0", "1.fluid.0", colnames(rarefaction.data))
colnames(rarefaction.data) <- gsub("3[.]0", "3.fluid.0", colnames(rarefaction.data))
colnames(rarefaction.data) <- gsub("6[.]0", "6.fluid.0", colnames(rarefaction.data))
colnames(rarefaction.data) <- gsub("DeMMO", "D", colnames(rarefaction.data))
colnames(rarefaction.data)[1] <- "rarefaction.id"
rarefaction.data$rarefaction.id <- gsub("_0.txt|_1.txt|_2.txt|_3.txt|_4.txt|_5.txt|_6.txt|_7.txt|_8.txt|_9.txt", "", rarefaction.data$rarefaction.id)

#remove any rows that sum to zero
rarefaction.data = rarefaction.data[ rowSums(rarefaction.data[3:length(colnames(rarefaction.data))], na.rm = TRUE)!=0, ] 

#aggregate by rarefaction id
rarefaction.summary <- rarefaction.data %>% 
  group_by(as.character(rarefaction.id)) %>% 
  summarize_at(vars(-rarefaction.id,-sequences.per.sample,-iteration), dplyr::funs(mean))

#reformat rarefaction summary for gathering
rarefaction.summary <- as.data.frame(rarefaction.summary)
rownames(rarefaction.summary) <- rarefaction.summary$`as.character(rarefaction.id)`
rarefaction.summary <- rarefaction.summary[-1]

rarefaction.summary <- as.data.frame(t(rarefaction.summary))
rarefaction.summary <- cbind(sample.id = rownames(rarefaction.summary), rarefaction.summary)

rarefaction.summary <- merge(data[,1:5], rarefaction.summary, by="row.names", all=TRUE)

#gather rarefaction data into long format for plotting
rarefaction.summary <- gather(rarefaction.summary, rarefaction.depth, OTU.abundance, 
                              colnames(rarefaction.summary[8]):colnames(rarefaction.summary[length(colnames(rarefaction.summary))]), factor_key=TRUE)

rarefaction.summary$rarefaction.depth <- gsub("alpha_rarefaction_", "", rarefaction.summary$rarefaction.depth)
rarefaction.summary$rarefaction.depth <- as.numeric(rarefaction.summary$rarefaction.depth)
rarefaction.summary$Label <- paste0(rownames(rarefaction.summary), ".", rarefaction.summary$Sample.Type)

rarefaction.summary <- na.omit(rarefaction.summary)
#plot rarefaction curves 
rarefaction.plot <- ggplot(rarefaction.summary, aes(rarefaction.depth, OTU.abundance, group=sample.id, color=Site.experiment)) +
  geom_line() + 
  theme(legend.position = "none") +
  scale_colour_discrete(guide = 'none') +
  #theme(legend.position = c(.1, .84), legend.text=element_text(size=6), legend.title = element_text(size=8, face="bold")) +
  #theme(legend.key.size =  unit(0.1, "in")) +
  #scale_x_continuous(expand = c(0.15, 0)) +
  #geom_dl(aes(label = sample.id), method = list(dl.combine("first.points", "last.points"), cex = 0.8)) +
  theme_grey()


#ranking diversity

#set diversity depth to 3910
diversity.ranking <- rarefaction.summary %>% 
  filter(rarefaction.depth == 3910)

#set diversity depth to 9760
diversity.ranking <- rarefaction.summary %>% 
  filter(rarefaction.depth == 9760)

diversity.ranking$site.type <- paste0(diversity.ranking$Site, ".",diversity.ranking$Sample.Type)
diversity.ranking$rarefaction.depth <- as.character(diversity.ranking$rarefaction.depth)
diversity.ranking$Site.experiment <- gsub("800.cont.control|D3.cont.control|4100L.fluid|4800.cont.control", "ambient.control", diversity.ranking$Site.experiment)

diversity.plot <- ggplot(diversity.ranking, aes(x=reorder(Site.experiment, -OTU.abundance), y=OTU.abundance, fill=Site.experiment)) + 
  geom_violin() +
  geom_boxplot(width=0.1) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  coord_flip() +
  theme_grey() +
  theme(legend.position = "none") 