###=========================###
NMDS using OTU table from qiime
###=========================###

library(ggplot2)
library(ggfortify)
library(vegan)
library(grid)
library(RColorBrewer)

#https://stackoverflow.com/questions/14711470/plotting-envfit-vectors-vegan-package-in-ggplot2

#calculate distance for NMDS
family_all<- read.csv("taxa_family_all_D3MayRemoved_unaltered.csv", row.name=1)

#include the outliers
family_all_NMDS <- metaMDS(family_all[,7:851],k=2)

#take out the outliers
#family_all<-family_all[1:101,]
#family_all_NMDS <- metaMDS(family_all[,6:850],k=2)

#save results in data.frame
NMDS.frame = data.frame(MDS1 = family_all_NMDS$points[,1], MDS2 = family_all_NMDS$points[,2])

#merge columns
merged_NMDS <- merge(NMDS.frame, family_all[,1:6], by="row.names", all=TRUE)
merged_NMDS$grp <- paste(merged_NMDS$Substrate,merged_NMDS$Date)

#fit vectors to family
vec.sp<-envfit(family_all_NMDS$points, family_all[,7:851], perm=1000)
vec.sp.df<-as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$pval <- vec.sp$vectors$pvals

#sort vectors by smallest p value
sig_vectors <- subset(vec.sp.df, pval <= 0.05, select=c(MDS1, MDS2, pval))
sig_vectors <- sig_vectors[order(sig_vectors$pval),] 
sig_vectors$species<-rownames(sig_vectors)
write.csv(sig_vectors, "sig_vectors_30Sep2018.csv")
sig_vector_to_keep<- read.csv("merge_sig_vec.csv", row.names=1,header = TRUE)
sig_vec_test<- merge(sig_vectors, sig_vector_to_keep, by="row.names")


#Now, plot them like a badass
NMDS_plot <- ggplot(merged_NMDS, aes(x=MDS1, y=MDS2, color=Site)) + 
  geom_point(size=4, stroke=2) +
  geom_text(aes(label=Label,hjust = 1, vjust = 1),  size=4, color="black") +
  geom_segment(data=sig_vec_test,inherit.aes = FALSE, aes(x=0,xend=MDS1,y=0,yend=MDS2), arrow = arrow(length = unit(0.2, "cm")),colour="grey")+
  geom_text(data=sig_vec_test,inherit.aes = FALSE,aes(x=MDS1,y=MDS2,label=species_name),size=4)+
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) + 
  stat_ellipse()


NMDS_plot_bf_vs_fluid <- ggplot(merged_NMDS, aes(x=MDS1, y=MDS2, color=Site, shape=type)) + 
  geom_point(size=4, stroke=2) +
  geom_text(aes(label=Label,hjust = 1, vjust = 1),  size=4, color="black") +
  geom_segment(data=sig_vec_test,inherit.aes = FALSE, aes(x=0,xend=MDS1,y=0,yend=MDS2), arrow = arrow(length = unit(0.2, "cm")),colour="grey")+
  geom_text(data=sig_vec_test,inherit.aes = FALSE,aes(x=MDS1,y=MDS2,label=species_name),size=4)+
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) + 
  stat_ellipse()