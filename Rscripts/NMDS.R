#load dependencies 
pacman::p_load(plotly, vegan, tidyverse, readr, plyr,lubridate)

#import otu table, metadata, taxonomy
otu_table <- read_delim("../orig_data/DeMMO136_Dec2015toApril2018_noChimera_otuTable_withTaxa_d10000.txt", delim = '\t', comment = "# ")
metadata <- read_csv("../orig_data/metadata.csv") 

#create taxonomy table updated with GTDB 
gtdb_taxonomy <- read_csv("../orig_data/Unassigned_assignments.csv") %>%
  select(Query_Seq_ID, GTDBtk_taxonomy) %>%
  rename(`#OTU ID` = Query_Seq_ID) %>%
  rename(taxonomy = GTDBtk_taxonomy)

taxonomy <- otu_table %>%
  select(`#OTU ID`, taxonomy) %>%
  filter(!`#OTU ID` %in% gtdb_taxonomy$`#OTU ID`) %>%
  bind_rows(gtdb_taxonomy) %>%
  mutate(tax = gsub("Gammaproteobacteria; D_3__Betaproteobacteriales", "Betaproteobacteria; D_3__Betaproteobacteriales", taxonomy), #fix taxonomy for Beta's,
         taxonomy = str_remove_all(tax, "D_0__| D_1__| D_2__| D_3__| D_4__| D_5__| D_6__")) %>%
  separate(taxonomy ,sep=';',c("domain", "phylum", "class", "order", "family", "genus", "species"))


# NMDS on OTU data with vectors representing phyla  -----------------------
otu_norm <- otu_table %>%
  select(-taxonomy) %>%
  mutate_at(vars(-`#OTU ID`), funs(./sum(.)*100)) %>%
  gather(sample_id, abundance, `7.DeMMO1.Steri.050917`:`18.800.DitchFluid.041818`) %>% 
  spread(key = `#OTU ID`,value = 'abundance') %>%
  right_join(metadata %>% select(sample_id)) %>%
  column_to_rownames("sample_id")
  

NMDS_ord <- otu_norm %>%
  metaMDS(k=2)

#pull out ordination and vector coordinates for plotting
NMDS_coords <- NMDS_ord[["points"]] %>%
  as_tibble(rownames = "sample_id") %>%
  left_join(metadata)

#make shape dictionary for ploting 
shape_dict <- c(0, 15, 15, 1, 19, 19, 2, 17, 17, 5, 5, 5, 5)
names(shape_dict) <- c("D1.fluid", "D1.inert.control", "D1.mineral", "D3.fluid", "D3.inert.control", "D3.mineral", "D6.fluid", "D6.inert.control", "D6.mineral","D3.cont.control", "4800.cont.control", "800.cont.control", "4100L.fluid")


#NMDS plot with controls 
NMDS_plot <- NMDS_coords %>%
  ggplot(aes(MDS1, MDS2)) +
  geom_point(size=2, alpha=0.8, aes(shape=site.type, color=site.type, label = sample_id)) +
  ggplot2::scale_shape_manual(values=shape_dict) +
  stat_ellipse(data=NMDS_coords, aes(color=site)) +
  #geom_text(aes(label=paste0(site, '.', date),hjust = 1, vjust = 1),  size=2, color="black") +  stat_ellipse(data=NMDS_coords, aes(color=site)) +
  theme(legend.key.size = unit(.5, "cm"))


#visualize interactive plot
plotly::ggplotly(NMDS_plot)




# NMDS plot without controls  ---------------------------------------------


otu_norm <- otu_table %>%
  select(-taxonomy) %>%
  mutate_at(vars(-`#OTU ID`), funs(./sum(.)*100)) %>%
  gather(sample_id, abundance, `7.DeMMO1.Steri.050917`:`18.800.DitchFluid.041818`) %>% 
  spread(key = `#OTU ID`,value = 'abundance') %>%
  right_join(metadata %>% select(sample_id, site)) %>%
  filter(!site == "ambient.control") %>%
  select(-site) %>%
  column_to_rownames("sample_id")

NMDS_ord <- otu_norm %>%
  metaMDS(k=2)

#pull out ordination and vector coordinates for plotting
NMDS_coords <- NMDS_ord[["points"]] %>%
  as_tibble(rownames = "sample_id") %>%
  left_join(metadata)


#NMDS plot with no controls 
NMDS_plot <- NMDS_coords %>%
  ggplot(aes(MDS1, MDS2)) +
  geom_point(size=2, alpha=0.8, aes(shape=site.type, color=site.type, label = sample_id)) +
  ggplot2::scale_shape_manual(values=shape_dict) +
  stat_ellipse(data=NMDS_coords, aes(color=site)) +
  #geom_text(aes(label=paste0(site, '.', date),hjust = 1, vjust = 1),  size=2, color="black") +  stat_ellipse(data=NMDS_coords, aes(color=site)) +
  theme(legend.key.size = unit(.5, "cm"))


#visualize interactive plot
plotly::ggplotly(NMDS_plot)

