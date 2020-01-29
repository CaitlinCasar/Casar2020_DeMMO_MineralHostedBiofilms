#load dependencies 
pacman::p_load(tidyverse, readr, plyr, plotly, lubridate)

#import otu table, metadata, taxonomy
otu_table <- read_delim("../orig_data/DeMMO136_Dec2015toApril2018_noChimera_otuTable_withTaxa_d10000.txt", delim = '\t', comment = "# ")
metadata <- read_csv("../orig_data/metadata.csv") 

#create taxonomy table updated with GTDB 
gtdb_taxonomy <- read_csv("../orig_data/Unassigned_assignments.csv") %>%
  select(Query_Seq_ID, GTDBtk_taxonomy) %>%
  rename(`#OTU ID` = Query_Seq_ID) %>%
  rename(taxonomy = GTDBtk_taxonomy)


barplot_samples <- c("12.DeMMO1.steri.041818", 
                     "26.DeMMO1.SC1.top.041818", 
                     "22.DeMMO1.C.top.041818",
                     "23.DeMMO1.D.top.041818", 
                     "24.DeMMO1.E.top.041818",
                     "27.DeMMO1.SC2.top.041818", 
                     "34.DeMMO1.SC10.top.041818", 
                     "45.DeMMO1.7.top.041818",
                     "46.DeMMO1.8.top.041818", 
                     "47.DeMMO1.9.top.041818",
                     "14.DeMMO3.steri.041718", 
                     "51.DeMMO3.A.top.041718",
                     "27.DeMMO3.T8.top.051117", 
                     "53.DeMMO3.C.top.041718", 
                     "54.DeMMO3.D.top.041718",
                     "55.DeMMO3.E.top.041718", 
                     "39.DeMMO3.1.top.041718", 
                     "40.DeMMO3.2.top.041718",
                     "41.DeMMO3.3.top.041718", 
                     "56.DeMMO3.F.top.041718",
                     "12.DeMMO6.Steri#2.051017", 
                     "13.DeMMO6.T1.top.051017",
                     "15.DeMMO6.T2.top.051017", 
                     "17.DeMMO6.T3.top.051017",
                     "19.DeMMO6.T4.top.051017",
                     "21.DeMMO6.T5.top.051017", 
                     "24.DeMMO6.T6.bottom.051017") 

abundance_table <- otu_table %>%
  select(-taxonomy) %>%
  mutate_at(vars(-`#OTU ID`), funs(./sum(.)*100)) %>% #normalize to relative abundance 
  gather(sample_id, abundance, `7.DeMMO1.Steri.050917`:`18.800.DitchFluid.041818`) 

#double check for date duplicates
test <- metadata %>%
  filter(sample_id %in% barplot_samples) %>% #remove abient controls 
  mutate(date = paste(month(mdy(date), label=TRUE), year(mdy(date)), sep = '.')) %>%
  select(site, date, experiment.type) %>%
  filter(duplicated(.))

taxon_abundance <- function(level, name){
  otu_table %>%
    select(`#OTU ID`, taxonomy) %>%
    filter(!`#OTU ID` %in% gtdb_taxonomy$`#OTU ID`) %>%
    bind_rows(gtdb_taxonomy) %>%
    mutate(taxonomy = gsub("Gammaproteobacteria; D_3__Betaproteobacteriales", "Betaproteobacteria; D_3__Betaproteobacteriales", taxonomy),
           taxa = str_extract(taxonomy, level),
           taxa = if_else(is.na(taxa), taxonomy, taxa)) %>%
    right_join(abundance_table) %>%
    ungroup() %>%
    group_by(sample_id, taxa) %>%
    summarise(abundance = sum(abundance)) %>%
    left_join(metadata) %>% #add metadata
    filter(sample_id %in% barplot_samples) %>% #remove unwanted samples 
    group_by(site, experiment.type, taxa) %>% 
    summarise(abundance = sum(abundance))
}

family_level <- "(.*)(?=; D_5__)"
class_level <- "(.*)(?=; D_3__)"
phylum_level <- "(.*)(?=; D_2__)"

less_abundant_taxa <- taxon_abundance(family_level, "family") %>%
  group_by(taxa) %>%
  filter(max(abundance) < 5) %>% #filter for families that represent less than 5% 
  group_by(site, experiment.type) %>%
  summarise(abundance = sum(abundance)) %>%
  mutate(taxa = 'Less abundant taxa')


taxon_abundance_table <- taxon_abundance(family_level, "family") %>%
  group_by(taxa) %>%
  filter(max(abundance) >= 5) %>%
  bind_rows(less_abundant_taxa)


#create color dictionary for figure
family_color <- read_csv("../orig_data/familycolors.csv")

family_color_dict <- family_color$color
names(family_color_dict) <- family_color$family

#bar plot figure 
bar_plot <- taxon_abundance_table %>%
  ungroup() %>%
  mutate(taxonomy = str_remove_all(taxa, "D_0__| D_1__| D_2__| D_3__| D_4__")) %>%
  separate(taxonomy ,sep=';',c("domain", "phylum", "class", "order", "family")) %>%
  mutate(family = if_else(is.na(family) | str_detect(family, "uncultured"), order, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), class, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), phylum, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), domain, family),
         family = if_else(family == "GWF2-40-263", order, family),
         experiment.type = factor(experiment.type, levels = c("fluid", "inert.control", "pyrolusite", "pyrite","hematite","magnetite","siderite","gypsum","muscovite","calcite")),
         family = factor(family, levels = family_color$family)) %>% 
  group_by(site, experiment.type, family) %>%
  summarise(abundance = sum(abundance)) %>%
  ggplot(aes(fill=family, y=abundance, x=experiment.type)) +
  geom_bar(stat='identity', position='fill') +
  scale_fill_manual(values=family_color_dict) +
  coord_flip() + 
  theme(axis.title = element_blank(),
        legend.title = ggplot2::element_blank(), 
        legend.text = ggplot2::element_text(size = 8),
        legend.key.size = unit(0.5, "cm")) +
  facet_grid(cols = vars(site), switch = "y") + 
  guides(fill = guide_legend(ncol = 1))
