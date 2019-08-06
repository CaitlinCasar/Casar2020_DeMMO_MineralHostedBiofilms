
library(demmogorgon)

#import OTU table binned at family level, normalize to relative abundances 
data <- import_otu_table(add_metadata = FALSE, format_date = FALSE, normalize=TRUE)
#remove samples that don't match metadata identifiers 
data <- create_metadata(data, date = FALSE, return_clean_data = TRUE)
#import metadata 
metadata <- create_metadata(data, date = FALSE)

taxa_dict_data <- read.csv("taxa_dict.csv", header = TRUE, stringsAsFactors = FALSE)
gtdbk_taxa_dict_data <- read.csv("Unassigned_assignments.csv", header = TRUE, stringsAsFactors = FALSE)

assign_taxa <- function(dict){
  if(dict[1] %in% gtdbk_taxa_dict_data$Query_Seq_ID){
    #print(paste0(dict[1], ',',match(dict[1], gtdbk_taxa_dict_data$Query_Seq_ID), '.', gtdbk_taxa_dict_data[match(dict[1], gtdbk_taxa_dict_data$Query_Seq_ID),4] ))
    return(gtdbk_taxa_dict_data[match(dict[1], gtdbk_taxa_dict_data$Query_Seq_ID),4])
    
  }else{
    return(dict[2])
  }
}

taxa_dict_data$taxonomy <- apply(taxa_dict_data,1, assign_taxa)

taxa_dict <- taxa_dict_data$taxonomy
names(taxa_dict) <- taxa_dict_data$OTU.ID

taxa_dict <- taxa_dict_data$taxonomy
names(taxa_dict) <- taxa_dict_data$OTU.ID


color_dict <- create_color_dict(data, metadata)

dendro_plot <- plot_communities(data=data, metadata=metadata, dendro_bar = TRUE, legend = TRUE, threshold_abundance = 2)

site_simper <- simper_test(data=data, metadata=metadata, comparison_index = 10, num_perm = 1000, p_threshold = 0.05, sim_plot=TRUE, phylum = FALSE, taxa_level="family")

###run simper using Vegan
# run_simper <- function(meta){
#   data <- data[!is.na(match(rownames(data), meta[,1])),]
#   data = data[ rowSums(data[1:ncol(data)])!=0, ]
#   simper_test_data <- vegan::simper(data, meta[,7], permutations = 1000)
#   test <- as.data.frame(do.call(cbind, purrr::flatten(simper_test_data)))
#   test <- sapply(unique(names(test)), function(x) unname(unlist(test[,names(test)==x])))
#   test <- cbind(comparisons=names(rep(simper_test_data, each=ncol(data))), test)
#   test <- as.data.frame(test, stringsAsFactors = FALSE) %>% dplyr::filter(p<=0.05)
# }
# D1_simper <- run_simper(D1_samples)
# D3_simper <- run_simper(D3_samples)
# D6_simper <- run_simper(D6_samples)
# test <- rbind(D1_simper, D3_simper, D6_simper)
# 
# level1 <- "D_1"
# level2 <- "D_5"
# taxa_level <- paste0(".*", level1, "__\\s*|; ", level2, "__.*")
# 
# 
# test$species <- plyr::mapvalues(test$species, from=names(taxa_dict), to=taxa_dict)
# test$taxa <- gsub(taxa_level, "", test$species)
# test$taxa <- gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__", "", test$taxa)
# otu_data <- test %>%
#   dplyr::group_by_(.dots = list(colnames(test)[1], colnames(test[ncol(test)]))) %>%
#   dplyr::summarize(avg_dissim=sum(as.numeric(average)))
# otu_data$site <- gsub( "[.].*$", "", otu_data$comparisons)
# otu_data$comparisons <- gsub("D1.|D3.|D6.", "", otu_data$comparisons)
# 
# test$taxa <- gsub("; uncultured bacterium|; uncultured archaeon|; uncultured", "", test$taxa)
# 
# otu_data_reduced <- test
# otu_data_reduced$ava <- as.numeric(otu_data_reduced$ava)
# otu_data_reduced$avb <- as.numeric(otu_data_reduced$avb)
# otu_data_reduced$ab_diff <- otu_data_reduced$avb-otu_data_reduced$ava
# otu_data_reduced <- otu_data_reduced %>%
#   dplyr::group_by(comparisons, taxa) %>%
#   dplyr::summarize(avg_abdiff=mean(ab_diff))
# otu_data_reduced$site <- gsub( "[.].*$", "", otu_data_reduced$comparisons)
# otu_data_reduced$comparisons <- gsub("D1.|D3.|D6.", "", otu_data_reduced$comparisons)
# 
# 
# test_plot <- ggplot2::ggplot(otu_data_reduced , ggplot2::aes(fill=taxa, y=avg_abdiff, x=comparisons)) +
#   ggplot2::theme_gray() +
#   ggplot2::geom_bar(stat='identity') +
#   ggplot2::geom_hline(yintercept=0, linetype="dotted", color="black", size=0.5) +
#   ggplot2::coord_flip() +
#   ggplot2::theme(legend.position = "none") +
#   ggplot2::scale_fill_manual(values=fam_palette) + 
#   ggplot2::facet_grid(cols = dplyr::vars(site))
# 
# #install scripts for SIMPER and kruskal-wallis tests
# source_url("https://raw.githubusercontent.com/kdillmcfarland/workshops_UW_Madison/master/Microbiota_analysis_R/Steinberger_scripts/simper_pretty.R")
# source_url("https://raw.githubusercontent.com/kdillmcfarland/workshops_UW_Madison/master/Microbiota_analysis_R/Steinberger_scripts/R_krusk.R")
# #https://rstudio-pubs-static.s3.amazonaws.com/268156_d3ea37937f4f4469839ab6fa2c483842.html
# #SIMPER
# 
# site_krusk <- read.csv("experiment_krusk_simper.csv", header=TRUE, row.names=1)
# site_krusk <- site_krusk %>% dplyr::filter(fdr_krusk_p.val <= 0.05)
# site_krusk$avg_abdiff <- paste0(site_krusk$Left.mean.abund - site_krusk$Right.mean.abund)
# 
# data_no_control <- as.data.frame(t(data[1:73,]))
# 
# metadata_sites <- metadata %>% dplyr::filter(!site.type %in% "ambient.control")
# metadata_sites$site.type <- as.factor(metadata_sites$site.type)
# rownames(metadata_sites) <- metadata_sites$sample_ids
# simper.pretty(data_no_control, metadata_sites, c('site.type'), perc_cutoff=1, low_cutoff = 'n', output_name = 'simper.site.type.otu10000')
# 
# simper.results = data.frame(read.csv("simper.site.type.otu10000_clean_simper.csv"))
# kruskal.pretty(data_no_control, metadata_sites, simper.results, c('site.type'), 'site.type.krusk')
# 
# #Import
# KW.results = data.frame(read.csv("site.type.krusk_krusk_simper.csv"))
# #Remove non-significant
# KW.results.signif = KW.results[KW.results$fdr_krusk_p.val < 0.05,]
# #Order by OTU#
# KW.results.signif = KW.results.signif[with(KW.results.signif, order(OTU)),]
# head(KW.results.signif)
# KW.results.signif <- na.omit(KW.results.signif)
# KW.results.signif$Taxonomy <- plyr::mapvalues(KW.results.signif$OTU, from=names(taxa_dict), to=taxa_dict)
# 
# KW.results.signif$family <- gsub(taxa_level, "", KW.results.signif$Taxonomy)
# KW.results.signif$family <- gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__", "", KW.results.signif$family)
# KW.results.signif$family <- gsub("; uncultured bacterium|; uncultured archaeon|; uncultured", "", KW.results.signif$family)
# 
# kruskal_data <- KW.results.signif %>%
#   dplyr::group_by(Comparison, family)
# kruskal_data <- kruskal_data %>% 
#   dplyr::filter(Comparison %in% c("D1.fluid_D1.mineral", "D1.fluid_D1.inert.control","D1.inert.control_D1.mineral",
#                                   "D3.fluid_D3.mineral","D3.inert.control_D3.mineral","D3.fluid_D3.inert.control",
#                                   "D6.fluid_D6.mineral", "D6.fluid_D6.inert.control","D6.inert.control_D6.mineral"))
#  
# kruskal_data$abdiff <-  paste0(kruskal_data$Right.mean.abund - kruskal_data$Left.mean.abund)         
# kruskal_data <- kruskal_data %>%
#   dplyr::group_by(Comparison, family) %>% 
#   dplyr::summarize(avg_abdiff=mean(as.numeric(abdiff)))
# kruskal_data$site <- gsub( "[.].*$", "", kruskal_data$Comparison)
# kruskal_data$Comparison <- gsub("D1.|D3.|D6.", "", kruskal_data$Comparison)
# kruskal_data$site.fam <- paste0(kruskal_data$site,'.',kruskal_data$Comparison, '.', kruskal_data$family)
# 
# 
# simp_data <- KW.results.signif %>%
#   dplyr::filter(Comparison %in% c("D1.fluid_D1.mineral", "D1.fluid_D1.inert.control","D1.inert.control_D1.mineral",
#                                   "D3.fluid_D3.mineral","D3.inert.control_D3.mineral","D3.fluid_D3.inert.control",
#                                   "D6.fluid_D6.mineral", "D6.fluid_D6.inert.control","D6.inert.control_D6.mineral")) %>%
#   dplyr::group_by(Comparison, family) %>% 
#   dplyr::summarize(sum_simp=sum(as.numeric(SIMPER)))
# simp_data$site <- gsub( "[.].*$", "",simp_data$Comparison)
# simp_data$Comparison <- gsub("D1.|D3.|D6.", "", simp_data$Comparison)
# simp_data$site.fam <- paste0(simp_data$site,'.', simp_data$Comparison, '.', simp_data$family)
# 
# 
# krusk_simp_summary <- merge(kruskal_data, simp_data, by="site.fam")
# 
# 
# krusk_plot <- ggplot2::ggplot(krusk_simp_summary, ggplot2::aes(fill=family.x, y=avg_abdiff, x=Comparison.x)) +
#   ggplot2::theme_gray() +
#   ggplot2::geom_bar(stat='identity') +
#   ggplot2::geom_hline(yintercept=0, linetype="dotted", color="black", size=0.5) + 
#   ggplot2::scale_fill_manual(values=fam_palette) +
#   ggplot2::coord_flip() +
#   ggplot2::scale_y_continuous(limits = c(-0.3, 0.3)) +
#   ggplot2::theme(legend.position = "none") +
#   ggplot2::facet_grid(cols = dplyr::vars(site.x))
# 
# plotly::ggplotly(krusk_plot)
# 
# n <- length(unique(c(otu_data_reduced$taxa, krusk_simp_summary$family.x)))
# fam_palette <- randomcoloR::distinctColorPalette(n)
# names(fam_palette) <- unique(c(otu_data_reduced$taxa, krusk_simp_summary$family.x))
# 
# cowplot::plot_grid(test_plot, krusk_plot, ncol = 1, align = 'v', axis = 'b')

adonis2(data ~ site*date, data = metadata, permutations = 1000)
#define shapes palette for nmds plot
shape_dict <- c(0, 15, 15, 1, 19, 19, 2, 17, 17, 5, 5, 5, 5)
names(shape_dict) <- c("D1.fluid", "D1.inert.control", "D1.mineral", "D3.fluid", "D3.inert.control", "D3.mineral", "D6.fluid", "D6.inert.control", "D6.mineral","D3.cont.control", "4800.cont.control", "800.cont.control", "4100L.fluid")

### raw OTU table data 

#convert OTU ID's to taxonomy
taxa_dict_data <- read.csv("taxa_dict.csv", header = TRUE, stringsAsFactors = FALSE)
taxa_dict <- taxa_dict_data$taxonomy
names(taxa_dict) <- taxa_dict_data$OTU.ID

test <- plyr::mapvalues(colnames(data), from=names(taxa_dict), to=taxa_dict)




# editing section ---------------------------------------------------------
###bar plots clustered at OTU level, display family level
`%>%` <- magrittr::`%>%`
barplot_samples <- c("12.D1.fluid.041818", "26.D1.pyrolusite.041818", "22.D1.pyrite.041818","23.D1.hematite.041818", "24.D1.magnetite.041818",
  "27.D1.siderite.041818", "34.D1.sand.041818", "45.D1.calcite.041818","46.D1.gypsum.041818", "47.D1.muscovite.041818",
  "14.D3.fluid.041718", "51.D3.pyrolusite.041718","27.D3.siderite.051117", "53.D3.pyrite.041718", "54.D3.hematite.041718",
  "55.D3.magnetite.041718", "39.D3.calcite.041718", "40.D3.gypsum.041718","41.D3.muscovite.041718", "56.D3.sand.041718",
  "12.D6.fluid.051017", "13.D6.pyrolusite.051017","15.D6.siderite.051017", "17.D6.pyrite.051017" ,"19.D6.hematite.051017",
  "21.D6.magnetite.051017", "24.D6.sand.051017" ) 
barplot_samples <- plyr::mapvalues(barplot_samples, from=metadata$clean_id, to=metadata$sample_ids)
barplot_samples <- which(rownames(data) %in% barplot_samples)
barplot_data <- data[barplot_samples,]
barplot_data <- data

dendrogram_data <- as.dendrogram(hclust(ecodist::bcdist(barplot_data)))
dendro.plot <- ggplot2::ggplot(dendextend::as.ggdend(dendrogram_data), horiz = T) + ggplot2::theme_gray()

barplot_data <- merge(metadata, barplot_data, by.x = colnames(metadata[1]), by.y = "row.names")
barplot_data <- tidyr::gather(barplot_data, OTU, Abundance, colnames(barplot_data[ncol(metadata)+1]):colnames(barplot_data[ncol(barplot_data)]), factor_key=TRUE)
barplot_data$taxonomy <- plyr::mapvalues(barplot_data$OTU, from=names(taxa_dict), to=taxa_dict)
barplot_data$taxonomy <- gsub("Gammaproteobacteria; D_3__Betaproteobacteriales", "Betaproteobacteria; D_3__Betaproteobacteriales", barplot_data$taxonomy)

level1 <- "D_1"
level2 <- "D_2"
level3 <- "D_3"
level5 <- "D_5"
phylum_level <- paste0(".*", level1, "__\\s*|; ", level2, "__.*")
barplot_data$phylum <- gsub(phylum_level, "", barplot_data$taxonomy)
class_level <- paste0(".*", level2, "__\\s*|; ", level3, "__.*")
barplot_data$class <- gsub(class_level, "", barplot_data$taxonomy)
proteobacteria_fun <- function(data){
  if(data[14] == "Proteobacteria"){
    return(data[15])
  }else{
    return(data[14])
  }}
barplot_data$phylum <- apply(barplot_data, 1, proteobacteria_fun)
family_level <- paste0(".*", level1, "__\\s*|; ", level5, "__.*")
barplot_data$family <- gsub(family_level, "", barplot_data$taxonomy)
barplot_data$family <- gsub("; D_2__|; D_3__|; D_4__", ';', barplot_data$family)
barplot_data$family <- gsub(";uncultured.*|; Ambiguous.*|D_0__;", "", barplot_data$family)

##pull out Thermodesulfovibrionia and Desulfobulbaceae OTUs for BLAST
Thermodesulfovibrionia_Otus <- barplot_data %>% dplyr::filter(grepl("Thermodesulfovibrionia", taxonomy))
Thermodesulfovibrionia_Otus <- Thermodesulfovibrionia_Otus %>% dplyr::group_by(site, sample.type, OTU) %>% dplyr::summarize(abundance=sum(as.numeric(Abundance)))
Thermodesulfovibrionia_Otus <- Thermodesulfovibrionia_Otus %>% dplyr::filter(sample.type == "pyrolusite") %>% 
  dplyr::group_by(site, OTU) %>% dplyr::summarize(abundance=sum(as.numeric(abundance)))


Desulfobulbaceae_Otus <- barplot_data %>% dplyr::filter(grepl("Desulfobulbaceae", taxonomy))
Desulfobulbaceae_Otus <- Desulfobulbaceae_Otus %>% dplyr::group_by(site, sample.type, OTU) %>% dplyr::summarize(abundance=sum(as.numeric(Abundance)))
Desulfobulbaceae_Otus <- Desulfobulbaceae_Otus %>% dplyr::filter(sample.type == "pyrolusite") %>% 
  dplyr::group_by(site, OTU) %>% dplyr::summarize(abundance=sum(as.numeric(abundance)))



barplot_data <- barplot_data %>% dplyr::group_by(site, sample.type, family) %>%
  dplyr::summarize(abundance=sum(as.numeric(Abundance)))

# barplot_data <- barplot_data %>% dplyr::group_by(clean_id, site, sample.type, date, phylum) %>%
#   dplyr::summarize(abundance=sum(as.numeric(Abundance)))



taxa_to_filter <- c("Desulfomicrobiaceae|Unassigned|Gallionellaceae|Pseudomonadaceae|Acetothermiia|Sulfuricellaceae|Desulfobulbaceae|Thermodesulfovibrionia|Desulfobacteraceae|^Firmicutes;Clostridia$|Latescibacteraceae|Zixibacteria|Woesearchaeia|Lentimicrobiaceae|Rhodocyclaceae|Anaerolineaceae|Omnitrophicaeota|Alphaproteobacteria|Ignavibacteriales|GWF2-40-263|Desulfarculaceae|D8A-2|Syntrophomonadaceae")
filtered_data <- barplot_data %>% dplyr::filter(grepl(taxa_to_filter,family))

filtered_data$family <- gsub(".*Desulfobulbaceae.*", "Deltaproteobacteria;Desulfobacterales;Desulfobulbaceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Thermodesulfovibrionia.*", "Nitrospirae;Thermodesulfovibrionia",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Gallionellaceae.*", "Betaproteobacteria;Betaproteobacteriales;Gallionellaceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Rhodocyclaceae.*", "Betaproteobacteria;Betaproteobacteriales;Rhodocyclaceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Sulfuricellaceae.*", "Betaproteobacteria;Betaproteobacteriales;Sulfuricellaceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Desulfomicrobiaceae.*", "Deltaproteobacteria;Desulfovibrionales;Desulfomicrobiaceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Pseudomonadaceae.*", "Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Acetothermiia.*", "Acetothermia;Acetothermiia",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Desulfobacteraceae.*", "Deltaproteobacteria;Desulfobacterales;Desulfobacteraceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*D8A-2.*", "Firmicutes;Clostridia;D8A-2",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Syntrophomonadaceae.*", "Firmicutes;Clostridia;Clostridiales;Syntrophomonadaceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Firmicutes;Clostridia", "Firmicutes;Clostridia",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Latescibacteraceae.*", "Latescibacteria;Latescibacteria;Latescibacterales;Latescibacteraceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Zixibacteria.*", "Zixibacteria",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Woesearchaeia.*", "Nanoarchaeaeota;Woesearchaeia",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Lentimicrobiaceae.*", "Bacteroidetes;Bacteroidia;Sphingobacteriales;Lentimicrobiaceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Anaerolineaceae.*", "Chloroflexi;Anaerolineae;Anaerolineales;Anaerolineaceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Omnitrophicaeota.*", "Omnitrophicaeota",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Alphaproteobacteria.*", "Alphaproteobacteria",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Ignavibacteriales.*", "Bacteroidetes;Ignavibacteria;Ignavibacteriales",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Uhrbacterales;GWF2-40-263.*", "Patescibacteria;ABY1;Uhrbacterales;GWF2-40-263",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Desulfarculaceae.*", "Deltaproteobacteria;Desulfarculales;Desulfarculaceae",filtered_data$family, perl=TRUE) 



filtered_data <- filtered_data %>% dplyr::group_by(site, sample.type, family) %>%
  dplyr::summarize(abundance=sum(as.numeric(abundance)))

barplot_data_sums <- filtered_data %>% dplyr::group_by(site, sample.type) %>%
  dplyr::summarize(abundance=sum(as.numeric(abundance)))

# filtered_data <-  barplot_data %>%
#    dplyr::filter(abundance > 5)

# barplot_data_sums <- filtered_data %>% dplyr::group_by(clean_id, site, sample.type, date, phylum) %>%
#   dplyr::summarize(abundance=sum(as.numeric(abundance)))

barplot_data_sums$abundance <- 100 - barplot_data_sums$abundance

barplot_data_sums <- tibble::add_column(barplot_data_sums, family="Less.abundant.taxa", .after = "sample.type")
filtered_data <- rbind(filtered_data, barplot_data_sums)


sample.type_order <- c("fluid","sand","pyrolusite", "pyrite","hematite","magnetite","siderite",
                       "gypsum","muscovite","calcite")
filtered_data$sample.type <- factor(filtered_data$sample.type, levels = sample.type_order)

# filtered_data <- filtered_data %>% dplyr::filter(!site %in% 'ambient.control')
# filtered_data <- filtered_data %>% 
#   dplyr::filter(site %in% c("25.D1.sand.Sep2016","20.D1.pyrolusite.041818", "22.D1.pyrite.041818",
#                             "23.D1.hematite.041818","24.D1.magnetite.041818", "27.D1.siderite.041818",
#                             ))
phylum_color_data <- read.csv("phylumcolors_all_LM_MRO.csv", header=TRUE, stringsAsFactors = FALSE)
phylum_color_data$hex.color <- gsub("\\<0\\>", "000000", phylum_color_data$hex.color)
phylum_color_data$hex.color <- paste0('#', phylum_color_data$hex.color)
phylum_color_dict <- phylum_color_data$hex.color
names(phylum_color_dict) <- phylum_color_data$full.name

bar_plot <- ggplot2::ggplot(filtered_data, ggplot2::aes(fill=phylum, y=abundance, x=date)) +
  #ggplot2::theme_gray() +
  ggplot2::geom_bar(stat='identity', position='fill') +
  #ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values=phylum_color_dict) +
  #ggplot2::theme_classic() +
  ggplot2::theme(
    legend.position = "none", 
    axis.title = ggplot2::element_blank()) +
  #ggplot2::guides(ncol=1, shape = ggplot2::guide_legend(override.aes = list(size = 0.5))) +
  #color = ggplot2::guide_legend(override.aes = list(size = 3)))
  #ggplot2::guides(col = ggplot2::guide_legend(ncol = 1)) +
  ggplot2::theme(legend.title = ggplot2::element_blank(), 
                 legend.text = ggplot2::element_text(size = 8)) 

bar_plot + ggplot2::facet_grid(rows = dplyr::vars(clean_id), switch = "y")


n <- length(unique(filtered_data$family))
palette <- randomcoloR::distinctColorPalette(n)
names(palette) <- unique(filtered_data$family)

family_dict <- read.csv("familycolors.csv", header = TRUE)
family_dict$color <- paste0('#', family_dict$color)
family_palette <- family_dict$color
names(family_palette) <- family_dict$Family

filtered_data$family <- factor(filtered_data$family, levels = rev(family_dict$Family))

bar_plot <- ggplot2::ggplot(filtered_data, ggplot2::aes(fill=family, y=abundance, x=sample.type)) +
  ggplot2::theme_gray() +
  ggplot2::geom_bar(stat='identity', position='fill') +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values=family_palette) +
  ggplot2::theme(
    #legend.position = "none", 
    axis.title = ggplot2::element_blank()) +
  ggplot2::theme(legend.title = ggplot2::element_blank(), 
                 legend.text = ggplot2::element_text(size = 5))

bar_plot + ggplot2::facet_grid(cols = dplyr::vars(site))



# 
# sample_order <- c("D1.fluid", "D1.inert.control", "D1.calcite","D1.muscovite","D1.gypsum",
#                   "D1.siderite","D1.magnetite","D1.hematite","D1.pyrite", "D1.pyrolusite",
#                   "D3.fluid", "D3.inert.control", "D3.calcite","D3.muscovite",  "D3.gypsum", 
#                   "D3.siderite",  "D3.magnetite" , "D3.hematite","D3.pyrite", "D3.pyrolusite",
#                   "D6.fluid",  "D6.sand", "D6.siderite", "D6.magnetite" ,
#                   "D6.hematite" ,   "D6.pyrite",  "D6.pyrolusite") 
# 
# filtered_data$site.experiment <- factor(filtered_data$site.experiment, levels = sample_order)



# 
# taxa_to_filter <- taxa_dict[grep(taxa_to_filter, taxa_dict)]
# 
# filtered_data <-  barplot_data[,which(colnames(barplot_data) %in% names(taxa_to_filter))]

#filtered_data[,1] <- factor(filtered_data[,1], levels = labels(dendrogram_data))

  



filtered_data$family <- gsub(".*Desulfobulbaceae.*", "Deltaproteobacteria.Desulfobacterales.Desulfobulbaceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Thermodesulfovibrionia.*", "Nitrospirae.Thermodesulfovibrionia",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Gallionellaceae.*", "Betaproteobacteria.Betaproteobacteriales.Gallionellaceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Rhodocyclaceae.*", "Betaproteobacteria.Betaproteobacteriales.Rhodocyclaceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Sulfuricellaceae.*", "Betaproteobacteria.Betaproteobacteriales.Sulfuricellaceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Desulfomicrobiaceae.*", "Deltaproteobacteria.Desulfovibrionales.Desulfomicrobiaceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Pseudomonadaceae.*", "Gammaproteobacteria.Pseudomonadales.Pseudomonadaceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Acetothermiia.*", "Acetothermia.Acetothermiia",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Desulfobacteraceae.*", "Deltaproteobacteria.Desulfobacterales.Desulfobacteraceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Clostridia.*", "Firmicutes.Clostridia",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Latescibacteraceae.*", "Latescibacteria.Latescibacteria.Latescibacterales.Latescibacteraceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Zixibacteria.*", "Zixibacteria",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Woesearchaeia.*", "Nanoarchaeaeota.Woesearchaeia",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Lentimicrobiaceae.*", "Bacteroidetes.Bacteroidia.Sphingobacteriales.Lentimicrobiaceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Anaerolineaceae.*", "Chloroflexi.Anaerolineae.Anaerolineales.Anaerolineaceae",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Omnitrophicaeota.*", "Omnitrophicaeota",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Alphaproteobacteria.*", "Alphaproteobacteria",filtered_data$family, perl=TRUE) 
filtered_data$family <- gsub(".*Ignavibacteriales.*", "Bacteroidetes.Ignavibacteria.Ignavibacteriales",filtered_data$family, perl=TRUE) 


taxa_order <- c("Deltaproteobacteria.Desulfobacterales.Desulfobulbaceae","Nitrospirae.Thermodesulfovibrionia","Unassigned", 
                "Betaproteobacteria.Betaproteobacteriales.Gallionellaceae",
                "Betaproteobacteria.Betaproteobacteriales.Rhodocyclaceae","Betaproteobacteria.Betaproteobacteriales.Sulfuricellaceae",
                "Deltaproteobacteria.Desulfovibrionales.Desulfomicrobiaceae", 
                "Gammaproteobacteria.Pseudomonadales.Pseudomonadaceae","Acetothermia.Acetothermiia",
                "Deltaproteobacteria.Desulfobacterales.Desulfobacteraceae",
                "Firmicutes.Clostridia","Latescibacteria.Latescibacteria.Latescibacterales.Latescibacteraceae",
                "Zixibacteria","Nanoarchaeaeota.Woesearchaeia","Bacteroidetes.Bacteroidia.Sphingobacteriales.Lentimicrobiaceae",
                "Chloroflexi.Anaerolineae.Anaerolineales.Anaerolineaceae",
                "Omnitrophicaeota","Alphaproteobacteria","Bacteroidetes.Ignavibacteria.Ignavibacteriales",
                "Less.Abundant.Taxa")
filtered_data$family <- factor(filtered_data$family, levels = taxa_order)




unique(filtered_data$family)[which(!unique(filtered_data$family) %in% names(family_palette))]

cowplot::plot_grid(dendro.plot,bar_plot, align = "h")

###NMDS plot with OTU data, no vectors 
NMDS <- vegan::metaMDS(data,k=2)
NMDS.frame = data.frame(MDS1 = NMDS$points[,1], MDS2 = NMDS$points[,2])
merged_NMDS <- merge(NMDS.frame, metadata, by.x = "row.names", by.y ="sample_ids", all=TRUE)
merged_NMDS <- merged_NMDS %>% dplyr::filter(!site.type %in% "ambient.control")

plot <- ggplot2::ggplot(merged_NMDS, ggplot2::aes(x=MDS2, y=MDS1)) + 
  ggplot2::geom_point(ggplot2::aes(shape=site.type, color=site.type),size=2, alpha=0.8) + 
  ggplot2::scale_shape_manual(values=shape_dict) +
  #theme(legend.position="bottom", legend.box = "horizontal") +
  ggplot2::stat_ellipse(data=merged_NMDS, ggplot2::aes(color=site.type)) +
  ggplot2::stat_ellipse(data=merged_NMDS, ggplot2::aes(color=site)) +
  ggplot2::xlim(-4.5, 4.5) +
  ggplot2::ylim(-4.5, 4.5) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::theme_grey()  

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



# alpha_div ---------------------------------------------------------------
data_path <- "collated_alpha"   # path to the data
files <- dir(data_path, pattern = "*.csv") # get file names
files <- paste(data_path, '/', files, sep="")

alpha_div = tibble::tibble(File = files) %>%
  tidyr::extract(File, "method", "([^.]+)", remove = FALSE) %>%
  dplyr::mutate(Data = lapply(File, readr::read_csv)) %>%
  tidyr::unnest(Data) %>%
  dplyr::select(-File)

format_data <- function(dataframe) {
  stringr::str_replace_all(dataframe, c("collated_alpha/" = "", "n/a" = NA))
}

alpha_div <- as.data.frame(lapply(alpha_div,format_data), stringsAsFactors=FALSE)
alpha_div[3:ncol(alpha_div)] <-lapply(alpha_div[3:ncol(alpha_div)], as.numeric)


#aggregate by rarefaction id
alpha_div <- alpha_div %>% 
  dplyr::group_by(as.character(method), as.character(sequences.per.sample)) %>% 
  dplyr::summarize_at(dplyr::vars(-method, -X1, -sequences.per.sample,-iteration), dplyr::funs(mean))

#gather rarefaction data into long format for plotting
alpha_div <- tidyr::gather(alpha_div, sample_id, abundance,colnames(alpha_div[3]):colnames(alpha_div[ncol(alpha_div)]), factor_key=TRUE)
colnames(alpha_div)[1:2] <- c("method","rarefaction_depth")
alpha_div$rarefaction_depth <- as.numeric(alpha_div$rarefaction_depth)
alpha_div$sample_id <- gsub("X", "", alpha_div$sample_id)
alpha_to_remove <- as.vector(unique(alpha_div$sample_id[(!alpha_div$sample_id %in% unlist(metadata[1]))]))
alpha_to_remove <- which(alpha_div$sample_id %in% alpha_to_remove)
alpha_div <- alpha_div[-alpha_to_remove,]
alpha_div$Site <- plyr::mapvalues(alpha_div$sample_id, from=metadata$sample_ids, to=metadata$site)
alpha_div$site.experiment <- plyr::mapvalues(alpha_div$sample_id, from=metadata$sample_ids, to=metadata$site.experiment)
alpha_div$experiment.type <- plyr::mapvalues(alpha_div$sample_id, from=metadata$sample_ids, to=metadata$experiment.type)
alpha_div$site.experiment <- gsub("sand", "inert.control", alpha_div$site.experiment)
alpha_div$site.type <- plyr::mapvalues(alpha_div$sample_id, from=metadata$sample_ids, to=metadata$site.type)

#plot rarefaction curves 
alpha_plot <- ggplot2::ggplot(alpha_div, ggplot2::aes(rarefaction_depth, abundance, group=sample_id, color=site.type)) +
  ggplot2::geom_line() + 
  #theme(legend.position = "none") +
  ggplot2::scale_colour_discrete(guide = 'none') +
  #theme(legend.position = c(.1, .84), legend.text=element_text(size=6), legend.title = element_text(size=8, face="bold")) +
  #theme(legend.key.size =  unit(0.1, "in")) +
  #scale_x_continuous(expand = c(0.15, 0)) +
  #geom_dl(aes(label = sample.id), method = list(dl.combine("first.points", "last.points"), cex = 0.8)) +
  ggplot2::theme_grey() +
  ggplot2::facet_wrap( ~ method, ncol=2, scales = "free_y") +
  ggplot2::guides(colour = ggplot2::guide_legend(title = "Site")) #+
  #theme(legend.position="none") 

rarefied_alpha_data <- alpha_div %>% dplyr::filter(rarefaction_depth == 10000) %>% 
  dplyr::filter(method %in% "PD_whole_tree") %>%
  dplyr::filter(!Site %in% "ambient.control") %>%
  #dplyr::filter(site.experiment %in% c("D1.inert.control", "D3.inert.control","D6.inert.control",
                                       #"D1.fluid", "D3.fluid", "D6.fluid", "D1.pyrolusite",
                                       #"D3.pyrolusite","D6.pyrolusite")) %>%
  dplyr::group_by(method, Site, site.experiment, experiment.type) %>%
  dplyr::summarise(avg.abundance = mean(abundance))

rarefied_alpha_data$experiment.type <- factor(rarefied_alpha_data$experiment.type, levels=c("fluid", "inert.control",
                                                                                            "pyrolusite", "siderite",
                                                                                            "pyrite", "magnetite", "hematite",
                                                                                            "gypsum", "muscovite", "calcite"))

rarefied.diversity.plot <- ggplot2::ggplot(rarefied_alpha_data, ggplot2::aes(x=experiment.type, y=avg.abundance))+#, shape=site.experiment)) + 
  #geom_violin() +
  ggplot2::geom_point(ggplot2::aes(color=experiment.type), size=4) +
  #geom_boxplot(width=0.1) +
  #stat_summary(fun.y=mean, geom="point", size=1, color="white") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1) +
  #geom_jitter(shape=16, position=position_jitter(0.2)) +
  ggplot2::coord_flip() +
  ggplot2::theme_grey() +
  ggplot2::theme(legend.position = "none", 
        axis.title = ggplot2::element_blank()) +
  ggplot2::facet_wrap( ~ Site, ncol=3) +
  ggplot2::facet_wrap( ~ Site, ncol=3, scales = "free_x")

  rarefied_alpha_data_otus <- alpha_div %>% dplyr::filter(rarefaction_depth == 10000) %>% 
    dplyr::filter(method %in% "observed_otus") %>%
    dplyr::filter(!Site %in% "ambient.control") %>%
    #dplyr::filter(site.experiment %in% c("D1.inert.control", "D3.inert.control","D6.inert.control",
    #"D1.fluid", "D3.fluid", "D6.fluid", "D1.pyrolusite",
    #"D3.pyrolusite","D6.pyrolusite")) %>%
    dplyr::group_by(method, Site, site.experiment, experiment.type) %>%
    dplyr::summarise(avg.abundance = mean(abundance))

pd_otus_alpha_div <- merge(rarefied_alpha_data, rarefied_alpha_data_otus, by="site.experiment")
pd_otus_alpha_div$experiment <- plyr::mapvalues(pd_otus_alpha_div$experiment.type.x, from = metadata$experiment.type, to=metadata$experiment)
pd_otus_alpha_div$site.type <- paste0(pd_otus_alpha_div$Site.x, '.', pd_otus_alpha_div$experiment)

otus_vs_pd_alpha_plot <- ggplot2::ggplot(pd_otus_alpha_div, ggplot2::aes(avg.abundance.x, avg.abundance.y, color=Site.x, shape=experiment)) +
  ggplot2::geom_point(size=4)

moser_data <- read.csv("moser_data.csv", header = TRUE)
moser_data <- tidyr::gather(moser_data, Site, Abundance, BH1:BF, factor_key=TRUE)
moser_data <- moser_data %>% dplyr::filter(Site %in% c("BH2", "SC"))
moser_data$Site <- factor(moser_data$Site, levels = c("SC", "BH2"))

moser_plot <- ggplot2::ggplot(moser_data, ggplot2::aes_string(fill="clone", y="Abundance", x="Site")) +
  ggplot2::theme_gray() +
  ggplot2::geom_bar(stat='identity', position='fill') +
  #ggplot2::coord_flip() + 
  ggplot2::theme(legend.position = "none")

