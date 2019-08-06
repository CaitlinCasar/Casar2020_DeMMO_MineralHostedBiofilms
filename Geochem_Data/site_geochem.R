geochem_data <- read.csv("site_geochem.csv")
avg_geochem_data <- geochem_data %>% dplyr::group_by(Site) %>%
  #dplyr::summarise_all(dplyr::funs(mean), na.rm=TRUE)
  dplyr::summarise_at(.vars = colnames(geochem_data[3:ncol(geochem_data)]),
             .funs = c(Mean="mean", Sd="sd"), na.rm=TRUE) %>%
  tidyr::gather(solute, mean, methane_Mean:DOC_Mean) %>%
  tidyr::gather(solute_sd, std_dev, methane_Sd:DOC_Sd) %>%
  dplyr::mutate(solute = gsub("_Mean", "", solute)) %>%
  dplyr::select(-solute_sd)
  
avg_geochem_data$solute <- factor(avg_geochem_data$solute, levels=c("hydrogen", "CO", "methane", "ammonia", "nitrate", "DOC", "ferrous_iron", "sulfide", "sulfate"))


geochem_plot <- ggplot2::ggplot(avg_geochem_data, ggplot2::aes(solute, mean, color=solute)) +
  ggplot2::geom_point(ggplot2::aes(color=solute),size=5) + 
  ggplot2::scale_y_log10() + 
  ggplot2::geom_errorbar(ggplot2::aes(ymin=mean-std_dev, ymax=mean+std_dev), width=0.2) + 
  ggplot2::coord_flip() +
  ggplot2::geom_line(ggplot2::aes(group=solute), color="gray", linetype="dotted") 
  #ggplot2::geom_path(ggplot2::aes(group=solute), color="gray", linetype="dotted") +
  #ggplot2::scale_y_reverse() +
  #ggplot2::facet_wrap( ~ Site, scales = "free_x") 
geochem_plot + ggplot2::facet_grid(cols = dplyr::vars(Site)) +
  ggplot2::theme(legend.position = "none", 
                 axis.title.y=ggplot2::element_blank()) +
  ggplot2::ylab("Concentration (mg/L)") +
  ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) 
