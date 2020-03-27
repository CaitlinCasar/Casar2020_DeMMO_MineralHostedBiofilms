#load dependencies 
pacman::p_load(tidyverse, readr, plyr, plotly, lubridate)

#load metadata
metadata <- read_csv("../orig_data/metadata.csv") 

# load collated alpha div data
data_path <- "../orig_data/collated_alpha"   # path to the data
files <- dir(data_path, pattern = "*.csv") # get file names
files <- paste(data_path, '/', files, sep="")

alpha_div = tibble::tibble(File = files) %>%
  tidyr::extract(File, "method", "(?<=/orig_data/collated_alpha/)(.*)(?=[.]csv)", remove = FALSE) %>%
  mutate(Data = lapply(File, readr::read_csv)) %>%
  tidyr::unnest(Data) %>%
  select(-File, -X1, -iteration) %>%
  gather(sample_id, value, `NU.39.DeMMO6.5mmPyrex.Sep2016`:`23.DeMMO1.D.top.041818`) %>%
  mutate(value = as.numeric(value)) %>%
  group_by(method, `sequences per sample`, sample_id) %>%
  summarise(value = mean(value)) %>%
  inner_join(metadata) 

#plot rarefaction curves 
alpha_plot <- alpha_div %>%
  ggplot(aes(`sequences per sample`, value, group=sample_id, color=site.type)) +
  geom_line() + 
  scale_colour_discrete(guide = 'none') +
  theme(legend.position="bottom") +
  theme_grey() +
  facet_wrap( ~ method, ncol=2, scales = "free_y") +
  guides(colour = guide_legend(title = "Site")) 


rarefied_alpha_plot <- alpha_div %>%
  filter(`sequences per sample` == 10000) %>%
  ggplot(aes(site, value, fill=site.type)) + 
  geom_violin() +
  stat_summary(fun.y=mean, geom="point", size=1, color="white") +
  coord_flip() +
  theme(legend.position = "none", 
        axis.title = element_blank()) +
  facet_wrap( ~ method, ncol=2, scales = "free_x")
