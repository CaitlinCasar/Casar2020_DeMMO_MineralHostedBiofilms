#load dependencies 
pacman::p_load(tidyverse, readr, plyr)

metadata <- read_csv("../orig_data/metadata.csv") 
biosample_accession <- tibble(site = unique(metadata$site), 
                              biosample_accession = c("SAMN12684770", "SAMN12684819", "SAMN12684862", "SAMN12768991"))


biosample_assignments <- read_delim("../orig_data/biosample_assignment.tsv", delim = '\t') %>%
  select(Sequence_ID) %>%
  mutate(sample_id = stringr::str_replace(Sequence_ID, "_.*", "")) %>%
  left_join(metadata) %>%
  left_join(biosample_accession) %>%
  select(Sequence_ID, biosample_accession) %>%
  write.table("../orig_data/biosample_assignment.tsv", sep = "\t", row.names=FALSE, quote=FALSE)

ncbi_chimeras <- read.table("../MicrobiomePub_16s_QIIME/ncbi_chimeric_seqs.txt",comment.char = "")
biosample_assignments <- biosample_assignments %>%
  dplyr::filter(!Sequence_ID %in% ncbi_chimeras$V1)

write.table(biosample_assignments, "biosample_assignments.tsv", sep = "\t", row.names=FALSE, quote=FALSE)



