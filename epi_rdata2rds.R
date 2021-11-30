## Load epigenetic features
library(tidyverse)
mainDir <- getwd()
subDir <- 'epi_rds'
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)

epigenetic_features_path <- "/home/ubuntu/data/epidata/hg38/" 
epigenetic_features_files <- list.files(epigenetic_features_path) %>%
  grep("HeLa", ., invert = TRUE, value = TRUE) %>%
  grep("CD133", ., invert = TRUE, value = TRUE) %>%
  grep("HEK293T", ., invert = TRUE, value = TRUE) %>%
  grep("Jurkat", ., invert = TRUE, value = TRUE) %>%
  file.path(epigenetic_features_path, .)


epi_names <- str_split(epigenetic_features_files, "/")
epi_names <- sapply(epi_names, tail, n = 1L)
epi_names <- str_remove(epi_names, ".RData") %>% str_replace_all("-", "_")

#path <- '/home/ubuntu/data/epidata/hg38//Act-CD4-HDAC6.RData'
lapply(epigenetic_features_files, function(path){
  load(path)
  name <- str_split(path, "/") %>% sapply(., tail, n = 1L) %>% str_remove(., ".RData") %>% str_replace_all("-", "_")
  saveRDS(epigenData,file.path(mainDir, subDir,paste0(name,'.rds')))
})


