library(RMySQL)
library(tidyverse)
library(hiAnnotator)
#library(microbenchmark)
library(doParallel)
library(furrr)

dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples <- dbGetQuery(dbConn,'select * from gtsp where Trial="CART19_ALL" or Trial="CART19_CLL"')                      



if( !file.exists('intSites_full_ALL_CLL.rds') ) {
  intSites <- getDBgenomicFragments(samples$SpecimenAccNum, 
                                    'specimen_management', 'intsites_miseq') %>%
    GenomicRanges::as.data.frame() %>% filter(refGenome == 'hg38') %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE) %>%
    stdIntSiteFragments(CPUs = numCores ) %>%
    collapseReplicatesCalcAbunds() %>%
    annotateIntSites(CPUs = numCores)
  saveRDS(intSites, 'intSites_full_ALL_CLL.rds')
} else {
  intSites <- readRDS('intSites_full_ALL_CLL.rds')
}

# test_sample <- intSites %>% as.data.frame() %>%
#   filter(GTSP=='GTSP0567') %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)


epigenetic_features_path <- "epi_rds" 
epi_files <- list.files(epigenetic_features_path)
names(epi_files) <- epi_files %>% str_remove(., ".rds")

plan(multisession, workers = 6)
kk <- future_imap(epi_files, function(x,name){
  epi_curr <- readRDS(file.path("epi_rds", x))
  p_kk <- getFeatureCounts(intSites, epi_curr, name) %>%
    as.data.frame() %>%
    select(!colnames(as.data.frame(intSites)))
  print(paste0('done with',name))
  return(p_kk)
})
epi_field <- Reduce(cbind,kk)
plan(sequential)
final_data <- cbind(test_sample %>% as.data.frame() ,epi_field)
saveRDS(final_data,'intSite_plus_epi_data.rds')