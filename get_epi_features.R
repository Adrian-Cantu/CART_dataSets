library(RMySQL)
library(tidyverse)
library(hiAnnotator)
#library(microbenchmark)
#library(doParallel)
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
intSites <- intSites %>% as.data.frame() %>%
  mutate(GTPSposID=paste0(GTSP,posid  )) %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)

 test_sample <- intSites %>% as.data.frame() %>%
   filter(GTSP=='GTSP0567') %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)


epigenetic_features_path <- "epi_rds" 
epi_files <- list.files(epigenetic_features_path)
names(epi_files) <- epi_files %>% str_remove(., ".rds")

## tested

all_names <- unique(intSites$GTSP)
# all_names <- all_names[1:3]
# 
# plan(sequential)
# plan(multisession, workers = 20)
# full_table <- lapply(all_names,function(c_gtsp){
#   c_sample <- intSites %>% as.data.frame() %>%
#     filter(GTSP==c_gtsp) %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)
#   print(paste0('starting work on ',c_gtsp))
#   kk <- future_imap(epi_files, function(x,name){
#     print(paste0('--starting with',name))
#     epi_curr <- readRDS(file.path("epi_rds", x))
#     p_kk <- getFeatureCounts(c_sample, epi_curr, name) %>%
#       as.data.frame() %>%
#       select(!colnames(as.data.frame(c_sample)))
#     rm(epi_curr)
#     print(paste0('--done with',name))
#     return(p_kk)
#   })
# #  plan(sequential)
#   epi_field <- Reduce(cbind,kk)
#   epi_field$GTPSposID <- c_sample$GTPSposID
#   return(epi_field)
# })
# 
# epi_final_data <- Reduce(rbind,full_table)
# final_final_data <- left_join(intSites %>% as.data.frame(),epi_final_data,by='GTPSposID')


## to test

all_epi <- lapply(epi_files,function(x){readRDS(file.path("epi_rds", x))})
#kk_test <- getFeatureCounts(intSites, all_epi[[1]], 'test')

plan(sequential)
plan(multisession, workers = 30)
full_table2 <- lapply(all_names,function(c_gtsp){
  c_sample <- intSites %>% as.data.frame() %>%
    filter(GTSP==c_gtsp) %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)
  print(paste0('starting work on ',c_gtsp))
  kk <- future_imap(all_epi, function(x,name){
    print(paste0('- - - - starting with',name))
    #epi_curr <- readRDS(file.path("epi_rds", x))
    p_kk <- getFeatureCounts(c_sample, x, name) %>%
      as.data.frame() %>%
      select(!colnames(as.data.frame(c_sample)))
    #print(paste0('--done with',name))
    return(p_kk)
  })
  #  plan(sequential)
  epi_field <- Reduce(cbind,kk)
  epi_field$GTPSposID <- c_sample$GTPSposID
  return(epi_field)
})

epi_final_data <- Reduce(rbind,full_table2)
final_final_data <- left_join(intSites %>% as.data.frame(),epi_final_data,by='GTPSposID')
saveRDS(epi_final_data,file='epi_data.rds')

