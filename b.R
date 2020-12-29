library(tidyverse)
library(RMySQL)
library(gt23)
library(GenomicRanges)

invisible(sapply(dbListConnections(MySQL()), dbDisconnect))

dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples_all <- dbGetQuery(dbConn, 'select * from gtsp')
samples <- dbGetQuery(dbConn, 'select * from gtsp where Trial in ("CART19_ALL", "CART19_CLL")')

CART19_samples <- dbGetQuery(dbConn, 'select * from gtsp where Trial = "CART19"')

dbConn  <- dbConnect(MySQL(), group='intsites_miseq')
intSitesamples <- unname(unlist(dbGetQuery(dbConn, 'select sampleName from samples where sampleName like "%GTSP%"')))
intSitesamples <- unique(gsub('\\-\\d+$', '', intSitesamples))

samples <- dplyr::filter(samples, SpecimenAccNum %in% intSitesamples) %>%
           dplyr:::select(SpecimenAccNum, Trial, Patient, CellType, Timepoint, VCN)

j <- unlist(lapply(readLines('JCI_sample_table'), function(x) unlist(stringr::str_extract_all(x, 'GTSP\\d+'))))

j <- j[! j %in% samples$SpecimenAccNum]
j <- dplyr:::select(subset(samples_all, SpecimenAccNum %in% j), SpecimenAccNum, Trial, Patient, CellType, Timepoint, VCN)

samples <- bind_rows(j, samples)

p <- read.table('JCI_patient_table', sep = '\t', header = TRUE, check.names = TRUE) %>%
     dplyr::select('Patient.ID', Disease, Response, Resp...Class, Clin...Trial)
names(p) <- c('Patient', 'Disease', 'Response', 'Response_class', 'Trial_ID')

samples <- left_join(samples, p, by = 'Patient')

v <- read.table('VCN_data', sep = '\t', header = TRUE)
samples <- left_join(samples, dplyr::select(v, SpecimenAccNum, VCN2), by = "SpecimenAccNum")
samples$VCN <- ifelse(is.na(samples$VCN2), samples$VCN, samples$VCN2)
samples$VCN2 <- NULL 
samples$VCN <- ifelse(samples$VCN == 0, NA, samples$VCN)

intSites <- getDBgenomicFragments(samples$SpecimenAccNum, 'specimen_management', 'intsites_miseq') %>%
            stdIntSiteFragments() %>%
            collapseReplicatesCalcAbunds()

sampleIntSites <- dplyr::group_by(data.frame(intSites), GTSP) %>%
                  dplyr::summarise(nSites = n_distinct(posid)) %>%
                  dplyr::ungroup()

samples <- left_join(samples, sampleIntSites, by = c("SpecimenAccNum" = "GTSP"))


# Relationship between VCN and number of sites recovered.
ggplot(subset(samples, ! is.na(VCN) & nSites > 0 & nSites <= 250), aes(nSites, VCN)) +
theme_bw() + geom_point() + geom_smooth(method = 'lm')

ggplot(subset(samples, ! is.na(VCN) & nSites > 0), aes(log10(nSites), log10(VCN))) + 
theme_bw() + geom_point() + geom_smooth(method = 'lm')


openxlsx::write.xlsx(samples, file = 'CART_samples_201106.xlsx')
