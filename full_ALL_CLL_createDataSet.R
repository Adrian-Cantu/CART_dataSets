options(stringsAsFactors = FALSE)
library(dplyr)
library(gt23)
library(stringr)
library(hiAnnotator)
library(GCcontent)
library(BSgenome)
library(gintools)
library(RMySQL)
library(geneRxCluster)
library(ggplot2)
numCores <- 7
source('utils.R')
source('supporting_functions.R')
source('clusters.R')
utilsDir <- 'utils'


# # Create a list of all GTSPs that passed through the INSPIIRED pipeline.
# dbConn  <- dbConnect(MySQL(), group='intsites_miseq')
# intSitesamples <- unname(unlist(dbGetQuery(dbConn, 'select sampleName from samples where sampleName like "%GTSP%"')))
# #DBI::dbDisconnect(dbConn)
# 
# intSitesamples <- unique(gsub('\\-\\d+$', '', intSitesamples))
# 
# sampleData <- read.table('sampleData.tsv', sep = '\t', header = TRUE) # this seems to be a table that includes all ALL
# sampleData$INSPIIRED_processed <- sampleData$SpecimenAccNum %in% intSitesamples
# 
# samples  <- subset(sampleData, INSPIIRED_processed == TRUE)$SpecimenAccNum

# #-----
# dbConn  <- dbConnect(MySQL(), group='specimen_management')
# responses <- read.table('sampleUpdates/patientResponses.tsv', sep = '\t', header = TRUE, fill = TRUE, comment.char = '', check.names = FALSE)
# gtsp <- dbGetQuery(dbConn, "select * from gtsp")
# gtsp <- subset(gtsp, Patient %in% responses$PatietID)
# gtsp <- subset(gtsp, SpecimenAccNum %in% intSitesamples)
# table(gtsp$SpecimenAccNum %in% samples)
# #----

dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples <- dbGetQuery(dbConn, 
                      'select * from gtsp where Trial="CART19_ALL" or Trial="CART19_CLL"')                      



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




# Recreate columns used in published data.
intSites$specimen <- intSites$GTSP
intSites$celltype <- intSites$cellType
intSites$timePoint <- tolower(intSites$timePoint)
intSites$timepoint <- intSites$timePoint
intSites$in_gene <- ifelse(intSites$nearestFeatureDist == 0, intSites$nearestFeature, FALSE)
intSites$in_geneOrt <- ifelse(intSites$nearestFeatureDist == 0, intSites$nearestFeatureStrand, NA)
intSites$nearest_geneDist <- intSites$nearestFeatureDist
intSites$nearest_gene <- intSites$nearestFeature
intSites$nearest_geneOrt <- intSites$nearestFeatureStrand
intSites$gene_id_wo_annot <- intSites$nearestFeature


# Generate summary stats for each sample --------------------------------------


# Load genomic and epigenetic features 
genome_sequence <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
genome_sequence@user_seqnames <- genome_sequence@user_seqnames[genome_sequence@user_seqnames %in% paste0("chr", c(1:22, "X", "Y", "M"))]
genome_sequence@seqinfo <- genome_sequence@seqinfo[paste0("chr", c(1:22, "X", "Y", "M"))]


cond_uniq_sites <- intSites[
  seqnames(intSites) %in% paste0("chr", c(1:22, "X", "Y", "M"))]
seqlevels(cond_uniq_sites) <- paste0("chr", c(1:22, "X", "Y", "M"))
seqinfo(cond_uniq_sites) <- seqinfo(genome_sequence)

refGenes <- readRDS(file.path(utilsDir, "hg38.refSeq.rds"))
refGenes <- refGenes[seqnames(refGenes) %in% paste0("chr", c(1:22, "X", "Y", "M"))]

rm(intSites)

#trial <- "CART19"
#std_clin_trials <- c("959", "04409", "03712")

# Split unique sites to pre and post transduction sets
tdn_sites <- cond_uniq_sites[cond_uniq_sites$timePoint == "d0"]
### tdn_sites <- tdn_sites[tdn_sites$patient %in% std_clin_patients]
timePoint_sites <- cond_uniq_sites[cond_uniq_sites$timePoint != "d0"]
### timePoint_sites <- timePoint_sites[timePoint_sites$patient %in% std_clin_patients]


o <- createClusters()
enriched_tdn_clus <- o[[1]]
enriched_timePoint_clus <- o[[2]]
composite_clusters <- o[[3]]
clusters <- o[[4]]
cart19_clusters <- o[[5]]
red_clusters <- o[[6]]
red_df <- o[[7]]
rm(o)



#for (thing in ls()) { message(thing); print(object.size(get(thing)), units='auto') }

#objSizes()


CpG_data <- cpg <- getUCSCtable("cpgIslandExt", "CpG Islands", freeze = "hg38") %>% 
                   dplyr::filter(chrom %in% paste0("chr", c(1:22, "X", "Y", "M")))

#getUCSCtable("cpgIslandExt", "CpG Islands", freeze = "hg38")


CpG_islands <- GenomicRanges::GRanges(
  seqnames = CpG_data$chrom,
  ranges = IRanges::IRanges(
    start = CpG_data$chromStart, end = CpG_data$chromEnd
  ),
  strand = "*",
  seqinfo = GenomeInfoDb::seqinfo(genome_sequence)
)

mcols(CpG_islands) <- CpG_data


DNaseI_data <- getUCSCtable("wgEncodeRegDnaseClustered", "DNase Clusters", freeze = "hg38") %>%
  dplyr::filter(chrom %in% paste0("chr", c(1:22, "X", "Y", "M")))

DNaseI <- GenomicRanges::GRanges(
  seqnames = DNaseI_data$chrom,
  ranges = IRanges::IRanges(start = DNaseI_data$chromStart, end = DNaseI_data$chromEnd),
  strand = "*",
  seqinfo = GenomeInfoDb::seqinfo(genome_sequence))

mcols(DNaseI) <- DNaseI_data

## windows
window_size_refSeq <- c("10k"=1e4, "100k"=1e5, "1M"=1e6)
window_size_CpG_counts <- c("1k"=1e3, "10k"=1e4)
window_size_CpG_density <- c("10k"=1e4, "100k"=1e5, "1M"=1e6)
window_size_GC <- c("100"=100, "1k"=1000, "10k"=1e4, "100k"=1e5, "1M"=1e6)
window_size_DNaseI <- c("1k"=1e3, "10k"=1e4, "100k"=1e5, "1M"=1e6)
window_size_epi <- c("10k"=1e4)



## Load epigenetic features
epigenetic_features_path <- "/home/ubuntu/data/epidata/hg38/" 
epigenetic_features_files <- list.files(epigenetic_features_path) %>%
  grep("HeLa", ., invert = TRUE, value = TRUE) %>%
  grep("CD133", ., invert = TRUE, value = TRUE) %>%
  grep("HEK293T", ., invert = TRUE, value = TRUE) %>%
  grep("Jurkat", ., invert = TRUE, value = TRUE) %>%
  file.path(epigenetic_features_path, .)


epi_names <- str_split(epigenetic_features_files, "/")
epi_names <- sapply(epi_names, "[[", 6)
epi_names <- str_remove(epi_names, ".RData") %>% str_replace_all("-", "_")

epi_env <- new.env()
epi_env$epi_features <- list()

null <- lapply(epigenetic_features_files, function(path){
  load(path)
  epi_env$epi_features <- c(epi_env$epi_features, list(epigenData))
})

epi_features <- epi_env$epi_features
names(epi_features) <- epi_names


oncoGenesData <- read.delim(file.path(utilsDir, "allOnco.human.v3.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

oncoGenes <- unique(oncoGenesData[,"symbol"])
nonOncoGenes <- unique(refGenes$name2[!refGenes$name2 %in% oncoGenes])

badActors <- read.delim(file.path(utilsDir, "humanLymph.v1.list"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1]


# Import HGNC reference data for annotation and consistency ----
hgnc_complete <- data.table::fread(
  paste0("zcat ", file.path(utilsDir, "hgnc_complete_set.180207.txt.gz")),
  sep = "\t", header = TRUE,
  select = c(
    "HGNC ID", "Approved Symbol", "Approved Name", "Locus Group", "Locus Type",
    "Synonyms", "Previous Symbols", "Entrez Gene ID", "Ensembl Gene ID",
    "RefSeq (supplied by NCBI)", "UniProt ID (supplied by UniProt)"),
  data.table = FALSE
)

names(hgnc_complete) <- c(
  "hgnc_id", "symbol", "name", "locus_group", "locus_type", "alias_symbol",
  "prev_symbol", "entrez_id", "ensembl_gene_id", "refseq_accession",
  "uniprot_ids"
)

hgnc_complete <- dplyr::filter(hgnc_complete, !grepl("withdrawn", symbol)) %>%
  dplyr::mutate(
    kegg_id = paste0("hsa:", entrez_id),
    entrez_id = paste0(entrez_id, ":EZID"),
    alias_symbol = gsub(", ", "|", alias_symbol),
    prev_symbol = gsub(", ", "|", prev_symbol),
    extended_alias = paste0(
      alias_symbol, "|", prev_symbol, "|", ensembl_gene_id, "|",
      refseq_accession, "|", uniprot_ids)
  )


rm(null)


p <- cond_uniq_sites$posid
g <- cond_uniq_sites$GTSP
w <- cond_uniq_sites$nearestFeature
o <- cond_uniq_sites$inFeatureSameOrt
a <- cond_uniq_sites$estAbund
d <- cond_uniq_sites$nearestFeatureDist

mcols(cond_uniq_sites)      <- NULL
cond_uniq_sites$posid       <- p
cond_uniq_sites$specimen    <- g
cond_uniq_sites$nearest_geneDist <- d
cond_uniq_sites$within_gene <- ifelse(cond_uniq_sites$nearest_geneDist == 0, TRUE, FALSE)
cond_uniq_sites$same_ort    <- o
cond_uniq_sites$estAbund    <- a
cond_uniq_sites$in_gene     <- ifelse(cond_uniq_sites$nearest_geneDist == 0, w, FALSE)



#o <- dplyr::select(mcols(cond_uniq_sites), patient, GTSP, cellType, timePoint, estAbund, relAbund)
#mcols(cond_uniq_sites) <- NULL

cond_uniq_list <- split(
  cond_uniq_sites, 
  ceiling(seq_along(cond_uniq_sites) / 
            (length(cond_uniq_sites)/length(unique(cond_uniq_sites$specimen))))
)


buster <- parallel::makeCluster(numCores) 

parallel::clusterExport(
  buster, 
  varlist = c(
    "from_counts_to_density", "refGenes", "CpG_islands", "DNaseI", 
    "window_size_refSeq", "window_size_GC", "window_size_CpG_counts", 
    "window_size_CpG_density", "window_size_DNaseI", "genome_sequence"
  ), 
  envir = environment()
)
#length(unique(cond_uniq_list@unlistData@elementMetadata@listData$specimen))
cond_uniq_annot_gen <- unname(unlist(GRangesList(parallel::parLapply(
  buster, 
  cond_uniq_list,
  function(gr){
    
    library(magrittr)
    library(GenomicRanges)
    
    gr %>%
      hiAnnotator::getFeatureCounts(
        refGenes, "refSeq_counts", width = window_size_refSeq) %>%
      GCcontent::getGCpercentage(
        "GC", window_size_GC, genome_sequence) %>%
      hiAnnotator::getFeatureCounts(
        CpG_islands, "CpG_counts", width = window_size_CpG_counts) %>%
      hiAnnotator::getFeatureCounts(
        CpG_islands, "CpG_density", width = window_size_CpG_density) %>%
      from_counts_to_density(
        "CpG_density", window_size_CpG_density) %>%
      hiAnnotator::getFeatureCounts(
        DNaseI, "DNaseI_count", width = window_size_DNaseI)
    
  }
))))

parallel::stopCluster(buster)

buster <- parallel::makeCluster(numCores) 
# Memory requirement ~2.5 GB per core

parallel::clusterExport(
  buster, 
  varlist = c("cond_uniq_sites", "window_size_epi"), 
  envir = environment()
)

epi_annots <- dplyr::bind_cols(parallel::clusterMap(
  buster,
  function(epi, name){
    
    library(GenomicRanges)
    
    annot_gr <- hiAnnotator::getFeatureCounts(
      cond_uniq_sites, epi, name, width = window_size_epi
    )
    
    feat_cols <- names(mcols(annot_gr))[
      !names(mcols(annot_gr)) %in% names(mcols(cond_uniq_sites))
      ]
    
    as.data.frame(mcols(annot_gr)[, feat_cols, drop = FALSE])
    
  },
  epi = epi_features,
  name = names(epi_features)
))


parallel::stopCluster(buster)

rm(epi_features)

## Combine and summarise features
cond_uniq_annot <- cond_uniq_annot_gen
mcols(cond_uniq_annot) <- bind_cols(
  as.data.frame(mcols(cond_uniq_annot)), 
  epi_annots
)

if( !all(cond_uniq_annot$posid == cond_uniq_sites$posid) ){
  stop("Indexing error occured during parallel processing.\n")
}



gen_epi_stats <- as.data.frame(cond_uniq_annot, row.names = NULL) %>%
  dplyr::select(-seqnames, -start, -end, -width)

cond_uniq_df <- as.data.frame(cond_uniq_sites, row.names = NULL)

stats <- cond_uniq_df %>%
  dplyr::group_by(specimen) %>%
  
  dplyr::mutate(relAbund = estAbund / sum(estAbund)) %>%
  
  dplyr::summarise(
    "numUniqSites" = n(), 
    "maxRelAbund" = max(relAbund),
    "numClones" = sum(estAbund), 
    "ShannonIndex" = pop_calcs(estAbund, calc = "shannon"),
    "GiniIndex" = pop_calcs(estAbund, calc = "gini"),
    "Chao1" = calculateChao(as.character(Rle(posid, estAbund))),
    "UC50" = pop_calcs(estAbund, calc = "uc50"),
    "pctTxnUnit" = 100 * sum(within_gene, na.rm = TRUE)/n(),
    "pctSameOrt" = 100 * sum(same_ort, na.rm = TRUE) / 
      sum(within_gene, na.rm = TRUE),
    "pctNearTxnUn" = 100 * sum(
      in_gene == FALSE & abs(nearest_geneDist) <= 5000, na.rm = TRUE) / 
      sum(in_gene == FALSE, na.rm = TRUE),
    "pctInOnco" = 100 * sum(in_gene %in% oncoGenes, na.rm = TRUE) / n())

gen_epi_summary <- gen_epi_stats %>%
  dplyr::select(-posid) %>%
  dplyr::group_by(specimen) %>%
  dplyr::summarise_all(mean, na.rm = TRUE)

result <- left_join(stats, gen_epi_summary, by = 'specimen')
result <- select(result, -strand, -in_gene)

samples_ML <- samples %>% 
  select(c('SpecimenAccNum',"Trial","CellType","Patient","Timepoint",'VCN')) %>%
  rename(SpecimenAccNum='specimen') %>%
  mutate(Trial=sub('^CART19_','',Trial)) %>%
  rename(Trial='cancer')

results_full <- left_join(result,samples_ML,by = 'specimen')

saveRDS(results_full, file = 'result_full_ALL_CLL.rds')
# 
# r <- left_join(select(sampleData, SpecimenAccNum, Trial, Patient, CellType, Timepoint, INSPIIRED_processed), result, by = c('SpecimenAccNum' = 'specimen'))
# 
# 
# 
# 
# 
# # > head(patient_data)
# # patient disease            response clin_trial simple_response general_response determinant_response
# # 1 p959-103    pALL                None        959            None    Non-responder                PR_NR
# # 2 p959-117    pALL                None        959            None    Non-responder                PR_NR
# # 3 p959-100    pALL            Complete        959        Complete        Responder              CR_PRtd
# # 4 p959-112    pALL            Complete        959        Complete        Responder              CR_PRtd
# # 5 p959-107    pALL Complete\nw relapse        959        Complete        Responder              CR_PRtd
# # 6 p959-101    pALL Complete\nw relapse        959        Complete        Responder              CR_PRtd
# 
# # > head(specimen_data)
# # trial   patient        celltype timepoint specimenaccnum
# # 1 CART19 p04409-10 Tcells:CAR+CD8+      d121       GTSP0560
# # 2 CART19 p04409-01 Tcells:CAR+CD8-       d14       GTSP0561
# # 3 CART19 p04409-01 Tcells:CAR+CD8+       d14       GTSP0562
# # 4 CART19 p04409-02 Tcells:CAR+CD8-       d28       GTSP0563
# # 5 CART19 p04409-02 Tcells:CAR+CD8+       d28       GTSP0564
# # 6 CART19 p04409-09 Tcells:CAR+CD8-       d14       GTSP0565
# 
# 
# 
# 
# # 
# # 
# # ## Specimen Summary ------------------------------------------------------------
# # 
# # specimen_data <- patient_data %>%
# #   dplyr::full_join(specimen_data, by = "patient") %>%
# #   dplyr::rename(specimen = specimenaccnum) %>%  
# #   dplyr::inner_join(stats, by = "specimen") %>%
# #   dplyr::left_join(gen_epi_summary, by = "specimen")
# # 
# # specimen_data$clustersRepresented <- sapply(
# #   specimen_data$specimen, function(x){
# #     sites <- cond_uniq_sites[cond_uniq_sites$specimen == x]
# #     hits <- findOverlaps(sites, red_clusters)
# #     length(unique(subjectHits(hits)))
# #   })
# 
# r$clustersRepresented <- sapply(r$SpecimenAccNum, function(x){
#   sites <- cond_uniq_sites[cond_uniq_sites$specimen == x]
#   hits <- findOverlaps(sites, red_clusters)
#   length(unique(subjectHits(hits)))
# })
# 
# r[is.na(r$numUniqSites),]$clustersRepresented <- NA
# 
# 
# # specimen_data$numSitesInClusters <- sapply(
# #   specimen_data$specimen, function(x){
# #     sites <- cond_uniq_sites[cond_uniq_sites$specimen == x]
# #     hits <- findOverlaps(sites, red_clusters)
# #     length(unique(queryHits(hits)))
# #   })
# 
# 
# r$numSitesInClusters <- sapply(r$SpecimenAccNum, function(x){
#        sites <- cond_uniq_sites[cond_uniq_sites$specimen == x]
#        hits <- findOverlaps(sites, red_clusters)
#        length(unique(queryHits(hits)))
#      })
#   
# r[is.na(r$numUniqSites),]$numSitesInClusters <- NA
# 
#  
# # specimen_data$abundInClusters <- sapply(
# #   specimen_data$specimen, function(x){
# #     sites <- cond_uniq_sites[cond_uniq_sites$specimen == x]
# #     hits <- findOverlaps(sites, red_clusters)
# #     sum(sites[queryHits(hits)]$estAbund)
# #   })
# # 
# 
# r$abundInClusters <- sapply(r$SpecimenAccNum, function(x){
#       sites <- cond_uniq_sites[cond_uniq_sites$specimen == x]
#       hits <- findOverlaps(sites, red_clusters)
#       sum(sites[queryHits(hits)]$estAbund)
#     })
# 
# r[is.na(r$numUniqSites),]$abundInClusters <- NA
# 
# 
# 
# 
# 
# openxlsx::write.xlsx(r, file = 'CART.xlsx')
# 
# 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
