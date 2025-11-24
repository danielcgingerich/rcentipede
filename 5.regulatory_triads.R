cond <- grepl('/hpc', getwd())
if (cond){
  source('/hpc/group/adrc/dcg27/african_american_multiome/scripts/config.R')
} else {
  source('C:/Users/danie/Desktop/african_american_multiome/scripts/config.R')
} ; rm(cond)
set.seed(1)

slurm_array_task_id <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric()

# active binding site locations from CENTIPEDE
setwd('C:/Users/danie/Desktop/african_american_multiome/data/15.centipede3/')
binding_sites <- list.files(pattern = '4.centipede___.*.rds')
cluster_id <- gsub('.*___|_[[:upper:]].*.rds', '', binding_sites)
cluster_id <- unique(cluster_id)[slurm_array_task_id]
binding_sites <- paste0('4.centipede3___', cluster_id, '.*.rds')
binding_sites <- list.files(pattern = binding_sites)
binding_sites <- lapply(binding_sites, readRDS)
binding_sites <- lapply(binding_sites, function(x) x$Z)

binding_sites <- do.call(what = rbind, args = binding_sites)
ix <- binding_sites[, 1] >= 0.95
binding_sites <- rownames(binding_sites)[ix]
binding_sites <- data.frame(tfbs = gsub('_.*', '', binding_sites), 
                            motif = gsub('chr.*-[[:digit:]]+-[[:digit:]]+_|_.$', '', binding_sites))

# library(TFBSTools)
# library(JASPAR2020)
# pfm <- getMatrixSet(
#   x = JASPAR2020,
#   opts = list(collection = "CORE", species = 9606, all_versions = FALSE)
# )
# motif_names <- sapply(pfm, function(x) x@name)
# setwd(centipede_objects)
# saveRDS(motif_names, '5.regulatory_triads___motif_names.rds')
setwd(centipede_objects)
motif_names <- readRDS('5.regulatory_triads___motif_names.rds')
binding_sites$tf <- motif_names[binding_sites$motif]

# linger output
trans_reg <- paste0('1.format_linger_output___', cluster_id, '_trans_net_norm.rds')
trans_reg <- readRDS(trans_reg)

trans_cutoff <- as.numeric(trans_reg)
trans_cutoff <- quantile(trans_cutoff, probs = 0.8)
ix <- which(trans_reg >= trans_cutoff, arr.ind = TRUE)
tg <- rownames(trans_reg)[ix[, 'row']]
tf <- colnames(trans_reg)[ix[, 'col']]
score <- trans_reg[ix]
trans_reg <- data.frame(tg = tg, tf = tf, normalized_trans_score = score)

cis_reg <- paste0('1.format_linger_output___', cluster_id, '_cis_net_norm.rds')
cis_reg <- readRDS(cis_reg)
cis_cutoff <- quantile(cis_reg$norm2, probs = 0.8)
ix <- cis_reg$norm2 >= cis_cutoff
cis_reg <- cis_reg[ix, ]
colnames(cis_reg) <- c('peak', 'gene', 'cis_score', 'normalized_cis_score1', 'normalized_cis_score2')

# find overlaps
tfbs_gr <- StringToGRanges(unique(binding_sites$tfbs))
peaks_gr <- StringToGRanges(unique(cis_reg$peak))
ix <- findOverlaps(query = tfbs_gr, subject = peaks_gr)
df <- data.frame(tfbs = GRangesToString(tfbs_gr)[queryHits(ix)],
                 peak = GRangesToString(peaks_gr)[subjectHits(ix)])


triads <- inner_join(binding_sites, df, by = c('tfbs' = 'tfbs'))
triads <- inner_join(triads, cis_reg, by = c('peak' = 'peak'), 
                     relationship = 'many-to-many')

trans_reg$id <- paste0(trans_reg$tf, '_', trans_reg$tg)
triads$id <- paste0(triads$tf, '_', triads$gene)
triads <- inner_join(triads, trans_reg[, c('id', 'normalized_trans_score')], by = 'id',
                     relationship = 'many-to-many')

triads <- triads[, c('tf', 'peak', 'gene', 'tfbs', 'motif', 'normalized_cis_score2', 'normalized_trans_score')]
colnames(triads) <- c('tf', 'peak', 'gene', 'tfbs', 'motif', 'normalized_peak_gene', 'normalized_tf_gene')

setwd(centipede_objects)
saveRDS(triads, paste0('5.regulatory_triads___', cluster_id, '.rds'))

sesh <- capture.output(sessionInfo())
print(sesh)