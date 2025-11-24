cond <- grepl('/hpc', getwd())
if (cond){
  source('/path/to/config.R')
} else {
  source('path/to/config.R')
} ; rm(cond)
set.seed(1)

library(motifmatchr)
library(TFBSTools)
library(JASPAR2020)
library(BSgenome.Hsapiens.NCBI.GRCh38)

slurm_array_task_id <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric()

setwd(ref_genome_path)
hg38 <- rtracklayer::import('gencode.v32.primary_assembly.annotation.gtf')

setwd(centipede_objects)
cis_net <- list.files(pattern = '1.format_linger_output___.*_cis_net_norm.rds')[slurm_array_task_id]
cluster_id <- gsub('.*___|_cis_net_norm.rds', '', cis_net)
trans_net <- paste0('1.format_linger_output___', cluster_id, '_trans_net_norm.rds')

cis_net <- readRDS(cis_net)
trans_net <- readRDS(trans_net)

tg <- rownames(trans_net)
tf <- colnames(trans_net)
df <- expand.grid(tg, tf)
colnames(df) <- c('tg', 'tf')

# load motif pwms
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", species = 9606, all_versions = FALSE)
)
pwm <- toPWM(pfm)



# find motifs
tfs <- unique(tf)
tfs <- paste(tfs, collapse = '|')
pwm_names <- sapply(pwm, function(x) x@name)
ix <- grep(pattern = tfs, pwm_names)
pwm <- pwm[ix]


genes <- unique(df$tg)
genes <- hg38[hg38$gene_name %in% genes & hg38$type == 'gene', ]
genes <- Extend(genes, 1e6, 1e6)
genes <- reduce(genes)

peaks <- unique(cis_net$peak)
peaks <- StringToGRanges(peaks)
ix <- findOverlaps(query = peaks, subject = genes)
peaks <- peaks[queryHits(ix)]
peaks <- unique(peaks)


v <- gsub('chr', '', seqlevels(peaks))
peaks <- renameSeqlevels(x = peaks, value = v)
pos <- matchMotifs(pwms = pwm, subject = peaks, 
                     genome = BSgenome.Hsapiens.NCBI.GRCh38,
                     out = 'positions')

motifs <- lapply(pos, length)
motifs <- lapply(1:length(motifs), 
                 function(x){
                   rep(names(motifs)[x], motifs[x])
                 })
motifs <- unlist(motifs)
pos <- unlist(pos)
pos$motif <- motifs
tf_key <- sapply(1:length(pwm), 
                 function(x){
                   motif <- names(pwm)[x]
                   tf <- pwm[[x]]@name
                   names(tf) <- motif
                   return(tf)
                 })
pos$tf <- tf_key[pos$motif]

setwd(centipede_objects)
setwd('/hpc/group/adrc/dcg27/african_american_multiome/data/15.centipede3')
saveRDS(pos, paste0('2.motifmatchr3___', cluster_id, '.rds'))

sesh <- capture.output(sessionInfo())
print(sesh)

