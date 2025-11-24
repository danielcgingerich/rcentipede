# srun --mem=100GB --pty bash -i
# export PATH="/hpc/group/adrc/dcg27/software/anaconda3/bin/conda:$PATH"
# conda init
# conda activate /hpc/group/adrc/dcg27/software/anaconda3/envs/r_conda_installation___2024_04_01
# R
cond <- grepl('/hpc', getwd())
if (cond){
  source('/hpc/group/adrc/dcg27/african_american_multiome/scripts/config.R')
} else {
  source('C:/Users/danie/Desktop/african_american_multiome/scripts/config.R')
} ; rm(cond)
set.seed(1)

library(Rsamtools)
library(rtracklayer)

slurm_array_task_id <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric()

setwd(centipede_objects)
setwd('/hpc/group/adrc/dcg27/african_american_multiome/data/15.centipede3')
pos <- list.files(pattern = '2.motifmatchr3.*.rds')
cluster_id <- gsub('2.motifmatchr3___|.rds', '', pos)
chr_id <- c(paste0('chr', 1:22), 'chrX', 'chrY')
df <- expand.grid(cluster_id, chr_id)
cluster_id <- df$Var1[slurm_array_task_id]
chr_id <- df$Var2[slurm_array_task_id]
pos <- paste0('2.motifmatchr3___', cluster_id, '.rds')
pos <- readRDS(pos)
pos <- renameSeqlevels(pos, value = paste0('chr', seqlevels(pos)))
ix <- as.character(seqnames(pos)) == chr_id
pos <- pos[ix]

setwd(centipede_objects)
cis_net <- paste0('1.format_linger_output___', cluster_id, '_cis_net_norm.rds')
trans_net <- paste0('1.format_linger_output___', cluster_id, '_trans_net_norm.rds')
cis_net <- readRDS(cis_net)
trans_net <- readRDS(trans_net)

# for each peak, take maximum cis regulatory potential
cis_net <- cis_net[order(- cis_net$norm2), ]
ix <- duplicated(cis_net$peak)
cis_net <- cis_net[!ix, ]

peaks <- cis_net$peak
peaks <- StringToGRanges(peaks)
ix <- findOverlaps(query = pos, subject = peaks)
pos$peak <- GRangesToString(peaks)[subjectHits(ix)]
pos$tg <- cis_net$gene[subjectHits(ix)]
pos$crp_max <- cis_net$norm2[subjectHits(ix)]
peaks <- pos$peak
pos$tf <- gsub('\\(.*\\)', '', pos$tf)
ix <- pos$tf %in% colnames(trans_net)
pos <- pos[ix, ]
ix <- cbind(pos$tg, pos$tf)
pos$trp <- trans_net[ix]

# strand(pos) <- '+'
# pos <- unique(pos)
motif_windows <- pos
motif_windows$id <- GRangesToString(motif_windows)
motif_windows$id <- paste0(motif_windows$id, '_',  motif_windows$motif, '_', strand(motif_windows))
motif_windows <- Extend(motif_windows, 100, 100)
motif_windows$window_start <- start(motif_windows)
motif_windows$window_end <- end(motif_windows)
# motif_windows <- renameSeqlevels(motif_windows, paste0('chr', seqlevels(motif_windows)))
motif_windows <- motif_windows[!duplicated(motif_windows$id)]

motif_index <- data.table(id = motif_windows$id, row_index = as.numeric(factor(motif_windows$id)))
motif_index <- unique(motif_index)
motif_index <- motif_index[order(row_index), ]

regions <- StringToGRanges(unique(peaks))

setwd(multimodal_clustering_objects)
cells <- readRDS('6.final_clusters___cluster_labels.rds')
ix <- cells$cluster == cluster_id
cells <- cells[ix, , drop = FALSE]
cells$sample_id <- gsub('_.*', '', rownames(cells))
sample_id <- unique(cells$sample_id)

for (i in 1:length(sample_id)){
  message(i, '/', length(sample_id))
  path <- paste0(cellranger_arc_objects, '/8.cellranger_arc_count/', sample_id[i], '/outs')
  setwd(path)

  fragment_file <- 'atac_fragments.tsv.gz'
  fragment_index <- 'atac_fragments.tsv.gz.tbi'
  
  tabix_file <- TabixFile(fragment_file, index = fragment_index)
  fragments_region <- scanTabix(tabix_file, param = regions)

  fragments_region <- unlist(fragments_region)
  fragments_region <- strsplit(fragments_region, split = '\t')
  fragments_region <- unlist(fragments_region)
  fragments_region <- matrix(fragments_region, ncol = 5, byrow = TRUE)
  fragments_region <- as.data.table(fragments_region)
  
  ix <- cells$sample_id == sample_id[i]
  cells_use <- rownames(cells)[ix]
  cells_use <- gsub('.*_', '', cells_use)
  colnames(fragments_region) <- c('chr', 'start', 'end', 'barcode', 'counts')
  ix <- fragments_region$barcode %in% cells_use
  fragments_region <- fragments_region[ix, ]
  fragments_region <- unique(fragments_region)
  
  if (nrow(fragments_region) > 0){
    cutsites_left <- fragments_region[, .(chr, start)]
    cutsites_right <- fragments_region[, .(chr, end)]
    colnames(cutsites_right) <- colnames(cutsites_left) <- c('chr', 'cut')
    cutsites <- rbind(cutsites_left, cutsites_right)
    cutsites$pos <- paste0(cutsites$chr, '-', cutsites$cut, '-', cutsites$cut)
    cutsites <- cutsites[, .N, by = pos]
    tmp <- cutsites
    cutsites <- StringToGRanges(cutsites$pos)
    mcols(cutsites) <- tmp
    
    ovrlp <- findOverlaps(query = cutsites, subject = motif_windows)
    cutsites <- cutsites[queryHits(ovrlp)]
    cutsites$window_start <- motif_windows$window_start[subjectHits(ovrlp)]
    cutsites$window_end <- motif_windows$window_end[subjectHits(ovrlp)]
    cutsites$id <- motif_windows$id[subjectHits(ovrlp)]
    cutsites$p1 <- start(cutsites) - cutsites$window_start
    cutsites$p2 <- cutsites$window_end - start(cutsites)
    
    ix <- cutsites$p1 < cutsites$p2
    cutsites$relative_position <- 0
    cutsites$relative_position[ix] <- start(cutsites)[ix] - cutsites$window_start[ix] + 1
    cutsites$relative_position[!ix] <- 200 - (cutsites$window_end[!ix] - start(cutsites)[!ix])
    
    cutsites <- mcols(cutsites)
    cutsites <- as.data.table(cutsites)
    
    cutsites <- left_join(x = cutsites, y = motif_index, by = 'id')
    n <- length(unique(motif_index$id))
    V <- matrix(0, nrow = n, ncol = 200)
    V <- as(V, 'dgCMatrix')
    ix <- cbind(cutsites$row_index, cutsites$relative_position)
    V[ix] <- cutsites$N
    
    rownames(V) <- motif_index$id
    
    if (i == 1){
      master <- V
    } else {
      master <- master[rownames(V), ] + V
    }
  }
  gc()
}
master <- list(counts = master, scores = mcols(motif_windows))
setwd(centipede_objects)
setwd('/hpc/group/adrc/dcg27/african_american_multiome/data/15.centipede3')
saveRDS(master, paste0('3.get_cutsites3___', cluster_id, '_', chr_id, '.rds'))

sesh <- capture.output(sessionInfo())
print(sesh)