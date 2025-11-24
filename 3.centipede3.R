cond <- grepl('/hpc', getwd())
if (cond){
  source('/path/to/config.R')
} else {
  source('C:/path/to/config.R')
} ; rm(cond)
set.seed(1)

library(nnet)
library(MASS)
library(mixtools)

slurm_array_task_id = Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric()

# prepare data
setwd(centipede_objects)
setwd('/hpc/group/adrc/dcg27/african_american_multiome/data/15.centipede3')
# setwd('C:/Users/danie/Desktop/african_american_multiome/data/15.centipede3/')
data <- list.files(pattern = '3.get_cutsites3')
cluster_id <- gsub('.*___|_chr.*.rds', '', data)
cluster_id <- unique(cluster_id)
cluster_id <- cluster_id[slurm_array_task_id]
data <- list.files(pattern = paste0('3.get_cutsites3___', cluster_id))
data <- lapply(data, readRDS)

# format data objects
P_all <- lapply(data, function(x) x$counts) # %>% as.matrix()
P_all <- do.call(what = rbind, args = P_all)
scores <- lapply(data, function(x) x$scores)
scores <- do.call(what = rbind, args = scores)
tfs <- scores$tf
thresh <- quantile(scores$trp, probs = 0.8)
scores$trp <- ifelse(scores$trp >= thresh, 1, 0)
X_all <- cbind(1, scores$score, scores$trp)
rownames(X_all) <- names(tfs) <- scores$id
colnames(X_all) <- c('intercept', 'score', 'trp')
R_all <- rowSums(P_all)

# remove motifs with 0 cutsites
ix <- R_all > 0  
ix <- names(R_all)[ix]
P_all <- P_all[ix, ]
R_all <- R_all[ix]
X_all <- X_all[ix, ]
tfs <- tfs[ix]

fix_linger = FALSE
for (tf_i in unique(tfs)){
  try(
    expr = {
      # subset the tf
      message(tf_i)
      ix <- tfs == tf_i
      ix <- names(tfs)[ix]
      P <- P_all[ix, ]
      R <- R_all[ix]
      X <- X_all[ix, ]
      
      # reverse the minus strand
      ix <- grepl('-$', rownames(P))
      Prev <- P[ix, , drop = FALSE]
      Pfwd <- P[!ix, , drop = FALSE]
      Prev <- Prev[, 200:1]
      P <- rbind(Pfwd, Prev)
      
      # impose symmetry on P
      P <- P[, 1:100] + P[, 200:101]
      
      # make sure R, P and X are ordered
      ix <- rownames(P)
      R <- R[ix]
      X <- X[ix, ]
      
      # scale X
      # X[, 2:ncol(X)] <- scale(X[, 2:ncol(X)])
      
      # initialize with k means clustering
      clust <- kmeans(x = R, centers = 2)
      ix <- which.max(clust$centers)
      ix <- clust$cluster == ix
      R1 <- R[ix]
      R0 <- R[!ix]
      P1 <- P[ix, ] %>% as.matrix()
      P0 <- P[!ix, ] %>% as.matrix()
      
      lambda1 <- colSums(P1) / sum(P1)
      
      B = c(0, 0, 0)
      
      h1 <- glm.nb(R1 ~ 1)
      u1 <- exp(coef(h1))
      phi1 <- h1$theta
      
      h0 <- glm.nb(R0 ~ 1)
      u0 <- exp(coef(h0))
      phi0 <- h0$theta
      
      params = c(B, lambda1, u0, u1, phi0, phi1)
      names(params) <- c(paste0('B', 1:length(B)), 
                         paste0('lambda', 1:ncol(P)), 
                         'u0', 'u1', 'phi0', 'phi1')
      
      i = 1
      epsilon = 1
      while (epsilon > 1e-3){
        
        # e step:
        # initial params
        B = params[paste0('B', 1:length(B))]
        
        lambda0 = rep(1/ncol(P), ncol(P))
        ix <- grep('lambda', names(params))
        lambda1 = params[ix]
        
        u0 = params['u0']
        phi0 = params['phi0']
        
        u1 = params['u1']
        phi1 = params['phi1']
        
        # calculate posterior
        message('updating posterior')
        b = c(B[1], exp(B[2:length(B)]))
        pi1 = plogis(X %*% b, log.p = TRUE)
        pi0 = plogis(-X %*% b, log.p = TRUE)
        
        nb1 = dnbinom(x = R, mu = u1, size = phi1, log = TRUE)
        nb0 = dnbinom(x = R, mu = u0, size = phi0, log = TRUE)
        
        mn1 = P %*% log(lambda1)
        mn0 = P %*% log(lambda0)
        
        scale_factors = -1 * pmin( pi1 + mn1 + nb1, pi0 + mn0 + nb0)
        Z = exp(pi1 + mn1 + nb1 + scale_factors) / ( exp(pi1 + mn1 + nb1 + scale_factors) + exp(pi0 + mn0 + nb0 + scale_factors))
        
        ix <- is.na(Z)
        if (any(ix)){
          message('NA values found')
          z1 <- pi1 + mn1 + nb1 + scale_factors
          z0 <- pi0 + mn0 + nb0 + scale_factors
          z <- z1[ix] > z0[ix]
          z <- ifelse(z, 1, 0)
          Z[ix] <- z
        }
        
        # pre update ll
        loglikelihood0 = sum(
          Z * pi1 + Z * mn1 + Z * nb1 + (1-Z) * pi0 + (1-Z) * mn0 + (1-Z) * nb0 
        )
        
        # m step:
        message('updating params')
        # optimize multinomial
        Z_matrix = lapply(1:ncol(P), function(x){Z})
        Z_matrix <- do.call(cbind, args = Z_matrix)
        M <- Z_matrix * P
        lambda1 = colSums(M) / sum(M)
        
        # optimize negative binom
        h1 <- glm.nb(R ~ 1, weights = as.numeric(Z))
        u1 <- exp(coef(h1))
        phi1 <- h1$theta
        
        h0 <- glm.nb(R ~ 1, weights = 1 - as.numeric(Z))
        u0 <- exp(coef(h0))
        phi0 <- h0$theta
        
        # optimize pi
        if (fix_linger){
          B_ = B[1:2]
          B_ = optim(B_, fn = function(b){
            beta = c(b[1], exp(b[2:length(b)]), 1)
            -1 * sum( Z * plogis(X %*% beta, log.p = TRUE) + (1 - Z) * plogis(- X %*% beta, log.p = TRUE) )
          })$par
          B[1:2] <- B_
        } else {
          B = optim(B, fn = function(b){
            beta = c(b[1], exp(b[2:length(b)]))
            -1 * sum( Z * plogis(X %*% beta, log.p = TRUE) + (1 - Z) * plogis(- X %*% beta, log.p = TRUE) )
          })$par
        }
        
        
        # update parameters
        new_params = c(B, lambda1, u0, u1, phi0, phi1)
        names(new_params) <- c(paste0('B', 1:length(B)), 
                               paste0('lambda', 1:ncol(P)), 
                               'u0', 'u1', 'phi0', 'phi1')
        params <- new_params
        
        # post update ll
        B = params[paste0('B', 1:length(B))]
        
        lambda0 = rep(1/ncol(P), ncol(P))
        ix <- grep('lambda', names(params))
        lambda1 = params[ix]
        
        u0 = params['u0']
        phi0 = params['phi0']
        
        u1 = params['u1']
        phi1 = params['phi1']
        
        b = c(B[1], exp(B[2:length(B)]))
        pi1 = plogis(X %*% b, log.p = TRUE)
        pi0 = plogis(-X %*% b, log.p = TRUE)
        
        mn1 = P %*% log(lambda1)
        mn0 = P %*% log(lambda0)
        
        nb1 = dnbinom(x = R, mu = u1, size = phi1, log = TRUE)
        nb0 = dnbinom(x = R, mu = u0, size = phi0, log = TRUE)
        
        loglikelihood1 = sum(
          Z * pi1 + Z * mn1 + Z * nb1 + (1-Z) * pi0 + (1-Z) * mn0 + (1-Z) * nb0 
        )
        epsilon = abs(loglikelihood0 - loglikelihood1)
        
        message('EM ', i, ': log likelihood = ', round(loglikelihood1))
        message('epsilon = ', round(epsilon, digits = 3))
        i = i + 1
        
        # save progress at each iteration
        em <- list(params = params, Z = Z)
        setwd(centipede_objects)
        filename <- paste0('4.centipede3___', cluster_id, '_', tf_i, '.rds')
        filename <- gsub(':', '-', filename)
        # setwd('/hpc/group/adrc/dcg27/african_american_multiome/data/15.centipede3')
        # saveRDS(em, file = filename)
      }
      em <- list(params = params, Z = Z)
      
      # param inspection
      # if (fix_linger){
      #   Z_linger_is_offset = em$Z
      # } else {
      #   Z_linger_is_parameter = em$Z
      # }
      # ix <- as.numeric(Z_linger_is_parameter) >= 0.95
      # mean(X[ix, c('trp')])
      # mean(X[!ix, c('trp')])
      # 
      # ix <- as.numeric(Z_linger_is_offset) >= 0.95
      # mean(X[ix, c('trp')])
      # mean(X[!ix, c('trp')])
      # 
      # cor(Z_linger_is_offset[, 1], Z_linger_is_parameter[, 1])
      # 
      # ix_parameter = as.numeric(Z_linger_is_parameter) >= 0.95
      # ix_offset <- as.numeric(Z_linger_is_offset) >= 0.95
      # ix = which(ix_parameter != ix_offset)
      # ix = rownames(Z_linger_is_offset)[ix]
      # 
      # P_diff = P[ix, ]
      # plot(colSums(P_diff))
      # 
      # exp(em$params['B3'])

      # visual inspection
      # ix <- tfs == tf_i
      # ix <- names(tfs)[ix]
      # P <- P_all[ix, ]
      # ix <- grepl('-$', rownames(P))
      # Prev <- P[ix, ]
      # Pfwd <- P[!ix, ]
      # Prev <- Prev[, 200:1]
      # P <- rbind(Pfwd, Prev)
      # 
      # bound <- as.numeric(Z_linger_is_parameter) >= 0.95
      # bound <- rownames(em$Z)[bound]
      # bound <- colSums(P[bound, , drop = FALSE])
      # 
      # unbound <- as.numeric(em$Z) <= 0.95
      # unbound <- rownames(em$Z)[unbound]
      # unbound <- colSums(P[unbound, , drop = FALSE])
      # 
      # plot(x = 1:200, y = bound)
      # plot(x = 1:200, y = unbound)

      setwd(centipede_objects)
      setwd('/hpc/group/adrc/dcg27/african_american_multiome/data/15.centipede3')
      saveRDS(em, file = filename)
    })
}

sesh <- capture.output(sessionInfo())

print(sesh)
