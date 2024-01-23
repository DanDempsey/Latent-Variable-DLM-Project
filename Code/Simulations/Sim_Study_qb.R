##### Simulation Study
##### Daniel Dempsey

### Load in libraries
library( dplyr )
library( tictoc )
library( dlnm )
library( MASS )
library( parallel )
library( mice )
source( 'Code/Main_Software/DLM_MCMC_qb.R' )
RNGkind( "L'Ecuyer-CMRG" ) # Sets a RNG stream that carries across parallel computing, so that code is reproducible

wd <- 'Output/Simulations/QB/'
dir.create( wd )
setwd( wd )

### Complete data
set.seed( 1001 )
X <- complete( mice(chicagoNMMAPS) )

### Simulate data
nealmon <- function( delta_t, theta ) {
  
  p <- length( theta )
  lw <- matrix(rep(delta_t, each = p) ^ (1:p), ncol = p, byrow = TRUE) %*% theta
  exp(lw - matrixStats::logSumExp(lw))
  
}

# Static variables
num_static <- 2
static_dat <- matrix( abs(rnorm(nrow(X)*num_static)), ncol = num_static ) %>% apply( 2, scale )
colnames( static_dat ) <- paste0( 'static', 1:num_static ) 

# Dynamic variable
lag_len <- 40
dynamic_dat <- apply( X[, 12:14], 2, scale )
lag_profile1 <- nealmon( 0:lag_len, c(0.32, -0.02) )
lag_profile2_raw <- lag_len - (0:lag_len)
lag_profile2 <- lag_profile2_raw / sum(lag_profile2_raw)
lag_profile3 <- 0.7*nealmon( 0:lag_len, c(0.2, -0.03) ) + 0.3*nealmon( 0:lag_len, c(6, -0.1) )

par( mfrow = c(1, 3) )
plot( lag_profile1, type = 'l', xlab = 'Lag', ylab = '', main = 'Lag Profile 1' ) # Delayed lag effect
grid()
plot( lag_profile2, type = 'l', xlab = 'Lag', ylab = '', main = 'Lag Profile 2' ) # Monotonic decreasing lag effect
grid()
plot( lag_profile3, type = 'l', xlab = 'Lag', ylab = '', main = 'Lag Profile 3' ) # Bimodal lag effect
grid()
par( mfrow = c(1, 1) )

rhum_lag <- tsModel::Lag( dynamic_dat[, 1], 0:lag_len )
pm10_lag <- tsModel::Lag( dynamic_dat[, 2], 0:lag_len )
o3_lag <- tsModel::Lag( dynamic_dat[, 3], 0:lag_len )

rhum_dlf <- rhum_lag %*% lag_profile1
pm10_dlf <- pm10_lag %*% lag_profile2
o3_dlf <- o3_lag %*% lag_profile3

dynamic_dlf <- data.frame( rhum = rhum_dlf, pm10 = pm10_dlf, o3 = o3_dlf )

# Combine
covars <- cbind( temp = static_dat, rhum = dynamic_dat )
covars_expanded <- cbind( 'Intercept' = 1, cbind(static_dat, dynamic_dlf) ) %>% na.omit %>% as.matrix

# True parameter values
N <- nrow( covars_expanded )
true_bet <- c( 0, -0.5, 0.01, # static effects
               0.5, 0.01, -0.5 ) # dynamic effects
quant <- 0.9

run_dlm <- function( x, dat, group_labs = NULL, nsamp = 25000, nburn = 25000, thin = 10 ) {
  DLM_MCMC_qb( y ~ ., data = dat, groups = group_labs, quantile = quant,
               nsamp = nsamp, nburn = nburn, thin = thin )
}

### Loop over simulations
sim_res <- list()
outlen <- 4
inlen <- 30

tic()
for ( i in 1:outlen ) {
  
  print( i )
  
  sim_res[[i]] <- list()
  
  ystar <- rALD( nrow(covars_expanded), mu = covars_expanded %*% true_bet, p = quant )
  y <- ifelse( ystar > 0, 1, 0 )
  
  # Combine
  sim_dat <- cbind( 'y' = c(rep(0, lag_len), y), covars ) %>% as.data.frame
  
  # Set degrees of freedom for b-splines
  bs_df <- 7
  
  ### Frequentist fit
  X <- dplyr::select( sim_dat, -y )
  X_static <- dplyr::select( X, paste0('static', 1:2) )
  X_dynamic_raw <- dplyr::select( X, 'rhum', 'pm10', 'o3' ) 
  X_dynamic <- lapply( X_dynamic_raw, crossbasis, lag = lag_len, arglag = list(fun = 'bs', df = bs_df),
                       argvar = list(fun = 'identity') ) %>% 
    do.call( what = 'cbind' )# %>% apply( 2, scale )
  colnames( X_dynamic ) <- paste0( rep( colnames(X_dynamic_raw), each = bs_df ), c('', as.character(2:bs_df)) )
  
  full_dat <- cbind( y = sim_dat$y, cbind(X_static, X_dynamic) ) %>% na.omit %>% data.frame
  
  len_static <- NCOL( X_static ) + 1
  group_labs_dyn <- rep( 1:ncol(X_dynamic_raw), each = bs_df ) + len_static
  group_labs <- c( 1:len_static, group_labs_dyn )
  
  ### Bayesian fit (with variable selection)
  sim_res[[i]]$MCMC <- mclapply( 1:inlen, run_dlm, dat = full_dat, mc.cores = 6 )
  
}
toc() # 5475.067 sec elapsed

### Save results
save( sim_res, file = 'Simulation_Results.Rdata' )

#### Smaller set
X <- X[1:290, ]
covars <- covars[1:290, ] 
covars_expanded <- covars_expanded[1:250, ]

# True parameter values
N <- nrow( covars_expanded )

### Loop over simulations
small_sim_res <- list()

set.seed( 2002 )
tic()
for ( i in 1:outlen ) {
  
  print( i )
  
  small_sim_res[[i]] <- list()
  
  ystar <- rALD( nrow(covars_expanded), mu = covars_expanded %*% true_bet, p = quant )
  y <- ifelse( ystar > 0, 1, 0 )
  
  # Combine
  sim_dat <- cbind( 'y' = c(rep(0, lag_len), y), covars ) %>% as.data.frame
  
  # Set degrees of freedom for b-splines
  bs_df <- 7
  
  ### Frequentist fit
  X <- dplyr::select( sim_dat, -y )
  X_static <- dplyr::select( X, paste0('static', 1:2) )
  X_dynamic_raw <- dplyr::select( X, 'rhum', 'pm10', 'o3' ) 
  X_dynamic <- lapply( X_dynamic_raw, crossbasis, lag = lag_len, arglag = list(fun = 'bs', df = bs_df),
                       argvar = list(fun = 'identity') ) %>% 
    do.call( what = 'cbind' )# %>% apply( 2, scale )
  colnames( X_dynamic ) <- paste0( rep( colnames(X_dynamic_raw), each = bs_df ), c('', as.character(2:bs_df)) )
  
  full_dat <- cbind( y = sim_dat$y, cbind(X_static, X_dynamic) ) %>% na.omit %>% data.frame
  
  small_sim_res[[i]]$MLE <- glm.nb( y ~ ., data = full_dat )
  
  len_static <- NCOL( X_static ) + 1
  group_labs_dyn <- rep( 1:ncol(X_dynamic_raw), each = bs_df ) + len_static
  group_labs <- c( 1:len_static, group_labs_dyn )
  
  ### Bayesian fit (with variable selection)
  small_sim_res[[i]]$MCMC <- mclapply( 1:inlen, run_dlm, dat = full_dat, mc.cores = 6 )
  
}
toc() # 742.69 sec elapsed

### Save results
save( small_sim_res, file = 'Small_Simulation_Results.Rdata' )

