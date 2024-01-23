###### Chicago Air Quality example analysis
###### Daniel Dempsey

### Load in libraries and set working directory
library( dlnm ) # Chicago dataset
library( magrittr ) # For the piping operator %>%
library( caret ) # One-hot encoding
library( mice ) # Missing value imputation
library( parallel ) # Parallel computation
library( tictoc ) # Time code
library( MASS ) # For Negative Binomial MLE
library( corrplot ) # Correlation matrix plots
source( 'Code/Main_Software/DLM_MCMC_qb.R' )
RNGkind("L'Ecuyer-CMRG") # Sets a RNG stream that carries across parallel computing, making code reproducible

wd <- 'Output/Application/QB/'
dir.create( wd )
setwd( wd )

### Impute missing values
set.seed( 1001 )
X_full <- complete( mice(chicagoNMMAPS) )
N <- nrow( X_full )

### Create desired variables
X_full$dow <- as.factor( X_full$dow )
seasons <- c( rep('Winter', 2), rep('Spring', 3), rep('Summer', 3), rep('Autumn', 3), 'Winter' )
X_full$season <- seasons[X_full$month] %>% as.factor

### Use a smaller set
train_prop <- 0.8
train_inds <- ceiling( N*train_prop ) %>% seq_len
test_inds <- setdiff( 1:N, train_inds )

X <- X_full[train_inds, ]
X_test <- X_full[test_inds, ]

### One-hot encode days of the week and set static and dynamic variables
n_df <- 7
lag_len <- 40
X_static <- dplyr::select( X, c('dow', 'season') )
X_dynamic_raw <- dplyr::select( X, 'temp', 'dptp', 'rhum', 'pm10', 'o3' ) 
X_dynamic <- lapply( X_dynamic_raw, crossbasis, lag = lag_len, arglag = list(fun = 'bs', df = n_df) ) %>% 
  do.call( what = 'cbind' ) %>% apply( 2, scale )
colnames( X_dynamic ) <- paste0( rep( colnames(X_dynamic_raw), each = n_df ), c('', as.character(2:n_df)) )

len_static <- sum( sapply( X_static, function(x) { length(levels(x)) - 1 } ) ) + 1
group_labs_dyn <- rep( 1:ncol(X_dynamic_raw), each = n_df ) + len_static
group_labs <- c( 1:len_static, group_labs_dyn )

# Bin the response variable
quant <- 0.9
full_dat <- cbind( death = ifelse( X$death >= quantile(X$death, quant), 1, 0 ),
                   cbind(X_static), X_dynamic ) %>% na.omit

### Set 4 different starting positions for each parameter
nvar <- length( group_labs )

gam1 <- c( TRUE, rep(FALSE, nvar - 1) ) # Sparse start
gam2 <- rep( TRUE, nvar ) # Saturated start
gam3 <- c( TRUE, as.logical( rbinom(nvar - 1, 1, 0.5) ) ) # Random start 1
gam4 <- c( TRUE, as.logical( rbinom(nvar - 1, 1, 0.5) ) ) # Random start 2

init_all <- list()
for ( i in 1:4 ) {
  init_all[[i]] <- list( beta = rnorm(nvar), gamma = get(paste0('gam', i)), xi = 1, pi = 0.5 )
}

### Run the model for each starting position
run_dlm <- function( x, nsamp = 50000, nburn = 50000, thin = 10 ) {
  DLM_MCMC_qb( death ~ ., data = full_dat, quantile = quant,
               init = x, nsamp = nsamp, nburn = nburn, thin = thin )
}

set.seed( 2023 )
tic()
mcmc_fit <- mclapply( init_all, run_dlm, mc.cores = 4 )
toc() # 361.244 sec elapsed

### Write results to file
all_res <- list( MCMC = mcmc_fit )
save( all_res, file = 'all_res.Rdata' )
#load( file = 'all_res.Rdata' )

### Visualisations
res <- all_res$MCMC

# Gamma
gamma_mean <- function( i ) {
  apply( res[[i]]$gamma, 2, mean )
}

gam_dat <- do.call( 'rbind', lapply( 1:4, gamma_mean ) )[, -1] * 100

pdf( 'Chicago_Gamma_QB.pdf', width = 12 )
plot( 0, xlim = c(1, ncol(gam_dat)), ylim = c(0, 100), xaxt = 'n', yaxt = 'n', type = 'n', xlab = '', ylab = '',
      main = 'MCMC Proportion of Variables Selected' )
abline( h = seq(0, 100, 25), lty = 3, col = adjustcolor('grey', 0.75) )
boxplot( gam_dat, xaxt = 'n', yaxt = 'n', add = TRUE, pch = 20, col = 'blue' )
axis( 1, at = c(2.5, 8, 13, 20, 27, 34, 41), tick = FALSE,
      labels = c('DOW', 'Season', 'Temp', 'Dew Temp', 'Humidity', 'PM10', 'O3') )
axis( 2, at = seq(0, 100, 25), labels = TRUE, las = 2 )
abline( v = c(6.5, 9.5, 16.5, 23.5, 30.5, 37.5), lty = 2, col = 'blue' )
dev.off()

# DLFs 
mycols <- c( '#56B4E9', '#E69F00', '#009E73', '#F0E442' ) %>% rev

beta_extract <- function( i, dat ) {
  dat[[i]]$beta
}

group_extract <- function( i, dat ) {
  lapply( dat, '[', i = , j = i )
} 

plot_fun <- function( x, nm, ls, use_col = mycols ) {
  
  fitted_weights_mean <- list()
  for ( i in 1:length(x) ) {
    fitted_weights_raw <- ls %*% t( x[[i]] )
    fitted_weights_raw_mean <- apply( fitted_weights_raw, 1, mean )
    mean_effect <- sum( fitted_weights_raw_mean )
    fitted_weights_mean[[i]] <- fitted_weights_raw_mean / mean_effect
    cat( paste0(nm, ': ', round(mean_effect, 2), '\n') )
  }
  
  ylim <- range( fitted_weights_mean, na.rm = TRUE )
  plot( 0, ylim = ylim, xlim = c(0, lag_len), type = 'n', ylab = '', xlab = 'Lag', main = nm )
  grid()
  
  for( i in 1:length(x) ) {
    lines( fitted_weights_mean[[i]], lwd = 3, col = use_col[i] )
  }
  
}

DLF_plot <- function( dat ) {
  
  temp_inds <- 11:17
  dptp_inds <- 18:24
  rhum_inds <- 25:31
  pm10_inds <- 32:38
  o3_inds <- 39:45
  
  lag_splines <- splines::bs( 0:lag_len, df = 7, intercept = TRUE )
  all_bet <- lapply( 1:4, beta_extract, dat )
  temp_bet <- group_extract( temp_inds, all_bet )
  dptp_bet <- group_extract( dptp_inds, all_bet )
  rhum_bet <- group_extract( rhum_inds, all_bet )
  pm10_bet <- group_extract( pm10_inds, all_bet )
  o3_bet <- group_extract( o3_inds, all_bet )
  
  par( mfrow = c(2, 3) )
  plot_fun( temp_bet, 'Temperature', lag_splines )
  plot_fun( dptp_bet, 'Dew Point Temperature', lag_splines )
  plot_fun( rhum_bet, 'Relative Humidity', lag_splines )
  plot_fun( pm10_bet, 'PM10', lag_splines )
  plot_fun( o3_bet, 'O3', lag_splines )
  par( mfrow = c(1, 1) )
  
}

pdf( 'Chicago_DLF_QB.pdf', width = 12 )
DLF_plot( res )
dev.off()

