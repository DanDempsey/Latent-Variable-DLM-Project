###### Visualising Simulation Study Results
###### Daniel Dempsey

### Preamble
library( dplyr )
setwd( 'Output/Simulations/' )

### Output function
output_fun <- function( wd ) {
  
  ### Extract results
  setwd( wd )
  load( file = 'Simulation_Results.Rdata' )
  load( file = 'Small_Simulation_Results.Rdata' )
  
  ### Gamma
  gamma_extract <- function( i, dat ) {
    x <- lapply( dat[[i]]$MCMC, '[[', 'gamma' )
    y <- lapply( x, apply, 2, mean )
    do.call( 'rbind', y )
  }
  
  gamma_plot <- function( i, dat ) {
    plot( 0, xlim = c(1, 24), ylim = c(0, 100), xaxt = 'n', yaxt = 'n', type = 'n', xlab = '', ylab = '',
          main = 'MCMC Proportion of Variables Selected' )
    abline( h = seq(0, 100, 25), lty = 3, col = adjustcolor('grey', 0.75) )
    boxplot( gamma_extract(i, dat) * 100, xaxt = 'n', yaxt = 'n', add = TRUE, pch = 20, col = 'blue' )
    axis( 1, at = c(1.5, 7, 14, 21.5), labels = c('Static', 'Relative Humidity', 'PM10', 'O3'),
          tick = FALSE )
    axis( 2, at = seq(0, 100, 25), labels = TRUE, las = 2 )
    abline( v = c(3.5, 10.5, 17.5), lty = 2, col = 'blue' )
  }
  
  pdf( paste0('gamma_sims_', wd, '.pdf'), width = 12 )
  gamma_plot( 1, dat = sim_res )
  dev.off()
  
  pdf( paste0('all_gamma_sims_', wd, '.pdf'), width = 12, height = 12 )
  par( mfrow = c(2, 2) )
  gamma_plot( 1, dat = sim_res )
  gamma_plot( 2, dat = sim_res )
  gamma_plot( 3, dat = sim_res )
  gamma_plot( 4, dat = sim_res )
  par( mfrow = c(1, 1) )
  dev.off()
  
  pdf( paste0('gamma_sims_small_', wd, '.pdf'), width = 12 )
  gamma_plot( 1, dat = small_sim_res )
  dev.off()
  
  pdf( paste0('all_gamma_sims_small_', wd, '.pdf'), width = 12, height = 12 )
  par( mfrow = c(2, 2) )
  gamma_plot( 1, dat = small_sim_res )
  gamma_plot( 2, dat = small_sim_res )
  gamma_plot( 3, dat = small_sim_res )
  gamma_plot( 4, dat = small_sim_res )
  par( mfrow = c(1, 1) )
  dev.off()
  
  ### Beta
  if ( wd == 'NB' ) {
    true_beta <- c( -1, -0.5, 0.01, 0.5, 0.01, -0.5 )
  }
  else {
    true_beta <- c( 0, -0.5, 0.01, 0.5, 0.01, -0.5 )
  }
  
  beta_extract_within <- function( x ) {
    nvar <- ncol( x$beta )
    res <- numeric( nvar )
    for ( j in 1:nvar ) {
      gam <- x$gamma[, j]
      if ( sum(gam) == 0 ) {
        res[j] <- 0
      }
      else {
        res[j] <- mean( x$beta[gam, j] )
      }
    }
    res
  }
  
  beta_extract <- function( i, dat ) {
    
    all_beta <- lapply( dat[[i]]$MCMC, beta_extract_within ) %>% do.call( what = 'rbind' )
    static_diff <- t( t(all_beta[, 1:3]) - true_beta[1:3] )
    
    rhum_inds <- 4:10
    pm10_inds <- 11:17
    o3_inds <- 18:24
    
    lag_splines <- splines::bs( 0:40, df = 7, intercept = TRUE )
    
    fitted_weights_raw <- lag_splines %*% t( all_beta[, rhum_inds] )
    rhum_diff <- colSums( fitted_weights_raw ) - true_beta[4]
    fitted_weights_raw <- lag_splines %*% t( all_beta[, pm10_inds] )
    pm10_diff <- colSums( fitted_weights_raw ) - true_beta[5]
    fitted_weights_raw <- lag_splines %*% t( all_beta[, o3_inds] )
    o3_diff <- colSums( fitted_weights_raw ) - true_beta[6]
    
    all_diffs <- data.frame( static_diff, rhum_diff, pm10_diff, o3_diff )
    colnames( all_diffs ) <- c( 'Intercept', 'Static1', 'Static2', 'Humidity', 
                                'PM10', 'O3' )
    all_diffs
    
  }
  
  beta_plot <- function( i, dat ) {
    boxplot( beta_extract( i, dat ), col = 'blue' )
    abline( h = 0, lty = 2, col = 'orange' )
  }
  
  pdf( paste0('beta_sims_', wd, '.pdf'), width = 12 )
  beta_plot( 1, sim_res )
  dev.off()
  
  pdf( paste0('all_beta_sims_', wd, '.pdf'), width = 12, height = 12 )
  par( mfrow = c(2, 2) )
  beta_plot( 1, sim_res )
  beta_plot( 2, sim_res )
  beta_plot( 3, sim_res )
  beta_plot( 4, sim_res )
  par( mfrow = c(1, 1) )
  dev.off()
  
  pdf( paste0('beta_sims_small_', wd, '.pdf'), width = 12 )
  beta_plot( 1, small_sim_res )
  dev.off()
  
  pdf( paste0('all_beta_sims_small_', wd, '.pdf'), width = 12, height = 12 )
  par( mfrow = c(2, 2) )
  beta_plot( 1, small_sim_res )
  beta_plot( 2, small_sim_res )
  beta_plot( 3, small_sim_res )
  beta_plot( 4, small_sim_res )
  par( mfrow = c(1, 1) )
  dev.off()
  
  ### DLFs
  beta_extract2 <- function( i, dat ) {
    all_beta <- lapply( dat[[i]]$MCMC, beta_extract_within ) %>% do.call( what = 'rbind' )
  }
  
  DLF_plot <- function( i, dat ) {
    
    nealmon <- function( delta_t, theta ) {
      
      p <- length( theta )
      lw <- matrix(rep(delta_t, each = p) ^ (1:p), ncol = p, byrow = TRUE) %*% theta
      exp(lw - matrixStats::logSumExp(lw))
      
    }
    
    lag_len <- 40
    lag_profile1 <- nealmon( 0:lag_len, c(0.32, -0.02) )
    lag_profile2_raw <- lag_len - (0:lag_len)
    lag_profile2 <- lag_profile2_raw / sum(lag_profile2_raw)
    lag_profile3 <- 0.7*nealmon( 0:lag_len, c(0.2, -0.03) ) + 0.3*nealmon( 0:lag_len, c(6, -0.1) )
    
    rhum_inds <- 4:10
    pm10_inds <- 11:17
    o3_inds <- 18:24
    
    lag_splines <- splines::bs( 0:lag_len, df = 7, intercept = TRUE )
    all_bet <- beta_extract2( i, dat )
    rhum_bet <- all_bet[, rhum_inds]
    pm10_bet <- all_bet[, pm10_inds]
    o3_bet <- all_bet[, o3_inds]
    
    plot_fun <- function( x, lp, nm, profile_lim = FALSE ) {
      
      fitted_weights_raw <- lag_splines %*% t( x )
      fitted_weights <- t( t(fitted_weights_raw) / colSums( fitted_weights_raw ) )
      fitted_weights_raw_mean <- apply( fitted_weights_raw, 1, mean )
      mean_effect <- sum( fitted_weights_raw_mean )
      fitted_weights_mean <- fitted_weights_raw_mean / mean_effect
      
      cat( paste0(nm, ': ', round(mean_effect, 2), '\n') )
      
      if ( profile_lim ) {
        ylim <- range( lp )
      }
      else {
        ylim <- range( fitted_weights, lp, na.rm = TRUE )
      }
      
      plot( 0, ylim = ylim, xlim = c(0, lag_len), type = 'n', ylab = '', xlab = 'Lag', main = nm )
      grid()
      for ( i in 1:ncol(fitted_weights) ) {
        lines( fitted_weights[, i], col = 'grey' )
      } 
      lines( fitted_weights_mean, lwd = 3 )
      lines( lp, col = 'blue', lwd = 2, lty = 2 )
      
    }
    
    par( mfrow = c(1, 3) )
    plot_fun( rhum_bet, lag_profile1, 'Relative Humidity' )
    plot_fun( pm10_bet, lag_profile2, 'PM10', profile_lim = TRUE )
    plot_fun( o3_bet, lag_profile3, 'O3' )
    par( mfrow = c(1, 1) )
    
  }
  
  for ( i in 1:4 ) {
    pdf( paste0('DLF_sims_large', i, '_', wd, '.pdf'), width = 12 )
    DLF_plot( i, sim_res )
    dev.off()
  }
  
  for ( i in 1:4 ) {
    pdf( paste0('DLF_sims_small', i, '_', wd, '.pdf'), width = 12 )
    DLF_plot( i, small_sim_res )
    dev.off()
  }
  
  if( wd == 'QB' ) {
    setwd( '..' )
    return()
  }
  
  ### Xi
  xi_extract <- function( i, dat ) {
    lapply( dat[[i]]$MCMC, '[[', 'xi' )
  }
  
  xi_plot <- function( i, dat ) {
    xi_dat <- xi_extract( i, dat ) %>% do.call( what = 'rbind' ) %>% t %>% as.data.frame
    x_range <- range( xi_dat )
    xi_densities <- lapply( xi_dat, density )
    y_range <- lapply( xi_densities, '[[', 'y' ) %>% range
    
    plot( 0, type = 'n', xlim = x_range, ylim = y_range, yaxt = 'n',
          xlab = 'Stopping Parameter', ylab = '', main = 'Stopping Parameter Posterior Density' )
    for ( i in 1:length(xi_densities) ) {
      lines( xi_densities[[i]], col = adjustcolor('grey', 0.75) )
    }
    abline( v = 50, col = 'blue', lty = 2 )
  }
  
  pdf( paste0('xi_sims_', wd, '.pdf'), width = 12 )
  xi_plot( 1, sim_res )
  dev.off()
  
  pdf( paste0('all_xi_sims_', wd, '.pdf'), width = 12, height = 12 )
  par( mfrow = c(2, 2) )
  xi_plot( 1, sim_res )
  xi_plot( 2, sim_res )
  xi_plot( 3, sim_res )
  xi_plot( 4, sim_res )
  par( mfrow = c(1, 1) )
  dev.off()
  
  pdf( paste0('xi_sims_small_', wd, '.pdf'), width = 12 )
  xi_plot( 1, small_sim_res )
  dev.off()
  
  pdf( paste0('all_xi_sims_small_', wd, '.pdf'), width = 12, height = 12 )
  par( mfrow = c(2, 2) )
  xi_plot( 1, small_sim_res )
  xi_plot( 2, small_sim_res )
  xi_plot( 3, small_sim_res )
  xi_plot( 4, small_sim_res )
  par( mfrow = c(1, 1) )
  dev.off()
  
  setwd( '..' )
  
}

### Run output
output_fun( 'NB' )
output_fun( 'QB' )

