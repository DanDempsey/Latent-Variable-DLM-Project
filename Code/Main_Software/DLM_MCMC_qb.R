##### Quantile Binary DLM with Reversible-Jump-like Variable Selection 
##### Daniel Dempsey

### Load libraries
library( dlnm ) # For extra DLM functionality 
library( magrittr ) # For piping operator %>%
library( statmod ) # For Inverse Gaussian
library( mvtnorm ) # For mu
library( pROC ) # For ROC curves
source( 'Code/Main_Software/ALD.R' ) # ALD functions

### Alternative to base sample function; doesn't treat scalars differently to vectors, unlike original
sample2 <- function( x, size, replace = FALSE, prob = NULL ) {
  if (missing(size)) {
    size <- length(x)
  } 
  x[sample.int(length(x), size, replace, prob)]
}

### Function for creating default priors
defaultPrior_qb <- function( X ) {
  
  # Setup
  nvar <- ncol( X )
  prior <- list()
  
  # beta
  prior$beta$b <- rep( 0, nvar )
  prior$beta$v <- diag( 100, nvar )
  
  # Gamma
  prior$gamma <- c( 1, rep(0.5, nvar - 1) )
  
  # Return result
  prior
  
}

### Function for creating default starting values for the MCMC algorithm
defaultInit_qb <- function( X, group_len ) {
  
  # Setup
  nvar <- ncol( X )
  init <- list()
  
  # beta
  init$beta <- rep( 0, nvar )
  
  # gamma
  init$gamma <- c( TRUE, rep(FALSE, group_len - 1) )
  
  # Return result
  init
  
}

### Main wrapper function for performing MCMC-based inference for the Negative Binomial DLM
DLM_MCMC_qb <- function( formula, data = NULL, quantile = 0.5, groups = NULL, 
                         prior = NULL, init = NULL, nsamp = 5000, nburn = 5000, thin = 1 ) {
  
  ### Initialize
  cat( "Initializing MCMC algorithm...\n" )
  MCMC_length <- nsamp + nburn
  
  # Set dataset and groups
  X_full <- model.matrix( formula, data )
  if ( is.null(groups) ) { groups <- 1:ncol(X_full) }
  group_len <- max( groups )
  
  # Prepare priors and associated statistics used in the algorithm
  if ( is.null(prior) ) {
    prior <- defaultPrior_qb( X_full )
  }
  
  prior_gamma <- prior$gamma
  
  V0i_full <- chol2inv( chol(prior$beta$v) )
  V0ib0_full <- V0i_full %*% prior$beta$b
  
  # Prepare starting values
  if ( is.null(init) ) {
    init <- defaultInit_qb( X_full, group_len )
  }
  if ( !init$gamma[1] ) {
    init$gamma[1] <- TRUE
    warning( 'The first element of the gamma starting value must be TRUE. This has been corrected.' )
  }
  
  # Parameter initialization
  nvar <- ncol( X_full )
  betares <- matrix( 0, ncol = nvar, nrow = MCMC_length ) # Regression slopes
  betares[1, ] <- init$beta 
  gammares <- matrix( FALSE, ncol = nvar, nrow = MCMC_length ) # Variable inclusion indicator (starting value set in the next block of code)
  colnames( betares ) <- colnames( gammares ) <- colnames( X_full )
  
  # Variable index to link the dynamic variable gammas and apply starting value
  var_split <- split( 1:nvar, groups )
  var_index_unique <- unique( groups )[-1] # Used when proposing new gamma values
  gammares[1, ][unlist(var_split[init$gamma])] <- TRUE
  gam <- gammares[1, ]
  #var_sum <- length( unique(groups[gam]) )
  
  X <- X_full[, gam, drop = FALSE]
  Xb <- X %*% init$beta[gam]
  V0i <- V0i_full[gam, gam]
  V0ib0 <- V0ib0_full[gam]
  
  y <- data[[all.vars(formula)[1]]]
  y_len <- length( y )
  y_max <- max( y ) + 1
  rtrunc <- ifelse( y, TRUE, FALSE )
  
  # Starting values for latent variables
  psi <- ( 1 - 2*quantile ) / ( quantile * (1 - quantile) )
  phi <- 2 / ( quantile * (1 - quantile) )
  delta <- 2 + ( ( psi^2 ) / phi )
  
  ystar <- rTALD( n = y_len, rtrunc = rtrunc, mu = Xb, sigma = 1, p = quantile )
  chi <- (ystar - Xb)^2 / phi
  nu <- 1/rinvgauss( y_len, mean = sqrt( delta/chi ), shape = delta )
  lambda <- ystar - (psi * nu)
  
  omega <- 1 / (phi * nu)
  XtO <- t(X * omega)
  
  V <- chol2inv( chol(V0i + XtO%*%X) )
  B <- V%*%( V0ib0 + XtO%*%lambda )
  
  ### Main Loop
  cat( 'Initialization complete. Running algorithm...\n' )
  for ( i in 2:MCMC_length ) {
    
    ### Update gamma
    # Propose change to parameter inclusion set
    change_ind <- var_split[[ sample2( var_index_unique, 1 ) ]]
    gam_star <- gam
    gam_star[change_ind] <- !gam_star[change_ind]
    
    # Construct log acceptance ratio for proposed move
    X_star <- X_full[, gam_star, drop = FALSE]
    V0i_star <- V0i_full[gam_star, gam_star]
    V0ib0_star <- V0ib0_full[gam_star]
    XtO_star <- t( X_star * omega )
    
    V_star <- chol2inv( chol(V0i_star + XtO_star%*%X_star) )
    B_star <- V_star%*%( V0ib0_star + XtO_star%*%lambda )
    
    ldet_V <- sum( log(diag(chol(V))) )
    ldet_V_star <- sum( log(diag(chol(V_star))) )
    ldet_V0i <- sum( log( diag(chol(V0i)) ) )
    ldet_V0i_star <- sum( log(diag(chol(V0i_star))) )
    
    gam_lprior <- dbinom( gam, 1, prior_gamma, log = TRUE )
    gam_lprior_star <- dbinom( gam_star, 1, prior_gamma, log = TRUE )
    
    lkernel <- crossprod(B, chol2inv(chol(V)))%*%B/2
    lkernel_star <- crossprod(B_star, chol2inv(chol(V_star)))%*%B_star/2
    
    lnum <- sum( ldet_V_star, ldet_V0i_star, lkernel_star, gam_lprior_star )
    ldenom <- sum( ldet_V, ldet_V0i, lkernel, gam_lprior )
    
    # Accept proposed update to gamma with M-H acceptance probability
    if ( (lnum - ldenom) > log(runif(1)) ) { # accept
      gammares[i, ] <- gam <- gam_star
      #var_sum <- ifelse( gam[change_ind], var_sum + 1, var_sum - 1 )
      X <- X_star
      V0i <- V0i_star
      V0ib0 <- V0ib0_star
      Xb <- X %*% betares[i-1, gam]
    } else { # reject
      gammares[i, ] <- gammares[i-1, ]
    }
    
    ### Update latent parameters
    # ystar
    ystar <- rTALD( n = y_len, rtrunc = rtrunc, mu = Xb, sigma = 1, p = quantile )
    
    # nu
    chi <- ( ystar - Xb )^2 / phi
    nu <- 1/rinvgauss( y_len, mean = sqrt(delta/chi), shape = delta )
    
    ### Update beta
    lambda <- ystar - ( psi * nu )
    omega <- 1 / ( phi * nu )
    XtO <- t( X * omega )
    V <- chol2inv( chol(V0i + XtO%*%X) )
    B <- V%*%( V0ib0 + XtO%*%lambda )
    betares[i, gam] <- B + t(chol(V)) %*% rnorm(sum(gam))
    Xb <- X %*% betares[i, gam]
    
    ### Print a progress report every 500 iterations
    if ( !i%%1000 ) {
      cat( paste0( 'Completed iteration ', i, '...\n' ) )
    }
    
  }
  
  ### Filter Markov chain and return result
  cat( 'Algorithm complete. Returning result.\n' )
  
  keep <- seq( nburn + 1, MCMC_length, thin )
  col_inds <- sapply(var_split, '[', 1)
  gamma_trunc <- gammares[keep, col_inds]
  
  list( beta = betares[keep, ], gamma = gammares[keep, ], quantile = quantile, 
        init = init, prior = prior, groups = groups, X = X_full, y = y )
  
}

### Post predictive distribution that binary response equals 1 given new data
PPD_qr <- function( formula, data, res ) {
  
  browser()
  ### Design matrices
  ynew <- data[, all.vars(formula)[[1]]]
  X_full <- model.matrix( formula, data )
  Xb_test <- X_full %*% t( res$beta )
  
  ### Compute posterior predictive probabilities and ROC
  test_p <- pALD( Xb_test, p = res$quantile )
  ppp <- apply( test_p, 1, mean )
  ROC <- roc( ynew, ppp )
  
  ### Return results
  list( predicted = ppp, actual = ynew, ROC = ROC )
  
}


