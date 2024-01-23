##### Negative-Binomial DLM with Reversible-Jump-like Variable Selection 
##### Daniel Dempsey

### Load libraries
library( BayesLogit ) # For sampling Polya-Gamma random variables
library( dlnm ) # For extra DLM functionality 
library( magrittr ) # For piping operator %>%
dyn.load( 'Code/Main_Software/C_funs/PSI_FUN.so' ) # C function used when updating the NB stopping parameter 

### Alternative to base sample function; doesn't treat scalars differently to vectors, unlike original
sample2 <- function( x, size, replace = FALSE, prob = NULL ) {
  if (missing(size)) {
    size <- length(x)
  } 
  x[sample.int(length(x), size, replace, prob)]
}

### Function for creating default priors
defaultPrior_nb <- function( X ) {
  
  # Setup
  nvar <- ncol( X )
  prior <- list()
  
  # beta
  prior$beta$b <- rep( 0, nvar )
  prior$beta$v <- diag( 100, nvar )
  
  # Gamma
  prior$gamma <- c( 1, rep(0.5, nvar - 1) )
  
  # xi
  prior$xi$a0 <- 2
  prior$xi$b0 <- 1/50
  
  # Return result
  prior
  
}

### Function for creating default starting values for the MCMC algorithm
defaultInit_nb <- function( X, group_len ) {
  
  # Setup
  nvar <- ncol( X )
  init <- list()
  
  # beta
  init$beta <- rep( 0, nvar )
  
  # xi
  init$xi <- 1
  
  # gamma
  init$gamma <- c( TRUE, rep(FALSE, group_len - 1) )
  
  # Return result
  init
  
}

### Main wrapper function for performing MCMC-based inference for the Negative Binomial DLM
DLM_MCMC_nb <- function( formula, data = NULL, groups = NULL, prior = NULL, init = NULL, 
                         nsamp = 5000, nburn = 5000, thin = 1 ) {
  
  ### Initialize
  cat( "Initializing MCMC algorithm...\n" )
  MCMC_length <- nsamp + nburn
  
  # Set dataset and groups
  X_full <- model.matrix( formula, data )
  if ( is.null(groups) ) { groups <- 1:ncol(X_full) }
  group_len <- max( groups )
  
  # Prepare priors and associated statistics used in the algorithm
  if ( is.null(prior) ) {
    prior <- defaultPrior_nb( X_full )
  }
  
  prior_xi_a0 <- prior$xi$a0
  prior_xi_b0 <- prior$xi$b0
  
  prior_gamma <- prior$gamma
  
  V0i_full <- chol2inv( chol(prior$beta$v) )
  V0ib0_full <- V0i_full %*% prior$beta$b
  
  # Prepare starting values
  if ( is.null(init) ) {
    init <- defaultInit_nb( X_full, group_len )
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
  xires <- numeric( MCMC_length ) # Negative-Binomial overdispersion parameter
  xires[1] <- xi <- init$xi
  colnames( betares ) <- colnames( gammares ) <- colnames( X_full )
  
  # Variable index to link the dynamic variable gammas and apply starting value
  var_split <- split( 1:nvar, groups )
  var_index_unique <- unique( groups )[-1] # Used when proposing new gamma values
  gammares[1, ][unlist(var_split[init$gamma])] <- TRUE
  gam <- gammares[1, ]
  var_sum <- length( unique(groups[gam]) )
  
  X <- X_full[, gam, drop = FALSE]
  Xb <- X %*% init$beta[gam]
  V0i <- V0i_full[gam, gam]
  V0ib0 <- V0ib0_full[gam]
  
  y <- data[[all.vars(formula)[1]]]
  y_len <- length( y )
  y_max <- max( y ) + 1
  
  # Starting value for Poyla-Gamma latent parameter and other associated statistics
  omega <- rpg( y_len, y + xi, Xb )
  XtO <- t( X * omega )
  lambda <- (y - xi) / (2 * omega)
  V <- chol2inv( chol(V0i + XtO%*%X) )
  B <- V%*%( V0ib0 + XtO%*%lambda )
  
  # Initialise matrix of pmfs for psi inference
  R <- matrix( 0, nrow = y_max, ncol = y_max )
  R[1, 1] <- R[2, 2] <- 1
  
  ### Main Loop
  cat( 'Initialization complete. Running algorithm...\n' )
  for ( i in 2:MCMC_length ) {
    
    ### Update xi
    # Poisson latent parameter: psi
    sum_psi <- .C( 'PSI_FUN', xi, as.integer(y), as.integer(y_len), as.integer(0:(y_max-1)), 
                   as.integer(y_max), R, as.integer(0) )[[7]]
    
    # Negative Binomial overdispersion parameter: xi
    xires[i] <- xi <- rgamma( 1, prior_xi_a0 + sum_psi, prior_xi_b0 + sum(log(1 + exp(Xb))) )
    
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
      var_sum <- ifelse( gam[change_ind], var_sum + 1, var_sum - 1 )
      X <- X_star
      V0i <- V0i_star
      V0ib0 <- V0ib0_star
      Xb <- X %*% betares[i-1, gam]
    } else { # reject
      gammares[i, ] <- gammares[i-1, ]
    }
    
    ### Update beta
    # Polya-Gamma latent parameter: omega
    omega <- rpg( y_len, y + xi, Xb )
    
    # Regression slopes: beta
    XtO <- t(X * omega)
    lambda <- (y - xi) / (2 * omega)
    V <- chol2inv(chol(V0i + XtO%*%X))
    B <- V%*%(V0ib0 + XtO%*%lambda)
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
  
  list( beta = betares[keep, ], gamma = gammares[keep, ], xi = xires[keep], 
        init = init, prior = prior, groups = groups, X = X_full, y = y )
  
}

