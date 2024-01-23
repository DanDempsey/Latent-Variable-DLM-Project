##### Sampler for a truncated ALD, based on ald package
##### Daniel Dempsey

pALD <- function(q, mu = 0, sigma = 1, p = 0.5) {
  ifelse( test = q < mu, 
          yes = p * exp((1 - p) * (q - mu)/sigma), 
          no = 1 - (1 - p) * exp(-p * (q - mu)/sigma) )
}

qALD <- function(prob, mu = 0, sigma = 1, p = 0.5)  {
  ifelse( test = prob < p, 
          yes = mu + (sigma*log(prob/p))/(1 - p), 
          no = mu - sigma*log((1 - prob)/(1 - p))/p )
}

# Used elsewhere
rALD <- function(n, mu = 0, sigma = 1, p = 0.5) {
  u <- runif(n)
  mapply(qALD, prob = u, mu = mu, sigma = sigma, p = p)
}

# Not used anywhere, but included here for the sake of completeness (based on ald package)
dALD <- function(y, mu = 0, sigma = 1, p = 0.5) {
  ifelse(test = y < mu, 
         yes = (p * (1 - p)/sigma) * exp((1 - p) * (y - mu)/sigma), 
         no = (p * (1 - p)/sigma) * exp(-p * (y - mu)/sigma))
}

rTALD <- function(n, rtrunc, mu = 0, sigma = 1, p = 0.5) {
  bound <- pALD( 0, mu = mu, sigma = sigma, p = p )
  u <- runif( n, min = ifelse( rtrunc, bound, 0 ), max = ifelse( rtrunc, 1, bound ) )
  qALD( prob = u, mu = mu, sigma = sigma, p = p )
}

