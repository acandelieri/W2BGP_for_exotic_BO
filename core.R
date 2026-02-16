library(lhs)
library(DiceKriging)






next_federated_query <- function( observations, agents_kernels, agent_weight, acq, acq_par, nStarts=5 ) {
  
  if( acq!="lcb" ) stop("NOT YET IMPLEMENTED!")
  stopifnot( is.vector(acq_par) & length(acq_par)==1 & acq_par>=0 )
  
  M <- length(agents_kernels)
  # d <- ncol(observations) - 6 # 'd' is retrieved by the number of columns in observations!
  # offset <- 4
  d <- ncol(observations) - 7 # 'd' is retrieved by the number of columns in observations!
  offset <- 5

  gps <- list()
  nuggets <- rep(F,M)
  
  for( m in 1:M ) {
    
    X <- observations[observations$agent==m & !is.na(observations$iter),offset+(0:(d-1)),drop=F]
    y <- observations$y[observations$agent==m & !is.na(observations$iter)]
    
    gp <- NULL
    # first we try to fit a GP without noise (because test functions are deterministic)...
    try( expr=(gp <- km( design=X, response=y, covtype=agents_kernels[m],
                         nugget.estim=F, control=list(trace=0) ) ),
         silent=T )
    # ...if it is not possible, then we estimate a nugget effect
    while( is.null(gp) ) {
      nuggets[m] <- T
      try( expr=(gp <- km( design=X, response=y, covtype=agents_kernels[m],
                           nugget.estim=T, control=list(trace=0) ) ),
           silent=T ) 
    }
    gps[[m]] <- gp
  }
  
  X <- NULL
  if( agent_weight != 1/M ) {
    for( m in 1:M ) {
      lambdas <- rep( (1-agent_weight)/(M-1), M )
      lambdas[m] <- agent_weight
      x_ <- solve_w2bgp_acquisition( gps=gps, lambdas=lambdas, acq=acq, acq_par=acq_par, nStarts=nStarts )
      X <- rbind( X, x_ )
    }
  } else { # just once!
    lambdas <- rep(1/M,M)
    x_ <- solve_w2bgp_acquisition( gps=gps, lambdas=lambdas, acq=acq, acq_par=acq_par, nStarts=nStarts )
    for( m in 1:M ) 
      X <- rbind( X, x_ )
  }
  rownames(X) <- NULL
  
  return( list( X=X, nuggets=nuggets ) )
}



next_batch_query <- function( observations, agents_kernels, agent_weight, acq, acq_par, nStarts=5 ) {
  
  if( acq!="lcb" ) stop("NOT YET IMPLEMENTED!")
  stopifnot( is.vector(acq_par) & length(acq_par)==1 & acq_par>=0 )
  
  M <- length(agents_kernels)
  # d <- ncol(observations) - 5 # 'd' is retrieved by the number of columns in observations!
  # offset <- 4
  d <- ncol(observations) - 7 # 'd' is retrieved by the number of columns in observations!
  offset <- 5
  
  tmp <- cbind( observations[!is.na(observations$iter),offset+(0:(d-1)),drop=F], y=observations$y[!is.na(observations$iter)] )
  tmp <- as.matrix(unique(tmp))
  colnames(tmp) <- NULL
  
  
  gps <- list()
  nuggets <- rep(F,M)
  
  for( m in 1:M ) {
    gp <- NULL
    # first we try to fit a GP without noise (because test functions are deterministic)...
    try( expr=(gp <- km( design=data.frame(x=tmp[,-ncol(tmp),drop=F]), response=tmp[,ncol(tmp)], covtype=agents_kernels[m],
                         nugget.estim=F, control=list(trace=0) ) ),
         silent=T )
    # ...if it is not possible, then we estimate a nugget effect
    while( is.null(gp) ) {
      nuggets[m] <- T
      try( expr=(gp <- km( design=data.frame(x=tmp[,-ncol(tmp),drop=F]), response=tmp[,ncol(tmp)], covtype=agents_kernels[m],
                           nugget.estim=T, control=list(trace=0) ) ),
           silent=T ) 
    }
    gps[[m]] <- gp
  }
  
  
  X <- NULL
  if( agent_weight != 1/M ) {
    for( m in 1:M ) {
      lambdas <- rep( (1-agent_weight)/(M-1), M )
      lambdas[m] <- agent_weight
      x_ <- solve_w2bgp_acquisition( gps=gps, lambdas=lambdas, acq=acq, acq_par=acq_par, nStarts=nStarts )
      X <- rbind( X, x_ )
    }
  } else { # just once!
    lambdas <- rep(1/M,M)
    x_ <- solve_w2bgp_acquisition( gps=gps, lambdas=lambdas, acq=acq, acq_par=acq_par, nStarts=nStarts )
    for( m in 1:M ) 
      X <- rbind( X, x_ )
  }
  rownames(X) <- NULL
  
  return( list( X=X, nuggets=nuggets ) )
}



# Multi-fidelity
next_multifidelity_query <- function( observations, kernel, fidelityCosts, acq_par, nStarts=5 ) {
  
  stopifnot( is.vector(acq_par) & length(acq_par)==1 & acq_par>=0 )
  
  M <- length(fidelityCosts)
  
  d <- ncol(observations) - 6 # 'd' is retrieved by the number of columns in observations!
  offset <- 4
  
  gps <- list()
  nuggets <- rep(F,M)
  
    
  for( m in 1:M ) {
    
    X <- observations[observations$source==m & !is.na(observations$iter),offset+(0:(d-1)),drop=F]
    y <- observations$y[observations$source==m & !is.na(observations$iter)]
    
    gp <- NULL
    # first we try to fit a GP without noise (because test functions are deterministic)...
    try( expr=(gp <- km( design=X, response=y, covtype=kernel,
                         nugget.estim=F, control=list(trace=0) ) ),
         silent=T )
    # ...if it is not possible, then we estimate a nugget effect
    while( is.null(gp) ) {
      nuggets[m] <- T
      try( expr=(gp <- km( design=X, response=y, covtype=kernel,
                           nugget.estim=T, control=list(trace=0) ) ),
           silent=T ) 
    }
    gps[[m]] <- gp
  }
    

  sol <- solve_w2bMFlcb_acquisition( gps=gps, fidelityCosts=fidelityCosts, acq_par=acq_par, nStarts=nStarts )
  
  return( list( x=sol$x, source_ix=sol$s, nuggets=nuggets ) )
}



solve_w2bgp_acquisition <- function( gps, lambdas, acq, acq_par, nStarts ) {
  d <- ncol(gps[[1]]@X)
  bestSol <- NULL
  for( i in 1:nStarts ) {
    par0 <- runif(d)
    resOpt <- optim( par=par0, fn=wlcb, gr=NULL, method="L-BFGS-B", lower=0, upper=1, control=list(trace=0),
                     gps=gps, lambdas=lambdas, beta=acq_par )
    if( is.null(bestSol) || resOpt$value<bestSol$value )
      bestSol <- resOpt
  }
  return( bestSol$par )
}



# Multi-fidelity
solve_w2bMFlcb_acquisition <- function( gps, fidelityCosts, acq_par, nStarts ) {
  
  d <- ncol(gps[[1]]@X)
  
  # the groudtruth cannot be queried less than all the others!!!
  queries_S1 <- nrow(gps[[1]]@X)
  queries_S2 <- 0
  for( i in 2:length(gps) ) {
    queries_S2 <- queries_S2 + nrow(gps[[i]]@X)
  }
  nSources <- ifelse( queries_S1 <= queries_S2, 1, length(gps) ) 
    
  bestSol <- bestSource <- NULL
  for( i in 1:nStarts ) {
    par0 <- runif(d)
    for( s in 1:nSources ) {
      resOpt <- optim( par=par0, fn=mfwlcb, gr=NULL, method="L-BFGS-B", lower=0, upper=1, control=list(trace=0),
                       gps=gps, fidelityCosts=fidelityCosts, beta=acq_par, source_ix=s )
      if( is.null(bestSol) || resOpt$value<bestSol$value ) {
        bestSol <- resOpt
        bestSource <- s
      }
    }
  }
  return( list( x=bestSol$par, s=bestSource ) )
}



wgp.predict <- function( x, gps, lambdas ) {
  
  stopifnot( is.vector(x) )
  
  if( length(x)>1 )
    x <- t(x)
    
  if( !is.list(gps) )
    gps <- list(gps)
  M <- length(gps)
  
  stopifnot( length(lambdas)==M )
  stopifnot( sum(lambdas)==1 )

  mu <- sg <- 0
  for( m in 1:M ) {
    pred <- predict( gps[[m]], data.frame(x=x), "UK" )
    mu <- mu + lambdas[m] * pred$mean
    sg <- sg + lambdas[m] * pred$sd
  }
  
  mu <- mu/M
  sg <- sg/M
  
  return( list( mean=mu, sd=sg ) )
}



wlcb <- function( x, gps, lambdas, beta ) {
  
  pred <- wgp.predict( x=x, gps=gps, lambdas=lambdas )
  
  return( pred$mean - beta * pred$sd )
}



# Multi-fidelity
mfwlcb <- function( x, gps, fidelityCosts, beta, source_ix ) {
  
  # WBGP
  lambdas <- fidelityCosts/sum(fidelityCosts)
  pred <- wgp.predict( x=x, gps=gps, lambdas=lambdas )
  
  # Acquisition for source "source_ix"
  y.best <- min(gps[[1]]@y) # always from ground-truth...
  numerator <- (y.best - (pred$mean - beta * pred$sd))
  pred_s <- predict( gps[[source_ix]], data.frame(x=t(x)), "UK" )
  value <- numerator / ( fidelityCosts[source_ix] * ((pred_s$mean-pred$mean)^2 - (pred_s$sd-pred$sd)^2) )

  return( value )
}