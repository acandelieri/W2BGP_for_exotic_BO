
# Forrester function: 1D, 4 fidelities
forrester <- function( x, fidelity ) {
  
  stopifnot( fidelity %in% 1:4 )
  
  # original search space: [0,1]
  
  if( fidelity == 2 ) {
    return( ( 5.5*x - 2.5 )^2 * sin( 12*x - 4 ) ) # fidelity == 2
  } else {
    res <- ( 6*x - 2 )^2 * sin( 12*x - 4 ) # fidelity == 1
    if( fidelity==3 ) {
      res <- 0.75*res + 5*(x-0.5) - 2 # fidelity == 3 
    } else {
      if( fidelity == 4 ) {
        res <- 0.5*res + 10*(x-0.5) - 5 # fidelity == 4
      }
    }
    return( res )
  }
}


# Rosenbrock function: nD with n>=2, 3 fidelities
rosenbrock <- function( x, fidelity ) {
  
  stopifnot( length(x)>=2 )
  stopifnot( fidelity %in% 1:3 )
  
  # rescaling to the original search space, that is [-2,2]^d
  x <- x * 4 - 2
  
  if( fidelity == 2 ) {
    return( sum(50*(x[-1] - x[-length(x)]^2)^2) + sum((-2-x)^2) - sum(0.5*x) ) # fidelity == 2
  } else {
    res <- sum(100*(x[-1]-x[-length(x)]^2)^2) + sum((1-x)^2) # fidelity == 1
    if( fidelity == 3 ) 
      res <- (res - 4 - sum(0.5*x)) / (10 + sum(0.25*x)) # fidelity == 3
    return( res )
  }
}

# Shifted-Rotated Rastriging: any d, 3 fidelities
shiftedRotatedRastrigin <- function( x, fidelity, phi1=10000, phi2=5000, phi3=2500 ) {
  
  stopifnot( fidelity %in% 1:3 )
  
  if( fidelity==1 ) {
    phi <- phi1
  } else {
    if( fidelity==2 ) {
      phi <- phi2
    } else {
      phi <- phi3
    }
  }
  
  # rescaling to the original search space, that is [-0.1,0.2]^d
  x <- x * 0.3 - 0.1

  z <- matrix( c(cos(0.2),-sin(0.2),sin(0.2),cos(0.2)),2,2,byrow=T ) %*% ( x - rep(0.1,length(x)) )
  a <- 1 - 0.0001*phi
  er <- sum( a * (cos( 10*pi*a*z + 0.5*pi*a + pi ))^2 )
  
  return( sum( z^2 + 1 - cos(10*pi*z) ) + er )
  
}


# Heterogeneous function: nD with n>=2, 2 fidelities
heterogeneous <- function( x, fidelity ) {
  
  stopifnot( fidelity %in% 1:2 )
  
  # original search space: [0,1]^d
  
  if( length(x)==1 ) {
    # 1D
    res <- sin(30*(x-0.9)^4) * cos(2*(x-0.9)) + 0.5*(x-0.9) # fidelity == 1
    if( fidelity == 2 )
      res <- (res-1+x) / (1+0.25*x) # fidelity == 2 
  } else {
    # 2D
    res <- sin(21*(x[1]-0.9)^4) * cos(2*(x[1]-0.9)) + 0.5*(x[1]-0.7) + (x[1]*sin(x[1])+2*x[2]^2*sin(prod(x)))
    if( fidelity== 2 )
      res <- (res - 2 + sum(x)) / (5 + sum(0.25*x) )
  }
  return(res)
}


# Paciorek function: nD with n>=2, 2 fidelities
paciorek <- function( x, fidelity, A=0.5, sig1=0.00125, sig2=0.075 ) {
  
  stopifnot( fidelity %in% 1:2 )
  
  # rescaling to the original search space, that is [0.3,1]^d
  x <- x * 0.7 + 0.3
  
  res <- sin(1/prod(x)) + rnorm(1,0,sig1) # fidelity == 1
  if( fidelity == 2 )
    res <- res - (9*A^2)*cos(1/prod(x)) + rnorm(1,0,sig2) # fidelity 2
  
  return( res )
}