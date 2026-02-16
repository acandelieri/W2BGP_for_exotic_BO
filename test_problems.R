#**************************************************************************************
# 1D test functions (defined in [0,1])
#**************************************************************************************

problem_02 <- function(x) {
  x <- x * (7.5-2.7) + 2.7
  f <- sin(x) + sin( x*10/3)
  return( f )
}

problem_03 <- function(x) {
  x <- x * 20 - 10
  f <- 0
  for( k in 0:5 )
    f <- f + ( k * sin( (k+1)*x + k ) )
  return( -f ) 
}

problem_05 <- function(x) {
  x <- x * 1.2
  f <- -(1.4 - 3*x)*sin(18*x)
  return( f ) 
}

problem_07 <- function(x) {
  x <- x * (7.5-2.7) + 2.7
  f <- sin(x ) + sin(x*10/3) + log(x) - 0.84*x + 3
  return( f )
}

problem_11 <- function(x) {
  x <- x * ( 2*pi + pi/2) -pi/2
  f <- 2*cos(x) + cos(2*x)
  return( f )
}

problem_14 <- function(x) {
  x <- x*4
  f <- -exp(-x) * sin(2*pi*x)
  return( f )
}

problem_15 <- function(x) {
  x <- x * (5 + 5) - 5
  f <- (x^2 - 5*x + 6)/(x^2 + 1)
  return( f )
}

problem_22 <- function(x) {
  x <- x * 20
  f <- exp( -3*x ) - (sin(x))^3
  return(f)
}


#**************************************************************************************
# 2D test functions (defined in [0,1]x[0,1])
#**************************************************************************************

bird <- function(x) {
  x <- x*(2*pi-(-2*pi)) - 2*pi
  y <- (x[1]-x[2])^2 + exp( (1-sin(x[1]))^2 ) * cos(x[2]) + exp( (1-cos(x[2]))^2 ) * sin(x[1])
  return(y)
}

michalewicz <- function(x) {
  m <- 10
  x <- x*pi
  ii <- c(1:length(x))
  y <- sum(sin(x) * (sin((ii*x^2)/pi))^(2*m))
  return( -y )
}

ursem03 <- function(x){
  x[1] <- x[1]*(2-(-2)) - 2
  x[2] <- x[2]*(1.5-(-1.5)) - 1.5
  y <- -sin(2.2*pi*x[1]+0.5*pi)*((2-abs(x[1]))/2)*((3-abs(x[1]))/2) - sin(2.2*pi*x[2]+0.5*pi)*((2-abs(x[2]))/2)*((3-abs(x[2]))/2)
  return( y )
}

# da Alroomi Website: https://al-roomi.org/benchmarks/unconstrained/2-dimensions/132-ursem-wave-function
ursemWaves <- function(x){
  x[1] <- x[1]*(1.2-(-0.9)) - 0.9
  x[2] <- x[2]*(1.2-(-1.2)) - 1.2
  y <- -(0.3*x[1])^3 + (x[2]^2 - 4.5*x[2]^2) * x[1]*x[2] + 4.7 * cos(3*x[1] - x[2]^2*(2+x[1])) * sin(2.5*pi*x[1])
  return( y )
}



#**************************************************************************************
# D>2 test functions (defined in [0,1]^d)
#**************************************************************************************

hartmann3 <- function(x) {
  
  alpha <- c(1.0, 1.2, 3.0, 3.2)
  
  A <- c(3.0, 10, 30,
         0.1, 10, 35,
         3.0, 10, 30,
         0.1, 10, 35)
  A <- matrix(A, 4, 3, byrow=TRUE)
  
  P <- 10^(-4) * c(3689, 1170, 2673,
                   4699, 4387, 7470,
                   1091, 8732, 5547,
                   381, 5743, 8828)
  P <- matrix(P, 4, 3, byrow=TRUE)
  
  xxmat <- matrix(rep(x,times=4), 4, 3, byrow=TRUE)
  inner <- rowSums(A[,1:3]*(xxmat-P[,1:3])^2)
  outer <- sum(alpha * exp(-inner))
  
  return( - outer )
}


hartmann6 <- function(x) {
  
  alpha <- c(1.0, 1.2, 3.0, 3.2)
  
  A <- c(10, 3, 17, 3.5, 1.7, 8,
         0.05, 10, 17, 0.1, 8, 14,
         3, 3.5, 1.7, 10, 17, 8,
         17, 8, 0.05, 10, 0.1, 14)
  A <- matrix(A, 4, 6, byrow=TRUE)

  P <- 10^(-4) * c(1312, 1696, 5569, 124, 8283, 5886,
                   2329, 4135, 8307, 3736, 1004, 9991,
                   2348, 1451, 3522, 2883, 3047, 6650,
                   4047, 8828, 8732, 5743, 1091, 381)
  P <- matrix(P, 4, 6, byrow=TRUE)
  
  xxmat <- matrix(rep(x,times=4), 4, 6, byrow=TRUE)
  inner <- rowSums(A[,1:6]*(xxmat-P[,1:6])^2)
  outer <- sum(alpha * exp(-inner))
  
  y <- -(2.58 + outer) / 1.94
  
  return( y )
}


alpine01 <- function(x) {
  x <- x*(10-(-10)) - 10
  y <- sum(abs(x * sin(x) + 0.1 * x))
  return( y )
}


styblinskiTang <- function(x){
  x <- x*(5-(-5)) - 5
  y <- 1/2 * sum(x^4 - 16*x^2 + 5*x)
  return( y )
}



# *********************************************************************************
# Utility functions
# *********************************************************************************

getDim <- function( test_problem ) {
  
  if( test_problem %in% paste0("problem_",c("02","03","05","07","11","14","15","22")) )
    return( 1 )
  if( test_problem %in% c("bird","michalewicz","ursem03","ursemWaves") )
    return( 2 )
  if( test_problem == "hartmann3" )
    return( 3 )
  if( test_problem == "hartmann6" )
    return( 6 )
  
  return( NA )
}


getOptY <- function( test_problem ) {
  
  # 1D test functions
  if( test_problem == "problem_02" ) {
    return(-1.8996)
  }
  if( test_problem == "problem_03" ) {
    return(-12.0312)
  }
  if( test_problem == "problem_05" ) {
    return(-1.4891)
  }
  if( test_problem == "problem_07" ) {
    return(-1.6013)
  }
  if( test_problem == "problem_11" ) {
    return(-1.5)
  }
  if( test_problem == "problem_14" ) {
    return(-0.7887)
  }
  if( test_problem == "problem_15" ) {
    return(-0.03553391)
  }
  if( test_problem == "problem_22" ) {
    return(exp(-27*pi/2)-1)
  }             
  
  # 2D test functions
  if( test_problem == "bird" )
    return(-106.7645)
  if( test_problem == "michalewicz" )
    return(-1.8013)
  if(  strsplit(test_problem,"_",fixed=T)[[1]][1] == "styblinskiTang" ) {
    # return(-78.33233)
    return( as.numeric(gsub("d","",strsplit(test_problem,"_",fixed=T)[[1]][2],fixed=T))*(-39.16599) )
  }
  if( test_problem == "ursem03" )
    return(-3)
  if( test_problem == "ursemWaves" )
    return(-7.306999)
  
  # D>2 test functions
  if( test_problem == "hartmann3" )
    return(-3.862782145)
  if( test_problem == "hartmann6" )
    return(-3.32237)
  if( strsplit(test_problem,"_",fixed=T)[[1]][1] %in% c("alpine01") )
    return(0)
  
}