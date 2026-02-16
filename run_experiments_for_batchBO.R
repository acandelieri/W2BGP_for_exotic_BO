rm(list=ls()); graphics.off(); cat("\014")

source("core.R")
source("test_problems.R")


#*****************************************************************************************************
# Experiment's setup
#*****************************************************************************************************

# [ 1D test problems ] *******************************************************************************
# d <- 1
# test_problem <- "problem_22" # 02, 03, 05, (06), 07, 11, 14, (15), 22


# [ 2D test problems ] *******************************************************************************
# d <- 2
# test_problem <- "bird" # bird, michalewicz, (rosenbrockModified), styblinskiTang, ursem03, ursemWaves


# [ D>2 test problems ] ******************************************************************************

# ***** fixed D *****
# test_problem <- "hartmann6" # hartmann3, hartmann6, sixHumpCamel
# d <- getDim(test_problem)

# ***** any D *****
d <- 20; test_problem <- "alpine01" # alpine01, styblinskiTang

# ****************************************************************************************************


# kernels of the agents (i.e., GPs) in the pool
agents_kernels <- c("exp","matern3_2","matern5_2","gauss")

# seeds for independent runs
start.seed <- 1; nSeeds <- 30


# check!
stopifnot( ( is.na(getDim(test_problem)) && d>1 ) || (d == getDim(test_problem)) )

# number of initial random queries
n0 <- max(d+1, min(2*d,10) )

# number of queries overall (including the initial random ones)
N <- min(30*d,150)


nStarts <- 1 # re-starts for L-BFGS-B in optimizing the acquisition function

#*****************************************************************************************************







#*****************************************************************************************************
# MAIN
#*****************************************************************************************************

testFun <- getFunction(test_problem)
M <- length(agents_kernels)

cat("|| * * * * * * * * * * [",test_problem,"d =",d,"] * * * * * * * * * * ||\n")

# 1 for 'self-confident' W2BGPBO, 2 for 'equally weighted' W2BGPBO, 3 for vanilla GPBO
RES1 <- RES2 <- RES3 <- NULL
times_1 <- times_2 <- times_3 <- NULL
  
for( seed in start.seed+(0:(nSeeds-1)) ) {
  
  # ******************************************************************************
  # Initialization (i.e., the same for all the methods)
  # ******************************************************************************
  
  cat("\n[ Experiment with seed =",seed,"]\n")
  set.seed(seed)
  
  cat("> Collecting initial random observations...")
  w2bgpbo <- w2bgpbo_ew <- gpbo <- NULL
  X <- as.matrix( maximinLHS(n0,d) )
  y <- apply(X,1,testFun)
  w2bgpbo <- data.frame( seed=rep(seed,n0),
                         iter=rep(0,n0),
                         agent=rep(NA,n0),
                         kernel=rep(NA,n0),
                         x=X,
                         y=y,
                         nugget=rep(NA,n0),
                         acquisition=rep("init",n0),
                         stringsAsFactors=F )
  
  w2bgpbo <- rbind( w2bgpbo, data.frame( seed=rep(NA,(N-n0)*M),
                                         iter=rep(NA,(N-n0)*M),
                                         agent=rep(NA,(N-n0)*M),
                                         kernel=rep(NA,(N-n0)*M),
                                         x=matrix(NA,(N-n0)*M,d),
                                         y=rep(NA,(N-n0)*M),
                                         nugget=rep(NA,(N-n0)*M),
                                         acquisition=rep(NA,(N-n0)*M),
                                         stringsAsFactors=F ) )

  w2bgpbo_ew <- gpbo <- w2bgpbo # are all the same!
  cat("Done!\n")
  
  
  
  
  # ***************************************************************************************
  # self-confident 
  # ***************************************************************************************
  
  agent_weight <- 0.5
  set.seed(seed)
  
  cat("> The 'self-confident' W2BGP-batchBO started:\n  [")
  t0 <- Sys.time()
  queryCount <- n0
  while( queryCount<N ) {
    cat("=")
    res <- next_batch_query( observations=w2bgpbo, agents_kernels=agents_kernels, agent_weight=agent_weight, acq="lcb", acq_par=1, nStarts=nStarts ) 
    y <- as.numeric(apply( res$X, 1, testFun ))
    # w2bgpbo <- rbind( w2bgpbo, data.frame( seed=seed,
    #                                        iter=rep(queryCount+1,length(y)),
    #                                        x=res$X,
    #                                        y=y,
    #                                        nugget=res$nuggets,
    #                                        acquisition="lcb",
    #                                        stringsAsFactors=F ) )
    ixs <- n0 + ( (M*(queryCount-n0)+1):(M*(queryCount-n0)+M) )
    w2bgpbo[ixs,] <- data.frame( seed=rep(seed,M),
                                 iter=rep(queryCount-n0+1,M),
                                 agent=1:M,
                                 kernel=agents_kernels,
                                 x=res$X,
                                 y=y,
                                 nugget=res$nuggets,
                                 acquisition=rep("lcb",M),
                                 stringsAsFactors=F )
    
    queryCount <- queryCount+1
  }
  cat("]\n")
  times_1 <- rbind( times_1, data.frame( seed=seed, 
                                         time=as.numeric( difftime( Sys.time(), t0, units="secs" ) ) ) )

  
  
  
  # ***************************************************************************************
  # Equally weighted
  # ***************************************************************************************
  
  agent_weight <- 1/M
  set.seed(seed)
  
  cat("> The 'equally' W2BGP-batchBO started:\n  [")
  t0 <- Sys.time()
  queryCount <- n0
  while( queryCount<N ) {
    cat("=")
    res <- next_batch_query( observations=w2bgpbo_ew, agents_kernels=agents_kernels, agent_weight=agent_weight, acq="lcb", acq_par=1, nStarts=nStarts ) 
    y <- as.numeric(apply( res$X, 1, testFun ))
    # w2bgpbo_ew <- rbind( w2bgpbo_ew, data.frame( seed=seed,
    #                                              iter=rep(queryCount+1,length(y)),
    #                                              x=res$X,
    #                                              y=y,
    #                                              nugget=res$nuggets,
    #                                              acquisition="lcb",
    #                                              stringsAsFactors=F ) )
    ixs <- n0 + ( (M*(queryCount-n0)+1):(M*(queryCount-n0)+M) )
    w2bgpbo_ew[ixs,] <- data.frame( seed=rep(seed,M),
                                    iter=rep(queryCount-n0+1,M),
                                    agent=1:M,
                                    kernel=agents_kernels,
                                    x=res$X,
                                    y=y,
                                    nugget=res$nuggets,
                                    acquisition=rep("lcb",M),
                                    stringsAsFactors=F )
    
    queryCount <- queryCount+1 
  }
  cat("]\n")
  times_2 <- rbind( times_2, data.frame( seed=seed, 
                                         time=as.numeric( difftime( Sys.time(), t0, units="secs" ) ) ) )
  
  
  
  
  # ***************************************************************************************
  # uncooperative
  # ***************************************************************************************
  
  agent_weight <- 1
  set.seed(seed)
  
  cat("> The 'uncooperative' GP-batchBO started:\n  [")
  t0 <- Sys.time()
  queryCount <- n0
  while( queryCount<N ) {
    cat("=")
    res <- next_batch_query( observations=gpbo, agents_kernels=agents_kernels, agent_weight=agent_weight, acq="lcb", acq_par=1, nStarts=nStarts ) 
    y <- as.numeric(apply( res$X, 1, testFun ))
    # gpbo <- rbind( gpbo, data.frame( seed=seed,
    #                                  iter=rep(queryCount+1,length(y)),
    #                                  x=res$X,
    #                                  y=y,
    #                                  nugget=res$nuggets,
    #                                  acquisition=rep("lcb",M),
    #                                  stringsAsFactors=F ) )
    ixs <- n0 + ( (M*(queryCount-n0)+1):(M*(queryCount-n0)+M) )
    gpbo[ixs,] <- data.frame( seed=rep(seed,M),
                              iter=rep(queryCount-n0+1,M),
                              agent=1:M,
                              kernel=agents_kernels,
                              x=res$X,
                              y=y,
                              nugget=res$nuggets,
                              acquisition=rep("lcb",M),
                              stringsAsFactors=F )    
    
    queryCount <- queryCount+1
  }
  cat("]\n")
  times_3 <- rbind( times_3, data.frame( seed=seed, 
                                         time=as.numeric( difftime( Sys.time(), t0, units="secs" ) ) ) )
  
  RES1 <- rbind( RES1, w2bgpbo )
  RES2 <- rbind( RES2, w2bgpbo_ew )
  RES3 <- rbind( RES3, gpbo )
}

cat("> Saving results...")
if( is.na(getDim(test_problem)) )
  test_problem <- paste0(test_problem,"_d",d)
today <- toString(Sys.Date())
if( !dir.exists( paste0("batchBO_results_",today,"/",test_problem) ) )
  dir.create( paste0("batchBO_results_",today,"/",test_problem), recursive=T )
saveRDS( RES1, paste0("batchBO_results_",today,"/",test_problem,"/W2BGPBO.RDS") )
saveRDS( RES2, paste0("batchBO_results_",today,"/",test_problem,"/W2BGPBO_ew.RDS") )
saveRDS( RES3, paste0("batchBO_results_",today,"/",test_problem,"/GPBO.RDS") )
saveRDS( times_1, paste0("batchBO_results_",today,"/",test_problem,"/W2BGPBO_times.RDS") )
saveRDS( times_2, paste0("batchBO_results_",today,"/",test_problem,"/W2BGPBO_ew_times.RDS") )
saveRDS( times_3, paste0("batchBO_results_",today,"/",test_problem,"/GPBO_times.RDS") )
cat(" Done!\n")
