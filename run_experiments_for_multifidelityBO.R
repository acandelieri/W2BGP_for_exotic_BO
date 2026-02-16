rm(list=ls()); graphics.off(); cat("\014")

source("core.R")
source("multi_fidelity_test_problems.R")


#*****************************************************************************************************
# Experiment's setup
#*****************************************************************************************************

# [ 1D test problems ] *******************************************************************************
# d <- 1
# test_problem <- "forrester"; fidelityCosts <- c(1,0.5,0.1,0.05)
# test_problem <- "heterogeneous"; fidelityCosts <- c(1,0.2)


# [ D>=2 test problems ] *******************************************************************************
d <- 5 # 2, 5, or 10
# test_problem <- "rosenbrock"; fidelityCosts <- c(1,0.5,0.1)
test_problem <- "shiftedRotatedRastrigin"; fidelityCosts <- c(1, 0.0625, 0.00390625)

# d <- 3 # 1, 2 or 3
# test_problem <- "heterogeneous"; fidelityCosts <- c(1,0.2)

# d <- 2 # only d=2
# test_problem <- "pacioreck"; fidelityCosts <- c(1, 0.2)

# ****************************************************************************************************


# kernel type for all the GPs
kernel <- "matern3_2"

# seeds for independent runs
start.seed <- 1; nSeeds <- 30


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
M <- length(fidelityCosts)

cat("|| * * * * * * * * * * [",test_problem,"d =",d,"] * * * * * * * * * * ||\n")

# 1 for 'self-confident' W2BGPBO, 2 for 'equally weighted' W2BGPBO
RES1 <- RES2 <- RES3 <- NULL
times_1 <- times_2 <- times_3 <- NULL
  
for( seed in start.seed+(0:(nSeeds-1)) ) {
  
  # ******************************************************************************
  # Initialization (i.e., the same for all the methods)
  # ******************************************************************************
  
  cat("\n[ Experiment with seed =",seed,"]\n")
  set.seed(seed)
  
  cat("> Initializing one GP model for every information source...")
  w2bgpbo <- w2bgpbo_pw <- w2bgpbo_ew <- NULL
  for( ns in 1:length(fidelityCosts) ) {
    X <- as.matrix( maximinLHS(n0,d) )
    y <- apply(X,1,testFun,fidelity=ns) # using default values for all the multi-fidelity test problems
    w2bgpbo <- rbind( w2bgpbo, data.frame( seed=rep(seed,n0),
                                           iter=numeric(n0),
                                           source=rep(ns,n0),
                                           x=X,
                                           y=y,
                                           nugget=rep(NA,n0),
                                           acquisition=rep("init",n0),
                                           stringsAsFactors=F ) )
  }
  w2bgpbo_ew <- w2bgpbo_pw <- w2bgpbo # are all the same!
  cat("Done!\n")
  
  
  
  
  # ***************************************************************************************
  # Multi-fidelity with WB's weights equal to fidelity
  # ***************************************************************************************
  
  set.seed(seed)
  
  cat("> The 'MF' W2BGPBO-based algorithm (with WB's weights equal to fidelities) started:\n  [")
  t0 <- Sys.time()
  queryCount <- n0
  while( queryCount<N ) {
    cat("=")
    res <- next_multifidelity_query( observations=w2bgpbo, kernel=kernel, fidelityCosts=fidelityCosts, acq_par=1, nStarts=nStarts ) 
    y <- testFun( res$x, fidelity=res$source_ix )
    w2bgpbo <- rbind( w2bgpbo, data.frame( seed=seed,
                                           iter=queryCount+1,
                                           source=res$source_ix,
                                           x=t(res$x),
                                           y=y,
                                           nugget=toString(res$nuggets),
                                           acquisition="mf-w2bgplcb",
                                           stringsAsFactors=F ) )
    queryCount <- queryCount+1
  }
  cat("]\n")
  times_1 <- rbind( times_1, data.frame( seed=seed, 
                                         time=as.numeric( difftime( Sys.time(), t0, units="secs" ) ) ) )

  
  
  
  # ***************************************************************************************
  # Multi-fidelity with WB's weights proportional to fidelity of each information source
  # ***************************************************************************************
  
  set.seed(seed)
  
  cat("> The 'MF' W2BGPBO-based algorithm (with WB's weights proportional but different from fidelities) started:\n  [")
  t0 <- Sys.time()
  queryCount <- n0
  propFidelityCosts <- numeric()
  aux <- 1
  for( i in 1:(M-1) ) {
    propFidelityCosts[i] <- 0.75*aux
    aux <- 0.25*aux
  }
  propFidelityCosts[M] <- aux
  while( queryCount<N ) {
    cat("=")
    res <- next_multifidelity_query( observations=w2bgpbo_pw, kernel=kernel, fidelityCosts=propFidelityCosts, acq_par=1, nStarts=nStarts ) 
    y <- testFun( res$x, fidelity=res$source_ix )
    w2bgpbo_pw <- rbind( w2bgpbo_pw, data.frame( seed=seed,
                                                 iter=queryCount+1,
                                                 source=res$source_ix,
                                                 x=t(res$x),
                                                 y=y,
                                                 nugget=toString(res$nuggets),
                                                 acquisition="mf-w2bgplcb",
                                                 stringsAsFactors=F ) )
    queryCount <- queryCount+1
  }
  cat("]\n")
  times_2 <- rbind( times_2, data.frame( seed=seed, 
                                         time=as.numeric( difftime( Sys.time(), t0, units="secs" ) ) ) )
  
  
  
  
  
  # ***************************************************************************************
  # Multi-fidelity with uniform WB's weights (i.e., equal fidelities)
  # ***************************************************************************************

  equalCosts <- rep(1,M) 
  set.seed(seed)

  cat("> The 'equally' W2BGPBO-based MF optimization:\n  [")
  t0 <- Sys.time()
  queryCount <- n0
  while( queryCount<N ) {
    cat("=")
    res <- next_multifidelity_query( observations=w2bgpbo_ew, kernel=kernel, fidelityCosts=equalCosts, acq_par=1, nStarts=nStarts ) 
    y <- testFun( res$x, fidelity=res$source_ix )
    w2bgpbo_ew <- rbind( w2bgpbo_ew, data.frame( seed=seed,
                                                 iter=queryCount+1,
                                                 source=res$source_ix,
                                                 x=t(res$x),
                                                 y=y,
                                                 nugget=toString(res$nuggets),
                                                 acquisition="mf-w2bgplcb",
                                                 stringsAsFactors=F ) )
    queryCount <- queryCount+1
  }
  cat("]\n")
  times_3 <- rbind( times_3, data.frame( seed=seed,
                                         time=as.numeric( difftime( Sys.time(), t0, units="secs" ) ) ) )
  
  
  
  RES1 <- rbind( RES1, w2bgpbo )
  RES2 <- rbind( RES2, w2bgpbo_pw )
  RES3 <- rbind( RES3, w2bgpbo_ew )
}

cat("> Saving results...")
test_problem <- paste0(test_problem,"_d",d)
today <- toString(Sys.Date())
if( !dir.exists( paste0("multifidelityBO_results_",today,"/",test_problem) ) )
  dir.create( paste0("multifidelityBO_results_",today,"/",test_problem), recursive=T )
saveRDS( RES1, paste0("multifidelityBO_results_",today,"/",test_problem,"/W2BGPBO_fw.RDS") )
saveRDS( RES2, paste0("multifidelityBO_results_",today,"/",test_problem,"/W2BGPBO_pw.RDS") )
saveRDS( RES3, paste0("multifidelityBO_results_",today,"/",test_problem,"/W2BGPBO_ew.RDS") )
saveRDS( times_1, paste0("multifidelityBO_results_",today,"/",test_problem,"/W2BGPBO_fw_times.RDS") )
saveRDS( times_2, paste0("multifidelityBO_results_",today,"/",test_problem,"/W2BGPBO_pw_times.RDS") )
saveRDS( times_3, paste0("multifidelityBO_results_",today,"/",test_problem,"/W2BGPBO_ew_times.RDS") )
cat(" Done!\n")
