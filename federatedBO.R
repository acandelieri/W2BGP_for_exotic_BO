rm(list=ls()); graphics.off(); cat("\014")

d <- 2

source("core.R")
source(paste0("test_functions_",d,"d.R"))


n0 <- 2*d
beta <- 1
N <- 30*d
M <- 4
w0 <- 0.5
covtypes <- c("exp","matern3_2","matern5_2","gauss")

if( M!=4 )
  covtypes <- sample( covtypes, M, replace=(M>4) )

# stopifnot(d==1); testFun <- "problem_22" # 02, 03, 05, 06, 07, 11, 14, 15, 22
stopifnot(d==2); testFun <- "ursemWaves"  # alpine01, bird, damavandi, levy03, michalewicz, qing,
                                                    # rosenbrockModified, styblinskiTang,ursem03, ursemWaves
seeds <- 1:30


RES <- NULL

for( seed in seeds ) {
  
  cat("\014")
  cat("********** seed =",seed,"/",length(seeds),"**********\n\n")

  cat("> Initializing agents: [")
  set.seed(seed)
  Xs <- ys <- list()
  for( m in 1:M ) {
    cat("=")
    Xs[[m]] <- as.matrix(maximinLHS( n0, d ))
    ys[[m]] <- apply(Xs[[m]],1,getFunction(testFun))
  }
  cat("]\n")
  
  
  cat("> Starting Federated BO:\n")
  
  while( length(ys[[1]])<N ) {
    
    cat("  - Fitting GPs: [")
    gps <- list()
    for( m in 1:M ) {
      cat("=")
      gp <- NULL
      try( expr=( gp <- km( design=data.frame(x=Xs[[m]]), response=ys[[m]], covtype=covtypes[[m]], nugget.estim=F, control=list(trace=0) )),
           silent=T )
      if( is.null(gp) )
        gp <- km( design=data.frame(x=Xs[[m]]), response=ys[[m]], covtype=covtypes[[m]], nugget.estim=T, control=list(trace=0) )
      gps[[m]] <- gp
    }
    cat("]\n")
    
    cat("  - weighted Wasserstein Barycenter GP-LCB and evaluation: [")
    par0 <- runif(d)
    for( m in 1:M ) {
      cat("=")
      ls <- rep((1-w0)/(M-1),M)
      ls[m] <- w0
      res <- optim( par=par0, fn=wlcb, gr=NULL, method="L-BFGS-B", lower=0, upper=1, control=list(trace=0),
                    gps=gps, ls=ls, beta=beta )
      x_ <- res$par
      y_ <- do.call(getFunction(testFun),list(x=x_))
      
      Xs[[m]] <- rbind( Xs[[m]], t(x_) )
      ys[[m]] <- c( ys[[m]], y_ )
    }
    cat("]\n")
    cat("> Performed queries:",length(ys[[1]]),"out of",N,"\n")
    
    
    #****************************************************************************
    # par(mfrow=c(M,2))
    # par(mar=c(4.1,4.1,1.1,1.1))
    # xx <- seq(0,1,by=0.001)
    # yy <- do.call(getFunction(testFun),list(x=xx))
    # clrs <- rainbow(M)
    # for( m in 1:M ) {
    #   plot( xx, yy, type="l", col="grey", lwd=2, ylim=2*range(yy) )
    #   pred <- predict( gps[[m]], data.frame(x=xx), "UK" )
    #   polygon( c(xx,rev(xx)), c(pred$mean+pred$sd,rev(pred$mean-pred$sd)),
    #            col=adjustcolor(clrs[m],alpha.f=0.2), border=F )
    #   lines( xx, pred$mean, col=clrs[m], lwd=2 )
    #   points( Xs[[m]], ys[[m]], pch=19, cex=2, col=clrs[m] )
    #   legend( "bottomleft", legend=covtypes[m], lwd=2, col=clrs[m], cex=1.5 )
    # 
    #   plot( xx, yy, type="l", col="grey", lwd=2, ylim=2*range(yy) )
    #   ls <- rep(0.5/(M-1),M); ls[m] <- 0.5
    #   pred <- wgp.predict( x=xx, gps=gps, ls=ls )
    #   polygon( c(xx,rev(xx)), c(pred$mean+pred$sd,rev(pred$mean-pred$sd)),
    #            col=adjustcolor(clrs[m],alpha.f=0.2), border=F )
    #   lines( xx, pred$mean, col=clrs[m], lwd=2 )
    #   points( Xs[[m]], ys[[m]], pch=19, cex=2, col=clrs[m] )
    #   points( Xs[[m]][length(ys[[m]])], ys[[m]][length(ys[[m]])], pch=25, cex=2, bg="orange" )
    # 
    # }
    # par(mfrow=c(1,1))
    # invisible( readline( "[RETURN]" ) )
    #****************************************************************************
  }
  
  X <- y <- NULL
  ags <- covs <- NULL
  for( m in 1:M ) {
    X <- rbind( X, Xs[[m]] )
    y <- c( y, ys[[m]] )
    ags <- c( ags, rep(m,N) )
    covs <- c( covs, rep(covtypes[m],N) ) 
  }
  
  RES <- rbind( RES, data.frame( seed=rep(seed,N*M),
                                 agent=ags,
                                 covtype=covs,
                                 acquisition=rep( c( rep("init",n0), rep("lcb",(N-n0)) ), M ),
                                 X=X, y=y ) )
  
}

folder <- paste0("federated-BO/test_function",d,"d/w2bgpbo/")
if( !dir.exists(folder) )
  dir.create( folder, recursive=T )
saveRDS( RES, paste0(folder,"/",testFun,".RDS") )
