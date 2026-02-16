rm(list=ls()); graphics.off(); cat("\014")

source("test_problems.R")

visualizeQuantiles <- F


# task <- "federatedBO"
task <- "batchBO"

tps <- c( paste0( "problem_", c("02","03","05","07","11","14","15","22") ),
          "alpine01_d2",
          "bird","levy03_d2","michalewicz","styblinskiTang_d2","ursem03","ursemWaves",
          "hartmann3","hartmann6",
          "alpine01_d5","alpine01_d10","alpine01_d20",
          "levy03_d5","levy03_d10","levy03_d20",
          "styblinskiTang_d5","styblinskiTang_d10","styblinskiTang_d20")


STATS <- STATS_bs <- STATS_rg <- NULL

for( tp in tps ) {
  
  # read results from files
  w2bgpbo <- readRDS( paste0(task,"_results/",tp,"/W2BGPBO.RDS") )
  w2bgpbo_ew <- readRDS( paste0(task,"_results/",tp,"/W2BGPBO_ew.RDS") )
  gpbo <- readRDS( paste0(task,"_results/",tp,"/GPBO.RDS") )
  
  
  # W2BGPBO
  tmp <- aggregate( w2bgpbo$y, by=list(w2bgpbo$iter,w2bgpbo$seed), min )
  names(tmp) <- c("iter","seed","y")
  bs <- rg <- NULL
  for( s in sort(unique(tmp$seed)) ) {
    bs <- rbind( bs, cummin(tmp$y[tmp$seed==s]) )
    rg <- rbind( rg, tmp$y[tmp$seed==s]-getOptY(tp) )
  }
  N <- max(tmp$iter); n0 <- N-length(unique(tmp$iter))+1
  
  gap <- t(apply( bs, 1, function(x) { (x[1]-x)/(x[1]-getOptY(tp)) } ))
  
  
  # W2BGPBO_ew
  tmp <- aggregate( w2bgpbo_ew$y, by=list(w2bgpbo_ew$iter,w2bgpbo_ew$seed), min )
  names(tmp) <- c("iter","seed","y")
  bs_ew <- rg_ew <- NULL
  for( s in sort(unique(tmp$seed)) ) {
    bs_ew <- rbind( bs_ew, cummin(tmp$y[tmp$seed==s]) )
    rg_ew <- rbind( rg_ew, tmp$y[tmp$seed==s]-getOptY(tp) )
  }
  N <- max(tmp$iter); n0 <- N-length(unique(tmp$iter))+1
  
  gap_ew <- t(apply( bs_ew, 1, function(x) { (x[1]-x)/(x[1]-getOptY(tp)) } ))
  
  
  # GPBO
  tmp <- aggregate( gpbo$y, by=list(gpbo$iter,gpbo$seed), min )
  names(tmp) <- c("iter","seed","y")
  bs_bo <- rg_bo <- NULL
  for( s in sort(unique(tmp$seed)) ) {
    bs_bo <- rbind( bs_bo, cummin(tmp$y[tmp$seed==s]) )
    rg_bo <- rbind( rg_bo, tmp$y[tmp$seed==s]-getOptY(tp) )
  }
  N <- max(tmp$iter); n0 <- N-length(unique(tmp$iter))+1
  
  gap_bo <- t(apply( bs_bo, 1, function(x) { (x[1]-x)/(x[1]-getOptY(tp)) } ))
  
  
  
  a_vs_b <- wilcox.test( round(apply(gap,1,mean),4), round(apply(gap_ew,1,mean),4), paired=T, exact=T )
  a_vs_c <- wilcox.test( round(apply(gap,1,mean),4), round(apply(gap_bo,1,mean),4), paired=T, exact=T )
    
  STATS <- rbind( STATS, data.frame( test_problem=tp,
                                     W2BGPBO=paste0(round(median(apply(gap,1,mean)),4),"(",round(sd(apply(gap,1,mean)),4),")"),
                                     W2BGPBO_ew=paste0(round(median(apply(gap_ew,1,mean)),4),"(",round(sd(apply(gap_ew,1,mean)),4),")"),
                                     GPBO=paste0(round(median(apply(gap_bo,1,mean)),4),"(",round(sd(apply(gap_bo,1,mean)),4),")"),
                                     pvalue_A_vs_B=round(a_vs_b$p.value,3),
                                     pvalue_A_vs_C=round(a_vs_c$p.value,3),
                                     stringsAsFactors=T ) )
  
  
  a_vs_b <- wilcox.test( round(bs[,ncol(bs)],4), round(bs_ew[,ncol(bs_ew)],4), paired=T, exact=T )
  a_vs_c <- wilcox.test( round(bs[,ncol(bs)],4), round(bs_bo[,ncol(bs_bo)],4), paired=T, exact=T )
  STATS_bs <- rbind( STATS_bs, data.frame( test_problem=tp,
                                           W2BGPBO=paste0(round(median(bs[,ncol(bs)]),4),"(",round(sd(bs[,ncol(bs)]),4),")"),
                                           W2BGPBO_ew=paste0(round(median(bs_ew[,ncol(bs_ew)]),4),"(",round(sd(bs_ew[,ncol(bs_ew)]),4),")"),
                                           GPBO=paste0(round(median(bs_bo[,ncol(bs_bo)]),4),"(",round(sd(bs_bo[,ncol(bs_bo)]),4),")"),
                                           pvalue_A_vs_B=round(a_vs_b$p.value,3),
                                           pvalue_A_vs_C=round(a_vs_c$p.value,3),
                                           stringsAsFactors=T ) )
  
  
  
  

  # ****************************************************************************************************************
  # Gap metric chart
  # ****************************************************************************************************************
  
  if( !visualizeQuantiles ) {
    md <- apply( gap, 2, mean ); lo <- md - apply( gap, 2, sd ); up <- md + apply( gap, 2, sd )
    up[up>1] <- 1; lo[lo<0] <- 0
  
    md_ew <- apply( gap_ew, 2, mean ); lo_ew <- md_ew - apply( gap_ew, 2, sd ); up_ew <- md_ew + apply( gap_ew, 2, sd )
    up_ew[up_ew>1] <- 1; lo_ew[lo_ew<0] <- 0
  
    md_bo <- apply( gap_bo, 2, mean ); lo_bo <- md_bo - apply( gap_bo, 2, sd ); up_bo <- md_bo + apply( gap_bo, 2, sd )
    up_bo[up_bo>1] <- 1; lo_bo[lo_bo<0] <- 0
  } else {
    tmp <- apply( gap, 2, quantile, probs=c(0.5,0.25,0.75) ); md <- tmp[1,]; lo <- tmp[2,]; up <- tmp[3,]
    tmp <- apply( gap_ew, 2, quantile, probs=c(0.5,0.25,0.75) ); md_ew <- tmp[1,]; lo_ew <- tmp[2,]; up_ew <- tmp[3,]
    tmp <- apply( gap_bo, 2, quantile, probs=c(0.5,0.25,0.75) ); md_bo <- tmp[1,]; lo_bo <- tmp[2,]; up_bo <- tmp[3,]
  }
  
  
  if( !dir.exists(paste0("gap_figures/",task)) )
    dir.create(paste0("gap_figures/",task),recursive=T)
  png( paste0("gap_figures/",task,"/GapMetric_",task,"_",ifelse(which(tps==tp)<10,"0",""),which(tps==tp),".png" ),
       width=500, height=500 ) 
  
  par(mar=c(5.1,5.1,4.1,1.1))
  plot( n0:N, md, type="l", ylim=0:1, main=tp, ylab="Gap metric", xlab="queries",
        cex.axis=2, cex.lab=2, cex.main=2 )
  abline( h=0 ); abline( h=1, lty=2, lwd=2 )

  
  
  # ****************************************************************************************************************
  # cumulative regret chart
  # ****************************************************************************************************************
  # 
  # tmp <- t(apply( rg, 1, cumsum ))
  # md <- apply( tmp, 2, median ); lo <- md - apply( tmp, 2, sd ); up <- md + apply( tmp, 2, sd )
  # 
  # tmp <- t(apply( rg_ew, 1, cumsum ))
  # md_ew <- apply( tmp, 2, median ); lo_ew <- md_ew - apply( tmp, 2, sd ); up_ew <- md_ew + apply( tmp, 2, sd )
  # 
  # tmp <- t(apply( rg_bo, 1, cumsum ))
  # md_bo <- apply( rg_bo, 2, median ); lo_bo <- md_bo - apply( rg_bo, 2, sd ); up_bo <- md_bo + apply( rg_bo, 2, sd )
  # 
  # 
  # plot( n0:N, md, type="l", ylim=range(lo,lo_ew,lo_bo,up,up_ew,up_bo), main=tp )
  # abline( h=0, lty=2, lwd=2 )
  # 
  # *****************************************************************************************************************
  
  
  polygon( c(n0:N,N:n0), c(lo_bo,rev(up_bo)), col=adjustcolor("red2",alpha.f=0.15), border=F )
  polygon( c(n0:N,N:n0), c(lo,rev(up)), col=adjustcolor("green3",alpha.f=0.15), border=F )
  polygon( c(n0:N,N:n0), c(lo_ew,rev(up_ew)), col=adjustcolor("deepskyblue",alpha.f=0.15), border=F )
  
  lines( n0:N, md_bo, col="red3", lwd=4 ); lines( n0:N, md, col="green4", lwd=4 ); lines( n0:N, md_ew, col="blue", lwd=4 )
  
  legend( "bottomright", legend=c("self-confident","equal weights","uncooperative"),
          col=c("green4","blue","red3"), lwd=4, cex=2  )
  
  dev.off()
  
}

cat("[AUGC]\n")
print(STATS)

cat("\n[Best Seen]\n")
print(STATS_bs)






