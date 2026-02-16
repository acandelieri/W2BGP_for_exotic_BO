rm(list=ls()); graphics.off(); cat("\014")


probs <- c(0.25,0.5,0.75)

test_problem_name <- "forrester"; d <- 1; ystar <- -6.020740; id <- 1
# test_problem_name <- "rosenbrock"; d <- 2; ystar <- 0; id <- 2
# test_problem_name <- "rosenbrock"; d <- 5; ystar <- 0; id <- 3
# test_problem_name <- "shiftedRotatedRastrigin"; d <- 2; ystar <- 0; id <- 4
# test_problem_name <- "heterogeneous"; d <- 1; ystar <-  -0.625; id <- 5
# test_problem_name <- "heterogeneous"; d <- 2; ystar <- -0.5627123; id <- 6
# test_problem_name <- "heterogeneous"; d <- 3; ystar <- -0.5627123; id <- 7
# test_problem_name <- "pacioreck"; d <- 2; ystar <- -1; id <- 8

RES1 <- readRDS( paste0("multifidelityBO_results/",test_problem_name,"_d",d,"/W2BGPBO_fw.RDS") )
RES2 <- readRDS( paste0("multifidelityBO_results/",test_problem_name,"_d",d,"/W2BGPBO_pw.RDS") )
RES3 <- readRDS( paste0("multifidelityBO_results/",test_problem_name,"_d",d,"/W2BGPBO_ew.RDS") )

seeds <- sort(unique(RES1$seed))
iters <- sort(unique(RES1$iter))
VALS1 <- VALS2 <- VALS3 <- matrix(Inf,length(seeds),max(RES1$iter,RES2$iter,RES3$iter)+1)

for( seed in seeds ) {
  
  tmp <- RES1[RES1$seed==seed & RES1$source==1, ]
  aux <- c( min(tmp$y[tmp$iter==0]), tmp$y[tmp$iter!=0] )
  VALS1[ seed, 1+sort(unique(tmp$iter)) ] <- aux
  
  tmp <- RES2[RES2$seed==seed & RES2$source==1, ]
  aux <- c( min(tmp$y[tmp$iter==0]), tmp$y[tmp$iter!=0] )
  VALS2[ seed, 1+sort(unique(tmp$iter)) ] <- aux
  
  tmp <- RES3[RES3$seed==seed & RES3$source==1, ]
  aux <- c( min(tmp$y[tmp$iter==0]), tmp$y[tmp$iter!=0] )
  VALS3[ seed, 1+sort(unique(tmp$iter)) ] <- aux
}

BS1 <- t( apply(VALS1,1,cummin) ); BS2 <- t( apply(VALS2,1,cummin) ); BS3 <- t( apply(VALS3,1,cummin) )

fBS1 <- BS1[,ncol(BS1)]; fBS2 <- BS2[,ncol(BS2)]; fBS3 <- BS3[,ncol(BS3)]
test13 <- wilcox.test( fBS1, fBS3, paired=T, exact=T )
test23 <- wilcox.test( fBS2, fBS3, paired=T, exact=T )

cat("***** [",test_problem_name,"] *****\n\n")
cat("final best seen:\n")
cat(" fidelities\t\trescaled\t\tequal\t\t\tp-value\t\tp-value\n")
cat(" ",round(median(fBS1),4)," (",round(sd(fBS1),4),")\t",
    round(median(fBS2),4)," (",round(sd(fBS2),4),")\t",
    round(median(fBS3),4)," (",round(sd(fBS3),4),")\t",
    round(test13$p.value,3),"\t\t",round(test23$p.value),"\n\n", sep="" )



qs1 <- apply(BS1,2,quantile,probs=probs); qs2 <- apply(BS2,2,quantile,probs=probs); qs3 <- apply(BS3,2,quantile,probs=probs)

xAxis <- 1:ncol(BS1) - 1

# plot( xAxis, qs1[2,], type="l", ylim=range(ystar,qs1,qs2,qs3) )
# polygon( c(xAxis,rev(xAxis)), c(qs1[1,],rev(qs1[3,])), col=adjustcolor("blue",alpha.f=0.1), border=F )
# polygon( c(xAxis,rev(xAxis)), c(qs2[1,],rev(qs2[3,])), col=adjustcolor("green",alpha.f=0.1), border=F )
# polygon( c(xAxis,rev(xAxis)), c(qs3[1,],rev(qs3[3,])), col=adjustcolor("red",alpha.f=0.1), border=F )
# lines( xAxis, qs1[2,], lwd=3, col=adjustcolor("blue",alpha.f=0.8) )
# lines( xAxis, qs2[2,], lwd=3, col=adjustcolor("green4",alpha.f=0.8) )
# lines( xAxis, qs3[2,], lwd=3, col=adjustcolor("red",alpha.f=0.8) )






# *************************************************************************************************************
# GAP metric and AUGC
# *************************************************************************************************************

G1 <- G2 <- G3 <- matrix(0,nrow(BS1),ncol(BS1))
for( i in 1:nrow(G1) ) {
  for( j in 2:ncol(G2) ) {
    G1[i,j] <- (BS1[i,1] - BS1[i,j])/(BS1[i,1] - ystar)
    G2[i,j] <- (BS2[i,1] - BS2[i,j])/(BS2[i,1] - ystar)
    G3[i,j] <- (BS3[i,1] - BS3[i,j])/(BS3[i,1] - ystar)
  }
}


# qs1 <- apply(G1,2,quantile,probs=probs); qs2 <- apply(G2,2,quantile,probs=probs); qs3 <- apply(G3,2,quantile,probs=probs)
qs1 <- t(apply(G1,2,mean)); qs1 <- rbind(qs1+apply(G1,2,sd),qs1); qs1 <- rbind(qs1,2*qs1[2,]-qs1[1,])
qs2 <- t(apply(G2,2,mean)); qs2 <- rbind(qs2+apply(G2,2,sd),qs2); qs2 <- rbind(qs2,2*qs2[2,]-qs2[1,])
qs3 <- t(apply(G3,2,mean)); qs3 <- rbind(qs3+apply(G3,2,sd),qs3); qs3 <- rbind(qs3,2*qs3[2,]-qs3[1,])
qs1[qs1<0] <- 0; qs2[qs2<0] <- 0; qs3[qs3<0] <- 0
qs1[qs1>1] <- 1; qs2[qs2>1] <- 1; qs3[qs3>1] <- 1

par(mar=c(5.1,5.1,4.1,1.1))

if( !dir.exists("gap_figures/mfBO") )
  dir.create("gap_figures/mfBO",recursive=T)
png( paste0("gap_figures/mfBO/GapMetric_mfBO_0",id,".png"), width=500, height=500 )
plot( xAxis, qs1[2,], type="l", ylim=0:1, main=paste0(test_problem_name,"_d",d),
      ylab="Gap metric", xlab="queries", cex.axis=2, cex.lab=2, cex.main=2 )
polygon( c(xAxis,rev(xAxis)), c(qs1[1,],rev(qs1[3,])), col=adjustcolor("blue",alpha.f=0.15), border=F )
polygon( c(xAxis,rev(xAxis)), c(qs2[1,],rev(qs2[3,])), col=adjustcolor("green",alpha.f=0.15), border=F )
polygon( c(xAxis,rev(xAxis)), c(qs3[1,],rev(qs3[3,])), col=adjustcolor("red",alpha.f=0.15), border=F )
lines( xAxis, qs1[2,], lwd=4, col=adjustcolor("blue",alpha.f=0.8) )
lines( xAxis, qs2[2,], lwd=4, col=adjustcolor("green4",alpha.f=0.8) )
lines( xAxis, qs3[2,], lwd=4, col=adjustcolor("red",alpha.f=0.8) )
abline( h=0 ); abline( h=1, lty=2, lwd=2 )

if( test_problem_name %in% c("forrester","shiftedRotatedRastrigin") ) {
  legend( "topleft", legend=c("fidelities as weights ","rescaled weights","equal weights"),
          col=c("green4","blue","red3"), lwd=4, cex=1.5  )  
} else {
  legend( "bottomright", legend=c("fidelities as weights ","rescaled weights","equal weights"),
          col=c("green4","blue","red3"), lwd=4, cex=1.5  )
}
dev.off()

AUGC1 <- apply(G1,1,sum); AUGC2 <- apply(G2,1,sum); AUGC3 <- apply(G3,1,sum)
test13 <- wilcox.test( AUGC1, AUGC3, paired=T, exact=T )
test23 <- wilcox.test( AUGC2, AUGC3, paired=T, exact=T )



cat("AUGC:\n")
cat(" fidelities\t\trescaled\t\tequal\t\t\tp-value\t\tp-value\n")
cat(" ",round(median(AUGC1),4)," (",round(sd(AUGC1),4),")\t",
    round(median(AUGC2),4)," (",round(sd(AUGC2),4),")\t",
    round(median(AUGC3),4)," (",round(sd(AUGC3),4),")\t",
    round(test13$p.value,3),"\t\t",round(test23$p.value),"\n\n", sep="" )


