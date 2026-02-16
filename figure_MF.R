rm(list=ls()); graphics.off(); cat("\014")

source("multi_fidelity_test_problems.R")
library(plot3D)


# Forrester
xs <- seq(0,1,by=0.001)
ys1 <- forrester( xs, fidelity=1 )
ys2 <- forrester( xs, fidelity=2 )
ys3 <- forrester( xs, fidelity=3 )
ys4 <- forrester( xs, fidelity=4 )

plot( xs, ys1, type="l", lwd=4, ylim=c(-10,16), main="Forrester", xlab="x", ylab="", cex.main=2, cex.lab=2, cex.axis=2 )
lines( xs, ys2, lwd=4, lty=2, col="green3" )
lines( xs, ys3, lwd=4, lty=4, col="red" )
lines( xs, ys4, lwd=4, lty=5, col="orange" )
legend( "topleft", legend=c(expression(f[1](x)),expression(f[2](x)),expression(f[3](x)),expression(f[4](x))),
        cex=2, lwd=4, lty=c(1,2,4,5), col=c("black","green3","red","orange") )



# Rosenbrock
gg <- 20
XX <- as.matrix( expand.grid( x1=seq(0,1,length.out=gg), x2=seq(0,1,length.out=gg) ) )
ys1 <- apply( XX, 1, rosenbrock, fidelity=1 )
ys2 <- apply( XX, 1, rosenbrock, fidelity=2 )
ys3 <- apply( XX, 1, rosenbrock, fidelity=3 )

par(mfrow=c(1,3))
par(mar=c(1,1,1,1))
persp3D( sort(unique(XX[,1])), sort(unique(XX[,2])), matrix(ys1,gg,gg), zlim=range(ys1),
         bty="g", col=adjustcolor("white", alpha.f=0), border="black", colkey=F,
         cex.lab=2, xlab="x1", ylab="x2", zlab="f(x)" )
image3D( z=0, x=sort(unique(XX[,1])), y=sort(unique(XX[,2])), colvar=matrix(ys1,gg,gg),
         clim=range(ys1), resfac=3, col=heat.colors(100,rev=T), add=T, colkey=F )

persp3D( sort(unique(XX[,1])), sort(unique(XX[,2])), matrix(ys2,gg,gg), zlim=range(ys1),
         lwd=2, bty="g", col=adjustcolor("white", alpha.f=0), border="green3", colkey=F,
         cex.lab=2, xlab="x1", ylab="x2", zlab="f(x)" )
image3D( z=0, x=sort(unique(XX[,1])), y=sort(unique(XX[,2])), colvar=matrix(ys2,gg,gg),
         clim=range(ys1), resfac=3, col=heat.colors(100,rev=T), add=T, colkey=F )

persp3D( sort(unique(XX[,1])), sort(unique(XX[,2])), matrix(ys3,gg,gg), zlim=range(ys1),
         lwd=2, bty="g", col=adjustcolor("white", alpha.f=0), border="red", colkey=F,
         cex.lab=2, xlab="x1", ylab="x2", zlab="f(x)" )
image3D( z=0, x=sort(unique(XX[,1])), y=sort(unique(XX[,2])), colvar=matrix(ys3,gg,gg),
         clim=range(ys1), resfac=3, col=heat.colors(100,rev=T), add=T, colkey=F )


# Shifted-rotated Rastrigin
gg <- 20
XX <- as.matrix( expand.grid( x1=seq(0,1,length.out=gg), x2=seq(0,1,length.out=gg) ) )
ys1 <- apply( XX, 1, shiftedRotatedRastrigin, fidelity=1 )
ys2 <- apply( XX, 1, shiftedRotatedRastrigin, fidelity=2 )
ys3 <- apply( XX, 1, shiftedRotatedRastrigin, fidelity=3 )

par(mfrow=c(1,3))
par(mar=c(1,1,1,1))
persp3D( sort(unique(XX[,1])), sort(unique(XX[,2])), matrix(ys1,gg,gg), zlim=range(0,ys1,ys2,ys3),
         bty="g", col=adjustcolor("white", alpha.f=0), border="black", colkey=F,
         cex.lab=2, xlab="x1", ylab="x2", zlab="f(x)" )
image3D( z=min(0,ys1,ys2,ys3), x=sort(unique(XX[,1])), y=sort(unique(XX[,2])), colvar=matrix(ys1,gg,gg),
         clim=range(ys1,ys2,ys3), resfac=3, col=heat.colors(100,rev=T), add=T, colkey=F )

persp3D( sort(unique(XX[,1])), sort(unique(XX[,2])), matrix(ys2,gg,gg), zlim=range(ys1,ys2,ys3),
         lwd=2, bty="g", col=adjustcolor("white", alpha.f=0), border="green3", colkey=F,
         cex.lab=2, xlab="x1", ylab="x2", zlab="f(x)" )
image3D( z=min(ys1,ys2,ys3), x=sort(unique(XX[,1])), y=sort(unique(XX[,2])), colvar=matrix(ys2,gg,gg),
         clim=range(ys1,ys2,ys3), resfac=3, col=heat.colors(100,rev=T), add=T, colkey=F )

persp3D( sort(unique(XX[,1])), sort(unique(XX[,2])), matrix(ys3,gg,gg), zlim=range(ys1,ys2,ys3),
         lwd=2, bty="g", col=adjustcolor("white", alpha.f=0), border="firebrick", colkey=F,
         cex.lab=2, xlab="x1", ylab="x2", zlab="f(x)" )
image3D( z=min(ys1,ys2,ys3), x=sort(unique(XX[,1])), y=sort(unique(XX[,2])), colvar=matrix(ys3,gg,gg),
         clim=range(ys3), resfac=3, col=heat.colors(100,rev=T), add=T, colkey=F )



# Heterogeneous 1D
xs <- seq(0,1,by=0.001)
ys1 <- sapply( xs, heterogeneous, fidelity=1 )
ys2 <- sapply( xs, heterogeneous, fidelity=2 )

par(mfrow=c(1,1))
par(mar=c(5.1,4.1,4.1,2.1))

plot( xs, ys1, type="l", lwd=4, ylim=range(ys1,ys2), main="heterogeneous-1D", xlab="x", ylab="", cex.main=2, cex.lab=2, cex.axis=2 )
lines( xs, ys2, lwd=4, lty=2, col="green3" )
legend( "bottomright", legend=c(expression(f[1](x)),expression(f[2](x))),
        cex=2, lwd=4, lty=1:2, col=c("black","green3") )

# Heterogeneous 2D
gg <- 20
XX <- as.matrix( expand.grid( x1=seq(0,1,length.out=gg), x2=seq(0,1,length.out=gg) ) )
ys1 <- apply( XX, 1, heterogeneous, fidelity=1 )
ys2 <- apply( XX, 1, heterogeneous, fidelity=2 )

par(mfrow=c(1,2))
par(mar=c(1,1,1,1))
persp3D( sort(unique(XX[,1])), sort(unique(XX[,2])), matrix(ys1,gg,gg), zlim=range(ys1),
         bty="g", col=adjustcolor("white", alpha.f=0), border="black", colkey=F,
         cex.lab=2, xlab="x1", ylab="x2", zlab="f(x)" )
image3D( z=-0.625, x=sort(unique(XX[,1])), y=sort(unique(XX[,2])), colvar=matrix(ys1,gg,gg),
         clim=range(ys1), resfac=3, col=heat.colors(100,rev=T), add=T, colkey=F )

persp3D( sort(unique(XX[,1])), sort(unique(XX[,2])), matrix(ys2,gg,gg), zlim=range(ys1),
         lwd=2, bty="g", col=adjustcolor("white", alpha.f=0), border="green3", colkey=F,
         cex.lab=2, xlab="x1", ylab="x2", zlab="f(x)" )
image3D( z=-0.625, x=sort(unique(XX[,1])), y=sort(unique(XX[,2])), colvar=matrix(ys2,gg,gg),
         clim=range(ys1), resfac=3, col=heat.colors(100,rev=T), add=T, colkey=F )


# Paciorek 1D
xs <- seq(0,1,by=0.001)
ys1 <- sapply( xs, paciorek, fidelity=1, A=0.5, sig1=0.0125, sig2=0.075 )
ys2 <- sapply( xs, paciorek, fidelity=2, A=0.5, sig1=0.0125, sig2=0.075 )

par(mfrow=c(1,1))
par(mar=c(5.1,4.1,4.1,2.1))

plot( xs*0.7+0.3, ys1, type="l", lwd=3, ylim=range(ys1,ys2), main="Paciorek-1D", xlab="x", ylab="", cex.main=2, cex.lab=2, cex.axis=2 )
lines( xs*0.7+0.3, ys2, lwd=3, col="green3" )
legend( "bottomleft", legend=c(expression(f[1](x)),expression(f[2](x))),
        cex=2, lwd=4, col=c("black","green3") )


# Paciorek 2D
gg <- 30
XX <- as.matrix( expand.grid( x1=seq(0,1,length.out=gg), x2=seq(0,1,length.out=gg) ) )
ys1 <- apply( XX, 1, paciorek, fidelity=1, A=0.5, sig1=0.0125, sig2=0.075 )
ys2 <- apply( XX, 1, paciorek, fidelity=2, A=0.5, sig1=0.0125, sig2=0.075 )

par(mfrow=c(1,2))
par(mar=c(1,1,1,1))
persp3D( sort(unique(XX[,1])), sort(unique(XX[,2])), matrix(ys1,gg,gg), zlim=range(ys1,ys2),
         bty="g", col=adjustcolor("white", alpha.f=0), border="black", colkey=F,
         cex.lab=2, xlab="x1", ylab="x2", zlab="f(x)" )
image3D( z=min(ys1,ys2), x=sort(unique(XX[,1])), y=sort(unique(XX[,2])), colvar=matrix(ys1,gg,gg),
         clim=range(ys1), resfac=3, col=heat.colors(100,rev=T), add=T, colkey=F )

persp3D( sort(unique(XX[,1])), sort(unique(XX[,2])), matrix(ys2,gg,gg), zlim=range(ys1,ys2),
         lwd=2, bty="g", col=adjustcolor("white", alpha.f=0), border="green3", colkey=F,
         cex.lab=2, xlab="x1", ylab="x2", zlab="f(x)" )
image3D( z=min(ys1,ys2), x=sort(unique(XX[,1])), y=sort(unique(XX[,2])), colvar=matrix(ys2,gg,gg),
         clim=range(ys2), resfac=3, col=heat.colors(100,rev=T), add=T, colkey=F )
