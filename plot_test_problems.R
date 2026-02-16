rm(list=ls()); graphics.off(); cat("\014")

source("test_problems.R")
source("multi_fidelity_test_problems.R")


# **********************************************************************
# 1D test problems
# **********************************************************************

par(mar=c(2.6,5.1,2.6,1.1) )
par(mfrow=c(4,2))
xs <- seq(0,1,by=0.001)

ys <- problem_02(xs)
plot(xs*(7.5-2.7)+2.7,ys,type="l",lwd=3,xlab="x",ylab="f(x)",main="problem 02",
     cex.lab=2,cex.axis=2,cex.main=2)
ys <- problem_03(xs)
plot(xs*20-10,ys,type="l",lwd=3,xlab="x",ylab="f(x)",main="problem 03",
     cex.lab=2,cex.axis=2,cex.main=2)
ys <- problem_05(xs)
plot(xs*1.2,ys,type="l",lwd=3,xlab="x",ylab="f(x)",main="problem 05",
     cex.lab=2,cex.axis=2,cex.main=2)
ys <- problem_07(xs)
plot(xs*(7.5-2.7) + 2.7,ys,type="l",lwd=3,xlab="x",ylab="f(x)",main="problem 07",
     cex.lab=2,cex.axis=2,cex.main=2)
ys <- problem_11(xs)
plot(xs*( 2*pi + pi/2) -pi/2,ys,type="l",lwd=3,xlab="x",ylab="f(x)",main="problem 11",
     cex.lab=2,cex.axis=2,cex.main=2)
ys <- problem_14(xs)
plot(xs*4,ys,type="l",lwd=3,xlab="x",ylab="f(x)",main="problem 14",
     cex.lab=2,cex.axis=2,cex.main=2)
ys <- problem_15(xs)
plot(xs*10-5,ys,type="l",lwd=3,xlab="x",ylab="f(x)",main="problem 15",
     cex.lab=2,cex.axis=2,cex.main=2)
ys <- problem_22(xs)
plot(xs*20,ys,type="l",lwd=3,xlab="x",ylab="f(x)",main="problem 22",
     cex.lab=2,cex.axis=2,cex.main=2)




# **********************************************************************
# 2D test problems
# **********************************************************************

library(plot3D)

par(mar=rep(1.5,4))
par(mfrow=c(2,3))

XX <- expand.grid( x1=seq(0,1,by=0.01), x2=seq(0,1,by=0.01) )

ys <- apply(XX,1,alpine01)
persp3D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
         xlab="x1", ylab="x2", zlab="f(x)", main="alpine01",
         resfac=2, cex.lab=2, cex.main=2, colkey=F )
ys <- apply(XX,1,bird)
persp3D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
         xlab="x1", ylab="x2", zlab="f(x)", main="bird",
         resfac=2, cex.lab=2, cex.main=2, colkey=F )
# ys <- apply(XX,1,levy03)
# persp3D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
#          xlab="x1", ylab="x2", zlab="f(x)", main="levy03",
#          resfac=2, cex.lab=2, cex.main=2, colkey=F )
ys <- apply(XX,1,michalewicz)
persp3D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
         xlab="x1", ylab="x2", zlab="f(x)", main="michalewicz",
         resfac=2, cex.lab=2, cex.main=2, colkey=F )
ys <- apply(XX,1,styblinskiTang)
persp3D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
         xlab="x1", ylab="x2", zlab="f(x)", main="styblinskiTang",
         resfac=2, cex.lab=2, cex.main=2, colkey=F )
ys <- apply(XX,1,ursem03)
persp3D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
         xlab="x1", ylab="x2", zlab="f(x)", main="ursem03",
         resfac=2, cex.lab=2, cex.main=2, colkey=F )
# plot.new()
ys <- apply(XX,1,ursemWaves)
persp3D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
         resfac=2, xlab="x1", ylab="x2", zlab="f(x)", main="ursemWaves",
         cex.lab=2, cex.main=2, colkey=F )
# plot.new()


par(mar=rep(2,4))

ys <- apply(XX,1,alpine01)
image2D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
         xlab="x1", ylab="x2", main="alpine01", xaxt="n", yaxt="n",
         cex.lab=2, cex.main=2, colkey=F )
contour2D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
         lwd=1.5, nlevels=7, cex.lab=2, cex.main=2, colkey=F, col="black", add=T )
ys <- apply(XX,1,bird)
image2D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
         xlab="x1", ylab="x2", main="bird", xaxt="n", yaxt="n",
         cex.lab=2, cex.main=2, colkey=F )
contour2D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
           lwd=1.5, nlevels=15, cex.lab=2, cex.main=2, colkey=F, col="black", add=T )
# ys <- apply(XX,1,levy03)
# image2D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
#          xlab="x1", ylab="x2", main="levy03", xaxt="n", yaxt="n",
#          cex.lab=2, cex.main=2, colkey=F )
contour2D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
           lwd=1.5, nlevels=11, cex.lab=2, cex.main=2, colkey=F, col="darkgrey", add=T )
ys <- apply(XX,1,michalewicz)
image2D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
         xlab="x1", ylab="x2", main="michalewicz", xaxt="n", yaxt="n",
         cex.lab=2, cex.main=2, colkey=F )
contour2D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
           lwd=1.5, nlevels=10, cex.lab=2, cex.main=2, colkey=F, col="black", add=T )
ys <- apply(XX,1,styblinskiTang)
image2D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
         xlab="x1", ylab="x2", main="styblinskyTang", xaxt="n", yaxt="n",
         cex.lab=2, cex.main=2, colkey=F )
contour2D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
           lwd=1.5, nlevels=15, cex.lab=2, cex.main=2, colkey=F, col="darkgrey", add=T )
ys <- apply(XX,1,ursem03)
image2D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
         xlab="x1", ylab="x2", main="ursem03", xaxt="n", yaxt="n",
         cex.lab=2, cex.main=2, colkey=F )
contour2D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
           lwd=1.5, nlevels=10, cex.lab=2, cex.main=2, colkey=F, col="black", add=T )
# plot.new()
ys <- apply(XX,1,ursemWaves)
image2D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
         xlab="x1", ylab="x2", main="ursemWaves", xaxt="n", yaxt="n",
         cex.lab=2, cex.main=2, colkey=F )
contour2D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys,length(unique(XX$x1))),
           lwd=1.5, nlevels=15, cex.lab=2, cex.main=2, colkey=F, col="black", add=T )
# plot.new()




# **********************************************************************
# multi-fidelity test problems
# **********************************************************************

par(mfrow=c(1,1))
par(mar=c(4.1,4.1,2.1,1.1))
ys1 <- forrester(xs,1)
ys2 <- forrester(xs,2)
ys3 <- forrester(xs,3)
ys4 <- forrester(xs,4)
plot( xs, ys1, type="l", lwd=4, col="black", ylim=range(ys1,ys2,ys3,ys4),
      xlab="x", ylab="", main="forrester", cex.lab=2, cex.axis=2, cex.main=2 )
lines( xs, ys2, col="green3", lwd=4, lty=2 )
lines( xs, ys3, col="orange", lwd=4, lty=4 )
lines( xs, ys4, col="red", lwd=4, lty=5 )
legend( "topleft", legend=c(expression(f[1](x)),
                            expression(f[2](x)),
                            expression(f[3](x)),
                            expression(f[4](x))),
        lwd=4, col=c("black","green3","orange","red"),
        lty=c(1,2,4,5), cex=1.5 )


invisible(readline("Reduce chart window's size and press [RETURN]"))


par(mar=rep(2,4))
par(mfrow=c(1,3))

ys1 <- apply(XX,1,rosenbrock,1)
ys2 <- apply(XX,1,rosenbrock,2)
ys3 <- apply(XX,1,rosenbrock,3)
persp3D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys1,length(unique(XX$x1))),
         xlab="x1", ylab="x2", zlab="f(x)", main="", contour=list(nlevels=30),
         zlim=range(-3000,ys1,ys2,ys3), col="white", border="blue", lwd=2, 
         resfac=0.15, cex.lab=2, cex.main=2, colkey=F )
legend("bottom",legend=expression(f[1](x)), col="blue", lwd=4, cex=2, horiz=T, bty="n" )
persp3D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys2,length(unique(XX$x1))),
         xlab="x1", ylab="x2", zlab="f(x)", main="", contour=list(nlevels=30),
         zlim=range(-3000,ys1,ys2,ys3), col="white", border="green3", lwd=2, 
         resfac=0.15, cex.lab=2, cex.main=2, colkey=F )
legend("bottom",legend=expression(f[2](x)), col="green3", lwd=4, cex=2, horiz=T, bty="n" )
# legend("top",legend="multi-fidelity rosenbrock",cex=2,bty="n")
persp3D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys1,length(unique(XX$x1))),
         xlab="x1", ylab="x2", zlab="f(x)", main="", contour=list(nlevels=30),
         zlim=range(-3000,ys1,ys2,ys3), col="white", border="red", lwd=2, 
         resfac=0.15, cex.lab=2, cex.main=2, colkey=F )
legend("bottom",legend=expression(f[3](x)), col="red", lwd=4, cex=2, horiz=T, bty="n" )



par(mar=rep(2,4))
par(mfrow=c(1,3))

ys1 <- apply(XX,1,shiftedRotatedRastrigin,1)
ys2 <- apply(XX,1,shiftedRotatedRastrigin,2)
ys3 <- apply(XX,1,shiftedRotatedRastrigin,3)
persp3D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys1,length(unique(XX$x1))),
         xlab="x1", ylab="x2", zlab="f(x)", main="", contour=list(nlevels=10),
         zlim=range(-5,ys1,ys2,ys3), col="white", border="blue", lwd=2, 
         resfac=0.15, cex.lab=2, cex.main=2, colkey=F )
legend("bottom",legend=expression(f[1](x)), col="blue", lwd=4, cex=2, horiz=T, bty="n" )
persp3D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys2,length(unique(XX$x1))),
         xlab="x1", ylab="x2", zlab="f(x)", main="", contour=list(nlevels=10),
         zlim=range(-5,ys1,ys2,ys3), col="white", border="green3", lwd=2, 
         resfac=0.15, cex.lab=2, cex.main=2, colkey=F )
legend("bottom",legend=expression(f[2](x)), col="green3", lwd=4, cex=2, horiz=T, bty="n" )
# legend("top",legend="shiftedRotatedRastrigin",cex=2,bty="n")
persp3D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys1,length(unique(XX$x1))),
         xlab="x1", ylab="x2", zlab="f(x)", main="", contour=list(nlevels=10),
         zlim=range(-5,ys1,ys2,ys3), col="white", border="red", lwd=2, 
         resfac=0.15, cex.lab=2, cex.main=2, colkey=F )
legend("bottom",legend=expression(f[3](x)), col="red", lwd=4, cex=2, horiz=T, bty="n" )



par(mar=rep(3,4))
par(mfrow=c(1,2))

ys1 <- apply(XX,1,heterogeneous,1)
ys2 <- apply(XX,1,heterogeneous,2)

persp3D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys1,length(unique(XX$x1))),
         xlab="x1", ylab="x2", zlab="f(x)", main="", contour=list(nlevels=10),
         zlim=range(-2,ys1,ys2), col="white", border="blue", lwd=2, 
         resfac=0.15, cex.lab=1.5, cex.main=1.5, colkey=F )
legend("topleft",legend=expression(f[1](x)), col="blue", lwd=4, cex=1.5 )
persp3D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys2,length(unique(XX$x1))),
         xlab="x1", ylab="x2", zlab="f(x)", main="", contour=list(nlevels=10),
         zlim=range(-2,ys1,ys2), col="white", border="green3", lwd=2, 
         resfac=0.15, cex.lab=1.5, cex.main=2, colkey=F )
legend("topleft",legend=expression(f[2](x)), col="green3", lwd=4, cex=1.5 )




par(mar=rep(3,4))
par(mfrow=c(1,2))

ys1 <- apply(XX,1,paciorek,1)
ys2 <- apply(XX,1,paciorek,2)

persp3D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys1,length(unique(XX$x1))),
         xlab="x1", ylab="x2", zlab="f(x)", main="", contour=list(nlevels=10),
         zlim=range(0,ys1,ys2), col="white", border="blue", lwd=2, 
         resfac=0.15, cex.lab=1.5, cex.main=1.5, colkey=F )
legend("topleft",legend=expression(f[1](x)), col="blue", lwd=4, cex=1.5 )
persp3D( x=sort(unique(XX$x1)), y=sort(unique(XX$x2)), z=matrix(ys2,length(unique(XX$x1))),
         xlab="x1", ylab="x2", zlab="f(x)", main="", contour=list(nlevels=10),
         zlim=range(0,ys1,ys2), col="white", border="green3", lwd=2, 
         resfac=0.15, cex.lab=1.5, cex.main=2, colkey=F )
legend("topleft",legend=expression(f[2](x)), col="green3", lwd=4, cex=1.5 )

