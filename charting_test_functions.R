rm(list=ls()); graphics.off(); cat("\014")

source("test_functions.R")
library(plot3D)

n.grid <- 100


for( tf in 1:length(test.functions) ) {
  Xs <- cbind(sort(rep(seq(test.functions[[tf]]$lower[1],test.functions[[tf]]$upper[1],length.out=n.grid), n.grid)),
              rep(seq(test.functions[[tf]]$lower[2],test.functions[[tf]]$upper[2],length.out=n.grid), n.grid))
  ys <- numeric(nrow(Xs))
  for( i in 1:nrow(Xs) )
    ys[i] <- test.functions[[tf]]$f(Xs[i,])
  image2D( z=matrix(ys,nrow=n.grid,byrow=T), x=sort(unique(Xs[,1])),
           y=sort(unique(Xs[,2])), col=cm.colors(100), colkey=F,
           xlab=expression(x[1]), ylab=expression(x[2]), cex.axis=1.5, cex.lab=1.5,
           main=test.functions[[tf]]$name )
  contour2D( z=matrix(ys,nrow=n.grid,byrow=T), x=sort(unique(Xs[,1])), y=sort(unique(Xs[,2])), col="grey", nlevels=20, add=T  )
  points2D( test.functions[[tf]]$x.star[1],test.functions[[tf]]$x.star[2], pch=8, lwd=2, col="red", add=T )
  
  persp3D( x=sort(unique(Xs[,1])), y=sort(unique(Xs[,2])),
           z=matrix(ys,nrow=n.grid,byrow=T), col=cm.colors(100), colkey=F,
           xlab="x1", ylab="x2", zlab="f(x)", cex.axis=1.5,
           phi=40, theta=40, 
           main=test.functions[[tf]]$name, border="grey")
  
}
