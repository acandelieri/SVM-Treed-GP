rm(list=ls()); graphics.off(); cat("\014")

source("test_functions.R")
source("svmtgp.R")

tf <- 4 # 4 - levy; 8 - ursem_waves
kern <- "gauss"
seed <- 1
n.grid <- 200
iters.chart <- c(0,30,60,90)



res <- readRDS( paste0("results/",test.functions[[tf]]$name,"_",kern,"_",seed,".RDS") ) 

par(mfrow=c(2,2))
# linSVMTGP
for( i in iters.chart ) {
  if( i == 0 ) {
    ixs <- which( res$iter==0 )
  } else {
    ixs <- which( res$method=="linSVMTGP" | res$iter<=i )  
  }
  X <- res[ixs,4:5]
  linSVMTGP <- svmtgp( X, res$y[ixs], thr=ncol(X), linear=T, covtype=kern )
  # plot2D.labels( linSVMTGP, lower=test.functions[[tf]]$lower, upper=test.functions[[tf]]$upper, n.grid=n.grid )
  plot2D.mean( linSVMTGP, lower=test.functions[[tf]]$lower, upper=test.functions[[tf]]$upper, n.grid=n.grid )
  title( paste0((i+10)," function evaluations" ), cex.main=1.5 )
  points2D( X[,1], X[,2], pch=20, col="blue", add=T )
  points2D( test.functions[[tf]]$x.star[1], test.functions[[tf]]$x.star[2], pch=8, lwd=2, col="red", add=T )
}

# Xs <- cbind( sort(rep(seq( test.functions[[tf]]$lower[1], test.functions[[tf]]$upper[1], length.out=n.grid),n.grid)),
#              rep(seq( test.functions[[tf]]$lower[2], test.functions[[tf]]$upper[2], length.out=n.grid),n.grid))
# z <- svmtgp.predict( linSVMTGP, Xs )
# z <- matrix( z$mean, nrow=n.grid, byrow=T )
# par(mfrow=c(1,1))
# persp3D( sort(unique(Xs[,1])), sort(unique(Xs[,2])), z, n.grid=n.grid, colkey=F, col=cm.colors(100) )
# par(mfrow=c(2,2))


# nolinSVMTGP
for( i in iters.chart ) {
  if( i == 0 ) {
    ixs <- which( res$iter==0 )
  } else {
    ixs <- which( res$method=="nolinSVMTGP" | res$iter<=i )  
  }
  
  X <- res[ixs,4:5]
  nolinSVMTGP <- svmtgp( X, res$y[ixs], thr=ncol(X), linear=F, covtype=kern )
  # plot2D.labels( nolinSVMTGP, lower=test.functions[[tf]]$lower, upper=test.functions[[tf]]$upper, n.grid=n.grid )
  plot2D.mean( nolinSVMTGP, lower=test.functions[[tf]]$lower, upper=test.functions[[tf]]$upper, n.grid=n.grid )
  title( paste0((i+10)," function evaluations" ), cex.main=1.5 )
  points2D( X[,1], X[,2], pch=20, col="blue", add=T )
  points2D( test.functions[[tf]]$x.star[1], test.functions[[tf]]$x.star[2], pch=8, lwd=2, col="red", add=T )
}

# Xs <- cbind( sort(rep(seq( test.functions[[tf]]$lower[1], test.functions[[tf]]$upper[1], length.out=n.grid),n.grid)),
#              rep(seq( test.functions[[tf]]$lower[2], test.functions[[tf]]$upper[2], length.out=n.grid),n.grid))
# z <- svmtgp.predict( nolinSVMTGP, Xs )
# z <- matrix( z$mean, nrow=n.grid, byrow=T )
# par(mfrow=c(1,1))
# persp3D( sort(unique(Xs[,1])), sort(unique(Xs[,2])), z, n.grid=n.grid, colkey=F, col=cm.colors(100) )
# par(mfrow=c(2,2))


par(mfrow=c(1,1))


