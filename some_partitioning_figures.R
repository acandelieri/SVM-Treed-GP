rm(list=ls()); graphics.off(); cat("\014")

source("test_functions.R")
source("svmtgp.R")

tf <- 4 # 4 - levy; 8 - ursem_waves
kern <- "matern3_2"
seed <- 1
n.grid <- 200
iters.chart <- c(0,30,60,90)



res <- readRDS( paste0("results/",test.functions[[tf]]$name,"_",kern,"_",seed,".RDS") ) 

par(mfrow=c(2,2))
# linSVMTGP
for( i in iters.chart ) {
  ixs <- which( res$method=="linSVMTGP" | res$iter<=i )
  X <- res[ixs,4:5]
  linSVMTGP <- svmtgp( X, res$y[ixs], thr=ncol(X), linear=T, covtype=kern )
  # plot2D.labels( linSVMTGP, lower=test.functions[[tf]]$lower, upper=test.functions[[tf]]$upper, n.grid=n.grid )
  plot2D.mean( linSVMTGP, lower=test.functions[[tf]]$lower, upper=test.functions[[tf]]$upper, n.grid=n.grid )
  title( paste0((i+10)," function evaluations" ))
  points2D( X[,1], X[,2], pch=19, col="blue", add=T )
  points2D( test.functions[[tf]]$x.star[1], test.functions[[tf]]$x.star[2], pch=8, lwd=2, col="red", add=T )
}

# nolinSVMTGP
for( i in iters.chart ) {
  ixs <- which( res$method=="linSVMTGP" | res$iter<=i )
  X <- res[ixs,4:5]
  nolinSVMTGP <- svmtgp( X, res$y[ixs], thr=ncol(X), linear=F, covtype=kern )
  # plot2D.labels( nolinSVMTGP, lower=test.functions[[tf]]$lower, upper=test.functions[[tf]]$upper, n.grid=n.grid )
  plot2D.mean( nolinSVMTGP, lower=test.functions[[tf]]$lower, upper=test.functions[[tf]]$upper, n.grid=n.grid )
  title( paste0((i+10)," function evaluations" ))
  points2D( X[,1], X[,2], pch=19, col="blue", add=T )
  points2D( test.functions[[tf]]$x.star[1], test.functions[[tf]]$x.star[2], pch=8, lwd=2, col="red", add=T )
}
par(mfrow=c(1,1))