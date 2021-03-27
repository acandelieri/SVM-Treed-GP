rm(list=ls()); graphics.off(); cat("\014")

library(DiceKriging)
source("test_functions.R")
source("svmtgp.R")

seeds <- 1:30
n0 <- 10
N <- 100

tf <- 1
gp_kernels <- c( "gauss", "exp", "powexp", "matern3_2", "matern5_2" )


cat("***",test.functions[[tf]]$name,"***\n\n")

for( gp_kernel in gp_kernels ) {
  
  for( seed in seeds ) {
  
    results <- data.frame( method = character(),
                           seed = numeric(),
                           iter = numeric(),
                           X = I(list()),
                           y = numeric(),
                           stringsAsFactors = F )
    
      
    cat("> analysing seed #",seed,"...\n",sep="")
  
    set.seed(seed)
    
    # initial design
    X <- lhs::maximinLHS(n0,length(test.functions[[tf]]$lower))
    for( i in 1:ncol(X) )
      X[,i] <- X[,i] * (test.functions[[tf]]$upper[i] - test.functions[[tf]]$lower[i] ) + test.functions[[tf]]$lower[i]
    
    # evaluating design
    y <- numeric(nrow(X))
    for( i in 1:length(y) )
      y[i] <- test.functions[[tf]]$f(X[i,])
    
    
    X.lin <- X
    y.lin <- y
    X.nl <- X
    y.nl <- y
    X.bo <- X
    y.bo <- y
    
    results <- rbind( results, data.frame( method = rep("init",n0),
                                           seed = rep(seed,n0), 
                                           iter = rep(0,n0),
                                           X = list(X),
                                           y = y,
                                           stringsAsFactors = F ) )
    
    
    while( nrow(X.lin)<N ) {
      
      # First linear
      SVMTGP <- NULL
      try( SVMTGP <- svmtgp( X.lin, y.lin, thr=ncol(X.lin), linear=T, covtype=gp_kernel ) )
      if( is.null(SVMTGP) ) {
        cat("Retrying fitting linSVMTGP...\n")
        SVMTGP <- svmtgp( X.lin, y.lin, thr=ncol(X.lin), linear=T, covtype=gp_kernel )
      }
      
      
      x_ <- next.query( SVMTGP, test.functions[[tf]]$lower, test.functions[[tf]]$upper,
                        N=length(test.functions[[tf]]$lower)*10L, eps=10^-8, M=2000L, rho=0.5 )
      
      X.lin <- rbind(X.lin,x_$par)
      y.lin <- c(y.lin,test.functions[[tf]]$f(x_$par))
      
  
      # then non-linear
      SVMTGP <- NULL
      try( SVMTGP <- svmtgp( X.nl, y.nl, thr=ncol(X), linear=F, covtype=gp_kernel ) )
      if( is.null(SVMTGP) ) {
        cat("Retrying fitting nolinSVMTGP...\n")
        SVMTGP <- svmtgp( X.nl, y.nl, thr=ncol(X), linear=F, covtype=gp_kernel )
      }
      
      x_ <- next.query( SVMTGP, test.functions[[tf]]$lower, test.functions[[tf]]$upper,
                        N=length(test.functions[[tf]]$lower)*10L, eps=10^-8, M=2000L, rho=0.5 )
      
      X.nl <- rbind(X.nl,x_$par)
      y.nl <- c(y.nl,test.functions[[tf]]$f(x_$par))
  
      
      # then non-linear
      gp <- NULL
      try( gp <- km( design=data.frame(X.bo), response=y.bo, covtype=gp_kernel, nugget.estim=T, control=list(trace=F) ) )
      if( is.null(gp) ) {
        cat("Retrying fitting GP...\n")
        gp <- km( design=data.frame(X.bo), response=y.bo, covtype=gp_kernel, nugget.estim=T, control=list(trace=F) )
      }
        
      SVMTGP <- list(gp)
      
      x_ <- next.query( SVMTGP, test.functions[[tf]]$lower, test.functions[[tf]]$upper,
                        N=length(test.functions[[tf]]$lower)*10L, eps=10^-8, M=2000L, rho=0.5 )
      
      X.bo <- rbind(X.bo,x_$par)
      y.bo <- c(y.bo,test.functions[[tf]]$f(x_$par))
          
          
      # # just for monitoring...
      # # gap metrics
      # plot( 0:(nrow(X.lin)-n0),
      #       (min(y.lin[1:n0]) - cummin( c( min(y.lin[1:n0]), y.lin[-c(1:n0)] )) ) / (min(y.lin[1:n0]) - test.functions[[tf]]$y.star ),
      #       type="o", lwd=3, col="blue", ylab="gap metric", xlab="BO iterations",
      #       cex.axis=1.5, cex.lab=1.5, cex.main=1.5,
      #       main=paste0(test.functions[[tf]]$name," (seed:",seed,") - gp.kern:",gp_kernel),
      #       ylim=0:1)
      # lines( 0:(nrow(X.nl)-n0),
      #        (min(y.nl[1:n0]) - cummin( c( min(y.nl[1:n0]), y.nl[-c(1:n0)] )) ) / (min(y.nl[1:n0]) - test.functions[[tf]]$y.star ),
      #        type="o",col="green3", lwd=3 )
      # lines( 0:(nrow(X.bo)-n0),
      #        (min(y.bo[1:n0]) - cummin( c( min(y.bo[1:n0]), y.bo[-c(1:n0)] )) ) / (min(y.bo[1:n0]) - test.functions[[tf]]$y.star ),
      #        type="o",col="orange3", lwd=3 )
      # abline( h=1.0, lty=2, col="red", lwd=2 )
      # if( nrow(X.lin)%%10 == 0 ) {
      #   graphics.off()
      # }
      
    }
    
    results <- rbind( results, data.frame( method = rep("linSVMTGP",N-n0),
                                           seed = rep(seed,N-n0), 
                                           iter = 1:(N-n0),
                                           X = list(X.lin[(n0+1):N,]),
                                           y = y.lin[(n0+1):N],
                                           stringsAsFactors = F ) )
    
    results <- rbind( results, data.frame( method = rep("nolinSVMTGP",N-n0),
                                           seed = rep(seed,N-n0), 
                                           iter = 1:(N-n0),
                                           X = list(X.nl[(n0+1):N,]),
                                           y = y.nl[(n0+1):N],
                                           stringsAsFactors = F ) )
    
    results <- rbind( results, data.frame( method = rep("BO",N-n0),
                                           seed = rep(seed,N-n0), 
                                           iter = 1:(N-n0),
                                           X = list(X.bo[(n0+1):N,]),
                                           y = y.bo[(n0+1):N],
                                           stringsAsFactors = F ) )
    
    
    if( !dir.exists("results") )
      dir.create("results")
    saveRDS( results, paste0("results/",test.functions[[tf]]$name,"_",gp_kernel,"_",seed,".RDS") )
  }
}




# tf <- as.integer(readline("Test function id:"))
# source("test_functions.R")
# tmp <- results[which(results$method=="init"),] 
# tmp <- aggregate( tmp$y, by=list(tmp$seed), min )
# 
# to.plot.lin <- tmp$x
# to.plot.nl <- tmp$x
# to.plot.bo <- tmp$x
# 
# for( i in min(results$iter):max(results$iter) ) {
#   to.plot.lin <- cbind( to.plot.lin, results$y[which(results$method=="linSVMTGP" & results$iter==i)] )
#   to.plot.nl <- cbind( to.plot.nl, results$y[which(results$method=="nolinSVMTGP" & results$iter==i)] )
#   to.plot.bo <- cbind( to.plot.bo, results$y[which(results$method=="BO" & results$iter==i)] )
# }
# 
# 
# lin.m <- apply( apply( to.plot.lin, 1, cummin ), 1, mean )
# lin.sd <- apply( apply( to.plot.lin, 1, cummin ), 1, sd )
# lin.m <- (lin.m[1] - lin.m) / (lin.m[1] - test.functions[[tf]]$y.star)
# lin.sd <- (lin.sd[1] - lin.sd) / (lin.sd[1] - test.functions[[tf]]$y.star)
# 
# nl.m <- apply( apply( to.plot.nl, 1, cummin ), 1, mean )
# nl.sd <- apply( apply( to.plot.nl, 1, cummin ), 1, sd )
# nl.m <- (nl.m[1] - nl.m) / (nl.m[1] - test.functions[[tf]]$y.star)
# nl.sd <- (nl.sd[1] - nl.sd) / (nl.sd[1] - test.functions[[tf]]$y.star)
# 
# bo.m <- apply( apply( to.plot.bo, 1, cummin ), 1, mean )
# bo.sd <- apply( apply( to.plot.bo, 1, cummin ), 1, sd )
# bo.m <- (bo.m[1] - bo.m) / (bo.m[1] - test.functions[[tf]]$y.star)
# bo.sd <- (bo.sd[1] - bo.sd) / (bo.sd[1] - test.functions[[tf]]$y.star)
# 
# plot( 0:(length(lin.m)-1), lin.m, type="l", lwd=3, col="blue", ylim=c(0,1.2),
#       ylab="Gap metric", xlab="iterations", cex.axis=1.5, cex.lab=1.5)
# polygon( c(0:(length(lin.m)-1),(length(lin.m)-1):0),
#          c( lin.m+lin.sd, rev(lin.m-lin.sd) ),
#          col=adjustcolor("blue",alpha.f=0.1), border=NA )
# polygon( c(0:(length(nl.m)-1),(length(nl.m)-1):0),
#          c( nl.m+nl.sd, rev(nl.m-nl.sd) ),
#          col=adjustcolor("green3",alpha.f=0.1), border=NA )
# polygon( c(0:(length(bo.m)-1),(length(bo.m)-1):0),
#          c( bo.m+bo.sd, rev(bo.m-bo.sd) ),
#          col=adjustcolor("orange3",alpha.f=0.1), border=NA )
# lines( 0:(length(lin.m)-1), lin.m, lwd=3, col="blue" )
# lines( 0:(length(nl.m)-1), nl.m, lwd=3, col="green3" )
# lines( 0:(length(bo.m)-1), bo.m, lwd=3, col="orange3" )
# abline(h=1,lwd=2,lty=2,col="red")




