rm(list=ls()); graphics.off(); cat("\014")

library(DiceKriging)
library(tgp)
source("test_functions.R")


next.query.tgp <- function( treedGP, lower, upper, N=length(lower)*10L, eps=10^-8, M=2000L, rho=0.5 ) {
  
  # center of the search space  
  xc <- (upper+lower)/2
  
  n <- 1
  l <- (upper-lower)
  
  
  while( n<=N & all(l>eps) ) {
    
    X_ <- NULL
    
    # replace with lhs?
    for( h in 1:length(upper) ) {
      if( max(lower[h],xc[h]-l/2) < min(upper[h],xc[h]+l/2) ) {
        X_ <- cbind(X_, runif(n=M,min=max(lower[h],xc[h]-l/2),max=min(upper[h],xc[h]+l/2)) )
      } else {
        stop("ERRORE GRAVE!!!!!")
      }
    }
    
    # X_ <- lhs::maximinLHS( n=M, k=length(lower) )
    # if( length(lower)==1 )
    #   X_ <- t(X_)
    # for( h in 1:ncol(X_) ) {
    #   dw <- max(lower[h],xc-l/2)
    #   up <- min(upper[h],xc+l/2)
    #   if( dw>up ) {
    #     dw <- xc
    #     up <- xc
    #   }
    #   X_[,h] <- X_[,h] * (up-dw) + up
    # }
    
    
    if( !is.null(X_) ) {
      
      y_ <- predict( treedGP, X_ )
      
      # lower confidence bound
      y_ <- y_$ZZ.km - y_$ZZ.ks2
      
      ix <- which.min(y_)
      xc <- X_[ix,]
      y <- y_[ix]
      
    }
    
    n <- n + 1
    l <- l * rho
    
  }
  
  return( list(par=xc, value=y) )
}



seeds <- 1:30
n0 <- 10
N <- 100

tfs <- 1:length(test.functions)
gp_kernels <- c( "gauss" ) # one is sufficient, it is just to retrieve initial design


for( tf in tfs ) {
  cat("\014")
  cat("***",test.functions[[tf]]$name,"***\n\n")
  
  for( gp_kernel in gp_kernels ) {
    
    for( seed in seeds ) {
    
      # results <- data.frame( method = character(),
      #                        seed = numeric(),
      #                        iter = numeric(),
      #                        X = I(list()),
      #                        y = numeric(),
      #                        stringsAsFactors = F )
      
        
      cat("> analysing seed #",seed,"...\n",sep="")
    
      set.seed(seed)
      
      # retrieving initial design
      retrieved <- readRDS( paste0("results/",test.functions[[tf]]$name,"_",gp_kernel,"_",seed,".RDS") )
      
      results <- retrieved[which(retrieved[,'method']=='init'),]
      
      X <- results[,4:5]
      y <- results[,'y']
      
      X.tgp <- X
      y.tgp <- y
      
      
      while( nrow(X.tgp)<N ) {
        
        if( gp_kernel=="exp" ) {
          kk <- "exp"
        } else {
          if( gp_kernel=="exp" ) {
            kk <- "exp"
          } else {
            
          }  
        }
        
        TreedGP <- btgp( X.tgp, y.tgp, corr="matern", nu=2, verb=0 )
        x_ <- next.query.tgp( TreedGP, test.functions[[tf]]$lower, test.functions[[tf]]$upper,
                          N=length(test.functions[[tf]]$lower)*10L, eps=10^-8, M=2000L, rho=0.5 )
        
        X.tgp <- rbind(X.tgp,x_$par)
        y.tgp <- c(y.tgp,test.functions[[tf]]$f(x_$par))
        
      }
      
      results <- rbind( results, data.frame( method = rep("treedGP",N-n0),
                                             seed = rep(seed,N-n0), 
                                             iter = 1:(N-n0),
                                             X = list(X.tgp[(n0+1):N,]),
                                             y = y.tgp[(n0+1):N],
                                             stringsAsFactors = F ) )
      
      saveRDS( results, paste0("results_tgp/",test.functions[[tf]]$name,"_",seed,".RDS") )
    }
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




