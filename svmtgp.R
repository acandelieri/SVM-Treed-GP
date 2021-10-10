library(DiceKriging)
library(e1071)
library(plot3D)

svmtgp <- function( X, y, thr=ncol(X), linear=T, covtype="gauss" ) {
  
  # TODO
  # args <- list(...) 
  
  if( thr<=ncol(X) ) {
    thr <- ncol(X)+1
    # TODO: print warning!
  }
  
  H <- list()
  L <- list()
  
  L[[1]] <- cbind(X,y)
  J <- c(1)
  
  
  while( length(L)>0 ) {
    
    D <- L[[1]]
    L <- L[-1]

    j <- J[1]
    J <- J[-1]
    
    m <- median( D[,ncol(D)] )
    l <- numeric(nrow(D))+1
    l[ which(D[,ncol(D)]<m) ] <- -1
    if( linear )
      h <- svm( x=D[,-ncol(D)], y=l, type="C-classification", kernel="linear", scale=F )
    else
      h <- svm( x=D[,-ncol(D)], y=l, type="C-classification", kernel="radial", scale=F )
    l.hat <- as.numeric(predict(h,D[,-ncol(D)]))
    
    
    if( length(which(l.hat==1))>thr & length(which(l.hat==2))>thr ) {
      
      H[[j]] <- h
      
      # left child
      D.neg <- D[which(l.hat==1),]
      J <- c(J,2*j)
      L[[length(L)+1]] <- D.neg
      
      # right child
      D.pos <- D[which(l.hat==2),]
      J <- c(J,2*j+1)
      L[[length(L)+1]] <- D.pos
      
    } else {
      
      H[[j]] <- km( design=data.frame(D[,-ncol(D)]), response=D[,ncol(D)],
                    covtype=covtype, nugget.estim=T,
                    control=list(trace=F))
      
    }
    
  }
  
  return( tree=H )
  
}

# svmtgp.predict <- function( svmtgp, x ) {
#   j <- 1
#   while( class(svmtgp[[j]])=="svm" ) {
#     y.hat <- as.numeric(predict( svmtgp[[j]], t(x) ))
#     if( y.hat==1 ) {
#       j <- 2*j  # class -1 --> left child
#     } else {
#       j <- 2*j+1 # class +1 --> right child
#     }
#   }
#   pred <- predict( svmtgp[[j]], data.frame(t(x)), "UK", checkNames=F )
#   return( list(leaf.id=j, gp=pred) )
# }


svmtgp.predict <- function( svmtgp, X ) {
  
  # preliminary check
  if( !is.matrix(X) ) {
    if( ((class(svmtgp[[1]])=="svm" ) && ncol(svmtgp[[1]]$SV)==1) || 
        ((class(svmtgp[[1]])=="km" ) && ncol(svmtgp[[1]]@X)==1) )
      X <- t(X)
  }

  # make predictions 
  lbls <- numeric(nrow(X))+1
  preds.mean <- numeric(nrow(X))+NA 
  preds.sd <- numeric(nrow(X))+NA
  
  to.pred <- 1
  
  while( length(to.pred)>0 ) {
    
    ixs <- which(lbls==to.pred[1])
    X_ <- X[ixs,]
    if( length(ixs)==1 )
      X_ <- t(X_)
    
    if( class(svmtgp[[to.pred[1]]])=="svm" ) {
      
      y.hat <- as.numeric(predict( svmtgp[[to.pred[1]]], X_ ))
      
      if( length(which(y.hat==1))>0 ) {
        # class -1 --> left child
        lbls[ixs[which(y.hat==1)]] <- 2*to.pred[1] 
        to.pred <- c(to.pred, 2*to.pred[1])  
      }
      if( length(which(y.hat==2))>0 ) {
        # class +1 --> right child
        lbls[ixs[which(y.hat==2)]] <- 2*to.pred[1]+1 
        to.pred <- c(to.pred, 2*to.pred[1]+1)  
      }
      
    } else {
      pred <- predict( svmtgp[[to.pred[1]]], data.frame(X_), "UK", checkNames=F )
      preds.mean[ixs] <- pred$mean
      preds.sd[ixs] <- pred$sd
    }
    
    to.pred <- to.pred[-1]
  }

  return( data.frame(leaf.id=lbls, mean=preds.mean, sd=preds.sd) )
}

next.query <- function( svmtgp, lower, upper,
                        N=length(lower)*10L,
                        eps=10^-8, M=2000L, rho=0.5 ) {
  
  
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
      
      y_ <- svmtgp.predict( svmtgp, X_ )
      
      # lower confidence bound
      y_ <- y_$mean - y_$sd
      
      ix <- which.min(y_)
      xc <- X_[ix,]
      y <- y_[ix]
      
    }
    
    n <- n + 1
    l <- l * rho
    
  }
  
  return( list(par=xc, value=y) )
}


plot2D.mean <- function( svmtgp, lower, upper, n.grid = 30 ) {
  
  stopifnot( length(lower)==2 && length(lower)==length(upper) )
  
  Xs <- cbind( sort(rep(seq(lower[1],upper[1],length.out=n.grid),n.grid)),
               rep(seq(lower[2],upper[2],length.out=n.grid),n.grid) )
  Ys.m <- numeric(nrow(Xs))
  Ys.sd <- numeric(nrow(Xs))
  leaves <- numeric(nrow(Xs))
  
  tmp <- svmtgp.predict(svmtgp,Xs)
  
  Ys.m <- matrix( tmp$mean, nrow=n.grid, byrow=T )
  Ys.sd <- matrix( tmp$sd, nrow=n.grid, byrow=T )
  
  image2D( x=seq(lower[1],upper[1],length.out=n.grid), y=seq(lower[2],upper[2],length.out=n.grid), z=Ys.m, colkey=F, xlab="x1", ylab="x2", col=cm.colors(100), cex.axis=2.4, cex.lab=2.4 )
  contour2D( x=seq(lower[1],upper[1],length.out=n.grid), y=seq(lower[2],upper[2],length.out=n.grid), z=Ys.m,
             add=T, col="grey", nlevels=30 )

}

plot2D.labels <- function( svmtgp, lower, upper, n.grid = 30 ) {
  
  stopifnot( length(lower)==2 && length(lower)==length(upper) )
  
  Xs <- cbind( sort(rep(seq(lower[1],upper[1],length.out=n.grid),n.grid)),
               rep(seq(lower[2],upper[2],length.out=n.grid),n.grid) )
  Ys.m <- numeric(nrow(Xs))
  Ys.sd <- numeric(nrow(Xs))
  leaves <- numeric(nrow(Xs))
  
  tmp <- svmtgp.predict(svmtgp,Xs)
  
  leaves <- matrix( tmp$leaf.id, nrow=n.grid, byrow=T )
  
  image2D( x=seq(lower[1],upper[1],length.out=n.grid), y=seq(lower[2],upper[2],length.out=n.grid), z=leaves, colkey=F, xlab="x1", ylab="x2", col=cm.colors(100), cex.axis=2.4, cex.lab=2.4 )
}