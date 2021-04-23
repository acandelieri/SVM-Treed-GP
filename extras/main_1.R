rm( list=ls() ); graphics.off(); cat("\014")

source("svmtgp.R")

# ******************************************
# Branin rescaled
# ******************************************
f <- function(xx) {
  x1 <- xx[1]
  x2 <- xx[2]

  x1bar <- 15*x1 - 5
  x2bar <- 15*x2

  term1 <- x2bar - 5.1*x1bar^2/(4*pi^2) + 5*x1bar/pi - 6
  term2 <- (10 - 10/(8*pi)) * cos(x1bar)

  y <- (term1^2 + term2 - 44.81) / 51.95
  return(y)
}
lower <- c(0,0)
upper <- c(1,1)
y.star <- f( c( (pi+5)/15, 2.275/15 ) )
# ******************************************

# ******************************************
# cosine mixture
# ******************************************
# f <- function( x ) {
#   f <- -0.1*sum(cos(5*pi*x))+sum(x^2)
#   return(f)
# }
# lower <- c(-1,-1)
# upper <- c(1,1)
# y.star <- f( c(0,0) )
# ******************************************

# ******************************************
# Rosenbrock modified
# ******************************************
# f <- function(x) {
#   f <- 74 + 100 * (x[2]-x[1]^2)^2 + (1-x[1])^2 - 400 * exp( -( (x[1]+1)^2 + (x[2]+1)^2 )/0.1 )
# }
# lower <- c(-2,-2)
# upper <- c(2,2)
# y.star <- f( c(-0.9,-0.95) )
# ******************************************

# ******************************************
# Levy
# ******************************************
# f <- function(xx) {
#   d <- length(xx)
#   w <- 1 + (xx - 1)/4
# 
#   term1 <- (sin(pi*w[1]))^2
#   term3 <- (w[d]-1)^2 * (1+1*(sin(2*pi*w[d]))^2)
# 
#   wi <- w[1:(d-1)]
#   sum <- sum((wi-1)^2 * (1+10*(sin(pi*wi+1))^2))
# 
#   y <- term1 + sum + term3
#   return(y)
# }
# lower <- c(-10,-10)
# upper <- c(10,10)
# y.star <- f( c(1,1) )
# ******************************************

# ******************************************
# Tripod
# ******************************************
# f <- function( x ) {
#   return( (x[2]>=0) + (1+(x[1]>=0)) + abs( x[1] +50*(x[2]>=0)*(1-2*(x[1]>=0)) ) + abs(x[2]+50*(1-2*(x[2]>=0))) )
# }
# lower <- c(-100,-100)
# upper <- c(100,100)
# y.star <- f( c(0,-50) )
# ******************************************


# ******************************************
# INPUT
# ******************************************
set.seed(42)
n0 <- 10
N <- n0+90
linear.kernel <- T

n.grid <- 200 # just for display
# ******************************************


# Initial design
X <- lhs::maximinLHS(n0,length(upper))
for( h in 1:length(upper) )
  X[,h] <- X[,h] * (upper[h]-lower[h]) + lower[h]

y <- numeric(n0)
for( i in 1:n0 )
  y[i] <- f(X[i,])

# SVM-Treed-GP based BO
perc <- 0
while( nrow(X)<N ) {
  
  if( 100*nrow(X)/N > perc+10 ) {
    perc <- perc+10
    cat("...",perc,"% ",sep="")
    if( perc == 100 )
      cat("\n")
  }
  
  model <- svmtgp( X, y, linear=linear.kernel )
  
  x.next <- next.query( svmtgp=model,
                        lower=lower,
                        upper=upper,
                        N=1000, #10*2,
                        eps=10^-8,
                        M=2000,
                        rho=0.5)
  
  X <- rbind( X, x.next$par )
  y <- c( y, f(x.next$par) )
  
  # par( mfrow=c(1,2) )
  # plot2D.mean( model, c(0,0), c(1,1), n.grid=30 )
  # points2D( X[1:n0,1], X[1:n0,2], pch=19, col="black", add=T )
  # points2D( X[-c(1:n0),1], X[-c(1:n0),2], pch=4, cex=1.3, lwd=2, col="black", add=T )
  # 
  # plot2D.labels( model, c(0,0), c(1,1), n.grid=30 )
  # points2D( X[1:n0,1], X[1:n0,2], pch=19, col="black", add=T )
  # points2D( X[-c(1:n0),1], X[-c(1:n0),2], pch=4, cex=1.3, lwd=2, col="black", add=T )
  # par( mfrow=c(1,1) )
  
}

model <- svmtgp( X, y, linear=linear.kernel )



par( mfrow=c(2,2) )
plot2D.labels( model, lower, upper, n.grid=n.grid )
points2D( X[1:n0,1], X[1:n0,2], pch=19, col="black", add=T )
points2D( X[-c(1:n0),1], X[-c(1:n0),2], pch=4, cex=1.3, lwd=2, col="black", add=T )

plot2D.mean( model, lower, upper, n.grid=n.grid )
points2D( X[1:n0,1], X[1:n0,2], pch=19, col="black", add=T )
points2D( X[-c(1:n0),1], X[-c(1:n0),2], pch=4, cex=1.3, lwd=2, col="black", add=T )

Xs <- cbind( sort(rep(seq(lower[1],upper[1],length.out=n.grid),n.grid)),
             rep(seq(lower[2],upper[2],length.out=n.grid),n.grid))
Ys <- numeric(nrow(Xs))
for( i in 1:nrow(Xs) )
  Ys[i] <- f(Xs[i,])
Ys <- matrix( Ys, nrow=n.grid, byrow=T )
contour2D( seq(lower[1],upper[1],length.out=n.grid), seq(lower[2],upper[2],length.out=n.grid), z=Ys, colkey=F, col="grey" )
points2D( X[1:n0,1], X[1:n0,2], pch=19, col="black", add=T )
points2D( X[-c(1:n0),1], X[-c(1:n0),2], pch=4, cex=1.3, lwd=2, col="black", add=T )



best.seen <- cummin(c(min(y[1:n0]),y[-c(1:n0)]))

# Plot best seen
# plot( 0:(N-n0), best.seen, type="l", lwd=3, col="blue",
#       ylim=c(y.star, max(best.seen)) ) 
# points( 0:(N-n0), cummin(c(min(y[1:n0]),y[-c(1:n0)])), pch=19, col="blue" )
# abline(h=y.star,lty=2,col="red")

plot( 0:(N-n0), (best.seen[1]-best.seen)/(best.seen[1]-y.star), type="l", lwd=3, col="blue",
      ylim=c(0,1) )
points( 0:(N-n0), (best.seen[1]-best.seen)/(best.seen[1]-y.star), pch=19, col="blue" )
abline(h=1,lty=2,col="red")

par( mfrow=c(1,1) )

# stop("END!")
# 
# ## Visualize
# 
# Xs <- cbind( sort(rep(seq(0,1,length.out=n.grid),n.grid)),
#              rep(seq(0,1,length.out=n.grid),n.grid) )
# Ys.m <- numeric(nrow(Xs))
# Ys.sd <- numeric(nrow(Xs))
# leaves <- numeric(nrow(Xs))
# 
# tmp <- svmtgp.predict(model,Xs)
# 
# Ys.m <- matrix( tmp$mean, nrow=n.grid, byrow=T )
# Ys.sd <- matrix( tmp$sd, nrow=n.grid, byrow=T )
# leaves <- matrix( tmp$leaf.id, nrow=n.grid, byrow=T )
# 
# 
# 
# 
# # contour( x=seq(0,1,length.out=n.grid), y=seq(0,1,length.out=n.grid), z=Ys )
# # persp( x=seq(0,1,length.out=n.grid), y=seq(0,1,length.out=n.grid), z=Ys,
# #        theta=55, phi=30 )
# 
# library(plot3D)
# # par(mfrow=c(2,2))
# image2D( x=seq(0,1,length.out=n.grid), y=seq(0,1,length.out=n.grid), z=Ys.m, colkey=F )
# contour2D( x=seq(0,1,length.out=n.grid), y=seq(0,1,length.out=n.grid), z=Ys.m,
#            add=T, col="white", nlevels=30 )
# points2D( X[1:n0,1], X[1:n0,2], pch=19, col="black", add=T )
# points2D( X[-c(1:n0),1], X[-c(1:n0),2], pch=4, cex=1.3, lwd=2, col="black", add=T )
# 
# # image2D( x=seq(0,1,length.out=n.grid), y=seq(0,1,length.out=n.grid), z=Ys.sd )
# # image2D( x=seq(0,1,length.out=n.grid), y=seq(0,1,length.out=n.grid), z=Ys.m-Ys.sd )
# 
# image2D( x=seq(0,1,length.out=n.grid), y=seq(0,1,length.out=n.grid), z=leaves, colkey=F )
# 
# 
# 
# # persp3D( x=seq(0,1,length.out=n.grid), y=seq(0,1,length.out=n.grid), z=Ys.m,
# #         theta=55, phi=30, colkey=F)
# # persp3D( x=seq(0,1,length.out=n.grid), y=seq(0,1,length.out=n.grid), z=Ys.sd,
#          # theta=55, phi=30, colkey=F )
# # par(mfrow=c(1,1))
# 
# 
# # plot( 0:(nrow(X)-n0), cummin( c( min(y[1:n0]),y[-c(1:n0)]) ), type="o", lwd=2 )
