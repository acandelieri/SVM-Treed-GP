rm( list=ls() ); graphics.off(); cat("\014")

source("svmtgp.R")

f <- function(xx) {
  x1 <- xx[1]
  x2 <- xx[2]
  
  x1bar <- 15*x1 - 5
  x2bar <- 15 * x2
  
  term1 <- x2bar - 5.1*x1bar^2/(4*pi^2) + 5*x1bar/pi - 6
  term2 <- (10 - 10/(8*pi)) * cos(x1bar)
  
  y <- (term1^2 + term2 - 44.81) / 51.95
  return(y)
}

set.seed(42)
N <- 100
n.grid <- 500
linear.kernel <- T

# X <- cbind( runif(N), runif(N) )
X <- lhs::maximinLHS(N,2)
y <- numeric(N)
for( i in 1:N )
  y[i] <- f(X[i,])


model <- svmtgp( X, y, linear=linear.kernel )

Xs <- cbind( sort(rep(seq(0,1,length.out=n.grid),n.grid)),
             rep(seq(0,1,length.out=n.grid),n.grid) )
Ys.m <- numeric(nrow(Xs))
Ys.sd <- numeric(nrow(Xs))
leaves <- numeric(nrow(Xs))

tmp <- svmtgp.predict(model,Xs)

Ys.m <- matrix( tmp$mean, nrow=n.grid, byrow=T )
Ys.sd <- matrix( tmp$sd, nrow=n.grid, byrow=T )
leaves <- matrix( tmp$leaf.id, nrow=n.grid, byrow=T )



par( mfrow=c(1,2) )
curr.mar <- par("mar")
par(mar=c(4.1,4.6,2.1,1.1))
plot2D.labels( model, c(0,0), c(1,1), n.grid=n.grid )
points2D( X[,1], X[,2], pch=19, col="blue", add=T )
plot2D.mean( model, c(0,0), c(1,1), n.grid=n.grid )
points2D( X[,1], X[,2], pch=19, col="blue", add=T )
par(mar=curr.mar)
par( mfrow=c(1,1) )

# par( mfrow=c(2,2) )
# persp3D( x=seq(0,1,length.out=n.grid), y=seq(0,1,length.out=n.grid), z=Ys.m,
#          theta=90, phi=30, colkey=F, xlab="x1", ylab="x2", zlab="f(x)")
# persp3D( x=seq(0,1,length.out=n.grid), y=seq(0,1,length.out=n.grid), z=Ys.m,
#        theta=60, phi=30, colkey=F, xlab="x1", ylab="x2", zlab="f(x)" )
# persp3D( x=seq(0,1,length.out=n.grid), y=seq(0,1,length.out=n.grid), z=Ys.m,
#          theta=30, phi=30, colkey=F, xlab="x1", ylab="x2", zlab="f(x)")
# persp3D( x=seq(0,1,length.out=n.grid), y=seq(0,1,length.out=n.grid), z=Ys.m,
#          theta=0, phi=30, colkey=F, xlab="x1", ylab="x2", zlab="f(x)")
# par( mfrow=c(1,1) )


n.grid <- 50
Xs <- cbind( sort(rep(seq(0,1,length.out=n.grid),n.grid)),
             rep(seq(0,1,length.out=n.grid),n.grid) )
Ys.m <- numeric(nrow(Xs))
Ys.sd <- numeric(nrow(Xs))
leaves <- numeric(nrow(Xs))

tmp <- svmtgp.predict(model,Xs)

Ys.m <- matrix( tmp$mean, nrow=n.grid, byrow=T )
Ys.sd <- matrix( tmp$sd, nrow=n.grid, byrow=T )
leaves <- matrix( tmp$leaf.id, nrow=n.grid, byrow=T )


fig <- persp3D( x=seq(0,1,length.out=n.grid), y=seq(0,1,length.out=n.grid), z=Ys.m, border="darkgrey",
                cex.axis=1.5,
                theta=60, phi=30, colkey=F, xlab="", ylab="", zlab="", col=cm.colors(100) )
text(trans3d(0.0,-0.2,2.0,fig), "f(x)", cex=2, col="black")
text(trans3d(0.7,-0.2,0.0,fig), label= "x1", cex=2, col="black")
text(trans3d(1.3,0.5,0.0,fig), label= "x2", cex=2, col="black", xpd=NA, pos=2)
