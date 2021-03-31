test.functions <- list()

test.functions[[1]] <- list(
  name = "branin_rescaled",
  f = function(xx) {
    x1 <- xx[1]
    x2 <- xx[2]

    x1bar <- 15*x1 - 5
    x2bar <- 15*x2

    term1 <- x2bar - 5.1*x1bar^2/(4*pi^2) + 5*x1bar/pi - 6
    term2 <- (10 - 10/(8*pi)) * cos(x1bar)

    y <- (term1^2 + term2 - 44.81) / 51.95
    return(y)
  },
  lower = c(0,0),
  upper = c(1,1),
  x.star = c( (pi+5)/15, 2.275/15 ),
  y.star = -1.047394
)

test.functions[[2]] <- list(
  name = "cosine_mixture",
  f = function( x ) {
    f <- -0.1*sum(cos(5*pi*x))+sum(x^2)
    return(f)
  },
  lower = c(-1,-1),
  upper = c(1,1),
  x.star = c(0,0),
  y.star = -0.2
)

test.functions[[3]] <- list(
  name = "rosenbrock_modified",
  f = function(x) {
    return( 74 + 100 * (x[2]-x[1]^2)^2 + (1-x[1])^2 - 400 * exp( -( (x[1]+1)^2 + (x[2]+1)^2 )/0.1 ) )
  },
  lower = c(-2,-2),
  upper = c(2,2),
  x.star = c(-0.9,-0.95),
  y.star = 34.37124
)

test.functions[[4]] <- list(
  name = "levy03",
  f = function(xx) {
    d <- length(xx)
    w <- 1 + (xx - 1)/4
    
    term1 <- (sin(pi*w[1]))^2
    term3 <- (w[d]-1)^2
    
    wi <- w[1:(d-1)]
    term2 <- sum((wi-1)^2 * (1+10*(sin(pi*wi[-1]))^2))
    
    y <- term1 + term2 + term3
    return(y)
  },
  lower = c(-10,-10),
  upper = c(10,10),
  x.star = c(1,1),
  y.star = 0
)

test.functions[[5]] <- list(
  name = "tripod",
  f = function( x ) {
    return( (x[2]>=0)*(1+(x[1]>=0)) + abs( x[1] +50*(x[2]>=0)*(1-2*(x[1]>=0)) ) + abs(x[2]+50*(x[1]>=0)*(1-2*(x[2]>=0))) )
  },
  lower = c(-100,-100),
  upper = c(100,100),
  x.star = c(0,-50),
  y.star = 0
)


test.functions[[6]] <- list(
  name = "qing",
  f = function( x ) {
    return( sum((x^2 - 1:2)^2) )
  },
  lower = c(-500,-500),
  upper = c(500,500),
  x.star = sqrt(1:2),
  y.star = 0
)


test.functions[[7]] <- list(
  name = "ursem01",
  f = function( x ) {
    return( -sin(2*x[1]-0.5*pi) -3*cos(x[2])-0.5*x[1] )
  },
  lower = c(-2.5,-2),
  upper = c(3,2),
  x.star = c(1.69714,0),
  y.star = -4.8168
)


test.functions[[8]] <- list(
  name = "ursem_waves",
  f = function( x ) {
    return( -0.9*x[1]^2 + (x[2]^2-4.5*x[2]^2)*x[1]*x[2] + 4.7*cos(3*x[1]-(x[2]^2)*(2+x[1]))*sin(2.5*pi*x[1]) )
  },
  lower = c(-0.9,-1.2),
  upper = c(1.2,1.2),
  x.star = c(1.2,1.2),
  y.star = -8.5536
)



# test.functions[[6]] <- list(
#   name = "step",
#   f = function( x ) {
#     return( sum( (floor(x) + 0.5)^2 ) )
#   },
#   lower = c(-100,-100),
#   upper = c(100,100),
#   x.star = c(0.5,0.5),
#   y.star = 0
# )