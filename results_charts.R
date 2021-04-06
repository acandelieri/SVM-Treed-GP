rm(list=ls()); graphics.off(); cat("\014")

source("test_functions.R")
kernels <- c("gauss")#,"exp","powexp","matern3_2","matern5_2")

avg.sd <- T

par(mfrow=c(3,3))
par(mar=c(4.1,3.1,3.1,1.1))

for( tf in 1:length(test.functions) ) {
#for( tf in c(1:5,8) ) {
  
  files <- list.files( path="results", pattern=test.functions[[tf]]$name, include.dirs=F )
  
  for( kern in kernels ) {
    bs.lin <- NULL
    bs.nl <- NULL
    bs.bo <- NULL
    ix <- grep( kern, files, fixed=T )
    if( length(ix)>0 ) {
      for( f in files[ix] ) {
        res <- readRDS( paste0("results/",f) ) 
        bs.lin <- rbind( bs.lin, 
                         c( min(res$y[which(res$method=="init")]), res$y[which(res$method=="linSVMTGP")] ) )
        bs.nl <- rbind( bs.nl, 
                         c( min(res$y[which(res$method=="init")]), res$y[which(res$method=="nolinSVMTGP")] ) )
        bs.bo <- rbind( bs.bo, 
                         c( min(res$y[which(res$method=="init")]), res$y[which(res$method=="BO")] ) )
      }
      bs.lin <- t(apply(bs.lin,1,cummin))
      bs.nl <- t(apply(bs.nl,1,cummin))
      bs.bo <- t(apply(bs.bo,1,cummin))
      
      tmp <- matrix(rep(bs.lin[,1],ncol(bs.lin)),ncol=ncol(bs.lin),byrow=F)
      bs.lin <- (tmp - bs.lin) / (tmp - test.functions[[tf]]$y.star)
      tmp <- matrix(rep(bs.nl[,1],ncol(bs.nl)),ncol=ncol(bs.nl),byrow=F)
      bs.nl <- (tmp - bs.nl) / (tmp - test.functions[[tf]]$y.star)
      tmp <- matrix(rep(bs.bo[,1],ncol(bs.bo)),ncol=ncol(bs.bo),byrow=F)
      bs.bo <- (tmp - bs.bo) / (tmp - test.functions[[tf]]$y.star)
      
      if( avg.sd ) {
        lin.avg <- apply(bs.lin,2,mean); lin.sd <- apply(bs.lin,2,sd)
        nl.avg <- apply(bs.nl,2,mean); nl.sd <- apply(bs.nl,2,sd)
        bo.avg <- apply(bs.bo,2,mean); bo.sd <- apply(bs.bo,2,sd)
        
        # plot( lin.avg, type="l", col=NA, main=paste0(test.functions[[tf]]$name,"\nGP's kernel:",kern), ylim=c(0,1.2),
        #       cex.axis=1.5, cex.lab=1.5, cex.main=1.5, xlab="BO iterations", ylab="Gap metric")
        
        plot( lin.avg, type="l", col=NA, main=test.functions[[tf]]$name, ylim=c(0,1),
              cex.axis=1.5, cex.lab=1.5, cex.main=1.5, xlab="BO iterations", ylab="Gap metric")

                dw <- lin.avg-lin.sd; up <- lin.avg+lin.sd
        dw[which(dw<0)] <- 0; up[which(up>1)] <- 1
        polygon( c(0:(length(lin.avg)-1),(length(lin.avg)-1):0), 
                 c(dw,rev(up)),
                 col=adjustcolor("blue",alpha.f=0.1), border=NA )
        dw <- nl.avg-nl.sd; up <- nl.avg+nl.sd
        dw[which(dw<0)] <- 0; up[which(up>1)] <- 1
        polygon( c(0:(length(nl.avg)-1),(length(nl.avg)-1):0), 
                 c(dw,rev(up)),
                 col=adjustcolor("green3",alpha.f=0.1), border=NA )
        dw <- bo.avg-bo.sd; up <- bo.avg+bo.sd
        dw[which(dw<0)] <- 0; up[which(up>1)] <- 1
        polygon( c(0:(length(bo.avg)-1),(length(bo.avg)-1):0), 
                 c(dw,rev(up)),
                 col=adjustcolor("orange",alpha.f=0.1), border=NA )
        lines( 0:(length(lin.avg)-1), lin.avg, col="blue", lwd=3 )
        lines( 0:(length(nl.avg)-1), nl.avg, col="green3", lwd=3 )
        lines( 0:(length(bo.avg)-1), bo.avg, col="orange", lwd=3 ) 
      } else {
        lin.med <- apply(bs.lin,2,median); lin.min <- apply(bs.lin,2,min); lin.max <- apply(bs.lin,2,max)
        nl.med <- apply(bs.nl,2,median); nl.min <- apply(bs.nl,2,min); nl.max <- apply(bs.nl,2,max)
        bo.med <- apply(bs.bo,2,median); bo.min <- apply(bs.bo,2,min); bo.max <- apply(bs.bo,2,max)
        
        plot( lin.med, type="l", col=NA, main=paste0(test.functions[[tf]]$name,"\nGP's kernel:",kern), ylim=c(0,1.2),
              cex.axis=1.5, cex.lab=1.5, cex.main=1.5, xlab="BO iterations", ylab="Gap metric")
        
        # lines( 0:(length(lin.min)-1), lin.min, col="blue", lwd=2, lty=3 )
        # lines( 0:(length(nl.min)-1), nl.min, col="green3", lwd=2, lty=3 )
        # lines( 0:(length(bo.min)-1), bo.min, col="orange", lwd=2, lty=3 )
        # lines( 0:(length(lin.max)-1), lin.max, col="blue", lwd=2, lty=3 )
        # lines( 0:(length(nl.max)-1), nl.max, col="green3", lwd=2, lty=3 )
        # lines( 0:(length(bo.max)-1), bo.max, col="orange", lwd=2, lty=3 )
        polygon( c(0:(length(lin.med)-1),(length(lin.med)-1):0),
                 c(lin.min,rev(lin.max)),
                 col=adjustcolor("blue",alpha.f=0.1), border=NA )
        polygon( c(0:(length(nl.med)-1),(length(nl.med)-1):0),
                 c(nl.min,rev(nl.max)),
                 col=adjustcolor("green3",alpha.f=0.1), border=NA )
        polygon( c(0:(length(bo.med)-1),(length(bo.med)-1):0),
                 c(bo.min,rev(bo.max)),
                 col=adjustcolor("orange",alpha.f=0.1), border=NA )
        lines( 0:(length(lin.med)-1), lin.med, col="blue", lwd=3 )
        lines( 0:(length(nl.med)-1), nl.med, col="green3", lwd=3 )
        lines( 0:(length(bo.med)-1), bo.med, col="orange", lwd=3 )
      }
    }
    # abline( h=1, lty=2, lwd=2, col="red" )
  }
  
  # abline( h=1, lty=2, lwd=2, col="red" )
  # plot.new()
  # legend( "topleft", legend=c("linSVMTGP","nolinSVMTGP","GP"),
  #         col=c("blue","green3","orange"), lwd=3, cex=1.5, bty="n")
  # legend( "top", legend=c("linSVMTGP","nolinSVMTGP","GP"),
  #        col=c("blue","green3","orange"), lwd=3, cex=1.5, horiz=T, bty="n")
}

plot.new()
legend( "topleft", legend=c("linSVMTGP","nolinSVMTGP","GP"),
        col=c("blue","green3","orange"), lwd=3, cex=1.5, bty="n")

par(mar=c(5.1,4.1,4.1,2.1))
par(mfrow=c(1,1))