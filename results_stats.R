rm(list=ls()); graphics.off(); cat("\014")

source("test_functions.R")
kernels <- c("gauss","exp","powexp","matern3_2","matern5_2")

for( tf in 1:length(test.functions) ) {
  
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
      
      final.lin <- round(bs.lin[,ncol(bs.lin)],4)
      final.nl <- round(bs.nl[,ncol(bs.nl)],4)
      final.bo <- round(bs.bo[,ncol(bs.bo)],4)
      
        
      test1.res <- wilcox.test( final.bo, final.lin, alternative="t", paired=F )
      test2.res <- wilcox.test( final.bo, final.nl, alternative="t", paired=F )
      test3.res <- wilcox.test( final.lin, final.nl, alternative="t", paired=F )
      
      cat(test.functions[[tf]]$name,";",kern,";",round(mean(final.bo),4)," (",round(sd(final.bo),5),") ",sep="")
      cat(";",round(mean(final.lin),4)," (",round(sd(final.lin),5),") ",sep="")
      cat(";",round(mean(final.nl),4)," (",round(sd(final.nl),5),") ",sep="")
      cat(";",round(test1.res$p.value,4),";",round(test2.res$p.value,4),";",round(test3.res$p.value,4),"\n",sep="")
      
    }

  }
  
}
