
plot.iNEXT2<-function(x, type=1, se=TRUE, SITE=NULL, show.legend=TRUE, show.main=TRUE, col=NULL, pch=NULL, add=F, ...){
  
  if(!inherits(x, "iNEXT"))
    stop("invalid object class")
  TYPE <-  c(1, 2, 3)
  # SPLIT <- c("none", "order", "site", "both")
  if(is.na(pmatch(type, TYPE)) | pmatch(type, TYPE) == -1)
    stop("invalid plot type")
  
  type <- pmatch(type, 1:3)
  
  y <- method <- site <- shape <- y.lwr <- y.upr <- NULL
  site <<- NULL
  
  z <- fortify(x, type=type)
  
  
  if("y.lwr" %in% names(z) == FALSE & se) {
    warning("invalid se setting, the iNEXT object do not consist confidence interval")
    se <- FALSE
  }else if("y.lwr" %in% names(z) & se) {
    se <- TRUE
  }else{
    se <- FALSE
  }
  
  if(type==1L) {
    #z$x <- z[,1]
    #z$y <- z$qD
    if(!is.null(xlab)) xlab <- ifelse(names(x$DataInfo)[2]=="n", "Number of individuals", "Number of sampling units")
    if(!is.null(ylab)) ylab <- "Species diversity (Hill coefficient)"
    # if(se){
    #   z$y.lwr <- z[,5]
    #   z$y.upr <- z[,6]
    # }
  }else if(type==2L){
    if(length(unique(z$Order.q))>1){
      # z <- subset(z, Order.q==unique(z$Order.q)[1])
      z <- z[z$Order.q==unique(z$Order.q)[1],]
    }
    # z$x <- z[,1]
    # z$y <- z$SC
    if(!is.null(xlab)) xlab <- ifelse(names(x$DataInfo)[2]=="n", "Number of individuals", "Number of sampling units")
    if(!is.null(ylab)) ylab <- "Sample coverage"
    # if(se){
    #   z$y.lwr <- z[,8]
    #   z$y.upr <- z[,9]
    # }
  }else if(type==3L){
    # z$x <- z$SC
    # z$y <- z$qD
    if(!is.null(xlab)) xlab <- "Sample coverage"
    if(!is.null(ylab)) ylab <- "Species diversity (Hill coefficient)"
    # if(se){
    #   z$y.lwr <- z[,5]
    #   z$y.upr <- z[,6]
    # }
  }
  
  conf.reg=function(x,LCL,UCL,...) {
    x.sort <- order(x)
    x <- x[x.sort]
    LCL <- LCL[x.sort]
    UCL <- UCL[x.sort]
    polygon(c(x,rev(x)),c(LCL,rev(UCL)), ...)
  }
  
  
  if (is.null(SITE)){
    SITE <- unique(z$Assemblage)    
  }
  
  ORDER <- unique(z$Order.q)
  
  if(is.null(col)){
    col <- gg_color_hue(length(SITE))
  }else{
    col <- rep(col,length(SITE))[1:length(SITE)]
  }
  if (is.null(pch)){
    pch <- (16+1:length(SITE))%%25
  }else{
    pch = rep(pch, length(SITE))
  }
  j=type
  # for(j in 1:length(ORDER)){
  if(se==TRUE){
    # tmp.sub <- subset(z, Order.q==ORDER[j])
    tmp.sub <- z[z$Order.q==ORDER[j],]
    tmp.j <- data.frame(Assemblage=tmp.sub$Assemblage,
                        Order.q=tmp.sub$Order.q,
                        Method=tmp.sub$Method, 
                        x=tmp.sub$x,
                        y=tmp.sub$y,
                        y.lwr=tmp.sub$y.lwr,
                        y.upr=tmp.sub$y.upr)
    if (add==F){
      plot(y.lwr~x, data=tmp.j[tmp.j$Method=="Extrapolation", ], type="n", xlab="", ylab="", ...)
      grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = 0.5)
      
    }
    
  }else{
    # tmp.sub <- subset(z, Order.q==ORDER[j])
    tmp.sub <- z[z$Order.q==ORDER[j],]
    
    tmp.j <- data.frame(Assemblage=tmp.sub$Assemblage, Order.q=tmp.sub$Order.q,
                        Method=tmp.sub$Method, 
                        x=tmp.sub$x, y=tmp.sub$y)
    if (add==F){
      plot(y~x, data=tmp.j, type="n", xlab="", ylab="", ...)
      grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = 0.5)
    }
  }
  
  for(i in 1:length(SITE)){
    # tmp <- subset(tmp.j, Assemblage==SITE[i])
    tmp <- tmp.j[tmp.j$Assemblage==SITE[i],]
    if(se==TRUE){
      conf.reg(x=tmp[tmp$Method=="Extrapolation",]$x,
               LCL=tmp[tmp$Method=="Extrapolation",]$y.lwr,
               UCL=tmp[tmp$Method=="Extrapolation",]$y.upr,
               border=NA, col=adjustcolor(col[i], 0.1))
    }
    # lines(y~x, data=subset(tmp, Method=="Rarefaction"), lty=1, lwd=2, col=col[i])
    # lines(y~x, data=tmp[tmp$Method=="Rarefaction",], lty=1, lwd=2, col=col[i])
    # lines(y~x, data=subset(tmp, Method=="Extrapolation"), lty=2, lwd=2, col=col[i])
    lines(y~x, data=tmp[tmp$Method=="Extrapolation",], lty=2, lwd=1, col=col[i])
    # points(y~x, data=subset(tmp, Method=="Observed"), pch=pch[i], cex=2, col=col[i])
    points(y~x, data=tmp[tmp$Method=="Observed",], pch=pch[i], cex=2, lwd=2, col=col[i])
    
  }
  if(show.legend==TRUE){
    if(type==3L){
      legend("topleft", legend=paste(SITE), col=col, lty=1, lwd=2, pch=pch, cex=1, bty="n")
    }else{
      legend("bottomright", legend=paste(SITE), col=col, lty=1, lwd=2, pch=pch, cex=1, bty="n")
    }
  }
  #title(ylab = ylab, line=2)
  #title(xlab = xlab, line=2.5)    
  if(show.main==TRUE) title(main=paste("Order q =", ORDER[j]))
  # }
}

gg_color_hue<-function(n){
    hues<-seq(15, 375, length = n+1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
