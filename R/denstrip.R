
denstrip <- function(x, dens, at, width, horiz=TRUE, colmax=par("fg"), scale=1, gamma=1, 
                     ticks=NULL, tlen=1.5, twd=par("lwd"), mticks=NULL, mlen=1.5, mwd=par("lwd")*2, 
                     ...)
{
    if (!is.numeric(x)) stop("\'x\' must be numeric")
    if (missing(dens)) {
        ## x assumed to be a sample from a distribution, density estimated
        de <- density(x, ...)
        x <- de$x; dens <- de$y
    }
    else {
        if (!is.numeric(dens)) stop("\'dens\' must be numeric")
        if (length(dens) != length(x)) stop("Lengths of \'dens\' and \'x\' must be the same")
        dens <- dens[order(x)]
        x <- sort(x)
    }
    dens <- dens / max(dens) * scale
    n <- length(x)
    rgbmax <- col2rgb(colmax, alpha=TRUE)
    if (gamma <= 0) stop("gamma must be greater than 0")
    p <- dens[1:(n-1)] ^ gamma
    cols <- rgb(p*rgbmax[1] + (1 - p)*255,
                p*rgbmax[2] + (1 - p)*255,
                p*rgbmax[3] + (1 - p)*255, alpha=rgbmax[4], maxColorValue=255)
    first.col <- c(TRUE, cols[2:(n-1)] != cols[1:(n-2)])
    next.col <- c(first.col,TRUE); next.col[1] <- FALSE    
    if (missing(width))
        width <- diff(par("usr")[if(horiz) 3:4 else 1:2]) / 30
    if (horiz) {
        xleft <- x[-n][first.col]; xright=x[next.col];
        ybottom <- at-width/2; ytop <- at+width/2
    }
    else {
        xleft <- at-width/2; xright <- at+width/2;
        ybottom <- x[-n][first.col]; ytop <- x[next.col]
    }
    rect(xleft=xleft, ybottom=ybottom, xright=xright, ytop=ytop,
         border=NA, col = cols[first.col])
    if (!is.null(ticks)){
        if (horiz) { tx0 <- tx1 <- ticks; ty0 <- at-width*tlen/2; ty1 <- at+width*tlen/2 }
        else { tx0 <- at-width*tlen/2; tx1 <- at+width*tlen/2; ty0 <- ty1 <- ticks }
        segments(tx0, ty0, tx1, ty1, lwd=twd)
    }
    if (!is.null(mticks)){
        if (horiz) { tmx0 <- tmx1 <- mticks; tmy0 <- at-width*mlen/2; tmy1 <- at+width*mlen/2  }
        else { tmx0 <- at-width*mlen/2; tmx1 <- at+width*mlen/2; tmy0 <- tmy1 <- mticks }
        segments(tmx0, tmy0, tmx1, tmy1, lwd=mwd)
    }
    invisible()
}

denstrip.normal <- function(mean, sd, log=FALSE, nx=1000, ...){    
    x <-
        if (log) qlnorm(seq(0,1,length=nx)[-c(1,nx)], mean, sd)
        else qnorm(seq(0,1,length=nx)[-c(1,nx)], mean, sd)
    dens <- if (log) dlnorm(x, mean, sd) else dnorm(x, mean, sd)
    denstrip(x=x, dens=dens, ...)
}

denstrip.legend <- function(x, # central x position 
                            y, # central y position 
                            width, len, colmax=par("fg"), gamma=1, horiz=FALSE,
                            max, nticks=5, value.adj=0, main="Density")
{
    if (missing(width))
        width <- diff(par("usr")[if(horiz) 3:4 else 1:2]) / 30
    if (missing(len))
        len <- diff(par("usr")[if(horiz) 1:2 else 3:4]) / 4
    if (horiz) {        
        pt <- x; at <- y
        xdim <- len; ydim <- width
    }
    else {
        pt <- y; at <- x
        xdim <- width; ydim <- len
    }
    npoints <- 1000 # number of distinct colors. no need to make this an argument
    dx <- seq(pt - len/2, pt + len/2, length=npoints)
    ddens <- seq(0, max, length=npoints)
    denstrip(x=dx, dens=ddens, at=at, width=width, colmax=colmax, gamma=gamma, horiz=horiz)    
    polygon(x = c(x - xdim/2, x + xdim/2, x + xdim/2, x - xdim/2), # draw box around strip
            y = c(y - ydim/2, y - ydim/2, y + ydim/2, y + ydim/2))
    ticks.at <- seq(pt - len/2, pt + len/2, length=nticks)
    tick.nos <- round(seq(0, max, length=nticks), 2)
    if (horiz) segments(ticks.at, at - width/2, ticks.at, at - width*0.75)
    else segments(at + width/2, ticks.at, at + width*0.75, ticks.at)
    if (horiz) text(ticks.at, at - width - value.adj, tick.nos, cex = par()$cex * 0.75, pos=1)
    else text(at + width + value.adj, ticks.at, tick.nos, cex = par()$cex * 0.75, pos=4)
    text(x, y + ydim/2, main, cex = par()$cex * 0.75, pos=3)
    invisible()
}

densregion <- function(x, ...) UseMethod("densregion")

densregion.default <- function(x, # times we have estimates for (vector)
                               y, # vector of values on y axis to show distinct densities at.  same number of distinct y values for each t.
                               z, # matrix of densities at t times y
                               pointwise=FALSE,
                               nlevels=100, # number of distinct densities to show
                               colmax=par("fg"),
                               scale=1,
                               gamma=1,
                               ...
                               )
{
    if (pointwise)
        z <- z/apply(z, 1, max)
    qz <- unique(quantile(z, seq(0,1,length=nlevels+1)))
    nlevels <- length(qz) - 1
    zq <- cut(z, qz, include.lowest=TRUE, labels=seq(length=nlevels))
    zq <- matrix(as.numeric(zq), nrow=nrow(z), ncol=ncol(z))
    dens <- tapply(z, zq, mean) # unique densities to plot 
    rgbmax <- col2rgb(colmax, alpha=TRUE)
    if (gamma <= 0) stop("gamma must be greater than 0")
    p <- (dens / max(dens)) ^ gamma
    cols <- rgb(p*rgbmax[1] + (1 - p)*255,
                p*rgbmax[2] + (1 - p)*255,
                p*rgbmax[3] + (1 - p)*255, alpha=rgbmax[4], maxColorValue=255)
    z <- z[order(x),order(y)]
    x <- sort(x)
    y <- sort(y)
    .Internal(filledcontour(as.double(x), as.double(y), scale*z/max(z), qz/max(z), col = cols))
    invisible()
}

densregion.survfit <- function(x, ny=20, ...) {
    if (!inherits(x, "survfit")) stop(deparse(substitute(x)), " should be a survfit object")
    if (is.null(x$conf.type))
        stop("No confidence intervals available in \"", deparse(substitute(x)), "\"")
    identity <- function(x)x # this is in base R from 2.7 
    tr <- switch(x$conf.type, "log-log"=function(x)log(-log(x)), "log" = log, plain=identity)
    invtr <- switch(x$conf.type, "log-log"=function(x)exp(-exp(x)), "log" = exp, plain=identity)
    drop <- is.na(x$surv) | is.na(x$lower)
    surv <- x$surv[!drop]; lower <- x$lower[!drop]; time <- x$time[!drop]
    lsurv <- tr(surv)
    lse <- abs((tr(surv) - tr(lower)) / qnorm(0.975))
    n <- length(time)
    if (n==0) stop("No observed events")
    y <- matrix(nrow=n, ncol=ny) 
    ## "ny" ordinates, based on normal quantiles, where density must be calculated for each event time
    for (i in 1:n)
        y[i,] <- invtr(qnorm(seq(0,1,length=ny+2)[-c(1,ny+2)], lsurv[i], lse[i]))
    ## but actually calculate density at all y for every time, defining a grid on the whole plot region
    yy <- sort(unique(y))
    z <- matrix(nrow=n, ncol=length(yy))
    for (i in 1:n)
        z[i,] <- dnorm(tr(yy), lsurv[i], lse[i])
    densregion.default(x=time, y=yy, z=z, ...)
    invisible()
}

densregion.normal <- function(x, mean, sd, ny=20, ...)
{
    n <- length(x) 
    if (n != length(mean))
        stop(deparse(substitute(mean)), " and ", deparse(substitute(x)), " should be the same length")
    if (n != length(sd))
        stop(deparse(substitute(sd)), " and ", deparse(substitute(x)), " should be the same length")
    y <- matrix(nrow=n, ncol=ny)
    for (i in 1:n)
        y[i,] <- qnorm(seq(0,1,length=ny+2)[-c(1,ny+2)], mean[i], sd[i])
    yy <- sort(unique(y))
    z <- matrix(nrow=n, ncol=length(yy))
    for (i in 1:n)
        z[i,] <- dnorm(yy, mean[i], sd[i])    
    densregion.default(x=x, y=yy, z=z, ...)
    invisible()
}

seqToIntervals <- function(x){
    x <- sort(unique(as.integer(x)))
    breaks <- x[c(1, diff(x)) > 1]
    groups <- cut(x, c(min(x)-1, breaks - 0.5, max(x)+1), labels = FALSE)
    ranges <- tapply(x, groups, range)
    res <- do.call("rbind", ranges)
    colnames(res) <- c("from","to")
    return(res)
}

sectioned.density <- function(x, dens, at, width, offset, ny,
                              method=c("kernel","frequency"), nx, horiz=TRUE,
                              colmax=par("fg"), colmin=par("bg"), gamma=1, ...)
{
    if (missing(width))
        width <- diff(par("usr")[if(horiz) 3:4 else 1:2]) / 20
    if (!is.numeric(x)) stop("\'x\' must be numeric")
    if (missing(offset)) offset <- width/3
    if (!missing(dens)) {
        if (!is.numeric(dens)) stop("\'dens\' must be numeric")
        if (length(dens) != length(x)) stop("Lengths of \'dens\' and \'x\' must be the same")
        de <-  list(x=sort(x), y=dens[order(x)])
    }
    else { 
        method <- match.arg(method)
        if (method=="frequency") {
            if (missing(nx)) nx <- nclass.Sturges(x)
            xcuts <- seq(min(x), max(x), length=nx)
            dens <- table(cut(x, xcuts)) / length(x)
            de <- list(x=sort(x), y=dens[findInterval(x, xcuts,rightmost.closed=TRUE)][order(x)])
        }
        else if (method=="kernel") de <- density(x, ...)
    }
    if (missing(ny)) ny <- nclass.Sturges(de$y)
    ycuts <- seq(0, max(de$y), length=ny+1)
    rgbmax <- col2rgb(colmax, alpha=TRUE)
    rgbmin <- col2rgb(colmin, alpha=TRUE)
    if (gamma <= 0) stop("gamma must be greater than 0")
    p <- seq(0,1,length=ny-1) ^ {1 / gamma}
    cols <- rgb(rgbmax[1] + p*(rgbmin[1] - rgbmax[1]),
                rgbmax[2] + p*(rgbmin[2] - rgbmax[2]),
                rgbmax[3] + p*(rgbmin[3] - rgbmax[3]), alpha=rgbmax[4], maxColorValue=255)
    for (i in 2:ny){
        ## draw rectangles for each region of x with density greather than cut-off
        ## one for each contiguous block in ind 
        ind <- which(de$y >= ycuts[i])
        ints <- seqToIntervals(ind) # from R.utils
        for (j in seq(length=nrow(ints))) {
            if (horiz) 
                rect(xleft=de$x[ints[j,1]],
                     xright=de$x[ints[j,2]],
                     ybottom=at + (i-1)*offset,
                     ytop=at + (i-1)*offset + width,
                     col=cols[i-1], border = (if(i==ny) cols[i-1] else NA) )
            else
                rect(ybottom=de$x[ints[j,1]],
                     ytop=de$x[ints[j,2]],
                     xleft=at - (i-1)*offset,
                     xright=at - (i-1)*offset - width,
                     col=cols[i-1], border = (if(i==ny) cols[i-1] else NA) )
        }
    }
    invisible()
}

cistrip <- function(x, at, d, horiz=TRUE, pch = 16, ...)
{
    if (missing(d)) d <- diff(par("usr")[if(horiz) 3:4 else 1:2]) / 60
    if (!is.numeric(x)) stop("\'x\' must be numeric")
    if (is.vector(x)) x <- matrix(x, ncol=3)
    n <- nrow(x)
    if (length(at) != n) stop("length of \'at\' should equal the number of estimates in \'x\'")
    for (i in 1:n) { 
        if (horiz) { 
            points(x[i,1], at[i], pch=pch, ...)
            segments(x[i,2], at[i], x[i,3], at[i], ...)
            segments(x[i,2], at[i]-d/2, x[i,2], at[i]+d/2, ...)
            segments(x[i,3], at[i]-d/2, x[i,3], at[i]+d/2, ...)
        }
        else {
            points(at[i], x[i,1], pch=pch, ...)
            segments(at[i], x[i,2], at[i], x[i,3], ...)
            segments(at[i]-d/2, x[i,2], at[i]+d/2, x[i,2], ...)
            segments(at[i]-d/2, x[i,3], at[i]+d/2, x[i,3], ...)
        }
        invisible()
    }
}

vwstrip <- function(x, dens, at, width, horiz=TRUE, scale=1, limits=c(-Inf,Inf),
                    col="gray", border=NULL, lwd=par("lwd"), lty=par("lty"),
                    ticks=NULL, tlen=1, twd=par("lwd"), tty=par("lty"), ...) 
{
    if (!is.numeric(x)) stop("\'x\' must be numeric")
    if (missing(dens)) {
        ## x assumed to be a sample from a distribution, density estimated
        de <- density(x, ...)
        x <- de$x; dens <- de$y
    }
    else {
        if (!is.numeric(dens)) stop("\'dens\' must be numeric")
        if (length(dens) != length(x)) stop("Lengths of \'dens\' and \'x\' must be the same")
        dens <- dens[order(x)]
        x <- sort(x)
    }
    dens <- dens[x>limits[1] & x<limits[2]]
    x <- x[x>limits[1] & x<limits[2]]
    if (missing(width))
        width <- diff(par("usr")[if(horiz) 3:4 else 1:2]) / 20
    dens <- dens / max(dens) * width/2 * scale
    xx <- c(x, rev(x))
    yy <- c(dens, -rev(dens)) + at 
    if (horiz) polygon(xx, yy, col=col, border=border, lty=lty, lwd=lwd)
    else polygon(yy, xx, col=col, border=border, lty=lty, lwd=lwd)
    if (!is.null(ticks)) {
        dticks <- dens[findInterval(ticks, x)]
        if (horiz) { tx0 <- tx1 <- ticks; ty0 <- at - dticks*tlen; ty1 <- at + dticks*tlen }
        else { ty0 <- ty1 <- ticks; tx0 <- at - dticks*tlen; tx1 <- at + dticks*tlen }
        segments(tx0, ty0, tx1, ty1, lwd=twd, lty=tty)
    }
    invisible()
}

vwstrip.normal <- function(mean, sd, log=FALSE, nx=1000, ...){    
    x <-
        if (log) qlnorm(seq(0,1,length=nx)[-c(1,nx)], mean, sd)
        else qnorm(seq(0,1,length=nx)[-c(1,nx)], mean, sd)
    dens <- if (log) dlnorm(x, mean, sd) else dnorm(x, mean, sd)
    vwstrip(x=x, dens=dens, ...)
}

bpstrip <- function(x, prob, at, width, horiz=TRUE, scale=1, limits=c(-Inf,Inf),
                    col="gray", border=NULL, lwd=par("lwd"), lty=par("lty"),
                    ticks=NULL, tlen=1, twd=par("lwd"), tty=par("lty"))
{
    if (!is.numeric(x)) stop("\'x\' must be numeric")
    if (missing(prob)) {
        x <- sort(x)
        prob <- ecdf(x)(x)
    }
    else {
        if (!is.numeric(prob)) stop("\'prob\' must be numeric")
        if (length(prob) != length(x)) stop("Lengths of \'prob\' and \'x\' must be the same")
        prob <- prob[order(x)]
        x <- sort(x)
    }
    prob <- pmin(prob, 1 - prob)
    prob <- prob[x>limits[1] & x<limits[2]]
    x <- x[x>limits[1] & x<limits[2]]
    if (missing(width))
        width <- diff(par("usr")[if(horiz) 3:4 else 1:2]) / 20
    prob <- prob / max(prob) * width/2 * scale
    xx <- c(x, rev(x))
    yy <- c(prob, -rev(prob)) + at 
    if (horiz) polygon(xx, yy, col=col, border=border, lty=lty, lwd=lwd)
    else polygon(yy, xx, col=col, border=border, lty=lty, lwd=lwd)
    if (!is.null(ticks)) {
        pticks <- prob[findInterval(ticks, x)]
        if (horiz) { tx0 <- tx1 <- ticks; ty0 <- at - pticks*tlen; ty1 <- at + pticks*tlen }
        else { ty0 <- ty1 <- ticks; tx0 <- at - pticks*tlen; tx1 <- at + pticks*tlen }
        segments(tx0, ty0, tx1, ty1, lwd=twd, lty=tty)
    }
    invisible()
}
