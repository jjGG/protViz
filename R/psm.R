#R
# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/psm.R $
# $Id: psm.R 6936 2014-11-07 08:15:20Z cpanse $
# $Date: 2014-11-07 09:15:20 +0100 (Fri, 07 Nov 2014) $


# TODO 
# compute score by sum of error div. by number of hits

psm<-function(sequence, 
    spec, 
    FUN=defaultIon,
    fi=fragmentIon(sequence, FUN=FUN)[[1]],
    fragmentIonError=0.6) { 

    n<-nchar(sequence)

    pim<-fi$y[nrow(fi)]

    # consider only b and y ions
    # by.mZ<-c(fi$b, fi$y)
    # by.label<-c(paste("b",1:n,sep=''), paste("y",1:n,sep=''))

    by.mZ<-numeric()
    by.label<-character()
    fi.names<-names(fi)

    for (i in 1:ncol(fi)){
        by.mZ <- c(by.mZ, fi[,i])
        by.label <- c(by.label, paste(fi.names[i],1:n,sep=''))
    }

    out <- .C("findNN_",
        nbyion=as.integer(length(by.mZ)),
        nmZ=as.integer(length(spec$mZ)),
        byion=as.double(by.mZ),
        mZ=as.double(spec$mZ),
        NN=as.integer(rep(-1, length(by.mZ))))


    mZ.error<-spec$mZ[out$NN+1] - by.mZ

    res <- list(mZ.Da.error=mZ.error, 
        mZ.error = mZ.error,
        mZ.ppm.error= 1E+6 * mZ.error / by.mZ,
        idx=out$NN+1,
        label=by.label, 
        score=-1, 
        pim=pim,
        fragmentIonError=fragmentIonError,
        sequence=sequence,
        fragmentIon=fi)

    class(res) <- "psm"

    return(res)
}


plot.psm <- function(x, ...){
        plot(x$mZ.error[ mZ.error.idx <- order(x$mZ.error) ],
            main=paste("Error of", x$sequence, "(parent ion mass =", round(x$pim,2) ,"Da)"),
            ylim=c(-5 * x$fragmentIonError, 5 * x$fragmentIonError),
            pch='o',
            sub=paste('The error cut-off is', 
                x$fragmentIonError, 'Da (grey line).'), ...)

        abline(h=x$fragmentIonError, col='grey')
        abline(h=-x$fragmentIonError, col='grey')
        abline(h=0, col='grey', lwd=2)

        text(1:length(x$label), 
            x$mZ.error[mZ.error.idx],  
            x$label[mZ.error.idx],
            cex=0.75, 
            pos=3) 

        hits <- (abs(x$mZ.error) < x$fragmentIonError)
        nHits <- sum(hits)

        sumMZerror <- round(sum(abs(x$mZ.error[hits])),2)

        avgMZerror <- round(sumMZerror / nHits, 2)
        cover <- round(nHits / (nrow(x$fragmentIon) * ncol(x$fragmentIon)), 2)

        legend("topleft", paste(c('nHits','sumMZerror','avgMZerror','cover'),
            as.character(c(nHits, sumMZerror, avgMZerror, cover)),sep='=')) 
}
