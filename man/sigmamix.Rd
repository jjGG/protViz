\docType{data}

\name{sigmamix}

\alias{sigmamix}

\title{MS dilution series of humen peptides.}

\description{

    The data set contains a dilution series of the MS Qual/Quant 
    QC Mix standard sample (SigmaAldrich, Buchs, Switzerland) 
    measured on five different mass spec platforms.

}

\details{

    The mass spectrometers are: 

    \itemize{

    \item \code{Qtrap} - {SRM on a Qtrap 5500}

    \item \code{TSQvantage} - {SRM on a TSQ Vantage}

    \item \code{Qexactive} - {PRM/t-MS2 on a Q-Exactive}

    \item \code{SWATH} - {SWATH on the TripleTOF 5600}

    \item \code{TripleTOF_MS1} - {MS1 quantification on the
    TripleTOF 5600 to evaluate similarities as well as differences in 
    quantitative output.}

    }
    
    For each platform there are 
    three technical replicates of a mixture of 14 SILAC 
    labled tryptic peptide pairs in known amount and a complex yeast background..
    The SILAC peptides span a concentration range over three orders
    of magnitude.


    \code{sigmamix$lf5} contains the absolute measured peak areas for
    each transition for all the SILAC peptide pairs (light and heavy), 
    all machines and
    all dilutions (in triplicates).
    The peak areas were computed using the Skyline software.

   \code{sigmamix$rep5x8} contains the 
    log2 ratios on all five plattforms from
    eight technical replicates.

}


\author{

    Asa Wahlander, Christian Panse, Jonas Grossmann 2014

}

\seealso{

        \url{http://www.sigmaaldrich.com/life-science/proteomics/mass-spectrometry/ms-qual-quant-qc-mix.html}, 2014

        Skyline: \url{http://www.ncbi.nlm.nih.gov/pubmed/20147306}

}

\examples{

    data(sigmamix)

    sigmamix$composition

    cm<-c(grey.colors(38), 'lightgreen', 'darkgreen')

    par(mfrow=c(1,1), mar=c(8,8,1,1));
    image(rep5x8.cor <- cor(sigmamix$rep5x8), 
        oldstyle=TRUE,
        main='cor(sigmamix$rep5x8)', col=cm, zlim=c(0,1),axes=FALSE); 
    axis(1, seq(0,1,length=nrow(rep5x8.cor)), row.names(rep5x8.cor),las=2);  
    axis(2, seq(0,1,length=ncol(rep5x8.cor)), colnames(rep5x8.cor),las=2);

}

\keyword{datasets}
