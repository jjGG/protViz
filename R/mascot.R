#R
# Christian Panse <cp@fgcz.ethz.ch> 20170525

.get_ <- function(obj, attribute = 'query_moverz') {
  as.vector(unlist(lapply(obj$queries, function(x){
    x[attribute]
  })))
}

.get_q_peptide <- function(obj, attribute = 'pep_seq'){
  rv <- unlist(lapply(obj$queries, function(x) { 
    if ("q_peptide" %in% names(x)){
      if (attribute %in% names(x$q_peptide)){
        res <- unlist(x$q_peptide[attribute])
        return(paste(res, collapse = ';'))
      }
    }
    NA
  }))
  as.vector(rv)
}

#' Convert a mascot nested list into a data.frame object
#'
#' @param a mascot object 
#'
#' @author Bernd Roschitzki, 2017
#' @return a data.frame
#' @export
#'
#' @examples
#'  lapply(list(F225712, F225714, F225715, F225716), 
#'     function(x){
#'        S <- as.data.frame.mascot(x); 
#'        plot(S$RTINSECONDS , S$moverz, pch=16, col=rgb(0.4,0.4,0.4,alpha=0.1))
#'     })
#'     
as.data.frame.mascot <- function(x, ...){
  # TODO
  # score 
  # %in%
  # shiny cut-off score
  # reformat charge into integer
  data.frame(RTINSECONDS = as.numeric(.get_(x, attribute = "RTINSECONDS")), 
             moverz = as.numeric(.get_(x, attribute = "query_moverz")), 
             query_charge = .get_(x, attribute = 'query_charge'),
             SCANS = .get_(x, attribute = 'SCANS'),
             TotalIonsIntensity = as.numeric(.get_(x, attribute = 'TotalIonsIntensity')),
             pep_exp_mz = .get_q_peptide(x, attribute = 'pep_exp_mz'),
             pep_exp_mr = .get_q_peptide(x, attribute = 'pep_exp_mr'),
             StringTitle = .get_q_peptide(x, attribute = 'StringTitle'),
             pep_exp_z = as.numeric(.get_q_peptide(x, attribute = 'pep_exp_z')),
             pep_calc_mr = .get_q_peptide(x, attribute = 'pep_calc_mr'),
             pep_delta = .get_q_peptide(x, attribute = 'pep_delta'),
             pep_miss = .get_q_peptide(x, attribute = 'pep_miss'),
             pep_local_mod_pos = .get_q_peptide(x, attribute = 'pep_local_mod_pos'),
             pep_scan_title = .get_q_peptide(x, attribute = 'pep_scan_title'),
             pep_seq = .get_q_peptide(x, attribute = 'pep_seq'),
             pep_expect = as.numeric(.get_q_peptide(x, attribute = 'pep_expect')),
             pep_score = as.numeric(.get_q_peptide(x, attribute = 'pep_score')),
             pep_var_mod = .get_q_peptide(x, attribute = 'pep_var_mod'),
             pep_var_mod_pos = .get_q_peptide(x, attribute = 'pep_var_mod_pos'),
             pep_summed_mod_pos = .get_q_peptide(x, attribute = 'pep_summed_mod_pos'),
             pep_var_mod_conf = .get_q_peptide(x, attribute = 'pep_var_mod_conf')
  )
}


# just define a generic S3 method
as.psm <- function(object, ...){
  UseMethod("as.psm")
}

as.psm.mascot_query <- function(query){
  L <-  .get_ms2(query)
  
  rv <- list(MonoisotopicAAmass = NA,
             charge = as.numeric(gsub("[+]", "", query$query_charge, perl=TRUE)),
             id = 0,
             intensity = L$intensity,
             mZ = L$mZ, 
             mascotScore = NA,
             modification = NA,
             pepmass = as.numeric(query$query_moverz),
             peptideSequence = NA,
             proteinInformation = NA,
             rtinseconds = as.numeric(query$RTINSECONDS),
             scans = query$SCANS, 
             searchEngine = "mascot",
             title = '',
             varModification = NA)

  # TODO: add proteinInformation
  if ('q_peptide' %in% names(query)){
    rv$title <- query$q_peptide$pep_scan_title
    rv$mascotScore <- as.numeric(query$q_peptide$pep_score)
    
    if ( is.null(query$q_peptide$pep_var_mod)){
      rv$modification <- NA
    }else{
      rv$modification <- query$q_peptide$pep_var_mod
    }
    
    rv$peptideSequence <- as.character(query$q_peptide$pep_seq)
    
    # TODO: $q_peptide$pep_var_mod_pos
    rv$varModification <- rep(0.0, nchar(query$q_peptide$pep_seq))
  }

  class(rv) <- c('psm', 'list')
  rv
}


# just define a generic S3 method
as.psmSet <- function(object, ...){
  UseMethod("as.psmSet")
}

#' transformas a mascot object into a psmSet
#'
#' @param mascot obj 
#'
#' @return a psmSet object
as.psmSet.mascot <- function(object, ...){
  
  if (is.mascot(object)){
    # todo(cp): class('psmSet')
    rv <- NULL
    #if(require(parallel)){
    #  rv <- mclapply(obj$queries, as.psm.mascot_query, ...)
    #} else{
    rv <- lapply(object$queries, as.psm.mascot_query)
    #}
    
    # assign the ``query number''
    for (idx in 1:(length(rv)-1))
      rv[[idx]]$id <- idx
    
    class(rv) <- c("psmSet", "list")
    return(rv)
  }
  
  NULL
}

.get_ms2 <- function(query){
  
  
  if (is.character(query$StringIons1)){
    S <- lapply(strsplit(query$StringIons1, ","), function(x){strsplit(x, ':')})[[1]]
  
    mZ <- as.numeric(sapply(S, function(x){x[1]}))
    intensity <- as.numeric(sapply(S, function(x){x[2]}))
    idx <- order(mZ)
    
    return(list(mZ=mZ[idx], intensity=intensity[idx]))
  }
  
  return(list(mZ=NULL, intensity=NULL))
}



summary.mascot <- function(object, ...){
  if (is.mascot(object)){
    cat("number of queries:\n")
    cat(paste("\t", length(object$queries), "\n"))
    cat("quantile pep_score:\n")
    quantile(.get_(object, "pep_score"), na.rm=TRUE)
  }
  #NextMethod('summary')
}

is.mascot <- function(obj){
  if ('mascot' %in% class(obj)){
    if (is.list(obj$queries) & sum(sapply(obj$queries, is.mascot_query)) == length(obj$queries)){
      return(TRUE)
    }
  }
  
  return(FALSE)
}

is.mascot_query <- function(obj){
  mascor_query.names <- c("query_charge", "query_moverz", "SCANS", "StringIons1",
                          "RTINSECONDS", "StringTitle")
  
  if ('mascot_query' %in% class(obj)){
      # all attributes should be present
      if(length(mascor_query.names) == sum(mascor_query.names %in% names(obj))){
        return(TRUE)
      }
  }
  return(FALSE)
}


plot.mascot <- function(x, ...){
  if (is.mascot(x)){
    df <- as.data.frame.mascot(x)
    
    # peptide scores versus e-value
    plot(df$score, 1 / log(df$pep_expect,10), 
         log = 'x', 
         pch = 16, 
         col = rgb(0.1, 0.1, 0.1, alpha = 0.2))
    
    # ssrc
    plot((rt.ssrc <-sapply(df$pep_seq, ssrc)) ~ (rt <- df$RTINSECONDS), 
         pch = 16, 
         col = rgb(0.1, 0.1, 0.1, alpha = 0.2))
    
    abline(lm(rt.ssrc ~ rt), 
           col = rgb(0.9, 0.1, 0.1, alpha=0.3), lwd=3)
    
    # LC-MS map
    lcmsmap(x, ...)
  }
}


plot.mascot_query <- function(x, obj = NULL, FUN=defaultIon, ...){
  if (!is.null(obj)){
    # get AA mass
    x.AA <- do.call('rbind', lapply(obj$masses, function(x){
      data.frame(letter1 = x$.attrs, 
                 Monoisotopic = as.numeric(x$text))
    }))
    
    x.filter <- x.AA$letter1  %in% c("A", "R", "N", "D", "C", "E", "Q", "G",
                                     "H", "I", "L", "K", "M", "F", "P", "S",
                                     "T", "W", "Y", "V") 
    x.AA <- x.AA[x.filter, ]
    # get fixed modification
    
    # get var modification
    x.variable_mods <- sapply(obj$variable_mods, function(x){as.numeric(x$delta)})
    names(x.variable_mods) <- sapply(obj$variable_mods, function(x){x$name})
    
    print(x.variable_mods)
  }
  
  if (is.mascot_query(x)){
    spec <- .get_ms2(x)
    pep_seq <- ""
    
    if ("q_peptide" %in% names(x)){
      if ("pep_seq" %in% names(x$q_peptide)){
        print(pep_seq)
        pep_seq <- x$q_peptide$pep_seq
        x.AA.mass <- aa2mass(pep_seq, letter1 = x.AA$letter1, mass = x.AA$Monoisotopic)[[1]]
        
        if (!is.null(x$q_peptide$pep_var_mod_pos)){
          x.pep_var_mod_pos <- strsplit(x$q_peptide$pep_var_mod_pos, "", perl=TRUE)[[1]]
          x.pep_var_mod_pos.n <- length(x.pep_var_mod_pos)
          x.pep_var_mod_pos <- as.integer(x.pep_var_mod_pos[-c(1, 2, x.pep_var_mod_pos.n - 1, x.pep_var_mod_pos.n)])
          x.pep_var_mod_pos.idx <- which(x.pep_var_mod_pos != 0)
          
          x.AA.mass[x.pep_var_mod_pos.idx] <- x.AA.mass[x.pep_var_mod_pos.idx] + x.variable_mods[x.pep_var_mod_pos[x.pep_var_mod_pos.idx]]
        }
        x.fi <- fragmentIon(x.AA.mass, FUN = FUN)[[1]]
        
        peakplot(peptideSequence = pep_seq, 
                 spec = spec, 
                 fi = x.fi, 
                 sub = x$StringTitle)
      } else {
        plot(spec$mZ, spec$intensity, type='h',  sub = x$StringTitle)
      }
    }
  }
}
  
  