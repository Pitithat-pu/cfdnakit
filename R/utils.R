
utils.unlist <- function (x){
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)) {
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}



util.get_sliding_windows <- function(binsize=1000, genome="hg19"){
  ### stop if given reference name is neither hg19 nor mm10
  if(! genome %in% c("hg19","hg38","mm10"))
    stop("Only hg19, hg38 or mm10 genome are possible")
  if(! binsize %in% c(100,500,1000))
    stop("The selected binsize (",binsize,
         ") is not available.\nAvailable binsize (kb) are 100, 500, 1000.")

  if(genome %in% c("hg19","hg38")){
    #### Reading in bin info from extdata if hg19 genome
    qdnaseq_sliding_windows_RDS =
      system.file("extdata",
                  paste0("AnnotationDataFrame_from_QDNAseq_",
                         genome,"_",
                         binsize,"k.rds"),
                  package = "cfdnakit")
    if (file.exists(qdnaseq_sliding_windows_RDS)) {
      bins = readRDS(qdnaseq_sliding_windows_RDS)

    } else {
      bins =  QDNAseq::getBinAnnotations(
        binSize=binsize,
        genome = genome)
    }
  } else if (genome=="mm10") {
    ### Loading bin information through QDNAseq package
    message("Loading ",genome, "(",binsize,")")
    bins = QDNAseq::getBinAnnotations(
      binSize=binsize, genome="mm10")
  }



  sliding_windows <- as.data.frame(bins@data)
  sliding_windows <- sliding_windows[
    which(sliding_windows$chromosome!="Y" &
            sliding_windows$mappability>=1),]
  sliding_windows_gr <- GenomicRanges::GRanges(
    seqnames = sliding_windows$chrom,
    ranges = IRanges::IRanges(
      start= sliding_windows$start,
      end = sliding_windows$end),
    gc=sliding_windows$gc,
    mappability=sliding_windows$mappability)

  return(sliding_windows_gr)
}

#' Correct GC Bias readcount
#'
#' @param readcount numeric
#' @param bias numeric
#'
#' @return numeric
#'
#' @importFrom stats loess predict median loess.control
util.bias_correct <- function(readcount, bias) {
  invalid_idx = which(readcount <= 0)


  i <- seq(min(bias, na.rm=TRUE),
           max(bias, na.rm=TRUE), by = 0.001)


  if(length(invalid_idx) == 0 ){
    readcount.trend <- loess(readcount ~ bias,
                             control=loess.control(surface="direct"))
  } else
    readcount.trend <- loess(readcount[-invalid_idx] ~ bias[-invalid_idx],
                             control=loess.control(surface="direct"))

  readcount.model <- loess(predict(readcount.trend, i) ~ i)
  readcount.pred <- predict(readcount.model, bias)
  readcount.corrected <-
    readcount - readcount.pred + median(readcount,na.rm = TRUE)
}



util.rowname_to_columns <- function(per_bin_profile){
  splited_list=lapply(strsplit(
    rownames(per_bin_profile),
    split = "[:-]"), function(x) {
      data.frame("chrom"=x[1],
                 "start"=as.numeric(x[2]),
                 "end"=as.numeric(x[3]))
    })
  cbind(do.call(rbind,splited_list),per_bin_profile)
}


util.get_chrLength_info <- function(chrLength_df){
  colnames(chrLength_df) = c("Chromosome", "Length")
  chrLength_df = chrLength_df[which(chrLength_df$Chromosome!="Y"),]
  chrNames = chrLength_df$Chromosome
  chrLength = chrLength_df$Length
  names(chrLength) = chrNames
  chroffsets = cumsum(as.numeric(chrLength))
  chroffsets <- c(0, chroffsets[0:(length(chroffsets)-1)])
  names(chroffsets) <- names(chrLength)
  chrMids <- cumsum(as.numeric(chrLength))
  chrMids <- (chrMids + chroffsets)/2
  names(chrMids) <- names(chrLength)

  chrLength_info =
    list("chrNames"=chrNames,
         "chrLength"=chrLength,
         "chroffsets"=chroffsets,
         "chrMids"=chrMids)
}
