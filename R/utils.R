utils.file_exists <- function(files_vector){
  all_exist = TRUE
  for (file in files_vector) {
    if(!file.exists(file)){
      paste(file,"do not existed",sep=" ")
      all_exist = FALSE
    }else{
      paste(file,"existed",sep=" ")
    }
  }
  return(all_exist)
}


utils.unlist <- function (x){
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)) {
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}



util.get_sliding_windows <- function(binsize=1000){
  qdnaseq_sliding_windows_RDS =
    system.file("extdata",
                paste0("AnnotationDataFrame_from_QDNAseq_",binsize,"k.rds"),
                package = "cfdnakit")
  if (utils.file_exists(qdnaseq_sliding_windows_RDS)) {
    bins = readRDS(qdnaseq_sliding_windows_RDS)
  } else {
    stop(paste0("The selected binsize (",binsize,") is not available.\nAvailable binsize (kb) are 1000, 500, 100."))

  }
  sliding_windows <- as.data.frame(bins@data)
  sliding_windows <- sliding_windows[which(sliding_windows$chromosome!="Y" &
                                             sliding_windows$mappability>=1),]
  sliding_windows_gr <- GRanges(seqnames = sliding_windows$chrom,
                                ranges = IRanges(
                                  start= sliding_windows$start,
                                  end = sliding_windows$end),
                                gc=sliding_windows$gc,
                                mappability=sliding_windows$mappability)

  return(sliding_windows_gr)
}

util.bias_correct <- function(readcount, bias) {
  i <- seq(min(bias, na.rm=TRUE),
           max(bias, na.rm=TRUE), by = 0.001)
  readcount.trend <- loess(readcount ~ bias)
  readcount.model <- loess(predict(readcount.trend, i) ~ i)
  readcount.pred <- predict(readcount.model, bias)
  readcount.corrected <- readcount - readcount.pred + median(readcount)
}
