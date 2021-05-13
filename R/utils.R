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

#' Correct GC Bias readcount
#'
#' @param readcount numeric
#' @param bias numeric
#'
#' @return numeric
#'
#' @examples
#' @importFrom stats loess predict median
util.bias_correct <- function(readcount, bias) {
  i <- seq(min(bias, na.rm=TRUE),
           max(bias, na.rm=TRUE), by = 0.001)
  readcount.trend <- loess(readcount ~ bias)
  readcount.model <- loess(predict(readcount.trend, i) ~ i)
  readcount.pred <- predict(readcount.model, bias)
  readcount.corrected <- readcount - readcount.pred + median(readcount)
}


#' Load control fragment profile
#'
#' @param control_rds Name of healthy fragment rds file in extdata directory
#' @param label Label of the sample
#'
#' @return List containing fragment length distribution table and all insert size
#' @export
#'
#' @examples
util.load_control_density_table <- function(
  control_rds = "healthycontrol.fragmentprofile.RDS",
  label="Healthy Control"){
  control_RDS_file =
    system.file("extdata",
                control_rds,
                package = "cfdnakit")
  # control_RDS_file = density_control_rds
  control_fragment_profile =
    readRDS(control_RDS_file)
  colnames(control_fragment_profile$distribution_table)=c("x","y")
  control_fragment_profile$distribution_table$Label = label
  return(list("control_density_df"=
                control_fragment_profile$distribution_table,
              "insert_size" =
                control_fragment_profile$insert_size))
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
