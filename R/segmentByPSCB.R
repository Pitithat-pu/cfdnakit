#' Segmentation data with PSCBS
#'
#' @param sample_transformed_sl dataframe of z-transformed SLRatio
#'
#' @return Dataframe of segmentation result
#' @export
#'
#' @examples
#' @importFrom PSCBS dropSegmentationOutliers findLargeGaps gapsToSegments segmentByCBS pruneByHClust getSegments
segmentByPSCB <- function(sample_transformed_sl){
  chrTotalLength_file= "hg19_chrTotalLength.tsv"
  chrLength_file =
    system.file("extdata",
                chrTotalLength_file,
                package = "cfdnakit")
  chrLength_df = read.table(file = chrLength_file,header=F, sep="\t")
  chrLength_info = util.get_chrLength_info(chrLength_df)

  data <- sample_transformed_sl[,c("chrom","start","size_based_zscore")]
  colnames(data) <- c("chromosome","x","y")
  data <- dropSegmentationOutliers(data)
  data[which(data$chromosome=="X"),"chromosome"] <- 23
  gaps <- findLargeGaps(data, minLength = 1e+07)
  knownSegments <- gapsToSegments(gaps)
  fit <- segmentByCBS(data, knownSegments = knownSegments,
                      joinSegments = TRUE, avg = "median")
  fitP <- pruneByHClust(fit, h = 0.50)

  # segment_df = getSegments(fit, simplify = TRUE)
  segment_df = getSegments(fitP, simplify = TRUE)

  # segment_df$pass_cutoff = segment_df$mean >= cutoff_zscore[1] | segment_df$mean <= cutoff_zscore[2]
  segment_df[which(segment_df$chromosome==23),"chromosome"] = "X"
  segment_df = dplyr::select(segment_df,
                           -c("sampleName"))
  segment_df = na.omit(segment_df)
  return(segment_df)
}
