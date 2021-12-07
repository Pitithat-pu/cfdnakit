#' Calculate PCA Score from Segmentation
#'
#' @param sample_segmentation Segmentation Dataframe
#'
#' @return Numeric; CPA score
#' @export
#'
#' @examples
calculate_CPA_score <- function(sample_segmentation){
  CPA_score_SLRatio =
    round(sum(abs(sample_segmentation$mean) * sample_segmentation$nbrOfLoci) /
            nrow(sample_segmentation),3)
  return(CPA_score_SLRatio)
}
