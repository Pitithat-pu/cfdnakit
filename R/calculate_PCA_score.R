#' Calculate PCA Score from Segmentation
#'
#' @param sample_segmentation Segmentation Dataframe
#'
#' @return Numeric; PCA score
#' @export
#'
#' @examples
calculate_PCA_score <- function(sample_segmentation){
  PCA_score_SLRatio =
    round(sum(abs(sample_segmentation$mean) * sample_segmentation$nbrOfLoci) /
            nrow(sample_segmentation),3)
  return(PCA_score_SLRatio)
}
