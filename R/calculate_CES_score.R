#' Calculate CES Score from Segmentation
#'
#' @param sample_segmentation Segmentation Dataframe
#'
#' @return Numeric; CES score
#' @export
#'
#' @examples
#' ### Loading example SampleBam file
#' example_file <-  system.file("extdata","example_patientcfDNA_SampleBam.RDS",package = "cfdnakit")
#' sample_bambin <- readRDS(example_file)
#' ### Example PoN
#' PoN_rdsfile <- system.file("extdata","ex.PoN.rds",package = "cfdnakit")
#' pon_profiles <- readRDS(PoN_rdsfile)
#' sample_profile <- get_fragment_profile(sample_bambin,sample_id = "Patient1")
#'
#' sample_zscore <- get_zscore_profile(sample_profile,pon_profiles)
#'
#' sample_zscore_segment <- segmentByPSCB(sample_zscore)
#'
#' calculate_CES_score(sample_zscore_segment)
#'
calculate_CES_score = function(sample_segmentation){
  CES_score_SLRatio <-
    round(sum(abs(sample_segmentation$mean) * sample_segmentation$nbrOfLoci) /
            nrow(sample_segmentation),3)
  return(CES_score_SLRatio)
}
