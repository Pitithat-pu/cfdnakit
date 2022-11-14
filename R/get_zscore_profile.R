#' Transform SLRatio with PoN Fragment profile
#'
#' @param fragment_profile Sample Profile
#' @param pon_profile PoN Profiles
#'
#' @return Dataframe of robust transformed SLratio
#' @export
#'
#' @examples
#' ### Loading example SampleBam file
#' example_file =  system.file("extdata","example_patientcfDNA.RDS",package = "cfdnakit")
#' sample_bambin <- readRDS(example_file)
#'
#' ### Example PoN
#' PoN_rdsfile = system.file("extdata","ex.PoN.rds",package = "cfdnakit")
#' pon_profiles = readRDS(PoN_rdsfile)
#' sample_profile <- get_fragment_profile(sample_bambin,sample_id = "Patient1")
#'
#' sample_zscore = get_zscore_profile(sample_profile,pon_profiles)
#'
#' sample_zscore_segment = segmentByPSCB(sample_zscore)

get_zscore_profile <- function(fragment_profile,
                               pon_profile){
  # main_functions.rbustz_transform(sample_SL_df,rownames(control_SL_ratio_df))
  sample_SLRatio_transform = zscore_transform(
    dplyr::select(fragment_profile$per_bin_profile,
                  .data$`S/L.Ratio.corrected`))

  SLRatio_zscore_df =
    calculate_SLRatio_zscore(sample_SLRatio_transform,
                             pon_profile)
  SLRatio_zscore_df =
    dplyr::mutate(SLRatio_zscore_df,
                  "SLRatio.corrected" =
                    dplyr::pull(fragment_profile$per_bin_profile,
                                .data$`S/L.Ratio.corrected`),
                  SLRatio_transform = sample_SLRatio_transform)

}


#' zscore_transform transforms SLRatio profile into z-score
#'
#' @param per_bin_profile SampleFragment from function get_fragment_profile
#'
#' @return dataframe of z-score per bin
#'
#'
#' @importFrom stats mad
zscore_transform <- function(per_bin_profile){
  per_bin_profile$region_id = rownames(per_bin_profile)
  ztransformed_SLRatio = apply(per_bin_profile,1,function(bin){
    sample_SLRatios = as.numeric(bin[1])
    reference_regions =
      unique(per_bin_profile$region_id[! per_bin_profile$region_id %in% bin["region_id"]])
    reference_SLRatio =
      as.numeric(per_bin_profile[reference_regions,1])
    SLRatio_zscore =
      ifelse(mad(reference_SLRatio,na.rm = TRUE)!=0,
             (sample_SLRatios - median(reference_SLRatio,na.rm = TRUE))/mad(reference_SLRatio,na.rm = TRUE),
             (sample_SLRatios - median(reference_SLRatio,na.rm = TRUE)+1)/(mad(reference_SLRatio,na.rm = TRUE)+1))

  })
  names(ztransformed_SLRatio) = rownames(per_bin_profile)
  return(ztransformed_SLRatio)

}



calculate_SLRatio_zscore <- function(sample_SLRatio_transform,
                                     control_SLratio_transform){
  all_region_name = names(sample_SLRatio_transform)
  SLRatio_zscore_df = data.frame(t(sapply(all_region_name, function(x){
    control_SLRatio = unlist(control_SLratio_transform[x,])
    chr_position = strsplit(x,"[:-]+")[[1]]
    control_SLRatio_median = median(control_SLRatio,na.rm = TRUE)
    control_SLRatio_mad = stats::mad(control_SLRatio,na.rm = TRUE)
    sample_SLRatio_transform =
      sample_SLRatio_transform[which(names(sample_SLRatio_transform)==x)]
    SL_zscore = ifelse(control_SLRatio_mad!=0,
                       (sample_SLRatio_transform - control_SLRatio_median)/control_SLRatio_mad,
                       (sample_SLRatio_transform - control_SLRatio_median+1)/(control_SLRatio_mad+1))

    c("chrom" = chr_position[1],
      "start" = as.numeric(chr_position[2]),
      "end" = as.numeric(chr_position[3]),
      "PoN_SLRatio_median"= control_SLRatio_median,
      "PoN_SLRatio_mad" = control_SLRatio_mad,
      "size_based_zscore" = SL_zscore
    )
  })),stringsAsFactors = FALSE)

  SLRatio_zscore_df =
    dplyr::mutate_at(SLRatio_zscore_df,
                     c("start", "end","PoN_SLRatio_median",
                       "PoN_SLRatio_mad",
                       "size_based_zscore"), as.numeric)

}
