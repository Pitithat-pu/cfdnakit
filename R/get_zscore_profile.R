#' Transform SLRatio with PoN Fragment profile
#'
#' @param fragment_profile Sample Profile
#' @param pon_profile_rds Path to PoN rds file
#'
#' @return Dataframe of robust transformed SLratio
#' @export
#'
#' @examples
get_zscore_profile <- function(fragment_profile,
                               pon_profile_rds){
  if(! file.exists(pon_profile_rds))
    stop(paste0("File ",basename(pon_profile_rds), " doesn't exist."))

  cat("Reading PoN profile",pon_profile_rds,"\n")
  control_SLratio_transform = readRDS(pon_profile_rds)
  # main_functions.rbustz_transform(sample_SL_df,rownames(control_SL_ratio_df))
  sample_SLRatio_transform = zscore_transform(
    dplyr::select(fragment_profile$per_bin_profile,
                  `S/L.Ratio.corrected`))

  SLRatio_zscore_df =
    calculate_SLRatio_zscore(sample_SLRatio_transform,
                             control_SLratio_transform)
  SLRatio_zscore_df =
    dplyr::mutate(SLRatio_zscore_df,
                  "SLRatio.corrected" =
                    dplyr::pull(fragment_profile$per_bin_profile,
                                  `S/L.Ratio.corrected`),
                  SLRatio_transform = sample_SLRatio_transform)

}


zscore_transform <- function(per_bin_profile){
  per_bin_profile$region_id = rownames(per_bin_profile)
  ztransformed_SLRatio = apply(per_bin_profile,1,function(bin){
    sample_SLRatios = as.numeric(bin[1])
    reference_regions =
      unique(per_bin_profile$region_id[! per_bin_profile$region_id %in% bin["region_id"]])
    reference_SLRatio =
      as.numeric(per_bin_profile[reference_regions,1])
    SLRatio_zscore =
      ifelse(stats::mad(reference_SLRatio,na.rm = TRUE)!=0,
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
  })))

  # SLRatio_zscore_df =
  #   dplyr::mutate_at(SLRatio_zscore_df,
  #                    c("start", "end","PoN_SLRatio_median",
  #                      "PoN_SLRatio_mad",
  #                      "size_based_zscore"), ~ as.numeric(levels(.))[.])
  SLRatio_zscore_df =
    dplyr::mutate(SLRatio_zscore_df,
                  start = as.numeric(levels(start))[start],
                  end = as.numeric(levels(end))[end],
                  PoN_SLRatio_median = as.numeric(levels(PoN_SLRatio_median))[PoN_SLRatio_median],
                  PoN_SLRatio_mad = as.numeric(levels(PoN_SLRatio_mad))[PoN_SLRatio_mad],
                  size_based_zscore = as.numeric(levels(size_based_zscore))[size_based_zscore])


}
