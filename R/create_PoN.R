#' Create Panel-of-Normal (PoN) Rdata
#'
#' @param list_rdsfiles a file contains paths to Profile.Rdata per line
#' @param output_dir Output directory
#' @param filename Output filename
#' @param overwrite Logical: Overwrite existing if result already existed
#'
#' @return
#' @export
#'
#' @examples
create_PoN <- function(list_rdsfiles,
                       output_dir,
                       filename = "PoN_fragment_profiles",
                       overwrite = TRUE){
  pon_profile_list = read_PoN_files(list_rdsfiles)
  control_SL_ratio_lst = lapply(pon_profile_list, function(pon_profile){
    per_bin_SLRatio = pon_profile$per_bin_profile$`S/L.Ratio.corrected`
    names(per_bin_SLRatio) = rownames(pon_profile$per_bin_profile)
    return(per_bin_SLRatio)
  })
  control_SLratio_df = do.call(cbind,control_SL_ratio_lst)
  control_SLratio_transform = apply(control_SLratio_df, 2, function(sample_SL_ratio){
    sample_SL_df = data.frame("SLRatio_corrected"=sample_SL_ratio)
    rownames(sample_SL_df) = rownames(control_SLratio_df)
    zscore_transform(sample_SL_df)
  })

  if(!dir.exists(output_dir)){
    stop("Output directory doesn't exist.")
  }

  output_file = paste0(output_dir,"/",filename,".rds")

  if(file.exists(output_file) ){
    if (is.logical(overwrite)) {
      if (overwrite) {
        file.remove(output_file)
      } else {
        stop(paste0("Output file already exist ",output_file,". Set overwrite=TRUE to overwrite."))
      }
    } else stop("Please provide param overwrite a logical value (TRUE or FALSE)")
  }
  paste0("Saving to ",output_file)
  saveRDS(control_SLratio_transform,
          file = output_file)
  paste0("Done")
}




#' Read Fragment Profile from a list of rds file
#'
#' @param list_rdsfiles path to file containing list of rds file
#'
#' @return list containing content of rds file
#'
#' @examples
read_PoN_files = function(list_rdsfiles){
  if(!file.exists(list_rdsfiles)){
    stop(list_rdsfiles," doesn't not existed")
  }

  pon_files = readLines(list_rdsfiles)

  pon_profiles_lst = list()
  for (pon_file in pon_files) {
    pon_profile = readRDS(pon_file)
    pon_profiles_lst[[pon_profile$Sample.ID]] = pon_profile
  }
  return(pon_profiles_lst)
}


check_PoNProfiles_mad <- function(pon_profiles_lst){
  pon_mad_vec = sapply(pon_profiles_lst, function(pon_profiles){
    return(pon_profiles$sample_profile$Mad)
  })

}
