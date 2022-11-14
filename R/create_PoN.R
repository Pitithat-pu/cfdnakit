#' Create Panel-of-Normal (PoN) object
#'
#' @param list_rdsfiles Character; a file contains paths to Profile.Rdata per line

#' @return Null
#' @export
#'
#' @examples
#' healthy.1 = system.file("extdata","ex.healthy1.rds",package = "cfdnakit")
#' healthy.2 = system.file("extdata","ex.healthy2.rds",package = "cfdnakit")
#'
#' path_to_PoN_txt = paste0(system.file("extdata",package = "cfdnakit"),"/temp.reference_healthy.listfile")
#' fileConn<-file(path_to_PoN_txt)
#' writeLines(c(healthy.1,healthy.2), fileConn)
#' close(fileConn)
#'
#'
#' PoN.profiles = create_PoN(path_to_PoN_txt)
#' file.remove(path_to_PoN_txt)

create_PoN <- function(list_rdsfiles){
  pon_profile_list = read_PoN_files(list_rdsfiles)
  control_SL_ratio_lst =
    lapply(pon_profile_list, function(pon_profile){
      per_bin_SLRatio = pon_profile$per_bin_profile$`S/L.Ratio.corrected`
      names(per_bin_SLRatio) = rownames(pon_profile$per_bin_profile)
      return(per_bin_SLRatio)
    })
  control_SLratio_df = do.call(cbind,control_SL_ratio_lst)
  control_SLratio_transform =
    apply(control_SLratio_df, 2, function(sample_SL_ratio){
      sample_SL_df = data.frame("SLRatio_corrected"=sample_SL_ratio)
      rownames(sample_SL_df) = rownames(control_SLratio_df)
      zscore_transform(sample_SL_df)
    })
  return(control_SLratio_transform)
}




#' Read Fragment Profile from a list of rds file
#'
#' @param list_rdsfiles path to file containing list of rds file
#'
#' @return list containing content of rds file
#'
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
