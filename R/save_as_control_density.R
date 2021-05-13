#' Save Fragment Profile for distribution comparison
#'
#' @param readbam_bin variable from read_bamfile
#' @param fragment_profile Fragment profile
#' @param file Filename
#'
#' @return
#' @export
#'
#' @examples
save_as_control_density <- function(
  readbam_bin,fragment_profile, file) {
  isize_vector = Biobase::subListExtract(readbam_bin, "isize")
  isize_vector = unlist(isize_vector)
  isize_vector = abs(isize_vector)[which(
    dplyr::between(abs(isize_vector),
                   left = fragment_profile$minimum_length,
                   right = fragment_profile$maximum_length))]
  isize_vector = unname(isize_vector)
  fragment_profile$insert_size = isize_vector

  saveRDS(fragment_profile,file = file)
  paste0("Done")
}
