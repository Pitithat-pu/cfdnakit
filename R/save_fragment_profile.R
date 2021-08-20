
#' Save Fragment Profile as RDS
#'
#' @param fragment_profile Fragment profile variable from function get_fragment_profile
#' @param output_dir Path to output directory
#' @param overwrite Boolean; TRUE to overwrite existing file
#'
#' @return Null
#' @export
#'
#' @examples
save_fragment_profile <- function(fragment_profile,
                                  output_dir,
                                  overwrite = TRUE){
  if(!dir.exists(output_dir)){
    stop("Output directory doesn't exist.")
  }
  output_file = paste0(output_dir,"/",fragment_profile$Sample.ID,".rds")

  if(file.exists(output_file) ){
    if (is.logical(overwrite)) {
      if (overwrite) {
        file.remove(output_file)
      } else {
        stop(paste0("Output file already exist ",output_file,". Set overwrite=TRUE to overwrite."))
      }
    } else stop("Please provide param overwrite a logical value (TRUE or FALSE)")

  }
  saveRDS(fragment_profile,
          file = output_file)
  paste0("Saving RDS : Done")
}

