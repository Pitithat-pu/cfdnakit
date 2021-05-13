#' Write solution table to file
#'
#' @param cnvcall cnvcalling result from function call_cnv.R
#' @param outputfile Character; Destination file name
#'
#' @export
#'
#' @examples
write_solution_table <- function(cnvcall, outputfile){
  if(!dir.exists(dirname(outputfile))) {
    stop("Output directory doesn't exist.")
  }
  solution_table = get_solution_table(cnvcall)
  write.table(solution_table, outputfile,sep="\t",quote = FALSE,
              row.names = FALSE,col.names = TRUE)
  paste0("Done writing table: ",basename(outputfile))
}

#' Get summarised table of cnv solutions
#'
#' @param cnv_solutions cnvcalling result from function call_cnv.R
#'
#' @return Dataframe of solution table
#' @export
#'
#' @examples
get_solution_table <- function(cnv_solutions) {
  solution_list = lapply(cnv_solutions, function(cnv_solution){
    cnv_solution[["cnv_solution"]]
  })
  as.data.frame(do.call(rbind,solution_list))
}
