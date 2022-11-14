#' Get summarised table of cnv solutions
#'
#' @param cnv_solutions cnvcalling result from function call_cnv.R
#'
#' @return Dataframe of solution table
#' @export
#'
#' @examples#'
#' ### Loading example SampleBam file
#' example_file =  system.file("extdata","example_patientcfDNA.RDS",package = "cfdnakit")
#' sample_bambin <- readRDS(example_file)
#' ### Example PoN
#' PoN_rdsfile = system.file("extdata","ex.PoN.rds",package = "cfdnakit")
#' pon_profiles = readRDS(PoN_rdsfile)
#' sample_profile <- get_fragment_profile(sample_bambin,sample_id = "Patient1")
#'
#' sample_zscore = get_zscore_profile(sample_profile,pon_profiles)
#'
#' sample_zscore_segment = segmentByPSCB(sample_zscore)
#'
#' sample_cnv = call_cnv(sample_zscore_segment,sample_zscore, tfs=c(0.1,0.3),ploidies=c(1.5,2), MaxCN=3)
#' get_solution_table(sample_cnv)

get_solution_table <- function(cnv_solutions) {
  solution_list = lapply(cnv_solutions, function(cnv_solution){
    cnv_solution[["cnv_solution"]]
  })
  as.data.frame(do.call(rbind,solution_list))
}
