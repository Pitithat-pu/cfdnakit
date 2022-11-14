#' Plot Distance Matrix from CNVCalling
#'
#' @param cnvcall cnvcalling result from function call_cnv.R
#'
#' @return ggplot object ; distance matrix per cnvcalling solution
#' @export
#' @importFrom ggplot2 .data
#'
#' @examples
#' ### Loading example SampleBam file
#' example_file =  system.file("extdata","example_patientcfDNA.RDS",package = "cfdnakit")
#' sample_bambin <- readRDS(example_file)
#' ### Example PoN
#' PoN_rdsfile = system.file("extdata","ex.PoN.rds",package = "cfdnakit")
#' pon_profiles = readRDS(PoN_rdsfile)
#' sample_profile <- get_fragment_profile(sample_bambin,sample_id = "Patient1")
#'
#' sample_zscore = get_zscore_profile(sample_profile,pon_profiles)
#' sample_zscore_segment = segmentByPSCB(sample_zscore)
#'
#' sample_cnv = call_cnv(sample_zscore_segment,sample_zscore, tfs=c(0.1,0.3),ploidies=c(1.5,2), MaxCN=3)
#' plot_distance_matrix(sample_cnv)
#'
#'
plot_distance_matrix <- function(cnvcall) {
  predicted_ployd_tf = as.data.frame(cnvcall[["ploidy_tf_mat"]])
  distance_per_ploidy = get_solution_table(cnvcall)
  # View(distance_per_ploidy)
  tile_plot = ggplot2::ggplot(predicted_ployd_tf,
                              ggplot2::aes_(~ tf,~ploidy)) +
    ggplot2::geom_tile(ggplot2::aes_(fill = ~average_distance))
  tile_plot = tile_plot +
    ggplot2::geom_text(
      data = distance_per_ploidy,
      ggplot2::aes(x=.data$TF, y= .data$ploidy,
                   label=paste0("*",.data$rank)),size=10)

  tile_plot = tile_plot +
    ggplot2::scale_fill_gradientn(
      colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF"),
      name="Average Distance")
  tile_plot = tile_plot +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0))
  tile_plot = tile_plot +
    ggplot2::theme_bw()
  tile_plot = tile_plot +
    ggplot2::theme(axis.text = ggplot2::element_text(size=14),
                   axis.title = ggplot2::element_text(size = 14),
                   legend.title = ggplot2::element_text(size=14)) +
    ggplot2::xlab("Tumor Fraction") + ggplot2::ylab("Ploidy")
  return(tile_plot)
}
