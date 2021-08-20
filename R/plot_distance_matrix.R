#' Plot Distance Matrix from CNVCalling
#'
#' @param cnvcall cnvcalling result from function call_cnv.R
#'
#' @return ggplot object ; distance matrix per cnvcalling solution
#' @export
#'
#' @examples
plot_distance_matrix <- function(cnvcall) {
  predicted_ployd_tf = as.data.frame(cnvcall[["ploidy_tf_mat"]])
  distance_per_ploidy = get_solution_table(cnvcall)
  # View(distance_per_ploidy)
  tile_plot = ggplot() +
    ggplot2::geom_tile(data=predicted_ployd_tf,
                       ggplot2::aes(tf,ploidy,fill = average_distance))
  tile_plot = tile_plot +
    ggplot2::geom_text(data = distance_per_ploidy,
                       ggplot2::aes(x=TF, y=ploidy,label=paste0("*",rank)),size=10)

  tile_plot = tile_plot +
    ggplot2::scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF"),
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

