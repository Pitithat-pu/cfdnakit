#' Plot z-tranformed Short/Long-fragment Ratio
#'
#' @param sample_transformed_sl Dataframe z-transformed SLRatio from get_zscore_profile
#' @param sample_segment_df Dataframe segmenation from segmentByPSCB
#' @param ylim plot y-axis limit
#'
#' @return Genome-wide plot of z-transformed SLRatio
#' @export
#'
#' @examples
plot_transformed_sl <- function(sample_transformed_sl,
                                sample_segment_df = NULL,
                                ylim=c(-20,20)){
  chrTotalLength_file= "hg19_chrTotalLength.tsv"
  chrLength_file =
    system.file("extdata",
                chrTotalLength_file,
                package = "cfdnakit")
  chrLength_df = read.table(file = chrLength_file,header=F, sep="\t")
  chrLength_info = util.get_chrLength_info(chrLength_df)


  sample_transformed_sl =
    dplyr::mutate(sample_transformed_sl,
                  start = as.numeric(start),
                  end = as.numeric(end),
                  size_based_zscore = as.numeric(size_based_zscore))
  per_bin_profile = sample_transformed_sl

  per_bin_profile =
    dplyr::mutate(per_bin_profile,
                  scaledPos = (start + end)/2 +
                    chrLength_info$chroffsets[chrom])

  sl_plot = ggplot2::ggplot(per_bin_profile,
                            ggplot2::aes(scaledPos,size_based_zscore))
  sl_plot = sl_plot +
    ggplot2::geom_point(fill="grey20", size = 1)
  sl_plot = sl_plot +
    ggplot2::geom_hline(yintercept =
                          median(per_bin_profile$size_based_zscore,
                                 na.rm = TRUE),
                        size=0.5)
  if (!is.null(sample_segment_df)){
    sample_segment_df =
      dplyr::mutate(sample_segment_df,
                    scaledPos_start = start +
                      chrLength_info$chroffsets[chromosome],
                    scaledPos_end= end +
                      chrLength_info$chroffsets[chromosome])
    sl_plot = sl_plot +
      ggplot2::geom_segment(data=sample_segment_df,
                            aes(y=mean, yend=mean,
                                x=scaledPos_start,xend=scaledPos_end),
                            size = 1.5, color="darkorange")
  }
  sl_plot = sl_plot +
    ggplot2::scale_x_continuous(breaks = chrLength_info$chrMids,
                                labels = chrLength_info$chrNames,
                                position = "bottom",expand=c(0,0))
  sl_plot = sl_plot +
    ggplot2::geom_vline(xintercept = chrLength_info$chroffsets,
                        linetype="dotted",size=0.5)
  sl_plot = sl_plot + ggplot2::xlab("Chromosome number") +
    ggplot2::ylab("S/LRatio z-transformed")
  y_upperbound =
    ifelse(ylim[2] <=
             median(per_bin_profile$size_based_zscore,na.rm = TRUE),
           median(per_bin_profile$size_based_zscore,na.rm = TRUE) + 1,
           ylim[2])
  sl_plot =
    sl_plot + ggplot2::scale_y_continuous(
      limits=c(ylim[1],y_upperbound))

  regionsOffTheChart <- per_bin_profile[
    per_bin_profile$size_based_zscore > y_upperbound,]
  sl_plot = sl_plot +
    ggplot2::geom_point(data=regionsOffTheChart,
                        ggplot2::aes(scaledPos,
                                     y_upperbound),
                        shape=2, size = 1)
  sl_plot <- sl_plot +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.title=ggplot2::element_text(size=14),
                   axis.text.x=ggplot2::element_text(size=14),
                   axis.text.y=ggplot2::element_text(size=14))
  return(sl_plot)
}
