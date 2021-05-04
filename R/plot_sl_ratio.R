
#' Plot Short/Long-fragment Ratio
#'
#' @param fragment_profile list
#' @param chrTotalLength_file character
#' @param ylim numeric vector
#' @return plot
#' @export
#'
#' @examples
#' @importFrom utils read.table
plot_sl_ratio <- function(fragment_profile,
                          chrTotalLength_file=
                            "hg19_chrTotalLength.tsv",
                          ylim=c(0,0.4)){
  chrLength_file =
    system.file("extdata",
                chrTotalLength_file,
                package = "cfdnakit")
  chrLength_df = read.table(file = chrLength_file,header=F, sep="\t")
  chrLength_info = get_chrLength_info(chrLength_df)
  per_bin_profile = rowname_to_columns(fragment_profile$per_bin_profile)
  per_bin_profile =
    dplyr::mutate(per_bin_profile,
                  scaledPos = (start + end)/2 +
                    chrLength_info$chroffsets[chrom])

  sl_plot = ggplot2::ggplot(per_bin_profile,
                            ggplot2::aes(scaledPos,`S/L.Ratio.corrected`))
  sl_plot = sl_plot +
    ggplot2::geom_point(fill="grey20", size = 1)
  sl_plot = sl_plot +
    ggplot2::geom_hline(yintercept =
                          median(per_bin_profile$`S/L.Ratio.corrected`,
                                 na.rm = TRUE),
                        size=0.5)
  sl_plot = sl_plot +
    ggplot2::scale_x_continuous(breaks = chrLength_info$chrMids,
                                labels = chrLength_info$chrNames,
                                position = "bottom",expand=c(0,0))
  sl_plot = sl_plot +
    ggplot2::geom_vline(xintercept = chrLength_info$chroffsets,
                        linetype="dotted",size=0.5)
  sl_plot = sl_plot + ggplot2::xlab("Chromosome number") +
    ggplot2::ylab("Short/Long-Fragment Ratio")
  y_upperbound = median(per_bin_profile$`S/L.Ratio.corrected`, na.rm = TRUE) +
    mad(per_bin_profile$`S/L.Ratio.corrected`, na.rm = TRUE)
  sl_plot =
    sl_plot + ggplot2::scale_y_continuous(
      limits=c(ylim[1],
                   ifelse(ylim[2] <= y_upperbound,
                          max(per_bin_profile$`S/L.Ratio.corrected`,na.rm = TRUE),
                          ylim[2])))

  regionsOffTheChart <- per_bin_profile[
    per_bin_profile$`S/L.Ratio.corrected` > ylim[2],]
  sl_plot = sl_plot +
    ggplot2::geom_point(data=regionsOffTheChart,
                        ggplot2::aes(scaledPos,
                                     ylim[2]),
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

rowname_to_columns <- function(per_bin_profile){
  splited_list=lapply(strsplit(
    rownames(per_bin_profile),
    split = "[:-]"), function(x) {
      data.frame("chrom"=x[1],
                 "start"=as.numeric(x[2]),
                 "end"=as.numeric(x[3]))
    })
  cbind(do.call(rbind,splited_list),per_bin_profile)
}

get_chrLength_info <- function(chrLength_df){
  colnames(chrLength_df) = c("Chromosome", "Length")
  chrLength_df = chrLength_df[which(chrLength_df$Chromosome!="Y"),]
  chrNames = chrLength_df$Chromosome
  chrLength = chrLength_df$Length
  names(chrLength) = chrNames
  chroffsets = cumsum(as.numeric(chrLength))
  chroffsets <- c(0, chroffsets[0:(length(chroffsets)-1)])
  names(chroffsets) <- names(chrLength)
  chrMids <- cumsum(as.numeric(chrLength))
  chrMids <- (chrMids + chroffsets)/2
  names(chrMids) <- names(chrLength)

  chrLength_info =
    list("chrNames"=chrNames,
         "chrLength"=chrLength,
         "chroffsets"=chroffsets,
         "chrMids"=chrMids)

}
