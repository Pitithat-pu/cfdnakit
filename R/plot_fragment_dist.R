#' Plot Fragment-length Distribution
#'
#' @param readbam_list List; A list containing SampleBam object/objects from the read_bamfile function
#' @param maximum_length Int; Maximum length of fragment. cfDNA fragment longer than this value will not be considered; Default 550
#' @param minimum_length Int; Minimum length of fragment. cfDNA fragment shorter than this value will not be considered;  Default 20
#'
#' @return distribution plot
#' @export
#'
#' @examples
#' @importFrom ggplot2 ggplot aes geom_line geom_text geom_vline
#' @importFrom dplyr mutate select
plot_fragment_dist <- function(readbam_list,
                               maximum_length = 550,
                               minimum_length = 20){

  density_lst = lapply(seq_along(readbam_list), function(i) {
    sample.name = names(readbam_list)[[i]]
    readbam_bin=readbam_list[[i]]
    isize_vector = extract_insert_size(readbam_bin,
                                       maximum_length,minimum_length)
    density_df = make_density_table(isize_vector,
                                    minimum_length,maximum_length)
    density_df$Label = sample.name
    return(density_df)
  })


  temp_density_df = do.call(rbind,density_lst)
  #>>>> Find the peak fragment-length
  peak_length <-
    sapply(split(temp_density_df, temp_density_df$Label),
           function(density_df){
             peak_length = density_df$x[
               which(density_df$y == max(density_df$y))]
             if(length(peak_length) > 1) #### if given more than 1 peak, return only one that close to 167
               peak_length = max(peak_length)[which.min(abs(max(peak_length - 167)))]
             peak_length = floor(peak_length)})

  peak_df = data.frame("Label" = unique(temp_density_df$Label),
                       "peak_length" = peak_length[
                         unique(temp_density_df$Label)]) %>%
    dplyr::mutate("Legend.Label" = paste0(Label," (",peak_length,")"))


  temp_density_df = temp_density_df %>%
    dplyr::left_join(peak_df %>% dplyr::select(Label,peak_length,Legend.Label),
              by = "Label")
  #<<<<<

  density_plot = ggplot2::ggplot(temp_density_df,
                                 aes(x,y))
  density_plot = density_plot +
    geom_line(aes(color=Legend.Label))
  density_plot = density_plot +
    ggplot2::scale_color_discrete(name="Peak Length")
  density_plot = density_plot +
    geom_vline(data=peak_df,
               aes(xintercept=peak_length,
                   color=Legend.Label),linetype=2,
               show.legend = FALSE)

  density_plot = density_plot +
    ggplot2::scale_x_continuous(limits = c(0,maximum_length),
                                breaks = seq(0,maximum_length,50), expand = c(0, 0),
                                name = "Fragment Length (bases)")
  density_plot = density_plot +
    ggplot2::scale_y_continuous(breaks = seq(0,max(temp_density_df$y),0.005),
                                expand = c(0, 0),
                                name = "Proportion")
  density_plot = density_plot +
    ggplot2::theme_linedraw()
  density_plot = density_plot +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = c(1, 1),
                   legend.justification = c(1, 1),
                   legend.title = ggplot2::element_text(size=14),
                   legend.text = ggplot2::element_text(size=14),
                   axis.title = ggplot2::element_text(size=14),
                   axis.text = ggplot2::element_text(size=12))
  return(density_plot)
}
