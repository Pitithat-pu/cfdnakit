#' Plot Fragment-length Distribution
#'
#' @param fragment_profile list
#' @param sample_name charactor
#' @param xlim numeric
#'
#' @return distribution plot
#' @export
#'
#' @examples
#' @importFrom ggplot2 ggplot aes geom_line geom_text geom_vline
plot_isize_distribution <- function(fragment_profile,
                                    sample_name = NULL,
                                    xlim = 500){
  if (is.null(fragment_profile[["distribution_table"]])){
    stop("Element distribution_table doen't exist in fragment_profile")
  }


  density_df = fragment_profile$distribution_table
  if (!is.null(sample_name)) {
    density_df$Label = sample_name
  } else density_df$Label = fragment_profile$Sample.ID

  control_fragment_profile = util.load_control_density_table()
  control_density_df = control_fragment_profile$control_density_df
  temp_density_df = rbind(control_density_df,density_df)

  #>>>>> extract KS test result comparing sample and control fragment distribution
  # ks_test = test_KolmogorovSmirnov(control_fragment_profile,
  #                                  fragment_profile)
  ks_test = dplyr::select(fragment_profile$sample_profile,
                          "K.S.p.value","K.S.stats")
  ks_test_text = paste0("KS test: stats ",
                        ks_test$K.S.stats)
  # ks_test_text = paste0("KS test: p-value ",
  #                       ks_test$K.S.p.value,"; stats ",
  #                       ks_test$K.S.stats)
  #<<<<<

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
                         unique(temp_density_df$Label)],
                       "hadjust"=c(1.5,-.4))
  #<<<<<

  density_plot = ggplot2::ggplot(temp_density_df,
                        aes(x,y))
  density_plot = density_plot +
    geom_line(aes(color=Label))
  density_plot = density_plot +
    geom_vline(data=peak_df,
               aes(xintercept=peak_length,
                   color=Label),linetype=2,
               show.legend = FALSE)
  density_plot = density_plot +
    geom_text(data=peak_df,
              aes(peak_length,0,label=peak_length,
                  color=Label,hjust=hadjust),vjust=-0.25,
              size=5,show.legend = FALSE)
  density_plot = density_plot +
    geom_text(x=Inf,y=Inf,hjust=1.1,vjust=1.4,
              label=ks_test_text,size=5, check_overlap = TRUE)
  density_plot = density_plot +
    ggplot2::scale_x_continuous(limits = c(0,xlim),
                       breaks = seq(0,xlim,50), expand = c(0, 0),
                       name = "Fragment Length (bases)")
  density_plot = density_plot +
    ggplot2::scale_y_continuous(breaks = seq(0,max(temp_density_df$y),0.005),
                       expand = c(0, 0),
                       name = "Proportion")
  density_plot = density_plot +
    ggplot2::theme_linedraw()
  density_plot = density_plot +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "top",
                   legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size=14),
                   axis.title = ggplot2::element_text(size=14),
                   axis.text = ggplot2::element_text(size=12))
  return(density_plot)
}
