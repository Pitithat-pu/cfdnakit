#' Plot Fragment-length Distribution
#'
#' @param fragment_profile list
#' @param xlim numeric
#'
#' @return
#' @export
#'
#' @examples
plot_isize_distribution <- function(fragment_profile,
                                    xlim = 500){
  if (is.null(fragment_profile[["distribution_table"]])){
    stop("Element distribution_table doen't exist in fragment_profile")
  }
  density_df = fragment_profile$distribution_table
  density_df$Label = "Patient cfDNA"
  control_fragment_profile = load_control_density_table()
  control_density_df = control_fragment_profile$control_density_df
  # isize_vec = control_fragment_profile$insert_size
  ks_test = test_KolmogorovSmirnov(control_fragment_profile,
                                   fragment_profile)
  temp_density_df = rbind(control_density_df,density_df)
  peak_length <-
    sapply(split(temp_density_df, temp_density_df$Label),
           function(density_df){
             peak_length = density_df$x[
               which(density_df$y == max(density_df$y))]
             peak_length = floor(peak_length)})
  peak_df = data.frame("Label" = unique(temp_density_df$Label),
                       "peak_length" = peak_length[
                         unique(temp_density_df$Label)],
                       "hadjust"=c(1.5,-.4))

  density_plot = ggplot(temp_density_df,
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
    scale_x_continuous(limits = c(0,xlim),
                       breaks = seq(0,xlim,50), expand = c(0, 0),
                       name = "Fragment Length (bases)")
  density_plot = density_plot +
    scale_y_continuous(breaks = seq(0,max(temp_density_df$y),0.005),
                       expand = c(0, 0),
                       name = "Proportion")
  density_plot = density_plot +
    theme_linedraw()
  density_plot = density_plot +
    theme(panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size=14),
          axis.title = element_text(size=14),
          axis.text = element_text(size=12))
  return(density_plot)
}

load_control_density_table <- function(){
  # control_density_file =
  #   system.file("extdata",
  #               "BH01_rmdup_paired_mapped_Fragment-length_report_50_insert_size_density.csv",
  #               package = "cfdnakit")
  # control_density_df = read.table(control_density_file, header = TRUE)
  # colnames(control_density_df)=c("x","y")
  # control_density_df$Label = "Healthy Control"
  # return(control_density_df)
  control_RDS_file =
    system.file("extdata",
                "BH01_rmdup_paired_mapped_0.01.fragmentprofile.RDS",
                package = "cfdnakit")
  control_fragment_profile =
    readRDS(control_RDS_file)
  colnames(control_fragment_profile$distribution_table)=c("x","y")
  control_fragment_profile$distribution_table$Label = "Healthy Control"
  return(list("control_density_df"=
         control_fragment_profile$distribution_table,
       "insert_size" =
         control_fragment_profile$insert_size))
}

test_KolmogorovSmirnov <- function(
  control_fragment_profile, sample_fragment_profile){
  ks_result = ks.test(control_fragment_profile$insert_size,
          sample_fragment_profile$insert_size)
  ks_result = signif(ks_result$p.value,3)
  return(ks_result)
}
#
# test_ecdf = ecdf(density_df$y)
# plot(test_ecdf)
