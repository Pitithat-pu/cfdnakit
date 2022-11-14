#' Plot Fragment-length profile with CNV calling result
#'
#' @param cnvcall solution results from call_cnv function
#' @param selected_solution solution rank to plot
#' @param genome Character; version of reference genome (default hg19)
#' @param ylim Vector of 2 Int; ylim of plot (default c(-20,20))
#'
#' @return ggplot object plot Genomics CNV profile of selected solution
#' @export
#' @importFrom rlang .data
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
#' plot_cnv_solution(sample_cnv,selected_solution = 1)
#'
plot_cnv_solution <- function(cnvcall,
                              selected_solution = 1,
                              genome="hg19",
                              ylim=c(-30,30)){
  cnv_col = c("AMP"="#00FF00","GAIN"="#008200",
              "NEUT"="grey30","HETD"="#820000",
              "HOMD"="#FF0000","Not available"="grey")

  if(! genome %in% c("hg19","hg38"))
    stop("Only hg19 or hg38 genome are possible")

  chrTotalLength_file= paste0(genome,"_chrTotalLength.tsv")
  chrLength_file =
    system.file("extdata",
                chrTotalLength_file,
                package = "cfdnakit")
  chrLength_df = read.table(file = chrLength_file,
                            header=FALSE, sep="\t")
  chrLength_info = util.get_chrLength_info(chrLength_df)


  per_bin_profile = cnvcall[["per_bin_profile"]]
  per_bin_profile =
    dplyr::mutate(per_bin_profile,
                  scaledPos = (.data$start + .data$end)/2 +
                    chrLength_info$chroffsets[.data$chrom])
  sample_segment_df =
    cnvcall[[selected_solution]]$solution_segmentation
  sample_segment_df =
    dplyr::mutate(sample_segment_df,
                  scaledPos_start = .data$start +
                    chrLength_info$chroffsets[.data$chromosome],
                  scaledPos_end= .data$end +
                    chrLength_info$chroffsets[.data$chromosome])
  sample_segment_df = sample_segment_df %>%
    dplyr::mutate(CNV.group = dplyr::case_when(
      predicted_CNV == 2 ~ "NEUT",
      predicted_CNV >= 4 ~ "AMP",
      predicted_CNV == 3 ~ "GAIN",
      predicted_CNV == 1 ~ "HETD",
      predicted_CNV == 0 ~ "HOMD",
      TRUE ~ "Not available"
    ))
  per_bin_profile =
    overlap_bin_with_segment(per_bin_profile,sample_segment_df)

  ##### Actually Create plot
  cnv_plot = ggplot2::ggplot(per_bin_profile,
                            ggplot2::aes_(~scaledPos,~size_based_zscore))
  cnv_plot = cnv_plot +
    ggplot2::geom_point( ggplot2::aes_(color=~CNV.group), size = 1)
  cnv_plot = cnv_plot +
    ggplot2::scale_colour_manual(values = cnv_col)
  cnv_plot = cnv_plot +
    ggplot2::geom_hline(yintercept = 0 ,
                        size=0.5)
  cnv_plot = cnv_plot +
    ggplot2::geom_segment(data=sample_segment_df,
                          ggplot2::aes_(y=~mean, yend=~mean,
                              x=~scaledPos_start,xend=~scaledPos_end),
                          size = 1.5, color="darkorange")
  cnv_plot = cnv_plot +
    ggplot2::scale_x_continuous(breaks = chrLength_info$chrMids,
                                labels = chrLength_info$chrNames,
                                position = "bottom",expand=c(0,0))
  cnv_plot = cnv_plot +
    ggplot2::geom_vline(xintercept = chrLength_info$chroffsets,
                        linetype="dotted",size=0.5)
  cnv_plot = cnv_plot +
    ggplot2::xlab("Chromosome number") +
    ggplot2::ylab("S/LRatio z-transformed")
  y_upperbound =
    ifelse(ylim[2] <=
             median(per_bin_profile$size_based_zscore,na.rm = TRUE),
           median(per_bin_profile$size_based_zscore,na.rm = TRUE) + 1,
           ylim[2])
  cnv_plot = cnv_plot +
    ggplot2::scale_y_continuous(
      limits=c(ylim[1],y_upperbound))

  regionsOffTheChart <- per_bin_profile[
    per_bin_profile$size_based_zscore > y_upperbound,]
  cnv_plot = cnv_plot +
    ggplot2::geom_point(data=regionsOffTheChart,
                        ggplot2::aes_(~scaledPos,
                                     ~y_upperbound),
                        shape=2, size = 1)
  cnv_plot = cnv_plot +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.title=ggplot2::element_text(size=14),
                   axis.text.x=ggplot2::element_text(size=14),
                   axis.text.y=ggplot2::element_text(size=14),
                   legend.position = "none")
  return(cnv_plot)
}


#' Overlap and merge bin data frame with segmentation dataframe
#'
#' @param per_bin_profile bin dataframe
#' @param sample_segmentation segmentation dataframe
#'
#' @return dataframe of overlapping bin and segmentation
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps pintersect
#' @importFrom IRanges width
#' @importFrom S4Vectors subjectHits queryHits
overlap_bin_with_segment <-
  function(per_bin_profile,sample_segmentation){
    segment_grange =
      makeGRangesFromDataFrame(sample_segmentation)

    per_bin_profile_grange =
      makeGRangesFromDataFrame(per_bin_profile)

    #### Only 18:78000001-78077248 has no hit (possible end to chromosome)
    hits = findOverlaps(segment_grange, per_bin_profile_grange)

    ### find overlapping segment
    ### If more than 1 hit, keep only highest percent overlapping
    overlapping_seg =
      cbind(per_bin_profile[subjectHits(hits),],
            sample_segmentation[queryHits(hits),
                                c("predicted_CNV","CNV.group")])
    overlaps = pintersect(segment_grange[queryHits(hits)],
                           per_bin_profile_grange[subjectHits(hits)])
    percentOverlap = width(overlaps) / width(per_bin_profile_grange[subjectHits(hits)])
    overlapping_seg$p.overlap = percentOverlap
    overlapping_seg = overlapping_seg[order(overlapping_seg$chrom,
                                overlapping_seg$start,
                                overlapping_seg$end,
                                -overlapping_seg$p.overlap ), ]
    overlapping_seg = overlapping_seg[ !duplicated(overlapping_seg[c("chrom","start","end")]), ]

    overlapping_seg = overlapping_seg %>%
      dplyr::select(-c("p.overlap"))
    rownames(overlapping_seg) = NULL
    return(overlapping_seg)

}
