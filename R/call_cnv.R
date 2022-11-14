#' Call Copy-number Variation from SLRatio and segmentation
#'
#' @param sample_segmentation segmentation dataframe from segmentByPSCBS
#' @param sample_zscore zscore dataframe
#' @param callChr chromosome to analysis : Default c(1:22)
#' @param tfs range of fitting tumor fraction : Default c(0,0.8)
#' @param ploidies range of fitting chromosomal ploidy : Default c(1.5,4)
#' @param MaxCN maximum copy-number : Default 4
#'
#' @return List of cnvcalling solutions
#' @export
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
#'
#' sample_zscore_segment = segmentByPSCB(sample_zscore)
#'
#' sample_cnv = call_cnv(sample_zscore_segment,sample_zscore, tfs=c(0.1,0.3),ploidies=c(1.5,2), MaxCN=3)
#' plot_cnv_solution(sample_cnv,selected_solution = 1)
#'
#' @importFrom IRanges mergeByOverlaps
#' @importFrom stats setNames
#' @importFrom rlang .data
call_cnv <- function(sample_segmentation, sample_zscore,
                     callChr=seq_len(22),
                     tfs = c(0,0.7),
                     ploidies = c(1.5,3),
                     MaxCN=4){

  tfs = seq(tfs[1],tfs[2],by = 0.01)
  ploidies = seq(ploidies[1],ploidies[2],by=0.05)
  CVs=seq(0,MaxCN)

  sample_segmentation =
    dplyr::filter(sample_segmentation,
                  .data$chromosome %in% callChr)
  sample_zscore = dplyr::filter(sample_zscore,
                                .data$chrom %in% callChr)

  segment_grange = GenomicRanges::GRanges(
    seqnames=sample_segmentation$chromosome,
    IRanges(start = sample_segmentation$start,
            end = sample_segmentation$end),
    segment_mean= sample_segmentation$mean)

  SLRatio_zscore_grange = GenomicRanges::GRanges(
    seqnames = sample_zscore$chrom,
    IRanges(start = sample_zscore$start,
            end = sample_zscore$end),
    SLRatio=sample_zscore$SLRatio.corrected)
  overlapping_seg =
    mergeByOverlaps(segment_grange, SLRatio_zscore_grange)
  overlapping_seg = as.data.frame(overlapping_seg)


  segmentation_SLRatio = overlapping_seg %>%
    dplyr::group_by(.data$segment_grange.seqnames,
                    .data$segment_grange.start,
                    .data$segment_grange.end,
                    .data$segment_grange.width,
                    .data$segment_mean) %>%
    dplyr::summarise(SLRatio_median=
                       median(.data$SLRatio,na.rm=TRUE))

  SLref = median(segmentation_SLRatio$SLRatio_median,na.rm=TRUE)



  SL_distance_tf = list()
  for (pred_tf in tfs) {

    ploidies.dominators =
      (pred_tf*ploidies + 2*(1-pred_tf))
    SL.expected.mat = sapply(ploidies.dominators, function(ploidies.dominator){
      SL.expected =
        SLref*(pred_tf*CVs + 2*(1-pred_tf))/ploidies.dominator
    })


    distances.list = apply(SL.expected.mat, 2, function(SL.expecteds){
      distances.lst =
        lapply(SL.expecteds, function(x){
          abs(x - segmentation_SLRatio$SLRatio_median)
        })
    })

    distance_mat = sapply(distances.list, function(x){
      distance.mat = do.call(rbind,x)
      min_distances_order = apply(distance.mat, 2, order)
      distance_vec = apply(distance.mat, 2, min)
    })

    predicted_CNV_mat = sapply(distances.list, function(x){
      distance.mat = do.call(rbind,x)
      min_distances_order = apply(distance.mat, 2, order)
      return(CVs[min_distances_order[1,]])
    })
    predicted_CNV_df = data.frame(predicted_CNV_mat)
    colnames(predicted_CNV_df) = ploidies

    SL_distance.df = apply(distance_mat, 2, function(distance_vec){
      predicted_CNV_df = data.frame(
        Chromosome = as.character(
          segmentation_SLRatio$segment_grange.seqnames),
        Start = as.numeric(
          segmentation_SLRatio$segment_grange.start),
        End = as.numeric(
          segmentation_SLRatio$segment_grange.end),
        Length=as.numeric(
          segmentation_SLRatio$segment_grange.width),
        tf=pred_tf,
        distance = distance_vec,
        stringsAsFactors = FALSE)
      # SLexp = SLexp[which(distances_corrected==min_distance)],stringsAsFactors = FALSE)
      return(predicted_CNV_df)
    })

    SL_distance.df = do.call(rbind,SL_distance.df)
    SL_distance.df$predicted_CNV =
      unlist(predicted_CNV_df[,c(as.character(ploidies))],
             use.names = FALSE)
    SL_distance.df$ploidy =
      sort(rep(ploidies,nrow(sample_segmentation)))

    SL_distance_tf[[as.character(pred_tf)]] = SL_distance.df
  }

  SL_distance_df = do.call(rbind, SL_distance_tf)
  SL_distance_df =
    dplyr::mutate(SL_distance_df,
                  dist_length=
                    as.numeric(.data$distance)*
                    as.numeric(.data$Length))

  ploidy_tf_mat =
    dplyr::group_by(SL_distance_df,
                    .data$ploidy,.data$tf) %>%
    dplyr::summarise(average_distance=
                       sum(.data$dist_length)/
                       sum(as.numeric(.data$Length))) %>%
    dplyr::arrange(.data$average_distance)
  ploidy_tf_mat =
    dplyr::mutate(ploidy_tf_mat,
                  rounded_ploidy = round(.data$ploidy))

  summarised_solution =
    dplyr::group_by(ploidy_tf_mat, .data$rounded_ploidy) %>%
    dplyr::summarise(distance=min(.data$average_distance)) %>%
    dplyr::inner_join(ploidy_tf_mat,
                      by = c("rounded_ploidy"="rounded_ploidy",
                             "distance"="average_distance"))
  summarised_solution$rank = rank(summarised_solution$distance)
  summarised_solution = summarised_solution %>%
    dplyr::arrange(rank)

  solution_list =
    setNames(split(summarised_solution,
                   seq(nrow(summarised_solution))),
             summarised_solution$rank)
  segment_by_solution_list =
    lapply(solution_list,get_segment_bysolution,
           sample_segmentation = sample_segmentation,
           SL_distance_df = SL_distance_df)

  segment_by_solution_list[["per_bin_profile"]] = sample_zscore
  segment_by_solution_list[["ploidy_tf_mat"]] = ploidy_tf_mat
  return(segment_by_solution_list)

}


#' Return CNV segmentation result from given all CNV solutions
#'
#' @param solution solution dataframe
#' @param sample_segmentation Segmeantion dataframe
#' @param SL_distance_df Distance matrix
#'
#' @return list of segmentation per solution
#'
#' @importFrom PSCBS ploidy
get_segment_bysolution <- function(solution,
                                   sample_segmentation,
                                   SL_distance_df){
  temp_solution_df =
    dplyr::filter(SL_distance_df, .data$tf==solution$tf,
             ploidy==solution$ploidy)
  solution_segment_df =
    dplyr::inner_join(sample_segmentation,
                      temp_solution_df,
                      by=c("chromosome"="Chromosome",
                           "start"="Start","end"="End")) %>%
    dplyr::select(-c("tf","ploidy"))
  solution_vec = c("TF"=solution$tf,
                   "rounded_ploidy"=solution$rounded_ploidy,
                   "ploidy"=solution$ploidy,
                   "distance"=solution$distance,
                   "rank"=solution$rank)
  solution_list = list(cnv_solution = solution_vec,
                       solution_segmentation = solution_segment_df)
  return(solution_list)
}
