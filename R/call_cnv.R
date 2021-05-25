#' Call Copy-number Variation from SLRatio and segmentation
#'
#' @param sample_segmentation segmentation dataframe from segmentByPSCBS
#' @param sample_zscore zscore dataframe
#' @param callChr Chromosome to analysis : Default c(1:22)
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom IRanges mergeByOverlaps
call_cnv <- function(sample_segmentation, sample_zscore, callChr=c(1:22)){
  sample_segmentation = dplyr::filter(sample_segmentation, chromosome %in% callChr)
  sample_zscore = dplyr::filter(sample_zscore, chrom %in% callChr)

  segment_grange = GenomicRanges::GRanges(seqnames=sample_segmentation$chromosome,
                           IRanges(start = sample_segmentation$start,
                                   end = sample_segmentation$end),
                           segment_mean= sample_segmentation$mean)

  SLRatio_zscore_grange = GenomicRanges::GRanges(seqnames = sample_zscore$chrom,
                                  IRanges(start = sample_zscore$start,
                                          end = sample_zscore$end),
                                  SLRatio=sample_zscore$SLRatio.corrected)
  overlapping_seg =
    mergeByOverlaps(segment_grange, SLRatio_zscore_grange)
  overlapping_seg = as.data.frame(overlapping_seg)


  segmentation_SLRatio = overlapping_seg %>%
    dplyr::group_by(segment_grange.seqnames,
             segment_grange.start,segment_grange.end,
             segment_grange.width,
             segment_mean) %>%
    dplyr::summarise(SLRatio_median=median(SLRatio,na.rm=TRUE))

  SLref = median(segmentation_SLRatio$SLRatio_median,na.rm=TRUE)
  ### CNV calling parameters
  tfs = seq(0.00,0.8,by = 0.01)
  ploidies = seq(1.5,4,by=0.05)
  CVs=seq(0,5)
  ####


  SL_distance = list()
  for (ploidy in ploidies) {
    for (pred_tf in tfs) {
      SLexp = SLref*(pred_tf*CVs + 2*(1-pred_tf))/(pred_tf*ploidy + 2*(1-pred_tf))
      distance_vec = apply(segmentation_SLRatio, 1, function(segment){
        distances = abs(SLexp - as.numeric(segment["SLRatio_median"]))
        # distances = abs(SLexp - as.numeric(segment["segment_SLRatio"]))

        ### Penalties
        # penalty = 0.2
        # distances_corrected = distances / (pred_tf^penalty)
        # min_distance = min(distances_corrected)
        # predicted_CNV = CVs[which(distances_corrected==min_distance)]
        ###
        min_distance = min(distances)
        predicted_CNV = CVs[which(distances==min_distance)]


        if(length(predicted_CNV)>1) {
          predicted_CNV = round(ploidy)
        }
        # predicted_CNV = mutate(expected_SL_df,
        #                       distance = abs(SLexp - as.numeric(segment["SLRatio_median"]))) %>%
        # arrange(distance) %>% slice(1)
        predicted_CNV_df = data.frame(
          Chromosome = as.character(segment["segment_grange.seqnames"]),
          Start = as.numeric(segment["segment_grange.start"]),
          End = as.numeric(segment["segment_grange.end"]),
          Length=as.numeric(segment["segment_grange.width"]),
          predicted_CNV = predicted_CNV,
          ploidy=ploidy,
          tf=pred_tf,
          distance = min_distance,
          SLexp = SLexp[which(distances==min_distance)],stringsAsFactors = FALSE)
        # SLexp = SLexp[which(distances_corrected==min_distance)],stringsAsFactors = FALSE)
        return(predicted_CNV_df)
      })
      # print(distance_vec)
      SL_distance[[paste(sep="_",ploidy,pred_tf)]] = distance_vec

    }
    # break
  }
  SL_distance_df = do.call(rbind,do.call(rbind, SL_distance))
  SL_distance_df =
    dplyr::mutate(SL_distance_df,
           dist_length=as.numeric(distance)*as.numeric(Length))

  ploidy_tf_mat =
    dplyr::group_by(SL_distance_df,
                    ploidy,tf) %>%
    dplyr::summarise(average_distance=
                sum(dist_length)/sum(as.numeric(Length))) %>%
    dplyr::arrange(average_distance)
  ploidy_tf_mat =
    dplyr::mutate(ploidy_tf_mat,
           rounded_ploidy = round(ploidy))

  summarised_solution =
    dplyr::group_by(ploidy_tf_mat, rounded_ploidy) %>%
    dplyr::summarise(distance=min(average_distance)) %>%
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


get_segment_bysolution <- function(solution,
                                   sample_segmentation,
                                   SL_distance_df){
  temp_solution_df =
    dplyr::filter(SL_distance_df, tf==solution$tf,
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
