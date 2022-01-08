#' Getting fragment-length information
#'
#' @param readbam_bin SampleBam Object
#' @param sample_id Character; Given sample ID
#' @param genome abbreviation of reference genome; namely hg19, mm10. default:hg19
#' @param short_range Vector of 2 Int; Range of fragment length to be defined as short fragment; Default c(100,150)
#' @param long_range Vector of 2 Int; Range of fragment length to be defined as long fragment; Default c(151,250)
#' @param maximum_length Int; Maximum length of fragment. cfDNA fragment longer than this value will not be considered; Default 600
#' @param minimum_length Int; Minimum length of fragment. cfDNA fragment shorter than this value will not be considered;  Default 20
#'
#' @return SampleFragment Object; Fragment length information for quality check and downstream analysis per bin and summary of sample
#' @export
#'
#' @examples
#' @importFrom stats sd
get_fragment_profile <- function(readbam_bin,
                                 sample_id,
                                 genome="hg19",
                                 short_range = c(100,150),
                                 long_range = c(151,250),
                                 maximum_length = 600,
                                 minimum_length = 20){
  if (sample_id=="" || is.na(sample_id)) {
    stop("Please specify sample id param (sample_id).")
  }
  binsize = imply_binsize(readbam_bin)/1000
  sliding_windows_gr = util.get_sliding_windows(binsize = binsize, genome=genome)
  bin_profile_df = as.data.frame(t(sapply(readbam_bin, function(bin){
    isize = bin$isize
    nfragment =
      length(isize[which(isize >= short_range[1] & isize <= long_range[2])])
    short =
      length(isize[which(isize >= short_range[1] & isize <= short_range[2])])
    long =
      length(isize[which(isize >= long_range[1] & isize <= long_range[2])])
    c("nfragment" = nfragment,
      "short" = short,
      "long" = long)
  })))
  # View(bin_profile_df)
  bin_profile_df =
    dplyr::mutate(bin_profile_df,
                  "S/L.Ratio" = short/long ,
                  total.corrected =
                    util.bias_correct(nfragment,sliding_windows_gr$gc/100),
                  short.corrected =
                    util.bias_correct(short,sliding_windows_gr$gc/100),
                  long.corrected =
                    util.bias_correct(long,sliding_windows_gr$gc/100))
  bin_profile_df =
    dplyr::mutate(bin_profile_df,
                  total.corrected =
                    util.bias_correct(total.corrected,
                                      sliding_windows_gr$mappability/100),
                  short.corrected =
                    util.bias_correct(short.corrected,
                                      sliding_windows_gr$mappability/100),
                  long.corrected =
                    util.bias_correct(long.corrected,
                                      sliding_windows_gr$mappability/100))
  bin_profile_df =
    dplyr::mutate(bin_profile_df,
                  "S/L.Ratio.corrected" =
                    short.corrected/long.corrected,
                  GC = sliding_windows_gr$gc,
                  mappability = sliding_windows_gr$mappability)
  bin_profile_df = dplyr::filter(bin_profile_df,
                                 nfragment>0)
  isize_vector = extract_insert_size(readbam_bin,
                                     maximum_length,minimum_length)

  #>>>>> Do KS test comparing sample and control fragment distribution
  # control_fragment_profile = util.load_control_density_table()
  # ks_test = test_isize_KolmogorovSmirnov(
  #   control_fragment_profile$insert_size, isize_vector)
  #<<<<<

  insert_info_df =
    data.frame("Sample.ID"=sample_id,
               "Total Fragments" = length(isize_vector),
               "Read Pairs in range" = sum(bin_profile_df$nfragment),
               "Read Pairs in range_corrected" = sum(bin_profile_df$total.corrected,na.rm=T),
               "Mode" = getmode(isize_vector),
               "Median" = median(isize_vector, na.rm = TRUE),
               "Mean" = round(mean(isize_vector, na.rm = TRUE), 2),
               "Mad" = round(mad(isize_vector, na.rm = TRUE),2),
               "short" = sum(bin_profile_df$short,na.rm=T),
               "long" = sum(bin_profile_df$long,na.rm=T),
               "S/L Ratio" = round(sum(bin_profile_df$short, na.rm=T) /
                                     sum(bin_profile_df$long, na.rm=T),2),
               "short_corrected" = sum(bin_profile_df$short.corrected,na.rm=T),
               "long_corrected" = sum(bin_profile_df$long.corrected,na.rm=T),
               "S/L Ratio_corrected" = round(sum(bin_profile_df$short.corrected, na.rm=T) /
                 sum(bin_profile_df$long.corrected, na.rm=T),2),
               # "K.S.p.value" = ks_test$p.value,
               # "K.S.stats" = ks_test$statistic,
               "Bin Size(KB)"=binsize)

  # density_table = make_density_table(isize_vector)

  fragment_profile = list("Sample.ID"=sample_id,
                          "per_bin_profile" = bin_profile_df,
                          "sample_profile" = insert_info_df,
                          # "distribution_table" = density_table,
                          "minimum_length"=minimum_length,
                          "maximum_length"=maximum_length)
  class(fragment_profile) = "SampleFragment"

  return(fragment_profile)
}

imply_binsize <- function(readbam_bin){
  first_binname = names(readbam_bin[1])
  splited_name = unlist(
    strsplit(first_binname, split = "[:-]+"))
  binsize =  as.numeric(splited_name[3]) -
    as.numeric(splited_name[2]) + 1
  return(binsize)
}

getmode <- function(v) {
  uniqv = unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#' Make Fragment-length density table
#'
#' @param isize_vector numeric
#' @param minimum_length numeric
#' @param maximum_length numeric
#'
#' @return data.frame
#'
#' @examples
#' @importFrom stats na.omit density
make_density_table = function(isize_vector,
                              minimum_length,
                              maximum_length){
  if(length(na.omit(isize_vector)) < 100){
    d = density(na.omit(isize_vector))
    density_df = data.frame(d[c("x","y")])
    density_df$x = round(density_df$x,digits = 0)
  } else {
    fraction_table = table(na.omit(isize_vector))
    temp_fraction_vec = fraction_table / sum(fraction_table)

    fraction_vec = rep(0,maximum_length-minimum_length+1)
    names(fraction_vec) = seq(minimum_length,maximum_length)
    fraction_vec[names(fraction_table)] = temp_fraction_vec
    density_df = data.frame("x"=as.numeric(names(fraction_vec)),
                            "y"=as.numeric(unname(fraction_vec)))
  }
  return(density_df)
}

extract_insert_size <- function(readbam_bin,
                                maximum_length = 600,
                                minimum_length = 20) {
  isize_vector = Biobase::subListExtract(readbam_bin, "isize")
  isize_vector = unlist(isize_vector)
  isize_vector = abs(isize_vector)[which(dplyr::between(abs(isize_vector),
                                                 left = minimum_length,
                                                 right = maximum_length))]
  isize_vector = unname(isize_vector)
}




#' KolmogorovSmirnov test for insert size
#'
#' @param control_insert_size Vector of insert size of a control sample
#' @param sample_insert_size Vector of insert size of a testing sample
#'
#' @return KS.Test result
#' @export
#'
#' @examples
test_isize_KolmogorovSmirnov <- function(
  control_insert_size, sample_insert_size){
  ks_result = suppressWarnings(ks.test(na.omit(control_insert_size),
                      na.omit(sample_insert_size)))
  ks_result$p.value = signif(ks_result$p.value,4)
  ks_result$statistic = round(ks_result$statistic,2)
  return(ks_result)
}


#' Get insert-size distribution table
#'
#' @param readbam_bin SampleBam Object from function read_bamfile
#' @param maximum_length Int; Maximum length of fragment. cfDNA fragment longer than this value will not be considered; Default 600
#' @param minimum_length Int; Minimum length of fragment. cfDNA fragment shorter than this value will not be considered;  Default 20
#'
#' @return Distribution table of fragment length
#' @export
#'
#' @examples
get_isize_density <- function(readbam_bin,
                              maximum_length = 600,
                              minimum_length = 20){
  isize_vector = extract_insert_size(readbam_bin,
                                     maximum_length,minimum_length)
  density_table = make_density_table(isize_vector)
}

