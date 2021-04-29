#' Getting fragment-length information
#' Function extract insert-size (fragment-length)
#'
#' @param readbam_bin list
#' @param short_range numeric vector
#' @param long_range numeric vector
#' @param maximum_length numeric
#' @param minimum_length numeric
#'
#' @return list
#' @export
#'
#' @examples
#'
get_fragment_profile <- function(readbam_bin,
                                 short_range = c(100,150),
                                 long_range = c(151,250),
                                 maximum_length = 600,
                                 minimum_length = 20){
  binsize = imply_binsize(readbam_bin)/1000
  sliding_windows_gr = util.get_sliding_windows(binsize)
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
                  "S/L.Ratio.corrected" =
                    short.corrected/long.corrected,
                  GC = sliding_windows_gr$gc,
                  mappability = sliding_windows_gr$mappability)

  # isize_vector = as.numeric(unlist(
  #   lapply(readbam_bin, function(bin) {
  #     abs(bin$isize)[which(dplyr::between(abs(bin$isize),
  #                                  left = minimum_length,
  #                                  right = maximum_length))]
  #   })))
  isize_vector = extract_insert_size(readbam_bin,
                                     maximum_length,minimum_length)
  insert_info_df =
    data.frame("Total Fragments" = length(isize_vector),
               "Read Pairs in range" = sum(bin_profile_df$nfragment),
               "Read Pairs in range_corrected" = sum(bin_profile_df$total.corrected,na.rm=T),
               "Mode" = getmode(isize_vector),
               "Median" = median(isize_vector, na.rm = TRUE),
               "Mean" = round(mean(isize_vector, na.rm = TRUE), 2),
               "SD" = round(sd(isize_vector, na.rm = TRUE),2),
               "short" = sum(bin_profile_df$short,na.rm=T),
               "long" = sum(bin_profile_df$long,na.rm=T),
               "short_corrected" = sum(bin_profile_df$short.corrected,na.rm=T),
               "long_corrected" = sum(bin_profile_df$long.corrected,na.rm=T),
               "S/L Ratio" = round(sum(bin_profile_df$short, na.rm=T) /
                                     sum(bin_profile_df$long, na.rm=T),2),
               "S/L Ratio_corrected" = round(sum(bin_profile_df$short.corrected, na.rm=T) /
                 sum(bin_profile_df$long.corrected, na.rm=T),2),
               "Bin Size(KB)"=binsize)

  density_table = make_density_table(isize_vector)

  fragment_profile = list("per_bin_profile" = bin_profile_df,
                          "sample_profile" = insert_info_df,
                          "distribution_table" = density_table,
                          "insert_size" = isize_vector)
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

make_density_table <- function(isize_vector){
  if(length(na.omit(isize_vector)) < 1e6){
    d = density(na.omit(isize_vector))
    density_df = d[c("x","y")]
  } else {
    fraction_table = table(na.omit(isize_vector))
    fraction_vec = fraction_table / sum(fraction_table)
    density_df = data.frame("x"=as.numeric(names(fraction_vec)),
                            "y"=as.numeric(unname(fraction_vec)))
  }
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
