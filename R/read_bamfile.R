
#' Read a bam file
#' Read a bam file from give path. Alignment and sequencing read information
#'  will be binned into non-overlapping size
#' @param bamfile_path Path to sample bamfile
#' @param binsize Size of non-overlapping windows in KB. Max = 1000; Default 1000
#' @param target_bedfile Path to exon/target bedfile; Default NULL
#'
#' @return List of fragment information per bin
#' @export
#'
#' @examples
read_bamfile <- function(bamfile_path, binsize=1000, target_bedfile=NULL){
  if(!utils.file_exists(bamfile_path)){
    stop("Bamfile doesn't exist. Please check if the path to bamfile is valid.")
  }
  if(!is.null(target_bedfile) & !utils.file_exists(target_bedfile)) {
    stop("Given target bedfile doesn't exist.")
  }
  which <- util.get_sliding_windows(binsize = binsize)
  flag <- Rsamtools::scanBamFlag(isPaired = TRUE,
                                 isUnmappedQuery = FALSE,
                                 isDuplicate = FALSE,
                                 isMinusStrand = FALSE,
                                 hasUnmappedMate = FALSE,
                                 isSecondaryAlignment = FALSE,
                                 isMateMinusStrand = TRUE,
                                 isNotPassingQualityControls = FALSE,
                                 isSupplementaryAlignment = FALSE)

  param <- Rsamtools::ScanBamParam(what = c("qname", "rname", "pos",
                                            "isize", "qwidth"),
                                   flag = flag,
                                   which = which,
                                   mapqFilter = 20)
  if(!if_exist_baifile(bamfile = bamfile_path)){
    print("BAM index is missing. Creating an index file")
    Rsamtools::indexBam(bamfile_path)
  }
  print("Reading bamfile")
  bam <- Rsamtools::scanBam(file = bamfile_path,
                            index = bamfile_path,
                            param=param)
    # bam <- Rsamtools::scanBam(file = bamfile_path,
    #                           param=param)
  print("Filtering-out read on blacklist")
  bam_blacklist_lst = filter_read_on_blacklist(bam)
  if(!is.null(target_bedfile)){
    print("Extracting on-target fragments")
    bam_blacklist_lst = extract_read_ontarget(bam_blacklist_lst,target_bedfile)
  }
  return(bam_blacklist_lst)
}

#' Check if bai file exist from given bam
#'
#' @param bamfile Path to bamfile
#'
#' @return Boolean if the bai file exist
#'
#' @examples
if_exist_baifile <- function(bamfile) {
  baifile_name = paste0(bamfile,".bai")
  utils.file_exists(baifile_name)
}


extract_read_ontarget <- function(sample_bin, target_bedfile){
  if(!utils.file_exists(target_bedfile)) stop("Given target bedfile doesn't exist.")
  if(if_gzfile(target_bedfile)){
    target_df <- as.data.frame(read.table(
      gzfile(target_bedfile),
      header = FALSE,stringsAsFactors = FALSE,
      sep = "\t"))[,1:3]
  } else {
    target_df <- as.data.frame(read.table(
      target_bedfile,
      header = FALSE,stringsAsFactors = FALSE,
      sep = "\t"))[,1:3]
  }

  target_gr=GenomicRanges::GRanges(
    seqnames = target_df[,1],
    ranges = IRanges::IRanges(start = as.numeric(target_df[,2]),
                              end=as.numeric(target_df[,3])))

  ontarget_bam_lst = lapply(sample_bin, function(region_lst){
    bin_gr =
      GenomicRanges::GRanges(seqnames =
                               as.character(region_lst$rname),
                             ranges = IRanges::IRanges(start = region_lst$pos,
                                                       end = region_lst$pos +
                                                         region_lst$qwidth),
                             qname=region_lst$qname,
                             rname=region_lst$rname,
                             pos=region_lst$pos,
                             qwidth=region_lst$qwidth,
                             isize = region_lst$isize)


    ontarget_gr <- GenomicRanges::findOverlaps(bin_gr,
                                               target_gr)
    if(length(ontarget_gr@from) == 0 )
      ontarget_bam_gr = bin_gr
    else{
      ontarget_bam_gr = bin_gr[ontarget_gr@from]
    }
    return_vec = list("qname" = ontarget_bam_gr$qname,
                      "rname" = ontarget_bam_gr$rname,
                      "pos" = ontarget_bam_gr$pos,
                      "qwidth" = ontarget_bam_gr$qwidth,
                      "isize" = ontarget_bam_gr$isize)
  })
  return(ontarget_bam_lst)
}


if_gzfile <- function(bedfile){
  if(summary( file(bedfile) )$class == "gzfile")
    TRUE
  else FALSE
}


#' Filter out reads on blacklist regions
#'
#' @param sample_bin Variable sample bin from read_bamfile function
#'
#' @return Sample read per bin after filtering out read on balck list regions
#'
#' @examples
filter_read_on_blacklist <- function(sample_bin){
  blacklist_targets_gr <- create_blacklist_gr()
  filtered_bam_lst = lapply(sample_bin, function(region_lst){
    bin_gr =
      GenomicRanges::GRanges(seqnames =
                               as.character(region_lst$rname),
                             ranges = IRanges::IRanges(start = region_lst$pos,
                                              end = region_lst$pos +
                                                region_lst$qwidth),
                             qname=region_lst$qname,
                             rname=region_lst$rname,
                             pos=region_lst$pos,
                             qwidth=region_lst$qwidth,
                             isize = region_lst$isize)


    filtered_gr <- GenomicRanges::findOverlaps(bin_gr,
                                               blacklist_targets_gr)
    if(length(filtered_gr@from) == 0 )
      filterd_bam_gr = bin_gr
    else{
      filterd_bam_gr = bin_gr[-filtered_gr@from]
    }
      return_vec = list("qname" = filterd_bam_gr$qname,
                      "rname" = filterd_bam_gr$rname,
                      "pos" = filterd_bam_gr$pos,
                      "qwidth" = filterd_bam_gr$qwidth,
                      "isize" = filterd_bam_gr$isize)
  })
  return(filtered_bam_lst)
}

#' Create Blacklist regions GRanges object
#'
#' @return GRanges object of blacklist regions
#'
#' @examples
create_blacklist_gr <- function(){
  centromere_region = system.file("extdata",
              "hg19_centromere.tsv.gz",
              package = "cfdnakit")
  dac_blacklist_region = system.file("extdata",
              "wgEncodeDacMapabilityConsensusExcludable.bed_GRCh37.gz",
              package = "cfdnakit")
  duke_blacklist_region = system.file("extdata",
              "wgEncodeDukeMapabilityRegionsExcludable.bed_GRCh37.gz",
              package = "cfdnakit")
  blacklist_regions = c(centromere_region,
                        dac_blacklist_region,
                        duke_blacklist_region)
  if(!utils.file_exists(blacklist_regions))
    stop("One of blacklist file doen't exist. Please check if the blacklist file exist in extdata directory.")
  blacklist_targets_gr = GenomicRanges::GRanges()
  for (blacklist_region in blacklist_regions) {
    blacklist_targets =
      as.data.frame(utils::read.table(gzfile(blacklist_region),
                               header = FALSE,sep = "\t",
                               stringsAsFactors = FALSE))[,1:3]
    colnames(blacklist_targets) = c("chromosome","start","end")
    blacklist_targets_gr_temp=
      GenomicRanges::GRanges(seqnames = blacklist_targets$chromosome,
              ranges = IRanges::IRanges(
                start = as.numeric(blacklist_targets$start),
                end=as.numeric(blacklist_targets$end)))

    blacklist_targets_gr = c(blacklist_targets_gr,blacklist_targets_gr_temp)
  }
  return(blacklist_targets_gr)
}

