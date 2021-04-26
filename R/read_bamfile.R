
#' Read a bam file
#' Read a bam file from give path. Alignment and sequencing read information
#'  will be binned into non-overlapping size
#' @param bamfile_path character
#' @param binsize numeric
#'
#' @return
#' @export
#'
#' @examples

read_bamfile <- function(bamfile_path, binsize){
  if(!utils.file_exists(bamfile_path)){
    paste0("Bamfile doesn't exist. Please check if the path to bamfile is valid.")
    quit(status=1)
  }
  blacklist_targets_gr <- create_blacklist_gr()
  which <- get_sliding_windows(binsize = binsize)
  which_filtered <-
    filter_blacklist_regions(which,blacklist_targets_gr)

  flag <- Rsamtools::scanBamFlag(isPaired = TRUE,
                   isUnmappedQuery = FALSE,
                   isDuplicate = FALSE,
                   hasUnmappedMate = FALSE,
                   isSecondaryAlignment = FALSE,
                   isSupplementaryAlignment = FALSE)

  param <- Rsamtools::ScanBamParam(what = c("qname", "rname", "pos",
                                            "isize", "qwidth"),
                                   flag = flag,
                                   which = which_filtered,
                                   mapqFilter = 20)

  if(if_exist_baifile(bamfile = bamfile_path)){
    bam <- Rsamtools::scanBam(file = bamfile_path,
                              index = bamfile_path,
                              param=param)
  } else
    bam <- Rsamtools::scanBam(file = bamfile_path,
                              param=param)

}

if_exist_baifile <- function(bamfile) {
  baifile_name = paste(bamfile,".bai")
  utils.file_exists(baifile_name)
}

filter_blacklist_regions <- function(which,blacklist_targets_gr){
  filtered_gr <- GenomicRanges::setdiff(which,blacklist_targets_gr)
  return(filtered_gr)
}

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
    q(save = 'no',status = 1)
  blacklist_targets_gr = GenomicRanges::GRanges()
  for (blacklist_region in blacklist_regions) {
    blacklist_targets =
      as.data.frame(utils::read.table(gzfile(blacklist_region),
                               header = FALSE,sep = "\t",
                               stringsAsFactors = FALSE))[,1:3]
    colnames(blacklist_targets) = c("chromosome","start","end")
    blacklist_targets_gr_temp=
      GRanges(seqnames = blacklist_targets$chromosome,
              ranges = IRanges(
                start = as.numeric(blacklist_targets$start),
                end=as.numeric(blacklist_targets$end)))

    blacklist_targets_gr = c(blacklist_targets_gr,blacklist_targets_gr_temp)
  }
  return(blacklist_targets_gr)
}

get_sliding_windows <- function(binsize){
  qdnaseq_sliding_windows_RDS =
    system.file("extdata",
                paste0("AnnotationDataFrame_from_QDNAseq_",binsize,"k.rds"),
                package = "cfdnakit")
  if (utils.file_exists(qdnaseq_sliding_windows_RDS)) {
    bins = readRDS(qdnaseq_sliding_windows_RDS)
  } else {
    paste0("The selected binsize (",binsize,") is not available.")
    paste0("Available binsize (kb) are 1000, 500, 100.")
    quit("no",status=1)
  }
  sliding_windows <- as.data.frame(bins@data)
  sliding_windows <- sliding_windows[which(sliding_windows$chromosome!="Y" &
                                            sliding_windows$mappability>=1),]
  sliding_windows_gr <- GRanges(seqnames = sliding_windows$chrom,
                                ranges = IRanges(
                                  start= sliding_windows$start,
                                  end = sliding_windows$end))

  return(sliding_windows_gr)
}
