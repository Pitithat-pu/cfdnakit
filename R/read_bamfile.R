
#' Read a bam file
#' Read a bam file from give path. Alignment and sequencing read information
#'  will be binned into non-overlapping size
#' @param bamfile_path Character; Path to sample bamfile
#' @param binsize Int; Size of non-overlapping windows in KB. Only 100,500 and 1000 is available; Default 1000
#' @param blacklist_files Character; Filepath to file containing blacklist regions
#' @param genome Character; abbreviation of reference genome; available genome: hg19,hg38, mm10. default:hg19
#' @param target_bedfile Character; Path to exon/target bedfile; Default NULL
#' @param min_mapq Int; minimum read mapping quality; Default 20
#' @param apply_blacklist Logical; To exclude read on the blacklist regions Default TRUE
#'
#' @return SampleBam Object; A list object containing read information from the BAM file.
#' @export
#'
#' @examples
read_bamfile <- function(bamfile_path, binsize=1000, blacklist_files=NULL ,
                         genome="hg19" ,target_bedfile=NULL,
                         min_mapq=20, apply_blacklist= TRUE){
  if(!utils.file_exists(bamfile_path)){
    stop("The bamfile doesn't exist. Please check if the path to bamfile is valid.")
  }
  if(!is.null(target_bedfile) & !utils.file_exists(target_bedfile)) {
    stop("The given target bedfile doesn't exist.")
  }


  which <- util.get_sliding_windows(binsize = binsize, genome=genome)
  if(if_ucsc_chrformat(bamfile_path)){
    which <- GRCh2UCSCGRanges(which)
  }

  flag <- Rsamtools::scanBamFlag(isPaired = TRUE,
                                 isUnmappedQuery = FALSE,
                                 isDuplicate = FALSE,
                                 isMinusStrand = FALSE,
                                 hasUnmappedMate = FALSE,
                                 isSecondaryAlignment = FALSE,
                                 isMateMinusStrand = TRUE)

  param <- Rsamtools::ScanBamParam(what = c("qname", "rname", "pos",
                                            "isize", "qwidth"),
                                   flag = flag,
                                   which = which,
                                   mapqFilter = min_mapq)
  if(!if_exist_baifile(bamfile = bamfile_path)){
    print("The BAM index file (.bai) is missing. Creating an index file")
    Rsamtools::indexBam(bamfile_path)
    print("Bam index file created.")
  }
  print("Reading bamfile")
  bam <- Rsamtools::scanBam(file = bamfile_path,
                            index = bamfile_path,
                            param=param)


  if (apply_blacklist) {
    print("Filtering-out read on the blacklist regions")
    bam = filter_read_on_blacklist(bam, blacklist_files , genome=genome)
  }
  if(!is.null(target_bedfile)){
    print("Extracting on-target fragments")
    bam = extract_read_ontarget(bam,target_bedfile)
  }
  class(bam)="SampleBam"

  if(if_ucsc_chrformat(bamfile_path)){
    bam <- UCSC2GRChSampleBam(bam)
  }
  return(bam)
}

#' Check if bai file exist from given bam
#'
#' @param bamfile Character; Path to sample bamfile
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

#' Check UCSC chromosomes format for input bam file
#'
#' @param bamfile_path Character; Path to sample bamfile
#'
#' @return Boolean; if the input bam file is UCSC format, chr prefix
#'
#' @examples
if_ucsc_chrformat <- function(bamfile_path){
  all_chrs = Rsamtools::idxstatsBam(bamfile_path)$seqnames
  return(any(all_chrs %in% c("chr1")))
}


#' Convert GRCh chromosome format to UCSC style
#'
#' @param which GRanges object;
#'
#' @return GRanges; GRanges after chromosome format conversion
#'
#' @examples
GRCh2UCSCGRanges<- function (which) {
  GenomeInfoDb::seqlevels(which)<-
    sub('chrM',
        'M',GenomeInfoDb::seqlevels(which))
  GenomeInfoDb::seqlevels(which)<-
    gsub('^(.*)$','chr\\1',GenomeInfoDb::seqlevels(which))
  return(which)
}

#' Convert UCSC chromosome format to GRCh style from a list of alignment information
#'
#' @param sample.bam list of alignment information from function read_bamfile
#'
#' @return List; list of alignment information after conversion
#'
#' @examples
UCSC2GRChSampleBam <- function (sample.bam) {
  names(sample.bam) = gsub(
    "^chr","", names(sample.bam)
  )

  return(sample.bam)
}
#' Filter out reads on blacklist regions
#'
#' @param sample_bin SampleBam; Object from function read_bamfile
#' @param blacklist_files Character; Filepath to file containing blacklist regions
#' @param genome Character; Abbreviation of reference genome; Either hg19 or mm10. default:hg19
#'
#' @return SampleBam after filtering out read on balck list regions
#'
#' @examples
filter_read_on_blacklist <- function(sample_bin, blacklist_files=NULL , genome="hg19"){
    ### Blacklist of hg19 is available as default without giving files
    if(is.null(blacklist_files) & genome=="hg19"){
      blacklist_files =
        c(system.file("extdata","hg19_centromere.tsv.gz",
                      package = "cfdnakit"),
          system.file("extdata", "wgEncodeDacMapabilityConsensusExcludable.bed_GRCh37.gz",
                      package = "cfdnakit"))
    } else if(is.null(blacklist_files) & genome!="hg19") {
      ### When no blacklist files is not provided, return the whole input.
      message("Blacklist files were not given.")
      return(sample_bin)
    }

  blacklist_targets_gr <- create_blacklist_gr(blacklist_files)
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


    filtered_gr <- suppressWarnings(GenomicRanges::findOverlaps(bin_gr,
                                               blacklist_targets_gr))
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
#' @param blacklist_files Character; Filepath to file containing blacklist regions
#' @return GRanges object of blacklist regions
#'
#' @examples
create_blacklist_gr <- function(blacklist_files){
  if(!utils.file_exists(blacklist_files))
    stop("One of blacklist file doen't exist. Please check if the blacklist file exist in extdata directory.")
  blacklist_targets_gr = GenomicRanges::GRanges()
  for (blacklist_region in blacklist_files) {
    if(if_gzfile(blacklist_region)){
      con = gzfile(blacklist_region)
      blacklist_targets =
        as.data.frame(utils::read.table(con,
                                        header = FALSE,sep = "\t",
                                        stringsAsFactors = FALSE))[,1:3]
      # close(con)
    } else {
      blacklist_targets =
        as.data.frame(utils::read.table(blacklist_region),
                      header = FALSE,sep = "\t",
                      stringsAsFactors = FALSE)[,1:3]
    }

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

