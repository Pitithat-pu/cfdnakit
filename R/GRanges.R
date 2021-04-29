GRanges <- function(seqnames, ranges, gc, mappability) {
  GenomicRanges::GRanges(seqnames = seqnames,
                         ranges = ranges,
                         gc=gc,
                         mappability=mappability)
}
