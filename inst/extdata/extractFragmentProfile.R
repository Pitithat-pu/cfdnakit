suppressPackageStartupMessages(library(cfdnakit))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-b", "--bamfile"),
              action = "store", type = "character",
              default = NULL, help = "Path to bamfile"),
  make_option(c("-i", "--sampleid"),
              action = "store", type = "character",
              default = NULL, help = "Sample id"),
  make_option(c("--outdir"), action = "store", type = "character",
              default = NULL,
              help = "Output directory: Path to output directory, Created if doesn't exist."),
  make_option(c("--binsize"), action = "store", type = "integer",
              default = 1000,
              help = "Size of non-overlapping bin in Kilobase. default=1000 Available options: 100, 500, 1000"),
  make_option(c("--fragmentinfo"), action = "store", type = "logical",
              default = TRUE,
              help = "TRUE : Print general fragment-length information into tab-separated file. default=TRUE"),
  make_option(c("--plot_dist"), action = "store", type = "logical",
              default = TRUE,
              help = "TRUE : Plot fragment-length distribution. default=TRUE"),
  make_option(c("--plot_SLRatio"), action = "store", type = "logical",
              default = TRUE,
              help = "TRUE : Plot genome-wide SLRatio . default=TRUE")
  )

opt <- parse_args(OptionParser(option_list = option_list))
# >>>>> Checking input arguments
if(! file.exists(opt$bamfile)){
  stop("Input bamfile doesn't exist.")
}
if(! dir.exists(opt$outdir)){

  print(paste0("Output directory ", basename(opt$outdir) ,
        " doesn't exist: Create output directory"))
  print(opt$outdir)
  dir.create(opt$outdir)
} else {
  print(paste0("Output directory ",opt$outdir))
}
if(! opt$binsize %in% c(100,500,1000)){
  stop(paste0("Suggested binsize argument ",
              opt$binsize ," is not possible."))
}
# <<<<<<<
print(paste0("Reading Bamfile and split into ",opt$binsize, "KB"))
sample_bambin = read_bamfile(opt$bamfile,binsize = opt$binsize)
print("Extracting fragment length profile")
sample_profile = get_fragment_profile(sample_bambin,
                                      sample_id = opt$sampleid)

print("Writing sample profile rds")
save_fragment_profile(sample_profile,
                      output_dir = opt$outdir,
                      overwrite = TRUE)

if(opt$fragmentinfo){
  print("Writing fragment info file")
  write.table(sample_profile$sample_profile,
              file = paste0(opt$outdir,"/",
                            sample_profile$Sample.ID,".fragmentinfo"),
              quote = FALSE,
              sep = "\t",row.names = FALSE,
              col.names = TRUE)
}


if(opt$plot_dist){
  print("Plotting fragment-length distribution")
  print(paste0(opt$outdir,"/",
               sample_profile$Sample.ID,".fragmentdist.png"))
  png(filename = paste0(opt$outdir,"/",
                    sample_profile$Sample.ID,".fragmentdist.png"),
      width = 700,height = 500)
  print(plot_isize_distribution(sample_profile,
                          sample_name = sample_profile$Sample.ID))
  dev.off()
}


if(opt$plot_SLRatio){
  print("Plotting genome-wide SLRatio")
  print(paste0(opt$outdir,"/",
               sample_profile$Sample.ID,".SLRatio.png"))
  png(filename = paste0(opt$outdir,"/",
                    sample_profile$Sample.ID,".SLRatio.png"),
      width = 1000,height = 350)
  print(plot_sl_ratio(fragment_profile = sample_profile))
  dev.off()
}
print("Done")
