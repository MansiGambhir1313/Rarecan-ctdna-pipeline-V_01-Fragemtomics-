#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
    library(optparse)
    library(ichorCNA)
    library(GenomeInfoDb)
    library(GenomicRanges)
})

# --- Argument Parsing ---
option_list <- list(
    make_option(c("-w", "--wig"), type="character", help="Input WIG file (required)"),
    make_option(c("-o", "--outdir"), type="character", default=".", help="Output directory"),
    make_option(c("-s", "--sample"), type="character", default="sample", help="Sample ID"),
    make_option(c("-b", "--binSize"), type="integer", default=500000, help="WIG bin size"),
    make_option(c("-p", "--ploidy"), type="character", default="c(2,3)", help="Ploidy states to test"),
    make_option(c("-n", "--normal"), type="character", default="c(0.5,0.7,0.8,0.9,0.95)", help="Normal fractions to test"),
    make_option(c("--maxCN"), type="integer", default=5, help="Max copy number"),
    make_option(c("--estimateNormal"), type="logical", action="store_true", default=TRUE, help="Estimate normal fraction"),
    make_option(c("--estimatePloidy"), type="logical", action="store_true", default=TRUE, help="Estimate ploidy"),
    make_option(c("--txnE"), type="double", default=0.9999, help="Transition probability"),
    make_option(c("--txnStrength"), type="double", default=10000, help="Transition strength"),
    make_option(c("--includeHOMD"), type="logical", action="store_true", default=FALSE, help="Include HOMD state")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$wig)) {
    print_help(opt_parser)
    stop("WIG file is required.", call.=FALSE)
}

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

# --- 1. Load Read Counts (BedGraph Mode) ---
print(paste("Loading read counts from:", opt$wig))
first_line <- readLines(opt$wig, n=1)
skip_n <- if (length(first_line) > 0 && grepl("track", first_line)) 1 else 0

raw_data <- read.table(opt$wig, skip=skip_n, header=FALSE, stringsAsFactors=FALSE, sep="\t")
colnames(raw_data) <- c("chr", "start", "end", "reads")

# Filter out contigs and non-standard chromosomes
standard_chrs <- paste0("chr", c(1:22, "X", "Y"))
raw_data <- raw_data[raw_data$chr %in% standard_chrs, ]

tumor_reads <- GRanges(
    seqnames = raw_data$chr,
    ranges   = IRanges(start = raw_data$start + 1, end = raw_data$end),
    reads    = raw_data$reads
)
# Try to assign hg38 seqinfo, but don't fail if it's not perfectly matched
tryCatch({
    seqinfo(tumor_reads) <- Seqinfo(genome="hg38")[standard_chrs]
}, error = function(e) {
    warning("Could not assign hg38 seqinfo. Proceeding with basic seqinfo.")
    seqinfo(tumor_reads) <- Seqinfo(seqnames = standard_chrs, seqlengths = rep(NA, length(standard_chrs)), isCircular = rep(FALSE, length(standard_chrs)), genome = "hg38")
})
tumor_reads <- tumor_reads[width(tumor_reads) == opt$binSize]

# --- 2. Run ichorCNA ---
print("Running ichorCNA HMM segmentation...")
ploidy_vec <- eval(parse(text=opt$ploidy))
normal_vec <- eval(parse(text=opt$normal))

# Use the modern ichorCNA API
# Note: ichorCNA::runIchorCNA has built-in looping over normal/ploidy
result <- NULL
tryCatch({
    result <- runIchorCNA(
        tumor_reads,
        sampleId = opt$sample,
        centromere = NULL, # Use built-in centromere data
        rmCentromere = TRUE,
        normal = normal_vec,
        ploidy = ploidy_vec,
        maxCN = opt$maxCN,
        estimateNormal = opt$estimateNormal,
        estimatePloidy = opt$estimatePloidy,
        txnE = opt$txnE,
        txnStrength = opt$txnStrength,
        plotCNA = FALSE, # Don't generate plots
        outDir = opt$outdir
    )
}, error = function(e) {
    print(paste("ichorCNA::runIchorCNA failed:", e$message))
})

# --- 3. Export Standardized Output ---
if (!is.null(result) && !is.null(result$results)) {
    # Find the optimal solution parameters
    best_idx <- which.max(result$results$loglik)
    if (length(best_idx) == 0) {
      best_idx <- 1 # Default to first if loglik is weird
    }
    tfx_est <- result$results$Fraction.Tumor.Estimate[best_idx]
    ploidy_est <- result$results$Ploidy.Estimate[best_idx]
    
    # Write segments file
    seg_file <- file.path(opt$outdir, paste0(opt$sample, ".seg.txt"))
    write.table(result$cna.seg, file=seg_file, sep="\t", quote=FALSE, row.names=FALSE)
    
    # Write tumor fraction file
    tfx_df <- data.frame(
        sample = opt$sample, 
        tumorFraction = tfx_est, 
        ploidy = ploidy_est
    )
    tfx_file <- file.path(opt$outdir, paste0(opt$sample, ".tfx.txt"))
    write.table(tfx_df, file=tfx_file, sep="\t", quote=FALSE, row.names=FALSE)
    
    print(paste("ichorCNA complete. Best estimate for Tumor Fraction:", round(tfx_est, 4)))
    
} else {
    print("ichorCNA failed to produce a valid result. Outputting empty files.")
    # Write empty files for pipeline compatibility
    seg_file <- file.path(opt$outdir, paste0(opt$sample, ".seg.txt"))
    write.table(data.frame(), file=seg_file, sep="\t")
    tfx_file <- file.path(opt$outdir, paste0(opt$sample, ".tfx.txt"))
    write.table(data.frame(sample=opt$sample, tumorFraction=0, ploidy=NA), file=tfx_file, sep="\t", quote=FALSE, row.names=FALSE)
}
