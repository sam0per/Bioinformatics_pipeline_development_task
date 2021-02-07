rm(list = ls())

pkgs <- c("ggplot2", "optparse", "patchwork")
# Install CRAN packages (if not already installed)
inst <- pkgs %in% installed.packages()
if(length(pkgs[!inst]) > 0) install.packages(pkgs[!inst])
# Load packages into session
invisible(lapply(pkgs, require, character.only=TRUE))

################################################################################################################
##### INPUT ####################################################################################################
option_list = list(
  make_option(c("-r", "--redundancy"), type = "character", default = NULL,
              help = "Coverage redundancy file.", metavar = "character"),
  make_option(c("-b", "--base-coverage"), type = "character", default = NULL,
              help = "Read coverage per base data (output DeepTools plotCoverage).", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Filename of the output figure.", metavar = "character"))

opt_parser = OptionParser(usage = paste("Rscript scripts/chr19_1_plot_coverage.R",
                                        "-r results/1_read_coverage/chr19_cov_redundancy.txt",
                                        "-b results/1_read_coverage/chr19_coverage.txt"),
                          option_list=option_list,
                          description = "Report a histogram that tells you how many bases are covered how many times.")
opt = parse_args(opt_parser)

if (is.null(opt$redundancy) | is.null(opt$`base-coverage`) | is.null(opt$output)){
  print_help(opt_parser)
  stop("All the arguments must be supplied.\n", call.=FALSE)
}

# Read coverage redundancy table
redu <- read.table(opt$redundancy)
# Reorder
redu <- redu[order(redu$V2), ]

# Fraction of bases covered
redu$prop <- redu$V1/sum(redu$V1)

# Select values of read coverage for the figure
brks <- c(0:10, 15, seq(20,100,length.out = 9), 125, seq(150, 350,by = 50))
red_bin <- redu[redu$V2 %in% brks, ]

# Plot fraction of bases against read coverage without the zero coverage
p_red <- ggplot(data = red_bin[-1,], aes(x = as.factor(V2), y = prop)) +
  geom_histogram(stat = "identity") +
  theme_bw() +
  labs(x = "Coverage (# reads per bp)", y = "Fraction of bases") +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 11)) +
  annotate("text", y = 0.038, x = 25, label="Mean = 14X", hjust=1, vjust=1, size = 6)

# Read coverage per base data (output DeepTools plotCoverage)
deepT <- read.table(opt$`base-coverage`)

# Plot a "reverse cumulative sum"
p_deep <- ggplot(deepT, aes(x = V4, y = 1 - ..y..)) +
  stat_ecdf(size = 1, pad = TRUE) +
  geom_segment(aes(x = 40, y = -Inf, xend = 40, yend = 0.1), linetype = "dashed") +
  geom_segment(aes(x = -Inf, y = 0.1, xend = 40, yend = 0.1), linetype = "dashed") +
  theme_bw() +
  labs(x = "Coverage (# reads per bp)", y = "Fraction of bases >= coverage") +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 11)) +
  scale_x_continuous(breaks = c(seq(0, 80, by = 20), 200, 400, 600),
                     labels = c(seq(0, 80, by = 20), 200, 400, 600)) +
  scale_y_continuous(breaks = c(0.1, seq(0, 1, by = 0.2)),
                     labels = c(0.1, seq(0, 1, by = 0.2))) +
  annotate("text", y = 0.2, x = 300, label="10% of the sampled bp have \nat least 40 overlapping reads",
           hjust=0.5, vjust=1, size = 6)

ggsave(filename = opt$output, plot = p_red + p_deep, width = 20, height = 7, dpi = "screen")
