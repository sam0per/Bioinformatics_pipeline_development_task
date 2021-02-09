rm(list = ls())

pkgs <- c("ggplot2", "optparse", "plotly", "purrr", "dplyr", "htmlwidgets", "MASS")
# Install CRAN packages (if not already installed)
inst <- pkgs %in% installed.packages()
if(length(pkgs[!inst]) > 0) install.packages(pkgs[!inst])
# Load packages into session
invisible(lapply(pkgs, require, character.only=TRUE))

################################################################################################################
##### INPUT ####################################################################################################
option_list = list(
  make_option(c("-x", "--overall"), type = "character", default = NULL,
              help = "Overall read coverage per base data (output DeepTools plotCoverage).", metavar = "character"),
  make_option(c("-t", "--target"), type = "character", default = NULL,
              help = "In- and off-target read coverage per base data (output DeepTools plotCoverage).", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Filename of the output figure without extension.", metavar = "character"))

opt_parser = OptionParser(usage = paste("Rscript scripts/chr19_1_plot_coverage.R",
                                        "-x results/1_read_coverage/chr19_coverage.txt",
                                        "-t results/1_read_coverage/chr19_coverage_target.txt",
                                        "-o figures/1_read_coverage/chr19_r_coverage"),
                          option_list=option_list,
                          description = "Report two figures that tell you how many bases are covered how many times.")
opt = parse_args(opt_parser)

if (is.null(opt$overall) | is.null(opt$target) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("All the arguments must be supplied.\n", call.=FALSE)
}

# Read coverage per base data (output DeepTools plotCoverage)
deepT <- read.table(opt$overall)
deepB <- read.table(opt$target)

# Prepare data for figures
deepT <- merge(x = deepT, y = deepB, by = c("V1", "V2", "V3"))
colnames(deepT) <- c("chr", "start", "end", "xxx", "intarget", "outarget")
# Count number of bases with the same depth
deepC <- apply(X = deepT[, -1:-3], MARGIN = 2, FUN = function(y) {
  aggregate(x = deepT$chr, by = list(y), length)
})

# Merge list of dataframes
ad <- deepC %>% reduce(left_join, by = "Group.1")
colnames(ad) <- c("depth", "xxx", "intarget", "outarget") 
# Replace NAs with 0
ad$outarget <- ifelse(test = is.na(ad$outarget), yes = 0, no = ad$outarget)

# Calculate "reverse cumulative sum"
calc_rev <- function(dat, coln, dep_val) {
  orow <- dat[dat[, 1] %in% dep_val, ]
  for (i in 1:nrow(orow)) {
    if (as.numeric(orow[i, 1])==0) {
      orow[i, paste(coln, "perc", sep = "_")] <- (orow[i, coln]/sum(dat[, coln])) * 100
    } else {
      orow[i, paste(coln, "perc", sep = "_")] <- (sum(dat[dat[,1]>=as.numeric(orow[i, 1]), coln])/sum(dat[, coln]))*100
    }
  }
  return(orow[, paste(coln, "perc", sep = "_")])
}
ad$xxx_perc <- calc_rev(dat = ad, coln = "xxx", dep_val = ad$depth)
ad$intarget_perc <- calc_rev(dat = ad, coln = "intarget", dep_val = ad$depth)
ad$outarget_perc <- calc_rev(dat = ad, coln = "outarget", dep_val = ad$depth)

# Select values of read coverage for the figure
brks <- c(0:10, 15, seq(20,100,length.out = 9), 125, seq(150, 350,by = 50))
red_bin <- ad[ad$depth %in% brks, 1:4]

# Calculate fraction of bases covered
red_bin$xxx_prop <- red_bin$xxx/sum(red_bin$xxx)
red_bin$in_prop <- red_bin$intarget/sum(red_bin$intarget)
red_bin$out_prop <- red_bin$outarget/sum(red_bin$outarget)

dtn <- data.frame(depth=rep(red_bin$depth,3),
                  prop=c(red_bin$xxx_prop, red_bin$in_prop, red_bin$out_prop),
                  pal=c(rep("xxx", nrow(red_bin)),
                        rep("intarget", nrow(red_bin)),
                        rep("outarget", nrow(red_bin))))
# Plot fraction of bases against read coverage without the zero coverage
p_deep <- ggplot(data = dtn[dtn$depth!=0, ], col = "black") +
  geom_histogram(aes(x = as.factor(depth), y = prop, fill = pal), stat = "identity",
                 position = "dodge") +
  theme_bw() +
  scale_fill_manual(values = c("#ff7f0e", "#2ca02c", "#1f77b4")) +
  labs(x = "Coverage (# reads per bp)", y = "Fraction of bases", fill = "") +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 9),
        legend.position = "none") +
  annotate("text", y = 0.038, x = 25, label="Overall mean = 13X", hjust=1, vjust=1, size = 6, col = "#1f77b4") +
  annotate("text", y = 0.031, x = 25, label="Target mean = 11X", hjust=1, vjust=1, size = 6, col = "#ff7f0e") +
  annotate("text", y = 0.024, x = 25, label="Off-target mean = 2X", hjust=1, vjust=1, size = 6,
           col = "#2ca02c")

ggsave(filename = paste0(opt$output, ".png"), plot = p_deep, width = 10, height = 7, dpi = "screen")

# Plot breadth against depth of coverage
fig <- plot_ly(ad[ad$depth<=200,], x = ~depth)
fig <- fig %>% add_trace(y = ~xxx_perc, name = 'Overall', type = 'scatter', mode = 'lines')
fig <- fig %>% add_trace(y = ~intarget_perc, name = 'Target', type = 'scatter', mode = 'lines')
fig <- fig %>% add_trace(y = ~outarget_perc, name = 'Off-target', type = 'scatter', mode = 'lines')
fig <- fig %>% layout(xaxis = list(title = "Coverage (# reads per bp)", titlefont = list(size = 18)),
                      yaxis = list(title = "Percentage of bases >= coverage", titlefont = list(size = 18)),
                      legend = list(font = list(size = 16)))

saveWidget(widget = fig, file = paste0(basename(opt$output), ".html"))
