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
              help = "In-target read coverage per base data (output DeepTools plotCoverage).",
              metavar = "character"),
  make_option(c("-f", "--oftarget"), type = "character", default = NULL,
              help = "Off-target read coverage per base data (output DeepTools plotCoverage).",
              metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Filename of the output figure.", metavar = "character"))

opt_parser = OptionParser(usage = paste("Rscript scripts/chr19_1_plot_coverage.R",
                                        "-x results/1_read_coverage/deep_coverage_chr19.txt",
                                        "-t results/1_read_coverage/deep_chr19_11_coverage.txt",
                                        "-f results/1_read_coverage/deep_chr19_00_coverage.txt",
                                        "-o figures/1_read_coverage/chr19_r_coverage.png"),
                          option_list=option_list,
                          description = "Report one figure that tells you how many bases are covered how many times.")
opt = parse_args(opt_parser)

if (is.null(opt$overall) | is.null(opt$target) | is.null(opt$oftarget) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("All the arguments must be supplied.\n", call.=FALSE)
}

# Read coverage per base data (output DeepTools plotCoverage)
# and sub-sample randomly without replacement
deepT <- read.table(opt$overall)
deepT$exp <- "overall"
idt <- deepT[sample(x = 1:nrow(deepT), size = 10000, replace = FALSE), ]

deepB <- read.table(opt$target)
deepB$exp <- "target"
idb <- deepB[sample(x = 1:nrow(deepB), size = 10000, replace = FALSE), ]

deepO <- read.table(opt$oftarget)
deepO$exp <- "oftarget"
ido <- deepO[sample(x = 1:nrow(deepO), size = 10000, replace = FALSE), ]

rm(list = c("deepB", "deepO", "deepT"))

idt <- idt[order(idt$V2),]
# mean(idt$V4)
idb <- idb[order(idb$V2),]
# mean(idb$V4)
ido <- ido[order(ido$V2),]
# mean(ido$V4)

# Plot depth against chr position
p_deep <- ggplot(data = idb) +
  geom_point(aes(x = V2, y = V4, col = "target"), alpha = 0.5) +
  geom_point(data = ido, aes(x = V2, y = V4, col = "off-target"), alpha = 0.5) +
  geom_point(data = idt, aes(x = V2, y = V4, col = "overall"), alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = c("#2ca02c", "#1f77b4", "#ff7f0e")) +
  labs(x = "portion of chromosome 19", y = "depth of coverage", col = "") +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14),
        legend.position = "top", legend.text = element_text(size = 16)) +
  annotate("text", y = 510, x = 2550000, label="off-target mean = 29X", hjust=0.5, vjust=1, size = 6, col = "#2ca02c") +
  annotate("text", y = 560, x = 7500000, label="overall mean = 13X", hjust=0.5, vjust=1, size = 6, col = "#1f77b4") +
  annotate("text", y = 510, x = 12550000, label="in-target mean = 144X", hjust=0.5, vjust=1, size = 6,
           col = "#ff7f0e") +
  guides(col = guide_legend(override.aes = list(size = 5)))
# p_deep

# Overdispersion and uniformity of coverage
# Negative binomial model
# Overall data
m1 <- glm.nb(V4 ~ 1, data = idt)
cat("\nOverall depth of coverage:")
uci <- round(1/(m1$theta - 2*m1$SE.theta), 3)
lci <- round(1/(m1$theta + 2*m1$SE.theta), 3)
cat("Dispersion parameter = ", round(1/m1$theta, 3), " [", lci, "-", uci, "]", "\n", sep = "")

# Target data
m2 <- glm.nb(V4 ~ 1, data = idb)
cat("\nTarget depth of coverage:")
uci <- round(1/(m2$theta - 2*m2$SE.theta), 3)
lci <- round(1/(m2$theta + 2*m2$SE.theta), 3)
cat("Dispersion parameter = ", round(1/m2$theta, 3), " [", lci, "-", uci, "]", "\n", sep = "")

# Off-target data
m3 <- glm.nb(V4 ~ 1, data = ido)
cat("\nOff-target depth of coverage:")
uci <- round(1/(m3$theta - 2*m3$SE.theta), 3)
lci <- round(1/(m3$theta + 2*m3$SE.theta), 3)
cat("Dispersion parameter = ", round(1/m3$theta, 3), " [", lci, "-", uci, "]", "\n", sep = "")

# Add dispersion results to figure
p_deep <- p_deep +
  geom_text(aes(x = 2550000, y = 510-30, label=paste0("dispersion = ", round(1/m3$theta, 3))),
            hjust=0.5, vjust=1, size = 5, col = "#2ca02c", data = data.frame()) +
  geom_text(aes(y = 560-30, x = 7500000), label=paste0("dispersion = ", round(1/m1$theta, 3)),
            hjust=0.5, vjust=1, size = 5, col = "#1f77b4", data = data.frame()) +
  geom_text(aes(y = 510-30, x = 12550000), label=paste0("dispersion = ", round(1/m2$theta, 3)),
            hjust=0.5, vjust=1, size = 5, col = "#ff7f0e", data = data.frame())

# Save figure
ggsave(filename = paste0(opt$output), plot = p_deep, width = 10, height = 7, dpi = "screen")

cat("...MISSION COMPLETE...\n")
