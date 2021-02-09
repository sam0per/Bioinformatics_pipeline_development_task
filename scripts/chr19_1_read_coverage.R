rm(list = ls())
library(ggplot2)

# Total number of bases of the chr19 region
length(60004:14992498)

# Read coverage redundancy table
redu <- read.table("Bioinformatics_pipeline_development_task/results/1_read_coverage/chr19_cov_redundancy.txt")
redu <- redu[order(redu$V2), ]
head(redu)

# Percentage of non-covered bases
(redu[1,1]/length(60004:14992498))*100

# Percentage of bases covered by at least one read
100-((redu[1,1]/length(60004:14992498))*100)
# Percentage of bases covered by at least 5 reads
(sum(redu[redu$V2>=5, 1])/length(60004:14992498))*100

# Plot coverage
redu$prop <- redu$V1/sum(redu$V1)
sum(redu$prop)
head(redu)
max(redu$V2)
brks <- c(0:10, 15, seq(20,100,length.out = 9), 125, seq(150, 350,by = 50))
red_bin <- redu[redu$V2 %in% brks, ]

head(red_bin)
ggplot(data = red_bin[-1,], aes(x = as.factor(V2), y = prop)) +
   geom_histogram(stat = "identity") +
   # scale_x_continuous(breaks=c(1,2,3,4,5,10,30), limits = c(0, 31)) +
   # scale_x_continuous(breaks=c(1,2,3,4,5,10,30,100,300,1000), trans="log1p", expand=c(0,0)) +
   theme_bw() +
   labs(x = "Coverage (# reads per bp)", y = "Fraction of bases") +
   theme(axis.title = element_text(size = 16), axis.text = element_text(size = 11)) +
   annotate("text", y = 0.038, x = 25, label="Mean = 14X", hjust=1, vjust=1, size = 6)
   # scale_x_continuous(breaks = c(1, seq(10, 50, by = 10)), limits = c(0, 51), labels = c(1, seq(10, 50, by = 10)))
   
head(redu)
redu$cums <- cumsum(redu$prop)
redu$icums <- round(c(redu$cums[1], 1-redu$cums[-1]), 4)

ggplot(data = redu) +
   # geom_line(aes(x = V2, y = cums)) +
   geom_line(aes(x = V2, y = icums))

x <- rnorm(12)
Fn <- ecdf(x)
Fn     # a *function*
Fn(x)  # returns the percentiles for x

deepT <- read.table("Bioinformatics_pipeline_development_task/results/1_read_coverage/chr19_coverage.txt")
head(deepT)

deepB <- read.table("Bioinformatics_pipeline_development_task/results/1_read_coverage/chr19_coverage_target.txt")
head(deepB)
identical(deepT$V2, deepB$V2)
identical(deepT$V3, deepB$V3)

deepT <- merge(x = deepT, y = deepB, by = c("V1", "V2", "V3"))
head(deepT)
colnames(deepT) <- c("chr", "start", "end", "xxx", "intarget", "outarget")
deepC <- apply(X = deepT[, -1:-3], MARGIN = 2, FUN = function(y) {
   aggregate(x = deepT$chr, by = list(y), length)
})

identical(deepC$xxx$Group.1,deepC$intarget$Group.1)
identical(deepC$xxx$Group.1,deepC$outarget$Group.1)

library(purrr)
library(dplyr)
ad <- deepC %>% reduce(left_join, by = "Group.1")
colnames(ad) <- c("depth", "oveall", "intarget", "outarget") 
head(ad)
ad$outarget <- ifelse(test = is.na(ad$outarget), yes = 0, no = ad$outarget)
summary(ad)

# Percentage of non-covered bases
(deepC[1,2]/sum(deepC$x))*100

# Percentage of bases covered by at least 10 reads
(sum(deepC[deepC$Group.1>=10, 2])/sum(deepC$x))*100
# Percentage of bases covered by at least 40 reads
(sum(deepC[deepC$Group.1>=40, 2])/sum(deepC$x))*100

# Plot a "reverse cumulative sum"
ggplot(deepT, aes(x = V4, y = 1 - ..y..)) +
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

library(plotly)
library(patchwork)
head(ad)
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
head(ad)

fig <- plot_ly(ad[ad$depth<=200,], x = ~depth)
fig <- fig %>% add_trace(y = ~xxx_perc, name = 'overall', type = 'scatter', mode = 'lines')
fig <- fig %>% add_trace(y = ~intarget_perc, name = 'intarget', type = 'scatter', mode = 'lines')
fig <- fig %>% add_trace(y = ~outarget_perc, name = 'outarget', type = 'scatter', mode = 'lines')
# pfont <- list(font = "Courier New, monospace", size = 18, color = "black")
fig <- fig %>% layout(xaxis = list(title = "Coverage (# reads per bp)", titlefont = list(size = 18)),
                      yaxis = list(title = "Fraction of bases >= coverage", titlefont = list(size = 18)),
                      legend = list(font = list(size = 16)))
fig
head(ad)
ad$xxx_prop <- ad$xxx/sum(ad$xxx)
sum(ad$xxx_prop)
figh <- ggplot(data = ad[ad$depth<=200,]) +
   geom_histogram(aes(x = depth, y = xxx_prop), stat = "identity")

library(htmlwidgets)

setwd(getwd())
saveWidget(widget = fig,
           file = "chr19_coverage_all.html")

# Overdispersion
deepT <- deepT[order(deepT$start), ]
head(deepT)
library("aod")
library(MASS)
summary(deepT)
summary(m1 <- glm.nb(xxx ~ 1, data = deepT))
summary(m2 <- glm.nb(intarget ~ 1, data = deepT))
m3 <- glm(xxx ~ 1, family = "poisson", data = deepT)
AIC(m3)-AIC(m1)
pchisq(2 * (logLik(m1) - logLik(m3)), df = 1, lower.tail = FALSE)
m4 <- glm(intarget ~ 1, family = "poisson", data = deepT)
summary(m4)
AIC(m4)-AIC(m2)
plot(deepT$start, deepT$intarget)

dt <- deepT[sample(x = 1:nrow(deepT), size = 10000, replace = FALSE), ]
summary(dt)
plot(dt$start, dt$intarget)
plot(dt$start, dt$xxx)
sum(is.na(dt$xxx))
hist(dt$xxx, breaks = 100)
hist(dt$intarget, breaks = 100)
summary(m2 <- glm.nb(intarget ~ 1, data = dt[dt$intarget!=0,], control = glm.control(maxit = 500),
                     init.theta = 1.0))
m4 <- glm(intarget ~ 1, family = "poisson", data = dt[dt$intarget!=0,])
pchisq(2 * (logLik(m2) - logLik(m4)), df = 1, lower.tail = FALSE)

summary(m1 <- glm.nb(xxx ~ 1, data = dt[dt$xxx!=0,]))
m3 <- glm(xxx ~ 1, family = "poisson", data = dt[dt$xxx!=0,])

summary(m5 <- glm.nb(outarget ~ 1, data = dt[dt$outarget!=0,], control = glm.control(maxit = 500),
                     init.theta = 1.0))

head(dt)
it <- data.frame(table(dt$intarget))
head(it)
it$prop <- it$Freq/sum(it$Freq)
ggplot(data = it[-1,]) +
   geom_histogram(aes(x = Var1, y = prop), stat = "identity")








