#!/usr/bin/env Rscript

# export R_LIBS=/gscmnt/gc2802/halllab/idas/jira/BIO-1834/vendor/R/rpkgs

# this can be considered as a one-off script, it'll probably change for each project
# see bash logs for 2016-10-03 for execution details

# L I B R A R I E S ###########################################################
library('optparse')
library('ggplot2')
library('GGally')

# F U N C T I O N S ###########################################################
"%w/o%" <- function(x, y) x[!x %in% y] # see ?match

get.samples <- function(string) {
  prune <- c()
  if (length(string) > 0) {
        prune <- strsplit(string, ",")[[1]]
    }
    return(prune)
}

get.df <- function(input_file, control_samples) {
  header <- read.table(input_file, skip=1, nrows=1, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  data   <- read.table(input_file, skip=2, header=FALSE, sep="\t")
  colnames(data) <- unlist(header)
  data$TYPE <- 'call-set'
  data[data$SAMPLE %in% control_samples, 'TYPE'] <- 'control'
  return(data)
}

ggpairs.pca.plot <- function(df, title="Top 4 Principal Components") {
  p <- ggpairs(df,
               columns=c('PC1', 'PC2', 'PC3', 'PC4'),
               mapping=aes(color=TYPE),
               axisLabels="show",
               title=title)
  return(p)
}

pca.plot <- function(df, hightlight.samples, title="Top 2 Principal Components") {
  p <- ggplot(df, aes(x=PC1, y=PC2, color=TYPE)) +
       geom_point() +
       ggtitle(title,
               subtitle="By LD pruning [window=50; step=5; r^2=0.3; geno=0] & eigenstat/smartpca on SNPs") + 
       xlab("PC 1") +
       ylab("PC 2")
  return(p)
}

pca.nstitziel.plot <- function(df, hightlight.samples, title="Top 2 Principal Components") {
  outlier.set1.df <- subset(df, SAMPLE == 'H_OS-7193-7193') 
  outlier.set2.df <- subset(df, SAMPLE == 'H_OS-7541-7541') 
  outlier.set3.df <- subset(df, SAMPLE == 'H_OS-10300-10300') 
  highlight.df <- df[df$SAMPLE %in% highlight.samples, ]
  p <- ggplot(df, aes(x=PC1, y=PC2, color=TYPE)) +
       geom_point() +
       geom_text(data=outlier.set1.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, color="black", nudge_x=0.0, nudge_y=0.045, angle=90) +
       geom_text(data=outlier.set2.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, color="black", nudge_x=-0.012, nudge_y=0.0) +
       geom_text(data=outlier.set3.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, color="black", nudge_x=-0.014, nudge_y=0.0) +
       geom_point(data=highlight.df, aes(x=PC1, y=PC2), color="red") +
       geom_text(data=highlight.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, color="red", nudge_x=0, nudge_y=-0.002) +
       ggtitle(title,
               subtitle="By LD pruning [window=50; step=5; r^2=0.3; geno=0] & eigenstat/smartpca on SNPs") + 
       xlab("PC 1") +
       ylab("PC 2")
  return(p)
}

pca.plot.ctl.outliers <- function(df, control.samples, title="Top 2 Principal Components") {
  ctl.df <- df[df$SAMPLE %in% control.samples, ]
  outlier.set1.df <- subset(ctl.df, PC2 > 0.4)
  outlier.set2.df <- subset(ctl.df, PC1 == -0.3200)
  outlier.set3.df <- subset(ctl.df, PC1 == -0.3160)
  outlier.set4.df <- subset(ctl.df, PC1 == -0.1881) 
  outlier.set5.df <- subset(ctl.df, PC1 == -0.1733) 
  outlier.set6.df <- subset(ctl.df, PC1 == -0.1688) 
  outlier.set7.df <- subset(ctl.df, PC1 == -0.1238) 
  outlier.set8.df <- subset(ctl.df, PC1 == -0.1134) 
  outlier.set9.df <- subset(ctl.df, PC1 == -0.1269) 
  outlier.set10.df <- subset(ctl.df, PC1 == -0.1054) 
  outlier.set11.df <- subset(ctl.df, PC1 == -0.0883) 
  outlier.set12.df <- subset(ctl.df, PC1 == -0.0832) 
  outlier.set13.df <- subset(ctl.df, PC1 == -0.0730) 
  outlier.set14.df <- subset(df, SAMPLE == 'H_OS-1336-1336') 
  p <- ggplot(df, aes(x=PC1, y=PC2, color=TYPE)) +
       geom_point() +
       geom_text(data=outlier.set1.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, color="black", nudge_x=0.06, nudge_y=0.0) +
       geom_text(data=outlier.set2.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, color="black", nudge_x=0.06, nudge_y=-0.007) +
       geom_text(data=outlier.set3.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, color="black", nudge_x=0.06, nudge_y=0.007) +
       geom_text(data=outlier.set4.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, color="black", nudge_x=0.06, nudge_y=0.0) +
       geom_text(data=outlier.set5.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, color="black", nudge_x=0.06, nudge_y=0.000) +
       geom_text(data=outlier.set6.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, color="black", nudge_x=-0.06, nudge_y=0.000) +
       geom_text(data=outlier.set7.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, color="black", nudge_x=-0.055, nudge_y=0.000) +
       geom_text(data=outlier.set8.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, color="black", nudge_x=-0.055, nudge_y=0.000) +
       geom_text(data=outlier.set9.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, color="black", nudge_x=-0.055, nudge_y=0.000) +
       geom_text(data=outlier.set10.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, color="black", nudge_x=-0.055, nudge_y=0.000) +
       geom_text(data=outlier.set11.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, color="black", nudge_x=-0.05, nudge_y=0.000) +
       geom_text(data=outlier.set12.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, color="black", nudge_x=0.055, nudge_y=-0.007) +
       geom_text(data=outlier.set13.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, color="black", nudge_x=0.055, nudge_y=0.000) +
       geom_text(data=outlier.set14.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, color="black", nudge_x=0.04, nudge_y=0.000) +
       ggtitle(title,
               subtitle="By LD pruning [window=50; step=5; r^2=0.3; geno=0] & eigenstat/smartpca on SNPs") + 
       xlab("PC 1") +
       ylab("PC 2")
  return(p)
}

pca.plot.exp.outliers <- function(df, hightlight.samples, title="Top 2 Principal Components") {
  outlier.set1.df <- subset(df, PC1 < -0.06)
  outlier.set2.df <- subset(df, PC2 > 0.03)
  outlier.set3.df <- subset(df, PC2 > -0.01 & PC1 < -0.04 & PC1 > -0.06)
  outlier.set4.df <- subset(df, PC2 < -0.04 & PC1 < -0.04 & PC1 > -0.045)
  outlier.set5.df <- subset(df, PC2 > -0.026 & PC1 < -0.049 & PC1 > -0.055)
  outlier.set6.df <- subset(df, PC2 > -0.025 & PC1 < -0.039 & PC1 > -0.041)
  outlier.set7.df <- subset(df, PC2 > -0.018 & PC1 < -0.037 & PC1 > -0.040)
  outlier.set8.df <- subset(df, PC2 < -0.040 & PC1 < -0.059 & PC1 > -0.061)
  outlier.set9.df <- subset(df, PC2 > 0.010 & PC1 < 0.000 & PC1 > -0.009)
  outlier.set10.df <- subset(df, PC2 < -0.025 & PC1 < -0.010 & PC1 > -0.015)
  highlight.df <- df[df$SAMPLE %in% highlight.samples, ]
  p <- ggplot(df, aes(x=PC1, y=PC2)) +
       geom_point() +
       geom_text(data=outlier.set1.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, nudge_x=0.011) +
       geom_text(data=outlier.set2.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, nudge_x=-0.011, nudge_y=0.002) +
       geom_text(data=outlier.set3.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, nudge_x=0, nudge_y=0.002) +
       geom_text(data=outlier.set4.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, nudge_x=0.011) +
       geom_text(data=outlier.set5.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, nudge_x=-0.011) +
       geom_text(data=outlier.set6.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, nudge_x=-0.020) +
       geom_text(data=outlier.set7.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, nudge_x=-0.020) +
       geom_text(data=outlier.set8.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, nudge_x=-0.020) +
       geom_text(data=outlier.set9.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, nudge_x=-0.018) +
       geom_text(data=outlier.set10.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, nudge_x=0.018) +
       geom_point(data=highlight.df, aes(x=PC1, y=PC2), color="green") +
       geom_text(data=highlight.df, aes(x=PC1, y=PC2, label=SAMPLE), size=2.5, color="red", nudge_x=0, nudge_y=-0.002) +
       ggtitle(title,
               subtitle="By LD pruning [window=50; step=5; r^2=0.3; geno=0] & eigenstat/smartpca on SNPs") + 
       xlab("PC 1") +
       ylab("PC 2")
  return(p)
}

# M A I N #####################################################################
option.list <- list(
    make_option(
        c("-i", "--input"),
        action  = "store",
        type    = "character",
        default = NULL,
        help    = "an eigenstrat eigenvector/evec converted tsv"
    ),

    make_option(
        c("-o", "--output"),
        action  = "store",
        type    = "character",
        default = "pca.pdf",
        help    = "name of the output plot file"
    ),

	make_option(
        c("-c", "--control"),
        action  = "store",
        type    = "character",
        default = NULL,
        help    = "sample names to label as 'controls' from the call set"
    ),

	make_option(
        c("-p", "--highlight"),
        action  = "store",
        type    = "character",
        default = NULL,
        help    = "sample names to highlight/label from the call set"
    )
);

opts <- parse_args(OptionParser(option_list=option.list))
cat( paste("input:   ", opts$input  , "\n") )
cat( paste("output:  ", opts$output , "\n") )
cat( paste("highlight:   ", opts$highlight , "\n") )
cat( paste("control:   ", opts$control , "\n") )

# controls: H_HF-CHM1htert-US-5A,H_IJ-HG00512-HG00512_1,H_IJ-HG00513-HG00513_1,H_IJ-HG00514-HG00514_1,H_IJ-HG00731-HG00731_2,H_IJ-HG00732-HG00732_1,H_IJ-HG00733-HG00733_2,H_IJ-NA12878-NA12878_K10,H_IJ-NA12891-NA12891_D2,H_IJ-NA12892-NA12892_E1,H_IJ-NA19238-NA19238_D3,H_IJ-NA19239-NA19239_B9,H_IJ-NA19240-NA19240_F1,H_PY-CHM13-CHM13h 

options(width=200)

highlight.samples <- get.samples(opts$highlight)
control.samples <- get.samples(opts$control)

df <- get.df(opts$input, control.samples)
minus.control.df <- df[df$SAMPLE %w/o% control.samples, ]

# pdf(file=opts$output)
# #plot <- ggpairs.pca.plot(df)
# cols <- character(nrow(df))
# cols[] <- "black"
# cols[df$TYPE == "call-set"] <- "orange"
# #pairs(df[2:5], main="Top 4 Principal Components", bg=c("black", "orange")[unclass(df$TYPE)], pch=21)
# pairs(df[2:5], main="Top 4 Principal Components", col=cols, pch=21)
# dev.off()

# input: BIO-2020/data/derived/manual/2-convert-data-frame/merged.eigenstrat.pca.evec.tsv
plot <- pca.nstitziel.plot(df, highlight.samples, title="Top 2 Principal Components w/o Control Samples in Eigenstrat")
ggsave(opts$output, plot)

# plot <- pca.plot.ctl.outliers(df, control.samples, title="Top 2 Principal Components with Control Samples")
# ggsave(opts$output, plot)

# plot <- pca.plot.exp.outliers(minus.control.df, hightlight.samples, title="Top 2 Principal Components w/o Control Samples")
# ggsave('experimental-outliers.pdf', plot)

#merged.stats <- merged.stats[merged.stats$Sample %w/o% prune,] # prune out selected samples

quit(save="no", status=0)

