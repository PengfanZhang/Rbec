#' Reference-based error correction of amplicon sequencing data
#'
#' @description
#' This function is designed for predicting the contaminated samples
#'
#' @details Ruben Garrido-Oter's group, Plant-Microbe interaction, Max Planck Institute for Plant Breeding Research
#' @author  Pengfan Zhang
#'
#' @param log_file the file contains a list of log files of each sample outputted with Rbec function
#' @param outdir output directory
#' @param outlier_constant the multiplier of variance to define the outlier
#'
#' @usage Contam_detect(log_file, outdir, outlier_constant=1.5)
#' @examples
#'
#' @import ggplot2
#' @import readr
#' 
#' @return Returns a plot showing the distribution of percentage of corrected reads across the whole sample set and a summary file recording which samples might be contaminated
#'


Contam_detect <- function(log_file, outdir, outlier_constant=1.5){

  log.file <- read_lines(log_file)

  samples <- c()
  ratio <- c()

  # read percentages of corrected reads for each sample
  for (i in seq(length(log.file))){
    corr_ratio <- read_lines(log.file[i], skip = 1)
    corr_ratio <- sub("% of reads were corrected!", "", corr_ratio)
    sample_name <- unlist(strsplit(log.file[i], "/"))
    sample_name <- sample_name[length(sample_name)-1]
    samples <- c(samples, sample_name)
    ratio <- c(ratio, corr_ratio)
  }

  ratio <- as.numeric(ratio)
  samples <- samples[order(ratio)]
  ratio <- sort(ratio)
  data4plot <- data.frame(sample = samples, ratio = ratio)
  # calculate the lower bound of percentage of corrected reads in clean samples
  lower_bd <- mean(ratio) - outlier_constant*var(ratio)

  data4plot$outlier <- "False"
  data4plot$outlier[which(data4plot$ratio < lower_bd)] <- "True"
  outlier_num <- length(which(data4plot$ratio < lower_bd))
  data4plot$sample <- factor(data4plot$sample, levels = as.character(data4plot$sample))

  # plot the distribution of percentages of corrected reads across the sample set
  p <- ggplot(data4plot, aes(sample, ratio, color = outlier)) + geom_point(size= 0.8) + geom_hline(yintercept = lower_bd, linetype = 3, size = 0.5, color = "grey") + geom_vline(xintercept = outlier_num+0.5, linetype = 3, size = 0.5, color = "grey") + theme(axis.text.x = element_text(angle=90, size=2, hjust = 1), axis.ticks.x = element_blank(), legend.position = "NA", panel.background = element_rect(colour = "black", fill="white"), panel.grid = element_blank()) + scale_color_manual(values = c("blue", "red")) + xlab("") + ylab("Percentage of corrected reads(%)")
  plot_out <- paste(outdir, "summary_plot.pdf", sep = "/")
  pdf(plot_out,width=10, height=6)
  print(p)
  dev.off()


  # output the statistic information of samples
  max_r <- round(max(ratio), digits = 2)
  min_r <- round(min(ratio), digits = 2)
  mean_r <- round(mean(ratio), digits = 2)
  output <- c(paste("The maximum percentage of corrected reads: ", max_r, sep=""), paste("The minimum percentage of corrected reads: ", min_r, sep=""), paste("Mean: ", mean_r, sep=""), " ", " ", "Potentially contaminated samples (<e2><89><a4> mean-1.5*variance):", paste(as.character(data4plot$sample[seq(outlier_num)]), round(data4plot$ratio[seq(outlier_num)], digits = 2), sep=' '))
  stat_out <- paste(outdir, "summary.log", sep = "/")
  write.table(as.data.frame(output), file=stat_out, sep="\t", quote=F, col.names = F, row.names = F)
}
