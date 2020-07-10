library(tidyr)

#' This function plots values from column 1 of a dataframe against the
#' corresponding values in column 2. The values can be numeric or character. 
#' 
#' @param df dataframe with two columns
#' @param logscale set to character string log10 or log2 to scale y-axis
#' @param axlabels vector of axis labels (xlabel, ylabel)
#' @param title character title of plot
#' @param xlabel.angle angle to rotate x-axis tick labels
#' @param lower.line set y-intercept for horizontal line
#' @param upper.line set y-intercept for horizontal line
#' @return Plot of column 1 values vs column 2 values
PlotXvalYval <- function(df, logscale = "", axlabels, title = "", xlabel.angle = 0, lower.line = F, 
                         upper.line = F) {
  colnames(df) <- c("x", "y")
  # Transform y values to log scale
  if (logscale=="log10"){
    df$y = log10(as.matrix(df$y + 1))
  }
  else if (logscale=="log2") {
    df$y = log2(as.matrix(df$y/10 + 1))
  }
  
  # Plot data
  p <- ggplot(data = df) + geom_point(mapping = aes(x = x, y = y)) + 
    labs(x = axlabels[1], y = axlabels[2], title = title) + 
    theme(axis.text.x = element_text(angle = xlabel.angle, hjust = 1))
    
  if (lower.line) {
    p <- p + geom_hline(aes(yintercept = lower.line), linetype = "dashed", color = "gray40")
  }
  if (upper.line) {
    p <- p + geom_hline(aes(yintercept = upper.line), linetype = "dashed", color = "gray40")
  }
  p
}

#' @param df dataframe with two columns where first column is the same string used 
#' as x-tick label and second column is a list of values for violin plot 

OneViolinPlot <- function(df, axlabels, logscale = "", jitter = FALSE, title ="") {
  # Transform y values to log scale
  colnames(df) <- c("x", "y")
  if (logscale=="log10"){
    df$y = log10(df$y + 1)
  }
  else if (logscale=="log2") {
    df$y = log2(df$y/10 + 1)
  }
  p <- ggplot(data = df, mapping = aes(x = x, y = y)) + 
    geom_violin(fill = "grey", color = "grey30") + 
    labs(x = axlabels[1], y = axlabels[2])
  if (jitter) {
    p <- p + geom_jitter(height = 0)
  }
  p
}

# Add option to set ylabels based on library prep (10x, SmartSeq)
PCPlots <- function(counts, gcdata, no.pcs=12, plotType='all', outDirectory="./figures", fileType = '.pdf', colors) { 

#' This function plots each PC score (PC1, PC2...) against cell count, reads, and # expressed
#' genes, and is colored by specified group. This allows you to see which PCs
#' are separating your groups, as well as whether different PCs encompass
#' likely technical variation. You need the "tidyr" pacakage installed.
#' 
#' @param counts counts matrix (gene x cell)
#' @param gcdata Seurat object
#' @param no.pcs number of principle components to look at
#' @param plotType set which plots to display ("density", "reads", "genes", "all") 
#' @param outDirectory output directory to save plots
#' @param fileType type of image file to save ('.png', '.pdf')
#' @param colors vector of colors, one for each group in unique(gc@ident)
#' @return Saved plots

  # calculate # reads and # expressed genes for genes that are retained in gcdata
  reads <- colSums(counts[rownames(gcdata@data), colnames(gcdata@data)])
  expressed.genes <- colSums(counts != 0)
  
  # set up data frame with cell PC scores, reads, expressed genes, and group info 
  df <- data.frame(gcdata@dr$pca@cell.embeddings[,1:no.pcs], reads = reads, expGenes = expressed.genes, group = gcdata@ident)
  #df <- data.frame(gcdata@pca.rot[,1:no.pcs], reads = reads, expGenes = expressed.genes, group = gcdata@ident)
  # Break the wide data frame to a long data frame with a column of all the PCs, for use with facet_wrap later 
  long.df <- gather(df, pc, pc_val, !! 1:no.pcs, factor_key = T) #PC indexing seems a little hacky
  # CHANGE THE ABOVE CODE ONCE TIDYR SWITCHES BACK TO EARLIER BEHAVIOR
  # tidyr 0.70 https://github.com/tidyverse/tidyr/releases
    
  n.groups <- length(unique(df$group))  # number of experimental conditions
  
  # calculate anova between groups for each PC as well as correlation, r^2, p-val between each PC and # reads and # expressed genes
  stat.data <- data.frame()
  for (i in 1:no.pcs) {
    # r^2 and correlation for reads and genes 
    fit.val.reads <- format(summary(lm(reads ~ df[, i], data = df))$adj.r.squared, digits=3)
    cor.val.reads <- format(cor(df$reads, df[, i]), digits = 3)
    p.val.reads <- format((cor.test(df$reads, df[, i]))$p.value, digits = 3)
    fit.val.genes <- format(summary(lm(expGenes ~ df[, i], data = df))$adj.r.squared, digits=3)
    cor.val.genes <- format(cor(df$expGenes, df[, i]), digits = 3)
    p.val.genes <- format((cor.test(df$expGenes, df[, i]))$p.value, digits = 3)
    
    # store results in stat.data 
    if (n.groups > 1) {
      # anova 
      aov.val <- format(summary(aov(df[,i] ~ group, data=df))[[1]]$`Pr(>F)`[1], digits=5, scientific=TRUE)
      aov.pc <- paste0('PC', i) # there is no reason this is called aov.pc other than to differentiate between big data frame pc
      stat.data <- rbind(stat.data, cbind(aov.val, aov.pc, fit.val.reads, cor.val.reads, 
                                          p.val.reads, fit.val.genes, cor.val.genes, 
                                          p.val.genes))
    } else {
      stat.data <- rbind(stat.data, cbind(fit.val.reads, cor.val.reads, p.val.reads, 
                                          fit.val.genes, cor.val.genes, p.val.genes))
    }
  }
  
  xrange <- range(long.df$pc_val)  # range of pc values
  pcs <- sapply(1:no.pcs, function(x) paste0("PC", x))  # names of PCs
  
  # PC density plots to compare PC values across multiple experimental conditions
  if ((plotType =='all' || plotType == 'density') && n.groups > 1) {
    gp <- ggplot(long.df, aes(x = pc_val, group = group, color = group)) +
      geom_density() +
      scale_color_manual(values = colors) +
      labs(y="Frequency of cells", x = "PC value") +
      theme(legend.title=element_blank(), legend.text = element_text(size = 8)) +
      theme(axis.text=element_text(size=10)) +
      scale_x_continuous(breaks=c(round(min(long.df$pc_val),1), 0, round(max(long.df$pc_val),1))) +
      facet_wrap(~pc)
    p <- ggplot_build(gp)
    text.max <- max(p$data[[1]]$y)
    gp <- gp + geom_text(data=data.frame(x=0, y=.95*text.max, label=paste0('p=',stat.data$aov.val), pc=stat.data$aov.pc),
                   aes(x=x,y=y,label=label), inherit.aes = FALSE, size=3)
    ggsave(paste0(paste0(outDirectory, 'density'), fileType), width=8, height=6, dpi=600)
    print(gp)
  }
  
  # reads vs. PCs 
  if (plotType =='all' || plotType == 'reads') {
    yrange <- range(long.df$reads)
    labs.r <- sapply(stat.data$fit.val.reads, function(x) paste0("R^2=", x))
    labs.c <- sapply(stat.data$cor.val.reads, function(x) paste0("corr=", x))
    labs.p <- sapply(stat.data$p.val.reads, function(x) paste0("p=", x))
    
    p <- ggplot(long.df, aes(x=pc_val, y=reads, group = group, color=group)) + 
      geom_point(size=.5) +
      scale_color_manual(values=colors) +
      labs(y="# UMI", x = "PC value") + 
      theme(legend.title=element_blank(), legend.text = element_text(size = 8)) + 
      theme(axis.text=element_text(size=10)) + 
      scale_x_continuous(breaks=c(round(min(long.df$pc_val),1), 0, round(max(long.df$pc_val),1))) + 
      facet_wrap(~pc) 
    p <- p + geom_text(data = data.frame(x = xrange[1]+0.4*(xrange[2]-xrange[1]), 
                                    y = yrange[2], pc=pcs, label=labs.r), 
                  aes(x,y, label=label), size = 3, inherit.aes=FALSE)
    p <- p + geom_text(data = data.frame(x = xrange[1]+0.4*(xrange[2]-xrange[1]), 
                                    y = yrange[2]-0.04*(yrange[2]-yrange[1]), pc=pcs, label=labs.c), 
                  aes(x,y, label=label), size = 3, inherit.aes=FALSE)
    p <- p + geom_text(data = data.frame(x = xrange[1]+0.4*(xrange[2]-xrange[1]), 
                                    y = yrange[2]-0.08*(yrange[2]-yrange[1]), pc=pcs, label=labs.p), 
                  aes(x,y, label=label), size = 3, inherit.aes=FALSE)
    ggsave(paste0(paste0(outDirectory, 'nUMI'), fileType), width=6, height=4, dpi=600)
    print (p)
  }
  # expressed genes vs. PCs 
  if (plotType =='all' || plotType == 'genes') {
    yrange <- range(long.df$expGenes)
    labs.r <- sapply(stat.data$fit.val.genes, function(x) paste0("R^2=", x))
    labs.c <- sapply(stat.data$cor.val.genes, function(x) paste0("corr=", x))
    labs.p <- sapply(stat.data$p.val.genes, function(x) paste0("p=", x))
    
    p <- ggplot(long.df, aes(x=pc_val, y=expGenes, group=group, color=group)) + 
      geom_point(size=.5) +
      scale_color_manual(values=colors) +
      labs(y="# expressed genes", x = "PC value") + 
      theme(legend.title=element_blank(), legend.text = element_text(size = 8)) + 
      theme(axis.text=element_text(size=10)) + 
      scale_x_continuous(breaks = c(round(min(long.df$pc_val),1), 0, round(max(long.df$pc_val),1))) + 
      facet_wrap(~pc)
    p <- p + geom_text(data = data.frame(x = xrange[1]+0.4*(xrange[2]-xrange[1]), 
                                         y = yrange[2], pc=pcs, label=labs.r), 
                       aes(x,y, label=label), size = 3, inherit.aes=FALSE)
    p <- p + geom_text(data = data.frame(x = xrange[1]+0.4*(xrange[2]-xrange[1]), 
                                         y = yrange[2]-0.04*(yrange[2]-yrange[1]), pc=pcs, label=labs.c), 
                       aes(x,y, label=label), size = 3, inherit.aes=FALSE)
    p <- p + geom_text(data = data.frame(x = xrange[1]+0.4*(xrange[2]-xrange[1]), 
                                         y = yrange[2]-0.08*(yrange[2]-yrange[1]), pc=pcs, label=labs.p), 
                       aes(x,y, label=label), size = 3, inherit.aes=FALSE)
    ggsave(paste0(paste0(outDirectory, 'genes'), fileType), width=6, height=4, dpi=600)
    print (p)
  }
}

notBoxPlot <- function(df, ylabel_string, xlabel_dir, title, ylims, lowerLine, upperLine) { 
  #' Code modified from Caroline Porter /ahg/regevdata/users/cporter/code/notBoxPlot_prettyBoxPlot.R
  #' 

  # calculate y limits 
  meanVal <- aggregate(values~vars, data=df, FUN=mean)
  sdVal <- aggregate(values~vars, data=df, FUN=sd)
  maxVal <- aggregate(values~vars, data=df, FUN=max)
  yMax <- max(max(meanVal[,2]+sdVal[,2]), max(maxVal[,2]))
  
  # plot data
  p <- ggplot(df, aes(vars, values)) + 
    stat_summary(fun.data=boxSTD, geom="boxplot", fill="steelblue", 
                 aes(width=0.7)) +
    stat_summary(fun.data=boxSEM, geom="boxplot", fill="lightblue", 
                 aes(width=0.7)) + 
    #geom_boxplot(fill="dodgerblue",outlier.colour = NA, 
    #             position = position_dodge(width=0.9), outlier.shape=NA, coef=0) +
    geom_jitter(width=0.15, shape=21,fill="gray", size=1.2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = xlabel_dir, hjust = 1)) + 
    labs(y = ylabel_string, x = "", title = title) + scale_y_log10(limits = ylims)
  #ylim(0, 20000) #+ 
  #heme(panel.border = element_blank(), panel.grid.major = element_blank(),
  #      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  # 
  if (lowerLine!=0){
    p <- p + geom_hline(aes(yintercept=lowerLine),linetype="dashed", color="gray40")
  }
  if (upperLine!=0){
    p <- p + geom_hline(aes(yintercept=upperLine),linetype="dashed", color="gray40")
  }
  p
}

boxSEM <- function(x) {
  v <- c(mean(x)-1.96*(sd(x)/sqrt(length(x))), mean(x)-1.96*(sd(x)/sqrt(length(x))), 
         mean(x), mean(x)+1.96*(sd(x)/sqrt(length(x))), mean(x)+1.96*(sd(x)/sqrt(length(x))))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}
boxSTD <- function(x) {
  v <- c(mean(x)-sd(x), mean(x)-sd(x), mean(x), mean(x)+sd(x), mean(x)+sd(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}
