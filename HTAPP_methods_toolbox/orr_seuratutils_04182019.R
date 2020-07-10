library(slam)  # NB.var.genes run apply function in sparse matrices for finding row variance
library(MASS)  # NB.var.genes uses fitdistr

WriteDataToPortal <- function(countData, outCluster, outDir='./sc_portal_prep', outStub='outPortal',
                              padCluster=T,writeDataTable=T) {
  #' This function takes a Seurat object and clustering information and prepares files that
  #' can be uploaded to the single cell portal at https://portals.broadinstitute.org/single_cell
  #' This code works with Seurat v2.0 and was modified from Caroline Porter seurat.utils_cbmp.R.
  #' 
  #' @param countData Seurat object with gene-cell expression matrix and tSNE data
  #' @param outCluster Vector of cluster assignments for each cell in countData. The cluster labels 
  #' will likely come from countData@ident, and can be either numeric or named subsets. 
  #' These labels are placed in both the tSNE and metadata text files. 
  #' @return Three files are generated. A zipped file containing the gene-cell 
  #' expression matrix. A txt file containing number of genes and number of UMIs/cell.
  #' A txt file containing tSNE representation and cluster assignments.
  
  require(data.table)
  require(R.utils)
  
  png("temp.png")
  
  options(datatable.showProgress=T)
  
  ct = outCluster
  
  # Write expression matrix
  if (writeDataTable) {
    out.ExpTableName = sprintf('%s/%s_expData.txt',outDir,outStub)
    dataMat = data.frame(Matrix::as.matrix(countData@data))
    data.table::fwrite(as.list(c("GENE",colnames(countData@data))),file=out.ExpTableName,showProgress = F,col.names = F)
    data.table::fwrite(dataMat,file=out.ExpTableName,row.names = T,col.names = F,append = T,showProgress = T)
    gzip(out.ExpTableName)
  }
  
  # Write meta data (nGene, nUMI, cluster)
  out.metaTableName = sprintf('%s/%s_metaData.txt',outDir,outStub)
  #metaData = data.table(countData@meta.data[,1:2],keep.rownames = T)
  #metaData = data.table(countData@meta.data[,c("nGene", "nUMI")], keep.rownames = T)
  #metaData = data.table(countData@meta.data[,c("nGene", "nUMI", "orig.ident")], keep.rownames = T)
  metaData = data.table(countData@meta.data[,c("nGene", "orig.ident", "experiment", "case", 
                                               "therapy", "clonal", "annotate")], keep.rownames = T)
  
  metaData[,"typeID"] = countData@ident
  
  if (padCluster) {
    metaData[,"SeuratClustering"] = sprintf("%02d",ct)
  } else {
    metaData[,"SeuratClustering"] = as.character(ct)
  }
  data.table::fwrite(as.list(c("NAME",colnames(metaData[,-1]))),file=out.metaTableName, sep='\t')
  #zTypes = as.list(gsub("character","group",gsub("integer|double","numeric",sapply(metaData[1,-1],typeof))));
  zTypes = as.list(c("numeric","group","group","group","group","numeric", "group", "numeric", "numeric"))
  data.table::fwrite(as.list(c("TYPE",zTypes)),file=out.metaTableName,append = T,sep='\t')
  data.table::fwrite(metaData,file=out.metaTableName,append = T,sep='\t')
  
  # Write tSNE data (tSNE_1, tSNE_2, cluster)
  out.clustTableName = sprintf('%s/%s_tSNE.txt',outDir,outStub)
  metaData = data.table(row.names = countData@cell.names,X=countData@dr$tsne@cell.embeddings[, "tSNE_1"],
                        Y=countData@dr$tsne@cell.embeddings[, "tSNE_2"])
  
  if (padCluster) {
    metaData[,"SeuratClustering"] = sprintf("%02d",ct)
  } else {
    metaData[,"SeuratClustering"] = as.character(ct)
  }
  data.table::fwrite(as.list(c("NAME",colnames(metaData[,-1]))),file=out.clustTableName, sep='\t')
  zTypes = as.list(gsub("character","group",gsub("integer|double","numeric",sapply(metaData[1,-1],typeof))));
  data.table::fwrite(as.list(c("TYPE",zTypes)),file=out.clustTableName ,append = T, sep='\t')
  data.table::fwrite(metaData,file=out.clustTableName,append = T, sep='\t')
}

# test[[1]] <- NB.var.genes(test[[1]], x.high.cutoff = 10, x.low.cutoff = 0.01, do.text = TRUE, num.sd = 1.3)
NB.var.genes <- function(object, cells.use=NULL, genes.use = NULL, do.plot=TRUE,set.var.genes=TRUE,
                         x.low.cutoff=0.005, x.high.cutoff=3, diffCV.cutoff=NULL,num.sd=NULL, 
                         cex.use=0.5,cex.text.use=0.5,do.spike=FALSE,pch.use=16, col.use="black", 
                         spike.col.use="red",do.ident=FALSE, do.text=TRUE, reads.use=FALSE, 
                         cut.quantile=1, sort.results=TRUE) {
  #' This function was written by Karthik Shekhar. I slightly modified this to work with Seurat2 objects
  #' I replaced object@count.data with object@raw.data (raw UMI) and I left object@data alone (log 
  #' trasformed UMI). I replaced object@mean.var with object@hvg.info. I named the dataframe columns in
  #' mv.df mean and dispersion.
  #'
  #'x.high.cutoff is the mean expression in raw data. For SmartSeq2 (TPM), try 100 as a cutoff.
  #' For 10x (UMI), 10 is a very high cutff.
 
  print("Identifying variable genes based on UMI Counts")
  cells.use=set.ifnull(cells.use,colnames(object@data))
  genes.use=set.ifnull(genes.use, rownames(object@data))
  if (!reads.use){
    count.data=object@raw.data[genes.use,cells.use]  # modified by Orr
  } else{
    count.data = object@reads.data[genes.use,cells.use]  # this option will not work with Seurat2
  }
  
  # Empirical mean, var and CV (modified by Orr Ashenberg to deal with sparse matrices)
  # mean_emp = apply(count.data, 1, mean)  # creates a dense matrix
  # var_emp = apply(count.data, 1, var)  # creates a dense matrix
  mean_emp = Matrix::rowMeans(count.data)
  var_emp = rowapply_simple_triplet_matrix(as.simple_triplet_matrix(count.data), FUN = var)
  
  genes.use=names(mean_emp)[mean_emp > 0]
  mean_emp = mean_emp[genes.use]
  var_emp = var_emp[genes.use]
  cv_emp = sqrt(var_emp) / mean_emp
  # NB sampling
  a=Matrix::colSums(count.data)
  a = a[a <= quantile(a,cut.quantile)]
  size_factor =  a/ mean(a)
  fit=fitdistr(size_factor, "Gamma")
  if (!reads.use){
    hist(size_factor, 50, probability=TRUE, xlab="N_UMI/<N_UMI>")
  } else {
    hist(size_factor, 50, probability=TRUE, xlab="N_Reads/<N_Reads>")
  }
  curve(dgamma(x, shape=fit$estimate[1], rate=fit$estimate[2]),from=0, to=quantile(size_factor, 0.95), add=TRUE, col="red",
        main="Gamma dist fit for size factor")
  text(5,0.6, paste("shape = ", round(fit$estimate[1],2)))
  text(5,0.5, paste("rate = ", round(fit$estimate[2],2)))
  
  # Gamma distributions of individual genes are just scaled versions. If X ~ Gamma(a,b)
  # then cX ~ Gamma(a, b/c)
  a_i = rep(fit$estimate[1], length(mean_emp)); names(a_i) = names(mean_emp)
  b_i = fit$estimate[2] / mean_emp; names(b_i) = names(mean_emp)
  mean_NB = a_i / b_i; var_NB = a_i*(1+b_i) / (b_i^2)
  cv_NB = sqrt(var_NB)/mean_NB
  diffCV = log(cv_emp) - log(cv_NB)
  
  hist(diffCV,500, main="Select a delta-logCV cutoff for variable gene: ", xlab="delta-logCV")
  
  if (!is.null(num.sd)){
    diffCV.cutoff = mean(diffCV) + num.sd*sd(diffCV)
  }
  if (is.null(diffCV.cutoff)){
    diffCV.cutoff = readline("Select a delta-logCV cutoff (genes with a higher value will be considered):")
    diffCV.cutoff = as.numeric(diffCV.cutoff)
  }
  
  
  print(paste0("Using diffCV = ", diffCV.cutoff, " as the cutoff"))
  abline(v=diffCV.cutoff)
  Sys.sleep(4)
  
  print(paste0("Considering only genes with mean counts less than ", x.high.cutoff, " and more than ", x.low.cutoff))
  pass.cutoff=names(diffCV)[which(diffCV > diffCV.cutoff & (mean_emp > x.low.cutoff & mean_emp < x.high.cutoff))]
  print(paste0("Found ", length(pass.cutoff), " variable genes"))
  mv.df=data.frame(gene.mean=mean_emp,gene.dispersion=cv_emp,gene.dispersion.scaled=0)  # modified by Orr
  rownames(mv.df)=names(mean_emp)
  object@hvg.info=mv.df  # modified by Orr
  
  if (do.spike) spike.genes=grep("^ERCC", rownames(count.data), value=TRUE)
  if (do.plot) {
    plot(mean_emp,cv_emp,pch=pch.use,cex=cex.use,col="black",xlab="Mean Counts",ylab="CV (counts)", log="xy")
    curve(sqrt(1/x), add=TRUE, col="red", log="xy", lty=2, lwd=2)
    or = order(mean_NB)
    lines(mean_NB[or], cv_NB[or], col="magenta", lwd=2)
    points(mean_emp[pass.cutoff], cv_emp[pass.cutoff], col="blue", pch=16, cex=cex.use)
    
    if (do.spike) points(mean_emp[spike.genes],cv_emp[spike.genes],pch=16,cex=cex.use,col=spike.col.use)
    if(do.text) text(mean_emp[pass.cutoff],cv_emp[pass.cutoff],pass.cutoff,cex=cex.text.use)
    
    if (do.ident) {
      identify(mean_emp,cv_emp,labels = names(mean_emp))
    }
  }
  if (set.var.genes) { 
    object@var.genes=pass.cutoff
    if (sort.results) {  # Orr added Seurat code
      object@hvg.info <- object@hvg.info[order(
        object@hvg.info$gene.dispersion,
        decreasing = TRUE
      ),]
    }
    return(object)
  }
  if (!set.var.genes) return(pass.cutoff)
}

set.ifnull=function(x,y) {
  if(is.null(x)) x=y
  return(x)
}

# Code from Chris Smillie to identify variable genes across batches.

MeanCVLOESS <- function(data, num_genes=1500, use_bins=TRUE, num_bins=20, window_size=100, do.plot=FALSE) {
  #' This function written by Chris Smillie fits a LOESS curve (local regression with low-degree polynomials)
  #' to the CV vs mean of gene counts for all genes, and returns the most variable genes. 
  #' 
  #' @param data 
  #' @param num_genes 
  #' @param use_bins 
  #' @param window_size 
  #' @param do.plot 
  #' @return A vector of variable genes.
  
  # calculate mean and cv
  u = apply(data, 1, mean)
  v = apply(data, 1, var)
  i = u > 0 & v > 0
  u = u[i]
  v = v[i]
  cv = sqrt(v)/u
  
  # fit loess curve
  l = loess(log(cv) ~ log(u), family='symmetric')
  d = log(cv) - l$fitted
  
  # get variable genes
  if(use_bins == TRUE) {
    
    # select variable genes from equal frequency bins
    library(Hmisc)
    k = as.integer(num_genes/num_bins)
    var_genes = as.character(unlist(tapply(d, cut2(u, g=num_bins), function(a) names(sort(a, decreasing=T)[1:k]))))
    
  } else {
    
    # select variable genes with rolling z-scores
    library(zoo)
    
    # re-order by mean expression
    D = d[order(u)]
    
    # rolling z-scores (only use unique values)
    ru = rollapply(D, window_size, function(a) mean(unique(a)), partial=T)
    rs = rollapply(D, window_size, function(a) sd(unique(a)), partial=T)
    rz = structure((D - ru)/rs, names=names(D))
    
    # select variable genes
    var_genes = names(sort(rz, decreasing=T)[1:num_genes])
  }
  
  # plot results
  if(do.plot == TRUE){
    colors = c('#cccccc', 'black')[as.integer(names(u) %in% var_genes) + 1]
    plot(log(u), log(cv), pch=16, col=colors, xlab='log(mean)', ylab='log(cv)')
    lines(l$x[order(u)], l$fitted[order(u)], col='red', lw=2)
  }
  
  return(var_genes)
}

GetCommonVarGenes <- function(data, ident=NULL, method='loess', do.tpm=F, num_genes=1500, min_cells=5, 
                              do.plot=F, prefix=NULL, do.flatten=T, n.cores=1, ...) {
  #' This function written by Chris Smillie finds the variable genes that are common across a group of
  #' datasets. For each dataset, this function calls another function to calculate the most highly
  #' variable genes. Then the function counts how often each gene is marked as variable across each of the 
  #' datasets. Those variable genes that are found most commonly across all datasets are returned. 
  #' This approach can help correct for batch effects, by not using variable genes in PCA and clustering 
  #' that are unique to only one or a few datasets. Genes that are variable across multiple datasets are 
  #' likely better at describing the common variation across the datasets.
  #' https://github.com/cssmillie/code/blob/master/single_cell/var_genes.r
  #' 
  #' @param data Seurat object if method == "seurat" or "karthik", or object@raw.data otherwise.
  #' @param ident Split cells by cell identity
  #' @param method 
  #' @param do.tpm 
  #' @param num_genes 
  #' @param min_cells 
  #' @param do.plot 
  #' @param prefix 
  #' @param do.flatten 
  #' @param n.cores 
  #' @return A vector of variable genes.
  
  #source('~/code/single_cell/parallel.r')  # Orr commented out
  
  # Split by cell identity
  if(is.null(ident)){ident = rep(1, ncol(counts))}
  ident = as.factor(as.character(ident))
  print(table(ident))
  
  var_genes = sapply(levels(ident), function(i) { 
    print(i)
    
    # Start plotting device
    if(!is.null(prefix)){pdf(paste(prefix, i, 'pdf', sep='.'), w=1000, h=800)}
    
    # Sample data to get individiual identity classes. Either Seurat object or matrix data.
    if (method == 'seurat' | method == 'karthik') {  # Seurat object
      data = SubsetData(data, cells.use = rownames(data@meta.data[ident == i, ]))
      print(dim(data@data))
    } else {  # matrix data
      data = data[,ident == i]
      if(do.tpm){data = calc_tpm(data=data)}
      genes.use = Matrix::rowSums(data > 0) >= min_cells	    
      data = data[genes.use,]
      print(dim(data))
    }
    
    # Variable gene selection methods
    if (method == 'loess') {
      vi = MeanCVLOESS(data, num_genes=num_genes, do.plot=do.plot, ...)
    }
    
    if(method == 'adam'){  # CAREFUL, I DO NOT HAVE THESE FUNCTIONS CURRENTLY
      source('~/dev/adam/rna_seq/r/samples.r')
      source('~/dev/adam/rna_seq/r/util.r')
      vi = get.variable.genes(data, do.plot=do.plot, ...)
      i = (!is.na(var_genes[,4])) & (var_genes[,4] <= .05)
      vi = vi[i,1]
    }
    if(method == 'karthik'){  
      # Chris original code for Karthik function 
      # vi = meanCVfit(data, diffCV.num_genes=num_genes, do.plot=do.plot, ...)
      # Orr code for Karthik function, which requires a Seurat object as input.
      vi = NB.var.genes(data, x.high.cutoff = 10, x.low.cutoff = 0.01, do.text = TRUE, 
                        num.sd = 0.7, do.plot=F, set.var.genes = F)
    }
    if(method == 'seurat'){
      vi = FindVariableGenes(data, mean.function = ExpMean, dispersion.function = LogVMR, 
                             x.low.cutoff = 0.1, x.high.cutoff = 7, y.cutoff = 0.5, do.plot = T, ...)
      vi = vi@var.genes
      print(length(vi))
    }
    
    # Stop plotting device
    if(!is.null(prefix)){dev.off()}
    
    return(vi)
  })
  
  # Find most common variable genes
  if(do.flatten == TRUE){
    a = sort(table(as.character(unlist(var_genes))), decreasing=T)
    num_genes = min(length(a), num_genes)
    k = a[num_genes]  # cutoff in number of occurrences of variable genes
    u = names(a)[a > k]  # genes above cutoff
    v = sample(names(a)[a == k], num_genes - length(u))  # genes at cutoff
    var_genes = c(u,v)
  }
  
  return(var_genes)
}


QCSS2 <- function(genome.mapping.file, transcriptome.mapping.file, qc.out.file, genic.out.file, sampleid, gcdata = NULL, 
                  exp.title = NULL)  {
  #' This function takes the genome_mapping_summary and transcriptome_mapping_summary files produced
  #' when running RNASeq-QC on SmartSeq2 data in the KCO RNAseq pipeline. It summarizes some of these 
  #' quality metrics related to exon content and number of genes per cell.
  #' If a Seurat object is provided, the function updates a few of the quality metrics. 
  #' The Seurat object contains the TPM gene expression matrix that has been filtered 
  #' for low-quality genes and cells. The function also plots the fraction of reads mapping
  #' to genic, exonic, intronic, and intergenic regions for each experimental replicate.
  #' 
  #' @param genome.mapping.file This file is generated by
  #' /seq/regev_genome_portal/SOFTWARE/KCO/RNASEQ_pipeline/generate_sample_summary_stats.pl.
  #' It contains information on how well reads aligned to the genome.
  #' @param transcriptome.mapping.file This file is generated by
  #' /seq/regev_genome_portal/SOFTWARE/KCO/RNASEQ_pipeline/util/summarize_rnaseqQC_results.pl.
  #' It contains information on how well genome-mapped reads align to known transcript features.
  #' @param qc.out.file File name for tsv file containing quality control metrics that have been 
  #' collected across cells.
  #' @param genic.out.file File name for plot of fraction of reads mapping to different genomic regions.
  #' @param sampleid Vector of sample names that should match the sample names in @param genome.mapping.file 
  #' and @param transcriptome.mapping.file. In addition if using @param gcdata, these names should match those
  #' in the gcdata@meta.data slot sample column. 
  #' @param gcdata A Seurat objet. The @meta.data slot should have a column named sample, which has 
  #' the sample name for each cell.
  #' @param exp.title Plot title for @param genic.out.file 
  #' @return A list with two elements: a dataframe with the quality control metrics, named table, and a
  #' ggplot object of the fraction of reads mapping to different genomic regions, named plot.genic .
  
  # Read genome_mapping_summary, taking specific columns
  qc.genome.vals <- read.table(genome.mapping.file, sep = "\t", stringsAsFactors = F)
  colnames(qc.genome.vals) <- strsplit(readLines(genome.mapping.file, n = 1), "\t")[[1]]
  qc.genome.vals <- qc.genome.vals[, c(1, 2, 5:8, 13)]
  colnames(qc.genome.vals)[1] <- "sample"  # used for merging qc.genome.vals and qc.transcriptome.vals
  colnames(qc.genome.vals)[7] <- "Number of Genes"
  
  # Read transcriptome_mapping_summary, taking specific columns
  qc.transcriptome.vals <- read.table(transcriptome.mapping.file, skip = 2, sep = "\t", stringsAsFactors = F)
  qc.colnames <- strsplit(readLines(transcriptome.mapping.file, n = 2)[2], "\t")[[1]][c(7:8, 10:11, 13, 28:31)]
  qc.transcriptome.vals <- qc.transcriptome.vals[, c(7:8, 10:11, 13, 28:31)]
  colnames(qc.transcriptome.vals) <- qc.colnames
  colnames(qc.transcriptome.vals)[1] <- "sample"  # used for merging qc.genome.vals and qc.transcriptome.vals
  
  # Merge genome and transcriptome mapping summaries based on the sample names, which are shared.
  qc <- merge(qc.genome.vals, qc.transcriptome.vals, by = "sample")
  qc <- qc[match(sampleid, qc$sample), ]  # sort rows using order of samples in sampleid
  
  # If Seurat object available, update a few quality metrics after the filtering.
  if (!is.null(gcdata)) {  
    qc$`Number of Genes` <- NULL  # no longer need the estimate from genome_mapping_summary
    qc <- merge(gcdata@meta.data[, c("sample", "nGene")], qc, by = "sample")
    colnames(qc)[2] <- "Number of Genes"
    qc <- qc[match(gcdata@meta.data$sample, qc$sample), ]  # sort rows using order of samples in Seurat object
    qc$sample <- as.character(qc$sample)  # do not use factor as this leads to warnings when adding the median values later
  }
  
  # Plot fraction of reads mapping to different genomic regions.
  # Store fraction of reads mapping to different genomic regions in dataframe.
  qc.colnames <- c("genic", "exonic", "intronic", "intergenic")
  qc.genic <- 100 * qc[ ,c("Intragenic Rate", "Exonic Rate", "Intronic Rate", "Intergenic Rate")]
  colnames(qc.genic) <- qc.colnames
  
  # Gather QC dataframe into genomic region : percent mapping key-value pairs for making violin plots. 
  qc.violin <- gather(qc.genic, "region", "percent", c(1:4))
  # Set order in which genic regions and samples are plotted in barplot
  qc.violin$region <- factor(qc.violin$region, levels = qc.colnames[1:4]) 
  p <- ggplot(data = qc.violin, mapping = aes(x = region, y = percent)) + geom_violin(fill = "grey", color = "grey30") + 
    geom_jitter(height = 0, size = 1) +
    scale_y_continuous(breaks=seq(0,100,10), limits = c(0, 100)) + 
    ylab("percent bases mapping") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(exp.title)
  ggsave(genic.out.file, width = 6, height = 4)
  
  # Add median value of each column, other than the first column which contains the names for each sample.
  median.vals <- apply(qc[-c(1)], 2, median)
  qc <- rbind(c("median", median.vals), qc)
  
  # Write qc table to file.
  write.table(t(qc), qc.out.file, sep = "\t", quote = F)
  
  # Return qc table and genic plot objects.
  list(table = qc, plot.genic = p)
}

QCDropSeq <- function(exon.files, gene.files, qc.out.file, genic.out.file, gcdata = NULL, 
                      exp.labels = NULL, exp.title = NULL)  {
  #' This function takes two files produced by the DropSeq pipeline and 
  #' extracts quality metrics related to exon content and number of genes per cell barcode.
  #' If a Seurat object is provided, the function updates a few of the quality metrics. 
  #' The Seurat object contains the digital gene expression matrix that has been filtered 
  #' for low-quality genes and cells. The function also plots the fraction of reads mapping
  #' to genic, exonic, intronic, and intergenic regions for each experimental replicate.
  #' 
  #' @param exon.files Vector of file names with percent of bases encoding introns and exons. Names 
  #' of vector elements should be the sample names.
  #' @param gene.files Vector of file names with total number of genes and transcript counts (total UMI) 
  #' per cell barcode. Names of vector elements should be the sample names.
  #' @param qc.out.file File name for tsv file containing quality control metrics that have been 
  #' concatenated across several samples.
  #' @param genic.out.file File name for plot of fraction of reads mapping to different genomic regions.
  #' @param gcdata A list of Seurat objects, one for each sample. This should follow the same order of
  #' samples as in @param exon.file and @param gene.file. The names of the gcdata list elements 
  #' should be the sample names.
  #' @param exp.labels Alternative labels to use in @param genic.out.file for the sample names. These
  #' should have the same length as the sample names in @param exon.files
  #' @param exp.title Plot title for @param genic.out.file 
  #' @return A list with two elements: a dataframe with the quality control metrics, named table, and a
  #' ggplot object of the fraction of reads mapping to different genomic regions, named plot.genic .
  qc <- data.frame()
  sampleid <- names(exon.files)  # names of experimental replicates
  for (sample in sampleid) {
    # Read exon information
    lines <- readLines(exon.files[sample]) 
    qc.exon.vals <- strsplit(lines[8], "\t")[[1]][11:16]
    qc.exon.vals <- round(as.numeric(qc.exon.vals), digits = 2)
    names(qc.exon.vals) <- strsplit(lines[7], "\t")[[1]][11:16]
    
    # Fraction of bases passing filter that aligned = PF_ALIGNED_BASES/PF_BASES
    qc.aligned <- strsplit(lines[8], "\t")[[1]][1:2] # PF_BASES PF_ALIGNED_BASES 
    qc.aligned <- as.numeric(qc.aligned)
    qc.aligned <- round(qc.aligned[2] / qc.aligned[1], digits = 2)
    names(qc.aligned) <- "PCT_GENIC_BASES"
    
    # Read total gene and UMI per barcode summary
    qc.gene.vals <- read.table(gene.files[sample], skip = 2, header = T)
    qc.gene.vals <- round(c(nrow(qc.gene.vals), median(qc.gene.vals$NUM_GENES), 
                            median(qc.gene.vals$NUM_TRANSCRIPTS)))
    names(qc.gene.vals) <- c("Estimated Number of Cells", "Median Genes per Cell", 
                             "Median UMI Counts per Cell")
    
    qc <- rbind(qc, t(data.frame(c(qc.gene.vals, qc.aligned, qc.exon.vals))))  # add qc for this sample to dataframe
  }
  rownames(qc) <- sampleid
  
  # If Seurat object available, update a few quality metrics after the filtering.
  if (!is.null(gcdata)) {  
    stopifnot(all.equal(names(gcdata), names(exon.files)))  # check the samples correspond to one another
    qc$`Estimated Number of Cells` <- sapply(gcdata, function(data) length(data@ident))
    qc$`Median Genes per Cell` <- sapply(gcdata, function(data) as.integer(median(data@meta.data$nGene)))
    qc$`Median UMI Counts per Cell` <- sapply(gcdata, function(data) as.integer(median(data@meta.data$nUMI)))
    qc$`Total Genes Detected` <- sapply(gcdata, function(data) nrow(data@data))
    qc <- qc[, c(1:3, 11, 4:10)]
  }
  
  # Write qc table to file
  write.table(t(qc), qc.out.file, sep = "\t", quote = F)
  
  # Plot fraction of reads mapping to different genomic regions.
  # Store fraction of reads mapping to different genomic regions in dataframe.
  qc.colnames <- c("genic", "exonic", "intronic", "intergenic", "ribosomal", "replicate")
  qc.genic <- 100 * qc[,c("PCT_GENIC_BASES", "PCT_CODING_BASES", "PCT_INTRONIC_BASES", 
                          "PCT_INTERGENIC_BASES", "PCT_RIBOSOMAL_BASES")]
  qc.genic$replicate <- sampleid
  colnames(qc.genic) <- qc.colnames
  
  # Gather QC dataframe into genomic region : percent mapping key-value pairs for making bar plots. 
  qc.bar <- gather(qc.genic, "region", "percent", c(1:5))
  # Set order in which genic regions and samples are plotted in barplot
  qc.bar$region <- factor(qc.bar$region, levels = qc.colnames[1:5])  
  qc.bar$replicate <- factor(qc.bar$replicate, levels = sampleid)
  if (is.null(exp.labels)) {exp.labels <- sampleid}
  p <- ggplot(data = qc.bar) + geom_bar(stat = "identity", mapping = aes(x = region, y = percent, 
                                                                         fill = replicate),
                                        position = position_dodge(), color = "black") + 
    scale_y_continuous(breaks=seq(0,100,10), limits = c(0, 100)) + 
    ylab("percent bases mapping") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    scale_fill_discrete(name="Experimental\nCondition", breaks=sampleid, labels=exp.labels) + 
    ggtitle(exp.title)
  ggsave(genic.out.file, width = 6, height = 4)
  
  list(table = qc, plot.genic = p)
}

QC10x <- function(qc.files, figures.dir, gcdata = NULL, exp.labels = NULL, exp.title = NULL) {
  #' This function takes the metrics_summary.csv files made from several runs of Cell Ranger count 
  #' and concantenates them into a single file. If a Seurat object is provided, the function 
  #' updates a few of the quality metrics. The Seurat object contains the digital gene expression matrix 
  #' that has been filtered for low-quality genes and cells. The function also plots the fraction of 
  #' reads mapping to genic, exonic, intronic, and intergenic regions for each experimental replicate.
  #' This code was written for Cell Ranger 2.1.0 and Seurat v2.3.2.
  #' 
  #' Each metrics_summary.csv file contains information on each sequenced molecule.
  #' [molecule_info.h5](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/molecule_info) 
  #' has information on each molecule from the 10x run. Each molecule is a unique (cell barcode, UMI, gene) 
  #' tuple and there is other associated information like number of reads and quality.
  #' @param qc.files Vector of metrics_summary.csv file names, one for each sample, created by 
  #' Cell Ranger count. The names of the vector elements should be the sample names.
  #' @param figures.dir Directory where to write figures.
  #' @param gcdata A list of Seurat objects, one for each sample. This should follow the same order of
  #' samples as in @param qc.files. The names of the gcdata list elements should be the sample names.
  #' The metadata slot of each Seurat object should have a nRead column generated by the function
  #' AddMoleculeInfoSeuratObject.
  #' @param exp.labels Alternative labels to use in @param qc.genic.file for the sample names. These
  #' should have the same length as the sample names in @param exon.files
  #' @param exp.title Plot title for qc.genic.file 
  #' @return A list with two elements: a dataframe with the quality control metrics, named table, and a
  #' ggplot object of the fraction of reads mapping to different genomic regions, named plot.genic .
  
  # Files to write containing 10x quality control metrics.
  # 1. tsv file containing quality control metrics that have been concatenated across several samples.
  # 2. Plot of fraction of reads mapping to different genomic regions.
  # 3. Plots of reads and genes per cell.
  qc.table.file <- paste0(figures.dir, "metrics_summary_all.tsv")
  qc.genic.file <- paste0(figures.dir, "qc_genic_regions.pdf")
  qc.plotstats.file <- paste0(figures.dir, "qc_plotstats.pdf")
  
  # Iterate over quality control metrics files, and collect the metrics for each sample into a dataframe.
  qc <- data.frame()
  sampleid <- names(qc.files)   # names of experimental replicates
  for (sample in sampleid) {
    qc.add <- read.table(qc.files[sample], header = T, sep = ",", 
                         check.names = F, stringsAsFactors = F)[,c(1:4, 6, 11:17, 19:20)]
    qc <- rbind(qc, qc.add)
  }
  qc$`Number of Reads` <- as.numeric(gsub(",", "", as.character(qc$`Number of Reads`)))  # remove commas
  rownames(qc) <- sampleid
  
  # If Seurat object available, update a few quality metrics after the filtering.
  if (!is.null(gcdata)) {  
    stopifnot(all.equal(names(gcdata), names(qc.files)))  # check the samples correspond to one another
    qc$`Estimated Number of Cells` <- sapply(gcdata, function(data) length(data@ident))
    qc$`Median Genes per Cell` <- sapply(gcdata, function(data) as.integer(median(data@meta.data$nGene)))
    qc$`Median UMI Counts per Cell` <- sapply(gcdata, function(data) as.integer(median(data@meta.data$nUMI)))
    qc$`Mean Reads per Cell` <- sapply(gcdata, function(data) as.integer(mean(data@meta.data$nRead)))
    qc$`Total Genes Detected` <- sapply(gcdata, function(data) nrow(data@data))
  }
  
  # Write quality metrics for all samples to file.
  write.table(t(qc), qc.table.file, sep = "\t", quote = F)
  
  # Plot fraction of reads mapping to different genomic regions.
  # For each sample, store fraction of reads mapping to different genomic regions in dataframe.
  qc.colnames <- c("genome", "exon", "intron", "intergene", "transcriptome", "sample")
  # Remove percent symbols from fraction mapping reads.
  # If there is only one sample, the result of apply is a vector, which then needs to be transposed.
  if (length(sampleid) > 1) {
    qc.genic <- data.frame(apply(qc[,c("Reads Mapped Confidently to Genome",
                                       "Reads Mapped Confidently to Exonic Regions", 
                                       "Reads Mapped Confidently to Intronic Regions", 
                                       "Reads Mapped Confidently to Intergenic Regions",
                                       "Reads Mapped Confidently to Transcriptome")], 
                                 2, function(x) as.numeric(sub("%", "", x)))) 
  } else {
    qc.genic <- data.frame(t(apply(qc[,c("Reads Mapped Confidently to Genome",
                                         "Reads Mapped Confidently to Exonic Regions", 
                                         "Reads Mapped Confidently to Intronic Regions", 
                                         "Reads Mapped Confidently to Intergenic Regions",
                                         "Reads Mapped Confidently to Transcriptome")], 
                                   2, function(x) as.numeric(sub("%", "", x))))) 
  }
  if (is.null(exp.labels)) {exp.labels <- sampleid}
  qc.genic$sample <- exp.labels  # mark what sample these statistics come from
  colnames(qc.genic) <- qc.colnames
  
  # Gather QC dataframe into genomic region : percent mapping key-value pairs for making bar plots. 
  qc.bar <- gather(qc.genic, "region", "percent", c(1:5))
  # Set order in which genic regions and samples are plotted in barplot
  qc.bar$region <- factor(qc.bar$region, levels = qc.colnames[1:5])  
  qc.bar$sample <- factor(qc.bar$sample, levels = exp.labels)
  plot.genic <- ggplot(data = qc.bar) + geom_bar(stat = "identity", mapping = aes(x = region, y = percent, 
                                                                                  fill = sample),
                                                 position = position_dodge(), color = "black") + 
    scale_y_continuous(breaks=seq(0,100,10), limits = c(0, 100)) + 
    ylab("percent bases mapping") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    scale_fill_discrete(name="Experimental\nCondition") + ggtitle(exp.title)
  ggsave(qc.genic.file, width = 6, height = 5)
  
  # Plots to show for each sample total number of sequencing reads, total number of cells, and
  # distributions for number of reads per cell and number of median genes per cell. 
  df.plot <- data.frame(nRead = qc$`Number of Reads`, nCell = qc$`Estimated Number of Cells`, 
                        sample = exp.labels, stringsAsFactors = F)
  df.plot$sample <- factor(df.plot$sample, levels = exp.labels)
  gcdata.merge <- MergeMultipleSeuratObjects(gcdata, project = "qc")
  df.violin <- data.frame(nRead = gcdata.merge@meta.data$nRead, nGene = gcdata.merge@meta.data$nGene,
                          ident = gcdata.merge@ident)
  df.violin$ident <- factor(df.violin$ident, levels = exp.labels)
  
  p1 <- ggplot(data = df.plot) + geom_bar(stat = "identity", mapping = aes(x = sample, y = nRead, fill = sample)) + 
    ylab("Number of Reads") + guides(fill=FALSE)
  p2 <- ggplot(data = df.plot) + geom_bar(stat = "identity", mapping = aes(x = sample, y = nCell, fill = sample)) + 
    ylab("Estimated Number of Cells") + guides(fill=FALSE)   
  p3 <- ggplot(data = df.violin, mapping = aes(x = ident, y = nRead)) + 
    geom_violin(aes(fill = factor(ident)), show.legend = F) +
    geom_boxplot(width = 0.1, outlier.shape = 1) + 
    labs(x = "sample", y = "Reads per Cell", title = "")
  p4 <- ggplot(data = df.violin, mapping = aes(x = ident, y = nGene)) + 
    geom_violin(aes(fill = factor(ident)), show.legend = F) +
    geom_boxplot(width = 0.1, outlier.shape = 1) + 
    labs(x = "sample", y = "Genes per Cell", title = "")
  plot.stat <- plot_grid(p1, p2, p3, p4, align = 'h', labels = c('A', 'B', 'C', 'D'), nrow = 2)
  save_plot(qc.plotstats.file, plot.stat, base_aspect_ratio = 1.8, base_height = 12)
  
  list(table = qc, plot.genic = plot.genic, plot.stat = plot.stat)
}

MergeMultipleSeuratObjects <- function(gcdata, project) {
  #' This function uses the MergeSeurat function to merge two or more Seurat objects.
  #' IMPORTANT note: do.normalize = F flag in Seurat::MergeSeurat() means that no normalization
  #' or scaling information (gcdata@data and gcdata@scale.data) exists in the merged Seurat object.
  #' gcdata@scale.data = NULL and gcdata@data = gcdata@raw.data. Normalization and scaling must be
  #' performed on the merged object before PCA.
  #' IMPORTANT: also note that Seurat chooses the initial identity class for each cell based on the 
  #' first field from the cell's column name, where "_" is used as a delimiter in the string split.
  #' 
  #' @param gcdata List of Seurat objects where the object elements are named by their sample.
  #' @param project Project name to give merged Seurat object.
  #' @return Merged Seurat object.
  if (length(gcdata) == 1) {  # if list of Seurat objects only has 1 element, no merge needs to be done
    return(gcdata)
  }
  
  # Iteratively merge Seurat objects two at a time.
  gcdata.merge <- gcdata[[1]]
  for(i in 2:length(gcdata)){
    gcdata.merge <- MergeSeurat(object1 = gcdata.merge, object2 = gcdata[[i]],
                                do.scale = F, do.center = F, do.normalize = F)
  }
  gcdata.merge
  
  # # CODE USED FOR VERSIONS OF SEURAT BEFORE v2.3
  # # Iteratively merge Seurat objects two at a time.
  # sampleid <- names(gcdata)  # sample names
  # gcdata.merge <- MergeSeurat(object1 = gcdata[[1]], object2 = gcdata[[2]],
  #                             add.cell.id1 = sampleid[1], add.cell.id2 = sampleid[2],
  #                             project = project, do.normalize = F)
  # for (i in seq(gcdata)[-c(1, length(gcdata))]) {  # skip first and last index
  #   gcdata.merge <- MergeSeurat(object1 = gcdata.merge, object2 = gcdata[[i+1]],
  #                               add.cell.id2 = sampleid[i+1],
  #                               project = project, do.normalize = F)
  # }
  # 
  # # Seurat has issues naming the identity class when there are underscores in the class names. 
  # # This is an issue in the @ident slot, where it will only keep part of the name. This can be 
  # # modified using names.field and names.delim but is annoying. For some reason, 
  # # @meta.data$orig.ident does not have this issue, so we use it to name the @ident class for each cell.
  # gcdata.merge <- SetIdent(gcdata.merge, ident.use = gcdata.merge@meta.data$orig.ident)
  # gcdata.merge
}

MultipleFeaturePlot <- function(gcdata, features="nGene", feature.type="meta", ncols=ceiling(sqrt(length(features))), pt.size=1, same.scale=FALSE) {
  #' This function takes a feature from the @meta.data of a Seurat object, or from the un-scaled gene expression 
  #' @data of a Seurat object, and plots it on a tSNE with a blue to yellow to red gradient of colors. 
  #' If more than one feature is supplied, the graphs are sub-plotted. Metadata features and genes cannot be plotted
  #' simultaneously. This extends the Seurat::FeaturePlot() function
  #' 
  #' @param gcdata Seurat object with @data and @meta.data slots.
  #' @param features Either a vector of features to be plotted from @meta.data, or a vector of gene names to be
  #' plotted from @data.
  #' @param feature.type Either "gene" or "meta" to indicate the type of feature, @data or @meta.data respectively.
  #' @param ncols Number of columns desired for facet wrap when making subplots.
  #' @param pt.size Size of points in plot.
  #' @param same.scale TRUE if color bar should be scaled based on entire dataset rather than just 
  #' the expression for plotted features. This is only used when @param feature.type = "gene".
  #' @return ggplot object.
  #' 
  #' Written by Caroline Porter and color functions provided by Sam Riesenfield in original function plot.tsne.feature. 
  #' I modified the code to better deal with plotting more than one @meta.data feature. 
  #' 
  #' TO DO: Add code to check that the gene or feature actually exists in the data set, throw error if not.
  #' Add option to scale colar bar baed on entire dataset, rather than the plotted sub-set
  
  # Load required libraries 
  library(tidyr)
  
  # Get color gradient (THIS COULD BE UPDATED TO ALLOW ADDITIONAL INPUTS PER SAM'S SCRIPT)
  # source("/ahg/regevdata/users/cporter/code/colrs.fromSamRiesen.R")
  # colors <- get.hmap.col()
  colors<-c("#191970","#121285","#0C0C9A","#0707B0","#0101C5","#0014CF","#0033D3","#0053D8","#0072DD","#0092E1","#00B2E6",
            "#00D1EB","#23E8CD","#7AF17B","#D2FA29","#FFEB00","#FFC300","#FF9B00","#FF8400","#FF7800","#FF6B00","#FF5F00","#FF5300",
            "#FF4700","#F73B00","#EF2E00","#E62300","#DD1700","#D50B00","#CD0000")
  
  # Error out if user does not specify correct feature type 
  if (feature.type!="meta" & feature.type!="gene") {
    stop("feature type must be 'meta' or 'gene'")
  }
  
  # If plotting gene expression, verify that the desired genes to plot exist in the @data slot
  if (feature.type == "gene") {
    features <- features[which(features %in% rownames(gcdata@data))]
  }
  if (!length(features)) {
    stop("None of the genes requested for plotting are in gcdata@data.")
  }
  
  # Collect feature info either from @meta.data or @data
  if (feature.type=="meta") { 
    feature.info <- as.matrix(gcdata@meta.data[,features])   # column format
    if (length(features) == 1) {
      colnames(feature.info) <- features
    }        
  } else if (feature.type=="gene") { 
    feature.info <- as.matrix(gcdata@data[features,])  # row format
    if (length(features) > 1) {
      feature.info <- t(feature.info)  # transpose into more convenient column format.
    }        
    colnames(feature.info) <- features  # only required when length of features is 1.
  }
  
  # Build data frame of feature info and tSNE coordinates 
  tmp.df <- data.frame(feature.info, gcdata@dr$tsne@cell.embeddings)
  plot.df <- gather(tmp.df, name, val, 1:length(features), factor_key=TRUE)
  
  # Set scales for color bar 
  if (same.scale==TRUE & feature.type=="gene") {
    lower=min(gcdata@data) 
    upper=max(gcdata@data)
  } else if (same.scale==TRUE & feature.type=="meta") {
    stop("same.scale can only be set to true when plotting gene features")
  } else {
    lower=min(plot.df$val) 
    upper=max(plot.df$val)
  }
  
  # Color tSNE plot by gene expression or by meta data
  p <- ggplot(plot.df, aes(x=tSNE_1, y=tSNE_2)) + geom_point(aes(color=val), alpha=0.8, shape=16, size=pt.size) + 
    theme(aspect.ratio = 1) + scale_color_gradientn(colors=colors, limits=c(lower, upper))
  p <- p + theme(aspect.ratio=1, text = element_text(size=10), axis.text=element_text(size=6), 
                 strip.text.x = element_text(margin = ggplot2::margin(.1, 0, .1, 0, "cm")),
                 strip.text = element_text(size=12)) + 
    facet_wrap( ~ name, ncol=ncols)
  return(p)
}

RankModuleScore <- function(gcdata, name.cluster, gs, plotdir, ctrl.size=10) {
  #' The purpose of this function is to help somewhat automate assigning cell types to
  #' cell clusters. For each cell cluster, the cluster is scored by annotated sets of genes
  #' using Seurat::AddModuleScore(). Each cell is individually scored, and then the average 
  #' score for all cells within a cluster is recorded. The gene sets that score the highest 
  #' for a given cluster can potentially reveal what cell type that cluster is. This function 
  #' requires a Seurat object containing gene expression values, cells that have been clustered, and a 
  #' tSNE dimensional reduction to display the gene set scores mapped onto the cells. 
  #' The function also requires a list of gene sets to score each cell cluster with. After
  #' running this function, it is good to visually inspect the created tSNE plots to see how specific
  #' the best scoring gene sets are for each cell cluster of interest. 
  #' 
  #' This function uses dplyr functions, so make sure the plyr package is detached before running.
  #' 
  #' @param gcdata Seurat object with @data and @meta.data slots. Clustering must have been performed
  #' and the cluster identities must be stored as a column of gcdata@meta.data. In addition,
  #' there must be a tSNE dimensional reduction slot for making plots.
  #' @param name.cluster Name of column in gcdata@meta.data with cluster labels. This defines the
  #' cell clusters that are scored by each gene set.
  #' @param gs This is a GSA (R pacckage) object containing the names of each gene set and a list
  #' of genes that define that corresponding gene set. The object can be created using  
  #' gs <- GSA.read.gmt(file.gmt). gmt files are further described 
  #' [here](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats).
  #' @param plotdir For each gene set, a tSNE plot is colored by the scores from that gene set and placed
  #' in this plotting directory.
  #' @param ctrl.size When calculating the average expression level of each gene set, Seurat::AddModuleScore()
  #' subtracts off expression from a control set of genes. This parameters defines how many control genes to
  #' take from each expression bin.
  #' @return The function returns the dataframe cluster.scores. The first column is the cluster names from
  #' @param name.cluster, and the subsequent columns are the mean gene set scores for all cells within 
  #' that cluster.
  
  # Create directory to store tSNE plots colored by mean gene set expression.
  dir.create(plotdir, showWarnings = FALSE) 
  
  # Create data frame where first column is cluster number and
  # subsequent columns are mean scores for a given gene set for all cells in given cluster.
  cluster.scores <- data.frame("clusters" = unique(gcdata@meta.data[,name.cluster]))
  names(cluster.scores)[1] <- name.cluster
  
  # File where each row of dataframe mouse.human is a pair of mouse-human gene homologs.
  # file.homolog <- "/Volumes/ahg_regevdata/users/orr/data/genelists/mouse_human_mapping.RData"
  # load(file.homolog)
  
  # Iterate over each gene set and score each cell.
  for (i in seq(gs$geneset.names)) {
    name <- gs$geneset.names[i]
    geneset <- gs$genesets[[i]]
    # geneset <- MouseHumanHomolog(mouse.human, geneset, "human")  # code to add to process homologs
    print(c(name, geneset))
    
    # Score cells by gene sets.
    # gcdata <- Seurat::AddModuleScore(gcdata, genes.list = list(geneset), ctrl.size = ctrl.size, 
    #                                  enrich.name = name)
    tryCatch({ 
      gcdata <- Seurat::AddModuleScore(gcdata, genes.list = list(geneset), ctrl.size = ctrl.size, enrich.name = name)
      #oldname <- paste0(name, "1")  # Remove 1 from end of name
      names(gcdata@meta.data)[names(gcdata@meta.data) == paste0(name, "1")] <- name
      
      # Rank clusters by average gene set score for cells in the cluster.
      # http://dplyr.tidyverse.org/articles/programming.html
      summ_name <- paste0("mean_", name)
      add.scores <- gcdata@meta.data %>% group_by(!!as.name(name.cluster)) %>% 
        summarise(!!summ_name := mean(!!as.name(name)))
      cluster.scores <- merge(cluster.scores, add.scores, by=name.cluster)  # need common name for cluster column
      
      # tSNE plot colored by gene set scores.
      # MultipleFeaturePlot(gcdata, features = name, feature.type = "meta", pt.size = 2)
      # ggsave(paste0(plotdir, "tSNE_", name, ".pdf"), width = 8, height = 8, dpi = 200)
      plot.umap.feature(gcdata, features = name, feature.type = "meta", pt.size = 2) #Edited by CBMP
      ggsave(paste0(plotdir, "UMAP_", name, ".pdf"), width = 8, height = 8, dpi = 200) #Edited by CBMP
    }, error=function(e) NA)
  }
  
  # Order clusters by their cluster number. 
  cluster.scores <- cluster.scores %>% dplyr::arrange(as.numeric(as.character(!!as.name(name.cluster))))
  return(cluster.scores)
}

SubclusterSeurat <- function(gcdata, figures.dir, prefix = "", num.sd = 0.6, res = 0.8, test = "wilcox", 
                             batchcorrect = F) {
  #' This function takes a Seurat object and clusters the cells. Cells and genes in this object have 
  #' already been filtered for quality, and the gene expression values have already been normalized 
  #' and log-transformed. This function's purpose is to do subclustering of specific cells taken from 
  #' the entire dataset, in order to find more detailed structure in the clusters. This function 
  #' identifies the variable genes for this cell subset, centers the gene expression, does PCA, 
  #' Louvain clustering, and calculates tSNE visualization. This code works with Seurat v2.0. 
  #' 
  #' @param gcdata Seurat object with gene-cell expression matrix and tSNE data. The gene expression values
  #' must already be normalized within each cell, and log-transformed, ie gcdata@data slot is filled.
  #' In addition, if doing batch corrrection then there must be a column in gcdata@meta.data named group
  #' that gives the batch grouping information by cell.
  #' @param figures.dir Directory where to write figures.
  #' @param prefix Name of subcluster to use in file names.
  #' @param num.sd Number of standard deviations away from background relation between gene expression
  #' and variability to classify a gene as variable.
  #' @param res This parameter is used in Seurat::FindClusters. It is a resolution parameter that sets 
  #' the granularity of the clustering, with increased values leading to a greater number of clusters.
  #' @param test Differential expression test to use in Seurat FindAllMarkers.
  #' @param batchcorrect Whether to run Seurat CCA batch correction. If batch correction is run, the
  #' Seurat object should have a column named group in gcdata@meta.data that specifies which batch
  #' group each cell comes from.
  #' @return List with Seurat object gcdata.sub and cluster gene markers gcdata.markers as the elements.
  
  if (batchcorrect) {
    # Take a Seurat object and divide it into a list of Seurat objects based on batch group they are from.
    gcdata.list <- vector("list", n_distinct(gcdata@meta.data$group))
    group <- unique(gcdata@meta.data$group)
    for (i in seq(group)) {
      cells <- rownames(gcdata@meta.data[gcdata@meta.data$group == group[i], ])
      gcdata.list[[i]] <- SubsetData(gcdata, cells.use = cells)
    }
    gcdata <- gcdata.list
    rm(gcdata.list)
    
    # Chris Smillie variable gene selection where we identify variable genes common across experiments.
    # Merge Seurat objects. Original sample identities are stored in object@meta.data$orig.ident.
    gcdata.merge <- MergeMultipleSeuratObjects(gcdata, project = "vargenes")
    var.genes <- GetCommonVarGenes(gcdata.merge, gcdata.merge@meta.data$group, method="karthik", 
                                   do.plot=T, num_genes=1500, prefix=paste0(figures.dir, prefix), 
                                   do.flatten=T)
    # var.genes <- GetCommonVarGenes(gcdata.merge@raw.data, gcdata.merge@meta.data$group, method="loess", 
    #                                do.plot=T, num_genes=1500, prefix=figures.dir, do.flatten=T)
    
    # Each gene in genes.use (RunCCA) must be present in each dataset in the gcdata list.
    names.list <- lapply(gcdata, function(data) rownames(data@data))  
    var.genes <- Reduce(intersect, c(list(var.genes), names.list))
    
    # Number of variable genes
    print(length(var.genes))
    
    # Within each individual batch, center the gene expression by subtracting from a gene its mean 
    # gene expression across all cells in the individual batch.
    gcdata <- sapply(gcdata, function(data) ScaleData(object = data, do.center = TRUE, do.scale = FALSE, genes.use = var.genes))
    
    # Run multi-set CCA to identify common sources of variation between two batches.
    # Each gene must be present in each dataset, and cells across batches must have unique names.
    # orig.ident slot in gcdata.CCA@meta.data is set by the first field in cell names.
    gcdata.CCA <- RunMultiCCA(gcdata, genes.use = var.genes, num.ccs = 20)
    # save(gcdata.CCA, sampleid, file = Rda.commonvargenes.CCA.path)  # save current progress
    
    # Visualize results of CCA plot CC1 versus CC2 and look at a violin plot of CC1 scores.
    p1 <- DimPlot(gcdata.CCA, reduction.use = "cca", group.by = "orig.ident", pt.size = 0.5, 
                  do.return = TRUE)
    p2 <- VlnPlot(gcdata.CCA, features.plot = "CC1", group.by = "orig.ident", do.return = TRUE)
    plot_grid(p1, p2)
    PrintDim(gcdata.CCA, reduction.type = "cca", dims.print = 1:2, genes.print = 10)
    
    # Identify how much each CCA component captures variation in cells and genes.
    MetageneBicorPlot(gcdata.CCA, grouping.var = "group", dims.eval = 1:15, display.progress = FALSE)
    ggsave(paste0(figures.dir, prefix, "commonvargenes_CCA_metagene_bicor.pdf"), width = 6, height = 6)
    DimHeatmap(object = gcdata.CCA, reduction.type = "cca", cells.use = 500, dim.use = 1:12, do.balanced = T, do.return = T)
    ggsave(paste0(figures.dir, prefix, "commonvargenes_CCA_heatmap.pdf"), width = 6, height = 6)
    
    # Align the CCA subspaces.
    gcdata.CCA <- AlignSubspace(gcdata.CCA, reduction.type = "cca", grouping.var = "group", dims.align = 1:12)
    #save(gcdata.CCA, sampleid, file = Rda.commonvargenes.CCA.path)  # save current progress
    p1 <- VlnPlot(gcdata.CCA, features.plot = "CC1", group.by = "orig.ident", do.return = TRUE)
    p2 <- VlnPlot(gcdata.CCA, features.plot = "CC2", group.by = "orig.ident", do.return = TRUE)
    p3 <- VlnPlot(gcdata.CCA, features.plot = "ACC1", group.by = "orig.ident", do.return = TRUE)
    p4 <- VlnPlot(gcdata.CCA, features.plot = "ACC2", group.by = "orig.ident", do.return = TRUE)
    plot_grid(p1, p2, p3, p4, ncol = 2)
    
    # Cluster and make tSNE after aligning batches. 
    gcdata.CCA <- FindClusters(gcdata.CCA, reduction.type = "cca.aligned", dims.use = 1:12, resolution = res, 
                               print.output = F, save.SNN = T, force.recalc = T)
    gcdata.CCA <- RunTSNE(gcdata.CCA, reduction.use = "cca.aligned", dims.use = 1:12, do.fast = TRUE, seed.use = 1)
    #gcdata.CCA <- RunUMAP(object = gcdata.CCA, reduction.use = "cca.aligned", dims.use = 1:12, seed.use = 1)
    
    # Find cluster biomarkers and store in list of dataframes
    gcdata.markers <- FindAllMarkers(gcdata.CCA, min.pct = 0.1, logfc.threshold = 0.5, test.use = test, 
                                     return.thresh = 0.1/nrow(gcdata.CCA@data), print.bar = T)
    
    return(list("gcdata" = gcdata.CCA, "markers" = gcdata.markers))
  } else {
    # Karthik method for identifying variable genes
    gcdata <- NB.var.genes(gcdata, x.high.cutoff = 10, x.low.cutoff = 0.01, do.text = TRUE, num.sd = num.sd)
    
    # Number of variable genes
    print(length(gcdata@var.genes))
    
    # Mean center gene expression: subtract from a gene its mean gene expression across all cells.
    gcdata <- ScaleData(gcdata, do.center = TRUE, do.scale = FALSE, genes.use = gcdata@var.genes)
    
    # Do PCA using variable genes or using all genes (pc.genes = rownames(object@data))
    gcdata <- RunPCA(gcdata, pc.genes = gcdata@var.genes, pcs.compute = 30, do.print = TRUE, pcs.print = 5, genes.print = 5)
    
    # For each gene (variable and non-variable), score correlation between its
    # expression across cells and its corresponding PCi score. This is done for all
    # PCi (i from 1 to n), so the final matrix dimension is # genes x # PCs.
    gcdata <- ProjectPCA(gcdata, genes.print = 5)
    
    # Jackstraw analysis
    # gcdata <- JackStraw(gcdata, num.replicate = 10, num.pc = 30)
    
    # PCA plots
    VizPCA(gcdata, pcs.use = 1:2)
    print(PCAPlot(gcdata, dim.1 = 1, dim.2 = 2, pt.size = 1, use.full = FALSE, do.return = TRUE))
    ggsave(paste0(figures.dir, prefix, "PCA.pdf"), width = 6, height = 6)
    
    # Identify significant principal components
    PrintPCA(gcdata, pcs.print = 1:5, genes.print = 10, use.full = FALSE)
    
    png(paste0(figures.dir, prefix, "heat.pdf"), width = 1800, height = 1800, res = 300)
    PCHeatmap(gcdata, pc.use = 1:6, cells.use = 100, do.balanced = TRUE)
    dev.off()
    
    print(PCElbowPlot(gcdata))
    ggsave(paste0(figures.dir, prefix, "screeplot.pdf"), width = 6, height = 6)
    
    # print(JackStrawPlot(gcdata, PCs = 1:12))
    # ggsave(paste0(figures.dir, prefix, "jackstraw.pdf"), width = 6, height = 6)
    
    # PC correlation with number of reads and number of expressed genes
    gcdata <- SetIdent(gcdata, ident.use = 1)
    counts <- as.data.frame(as.matrix(gcdata@data))
    PCPlots(counts, gcdata, no.pcs=6, plotType='all', outDirectory=figures.dir, fileType = '.pdf', colors = c("blue"))
    
    gcdata <- FindClusters(gcdata, reduction.type = "pca", dims.use = 1:20, resolution = res, print.output = F, save.SNN = T, force.recalc = T)
    PrintFindClustersParams(gcdata)
    
    gcdata.markers <- FindAllMarkers(gcdata, min.pct = 0.1, logfc.threshold = 0.5, test.use = test, 
                                     return.thresh = 0.1/nrow(gcdata@data), print.bar = T)
    
    # tSNE visualization
    gcdata <- RunTSNE(gcdata, dims.use = 1:20, do.fast = T, reduction.use = "pca", perplexity = 30)
    
    return(list("gcdata" = gcdata, "markers" = gcdata.markers))
  }
}

AddMoleculeInfoSeuratObject <- function(gcdata, h5.file) {
  #' This function adds information from the [molecule_info.h5]
  #' (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/molecule_info) 
  #' file generated by the Cell Ranger count script to the metadata slot of a Seurat object. 
  #' The function assumes that the names of cell barcodes in the molecule_info.h5 file and in the Seurat 
  #' object are identical. The function uses the 
  #' [DropletUtils library]
  #' (http://bioconductor.org/packages/release/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html).
  #' 
  #' @param gcdata Seurat object containing a metadata slot. Cell barocde names in this object should be 
  #' identical to those in the associated molecule_info.h5 file generated by Cell Ranger count.
  #' @param h5.file Path to molecule_info.h5 file associated with the Seurat object.
  #' @return Seurat object with additional metadata. Currently, the number of reads per cell is calculated 
  #' and added into the metadata slot.
  
  # Read molecule_info.h5 into DataFrame object using DropletUtils.
  mol.info <- read10xMolInfo(h5.file)
  
  # Get cell barcode names from Seurat object, and look up statistics for each of these barcodes within
  # molecule_info.h5. 
  cell.names <- rownames(gcdata@meta.data)
  # df.mol.info <- as.data.frame(mol.info$data[mol.info$data$cell %in% cell.names, ]) %>% group_by(cell) %>% 
  #   summarise(nGene=n_distinct(gene), nUMI=n(), nRead=sum(reads))
  df.mol.info <- as.data.frame(mol.info$data[mol.info$data$cell %in% cell.names, ]) %>% group_by(cell) %>%
    summarise(nRead=sum(reads))
  
  # Add the statistics into the metadata slot of the Seurat object.
  # Make sure that the cell barcode names in the metadata slot are in the same order as before the merge. 
  # This is important as otherwise the Seurat object will be unusable.
  gcdata@meta.data <- merge(x = gcdata@meta.data, y = df.mol.info, by.x = "row.names", by.y = "cell") 
  rownames(gcdata@meta.data) <- gcdata@meta.data$Row.names
  gcdata@meta.data$Row.names <- NULL
  if (!all.equal(rownames(gcdata@meta.data), cell.names))
    stop("The cell barcode names in the metadata slot are not the same as before, after adding molecule_info.")
  return(gcdata)
  
  # a <-AddMoleculeInfoSeuratObject(gcdata[[2]], h5.files[2])
  # mol.info <- read10xMolInfo(h5.files[2])
  # sum(a@meta.data$nReads)/sum(mol.info$data$reads)  # fraction of reads in cells
}

RunEmptyDrops <- function(counts, figures.dir, title = "", lower = NULL, bc.rank = 10000) {
  #' This function identifies empty droplets (10x) by comparing the gene expression distribution
  #' of each cell barcode to the gene expression distribution from ambient RNA. 
  #' The ambient RNA distribution is defined as all cell barcodes with total number of UMI below a 
  #' cutoff (usually 100). This should largely be RNA from highly expressed transcripts in lysed cells. 
  #' The function uses the [DropletUtils library]
  #' (http://bioconductor.org/packages/release/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html).
  #' 
  #' @param data Counts matrix of genes x cells read in from 10x Cell Ranger sparse matrix. 
  #' @param figures.dir Directory where to write figures.
  #' @param lower Total UMI cutoff below which cell barcodes are defined as being ambient RNA. 
  #' The default lower nUMI cutoff in emptydrops is 100, but if sequencing depth is higher, we could 
  #' see a higher number of UMI for the ambient RNA. As an alternative, given we know a certain number 
  #' of cells are loaded on each channel, there is a maximum number of expected cells to be recovered. 
  #' Typically 6000 cells are loaded on a channel, so to be conservative, especially when accounting 
  #' for differences in amount of cell RNA across cell types, we can set the cutoff as the barcode 
  #' with rank 10000 when barcodes are ranked by nUMI. If @param lower is specified, that is the value
  #' used for the cutoff. Otherwise, the function looks at the cell barcodes sorted by nUMI.
  #' @param bc.rank This is the rank of the barcodes sorted by nUMI that is used to set the lower nUMI
  #' cutoff for ambient RNA, if @param lower is not specified. The script looks up the nUMI for the 
  #' barcode with this rank and compares it to the default emptydrops cutoff of 100 nUMI, and takes the 
  #' maximum of those two values. The default is set to rank 10000, as we typically load 6000 cells in a
  #' 10x channel.
  #' @param title Sample name to use as title in plots of empty droplets.
  #' @return Names of cells that are not empty droplets are returned.
  #' 
  set.seed(100)  # EmptyDrops is written in R and we want reproducible runs
  
  # Ambient RNA is droplets containing only a barcoded bead and no cell.  I expect 10-100 nUMI in these 
  # empty drops. Below 5 or 10 nUMI, the counts are likely from droplets that contained a cell, but there
  # was more than one sequencing error in the cell barcode, so the molecule was unable to be assigned to 
  # the correct cell barcode. I remove all these instances, so that they do not contribute to the ambient RNA.
  bcerrors <- Matrix::colSums(counts) <= 5  # cell barcode errors should not have more than 5 nUMI
  counts <- counts[, !bcerrors]  # remove cell barcodes that are likely from sequencing errors
  
  # For low-quality samples, there may be fewer than 10000 droplets. Update bc.rank to use all available cells.
  bc.rank <- min(ncol(counts), bc.rank)  
  
  # Determine what to use for the lower nUMI cutoff to define ambient RNA. 
  br.out <- barcodeRanks(counts)
  # Look at cell barcodes ranked by total UMI if a lower nUMI cutoff for ambient RNA not provided.
  if(is.null(lower)) {
    bc.rank.max <- br.out$total[order(br.out$rank)][bc.rank]  # rank cell barcodes by nUMI
    lower <- max(bc.rank.max, 100)  # never have an ambient background cutoff <100 nUMI
  }
  
  # Determine what to use for upper nUMI cutoff above which all cell barcodes are retained.
  # Use either the knee of the distribution calculated by barcodeRanks(), or 5 x lower NUMI cutoff.
  # Occasionnally the calculated knee is very close to the lower cutoff, and then we use 5 x lower instead.
  retain <- max(br.out$knee, 5*lower)
  
  # Identify empty drops using ambient RNA lower nUMU cutoff and retain cutoff.
  e.out <- emptyDrops(counts, lower = lower, retain = retain, ignore = 200)  # always remove cells with < 200 nUMI
  # e.out <- emptyDrops(counts, lower = lower, retain = Inf, ignore = 200)  
  
  is.cell <- e.out$FDR <= 0.01  # T, F, or NA
  cells <- which(is.cell)  # get only cell indices with T (FDR <= 0.01), and remove NA and F
  
  # The empty drops that are removed are neither cells nor the ambient RNA below the nUMI cutoff.
  # It is useful to collect these cell barcodes to see whether they have any interesting features.
  ambient <- which(Matrix::colSums(counts) <= lower)  # cell bcs below nUMI cutoff that make up ambient RNA.
  emptydrops <- setdiff(seq(ncol(counts)), c(cells, ambient))  
  
  # Determine whether more permutations would lower the p-value and potentially make FDR <= 0.01.
  # These are entries where Limited is True and Significant is False. For now, just print the table.
  table.limited <- table(Limited=e.out$Limited, Significant=is.cell)
  if (table.limited[2, 1] > 0) {
    print(table.limited)
  }
  
  # Plot ranked cell barcodes vs nUMI with the lower cutoff overlaid.
  p1 <- ggplot(data = data.frame(rank=br.out$rank, total=br.out$total)) + 
    geom_point(mapping = aes(x = rank, y = total)) + 
    scale_x_log10(breaks = c(1e1, 1e2, 1e3, 1e4, 1e5)) + 
    scale_y_log10(breaks = c(1e1, 1e2, 1e3, 1e4, 5e5)) + geom_hline(yintercept = lower, color = "blue") +
    geom_hline(yintercept = retain, color = "red") + 
    labs(x="barcode rank", y="Total UMI count") + ggtitle(title)
  # geom_text(aes(10, lower, label = "lower", vjust = -1)) +
  # geom_text(aes(10, br.out$knee, label = "knee", vjust = -1)) + 
  # Plot nUMI vs log(Prob) for each cell barcode. Focus on the cell barcodes with low total UMI.
  df.plot <- data.frame(as.data.frame(e.out), color = is.cell) %>% filter(Total <= 5000) %>% drop_na()
  p2 <- ggplot(df.plot) + geom_point(mapping = aes(x = Total, y = -LogProb, color = color), shape = 1) + 
    labs(x="Total UMI count", y="-Log Probability") + 
    scale_color_manual(name="Cell", breaks=c(TRUE, FALSE), values=c("#ff0000", "#000000")) + ggtitle(title)
  # Plot number of genes per cell barcode for each cell barcode marked as being an empty drop.
  nGene <- Matrix::colSums(counts[, emptydrops]>0)
  p3 <- ggplot(data = data.frame(sample = "emptydrops", nGene = nGene), mapping = aes(x = sample, y = nGene)) + 
    geom_violin(fill = "grey", color = "grey30") + geom_boxplot(width = 0.1, outlier.shape = 1) +
    ggtitle(title)
  p <- plot_grid(p1, p2, p3, align = 'h', labels = c('A', 'B', 'C'), nrow = 1)
  save_plot(paste0(figures.dir, 'emptydrops.pdf'), p, base_aspect_ratio = 2.2, base_height = 6)
  
  # return(counts[, cells])
  return(colnames(counts)[cells])
}

MouseHumanHomolog <- function(mouse.human, genes, species) {
  #' This function takes a vector of mouse genes or human genes, and returns the corresponding human or mouse 
  #' homologs respectively.
  #' 
  #' @param mouse.human Data frame of mouse and human genes where each row contains a pair of homologous genes.
  #' The column names are "mouse" and "human".
  #' @param genes A vector of mouse or human genes to return homologs for.
  #' @param species The species of genes in @param genes.
  #' @return ggplot Vector of homologous genes.
  #' 
  if (species == "human") {
    as.character(mouse.human[mouse.human$human %in% genes, "mouse"])
  } else if (species == "mouse") {
    as.character(mouse.human[mouse.human$mouse %in% genes, "human"])
  }
  else {
    warning(paste0(species, " is invalid choice for species."))
  }
}


PlotHighestExpression <- function(counts, file.plot, n = 50, max.percent = NA) {
  #' This function takes a counts matrix and finds the most highly expressed genes. The function plots how much 
  #' each gene contributes to the total expression in each cell. The function also returns a vector of 
  #' the most highly expressed genes.
  #' 
  #' @param counts Counts matrix where rows are genes, columns are cells, and values are counts. Counts may be UMI
  #' or TPM.
  #' @param file.plot Filename where to write figures.
  #' @param n Total number of most highly expressed genes to plot.
  #' @param max.percent Upper limit of a gene's percent of total expression to display in plot. 
  #' This upper limit number ranges from 0 to 100.
  #' @return Vector of n most highly expressed genes.
  
  # Order genes in descending order by their total counts or expression across all cells.
  exp.genes <- Matrix::rowSums(counts)
  exp.genes <- names(exp.genes[order(exp.genes, decreasing = T)])
  
  # Calculate within each cell what fraction of total counts come from the most highly expressed genes.
  percent.gene.totalcounts <- 100 * sweep(as.matrix(counts[exp.genes[1:n], ]), 2, Matrix::colSums(counts), '/') 
  percent.totalcounts <- round(100 * sum(counts[exp.genes[1:n], ]) / sum(counts), digits = 1)
  if (is.na(max.percent)) {
    max.percent <- max(percent.gene.totalcounts)
  }
  
  # Plot fraction of total counts coming from each gene.
  df.plot <- gather(as.data.frame(t(percent.gene.totalcounts)), "gene", "percent", 1:n)
  df.plot$gene <- factor(df.plot$gene, levels = rev(exp.genes[1:n]))
  p <- ggplot(df.plot, mapping = aes(x = gene, y = percent)) + geom_boxplot(outlier.size = 0.5) + labs(y = "% of total counts") + ggtitle(paste0("Top ", n, " genes account for ", percent.totalcounts, "% of total counts")) + coord_flip(ylim = c(0, max.percent))
  ggsave(file.plot, width = 7, height = 10*n/50)
  
  exp.genes[1:n]  
}