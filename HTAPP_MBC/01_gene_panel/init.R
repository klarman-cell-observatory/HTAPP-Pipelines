library(reticulate)
library(data.table)
library(ggplot2)
library(gtools)
library(RColorBrewer)
library(simpleCache)
library(Seurat)
library(readxl)

#Jupyter lab settings
options(repr.matrix.max.rows=20, repr.matrix.max.cols=200)

#plot settings
theme_set(theme_bw(base_size=10))
theme_update(axis.text=element_text(color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white"))
rotate_labels=function(angle=60,vjust=1,hjust=1){return(theme(axis.text.x = element_text(angle = angle, vjust = vjust,hjust=hjust)))}
pdf.options(useDingbats=FALSE)


#Directories
projectDir="/ahg/regevdata/projects/HTAPP_MBC/"
baseDir=file.path(projectDir,"phase1")
codeDir="/ahg/regevdata/users/jklugham/projects/HTAPP_MBC"
metaDir=file.path(codeDir,"meta")
dataDir=file.path(baseDir,"pipeline/alignreads/")
cacheDir=file.path(baseDir,"RCache/")
analysisDir=file.path(baseDir,"analysis/")
statsDir=file.path(analysisDir,"00_stats")
statsSummaryDir=file.path(statsDir,"summary")
tmpDir=file.path(baseDir,"tmp/")
extDir=file.path(projectDir,"ext_data/")
dir.create(path=tmpDir,showWarnings = F,recursive = TRUE)

#Annotation
combined_annot_file=file.path(statsDir,"summary/sample_annot.tsv")
if (file.exists(combined_annot_file)){
    message("Loading combined annotation: annot")
    annot=fread(combined_annot_file,quote="")
}
sample_sheet=fread(file.path(codeDir,"meta/sampletracking_HTAPP_MBC.csv"))

#Settings
cell_types=c("Astrocyte","Neurons","Neuroepithelial_cell","Embryonic_stem_cells","iPS_cells","Chondrocytes","Osteoblasts","Smooth_muscle_cells","Epithelial_cells","Endothelial_cells","Fibroblasts","Hepatocytes","Macrophage","Monocyte","NK_cell","T_cells","B_cell")
genes=c("ESR1","PGR","ERBB2","EPCAM","KRT8","KRT18","CD74","PTPRC","CD3D","CD19","CD68","ACTB","GAPDH","EIF2A","ALB")
marker_genes=c("ESR1","PGR","ERBB2","KRT8","KRT18","EPCAM","ACTB","PTPRC","CD3D","CD8A","CD4","CD19" ,"CD163","CD68","LYZ","MPO","HLA-DRA","CD74","THY1","PECAM1","VWF","ICAM1","VIM","COL1A1","CYP1A1", "HNF4", "ABCC2", "ALB")


#Functions
#get expression with annotation
get_expression=function(s_obj, genes){
    exp_raw=as.data.table(t(as.matrix(s_obj@raw.data[genes,])),keep.rownames="cellid")
    exp_raw[,type:="raw",]
    exp_norm=as.data.table(t(as.matrix(s_obj@data[genes,])),keep.rownames="cellid")
    exp_norm[,type:="norm",]
    exp_scaled=as.data.table(t(as.matrix(s_obj@scale.data[genes,])),keep.rownames="cellid")
    exp_scaled[,type:="scaled",]
    all_exp=rbindlist(list(exp_raw,exp_norm,exp_scaled))
    
    annot=as.data.table(s_obj@meta.data,keep.rownames="cellid")
    stopifnot(exp_raw$cellid==exp_norm$cellid&exp_norm$cellid==exp_scaled$cellid&exp_scaled$cellid==annot$cellid)
    exp=merge(all_exp,annot,by = "cellid")

    return(exp)
}


#function to add annotations
add_annotations=function(s_obj){
    cond=as.data.table(s_obj@meta.data[,c("sampleID","labels"),drop=FALSE],keep.rownames = "cellid")
    #add new labels column
    cond[,labels_simpl:=ifelse(labels%in%cell_types,labels,"other"),]
    cond[,labels_simpl:=ifelse(labels%in%c("Astrocyte","Neurons","Neuroepithelial_cell"),"Epithelial_neuro",labels_simpl),]
    cond[,labels_simpl:=ifelse(labels%in%c("Embryonic_stem_cells","iPS_cells"),"Epithelial_stem",labels_simpl),]
    cond[,labels_simpl:=ifelse(labels%in%c("Chondrocytes","Osteoblasts","Smooth_muscle_cells"),"Mesenchymal",labels_simpl),]
    cond[,labels_simpl:=ifelse(grepl("B_cell",labels),"B_cell",labels_simpl),]
    cond[,labels_group:=ifelse(grepl("Epithelial",labels_simpl),"Malignant",ifelse(grepl("B_|T_|NK_",labels_simpl),"Lymphocytes",ifelse(grepl("Mono|Macro",labels_simpl),"Mono/Macro",labels_simpl))),]
    #clinical and technical annotations
    cond_an=merge(cond,unique(annot[,c("sampleid","site","receptors_biopsy","site_biopsy","inflammatory"),]),by.x="sampleID",by.y="sampleid")
    cond_an=merge(cond_an,unique(sample_sheet[,c("sampleid","mode","date","flowcell","site"),])[!duplicated(sampleid)],by.x="sampleID",by.y="sampleid")
    cond_an=data.frame(cond_an,row.names = which(colnames(cond_an)=="cellid"))
    s_obj=AddMetaData(s_obj,metadata = cond_an)
    return(s_obj)
}



#function to order genes according to their expression pattern
order_genes=function(obj,sname="labels",svals,ncells=500,genes){
    dat=SubsetData(obj,subset.name = sname,accept.value = svals,max.cells.per.ident = ncells)@scale.data[genes,]
    cl=hclust(dist(dat))
    genes_ordered=data.table(gene=cl$labels[cl$order],order=1:length(cl$order))
    return(genes_ordered[order(order)])  
}


FeaturePlot_annot=function(object,features, ncol,...){
  #calculate cluster centers for annotation in FeaturePlot
  tsne_clusters=merge(as.data.table(object@ident,keep.rownames=TRUE),as.data.table(object@dr$tsne@cell.embeddings,keep.rownames = TRUE),by.x="V1",by.y="rn")
  centers=tsne_clusters[,list(tSNE_1=mean(tSNE_1),tSNE_2=mean(tSNE_2)),by=V2]

  pl=FeaturePlot(object = object, features.plot = features, cols.use = c("lightgrey", "blue"), pt.size = 0.5,nCol = ncol,do.return = TRUE)
  pl_an=lapply(pl,function(x)x+geom_text(data = centers,aes(x=tSNE_1,y=tSNE_2,label=V2)))
  print(x = cowplot::plot_grid(plotlist = pl_an, ncol = ncol))
}


plot_expression=function(object,genes,expr_thres=1,groups=NULL){
  genes=genes[genes%in%row.names(object@data)]
  features_df=FetchData(object,vars.all=genes)
  annot=object@meta.data[,groups,drop=FALSE]
  comb=as.data.table(merge(features_df,annot,by=0))
  comb_long=melt(comb,id.vars=c("Row.names",groups))
  stats=comb_long[,.(Npos=sum(value>expr_thres),Ntotal=.N,frac_pos=sum(value>expr_thres)/.N,max_expr=max(value)),by=c("variable",groups)]
  
  if (!is.null(groups)){  
    pl=ggplot(comb_long,aes(x=variable,fill=get(groups),col=get(groups)))+geom_violin(aes(y=value),scale = "width",bw=0.1)+geom_hline(yintercept = expr_thres,col="red",lty=20)+geom_text(data = stats,position =position_dodge(width=1),aes(y=max_expr+0.5,col=get(groups),label=paste0(signif(frac_pos,2)*100,"%")))+xlab("")+ylab("Normalized expression")
      }else{
        pl=ggplot(comb_long,aes(x=variable))+geom_violin(aes(y=value),scale = "width",bw=0.1,fill="transparent")+geom_hline(yintercept = expr_thres,col="red",lty=20)+geom_text(data = stats, aes(y=max_expr+0.5,label=paste0(signif(frac_pos,2)*100,"%")))+xlab("")+ylab("Normalized expression")     
  }

  return(pl)
}
               
#### modified infercnv v0.8.2
library(infercnv)
infercnv_run=function (infercnv_obj, cutoff = 1, min_cells_per_gene = 3, out_dir = ".", 
    normalize_factor = NA, window_length = 101, num_ref_groups = NULL, 
    max_centered_threshold = NA, noise_filter = NA, sd_amplifier = 1.5, 
    cluster_by_groups = FALSE, k_obs_groups = 1, outlier_method_bound = "average_bound", 
    outlier_lower_bound = NA, outlier_upper_bound = NA, hclust_method = "complete", 
    anscombe_normalize = TRUE, use_zscores = FALSE, remove_genes_at_chr_ends = FALSE, 
    mask_nonDE_genes = FALSE, mask_nonDE_pval = 0.05, test.use = "wilcoxon", 
    plot_steps = FALSE, debug = FALSE, include.spike = FALSE, 
    spike_in_chrs = NULL, spike_in_multiplier_vec = NULL, pseudocount = 0,plot=FALSE) 
{
    if (debug) {
        flog.threshold(DEBUG)
    }
    flog.info(paste("::process_data:Start", sep = ""))
    if (out_dir != "." & !file.exists(out_dir)) {
        dir.create(out_dir)
    }
    step_count = 0
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %d: incoming data\n", step_count))
    if (plot_steps) {
        infercnv_obj_incoming_data <- infercnv_obj
        save("infercnv_obj_incoming_data", file = file.path(out_dir, 
            sprintf("%02d_incoming_data.infercnv_obj", step_count)))
    }
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: Removing lowly expressed genes\n", 
        step_count))
    infercnv_obj <- require_above_min_mean_expr_cutoff(infercnv_obj, 
        cutoff)
    infercnv_obj <- require_above_min_cells_ref(infercnv_obj, 
        min_cells_per_gene = min_cells_per_gene)
    if (plot_steps) {
        infercnv_obj_low_expr_genes_pruned <- infercnv_obj
        save("infercnv_obj_low_expr_genes_pruned", file = file.path(out_dir, 
            sprintf("%02d_reduced_by_cutoff.infercnv_obj", step_count)))
    }
    if (pseudocount != 0) {
        flog.info(sprintf("Adding pseudocount: %g", pseudocount))
        infercnv_obj <- add_pseudocount(infercnv_obj, pseudocount)
    }
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: normalization by sequencing depth\n", 
        step_count))
    infercnv_obj <- normalize_counts_by_seq_depth(infercnv_obj, 
        normalize_factor = normalize_factor)
    if (plot_steps) {
        infercnv_obj_normalize_by_depth <- infercnv_obj
        save("infercnv_obj_normalize_by_depth", file = file.path(out_dir, 
            sprintf("%02d_normalized_by_depth.infercnv_obj", 
                step_count)))
    }
    if (include.spike) {
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: Spiking in genes with variation added for tracking\n", 
            step_count))
        if (!(is.null(spike_in_chrs) && is.null(spike_in_multiplier_vec))) {
            infercnv_obj <- spike_in_variation_chrs(infercnv_obj, 
                spike_in_chrs, spike_in_multiplier_vec)
        }
        else {
            infercnv_obj <- spike_in_variation_chrs(infercnv_obj)
        }
        if (plot_steps) {
            infercnv_obj_spiked <- infercnv_obj
            save("infercnv_obj_spiked", file = file.path(out_dir, 
                sprintf("%02d_spiked.infercnv_obj", step_count)))
            plot_cnv(infercnv_obj = infercnv_obj, k_obs_groups = k_obs_groups, 
                cluster_by_groups = cluster_by_groups, out_dir = out_dir, 
                color_safe_pal = FALSE, x.center = mean(infercnv_obj@expr.data), 
                x.range = "auto", title = sprintf("%02d_spike_added", 
                  step_count), obs_title = "Observations (Cells)", 
                ref_title = "References (Cells)", output_filename = sprintf("infercnv.%02d_spike_added", 
                  step_count), write_expr_matrix = TRUE)
        }
    }
    if (anscombe_normalize) {
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: anscombe normalization\n", 
            step_count))
        infercnv_obj <- anscombe_transform(infercnv_obj)
        if (plot_steps) {
            infercnv_obj_anscombe_norm <- infercnv_obj
            save("infercnv_obj_anscombe_norm", file = file.path(out_dir, 
                sprintf("%02d_anscombe_normalization.infercnv_obj", 
                  step_count)))
            plot_cnv(infercnv_obj = infercnv_obj, k_obs_groups = k_obs_groups, 
                cluster_by_groups = cluster_by_groups, out_dir = out_dir, 
                color_safe_pal = FALSE, x.center = mean(infercnv_obj@expr.data), 
                x.range = "auto", title = sprintf("%02d_anscombe_norm", 
                  step_count), obs_title = "Observations (Cells)", 
                ref_title = "References (Cells)", output_filename = sprintf("infercnv.%02d_anscombe_norm", 
                  step_count), write_expr_matrix = TRUE)
        }
    }
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: log transformation of data\n", 
        step_count))
    infercnv_obj <- log2xplus1(infercnv_obj)
    if (plot_steps) {
        infercnv_obj_log_transformed <- infercnv_obj
        save("infercnv_obj_log_transformed", file = file.path(out_dir, 
            sprintf("%02d_logtransformed.infercnv_obj", step_count)))
        plot_cnv(infercnv_obj = infercnv_obj, k_obs_groups = k_obs_groups, 
            cluster_by_groups = cluster_by_groups, out_dir = out_dir, 
            color_safe_pal = FALSE, x.center = mean(infercnv_obj@expr.data), 
            x.range = "auto", title = sprintf("%02d_log_transformed_data", 
                step_count), obs_title = "Observations (Cells)", 
            ref_title = "References (Cells)", output_filename = sprintf("infercnv.%02d_log_transformed", 
                step_count), write_expr_matrix = TRUE)
    }
    if (use_zscores) {
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: Z-score transformation of data\n", 
            step_count))
        infercnv_obj <- transform_to_reference_based_Zscores(infercnv_obj)
        if (plot_steps) {
            infercnv_obj_zscores <- infercnv_obj
            save("infercnv_obj_zscores", file = file.path(out_dir, 
                sprintf("%02d_Z-scores.infercnv_obj", step_count)))
            plot_cnv(infercnv_obj = infercnv_obj, k_obs_groups = k_obs_groups, 
                cluster_by_groups = cluster_by_groups, out_dir = out_dir, 
                color_safe_pal = FALSE, x.center = 0, x.range = "auto", 
                title = sprintf("%02d_centering_gene_expr", step_count), 
                obs_title = "Observations (Cells)", ref_title = "References (Cells)", 
                output_filename = sprintf("infercnv.%02d_centering_gene_expr", 
                  step_count), write_expr_matrix = TRUE)
        }
    }
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: apply max centered expression threshold\n", 
        step_count))
    threshold = max_centered_threshold
    if (is.na(max_centered_threshold)) {
        threshold = mean(abs(get_average_bounds(infercnv_obj)))
    }
    infercnv_obj <- apply_max_threshold_bounds(infercnv_obj, 
        threshold = threshold)
    if (plot_steps) {
        infercnv_obj_max_centered_expr <- infercnv_obj
        save("infercnv_obj_max_centered_expr", file = file.path(out_dir, 
            sprintf("%02d_apply_max_centered_expr_threshold.infercnv_obj", 
                step_count)))
        plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups, cluster_by_groups = cluster_by_groups, 
            out_dir = out_dir, color_safe_pal = FALSE, x.center = mean(infercnv_obj@expr.data), 
            x.range = "auto", title = sprintf("%02d_apply_max_centered_expr_threshold", 
                step_count), obs_title = "Observations (Cells)", 
            ref_title = "References (Cells)", output_filename = sprintf("infercnv.%02d_apply_max_centred_expr_threshold", 
                step_count), write_expr_matrix = TRUE)
    }
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: Smoothing data per cell by chromosome\n", 
        step_count))
    infercnv_obj <- smooth_by_chromosome(infercnv_obj, window_length = window_length, 
        smooth_ends = TRUE)
    if (plot_steps) {
        infercnv_obj_smoothed_by_chr <- infercnv_obj
        save("infercnv_obj_smoothed_by_chr", file = file.path(out_dir, 
            sprintf("%02d_smoothed_by_chr.infercnv_obj", step_count)))
        plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups, cluster_by_groups = cluster_by_groups, 
            out_dir = out_dir, color_safe_pal = FALSE, x.center = mean(infercnv_obj@expr.data), 
            x.range = "auto", title = sprintf("%02d_smoothed_by_chr", 
                step_count), obs_title = "Observations (Cells)", 
            ref_title = "References (Cells)", output_filename = sprintf("infercnv.%02d_smoothed_by_chr", 
                step_count))
    }
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: re-centering data across chromosome after smoothing\n", 
        step_count))
    infercnv_obj <- center_cell_expr_across_chromosome(infercnv_obj, 
        method = "median")
    if (plot_steps) {
        infercnv_obj_cell_centered <- infercnv_obj
        save("infercnv_obj_cell_centered", file = file.path(out_dir, 
            sprintf("%02d_recentered_cells_by_chr.infercnv_obj", 
                step_count)))
        plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups, cluster_by_groups = cluster_by_groups, 
            out_dir = out_dir, color_safe_pal = FALSE, x.center = mean(infercnv_obj@expr.data), 
            x.range = "auto", title = sprintf("%02d_centering_of_smoothed", 
                step_count), obs_title = "Observations (Cells)", 
            ref_title = "References (Cells)", output_filename = sprintf("infercnv.%02d_centering_of_smoothed", 
                step_count))
    }
    if (!is.null(num_ref_groups)) {
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: splitting reference data into %d clusters\n", 
            step_count, num_ref_groups))
        infercnv_obj <- split_references(infercnv_obj, num_groups = num_ref_groups, 
            hclust_method = hclust_method)
    }
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: removing average of reference data\n", 
        step_count))
    infercnv_obj <- subtract_ref_expr_from_obs(infercnv_obj, 
        inv_log = TRUE)
    if (plot_steps) {
        infercnv_obj_subtract_ref <- infercnv_obj
        save("infercnv_obj_subtract_ref", file = file.path(out_dir, 
            sprintf("%02d_remove_ref_avg_from_obs.infercnv_obj", 
                step_count)))
        plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups, cluster_by_groups = cluster_by_groups, 
            out_dir = out_dir, color_safe_pal = FALSE, x.center = 0, 
            x.range = "auto", title = sprintf("%02d_remove_average", 
                step_count), obs_title = "Observations (Cells)", 
            ref_title = "References (Cells)", output_filename = sprintf("infercnv.%02d_remove_average", 
                step_count))
    }
    if (remove_genes_at_chr_ends == TRUE) {
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: removing genes at chr ends\n", 
            step_count))
        infercnv_obj <- remove_genes_at_ends_of_chromosomes(infercnv_obj, 
            window_length)
        if (plot_steps) {
            infercnv_obj_remove_chr_end_genes <- infercnv_obj
            save("infercnv_obj_remove_chr_end_genes", file = file.path(out_dir, 
                sprintf("%02d_remove_gene_at_chr_ends.infercnv_obj", 
                  step_count)))
            plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups, 
                cluster_by_groups = cluster_by_groups, out_dir = out_dir, 
                color_safe_pal = FALSE, x.center = 0, x.range = "auto", 
                title = sprintf("%02d_remove_genes_at_chr_ends", 
                  step_count), obs_title = "Observations (Cells)", 
                ref_title = "References (Cells)", output_filename = sprintf("infercnv.%02d_remove_genes_at_chr_ends", 
                  step_count), write_expr_matrix = TRUE)
        }
    }
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: invert log2(FC) to FC\n", 
        step_count))
    infercnv_obj <- invert_log2(infercnv_obj)
    if (plot_steps) {
        infercnv_obj_invert_log_transform <- infercnv_obj
        save("infercnv_obj_invert_log_transform", file = file.path(out_dir, 
            sprintf("%02d_invert_log_transform.infercnv_obj", 
                step_count)))
        plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups, cluster_by_groups = cluster_by_groups, 
            out_dir = out_dir, color_safe_pal = FALSE, x.center = 1, 
            x.range = "auto", title = sprintf("%02d_invert_log_transform log(FC)->FC", 
                step_count), obs_title = "Observations (Cells)", 
            ref_title = "References (Cells)", output_filename = sprintf("infercnv.%02d_invert_log_FC", 
                step_count), write_expr_matrix = TRUE)
    }
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: Denoising\n", step_count))
    if (!is.na(noise_filter)) {
        if (noise_filter > 0) {
            flog.info(paste("::process_data:Remove noise, noise threshold at: ", 
                noise_filter))
            infercnv_obj <- clear_noise(infercnv_obj, threshold = noise_filter)
        }
        else {
        }
    }
    else {
        flog.info(paste("::process_data:Remove noise, noise threshold defined via ref mean sd_amplifier: ", 
            sd_amplifier))
        infercnv_obj <- clear_noise_via_ref_mean_sd(infercnv_obj, 
            sd_amplifier = sd_amplifier)
    }
    if (plot_steps) {
        infercnv_obj_denoised <- infercnv_obj
        save("infercnv_obj_denoised", file = file.path(out_dir, 
            sprintf("%02d_denoise.infercnv_obj", step_count)))
        plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups, cluster_by_groups = cluster_by_groups, 
            out_dir = out_dir, color_safe_pal = FALSE, x.center = 1, 
            x.range = "auto", title = sprintf("%02d_denoised", 
                step_count), obs_title = "Observations (Cells)", 
            ref_title = "References (Cells)", output_filename = sprintf("infercnv.%02d_denoised", 
                step_count))
    }
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: Removing outliers\n", 
        step_count))
    infercnv_obj = remove_outliers_norm(infercnv_obj, out_method = outlier_method_bound, 
        lower_bound = outlier_lower_bound, upper_bound = outlier_upper_bound)
    if (plot_steps) {
        infercnv_obj_remove_outliers <- infercnv_obj
        save("infercnv_obj_remove_outliers", file = file.path(out_dir, 
            sprintf("%02d_remove_outlier.infercnv_obj", step_count)))
        plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups, cluster_by_groups = cluster_by_groups, 
            out_dir = out_dir, color_safe_pal = FALSE, x.center = 1, 
            x.range = "auto", title = sprintf("%02d_removed_outliers", 
                step_count), obs_title = "Observations (Cells)", 
            ref_title = "References (Cells)", output_filename = sprintf("infercnv.%02d_removed_outliers", 
                step_count))
    }
    plot_data = infercnv_obj@expr.data
    high_threshold = max(abs(quantile(plot_data[plot_data != 
        0], c(0.05, 0.95))))
    low_threshold = -1 * high_threshold
    if (include.spike) {
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: Scaling according to spike\n", 
            step_count))
        infercnv_obj <- scale_cnv_by_spike(infercnv_obj)
        low_threshold = 0
        high_threshold = 2
        if (plot_steps) {
            infercnv_obj_scaled_by_spike <- infercnv_obj
            save("infercnv_obj_scaled_by_spike", file = file.path(out_dir, 
                sprintf("%02d_scaled_by_spike.infercnv_obj", 
                  step_count)))
            plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups, 
                cluster_by_groups = cluster_by_groups, out_dir = out_dir, 
                color_safe_pal = FALSE, x.center = 1, x.range = c(low_threshold, 
                  high_threshold), title = sprintf("%02d_scaled_by_spike", 
                  step_count), obs_title = "Observations (Cells)", 
                ref_title = "References (Cells)", output_filename = sprintf("infercnv.%02d_scaled_by_spike", 
                  step_count))
        }
    }
    if (mask_nonDE_genes) {
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: Identify and mask non-DE genes\n", 
            step_count))
        infercnv_obj <- mask_non_DE_genes_basic(infercnv_obj, 
            test.use = test.use, center_val = mean(plot_data))
        if (plot_steps) {
            infercnv_obj_mask_nonDE <- infercnv_obj
            save("infercnv_obj_mask_nonDE", file = file.path(out_dir, 
                sprintf("%02d_mask_nonDE.infercnv_obj", step_count)))
            plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups, 
                cluster_by_groups = cluster_by_groups, out_dir = out_dir, 
                color_safe_pal = FALSE, x.center = 1, x.range = c(low_threshold, 
                  high_threshold), title = sprintf("%02d_mask_nonDE", 
                  step_count), obs_title = "Observations (Cells)", 
                ref_title = "References (Cells)", output_filename = sprintf("infercnv.%02d_mask_nonDE", 
                  step_count))
        }
    }
    if (include.spike) {
        infercnv_obj <- remove_spike(infercnv_obj)
    }
    save("infercnv_obj", file = file.path(out_dir, "run.final.infercnv_obj"))
    if (plot==TRUE){
    flog.info("Making the final infercnv heatmap")
    plot_cnv(infercnv_obj, k_obs_groups = k_obs_groups, cluster_by_groups = cluster_by_groups, 
        out_dir = out_dir, color_safe_pal = FALSE, x.center = 1, 
        x.range = c(low_threshold, high_threshold), title = "inferCNV", 
        obs_title = "Observations (Cells)", ref_title = "References (Cells)", 
        output_filename = "infercnv")}
    return(infercnv_obj)
}               
environment(infercnv_run) <- asNamespace('infercnv')               
               
               
#' Read 10X hdf5 file for Cell Ranger V3 data (Read10X_h5 taken from Seurat release 3 branch of github)
#' https://github.com/satijalab/seurat/blob/release/3.0/R/preprocessing.R
#'
#' Read count matrix from 10X CellRanger hdf5 file.
#' This can be used to read both scATAC-seq and scRNA-seq matrices.
#'
#' @param filename Path to h5 file
#' @param use.names Label row names with feature names rather than ID numbers.
#' @param unique.features Make feature names unique (default TRUE)
#'
#' @return Returns a sparse matrix with rows and columns labeled. If multiple
#' genomes are present, returns a list of sparse matrices (one per genome).
#'
#' @export
#'

Read10X_h5_Seurat3 <- function(filename, use.names = TRUE, unique.features = TRUE) {
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
  genomes <- names(x = infile)
  output <- list()
  if (!infile$attr_exists("PYTABLES_FORMAT_VERSION")) {
    # cellranger version 3
    if (use.names) {
      feature_slot <- 'features/name'
    } else {
      feature_slot <- 'features/id'
    }
  } else {
    if (use.names) {
      feature_slot <- 'gene_names'
    } else {
      feature_slot <- 'genes'
    }
  }
  for (genome in genomes) {
    counts <- infile[[paste0(genome, '/data')]]
    indices <- infile[[paste0(genome, '/indices')]]
    indptr <- infile[[paste0(genome, '/indptr')]]
    shp <- infile[[paste0(genome, '/shape')]]
    features <- infile[[paste0(genome, '/', feature_slot)]][]
    barcodes <- infile[[paste0(genome, '/barcodes')]]
    sparse.mat <- sparseMatrix(
      i = indices[] + 1,
      p = indptr[],
      x = as.numeric(x = counts[]),
      dims = shp[],
      giveCsparse = FALSE
    )
    if (unique.features) {
      features <- make.unique(names = features)
    }
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as(object = sparse.mat, Class = 'dgCMatrix')
    # Split v3 multimodal
    if (infile$exists(name = paste0(genome, '/features/feature_type'))) {
      types <- infile[[paste0(genome, '/features/feature_type')]][]
      types.unique <- unique(x = types)
      if (length(x = types.unique) > 1) {
        message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
        sparse.mat <- sapply(
          X = types.unique,
          FUN = function(x) {
            return(sparse.mat[which(x = types == x), ])
          },
          simplify = FALSE,
          USE.NAMES = TRUE
        )
      }
    }
    output[[genome]] <- sparse.mat
  }
  infile$close_all()
  if (length(x = output) == 1) {
    return(output[[genome]])
  } else{
    return(output)
  }
}

#' Read 10X molecule.info file for Cell Ranger 3 V3 data (read10xMolInfo taken from DropletUtils 1.3.12)
#' https://github.com/MarioniLab/DropletUtils/blob/master/R/read10xMolInfo.R
#' @export
#' @importFrom rhdf5 h5read
#' @importFrom S4Vectors DataFrame
read10xMolInfo <- function(sample, barcode.length=NULL, keep.unmapped=FALSE, 
                           get.cell=TRUE, get.umi=TRUE, get.gem=TRUE, get.gene=TRUE, get.reads=TRUE, 
                           version=c("auto", "2", "3"))
  # Utility function to read useful information from a 10X molecule information file.
  #
  # written by Aaron Lun
  # based on code from Jonathan Griffiths
  # created 20 December 2017    
{
  version <- match.arg(version)
  if (version=="auto") {
    available <- rhdf5::h5ls(sample, recursive=FALSE)
    version <- if ("barcode_idx" %in% available$name) "3" else "2"
  }
  
  data <- list()
  
  if (get.cell) {
    if (version=="3") {
      all.barcodes <- as.vector(rhdf5::h5read(sample, "/barcodes"))
      all.barcodes <- sub("-[0-9]+", "", all.barcodes) # removing GEM group.
      data$cell <- all.barcodes[as.vector(rhdf5::h5read(sample, "/barcode_idx")) + 1L]
    } else {
      data$cell <- .Call(cxx_get_cell_barcodes, sample, "barcode", barcode.length)
    }
  }
  
  if (get.umi) {
    data$umi <- as.vector(rhdf5::h5read(sample, "/umi"))
  }
  
  if (get.gem) {
    data$gem_group <- as.vector(rhdf5::h5read(sample, "/gem_group"))
  }
  
  if (get.gene || !keep.unmapped) {
    # Both of these are zero-indexed by default
    if (version=="3") {
      data$gene <- as.vector(rhdf5::h5read(sample, "/feature_idx")) + 1L 
    } else {
      data$gene <- as.vector(rhdf5::h5read(sample, "/gene")) + 1L
    }
  }
  
  if (get.reads || length(data)==0) {
    if (version=="3") {
      nreads <- as.vector(rhdf5::h5read(sample, "/count"))
    } else {
      nreads <- as.vector(h5read(rhdf5::sample, "/reads")) 
    }
    
    if (get.reads) {
      data$reads <- nreads
    } else {
      # Just to ensure we get the right number of rows,
      # if there were no other fields requested.
      data$reads <- matrix(0L, length(nreads), 0) 
    }
  }
  
  data <- do.call(DataFrame, data)
  
  # Defining the set of all genes, removing unassigned gene entries.
  if (version=="3") {
    gene.ids <- rhdf5::h5read(sample, "/features/id") 
  } else {
    gene.ids <- rhdf5::h5read(sample, "/gene_ids") 
  }
  if (!keep.unmapped) {
    keep <- data$gene <= length(gene.ids)
    if (!get.gene) {
      data$gene <- NULL
    }
    data <- data[keep,]
  }
  
  # Don't define the total cell pool here, as higher level functions may want to use gem_group.
  return(list(data=data, genes=gene.ids))
}


               
#modified SingleR makeSeurat function to only keep highly variable genes
               
SingleR.CreateSeurat_mod <- function(project.name,sc.data,min.genes = 200,
                                 min.cells = 2,regress.out = 'nUMI',
                                 npca = 10,resolution=0.8,temp.dir=NULL) {
  mtgenes = '^mt-'
  
  if (packageVersion('Seurat')>=3) {
    sc = CreateSeuratObject(sc.data, min.cells = min.cells, 
                            min.features = min.genes, project = project.name)
    percent.mito <- PercentageFeatureSet(object = sc, pattern = "^(?i)mt-")
    # mito.features <- grep(pattern = mtgenes, x = rownames(x = sc), value = TRUE,ignore.case=TRUE)
    #  percent.mito <- Matrix::colSums(x = GetAssayData(object = sc, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = sc, slot = 'counts'))
    sc <- AddMetaData(object = sc, metadata = percent.mito, 
                      col.name = "percent.mito")
  } else {
    sc = CreateSeuratObject(sc.data, min.cells = min.cells, 
                            min.genes = min.genes, project = project.name)
    mito.genes <- grep(pattern = mtgenes, x = rownames(x = sc@data), 
                       value = TRUE,ignore.case=TRUE)
    percent.mito <- colSums((sc.data[mito.genes, ]))/colSums(sc.data)
    sc <- AddMetaData(object = sc, metadata = percent.mito, 
                      col.name = "percent.mito")
    
    sc <- NormalizeData(object = sc, 
                        normalization.method = "LogNormalize", 
                        scale.factor = 10000)
  }
  
  if (packageVersion('Seurat')>=3) {
    sc <- SCTransform(object = sc, vars.to.regress = "percent.mito", verbose = FALSE,
                      do.correct.umi=T)
    
    sc <- RunPCA(object = sc,verbose = FALSE)
    sc <- FindNeighbors(object = sc, dims = 1:30)
    sc <- FindClusters(object = sc)
    if (ncol(sc@assays$RNA@data)<100) {
      sc <- RunTSNE(sc,perplexity=10,dims = 1:npca)
    } else {
      sc <- RunTSNE(sc,dims = 1:30)
    }
    sc <- RunUMAP(sc,dims = 1:30, verbose = FALSE)
  } else {
    sc <- FindVariableGenes(object = sc, mean.function = ExpMean, 
                            dispersion.function = LogVMR, 
                            x.low.cutoff = 0.0125, x.high.cutoff = 3, 
                            y.cutoff = 0.5, do.contour = F, do.plot = F)
      #added by JK
      sc@data=sc@data[sc@var.genes,]
      
    if (!is.null(regress.out)) {
      sc <- ScaleData(object = sc, vars.to.regress = regress.out)
    } else {
      sc <- ScaleData(object = sc)
    }
    sc <- RunPCA(object = sc, pc.genes = sc@var.genes, do.print = FALSE)
    #PCElbowPlot(object = sc)
    sc <- FindClusters(object = sc, reduction.type = "pca", 
                       dims.use = 1:npca,resolution = resolution, 
                       print.output = 0, save.SNN = F, 
                       temp.file.location = temp.dir)
    if (ncol(sc@data)<100) {
      sc <- RunTSNE(sc, dims.use = 1:npca, do.fast = T,perplexity=10  )
    } else {
      sc <- RunTSNE(sc, dims.use = 1:npca, do.fast = T,check_duplicates = FALSE)
      
    }
  }
  
  sc
}