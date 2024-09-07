try(library(reticulate,verbose = FALSE, quietly=TRUE),silent=TRUE)
library(data.table,verbose = FALSE, quietly=TRUE)
library(ggplot2,verbose = FALSE, quietly=TRUE)
library(gtools,verbose = FALSE, quietly=TRUE)
library(RColorBrewer,verbose = FALSE, quietly=TRUE)
try(library(simpleCache,verbose = FALSE, quietly=TRUE),silent=TRUE)
try(library(Seurat,verbose = FALSE, quietly=TRUE),silent=TRUE)
try(library(readxl,verbose = FALSE, quietly=TRUE),silent=TRUE)


#Jupyter lab settings
options(repr.matrix.max.rows=40, repr.matrix.max.cols=200)

#plot settings
theme_set(theme_bw(base_size=10))
theme_update(axis.text=element_text(color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white"))
rotate_labels=function(angle=60,vjust=1,hjust=1){return(theme(axis.text.x = element_text(angle = angle, vjust = vjust,hjust=hjust)))}
pdf.options(useDingbats=FALSE)

#use like this: stat_summary(fun.data = give.n,fun.args = c(y=0), geom = "text",size=4)
give.n <- function(x,y = NULL){
  return(data.frame(y = ifelse(is.null(y),0,y), label = paste0("N=",length(x)))) 
}


#Directories
projectDir="/scratch/projects/HTAPP_MBC/"
baseDir=file.path(projectDir,"phase2")
codeDir="/secure/projects/HTAPP_MBC"
metaDir=file.path(codeDir,"meta")
dataDir=file.path(baseDir,"pipeline/alignreads/")
cacheDir=file.path(baseDir,"RCache/")
analysisDir=file.path(baseDir,"analysis/")
statsDir=file.path(analysisDir,"00_stats")
statsSummaryDir=file.path(statsDir,"summary")
tmpDir=file.path(baseDir,"tmp/")
extDir=file.path(projectDir,"ext_data/")
spatialDir="/scratch/projects/HTAPP_MBC/spatial/"
dir.create(path=tmpDir,showWarnings = F,recursive = TRUE)
dir.create(path=statsDir,showWarnings = F,recursive = TRUE)
dir.create(path=analysisDir,showWarnings = F,recursive = TRUE)


#define colors
ct_colors=list()
ct_colors['MBC']='#73d56d'
ct_colors['MBC_stem-like']='#146c18'
ct_colors['MBC_neuronal']='#39a13c'
ct_colors['MBC_chondroid']='#003b00'

ct_colors['Endothelial']='#feb052'
ct_colors['Endothelial_sinusoidal']='#6d0000'
ct_colors['Endothelial_angiogenic']='#dc7014'
ct_colors['Endothelial_vascular']='#aa3700'

ct_colors['Fibroblast']='#ced208'
ct_colors['Chondrocyte']='#515900'
ct_colors['Smooth muscle_vascular']='#748000'
ct_colors['Stellate']='#323400'
ct_colors['Skeletal muscle']='#a0a800'

ct_colors['Adipocytes']='#bb6fc4'
ct_colors['Hepatocyte']='#f3a3f6'
ct_colors['Keratinocyte']='#53065f'
ct_colors['Neuron']='#873b92'

ct_colors['Macrophage']='#99d1fe'
ct_colors['Monocyte']='#387fb9'
ct_colors['Neutrophil']='#003365'
ct_colors['Erythrocyte']='#66a8dd'
ct_colors['Mast']='#1a588e'

ct_colors['B_plasma']='#f86652'
ct_colors['B']='#cf1917'
ct_colors['T']='#fbb2a1'
ct_colors['NK']='#860000'


#receptor colors
rc_colors=list()
rc_colors['HR+/HER2-'] = '#255668'
rc_colors['HR+/HER2+'] = '#008C80'
rc_colors['HR-/HER2+'] = '#61C074'
rc_colors['HR-/HER2-'] = '#EDEF5C'


#tissue colors
ti_colors=list()
#ti_colors['Axilla'] = '#7D1D67'
#ti_colors['Bone'] = '#A02074'
#ti_colors['Brain'] = '#C22B79'
#ti_colors['Breast'] = '#E13C77'
#ti_colors['Chest wall'] = '#F65A6D'
#ti_colors['Liver'] = '#FD806D'
#ti_colors['Lung'] = '#FFA076'
#ti_colors['Neck'] = '#FFBD88'
#ti_colors['Skin'] = '#FFD99F'

ti_colors['Axilla']='#7c0e6f'
ti_colors['Bone']='#c41296'
ti_colors['Brain']='#e86bbb'
ti_colors['Breast']='#ba100e'
ti_colors['Chest wall']='#fc5644'
ti_colors['Liver']='#fc9a90'
ti_colors['Lung']='#c16b18'
ti_colors['Neck']='#fcaf3f'
ti_colors['Skin']='#fcd7a2'



#Annotation
combined_annot_file=file.path(statsDir,"summary/sample_annot.tsv")
if (file.exists(combined_annot_file)){
    message("Loading combined annotation: annot")
    annot=fread(combined_annot_file,quote="")
}
sample_sheet=fread(file.path(codeDir,"meta/sampletracking_HTAPP_MBC_check.csv"))
sample_sheet[,channel_id:=paste0(sampleid,"_",condition,"_",replicate),]
sample_sheet[,channel_id_match:=paste0(sampleid,"-",condition,"-",replicate),]


sample_sheet[,site_group:=site,]
sample_sheet[site%in%c("axilla","axillary lymph node","lymph node", "right axilla","left axillary lymph node","left axilla"),site_group:="axilla/ln"]
sample_sheet[site%in%c("right breast skin punch","left lower abdomen skin punch"),site_group:="skin"]
sample_sheet[site%in%c("breast","left breast"),site_group:="breast"]
sample_sheet[site%in%c("left parietal mass","left parietal"),site_group:="brain"]
sample_sheet[site%in%c("chest wall","left peristernal soft tissue mass"),site_group:="chest"]

#source(file.path(metaDir,"mbc_signatures_ofir.R"))




#Functions

#simplify receptors
simplify_rec=function(annot){
    annot[,ER:=grepl("ER\\+",receptors_biopsy),]
    annot[,PR:=grepl("PR\\+",receptors_biopsy),]
    annot[,HER2:=grepl("HER2\\+",receptors_biopsy),]
    annot[,receptors_biopsy_simpl:=ifelse(grepl("ER\\+|PR\\+",receptors_biopsy)&grepl("HER2\\+",receptors_biopsy),"HR+/HER2+",
                                          ifelse(grepl("ER\\+|PR\\+",receptors_biopsy)&grepl("HER2\\-",receptors_biopsy),"HR+/HER2-",
                                                 ifelse(grepl("ER\\-",receptors_biopsy)&grepl("PR\\-",receptors_biopsy)&grepl("HER2\\+",receptors_biopsy),"HR-/HER2+",
                                                 ifelse(grepl("ER\\-",receptors_biopsy)&grepl("PR\\-",receptors_biopsy)&grepl("HER2\\-",receptors_biopsy),"HR-/HER2-",NA)))),]
return(annot)
}

#Function to create Seurat object and annotate with SingleR
cell_unification=list("NK"=c("NK_cell","NK cells"),"Endothelial"=c("Endothelial cells","Endothelial_cells"),"B"=c("B-cells","B cell","B_cell"),"Macrophage"=c("Macrophage","Macrophages"),
     "Chondrocyte"=c("Chondrocytes"),"Monocyte"=c("Monocyte","Monocytes"),"Neutrophil"=c("Neutrophils"),"Fibroblast"=c("Fibroblasts"),
     "Smooth muscle"=c("Smooth muscle","Smooth_muscle_cells"),"Epithelial"=c("Epithelial_cells","Epithelial cells"))
    

make_seurat_annot = function(cb,min.features=0,ft=TRUE,cellTypes=NULL){
    
    if(!exists("ref.se")){
    suppressMessages(expr=ref.se <- HumanPrimaryCellAtlasData())}
    
    if(!exists("ref.be")){
    suppressMessages(ref.be <- BlueprintEncodeData())}
    
    if (!is.null(cellTypes)){
        ref.se@colData <- ref.se@colData[ref.se@colData$label.main%in%select,]
        ref.se@assays@data$logcounts=ref.se@assays@data$logcounts[,row.names(ref.se@colData)]
        ref.be@colData <- ref.be@colData[ref.be@colData$label.main%in%select,]
        ref.be@assays@data$logcounts<-ref.be@assays@data$logcounts[,row.names(ref.be@colData)]
    }
    
    
    so <- CreateSeuratObject(counts = cb,min.features = min.features, min.cells = 3)
    so <- PercentageFeatureSet(so,pattern = "^MT-",col.name = "percent.mito")
    so <- NormalizeData(object = so)
    so <- FindVariableFeatures(object = so)
    so <- ScaleData(object = so,vars.to.regress = c("nCount_RNA","percent.mito"))
    so <- RunPCA(object = so)
    so <- FindNeighbors(object = so)
    so <- FindClusters(object = so,algorithm = 4)
    so <- RunTSNE(object = so,check_duplicates = FALSE,dims = 1:10)
    so <- RunUMAP(object = so,umap.method = 'umap-learn' , metric = 'correlation',dims=1:10)
    
    print("Seurat done.")
    
    
     #modify se labels to include plasma cells
    se_labels_main=ref.se$label.main
    se_labels_main[grepl("B_cell:Plasma_cell",ref.se$label.fine)] <- "Plasma_cell"
    
    be_labels_main=ref.be$label.main
    be_labels_main[grepl("Plasma cells",ref.be$label.fine)] <- "Plasma_cell"
    
    s_co=SingleR(test = GetAssayData(object = so, slot = "data"),ref= list(ref.se,ref.be), labels = list(se_labels_main,be_labels_main),fine.tune = ft)
    s_co_cl=SingleR(test = GetAssayData(object = so, slot = "data"),clusters = so@active.ident,ref= list(ref.se,ref.be),method = "cluster", labels = list(se_labels_main,be_labels_main),fine.tune = ft)
   
     
    #modify se labels to include plasma cells (original - missing the correct plasma cell annotation)
    #se_labels_main=ref.se$label.main
    #se_labels_main[grepl("B_cell:Plasma_cell",ref.se$label.fine)] <- "Plasma_cell"
    
    #s_co=SingleR(test = cb,ref= list(ref.se,ref.be), labels = list(se_labels_main,ref.be$label.main),fine.tune = ft)
    #s_co_cl=SingleR(test = GetAssayData(object = so, slot = "data"),clusters = so@active.ident,ref= list(ref.se,ref.be),method = "cluster", labels = list(se_labels_main,ref.be$label.main),fine.tune = ft)
    
    s_co_dt=as.data.table(s_co$scores)
    s_co_dt[,labels:=s_co$labels,]
    s_co_dt[,labels_score:=.SD[,get(labels)],by=1:nrow(s_co_dt)]
    s_co_dt[,labels_unif:=ifelse(labels%in%unlist(cell_unification),names(grep(labels,cell_unification,value = TRUE)),labels),by=1:nrow(s_co_dt)]
    s_co_df=as.data.frame(s_co_dt[,c("labels","labels_score","labels_unif")])
    row.names(s_co_df)=row.names(s_co)
    so=AddMetaData(so,s_co_df)

    labels_cl=so@active.ident
    levels(labels_cl)=s_co_cl$labels[order(as.numeric(row.names(s_co_cl)))]
    labels_cl_unif=sapply(as.character(labels_cl),function(x){return(ifelse(x%in%unlist(cell_unification),names(grep(x,cell_unification,value = TRUE)),x))})
    so=AddMetaData(so,as.character(labels_cl),col.name = "labels_cl")
    so=AddMetaData(so,as.character(labels_cl_unif),col.name = "labels_cl_unif")    
    return(so)   
}


subset_seurat = function(so,sel_cells,var.only=FALSE){
    
#    so <- subset(so,cells = sel_cells) # only used the features of the data slot even if counts had all features therefore changed to the line below
    so <- CreateSeuratObject(GetAssayData(so,slot="counts")[,sel_cells],meta.data=so@meta.data[sel_cells,])
    so <- NormalizeData(object = so)
    so <- FindVariableFeatures(object = so)
    if (var.only==TRUE){ so@assays$RNA@data=so@assays$RNA@data[so@assays$RNA@var.features,]}
    so <- ScaleData(object = so,vars.to.regress = c("nCount_RNA","percent.mito"))
    so <- RunPCA(object = so)
    so <- FindNeighbors(object = so)
    so <- FindClusters(object = so,algorithm = 4)
    so <- RunTSNE(object = so,check_duplicates = FALSE,dims = 1:10)
    so <- RunUMAP(object = so,umap.method = 'umap-learn' , metric = 'correlation',dims=1:10)
    
    return(so)
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

make_singleR_annot = function(so,ft=TRUE){
    
    if(!exists("ref.se")){
    suppressMessages(expr=ref.se <- HumanPrimaryCellAtlasData())}
    
    if(!exists("ref.be")){
    suppressMessages(ref.be <- BlueprintEncodeData())}
    
    
    #modify se labels to include plasma cells
    se_labels_main=ref.se$label.main
    se_labels_main[grepl("B_cell:Plasma_cell",ref.se$label.fine)] <- "Plasma_cell"
    
    be_labels_main=ref.be$label.main
    be_labels_main[grepl("Plasma cells",ref.be$label.fine)] <- "Plasma_cell"
    
    s_co=SingleR(test = GetAssayData(object = so, slot = "data"),ref= list(ref.se,ref.be), labels = list(se_labels_main,be_labels_main),fine.tune = ft)
    s_co_cl=SingleR(test = GetAssayData(object = so, slot = "data"),clusters = so@active.ident,ref= list(ref.se,ref.be),method = "cluster", labels = list(se_labels_main,be_labels_main),fine.tune = ft)
    
    s_co_dt=as.data.table(s_co$scores)
    s_co_dt[,labels:=s_co$labels,]
    s_co_dt[,labels_score:=.SD[,get(labels)],by=1:nrow(s_co_dt)]
    s_co_dt[,labels_unif:=ifelse(labels%in%unlist(cell_unification),names(grep(labels,cell_unification,value = TRUE)),labels),by=1:nrow(s_co_dt)]
    s_co_df=as.data.frame(s_co_dt[,c("labels","labels_score","labels_unif")])
    row.names(s_co_df)=row.names(s_co)
    so=AddMetaData(so,s_co_df)

    labels_cl=so@active.ident
    levels(labels_cl)=s_co_cl$labels[order(as.numeric(row.names(s_co_cl)))]
    labels_cl_unif=sapply(as.character(labels_cl),function(x){return(ifelse(x%in%unlist(cell_unification),names(grep(x,cell_unification,value = TRUE)),x))})
    so=AddMetaData(so,as.character(labels_cl),col.name = "labels_cl")
    so=AddMetaData(so,as.character(labels_cl_unif),col.name = "labels_cl_unif")    
    return(so)   
}
