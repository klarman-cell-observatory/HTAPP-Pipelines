# load packages
library(infercnv)
library(hdf5r)
library(Matrix) # must be loaded before SparseM to enable cast from SparseM to Matrix
library(SparseM)

# define function to read h5ad (on github.com/KlughammerLab/anndata2r)
read_matrix = function(file, name) {
if (name %in% list.datasets(file,recursive=FALSE)) {
print('dense')
newX = t(file$open(name)$read())
} else if ('X' %in% list.groups(file,recursive=FALSE)) {
groupX = openGroup(file, name)
groupX_encoding_type = h5attr(groupX,'encoding-type')
if (groupX_encoding_type == 'csr_matrix' || groupX_encoding_type == 'csc_matrix' ) {
ra = groupX$open('data')$read()
ja = as.integer(groupX$open('indices')$read()+1)
ia = as.integer(groupX$open('indptr')$read()+1)
dim = h5attr(groupX,'shape')
if (groupX_encoding_type == 'csr_matrix') {
print('csr')
newX = new("matrix.csr", ra=ra, ja=ja, ia=ia, dimension=dim)
} else if (groupX_encoding_type == 'csc_matrix') {
print('csc')
newX = new("matrix.csc", ra=ra, ja=ja, ia=ia, dimension=dim)
}
newX = as(newX,'dgCMatrix')
} else {
print('unkown encoding for X...')
}
} else {
print('unkown encoding for X...')
}
return(newX)
}

read_df = function(file, name) {
group = openGroup(file, name)
categories = NULL
if (group$exists('__categories')) {
categories = group$open('__categories')
categories_ds = list.datasets(categories)
}
df = data.frame(row.names=group$open(h5attr(group, '_index'))$read())
if (length(list.datasets(group)) > 1) { # catch case with only the index and no additional column
for (col in h5attr(group, 'column-order')) {
temp = group$open(col)$read()
if (!is.null(categories) && col %in% categories_ds) {
temp = categories$open(col)$read()[temp+1]
}
df[col] = temp
}
}
return(df)
}


read_adata = function(filename, transpose=FALSE) {
file = H5File$new(filename, mode = "r")
newX = read_matrix(file, 'X')
obs_df = read_df(file, 'obs')
var_df = read_df(file, 'var')
colnames(newX) = row.names(var_df)
row.names(newX) = row.names(obs_df)
if (transpose) {
return(list('X'=t(newX),'obs'=var_df,'var'=obs_df))
} else {
return(list('X'=newX,'obs'=obs_df,'var'=var_df))
}
}


# load h5ad
adata_file='frozen_selected_for_inferCNV.h5ad'
adata = read_adata(adata_file)

# quick first look on input
head(adata$var)
head(adata$obs)
head(adata$X[,1:10])
length(colnames(adata$X))



# create filtered dataframes from input and adjust rownnmes and colnames to match the required input format for inferCNV
head(adata$X)
anno_X <-as.data.frame(t(adata$X))
head(anno_X)
anno_file <- as.data.frame(cbind(rownames(adata$obs),adata$obs$celltype_by_patient)) # change if you want to access other column to distinguish e.g. patients when looking at tumours 
head(anno_file)
gene_o_file <- as.data.frame(cbind(adata$var$ensembl_gene_id,adata$var$chromosome,adata$var$start,adata$var$end))
gene_o_file[,3] <- as.numeric(as.character(gene_o_file[,3]))
gene_o_file[,4] <- as.numeric(as.character(gene_o_file[,4]))
colnames(gene_o_file) <- NULL
head(gene_o_file)
unique(adata$var$chromosome)


# check if rownames of the count matrix and indices of var are identical, i.e.  of same order and lentgh. Then assign gene IDs as new rownames for count matrix  
unique(rownames(t(adata$X))==rownames(adata$var)) # order of rownames seems ok, so transferring of gene names should be too
length(rownames(anno_X))
length(adata$var$ensembl_gene_id)
rownames(anno_X) <- adata$var$ensembl_gene_id

# delete adata here to save memory
rm(adata)
gc()

# further adjustment of dataframes to match required input format
colnames(anno_file) <- NULL
colnames(gene_o_file) <- NULL

head(gene_o_file)
rownames(gene_o_file) <- gene_o_file[,1]
gene_o_file[,1] <-NULL
head(gene_o_file)         

head(anno_file)
rownames(anno_file) <- anno_file[,1]
anno_file[,1] <-NULL
head(anno_file) 

# the final prepared input dataframes
head(gene_o_file)
head(anno_file)
head(anno_X[,1:10])

# create inferCNV object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=anno_X,
                                    annotations_file=anno_file,
                                    delim="\t",
                                    gene_order_file=gene_o_file,
                                    ref_group_names=c('T', 'Macrophage', 'Fibroblast', 'Endothelial_sinusoidal'),
                                    chr_exclude = c('chrGL000009.2', 'chrGL000194.1', 'chrGL000195.1', 'chrGL000213.1', 'chrGL000218.1', 'chrGL000219.1', 'chrKI270711.1', 'chrKI270713.1', 'chrKI270721.1', 'chrKI270726.1','chrKI270727.1', 'chrKI270728.1', 'chrKI270731.1', 'chrKI270734.1', 'chrnan', 'chrMT', 'chrY')
                                   ) 
any(is.na(infercnv_obj@count.data))
any(!is.finite(infercnv_obj@count.data))
any(is.na(infercnv_obj@expr.data))
any(!is.finite(infercnv_obj@expr.data))
lengths(infercnv_obj@reference_grouped_cell_indices)
ouput_dir_path <- "HTAPP_frozen_by_patient"
dir.create(ouput_dir_path)

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=ouput_dir_path, 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=FALSE, # otherwise not enough RAM or other problems that lead to crashes
			     useRaster=FALSE) # to avoid weird RAM error. Set to TRUE to speed up

new_gene_order = data.frame()
for (chr_name in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")) {
	new_gene_order = rbind(new_gene_order, infercnv_obj@gene_order[which(infercnv_obj@gene_order[["chr"]] == chr_name) , , drop=FALSE])
}
names(new_gene_order) <- c("chr", "start", "stop")
infercnv_obj@gene_order = new_gene_order
infercnv_obj@expr.data = infercnv_obj@expr.data[rownames(new_gene_order), , drop=FALSE]

plot_cnv(infercnv_obj,out_dir = ouput_dir_path,useRaster=FALSE)