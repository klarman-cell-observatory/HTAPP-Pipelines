import logging as logg
import scanpy as sc
import pandas as pd
import os
try:
    from plotnine import * 
    #plotnine
    theme_set(theme_bw(base_size=10))
    theme_update(axis_text=element_text(color="black"),panel_grid_major = element_blank(), panel_grid_minor = element_blank(),panel_background = element_rect(fill = "white"))
    def rotate_labels (angle=60,vjust=1,hjust=1):
        return theme(axis_text_x = element_text(angle = angle, vjust = vjust,hjust=hjust))
except:
    print("running without plotnine")

### Display settings
pd.set_option('display.max_columns', 100)
sc.settings.set_figure_params(dpi=80)




### Paths

baseDir="/scratch/projects/HTAPP_MBC/phase2/"
codeDir="/secure/projects/HTAPP_MBC"
metaDir=codeDir+"/meta"
dataDir=baseDir+"/pipeline/alignreads/"
out_dir=baseDir+"/PyCache"
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

de_dir=baseDir + "/differential_expression"
if not os.path.isdir(de_dir):
    os.mkdir(de_dir)

### Functions
# adapted from here: https://scanpy.discourse.group/t/heatmap-max-number-of-cells-per-cluster/202/2  
def downsample_to_smallest_category(
        adata,
        column="sample_short",
        random_state=42,
        min_cells=15,
        keep_small_categories=False
) -> sc.AnnData:
    """
    returns an annData object in which all categories in 'column' have
    the same size

    column
        column with the categories to downsample
    min_cells
        Minimum number of cells to downsample.
        Categories having less than `min_cells` are discarded unless
        keep_small_categories is True
    keep_small_categories
        Be default categories with less than min_cells are discarded.
        Set to true to keep them
    """
    counts = adata.obs[column].value_counts(sort=False)
    if len(counts[counts < min_cells]) > 0 and keep_small_categories is False:
        logg.warning(
            "The following categories have less than {} cells and will be "
            "ignored: {}".format(min_cells, dict(counts[counts < min_cells]))
        )
    min_size = min(counts[counts >= min_cells])
    sample_selection = None
    for sample, num_cells in counts.items():
        if num_cells <= min_cells:
            if keep_small_categories:
                sel = adata.obs.index.isin(
                    adata.obs[adata.obs[column] == sample].index)
            else:
                continue
        else:
            sel = adata.obs.index.isin(
                adata.obs[adata.obs[column] == sample]
                .sample(min_size, random_state=random_state)
                .index
            )
        if sample_selection is None:
            sample_selection = sel
        else:
            sample_selection |= sel
    logg.info(
        "The cells in category {!r} had been down-sampled to have each {} cells. "
        "The original counts where {}".format(column, min_size, dict(counts))
    )
    return adata[sample_selection].copy()


def get_sub_markers (adata_sel,types,col='named_cluster_split',rest_thres=0.1,mode="genes"):
    key="temp_sub_markers"
    sub=adata_sel[adata_sel.obs[col].isin(types)].copy()
    sc.tl.rank_genes_groups(sub, col,groups=types, reference='rest',key_added=key, method='t-test_overestim_var',pts=True)
    de_info=pd.concat([sc.get.rank_genes_groups_df(sub,key= key, group=g)[:200].set_index('names').assign(pts=sub.uns[key]['pts'][g],pts_rest=sub.uns[key]['pts_rest'][g],type=g).reset_index() for g in types])
    de_info=de_info.reset_index()
    if (mode == "genes"):
        return(de_info[de_info.pts_rest<rest_thres].groupby('type').apply(lambda x: x['names'][:10]))
    if (mode == "info"):
        return(de_info)



    
    
def get_sub_markers_big (adatas,subset=None,types=None,col='leiden',rest_thres=1,p_thres=0.05,mode="plot"):
    types_orig=types
    key="temp_sub_markers"
    if subset is not None:
        sub=adatas[adatas.obs[subset[0]].isin(subset[1])].copy()
    if types is not None:
        sub=adatas[adatas.obs[col].isin(types)].copy()
    if subset is None and types is None:
        sub=adatas.copy()
        
    types=list(sub.obs[col].dtype.categories)

    sc.tl.rank_genes_groups(sub, col,groups=types, reference='rest',key_added=key, method='t-test_overestim_var',pts=True)
    de_info=pd.concat([sc.get.rank_genes_groups_df(sub,key= key, group=g,pval_cutoff=p_thres)[:200].set_index('names').assign(pts=sub.uns[key]['pts'][g],pts_rest=sub.uns[key]['pts_rest'][g],type=g).reset_index() for g in types])
    de_info=de_info.reset_index()
    de_info["pts_rat"]=de_info.pts_rest/de_info.pts
    de_info=de_info.sort_values("pts_rat")
    if (mode == "genes"):
        return(de_info[de_info.pts_rest<rest_thres].groupby('type').apply(lambda x: x['names'][:10]))
    if (mode == "info"):
        return(de_info)
    if (mode=="plot"):
        markers=de_info[de_info.pts_rest<rest_thres].groupby('type').apply(lambda x: x['names'][:5])
        markers=dict(iter(markers.reset_index().groupby('type')['names']))
        sc.pl.stacked_violin(sub, markers, groupby=col, dendrogram=True,row_palette="pastel")
    types=types_orig

### colors
# basic colors: "#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33"

colors=pd.Series(dtype=object)
colors['MBC']='#73d56d'
colors['MBC_stem-like']='#146c18'
colors['MBC_neuronal']='#39a13c'
colors['MBC_chondroid']='#003b00'

colors['Endothelial']='#feb052'
colors['Endothelial_sinusoidal']='#6d0000'
colors['Endothelial_angiogenic']='#dc7014'
colors['Endothelial_vascular']='#aa3700'

colors['Fibroblast']='#ced208'
colors['Chondrocyte']='#515900'
colors['Smooth muscle_vascular']='#748000'
colors['Stellate']='#323400'
colors['Skeletal muscle']='#a0a800'

colors['Adipocytes']='#bb6fc4'
colors['Hepatocyte']='#f3a3f6'
colors['Keratinocyte']='#53065f'
colors['Neuron']='#873b92'

colors['Macrophage']='#99d1fe'
colors['Monocyte']='#387fb9'
colors['Neutrophil']='#003365'
colors['Erythrocyte']='#66a8dd'
colors['Mast']='#1a588e'

colors['B_plasma']='#f86652'
colors['B']='#cf1917'
colors['T']='#fbb2a1'
colors['NK']='#860000'


#receptor colors
rc_colors=pd.Series(dtype=object)
rc_colors['HR+/HER2-'] = '#255668'
rc_colors['HR+/HER2+'] = '#008C80'
rc_colors['HR-/HER2+'] = '#61C074'
rc_colors['HR-/HER2-'] = '#EDEF5C'


#tissue colors
ti_colors=pd.Series(dtype=object)
ti_colors['Axilla']='#7c0e6f'
ti_colors['Bone']='#c41296'
ti_colors['Brain']='#e86bbb'
ti_colors['Breast']='#ba100e'
ti_colors['Chest wall']='#fc5644'
ti_colors['Liver']='#fc9a90'
ti_colors['Lung']='#c16b18'
ti_colors['Neck']='#fcaf3f'
ti_colors['Skin']='#fcd7a2'