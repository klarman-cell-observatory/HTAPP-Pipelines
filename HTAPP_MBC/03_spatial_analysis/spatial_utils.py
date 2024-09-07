import pandas as pd
import scanpy as sc
import anndata
import os
import matplotlib.image as mpimg
import numpy as np
import re
from sklearn.cluster import DBSCAN
from skimage import transform
from skimage import color
from skimage import filters
import matplotlib.pyplot as plt
import matplotlib
import sys
import seaborn as sns
import PIL
import glob
import tacco as tc
import tacco

try:
    sys.path.insert(1, '/secure/projects/tacco/tacco')
    import spatial_colocalization
    tacco.spatial_colocalization = spatial_colocalization
except:
    print("Could not import spatial_colocalization")
    
#allows to load large images (otherwise decopression bomb warning)
PIL.Image.MAX_IMAGE_PIXELS=None


#paths
baseDir="/scratch/projects/HTAPP_MBC/phase2/"
projectDir="/scratch/projects/HTAPP_MBC/spatial/"
codeDir="/secure/projects/HTAPP_MBC"
data_dir=projectDir+"/coordination/"
out_dir_base=projectDir+"/01_typing/"
anno_dir=projectDir+"he_annot"

#plotting
pd.set_option('display.max_columns', 500)
sc.settings.set_figure_params(dpi=80)
high_dpi=300
low_dpi=80
matplotlib.rcParams['figure.dpi'] = high_dpi

# scaling
#MERFISH: 1 px is 1um --> 10 um --> 10 px
#Slide-seq: 3mm = 5000 px -->1 px = 0.6um --> 10um = 16.67px #NOTE: 1 px is actually 0.65um --> adjust upon next rerun
#CODEX: 1 px is 386 nm = 0.386um --> 10um = 25.91 px
#ExSeq: 1 px is 51 nm = 0.051um --> 10um = 196.07 px
microns={"slide_seq": 5000/3000,"exseq": 1/0.051,"exseq_bin": 1/0.051,"merfish":1,"merfish_bin":1,"codex":1/0.386}

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

#Mappigs
compartments=pd.Series(dtype=object)
compartments["Malignant"]=['MBC','MBC_stem-like','MBC_neuronal','MBC_chondroid']
compartments["Stromal"]=['Endothelial','Endothelial_sinusoidal','Endothelial_angiogenic','Endothelial_vascular','Fibroblast','Chondrocyte','Smooth muscle_vascular','Stellate','Skeletal muscle','Adipocytes','Hepatocyte','Keratinocyte','Neuron']
compartments["Myeloid"]=['Macrophage','Monocyte','Neutrophil','Erythrocyte','Mast']
compartments["Lymphoid"]=['B_plasma','B','T','NK']

#functions
################# Image registration #######################################
def rotate(coords, origin=None,shift=(0, 0), degrees=0,flip=None):
    coords=coords.copy()
    if origin is None:
        origin=(coords.x_orig.mean(),coords.y_orig.mean())
    
    if flip == "v":
        coords.x_orig =   abs(coords.x_orig - max(coords.x_orig))
    if flip == "h":
        coords.y_orig =   abs(coords.y_orig - max(coords.y_orig))
    
    angle = np.deg2rad(degrees)
    R = np.array([[np.cos(angle), -np.sin(angle)],
                  [np.sin(angle),  np.cos(angle)]])    
    coords-=origin
    coords@= R.T
    coords+=origin
    coords+=shift
    
    coords.rename(columns={0:"x",1:"y"},inplace=True)
    return coords

def scale(im, rat=0.5,crop_x=(0,0),crop_y=(0,0), fill_value=1.0):
    im=im[max(int(crop_y[0]),0):len(im)-int(crop_y[1]),max(int(crop_x[0]),0):len(im[0])-int(crop_x[1])]
    if (crop_x[0]<0) or (crop_y[0]<0):
        if crop_x[0]<0:
            add_left=abs(crop_x[0])
        else:
            add_left=0
        if crop_y[0]<0:
            add_bottom=abs(crop_y[0])
        else:
            add_bottom=0
        fill_value = fill_value * im.max()
        new=np.full(fill_value=fill_value,shape=(im.shape[0]+add_bottom,im.shape[1]+add_left,im.shape[2]),dtype=im.dtype)
        ##new=np.full_like(im,im.max(),shape=(im.shape[0]+add_bottom,im.shape[1]+add_left,im.shape[2]))
        ##new=np.full(fill_value=im.max(),shape=(im.shape[0]+add_bottom,im.shape[1]+add_left,im.shape[2]))
        #new=np.zeros(shape=(im.shape[0]+add_bottom,im.shape[1]+add_left,im.shape[2]),dtype=im.dtype)
        ##new[:,:,:3] = 0
        ##new[:,:,3:] = im.max()
        new[add_bottom:add_bottom+im.shape[0],add_left:add_left+im.shape[1]]=im
        im=new
    im=transform.rescale(im,rat,multichannel=True)
    return im

#alternative way of scaling
#    nR0 = len(im)     # source number of rows 
#    nC0 = len(im[0])  # source number of columns 
#    nR = int(nR0 * rat)
#    nC = int(nC0 * rat)
#    return im[np.arange(0, nR0, nR0 / nR).astype(int)][:,np.arange(0, nC0, nC0 / nC).astype(int)]
    return im
    

def aligner(adatas,img_idx,registration,dat_idxs,img_idx2=None,figsize=(12,12),dot_sizes=4.5,he_method="HE_highres",adjust=False,change=False,adjust_HE=True):
    if img_idx2 is not None: 
        fig, axs = plt.subplots(2,len(dat_idxs)+1,figsize=figsize,squeeze=False, gridspec_kw={'wspace':0,'hspace':0.3})
    else:
        fig, axs = plt.subplots(1,len(dat_idxs),figsize=figsize,squeeze=False)
    img=replicate_splitter(adatas[[he_method]])[img_idx]
    if change==False:
        downsampling = 20
    else:
        downsampling = 1
    img = img[::downsampling,::downsampling,:].copy()
    
    if img_idx2 is not None:
        img2=replicate_splitter(adatas[[he_method]])[img_idx2]
        img2 = img2[::downsampling,::downsampling,:].copy()
    
    if adjust==True and adjust_HE==True:
        img = transform.rotate(img,angle=registration.loc[img_idx,"degrees"],resize=True,mode="edge")
        img=scale(img,rat=registration.loc[img_idx,"scale"],
                      crop_x=(registration.loc[img_idx,"crop_left"]//downsampling,registration.loc[img_idx,"crop_right"]//downsampling),
                      crop_y=(registration.loc[img_idx,"crop_bottom"]//downsampling,registration.loc[img_idx,"crop_top"]//downsampling))
        img=np.clip(img,0,1)
        if img_idx2 is not None:
            img2 = transform.rotate(img2,angle=registration.loc[img_idx2,"degrees"],resize=True,mode="edge")
            img2 = scale(img2,rat=registration.loc[img_idx2,"scale"],
                      crop_x=(registration.loc[img_idx2,"crop_left"]//downsampling,registration.loc[img_idx2,"crop_right"]//downsampling),
                      crop_y=(registration.loc[img_idx2,"crop_bottom"]//downsampling,registration.loc[img_idx2,"crop_top"]//downsampling))
            img2=np.clip(img2,0,1)
        if change == True:
            adatas[he_method][img_idx]=img
            if img_idx2 is not None:
                adatas[he_method][img_idx2]=img2
    
    if img_idx2 is not None:
        from matplotlib.colors import rgb_to_hsv as hsv
        
        def greenify(img):
            maxv = 255 if np.issubdtype(img.dtype, np.integer) else 1
            green = np.zeros((*img.shape[:-1],4),dtype=img.dtype)
            green[...,1] = maxv
            green[...,3] = maxv - hsv(img)[...,2]
            return green
            
        
        green = greenify(img)
        green2 = greenify(img2)
        
        axs[0,0].imshow(img,origin='lower')
        axs[0,0].imshow(green2,origin='lower')
        axs[0,0].grid(False)
        axs[0,0].set_title(img_idx2 +" on "+ img_idx)
        axs[0,0].set_ylim(0, max(img.shape[0],img2.shape[0]))
        axs[0,0].set_xlim(0, max(img.shape[1],img2.shape[1]))
        axs[1,0].imshow(img2,origin='lower')
        axs[1,0].imshow(green,origin='lower')
        axs[1,0].grid(False)
        axs[1,0].set_title(img_idx +" on "+ img_idx2)
        axs[1,0].set_ylim(0, max(img.shape[0],img2.shape[0]))
        axs[1,0].set_xlim(0, max(img.shape[1],img2.shape[1]))
    
    for i,dat_idx in enumerate(dat_idxs):
        iorig=i
        if img_idx2 is not None:
            i=i+1
        dat_obs=replicate_splitter(adatas[[re.sub("[\s\d]*","",dat_idx)]])[dat_idx].obs
        
        if ("x" in dat_obs.columns) & ("y" in dat_obs.columns):
            dat=dat_obs[["x","y"]]
        else:
            dat=dat_obs[["x_orig","y_orig"]]
        
        
        if "region" in dat_obs.columns and adjust==True and len(dat_obs["region"].unique())>1:
            regions=dat_obs["region"].unique().tolist()
            dat=pd.DataFrame()
            for region in regions:
                dat_reg=dat_obs.loc[dat_obs.region==region,["x_orig","y_orig"]]
                dat_reg=rotate(dat_reg,shift=(eval(str(registration.loc[dat_idx,"shift_x"]))[region-1],eval(str(registration.loc[dat_idx,"shift_y"]))[region-1]),degrees=eval(str(registration.loc[dat_idx,"degrees"]))[region-1])
                dat=dat.append(dat_reg)
                if change == True:
                    adatas[re.sub("[\s\d]*","",dat_idx)].obs.loc[(adatas[re.sub("[\s\d]*","",dat_idx)].obs["replicate"]==re.sub("[^\d]*","",dat_idx)) & (adatas[re.sub("[\s\d]*","",dat_idx)].obs["region"]==region),"x"]=dat_reg.x
                    adatas[re.sub("[\s\d]*","",dat_idx)].obs.loc[(adatas[re.sub("[\s\d]*","",dat_idx)].obs["replicate"]==re.sub("[^\d]*","",dat_idx)) & (adatas[re.sub("[\s\d]*","",dat_idx)].obs["region"]==region),"y"]=dat_reg.y
     
        else:
            if adjust==True:
                dat=dat_obs[["x_orig","y_orig"]]
                if "flip" in registration.columns:
                    flip = registration.loc[dat_idx,"flip"]
                else:
                    flip=None
                dat=rotate(dat,shift=(registration.loc[dat_idx,"shift_x"],registration.loc[dat_idx,"shift_y"]),degrees=registration.loc[dat_idx,"degrees"],flip=flip)
                if change == True:
                    adatas[re.sub("[\s\d]*","",dat_idx)].obs.loc[adatas[re.sub("[\s\d]*","",dat_idx)].obs["replicate"]==re.sub("[^\d]*","",dat_idx),"x"]=dat.x
                    adatas[re.sub("[\s\d]*","",dat_idx)].obs.loc[adatas[re.sub("[\s\d]*","",dat_idx)].obs["replicate"]==re.sub("[^\d]*","",dat_idx),"y"]=dat.y
            else:
                if (not "x" in dat_obs.columns) | (not "y" in dat_obs.columns):
                    dat.rename(columns={"x_orig":"x","y_orig":"y"},inplace=True)
        
        
        dat_counts=replicate_splitter(adatas[[re.sub("[\s\d]*","",dat_idx)]])[dat_idx].obs["n_counts"]
        dat_counts /= dat_counts.max()
        colors = [ (0,1,0,a) for a in dat_counts ]
        
        if (type(dot_sizes) is list) and (len(dot_sizes)>1):
            dot_size=dot_sizes[iorig]
        else:
            dot_size=dot_sizes
        
        axs[0,i].imshow(img,origin='lower')
        axs[0,i].scatter(x=dat.x//downsampling,y=dat.y//downsampling, s=dot_size, color=colors)
        axs[0,i].grid(False)
        axs[0,i].set_title(dat_idx +" on "+ img_idx)
        if img_idx2 is not None:
            axs[1,i].imshow(img2,origin='lower')
            axs[1,i].scatter(x=dat.x//downsampling,y=dat.y//downsampling, s=dot_size, color=colors)
            axs[1,i].grid(False)
            axs[1,i].set_title(dat_idx +" on "+ img_idx2)
            
            
def create_mask(adatas, img_idxs,figsize,dat_idxs=["slide_seq"],he_method="HE_highres",dot_sizes=[0.5],s=2,t=0.8,long=True):
    if long == True:
        fig, axs = plt.subplots(len(dat_idxs),len(img_idxs),figsize=figsize,squeeze=False)
    else:
        fig, axs = plt.subplots(1,len(img_idxs)*len(dat_idxs),figsize=figsize,squeeze=False)
    count=0
    for i,img_idx in enumerate(img_idxs):
        img=replicate_splitter(adatas[[he_method]])[img_idx]
        # blur and grayscale before thresholding
        blur = color.rgb2gray(img)
        blur = filters.gaussian(blur, sigma=s)
        # perform inverse binary thresholding
        mask = blur < t
        # use the mask to select the "interesting" part of the image
        sel = np.zeros_like(img)
        sel[mask] = img[mask]
        for j,dat_idx in enumerate(dat_idxs):
            mask.shape
            ijs = adatas[dat_idx].obs[['x','y']].astype(int).to_numpy()
            in_mask = ((ijs < mask.T.shape) & (ijs >= 0)).all(axis=1)
            adatas[dat_idx].obs["ut_"+re.sub(" ","",img_idx)]=in_mask
            adatas[dat_idx].obs.loc[in_mask,"ut_"+re.sub(" ","",img_idx)]=mask[adatas[dat_idx].obs.y.astype(int)[in_mask],adatas[dat_idx].obs.x.astype(int)[in_mask]]
            
            if len(dot_sizes)==1:
                dot_size=dot_sizes[0]
            else:
                dot_size=dot_sizes[j]

            if long ==True:
                loc=i
            else:
                loc=count
                j=0
            
            # display the result
            axs[j,loc].imshow(sel,origin='lower')
    
            plotxy_t=adatas[dat_idx][adatas[dat_idx].obs["ut_"+re.sub(" ","",img_idx)]==True].obs[["x","y"]]
            plotxy_f=adatas[dat_idx][adatas[dat_idx].obs["ut_"+re.sub(" ","",img_idx)]==False].obs[["x","y"]]
            axs[j,loc].scatter(x=plotxy_t.x,y=plotxy_t.y, s=dot_size,alpha=0.2, color="green")
            axs[j,loc].scatter(x=plotxy_f.x,y=plotxy_f.y, s=dot_size,alpha=0.2, color="red")
            axs[j,loc].set_facecolor("black")
            axs[j,loc].set_title(dat_idx+" on "+img_idx)
            axs[j,loc].grid(False)
            count+=1
            
            

###########Adding annotations##############################
def add_annotations(adatas,sample,dat_idxs,registration):
    for replicate in [1,2]:
        replicate=str(replicate)
        sample_ID=sample + "[_HE_|-]*" + replicate
        img=anno_dir+'/'+sample+'/**/'+ sample_ID + "_[a-zA-Z]*_[Ff]ill.png"
        source_glob=glob.glob(img, recursive=True)
        for anno_file in source_glob:
            anno_name=anno_file.split("_")[-2]
            print(anno_name)
            img=mpimg.imread(anno_file)
            
            #adjust annotation HE
            img_idx="HE "+replicate+" annot"
            img = transform.rotate(img,angle=registration.loc[img_idx,"degrees"],resize=True,mode="edge")
            img=scale(img,rat=registration.loc[img_idx,"scale"],
                      crop_x=(registration.loc[img_idx,"crop_left"],registration.loc[img_idx,"crop_right"]),
                      crop_y=(registration.loc[img_idx,"crop_bottom"],registration.loc[img_idx,"crop_top"]),
                      fill_value=0)

            blur = img[:,:,3]
            mask = blur > 0
            
            for j,dat_idx in enumerate(dat_idxs):
                mask.shape
                ijs = adatas[dat_idx].obs[['x','y']].astype(int).to_numpy()
                in_mask = ((ijs < mask.T.shape) & (ijs >= 0)).all(axis=1)
                adatas[dat_idx].obs[anno_name+"_"+replicate]=in_mask
                adatas[dat_idx].obs.loc[in_mask,anno_name+"_"+replicate]=mask[adatas[dat_idx].obs.y.astype(int)[in_mask],adatas[dat_idx].obs.x.astype(int)[in_mask]]
                

def plot_annnotation(adatas,dat_idxs,annots,HE,figsize=[15,12],s=0.1):
    fig, axs= plt.subplots(len(dat_idxs),len(annots),figsize=figsize,sharex=True,sharey=True)
    for i,dat_idx in enumerate(dat_idxs):
        for j,annot in enumerate(annots):
            selection = adatas[dat_idx].obs.loc[adatas[dat_idx].obs[annot]==True]
            axs[i,j].imshow(adatas['HE_highres'][HE],origin='lower')
            axs[i,j].scatter(selection['x'],selection['y'],s=s,c='g')
            axs[i,j].grid("off")
            axs[i,j].set_title(dat_idx +"|"+annot)
            

################ Plotting ##################################        
def plot_spatial_obs(adatas,methods,axsize=(1.5,5),wspace=0.6,hspace=0.2,width_ratios=[3, 0.2],x_shifts=[0, -0.05],xlim=(None,None),ylim=(None,None)):
    n_rep=0
    for method in methods:
        replicates = adatas[method].obs.replicate.unique()
        n_rep += len(replicates)
    fig, axs = tc.pl.subplots(n_rep*2,2,axsize=axsize,wspace=wspace,hspace=hspace,width_ratios=width_ratios*n_rep,x_shifts=x_shifts*n_rep)
    counter=0
    for j,method in enumerate(methods):
        replicates = adatas[method].obs.replicate.unique()
        n_rep = len(replicates)
        for i,rep in enumerate(replicates):
            i=i+counter
            sub_data = adatas[method]
            sub_data = sub_data[sub_data.obs.replicate==rep]
            sc.pl.scatter(adatas[method][adatas[method].obs.replicate==rep],x="x",y="y",color="leiden",show=False, ax=axs[0,i*2]).set_aspect("equal")
            axs[0,i*2].set_title(method + ': replicate ' + rep)
            axs[0,i*2].grid(False)
            axs[0,i*2].set_xlabel('')
            axs[0,i*2].set_ylabel('')
            axs[0,i*2].set_xlim(xlim)
            axs[0,i*2].set_ylim(ylim)
            axs[0,i*2+1].set_axis_off()
            im=axs[1,i*2].scatter(x=sub_data.obs['x'],y=sub_data.obs['y'],s=0.5,alpha=1,c=np.minimum(np.log(sub_data.obs['n_genes'])/np.log(10),3),cmap='viridis')
            axs[1,i*2].set_aspect("equal")
            axs[1,i*2].set_title(method + ': replicate ' + rep)
            axs[1,i*2].set_xlim(xlim)
            axs[1,i*2].set_ylim(ylim)
            axs[1,i*2].grid(False)
            fig.colorbar(im, cax=axs[1,i*2+1])
            [left,bottom,width,height] = axs[1,i*2+1].get_position().bounds
            axs[1,i*2+1].set_aspect(10)   
        counter+=n_rep
        
def plot_spatial_genes(adatas,methods,genes,axsize=(1.5,5),wspace=0.4,hspace=0.2,width_ratios=[3, 0.2],x_shifts=[0, -0.05],s=0.5):
    for method in methods:
        n_genes = len(genes)
        replicates = adatas[method].obs.replicate.unique()
        n_rep = len(replicates)
        fig, axs = tc.pl.subplots(n_genes*2,n_rep,axsize=axsize,wspace=wspace,hspace=hspace,width_ratios=width_ratios*n_genes,x_shifts=x_shifts*n_genes)
        for r,rep in enumerate(replicates):
            for g in range(n_genes):
                sub_data = adatas[method]
                sub_data = sub_data[sub_data.obs.replicate==rep].raw.to_adata()
                im = axs[r,g*2].scatter(x=sub_data.obs['x'],y=sub_data.obs['y'],s=s,alpha=1,c=np.minimum(np.array(sub_data[:,genes[g]].X.sum(axis=1)).flatten(),6),cmap='coolwarm')
                axs[r,g*2].set_aspect("equal")
                axs[r,g*2].set_title('%s: %s %s' % (genes[g],method, rep))
                fig.colorbar(im, cax=axs[r,g*2+1])
                [left,bottom,width,height] = axs[r,g*2+1].get_position().bounds
                axs[r,g*2+1].set_aspect(10)
                axs[r,g*2].grid(False)
    return(fig)
        
def plot_variable_grid(adatas,methods,variables,gridspec_kw={'wspace':0.9,'hspace':0.4},width=6,height=4):
    n_methods=len(methods)
    n_variables=len(variables)
    fig,axs=plt.subplots(n_variables,n_methods,figsize=(width*n_methods,height*n_variables),gridspec_kw=gridspec_kw)
    for i,method in enumerate(methods):
        for j,variable in enumerate(variables): 
            sc.pl.umap(adatas[method], color=[variable], ax=axs[j,i], show=False,use_raw=True);
            axs[j,i].set_title(method+':'+variable)
            
def plot_variable_individual(adatas,method_vars,gridspec_kw={'wspace':0.9,'hspace':0.4}):
    n_method_vars=len(method_vars)
    fig,axs=plt.subplots(1,n_method_vars,figsize=(6*n_method_vars,4),gridspec_kw=gridspec_kw)
    for i,method_var in enumerate(method_vars):
        sc.pl.umap(adatas[method_var[0]], color=[method_var[1]], ax=axs[i], show=False);
        axs[i].set_title(method_var[0]+':'+method_var[1])            
            

def get_sub_markers (adatas,methods,subset=None,types=None,col='leiden',rest_thres=1,p_thres=0.05,mode="plot"):
    types_orig=types
    for method in methods:
        key="temp_sub_markers"
        if subset is not None:
            sub=adatas[method][adatas[method].obs[subset[0]].isin(subset[1])].copy()
        if types is not None:
            sub=adatas[method][adatas[method].obs[col].isin(types)].copy()
        if subset is None and types is None:
            sub=adatas[method].copy()
            
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
            sc.pl.stacked_violin(sub, markers,title=method, groupby=col, dendrogram=True,row_palette="pastel")
        types=types_orig

def cluster_ct_composition(adatas,methods,typings,sel_colors,ct_column="cell_type",split_replicates=False,legend=False):    
    all_weights=pd.DataFrame()
    for method in methods:
        sub=adatas[method]
        for typing in typings:
            combined=pd.DataFrame(sub.obs[["leiden","replicate"]])
            weights=sub.obsm[typing]
            weights=weights/weights.sum(axis=1).to_numpy()[:,None]
            combined=combined.merge(weights,left_index=True,right_index=True)
            combined_long=pd.melt(combined,id_vars=["leiden","replicate"],value_name="weight",var_name=ct_column)
            combined_long["method"]=method
            combined_long["typing"]=typing
            if split_replicates:
                combined_long["sample"]=combined_long["method"]+" "+combined_long["replicate"].astype(str)
            else:
                combined_long["sample"]=combined_long["method"]
            all_weights=all_weights.append(combined_long)
            
    if split_replicates:
        mean_weights=all_weights.groupby(['leiden',ct_column,'replicate','method','typing','sample']).mean('weight').reset_index()
    else:
        mean_weights=all_weights.groupby(['leiden',ct_column,'method','typing','sample']).mean('weight').reset_index()

    mean_weights['leiden']=mean_weights['leiden'].astype(pd.CategoricalDtype(categories=pd.Index(mean_weights['leiden'].unique()).astype(int).sort_values().astype(str),ordered=True))
    mean_weights[ct_column]=mean_weights[ct_column].astype(pd.CategoricalDtype(categories=sel_colors.index,ordered=True))
    
    mean_weights=mean_weights.sort_values(by=['leiden',ct_column,'sample',"typing"])
    
    fig=sns.relplot(data=mean_weights,x='leiden',y=ct_column,hue=ct_column,size="weight",
                hue_order=sel_colors.index.tolist(),
                height=3.5, aspect=1.8, sizes=(5,250),
                palette=sel_colors.values.tolist(),col='sample',row='typing',facet_kws={"sharex":False,"sharey":True},legend=legend);
    return(fig.fig)


def plot_confusion(adata,keys_1,keys_2,sel_colors,sort_cols=True, sort_rows=True,sharex=True,sharey=True,axsize=(8,3),vmax=1000):
    fig,axs=tc.pl.subplots(n_x=len(keys_1),n_y=len(keys_2),axsize=axsize,sharex=sharex,sharey=sharey)
    for k1,key1 in enumerate(keys_1):
        for k2,key2 in enumerate(keys_2):
            transitions = tc.tl.dataframe2anndata(adata.obs, obs_key=key1, var_key=key2)
            transitions = pd.DataFrame(transitions.X.A, index=transitions.obs.index, columns=transitions.var.index)
            if sort_cols:
                transitions=transitions.reindex(columns=sel_colors.index)
            if sort_rows:
                transitions=transitions.reindex(index=sel_colors.index)
            transitions=transitions.dropna(axis=0,how='all').dropna(axis=1,how='all')
            g=sns.heatmap(transitions.astype(int),xticklabels=True,yticklabels=True,annot=True,vmax=vmax,ax=axs[k2,k1],annot_kws={'size':8},fmt='d')
            g.set_xticklabels(g.get_xticklabels(), rotation = 45, ha='right')
            g.title.set_text(key2)
    return(fig)


############# Processing #######################################        
        
# function to process single-cell or soatial data with standard scanpy pipeline
def create_scanpy(adatas,replicates=None,var_genes=None,batch_key='replicate',redo=False,process=True,mode=None):
    if redo == True:
        if "counts" in list(adatas.obsm.keys()):
            adata_reset=anndata.AnnData(X=adatas.obsm["counts"],obs=adatas.obs,var=pd.DataFrame(index=adatas.uns["counts_var"]))
        else:
            adata_reset=adatas
   
        if ("ut_HE1" in adata_reset.obs.columns.to_list()) and ("ut_HE2" in adata_reset.obs.columns.to_list()):
            #TODO: let each replicate choose its own HE 
            covered={'ut_HE1':adata_reset.obs['ut_HE1'].sum(),'ut_HE2':adata_reset.obs['ut_HE2'].sum()}
            adata_reset=adata_reset[adata_reset.obs[max(covered, key=covered.get)]==True]
        elif "ut_HE1" in adata_reset.obs.columns.to_list():
            adata_reset=adata_reset[adata_reset.obs["ut_HE1"]==True]
        elif "ut_HE2" in adata_reset.obs.columns.to_list():
            adata_reset=adata_reset[adata_reset.obs["ut_HE2"]==True]
        else:
            print("No HE masks found. Not trimming data.")
        adata=adata_reset
    else:
        if replicates is None: 
            adata=adatas[0].copy()
        else:
            adata=adatas[0].concatenate(adatas[1:],batch_categories=replicates,batch_key = batch_key,join="outer")
        
    
    if mode !="codex":
        if (mode == "slide_seq") or (mode=="scRNAseq"):
            min_counts=30
            min_genes=30
            sc.pp.filter_cells(adata, min_counts=min_counts)
            sc.pp.filter_cells(adata, min_genes=min_genes)
            sc.pp.filter_genes(adata, min_cells=3)
    
            #make sure percentage of super low beads is below 0.35 (if higher clustering etc goes crazy)
            all_counts=adata.obs.n_counts
            cf=min_counts
            while (all_counts[all_counts>cf]<100).value_counts(normalize=True).loc[False] <=0.65:
                cf+=1
            sc.pp.filter_cells(adata, min_counts=cf)
    
        if mode == "merfish" or mode == "merfish_bin" or mode == "exseq" or mode == "exseq_bin":
            min_counts=20
            min_genes=1
            sc.pp.filter_cells(adata, min_counts=min_counts)
            sc.pp.filter_cells(adata, min_genes=min_genes)
            sc.pp.filter_genes(adata, min_cells=3)
            sc.pp.filter_cells(adata, min_counts=min_counts)
        
        if process == False:
            return adata
    
        counts = adata.X.copy()
        var = adata.var.index.to_numpy()
    
        try:
            sc.pp.highly_variable_genes(adata,flavor='seurat_v3',n_top_genes=5000, batch_key=batch_key)
        except:
            print("Batch aware HVG selection failed. Running without batch_key.")
            sc.pp.highly_variable_genes(adata,flavor='seurat_v3',n_top_genes=5000)
    
        adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata
        adata.obsm["counts"] = counts
        adata.uns["counts_var"] = var
        if (not var_genes is None):
            adata.var['highly_variable']= adata.var.index.isin(var_genes)
        
        adata = adata[:, adata.var['highly_variable']]
        sc.pp.regress_out(adata, ['n_counts', 'pct_counts_mt']);
    
    else:
        min_counts=1
        min_genes=1
        sc.pp.filter_cells(adata, min_counts=min_counts)
        sc.pp.filter_cells(adata, min_genes=min_genes)
        if process == False:
            return adata
        print("Normalization turned off.")
        counts = adata.X.copy()
        var = adata.var.index.to_numpy()
        sc.pp.log1p(adata)
        adata.raw = adata
        adata.obsm["counts"] = counts
        adata.uns["counts_var"] = var
#        sc.pp.regress_out(adata, ['n_counts']);
    
    
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack',random_state=42)
    sc.pp.neighbors(adata,random_state=42)
    sc.tl.leiden(adata,resolution=0.8,random_state=42)
    sc.tl.umap(adata,random_state=42)
    return adata

#load, process and cache all spatial and sc data
def get_and_process_data(sample,data_dir,out_dir,methods,scale_he=0.4,sc_method="scRNAseq",scratch=False,redo=False,process=True,no_write=False,exclude=(None,None),verbose=True):
    adatas=pd.Series(index=methods,dtype=object) 
    for method in methods:
        fname = out_dir+"/"+sample+"_"+method+"_processed.h5ad"
        if (os.path.isfile(fname)) and (scratch==False):
            
            if redo==False:
                print('reading ', method, ' from buffer') if verbose else None
                adatas[method]=sc.read_h5ad(fname)
            elif redo == True:
                print('recreating',method) if verbose else None
                adatas[method]=create_scanpy(sc.read_h5ad(fname),redo=True,process=process,mode=method)
                if no_write==False:
                    adatas[method].write(fname)
        else:
            if method == sc_method:
                print("Creating "+method) if verbose else None
                sc_raw=sc.read_h5ad(data_dir+"/"+sample+"/"+sample+"_"+method+".h5ad")
                sc_raw.obs["replicate"]=[re.sub("[a-zA-Z]*","",i) for i in sc_raw.obs["replicate"].astype(str)]
                adatas[method]=create_scanpy([sc_raw],process=process,mode=method)
                if no_write==False:
                    adatas[method].write(fname)
            elif "HE" in method:
                adatas[method]=pd.Series(dtype="object")
                for replicate in range(1,4):
                    replicate=str(replicate)
                    sample_ID=sample + "_" + replicate
                    img_file=data_dir+"/"+sample+"/"+sample_ID+"_"+method+".jpg"
                    img_file_processed=out_dir+"/"+sample+"_HE_"+replicate+"_processed.jpg"
                    if os.path.isfile(img_file_processed) and (scratch==False):
                        print("Loading processed "+method) if verbose else None
                        adatas[method]["HE "+replicate]=mpimg.imread(img_file_processed)
                    else:
                        if os.path.isfile(img_file):
                            print("Creating "+method) if verbose else None
                            if scale_he == 1:
                                adatas[method]["HE "+replicate]=mpimg.imread(img_file)
                            else:
                                mpimg.thumbnail(img_file,fname.replace("h5ad","jpg"),scale=scale_he)
                                adatas[method]["HE "+replicate]=mpimg.imread(fname.replace("h5ad","jpg"))
                                os.remove(fname.replace("h5ad","jpg"))
            else:
                rows=[]
                for replicate in range(1,4):
                    if (method==exclude[0]) and (replicate==exclude[1]):
                        next
                    replicate=str(replicate)
                    sample_ID=sample + "_" + replicate
                    file=data_dir+"/"+sample+"/"+sample_ID+"_"+method+".h5ad"
                    if os.path.isfile(file):
                        adata=sc.read_h5ad(file)
                        adata.obs.rename(columns={"xcoord": "x_orig", "ycoord": "y_orig"},inplace=True)
                        adata.obs["x_orig"]=(adata.obs["x_orig"]-min(adata.obs["x_orig"]))/microns[method]
                        adata.obs["y_orig"]=(adata.obs["y_orig"]-min(adata.obs["y_orig"]))/microns[method]
                        if method=="slide_seq":
                            clustering = DBSCAN(eps=50, min_samples=20).fit(adata.obs[["x_orig","y_orig"]])
                            adata=adata[np.array(clustering.labels_)!=-1]
                            adata.obs["x_orig"]=adata.obs["x_orig"]-min(adata.obs["x_orig"])
                            adata.obs["y_orig"]=adata.obs["y_orig"]-min(adata.obs["y_orig"])
                        rows.append((replicate,adata))
                if len(rows)>0:
                    print("Creating "+method) if verbose else None
                    sp_raws=pd.DataFrame(rows,columns=('replicate','data'))
                    adatas[method]=create_scanpy(sp_raws.data,sp_raws.replicate,process=process,mode=method)
                    if no_write==False:
                        adatas[method].write(fname)
                else:
                    print("No data for "+method) if verbose else None
    return(adatas.dropna())



def update_processed_adatas(adatas,sample,out_dir,methods,overwrite=False,process=False,exclude=(None,None)):
    for method in methods:
        if method not in adatas.index:
            print(method + " not available.")
            continue
        if "HE" in method:
            for img,idx in zip(adatas[method],adatas[method].index):
                mpimg.imsave(out_dir+"/"+sample+"_"+idx.replace(" ","_")+"_processed.jpg",img)
            
        else:
            fname = out_dir+"/"+sample+"_"+method+"_processed.h5ad"
            if (os.path.isfile(fname) & overwrite==True)|(not os.path.isfile(fname)):
                if process == True:
                    if (method==exclude[0]):
                        adatas[method]=adatas[method][adatas[method].obs['replicate']!=str(exclude[1])]
                    print("Processing " + method)
                    adatas[method]=create_scanpy(adatas[method],redo=True,process=True,mode=method)
                    adatas[method].write(fname)
                else:
                    print("Writing " + method)
                    adatas[method].write(fname)
            else:
                print("Exists and overwirte disabled. Not writing")
    return adatas
        
        


########### Cell typing ##############################

def run_typing(adatas,out_dir,sample,ct_column="cell_type",sc_method="scRNAseq",kick_genes=None, threshold=None,rerun=False, RCTD=True):
    def subset_genes(adata, kick_genes):
        if (kick_genes is None) and (threshold is None):
            return adata
        if isinstance(kick_genes, str):
            kick_genes = [kick_genes]
        if kick_genes is None:
            kick_genes=[]
        bad_genes = np.isin(adata.uns['counts_var'],kick_genes)
        filtered = adata.copy()
        filtered.uns['counts_var'] = filtered.uns['counts_var'][~bad_genes]
        filtered.obsm['counts'] = filtered.obsm['counts'][:,~bad_genes]
        
        if threshold is not None:
            filtered.obsm['counts'].data = np.maximum(filtered.obsm['counts'].data-threshold, 0)
        
        return filtered
        
    pipes = pd.Series(dtype=object)
    pipes['OT'] = lambda sp_raw: tc.tl.annotate(sp_raw, reference=adatas[sc_method], annotation_key=ct_column,method='OT',platform_iterations=0,normalize_to='reference',lamb=0.001,multi_center=4 ,bisections=4, bisection_divisor=3,counts_location=("obsm","counts"))
    pipes['OT_max'] = lambda sp_raw: tc.tl.annotate(subset_genes(sp_raw,kick_genes), reference=adatas[sc_method], annotation_key=ct_column,method='OT',max_annotation=1,platform_iterations=0,normalize_to='reference', lamb=0.001, multi_center=4,bisections=4, bisection_divisor=3,counts_location=("obsm","counts"))
    if RCTD is True:
        pipes['RCTD'] = lambda sp_raw: tc.tl.annotate(sp_raw, reference=adatas[sc_method], annotation_key=ct_column,method='RCTD',conda_env='/ahg/regevdata/users/smages/conda_envs/RCTD', platform_iterations=None, doublet=False, n_cores=6, min_ct=2, verbose=False, counts_location=("obsm","counts"))
    
    
    rows = []
    
    repl_idxs=[i for i in list(replicate_splitter(adatas).index) if not (('scRNAseq' in i)|('HE' in i))]
    for method, sample_ID, sp_raw in zip([re.sub("[\s\d]*","",i) for i in repl_idxs],[sample+"_"+re.sub("[^\d]*","",i) for i in repl_idxs],replicate_splitter(adatas)[repl_idxs]):
        print("++"+method+"++")
        reses = pipes.copy()
        for key, pipe in pipes.items():
            fname = out_dir+"/"+sample_ID+"_"+method+"_"+key+".csv"
            if (os.path.isfile(fname)) and (rerun==False):
                print('reading ', key, ' from buffer')
                res = pd.read_csv(fname)
                res.set_index(res.columns[0], inplace=True)
                res.index.rename(None, inplace=True)
            else:
                print('running ', key)
                try:
                    res = pipe(sp_raw)
                    res.index.rename(None, inplace=True)
                    res.to_csv(fname)
                except Exception as e:
                    res=None
                    print(e)
                    
            if res is not None:
                reses[key] = res
            
        rows.append(reses)
        
    sp_typing = pd.DataFrame(rows,columns=pipes.index,index=repl_idxs)
    return sp_typing

def add_typing_results(adatas,sp_typing):
    for idx,run in sp_typing.iterrows():
        method, replicate = run.name.split(" ")
        for typing in sp_typing.columns:
            try:
                data=sp_typing.loc[run.name,typing].copy()
            except Exception as e:
                continue
            if not typing in list(adatas[method].obsm.keys()):
                adatas[method].obsm[typing]=data.reindex(adatas[method].obs.index)
            else:
                indexintersection = data.index.intersection(adatas[method].obs.index)
                adatas[method].obsm[typing].loc[indexintersection]= data.loc[indexintersection]
            
##########Analysis#####################################
from scipy.cluster.hierarchy import linkage, leaves_list
def plot_spatial_coloc(adatas, methods,typings,sel_colors,cluster=False,npermut=100,radius=20):
    cm_neighbors = 50
    methodXreplicates= replicate_splitter(adatas[methods]).index
    morans={}
    pvals={}
    for methodXreplicate in methodXreplicates:
        sel=replicate_splitter(adatas)[methodXreplicate].copy()
        for typing in typings:
#            print(methodXreplicate +" "+ typing)
            freqs=sel.obsm[typing].dropna(1,how='all')
            freqs=freqs.loc[:, (freqs != 0).any(axis=0)]            
            freqs.fillna(0,inplace=True)
            cross_morans,pvals_tmp=tacco.spatial_colocalization.cross_moran_pvals_shuffle_cells(sel.obs[['x', 'y']], freqs.T, freqs.T, radius=radius, verbose=False, npermut=npermut, n_neighbors=cm_neighbors, exclude_self=False)
#            cross_morans = tacco.spatial_colocalization.cluster_dist_df(cross_morans) 
            
            morans[methodXreplicate+"|"+typing]=cross_morans
            pvals[methodXreplicate+"|"+typing]=pvals_tmp
            
    fig,axs=tc.pl.subplots(len(methodXreplicates),len(typings),hspace=0.7,wspace=0.6,axsize=(5,4))
    for moran,ax in zip(morans.keys(),axs.flatten()):
        _moran = morans[moran].reindex(columns=sel_colors.index,index=sel_colors.index)
        _pval = pvals[moran].reindex(columns=sel_colors.index,index=sel_colors.index)
        if cluster==True:
            reordering = leaves_list(linkage(_moran.fillna(-1))) # TODO: maybe a common clustered ordering for all heatmaps would be better for the comparison between
            _moran=_moran.iloc[reordering,reordering]
            _pval=_pval.iloc[reordering,reordering]
        sns.heatmap(_moran,annot=_pval,center=0,ax=ax)
        ax.set_title(moran)
    return fig
               
########### Convenience ###############################
#find well covered squares across all sections
def get_squares(adatas,stride=10,size=1000):
    
    if 'exseq' in adatas.index:
        sel=adatas['exseq'].copy()
    else:
        sel=adatas['slide_seq'].copy()
    max_x=sel.obs.x.max()
    max_y=sel.obs.y.max()
    
    x_min=0
    x_max=min(size,max_x-1)
    y_min=0
    y_max=min(size,max_y-1)
    
    res={"x_min":[],"x_max":[],"y_min":[],"y_max":[],"count":[]}
    while (x_max<max_x):
        while (y_max<max_y):
            count=len(sel.obs.loc[(sel.obs['x']>=x_min)&(sel.obs['x']<=x_max)&(sel.obs['y']>=y_min)&(sel.obs['y']<=y_max),:])
            res['x_min'].append(x_min)
            res['y_min'].append(y_min)
            res['x_max'].append(x_max)
            res['y_max'].append(y_max)
            res['count'].append(count)
            
            y_min=y_min+stride
            y_max=y_max+stride
            
        x_min=x_min+stride
        x_max=x_max+stride
        y_min=0
        y_max=size
    res=pd.DataFrame(res)
    res=res.sort_values('count',ascending=False)
    for method in adatas.index:
        if method=="scRNAseq":
            continue
        adatas[method].uns['squares_'+str(size)]=res.loc[res['count']>100,:]
    return adatas



#select regions to plot
def square(adatas,square='squares_1000',coords_loc=0,coords=None):
    for method in adatas.index:
        sel=adatas[method]
        if coords == None:
            x_min,x_max,y_min,y_max = sel.uns[square].iloc[coords_loc][['x_min','x_max','y_min','y_max']]
        else:
            x_min,x_max,y_min,y_max=coords
        adatas[method]=sel[(sel.obs['x']>=x_min)&(sel.obs['x']<=x_max)&(sel.obs['y']>=y_min)&(sel.obs['y']<=y_max)]
    return adatas


def replicate_splitter(adatas,repl_col="replicate",he_method="HE_highres"):
    indices=[]
    #get indices
    for key in adatas.keys():
        sel=adatas[key]
        if key in [he_method]:
            indices.append(sel.index.tolist())
        else:
            pre = key
            pre += ' % s'
            indices.append([pre % i for i in sel.obs[repl_col].unique().astype(str)])
    
    indices = [item for sublist in indices for item in sublist]
    split=pd.Series(index=indices,dtype=object)
    for key in adatas.keys():
        sel=adatas[key]
        if key in [he_method]:
            for key2 in sel.index.tolist():
                split[key2]=sel[key2]
        else:
            for repl in sel.obs[repl_col].unique():
                key2=key+" "+repl
                split[key2]=sel[sel.obs[repl_col]==repl].copy()
    return(split)



#generate combined single cell data
def generate_combined_reference(data_dir=data_dir,ct_column="cell_type",sc_method="scRNAseq"):
    all_sc_adatas=[]
    for sc_sample in os.listdir(data_dir):
        sc_adata=sc.read_h5ad(data_dir+"/"+sc_sample+"/"+sc_sample+"_"+sc_method+".h5ad")
        sc_adata.obs["cell_type_combined"]=sc_adata.obs[ct_column].astype(str)
        sc_adata.obs.loc[sc_adata.obs["compartments"]=="Malignant","cell_type_combined"]=sc_adata.obs["cell_type_combined"]+sc_sample.split("-")[1]
        sc_adata.obs["cell_type_combined"]=sc_adata.obs["cell_type_combined"]+sc_adata.obs["condition"].astype(str)
        all_sc_adatas.append((sc_sample,sc_adata))
    all_sc_adatas = pd.DataFrame(all_sc_adatas,columns=('index','data')).set_index('index')['data']
    sc_raw_all=all_sc_adatas.iloc[0].concatenate(all_sc_adatas.iloc[1:],batch_categories=all_sc_adatas.index,batch_key = "sample",join="outer")
    return sc_raw_all


#collect and concatenate plots per sample

import sys
from PIL import Image
def collect_plots(plots_dir,name):
    img=plots_dir+"HTAPP*_"+name+".png"
    source_glob=glob.glob(img, recursive=True)
    source_glob.sort()
    images = [Image.open(x) for x in source_glob]
    widths, heights = zip(*(i.size for i in images))
    
    total_height = sum(heights)
    max_width = max(widths)
    
    new_im = Image.new('RGB', (max_width,total_height), color='#FFF')
    
    y_offset = 0
    for im in images:
        new_im.paste(im, (0,y_offset))
        y_offset += im.size[1]
    
    new_im.save(plots_dir+'00_'+name+'.jpg')
   

#max assign annotations
def typing_to_obs(adata,typings):
    for typing in typings:
        adata.obs[typing]=list(adata.obsm[typing].apply(lambda x: adata.obsm[typing].columns[np.argmax(x)],axis=1))
        
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
    return adata[sample_selection].copy()