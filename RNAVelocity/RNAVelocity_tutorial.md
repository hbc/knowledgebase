## Running RNA Velocity with scVelo

Heather Wick

There are several ways to approach RNA velocity, once you have [generated your loom files](./velocyto_tutorial.md).

* One way is to analyze the loom files directly, as outlined in [Mary Piper's tutorial](https://github.com/hbc/tutorials/blob/7ff670c5e3b477da09b6c2e832e05bd43e25448f/scRNAseq/scRNAseq_analysis_tutorial/lessons/velocity.md). However, if you already have a Seurat object, you may not want to re-analyze your data.

* Another theoretical way to approach would be to merge your loom files and subset them based on a pre-existing, annotated Seurat file in R, also described by Mary [here](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/seurat_loom_subset_velocity.md). However, several packages, including Seurat, have changed enough that this can be highly challenging, especially if you are working with Seurat objects generated with different versions of Seurat, which makes it not possible to frankenstein the objects together in a way that scVelo will comprehend. I originally started with this method and abandoned it. I do not recommend.

* Instead, if you already have an annotated Seurat file I recommend merging your loom files and combining with your Seurat object in python, as outlined [here](https://uci-genpals.github.io/pseudotime/2021/02/09/scvelo-tutorial.html). For your convenience, I have outlined below the exact steps I took on O2 to do this. There are some small differences between the linked tutorial and what I executed.

**Please note these are real scripts linking to live files and folders; please do not overwrite**

### Step 1: create loom files from CellRanger output

This step is described in the [velocyto tutorial](./velocyto_tutorial.md)

### Step 2: create annotated, clustered Seurat object for your data using the methods of your choice

In my case, this was a Seurat v4 object saved as an RDS.

### Step 3: convert Seurat object into an anndata object

I always loaded my velocyto environment and modules from the velocyto tutorial first:

```
module load gcc/9.2.0 python/3.8.12
source ~/velocyto/bin/activate
```

You can run this interactively which may be useful for deciding what to label your UMAPs with.
**Note that unlike in R, where you can choose a default resolution at which you can plot your clusters, this may not be the case in python, thus, `seurat_clusters` may show unexpect results. Use the clusters at the specific resolution you want, for example, `SCT_snn_res.0.4` below:**

```
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

os.chdir("/n/data1/cores/bcbio/PIs/werner_neuhausser/neuhausser_scRNA-seq_human_embryo_hbc04528/multiome_combined-2023-10-16/RNAvelocity/scVelo/data/redo")

# load sparse matrix:
X = io.mmread("counts.mtx")

# create anndata object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)

# load cell metadata:
cell_meta = pd.read_csv("metadata.csv")

# load gene names:
with open("gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot UMAPs colored by elements of interest
sc.pl.umap(adata, color='Timepoint', frameon=False, legend_loc='on data', title='Timepoint', save='_timepoint.pdf')
sc.pl.umap(adata, color='celltypes', frameon=False, legend_loc='on data', title='Cell type annotation', save='_celltypes.pdf')
sc.pl.umap(adata, color='seurat_clusters', frameon=False, legend_loc='on data', title='Clusters', save='_clusters.pdf')
sc.pl.umap(adata, color='SCT_snn_res.0.4', frameon=False, legend_loc='on data', title='Clusters res 0.4', save='_clusters_res.0.4.pdf')
sc.pl.umap(adata, color='sample_id', frameon=False, legend_loc='on data', title='Samples', save='_sample_id.pdf')

# save dataset as anndata format
adata.write('my_data.h5ad')
```
### Step 4: concatenate loom and merge with seurat to create combined anndata object for RNA velocity analysis

I have the individual important aspects of this outlined below, but you can [click here to jump directly to the full scripts for step 4 and 5](#scripts-for-steps-4-and-5)

Depending on the size of your loom, you may not be able to run this interactively.

#### 4.1: Rename the cell ids in your loom files to match the cell ids in your Seurat object and make unique variable names

```
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad
import os

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2

os.chdir("/n/data1/cores/bcbio/PIs/werner_neuhausser/neuhausser_scRNA-seq_human_embryo_hbc04528/multiome_combined-2023-10-16/RNAvelocity/scVelo/data/redo")

adata = sc.read_h5ad('my_data.h5ad')

# load loom files for spliced/unspliced matrices for each sample:
ldata1 = scv.read('../../../velocyto_out/day7_old/gex_possorted_bam_0964Y.loom', cache=True)
ldata2 = scv.read('../../../velocyto_out/day7_new/gex_possorted_bam_O5WIP.loom', cache=True)
ldata3 = scv.read('../../../velocyto_out/day9/gex_possorted_bam_32XTM.loom', cache=True)
ldata4 = scv.read('../../../velocyto_out/day11/gex_possorted_bam_GUJM4.loom', cache=True)
ldata5 = scv.read('../../../velocyto_out/day11_new/gex_possorted_bam_PMY57.loom', cache=True)
ldata6 = scv.read('../../../velocyto_out/day13/gex_possorted_bam_3L88O.loom', cache=True)
#may get warnings about non-unique variable names; that's ok we'll fix later

# rename barcodes in order to merge (this mapping will be different for your experiment):
# format: gex_possorted_bam_sample:GATTACAx needs to be sample_id:GATTACA-1
# key for my samples:
# "gex_possorted_bam_0964Y:", "D07_"
# "gex_possorted_bam_O5WIP:", "D07new_"
# "gex_possorted_bam_32XTM:", "D09_"
# "gex_possorted_bam_GUJM4:", "D11old_"
# "gex_possorted_bam_PMY57:", "D11new_"
# "gex_possorted_bam_3L88O:", "D13_"
# "x$", "-1" #replace ending with -1

barcodes = [bc.split(':')[1][:-1] for bc in ldata1.obs.index.tolist()] #get just the nucleotide sequence 
barcodes = [bc + '-1' for bc in barcodes] #add -1 to the end
barcodes = ['D07_'  + bc for bc in barcodes] #add proper sample_id to beginning
ldata1.obs.index = barcodes

barcodes = [bc.split(':')[1][:-1] for bc in ldata2.obs.index.tolist()] #get just the nucleotide sequence 
barcodes = [bc + '-1' for bc in barcodes] #add -1 to the end
barcodes = ['D07new_'  + bc for bc in barcodes] #add proper sample_id to beginning
ldata2.obs.index = barcodes

barcodes = [bc.split(':')[1][:-1] for bc in ldata3.obs.index.tolist()] #get just the nucleotide sequence 
barcodes = [bc + '-1' for bc in barcodes] #add -1 to the end
barcodes = ['D09_'  + bc for bc in barcodes] #add proper sample_id to beginning
ldata3.obs.index = barcodes

barcodes = [bc.split(':')[1][:-1] for bc in ldata4.obs.index.tolist()] #get just the nucleotide sequence 
barcodes = [bc + '-1' for bc in barcodes] #add -1 to the end
barcodes = ['D11old_'  + bc for bc in barcodes] #add proper sample_id to beginning
ldata4.obs.index = barcodes

barcodes = [bc.split(':')[1][:-1] for bc in ldata5.obs.index.tolist()] #get just the nucleotide sequence 
barcodes = [bc + '-1' for bc in barcodes] #add -1 to the end
barcodes = ['D11new_'  + bc for bc in barcodes] #add proper sample_id to beginning
ldata5.obs.index = barcodes

barcodes = [bc.split(':')[1][:-1] for bc in ldata6.obs.index.tolist()] #get just the nucleotide sequence 
barcodes = [bc + '-1' for bc in barcodes] #add -1 to the end
barcodes = ['D13_'  + bc for bc in barcodes] #add proper sample_id to beginning
ldata6.obs.index = barcodes

# make variable names unique
ldata1.var_names_make_unique()
ldata2.var_names_make_unique()
ldata3.var_names_make_unique()
ldata4.var_names_make_unique()
ldata5.var_names_make_unique()
ldata6.var_names_make_unique()
```

#### 4.2 Concatenate the loom files, remove any extraneously added sample designations, and merge with the anndata object (saved from step 3 from your Seurat object)

Check the barcodes from your concatenated ldata object before merging with the anndata object. In my case, it added -0, -2, ... -5 to all of my samples to differentiate them, but that would have created a mismatch with the cell ID naming in my anndata object. This was not something mentioned in the linked tutorial.

Plotting UMAPs again after merging and comparing to the UMAPs generated from just the Seurat object proved to be a useful way to see if something is off in the mapping.

```
# concatenate the looms
ldata = ldata1.concatenate([ldata2, ldata3, ldata4, ldata5, ldata6])

# remove individual sample designations added during concatenation (-0, -1, ... -5)
barcodes = [bc[0:len(bc)-2] for bc in ldata.obs.index.tolist()]
ldata.obs.index = barcodes

# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)

# plot UMAPs to check for consistency with previous UMAPs
sc.pl.umap(adata, color='Timepoint', frameon=False, legend_loc='on data', title='Timepoint', save='_timepoint_postmerge.pdf')
sc.pl.umap(adata, color='celltypes', frameon=False, legend_loc='on data', title='Cell type annotation', save='_celltypes_postmerge.pdf')
sc.pl.umap(adata, color='seurat_clusters', frameon=False, legend_loc='on data', title='Clusters', save='_clusters_postmerge.pdf')
sc.pl.umap(adata, color='SCT_snn_res.0.4', frameon=False, legend_loc='on data', title='Clusters res 0.4', save='_clusters_res.0.4_postmerge.pdf')
sc.pl.umap(adata, color='sample_id', frameon=False, legend_loc='on data', title='Samples', save='_sample_id_postmerge.pdf')
```

### Step 5: Run RNA Velocity

This is perhaps the most straightforward part. Note that if you run the dynamical model you need an extra line of code, which I have commented out below. You will also need to request more memory and time

```
# Part II: RNA Velocity
scv.pl.proportions(adata, groupby='celltypes', save='celltypes.pdf')
scv.pl.proportions(adata, groupby='Timepoint', save='Timepoint.pdf')
#scv.pl.proportions(adata, groupby='seurat_clusters', save='seurat_clusters.pdf') #doesn't work; probably need to convert to string or factor first. These are recommended to plot
#scv.pl.proportions(adata, groupby='SCT_snn_res.0.4', save='SCT_snn_res.0.4.pdf') #doesn't work; probably need to convert to string or factor first. These are recommended to plot
scv.pl.proportions(adata, groupby='sample_id', save='sample_id.pdf')

# pre-process
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
# may get warnings about normalization -- that's ok

# compute velocity; stochastic
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

# compute velocity; dynamic
#scv.tl.recover_dynamics(adata) #for dynamic model only
#scv.tl.velocity(adata, mode='dynamical')
#scv.tl.velocity_graph(adata)

#visualize velocity fields
scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='umap_embedding.pdf')

scv.pl.velocity_embedding_grid(adata, basis='umap', color='celltypes', save='celltypes_embedding_grid.pdf', title='Cell type annotation', scale=0.25)
scv.pl.velocity_embedding_grid(adata, basis='umap', color='Timepoint', save='Timepoint_embedding_grid.pdf', title='Timepoint', scale=0.25)
scv.pl.velocity_embedding_grid(adata, basis='umap', color='seurat_clusters', save='Clusters_embedding_grid.pdf', title='Clusters', scale=0.25)
scv.pl.velocity_embedding_grid(adata, basis='umap', color='sample_id', save='sample_id_embedding_grid.pdf', title='Samples', scale=0.25)
scv.pl.velocity_embedding_grid(adata, basis='umap', color='SCT_snn_res.0.4', save='SCT_snn_res.0.4_embedding_grid.pdf', title='Clusters', scale=0.25)

scv.pl.velocity_embedding_stream(adata, basis='umap', color=['celltypes','Timepoint'], save='celltypes_Timepoint_embedding_stream.pdf', title=['Cell type annotation', 'Timepoint'])
#scv.pl.velocity_embedding_stream(adata, basis='umap', color='Timepoint', save='_Timepoint_embedding_grid.pdf', title='Timepoint', scale=0.25)
#scv.pl.velocity_embedding_stream(adata, basis='umap', color='seurat_clusters', save='_Clusters_embedding_grid.pdf', title='Clusters', scale=0.25)
#scv.pl.velocity_embedding_stream(adata, basis='umap', color='sample_id', save='_sample_id_embedding_grid.pdf', title='Samples', scale=0.25)

# plot velocity of a selected gene; you will need to change this to reflect gene(s) of your choice
scv.pl.velocity(adata, var_names=['FABP5'], color='celltypes', save='_celltypes_velocity.pdf')
scv.pl.velocity(adata, var_names=['FABP5'], color='Timepoint', save='_Timepoint_velocity.pdf')

#save object for downstream analysis
adata.write('merged_adata_ldata_6samples.h5ad')
```

### Step 6: Run downstream analysis

This involves looking at top genes, doing additional QC, looking at subsets of the data, etc. I will flesh this out as we do them. Running interactively would be ideal

## Scripts for Steps 4 and 5:

Submission script:
  
```
#!/bin/sh
#SBATCH -p short
#SBATCH -J scVelo
#SBATCH -o %x_%j.o
#SBATCH -e %x_%j.e
#SBATCH -t 0-00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --mail-type=ALL

#load before running
module load gcc/9.2.0 python/3.8.12 samtools
source ~/velocyto1/bin/activate

python mergeLoomAnndataRNAVelocity.py
```

Python script (Note that if you run the dynamical model you need an extra line of code, which I have commented out below. You will also need to request more memory and time):

```
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad
import os

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2

os.chdir("/n/data1/cores/bcbio/PIs/werner_neuhausser/neuhausser_scRNA-seq_human_embryo_hbc04528/multiome_combined-2023-10-16/RNAvelocity/scVelo/data/redo")

adata = sc.read_h5ad('my_data.h5ad')

# load loom files for spliced/unspliced matrices for each sample:
ldata1 = scv.read('../../../velocyto_out/day7_old/gex_possorted_bam_0964Y.loom', cache=True)
ldata2 = scv.read('../../../velocyto_out/day7_new/gex_possorted_bam_O5WIP.loom', cache=True)
ldata3 = scv.read('../../../velocyto_out/day9/gex_possorted_bam_32XTM.loom', cache=True)
ldata4 = scv.read('../../../velocyto_out/day11/gex_possorted_bam_GUJM4.loom', cache=True)
ldata5 = scv.read('../../../velocyto_out/day11_new/gex_possorted_bam_PMY57.loom', cache=True)
ldata6 = scv.read('../../../velocyto_out/day13/gex_possorted_bam_3L88O.loom', cache=True)
#may get warnings about non-unique variable names; that's ok we'll fix later

# rename barcodes in order to merge (this mapping will be different for your experiment):
# format: gex_possorted_bam_sample:GATTACAx needs to be sample_id:GATTACA-1
# key for my samples:
# "gex_possorted_bam_0964Y:", "D07_"
# "gex_possorted_bam_O5WIP:", "D07new_"
# "gex_possorted_bam_32XTM:", "D09_"
# "gex_possorted_bam_GUJM4:", "D11old_"
# "gex_possorted_bam_PMY57:", "D11new_"
# "gex_possorted_bam_3L88O:", "D13_"
# "x$", "-1" #replace ending with -1

barcodes = [bc.split(':')[1][:-1] for bc in ldata1.obs.index.tolist()] #get just the nucleotide sequence 
barcodes = [bc + '-1' for bc in barcodes] #add -1 to the end
barcodes = ['D07_'  + bc for bc in barcodes] #add proper sample_id to beginning
ldata1.obs.index = barcodes

barcodes = [bc.split(':')[1][:-1] for bc in ldata2.obs.index.tolist()] #get just the nucleotide sequence 
barcodes = [bc + '-1' for bc in barcodes] #add -1 to the end
barcodes = ['D07new_'  + bc for bc in barcodes] #add proper sample_id to beginning
ldata2.obs.index = barcodes

barcodes = [bc.split(':')[1][:-1] for bc in ldata3.obs.index.tolist()] #get just the nucleotide sequence 
barcodes = [bc + '-1' for bc in barcodes] #add -1 to the end
barcodes = ['D09_'  + bc for bc in barcodes] #add proper sample_id to beginning
ldata3.obs.index = barcodes

barcodes = [bc.split(':')[1][:-1] for bc in ldata4.obs.index.tolist()] #get just the nucleotide sequence 
barcodes = [bc + '-1' for bc in barcodes] #add -1 to the end
barcodes = ['D11old_'  + bc for bc in barcodes] #add proper sample_id to beginning
ldata4.obs.index = barcodes

barcodes = [bc.split(':')[1][:-1] for bc in ldata5.obs.index.tolist()] #get just the nucleotide sequence 
barcodes = [bc + '-1' for bc in barcodes] #add -1 to the end
barcodes = ['D11new_'  + bc for bc in barcodes] #add proper sample_id to beginning
ldata5.obs.index = barcodes

barcodes = [bc.split(':')[1][:-1] for bc in ldata6.obs.index.tolist()] #get just the nucleotide sequence 
barcodes = [bc + '-1' for bc in barcodes] #add -1 to the end
barcodes = ['D13_'  + bc for bc in barcodes] #add proper sample_id to beginning
ldata6.obs.index = barcodes

# make variable names unique
ldata1.var_names_make_unique()
ldata2.var_names_make_unique()
ldata3.var_names_make_unique()
ldata4.var_names_make_unique()
ldata5.var_names_make_unique()
ldata6.var_names_make_unique()

# concatenate the looms
ldata = ldata1.concatenate([ldata2, ldata3, ldata4, ldata5, ldata6])

# remove individual sample designations added during concatenation (-0, -1, ... -5)
barcodes = [bc[0:len(bc)-2] for bc in ldata.obs.index.tolist()]
ldata.obs.index = barcodes

# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)

# plot UMAPs to check for consistency with previous UMAPs
sc.pl.umap(adata, color='Timepoint', frameon=False, legend_loc='on data', title='Timepoint', save='_timepoint_postmerge.pdf')
sc.pl.umap(adata, color='celltypes', frameon=False, legend_loc='on data', title='Cell type annotation', save='_celltypes_postmerge.pdf')
sc.pl.umap(adata, color='seurat_clusters', frameon=False, legend_loc='on data', title='Clusters', save='_clusters_postmerge.pdf')
sc.pl.umap(adata, color='SCT_snn_res.0.4', frameon=False, legend_loc='on data', title='Clusters res 0.4', save='_clusters_res.0.4_postmerge.pdf')
sc.pl.umap(adata, color='sample_id', frameon=False, legend_loc='on data', title='Samples', save='_sample_id_postmerge.pdf')

# RNA Velocity
scv.pl.proportions(adata, groupby='celltypes', save='celltypes.pdf')
scv.pl.proportions(adata, groupby='Timepoint', save='Timepoint.pdf')
#scv.pl.proportions(adata, groupby='seurat_clusters', save='seurat_clusters.pdf') #doesn't work; probably need to convert to string or factor first. These are recommended to plot
#scv.pl.proportions(adata, groupby='SCT_snn_res.0.4', save='SCT_snn_res.0.4.pdf') #doesn't work; probably need to convert to string or factor first. These are recommended to plot
scv.pl.proportions(adata, groupby='sample_id', save='sample_id.pdf')

# pre-process
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
# may get warnings about normalization -- that's ok

# compute velocity; stochastic
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

# compute velocity; dynamic
#scv.tl.recover_dynamics(adata) #for dynamic model only
#scv.tl.velocity(adata, mode='dynamical')
#scv.tl.velocity_graph(adata)

#visualize velocity fields
scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='umap_embedding.pdf')

scv.pl.velocity_embedding_grid(adata, basis='umap', color='celltypes', save='celltypes_embedding_grid.pdf', title='Cell type annotation', scale=0.25)
scv.pl.velocity_embedding_grid(adata, basis='umap', color='Timepoint', save='Timepoint_embedding_grid.pdf', title='Timepoint', scale=0.25)
scv.pl.velocity_embedding_grid(adata, basis='umap', color='seurat_clusters', save='Clusters_embedding_grid.pdf', title='Clusters', scale=0.25)
scv.pl.velocity_embedding_grid(adata, basis='umap', color='sample_id', save='sample_id_embedding_grid.pdf', title='Samples', scale=0.25)
scv.pl.velocity_embedding_grid(adata, basis='umap', color='SCT_snn_res.0.4', save='SCT_snn_res.0.4_embedding_grid.pdf', title='Clusters', scale=0.25)

scv.pl.velocity_embedding_stream(adata, basis='umap', color=['celltypes','Timepoint'], save='celltypes_Timepoint_embedding_stream.pdf', title=['Cell type annotation', 'Timepoint'])
#scv.pl.velocity_embedding_stream(adata, basis='umap', color='Timepoint', save='_Timepoint_embedding_grid.pdf', title='Timepoint', scale=0.25)
#scv.pl.velocity_embedding_stream(adata, basis='umap', color='seurat_clusters', save='_Clusters_embedding_grid.pdf', title='Clusters', scale=0.25)
#scv.pl.velocity_embedding_stream(adata, basis='umap', color='sample_id', save='_sample_id_embedding_grid.pdf', title='Samples', scale=0.25)

# plot velocity of a selected gene; you will need to change this to reflect gene(s) of your choice
scv.pl.velocity(adata, var_names=['FABP5'], color='celltypes', save='_celltypes_velocity.pdf')
scv.pl.velocity(adata, var_names=['FABP5'], color='Timepoint', save='_Timepoint_velocity.pdf')

#save object for downstream analysis
adata.write('merged_adata_ldata_6samples.h5ad')
```








