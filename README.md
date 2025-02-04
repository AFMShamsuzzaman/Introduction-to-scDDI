# Introduction to scDDI (Darmanis Data)
A novel dropout detection and imputation method

# Summary

A novel approach leveraging Poisson-Negative Binomial mixture models, cell-to-cell similarity calculations, and decision tree regression for accurate detection and imputation of dropout events in scRNA-seq data.



# How to use scDDI

## Data Loading and Preprocessing


Load the Libraries

```
library(SingleCellExperiment)
library('Linnorm')
```

We will give the method demonstration on Darmanis Dataset. For more details about the data, see [A survey of human brain transcriptome diversity at the single cell level](https://www.pnas.org/content/112/23/7285#:~:text=Our%20results%20show%20that%20MHCI,as%20endothelial%20cells%20and%20microglia.)

Read the gene expression data using (SingleCellExperiment object), calculate CPM values and extract metadata.


```
rawdata <- readRDS("darmanis.rds")
data <- assay(rawdata)
```
For demonstration purposes, we apply a standard *Linnorm* normalization with minimum read count =5 in 10 percent cell. However any other normalization approach may be used.
Gene should be in row, Cells should be in coloumn


```
darmanis_process= normalized_data(data)

```

```
dim(darmanis_process) 
[1] 8994    466

preprocessedata[1:2,1:3]
      Brain    Brain    Brain
A2M  0.000000 4.953487 4.908761
AAAS 1.526881 0.000000 0.000000
```

A total of 466 cells and 8994 genes are remaining in the dataset after cell, gene filtering, and Normalization.

## Calculation of dropout probability matrix and cell-to-cell similarity matrix

Load the libraries

```
library(foreach)
library(doParallel)
library(scDoc)

```

Now, calculate dropout probabilty matrix and cell-to-cell similarity matrix as follows :

```
offsets_darmanis <- as.numeric(log(colSums(darmanis_process)))
dp_darmanis <- prob.dropout(input = darmanis_process, offsets = offsets_darmanis, mcore = 6)  ## dp_darmanis is the dropout probability matrix
sim_darmanis <- sim.calc(log2(count_darmanis+1), dp_darmanis)   ## sim_darmanis is the cell-to-cell similarity matrix

```


## Saving the results

Then, save the processed dataset, dropout probability matrix and cell-to-cell similarity matrix into csv file fromat. 

      write.csv(darmanis_process,"/home/zaman/New2/darmanis_process.csv",row.names = FALSE)
      write.csv(dp_darmanis,"/home/zaman/New2/dp_darmanis.csv",row.names = FALSE)
      write.csv(sim_darmanis,"/home/zaman/New2/sim_darmanis.csv",row.names = FALSE)

Or you can also simply run the Rscript file as follows:

      Rscript Imputation.R


## Usage of the Python functions 

Now, Run the following python code to impute the given dataset as follows:

      python3 imputation_scDDI.py


## Clustering of the imputed dataset and calculation of Adjusted Rand Index(ARI)
    
Then, to validate the results, we utilized clustering performance metrics, specifically the Adjusted Rand Index (ARI). We compared the ARI value for both unimputed dataset and as well as imputed dataset using scDDI.

import libraries in python and importing the data

```
import numpy as np
import pandas as pd
import scanpy as sc
adata=sc.read_csv('darmanis_imputed.csv',delimiter=',', first_column_names=None, dtype='float32')
      ```












      
Now,we visualized the clustering outcomes and analyzed the clusters to identify marker genes.

    jupyter notebook marker_gene.ipynb
```

Using PCA dimensionality reduction and Leiden clustering

```
sc.tl.pca(adata1, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata1,n_pcs=20,log=True)
```
<img src="./pca_darmanis.png">


```
#create neighborhood graph using 20 pcs 
sc.pp.neighbors(adata1, n_neighbors=15, n_pcs=20)
##dim reduction using umap
sc.tl.umap(adata1)
#Leiden clustering
import leidenalg
sc.tl.leiden(adata1)
##visualizing clusters
sc.pl.umap(adata1, color=['leiden'])
```
<img src="./cluster_darmanis.png">

save the clustering results

```
pd.DataFrame(adata1.obs).to_csv("darmanis_leiden.csv")
```


## Marker selection form identified clusters

Compute a ranking for the highly differential genes in each cluster using Wilcoxon-RankSum test

```
sc.tl.rank_genes_groups(adata1, 'leiden', method='wilcoxon',key_added = "wilcoxon")
sc.pl.rank_genes_groups(adata1, n_genes=30, sharey=False,key="wilcoxon")
```

<img src="./wilcox_darmanis.png">

Top 10 DE genes for each cluster using Wilcox-Ranksum Test

```
#save top 10 DE genes for each cluster using wilcox-ranksum test
result = adata1.uns['wilcoxon']
groups = result['names'].dtype.names
p=pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(10)
pd.DataFrame(p).to_csv("darmanis_marker.csv")
```
Visualizing top 5 DE genes for each cluster in a heatmap using wilcox results

```
sc.pl.rank_genes_groups_heatmap(adata1, n_genes=5, key="wilcoxon", groupby="leiden", show_gene_labels=True)
```
<img src="./heat_darmanis.png">

Visualizing top 5 DE genes for each cluster in a dotplot using t-test results. Here color of dot represents mean expression of the gene in those cell, dot size represents fraction of cells expressing a gene  

```
sc.pl.rank_genes_groups_dotplot(adata1, n_genes=5, key="wilcoxon", groupby="leiden")
```
<img src="./dotplot_darmanis.png">

Visualizing top 5 DE genes for each cluster in a stacked violin plot using t-test results 

```
sc.pl.rank_genes_groups_stacked_violin(adata1, n_genes=5, key="wilcoxon", groupby="leiden")
```
<img src="./violin_darmanis.png">

Visualizing top 5 DE genes for each cluster in a matrixplot using wilcox results. matrixplot represents mean expression of a gene in a cluster as a heatmap.

```
sc.pl.rank_genes_groups_matrixplot(adata1, n_genes=5, key="wilcoxon", groupby="leiden")
```
<img src="./heat2_darmanis.png">

Showing expression of some marker genes (e.g VIP,DCX) across Leiden groups

```
sc.pl.violin(adata1, ['VIP'], groupby='leiden')
```
<img src="./VIP_darmanis.png">

```
sc.pl.violin(adata1, ['DCX'], groupby='leiden')
```
<img src="./dcx_darmanis.png">





