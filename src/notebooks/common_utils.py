# Shared utility functions for the notebooks

import scanpy as sc
import anndata
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import harmonypy
import pickle
import numpy as np
import matplotlib as mpl
import matplotlib.font_manager
from matplotlib import font_manager
from matplotlib.font_manager import fontManager, FontProperties
import infercnvpy as cnv


def setup_dirs(outDir):
    """
    Setup directories for figures, data, and tables.

    Parameters
    ----------
    outDir : str
        Output directory
    
    Returns
    -------
    str, str, str
        Figures directory, data directory, tables directory
    
    Notes
    -----
    
    Creates the following directories:
    1. Figures directory: outDir/figures
    2. Data directory: outDir/data
    3. Tables directory: outDir/tables
    """
    figuresDir = os.path.join(outDir, 'figures')
    dataDir = os.path.join(outDir, 'data')
    tablesDir = os.path.join(outDir, 'tables')
    os.makedirs(figuresDir, exist_ok=True)
    os.makedirs(dataDir, exist_ok=True)
    os.makedirs(tablesDir, exist_ok=True)
    return figuresDir, dataDir, tablesDir

def local_arial_font(arial_font_path):
    """
    Add Arial font given its path.
    """
    font_manager.fontManager.addfont(arial_font_path)
    prop = font_manager.FontProperties(fname=arial_font_path)

# set the font
def find_arial_font():
    """
    Find Arial font and set it as the default font.
    """
    arial_font = None
    for font in font_manager.findSystemFonts():
        if "arial" in font.lower():
            arial_font = font
            break
        if arial_font:
            print("Found Arial font at: ", arial_font)
            prop = font_manager.FontProperties(fname=arial_font)
            sns.set(font=prop.get_name())
    if arial_font is None:
        print("Arial font not found")
        arial_font_path = PATH_TO_ARIAL_FONT
        local_arial_font(arial_font_path)


def filter_genes(adata):
    """
    Filtering the following genes to avoid the dominant effect of
    IG{}V, (Immunoglobulin variable)
    TR{}, (T cell receptor variable genes)
    linc, (Long intergenic non-coding),
    genes starting with RP (ribosomal protein),
    genes starting with MT- (mitochondrial genes)
    HLA genes
    """
    genes = [x for x in adata.var.index.tolist() if "MT-" not in x]
    genes = [x for x in genes if "." not in x]
    genes = [x for x in genes if not x.startswith("RP")]
    genes = [x for x in genes if "linc" not in x.lower()]
    genes = [x for x in genes if "TRA" not in x.upper()]
    genes = [x for x in genes if "TRB" not in x.upper()]
    genes = [x for x in genes if "TRG" not in x.upper()]
    genes = [x for x in genes if "TRD" not in x.upper()]
    genes = [x for x in genes if "IGKV" not in x.upper()]
    genes = [x for x in genes if "IGHV" not in x.upper()]
    genes = [x for x in genes if "IGLV" not in x.upper()]
    genes = [x for x in genes if "-" not in x.upper() and "HLA" not in x.upper()]
    adata = adata[:, genes].copy()
    return adata


def add_gene_binary_status(adata, gene, threshold=np.log1p(1), use_counts=False):
    """
    Find a cut-off for expressed vs not expressed

    Parameters
    ----------
    adata : AnnData
        AnnData object
    gene : str
        Gene to check
    threshold : float
        Threshold for expressed vs not expressed (default: np.log1p(1))
    use_counts : bool
        Use counts instead of normalized data (default: False)
    
    Returns
    -------
    AnnData
        AnnData with a new column {gene}_is_expressed

    Notes
    ----
    Very simple, expressed, if there is more than 1 count (i.e., log1p > log(1+1) = np.log1p(1))
    Add a column to adata.obs, {gene}_is_expressed. 
    Will overwrite the column if it exists.
    """    
    assert gene in adata.var_names, f"Gene {gene} not in adata.var_names..."
    if f'{gene}_is_expressed' in adata.obs.columns:
        adata.obs.drop(f'{gene}_is_expressed', axis=1, inplace=True)
    if f'{gene}_EXPR' in adata.obs.columns:
        adata.obs.drop(f'{gene}_EXPR', axis=1, inplace=True)
    if use_counts:
        adata.obs[f'{gene}_EXPR'] = adata.layers['counts'][:, (adata.var_names == gene)].toarray()
    else:
        adata.obs[f'{gene}_EXPR'] = adata[:, gene].X.A
    adata.obs[f'{gene}_is_expressed'] = adata.obs[f'{gene}_EXPR'] > threshold
    adata.obs[f'{gene}_is_expressed'] = adata.obs[f'{gene}_is_expressed'].astype('category')    
    return adata 


def mini_process(adata, use_harmony=False, do_scale=False, max_iter_harmony=100, harmony_column='sample'):
    """
    Minimal Processing on the top 2000 genes. 

    Parameters
    ----------
    adata : AnnData
        AnnData object
    use_harmony : bool
        Use Harmony for batch correction (default: False)
    do_scale : bool
        Scale the data (default: False)
    max_iter_harmony : int
        Maximum number of iterations for Harmony (default: 100)
    harmony_column : str
        Column to use for Harmony (default: 'sample')
    
    Returns
    -------
    AnnData
        Processed AnnData. See Notes for details.
     
    Notes
    -----
    1. Normalize the data
    2. Log1p the data
    3. Scale the data
    4. Run PCA
    5. Run Harmony
    6. Run UMAP
    """
    adata.X = adata.layers['counts'].copy()
    assert adata.X.max() > 100, "The data is not raw counts"
    assert adata.X.min() >= 0, "The data is not raw counts"
    # Normalize the data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    if do_scale:
        sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    use_rep = 'X_pca'
    if use_harmony:
        sc.external.pp.harmony_integrate(adata, key=harmony_column, max_iter_harmony=max_iter_harmony)
        use_rep = 'X_pca_harmony'
    sc.pp.neighbors(adata, use_rep=use_rep)
    sc.tl.umap(adata)
    return adata
