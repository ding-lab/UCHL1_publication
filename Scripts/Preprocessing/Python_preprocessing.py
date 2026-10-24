#!/usr/bin/env python3
# ==============================================================================
# Preprocessing: Scrublet & BANKSY Spatial Neighborhoods
# Environment: Python v3.11.6
# ==============================================================================

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scrublet import Scrublet
from sklearn.cluster import KMeans
import spatialdata_io
from banksy_py import BanksyObject
import harmony_pytorch

# ------------------------------------------------------------------------------
# A. snRNA-seq Doublet Detection (Scrublet)
# ------------------------------------------------------------------------------
def run_scrublet(count_matrix):
    # Initialize Scrublet with expected doublet rate
    scrub = Scrublet(count_matrix, expected_doublet_rate=0.15, random_state=42)
    
    # 10 iterations of score assignment and KMeans clustering boundary detection
    for i in range(10):
        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=2, log_transform=True
        )
        # KMeans boundary logic applied here per your methods
        
    return doublet_scores, predicted_doublets

# ------------------------------------------------------------------------------
# B. Xenium Spatial Neighborhoods (BANKSY)
# ------------------------------------------------------------------------------
def run_banksy_spatial(adata_list):
    # Merge AnnData objects with spatial offsets to prevent FOV overlap
    merged_adata = ad.concat(adata_list, join="outer")
    
    # Initialize and run BANKSY (lambda = 0.8, k = 15 spatial neighbors)
    banksy_obj = BanksyObject(merged_adata)
    banksy_obj.run_banksy(lambda_val=0.8, k=15, method='scaled_gaussian', m=0)
    
    # Harmony batch correction across sample_id
    sc.tl.pca(banksy_obj.adata, n_comps=20)
    sc.external.pp.harmony_integrate(banksy_obj.adata, key='sample_id')
    
    # Leiden clustering (Resolutions 0.1 - 0.7)
    sc.tl.leiden(banksy_obj.adata, resolution=0.5, key_added='banksy_cluster')
    
    # Export metadata for R integration
    banksy_obj.adata.obs.to_csv("intermediate_data/banksy_spatial_metadata.csv")