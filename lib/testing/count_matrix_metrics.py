import numpy as np
import scipy.sparse as sp_sparse
import h5py

cr_basedir = '/mnt/home/daniel.liu/cellranger/pbmc1k_star/outs/'
thermite_basedir = '/mnt/home/daniel.liu/cellranger/pbmc1k_thermite/outs/'

def load_h5_file(basedir):
    h5_file = basedir + "/raw_feature_bc_matrix.h5"
    h5 = h5py.File(h5_file, 'r')
    matrix= h5['matrix']
    barcodes_np = matrix.get('barcodes')[()]
    data_np = matrix.get('data')[()]
    indices_np = matrix.get('indices')[()]
    indptr_np = matrix.get('indptr')[()]
    shape = matrix.get('shape')[()]
    count_matrix = sp_sparse.csc_matrix((data_np, indices_np, indptr_np), shape=shape)
    features_np = matrix.get('features').get('id')[()]
    return count_matrix, features_np

thermite_count_matrix, thermite_features = load_h5_file(thermite_basedir)
cellranger_count_matrix, cr_features = load_h5_file(cr_basedir)
nonzero = cellranger_count_matrix.nonzero()
thermite = thermite_count_matrix[nonzero]
cellranger = cellranger_count_matrix[nonzero]

print("Thermite shape", thermite_count_matrix.shape)
print("Cellranger shape", cellranger_count_matrix.shape)

thermite_nz = thermite_count_matrix[thermite_count_matrix.nonzero()]
print("Thermite mean of elements", thermite_nz.mean())
print("Thermite min of elements", thermite_nz.min())
print("Thermite max of elements", thermite_nz.max())

cellranger_nz = cellranger_count_matrix[cellranger_count_matrix.nonzero()]
print("Cellranger mean of elements", cellranger_nz.mean())
print("Cellranger min of elements", cellranger_nz.min())
print("Cellranger max of elements", cellranger_nz.max())

abs_diff = np.abs(thermite - cellranger)
print("Mean abs diff of elements", abs_diff.mean())
print("Min abs diff of elements", abs_diff.min())
print("Max abs diff of elements", abs_diff.max())

rel_diff = abs_diff / np.maximum(thermite, cellranger)
print("Mean rel diff of elements", rel_diff.mean())
print("Min rel diff of elements", rel_diff.min())
print("Max rel diff of elements", rel_diff.max())

print("Thermite num nonzero", thermite_count_matrix.nnz)
print("Cellranger num nonzero", cellranger_count_matrix.nnz)
nne = (thermite_count_matrix != cellranger_count_matrix).nnz
print("Num not equal", nne)
print("Fraction not equal", nne / cellranger_count_matrix.nnz)

abs_diff_nz = abs_diff[abs_diff.nonzero()]
print("Mean abs diff of nonzero elements", abs_diff_nz.mean())
print("Min abs diff of nonzero elements", abs_diff_nz.min())
print("Max abs diff of nonzero elements", abs_diff_nz.max())

rel_diff_nz = abs_diff_nz / np.maximum(thermite, cellranger)[abs_diff.nonzero()]
print("Mean rel diff of nonzero elements", rel_diff_nz.mean())
print("Min rel diff of nonzero elements", rel_diff_nz.min())
print("Max rel diff of nonzero elements", rel_diff_nz.max())
