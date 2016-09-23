import dipy
import numpy as np
import nibabel as nib
from nibabel import trackvis as tv
from dipy.tracking.streamline import set_number_of_points
from dipy.segment.mask import median_otsu
from dipy.io import read_bvals_bvecs
from dipy.core.gradients import gradient_table
from dipy.reconst.dti import TensorModel
from dipy.reconst.dti import fractional_anisotropy
from dipy.reconst.dti import color_fa
from dipy.reconst.shm import CsaOdfModel
from dipy.data import get_sphere
from dipy.reconst.peaks import peaks_from_model
from dipy.tracking.eudx import EuDX
from dipy.tracking.utils import density_map
#from dipy.viz import fvtk
from dipy.tracking import utils
import matplotlib.pyplot as plt

fimg = "DTICAP_bet.nii.gz"
img = nib.load(fimg)
data = img.get_data()
affine = img.get_affine()
header = img.get_header() 
voxel_size = header.get_zooms()[:3]
mask, S0_mask = median_otsu(data[:, :, :, 0])
fbval = "../TXT/DTICAP.bval"
fbvec = "../TXT/DTICAP.bvec"

bvals, bvecs = read_bvals_bvecs(fbval, fbvec)
gtab = gradient_table(bvals, bvecs)
ten_model = TensorModel(gtab)
ten_fit = ten_model.fit(data, mask)

fa = fractional_anisotropy(ten_fit.evals)
cfa = color_fa(fa, ten_fit.evecs)
csamodel = CsaOdfModel(gtab, 6)
sphere = get_sphere('symmetric724')

pmd = peaks_from_model(model=csamodel,
                       data=data,
                       sphere=sphere,
                       relative_peak_threshold=.5,
                       min_separation_angle=25,
                       mask=mask,
                       return_odf=False)
#Deterministic tractography 
eu = EuDX(a=fa, ind=pmd.peak_indices[..., 0], seeds=2000000, odf_vertices=sphere.vertices, a_low=0.01)
affine = eu.affine
csd_streamlines= list(eu)

#Remove tracts shorter than 30mm
#print np.shape(csd_streamlines)
from dipy.tracking.utils import length 
csd_streamlines=[t for t in csd_streamlines if length(t)>30]
 

#Trackvis
hdr = nib.trackvis.empty_header()
hdr['voxel_size'] = img.get_header().get_zooms()[:3]
hdr['voxel_order'] = 'LAS'
hdr['dim'] = fa.shape
tensor_streamlines_trk = ((sl, None, None) for sl in csd_streamlines)
ten_sl_fname = 'tensor_streamlines.trk'
nib.trackvis.write(ten_sl_fname, tensor_streamlines_trk, hdr, points_space='voxel')

print np.shape(csd_streamlines)
atlas = nib.load('atlas_reg.nii.gz')
labels = atlas.get_data()
labelsint = labels.astype(int)
 
#M, grouping = utils.connectivity_matrix(csd_streamlines, labelsint, affine=affine,    return_mapping=True,  mapping_as_streamlines=True)
M = utils.connectivity_matrix(csd_streamlines, labelsint, affine=affine  )

#Remove background
M = M[1:,1:]
#Remove the last rows and columns since they are cerebellum and brainstem
M = M[:90,:90]
'''
#Reshuffle making all left areas first right areas
odd_odd = M[::2, ::2]
odd_even = M[::2, 1::2]
first = np.vstack((odd_odd,odd_even))
even_odd = M[1::2, ::2]
even_even= M[1::2, 1::2]
second = np.vstack((even_odd,even_even))
M = np.hstack((first,second))
'''
np.fill_diagonal(M,0)

print np.shape(M)
np.savetxt("foo.csv", M, delimiter=",")
#np.savetxt('connectome.txt', M) 
 
tem = M.sum(axis=1)
tem2 = tem.sum(axis=0)
print tem2

#plt.imshow(np.log1p(M), interpolation='nearest')
#plt.savefig("connectivity.png")
