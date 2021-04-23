These scripts aim at computing the effective connectivity by using functional and structural MRI data. 

The main function to compute the autoregressive model as in <a href="https://link.springer.com/chapter/10.1007/978-3-319-46720-7_17" target="_blank"> A.Crimi et al. "Effective Brain Connectivity Through a Constrained Autoregressive Model" MICCAI 2016
</a> 
The path is assumed to be of the data from the NKI 1000Connectome project (in the original path definition, to be used in BIDS structure, the path has to be changed)
Some interactive visualization is acceesible on the website
<a href="https://www.effectiveconnectivity.net" target="_blank">www.effectiveconnectivity.net
</a> 

Main function (structured_G_causality.m)
Input: 

name_session: the folder with the subjects for the same session (or folder with all subjects)

name_fold: main folder for single subject e.g.'MR_9421819_1328';

norm_opt: = 0 means no normalization, = 1 means binarize, = 2 means divide by the maximum
(Although convergence should be reached anyhow)

use_multistsp: = 0 only direct connections, 1 use direct and indirect connections

use_deconv: = 0 use normal BOLD signal, 1 perform the deconvolution of the HRF (this will require SPM12 which is not included in this repo)


Output:

M: the resulting effective connectivity matrix

etem: the reconstruction error over the iterations

![alt text](https://github.com/alecrimi/effective_connectivity_toolbox/blob/master/nft.jpg)


Data are loaded using the script of Jimmy Shen to load NIFTI and Analyze files, those are in the folder libs.

The code can be run on "pure" pre-processed BOLD series, or on the deconvoluted version which removes the hemodynamics response. In this latter case, the code related to  Wu et al, "A blind deconvolution approach to recover effective connectivity brain networks from resting state fMRI data, (doi: 10.1016/j.media.2013.01.003) is used, and for this installation  SPM12 is necessary.


**SIMULATIONS**

You can also test well known ground-truth BOLD simulations from <a href="https://www.fmrib.ox.ac.uk/datasets/netsim/" target="_blank"> Netsim (Smith et al. Neuroimage 2011) </a> 
In the folder *simulations* there are the basic scripts to load simulations and run the CMAR model.


**Data preprocessed and final matrices**
Data preprocessed and final matrices for part of the used database are available <a href="https://doi.org/10.5281/zenodo.4711994" target="_blank">here.</a>
