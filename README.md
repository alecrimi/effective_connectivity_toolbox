These scripts aim at computing the effective connectivity by using functional and structural MRI data. 

The main function to compute the autoregressive model as in <a href="https://link.springer.com/chapter/10.1007/978-3-319-46720-7_17" target="_blank"> A.Crimi et al. "Effective Brain Connectivity Through a Constrained Autoregressive Model" MICCAI 2016
</a> 
The path is assumed to be of the data from the NKI 1000Connectome project (in the original path definition, to be used in BIDS structure, the path has to be changed)

Input: 

name_session: the folder with the subjects for the same session (or folder with all subjects)

name_fold: main folder for single subject e.g.'MR_9421819_1328';

norm_opt: = 0 means no normalization, = 1 means binarize, = 2 means divide by the maximum
(Although convergence should be reached anyhow)

Output:

M: the resulting effective connectivity matrix

etem: the reconstruction error over the iterations

![alt text](https://github.com/alecrimi/effective_connectivity_toolbox/blob/master/nft.jpg)


Data are loaded using the script of Jimmy Shen to load NIFTI and Analyze files, those are in the folder libs.

The code can be run on "pure" pre-processed BOLD series, or on the deconvoluted version which removes the hemodynamics response. In this latter case, the code related to  Wu et al, "A blind deconvolution approach to recover effective connectivity brain networks from resting state fMRI data, (doi: 10.1016/j.media.2013.01.003) is used, and for this installation of SPM12 is necessary.


**SIMULATIONS**

You can also test well known ground-truth BOLD simulations from <a href="https://www.fmrib.ox.ac.uk/datasets/netsim/" target="_blank"> Netsim (Smith et al. Neuroimage 2011) </a> 
In the folder *simulations* there are the basic scripts to load simulations and run the CMAR model.
