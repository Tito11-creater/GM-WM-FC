![GM-WM FC](./GM-WM-FC.png)

# Grey Matter-White Matter Functional Connectome (GM-WM FC)

Welcome to the GM-WM FC package; a companion to our article "Mapping the functional connectome between grey matter and white matter".

Paper link: XXX

Corresponding authors:
Wei LIAO (weiliao.wl@gmail.com) and Jiao LI (jiaoli@uestc.edu.cn), The Clinical Hospital of Chengdu Brain Science Institute, MOE Key Laboratory for Neuroinformation, University of Electronic Science and Technology of China, Chengdu 611731, P.R. China. Fax: +86-28-61831626. Tel: +86-28-61831626

# Data Sources

HCP neuroimaging data (https://www.humanconnectome.org/software/connectomedb)

AD neuroimaging data (https://adni.loni.usc.edu/data-samples/adni-data)

ASD neuroimaging data (https://fcon_1000.projects.nitrc.org/indi/abide/abide_I.html)

MDD neuroimaging data (https://www.synapse.org/Synapse:syn22317081)

Neurotransmitter density PET data (https://github.com/netneurolab/hansen_receptors)

Neurosynth activation maps (https://www.neurosynth.org/)

Schaefer-400 atlas (https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal) for GM parcellation

MWMA-200 atlas (https://github.com/Tito11-creater/MWMA) for WM parcellation

# Software Dependencies
MATLAB (tested on R2020a version - https://www.mathworks.com/products/matlab.html)

R (tested on 2023.06.10-524 version - https://dailies.rstudio.com/version/2023.06.1+524/)

# Toolbox Dependencies

SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) for fMRI data preprocessing

DPABI (https://rfmri.org/DPABI) for fMRI data preprocessing

BrainSpace (https://brainspace.readthedocs.io/en/latest/) for functional gradients analysis

CPM (https://github.com/YaleMRRC/CPM) for prediction model analysis

NeuroSynth decoder (https://github.com/gpreti/GSP_StructuralDecouplingIndex/blob/master/Code_NCOMMS/Python/.ipynb_checkpoints)

relaimpo (https://cran.r-project.org/web/packages/relaimpo/index.html) for assessing relative importance in linear models

ITK-SNAP (https://www.itksnap.org/pmwiki/pmwiki.php) for visualizing 3D brain atlas images

slanCM (https://github.com/slandarer/slanColor) for visualizing FC matrices
