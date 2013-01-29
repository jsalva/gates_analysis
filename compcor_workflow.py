
import os
import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio          # input/output
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as util
from nipype.algorithms.misc import TSNR
import nipype.interfaces.fsl as fsl
import nipype.algorithms.rapidart as ra     # rapid artifact detection


def pickfirst(files):
    """Return first file from a list of files

    Parameters
    ----------
    files : list of filenames

    Returns
    -------
    file : returns the filename corresponding to the middle run
    """
    if isinstance(files, list):
        return files[0]
    else:
        return files

def extract_csf_mask():
    """Create a workflow to extract a mask of csf voxels
    
    Inputs
    ------
    inputspec.mean_file :
    inputspec.reg_file :
    inputspec.fsaseg_file :
    
    Outputs
    -------
    outputspec.csf_mask :
    
    Returns
    -------
    workflow : workflow that extracts mask of csf voxels
    """
    extract_csf = pe.Workflow(name='extract_csf_mask')
    inputspec = pe.Node(util.IdentityInterface(fields=['mean_file',
                                                       'reg_file',
                                                       'fsaseg_file']),
                        name='inputspec')

    bin = pe.Node(fs.Binarize(), name='binarize')
    bin.inputs.wm_ven_csf = True
    bin.inputs.match = [4, 5, 14, 15, 24, 31, 43, 44, 63]
    #erode 2 voxels into ventricles to make SURE no BG is being removed by compcorr
    bin.inputs.erode = 2
    
    extract_csf.connect(inputspec, 'fsaseg_file',
                        bin, "in_file")
    voltransform = pe.Node(fs.ApplyVolTransform(inverse=True),
                           name='inverse_transform')
    extract_csf.connect(bin, 'binary_file',
                        voltransform, 'target_file')
    extract_csf.connect(inputspec, 'reg_file',
                        voltransform, 'reg_file')
    extract_csf.connect(inputspec, 'mean_file',
                        voltransform, 'source_file')
    outputspec = pe.Node(util.IdentityInterface(fields=['csf_mask']),
                         name='outputspec')
    extract_csf.connect(voltransform, 'transformed_file',
                        outputspec, 'csf_mask')
    return extract_csf


def extract_noise_components(realigned_file, noise_mask_file, num_components,
                             csf_mask_file, selector):
    """Derive components most reflective of physiological noise
    
    Parameters
    ----------
    realigned_file :
    noise_mask_file :
    num_components :
    csf_mask_file :
    selector :
    
    Returns
    -------
    components_file :
    """

    import os
    from nibabel import load
    import numpy as np
    import scipy as sp
    from scipy.signal import detrend

    options = np.array([noise_mask_file, csf_mask_file])
    selector = np.array(selector)
    imgseries = load(realigned_file)
    if selector.all():  # both values of selector are true, need to concatenate
        tcomp = load(noise_mask_file)
        acomp = load(csf_mask_file)
        voxel_timecourses = imgseries.get_data()[np.nonzero(tcomp.get_data() +
                                                            acomp.get_data())]
    else:
        noise_mask_file = options[selector][0]
        noise_mask = load(noise_mask_file)
        voxel_timecourses = imgseries.get_data()[np.nonzero(noise_mask.get_data())]
    for timecourse in voxel_timecourses:
        timecourse[:] = detrend(timecourse, type='constant')
    voxel_timecourses = voxel_timecourses.byteswap().newbyteorder()
    voxel_timecourses[np.isnan(np.sum(voxel_timecourses,axis=1)),:] = 0
    _, _, v = sp.linalg.svd(voxel_timecourses, full_matrices=False)
    components_file = os.path.join(os.getcwd(), 'noise_components.txt')
    np.savetxt(components_file, v[:num_components, :].T)
    return components_file



def create_compcorr(name='CompCor'):
    """Workflow that implements (t and/or a) compcor method from 
    
    Behzadi et al[1]_.
    
    Parameters
    ----------
    name : name of workflow. Default = 'CompCor'
    
    Inputs
    ------
    inputspec.num_components :
    inputspec.realigned_file :
    inputspec.in_file :
    inputspec.reg_file :
    inputspec.fsaseg_file :
    inputspec.selector :
    
    Outputs
    -------
    outputspec.noise_components :
    outputspec.stddev_file :
    outputspec.tsnr_file :
    outputspec.csf_mask :
    
    References
    ----------
    .. [1] Behzadi Y, Restom K, Liau J, Liu TT. A component based\
           noise correction method (CompCor) for BOLD and perfusion\
           based fMRI. Neuroimage. 2007 Aug 1;37(1):90-101. DOI_.

    .. _DOI: http://dx.doi.org/10.1016/j.neuroimage.2007.04.042
    """
    compproc = pe.Workflow(name=name)
    inputspec = pe.Node(util.IdentityInterface(fields=['num_components',
                                                       'realigned_file',
                                                       'mean_file',
                                                       'reg_file',
                                                       'fsaseg_file',
                                                       'selector']),
                        name='inputspec')
    # selector input is bool list [True,True] where first is referring to
    # tCompcorr and second refers to aCompcorr
    outputspec = pe.Node(util.IdentityInterface(fields=['noise_components',
                                                        'stddev_file',
                                                        'tsnr_file',
                                                        'csf_mask',
                                                        'tsnr_detrended']),
                         name='outputspec')
    # extract the principal components of the noise
    tsnr = pe.MapNode(TSNR(regress_poly=2),  #SG: advanced parameter
                      name='tsnr',
                      iterfield=['in_file'])

    # additional information for the noise prin comps
    getthresh = pe.MapNode(interface=fsl.ImageStats(op_string='-p 98'),
                            name='getthreshold',
                            iterfield=['in_file'])

    # and a bit more...
    threshold_stddev = pe.MapNode(fsl.Threshold(),
                                  name='threshold',
                                  iterfield=['in_file', 'thresh'])

    acomp = extract_csf_mask()

    # compcor actually extracts the components
    compcor = pe.MapNode(util.Function(input_names=['realigned_file',
                                                    'noise_mask_file',
                                                    'num_components',
                                                    'csf_mask_file',
                                                    'selector'],
                                       output_names=['noise_components'],
                                       function=extract_noise_components),
                                       name='compcor_components',
                                       iterfield=['realigned_file',
                                                  'noise_mask_file'])
    # Make connections
    compproc.connect(inputspec, 'mean_file',
                     acomp, 'inputspec.mean_file')
    compproc.connect(inputspec, 'reg_file',
                     acomp, 'inputspec.reg_file')
    compproc.connect(inputspec, 'fsaseg_file',
                     acomp, 'inputspec.fsaseg_file')
    compproc.connect(inputspec, 'selector',
                     compcor, 'selector')
    compproc.connect(acomp, ('outputspec.csf_mask',pickfirst),
                     compcor, 'csf_mask_file')
    compproc.connect(inputspec, 'realigned_file',
                     tsnr, 'in_file')
    compproc.connect(inputspec, 'num_components',
                     compcor, 'num_components')
    compproc.connect(inputspec, 'realigned_file',
                     compcor, 'realigned_file')
    compproc.connect(getthresh, 'out_stat',
                     threshold_stddev, 'thresh')
    compproc.connect(threshold_stddev, 'out_file',
                     compcor, 'noise_mask_file')
    compproc.connect(tsnr, 'stddev_file',
                     threshold_stddev, 'in_file')
    compproc.connect(tsnr, 'stddev_file',
                     getthresh, 'in_file')
    compproc.connect(tsnr, 'stddev_file',
                     outputspec, 'stddev_file')
    compproc.connect(tsnr, 'tsnr_file',
                     outputspec, 'tsnr_file')
    compproc.connect(tsnr, 'detrended_file',
                     outputspec, 'tsnr_detrended')
    compproc.connect(compcor, 'noise_components',
                     outputspec, 'noise_components')
    return compproc
