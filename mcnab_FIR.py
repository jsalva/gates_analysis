import numpy as np
import matplotlib as mpl
import scipy as sp
import nipype.interfaces.io as nio
import nipype.pipeline.engine as pe 
import nipype.algorithms.modelgen as model
import nipype.interfaces.matlab as mlab
import nipype.interfaces.spm as spm
import nipype.interfaces.ants as ants
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as util
from nipype import config
from subject_info import *
import os
from fir_utils import *
study_dir = os.path.abspath('/mindhive/gablab/GATES/Analysis/MCNAB/')
data_dir = os.path.abspath('/mindhive/gablab/GATES/data/')
onsets_dir = os.path.join(data_dir, 'onsets/MCNAB/')
subjects_dir = os.path.abspath('/mindhive/xnat/surfaces/GATES/')
os.environ["SUBJECTS_DIR"] = subjects_dir
fs.FSCommand.set_default_subjects_dir(subjects_dir)
mlab.MatlabCommand.set_default_matlab_cmd("/software/matlab_versions/2010b/bin/matlab -nodesktop -nosplash")
mlab.MatlabCommand.set_default_paths('/software/spm8_4290')
config.set('execution','keep_inputs','true')
config.set('execution','remove_unnecessary_outputs','false')


def gen_peak_mask(anat_mask,func_activation):
    import nibabel as nb 
    import numpy as np
    from fir_utils import (binarize_peak, dilate_mask)
    import sys
    import os

    anatmask = nb.load(anat_mask)
    mask_data = np.asarray(anatmask.get_data())
    mask_header = anatmask.get_header()
    mask_affine = anatmask.get_affine()

    func = nb.load(func_activation)
    func_data = np.asarray(func.get_data())
    func_header = func.get_header()
    func_affine = func.get_affine()

    if np.array(func_data).shape != np.array(mask_data).shape:
        print np.array(func_data).shape
        print np.array(mask_data).shape
        sys.exit("functional and mask do not have the same 3d shape!")

    masked_activation = np.zeros(func_data.shape)
    masked_activation[np.nonzero(mask_data)] = func_data[np.nonzero(mask_data)]

    peak_mask = binarize_peak(masked_activation)
    dilated_peak_mask = dilate_mask(peak_mask)

    path = '%s/%s.nii'%(os.path.abspath(os.path.curdir),'peak_mask')

    out_image = nb.Nifti1Image(dilated_peak_mask,mask_affine,mask_header)
    out_image.to_filename(path)

    return path

def get_timeseries_from_mask(func, mask):
    """extracts the average timeseries from the mask"""
    import nibabel as nb
    import numpy as np
    import sys
    import os

    func_nb = nb.load(func)
    func_data = np.asarray(func_nb.get_data())
    func_affine = func_nb.get_affine()
    func_header = func_nb.get_header()

    mask_nb = nb.load(mask)
    mask_data = np.asarray(mask_nb.get_data())
    mask_affine = mask_nb.get_affine()
    mask_header = mask_nb.get_header()

    if np.array(func_data).shape[:3] != np.array(mask_data).shape:
        print np.array(func_data).shape[:3]
        print np.array(mask_data).shape
        sys.exit("functional and mask do not have the same 3d shape!")
    
    timeseries_matrix = func_data[np.nonzero(mask_data)]
    
    average_timeseries = np.average(timeseries_matrix,axis=0)

    if len(average_timeseries) != np.array(func_data).shape[3]:
        sys.exit('timeseries extracted from mask does not match that of the input functional!')

    path = '%s/%s.txt'%(os.path.abspath(os.path.curdir),os.path.split(func)[1].split('.')[0])
    np.savetxt(path,average_timeseries)

    return path


def extract_segmentation_roi(func_aparc_aseg, segmentation_value):
    """binarizes the volume at locations equal to segmentation_value"""
    import nibabel as nb
    import numpy as np
    import os
    from fir_utils import (functional_binarize, dilate_mask)
    import sys

    segmentation = nb.load(func_aparc_aseg)
    seg_data = segmentation.get_data()
    header = segmentation.get_header()
    affine = segmentation.get_affine()

    seg_mask = functional_binarize(seg_data,lambda x: x == segmentation_value)

    path = '%s/%s.nii'%(os.path.abspath(os.path.curdir),'mask_'+str(segmentation_value))

    dilated_mask = dilate_mask(np.array(seg_mask))

    out_image = nb.Nifti1Image(dilated_mask,affine,header)
    out_image.to_filename(path)

    return path


def model_fir_from_onsets(subject_id, modeled_event_name,timeseries_files, onsets_files, motion_param_files, art_outlier_files):
    from scipy.io.matlab import loadmat
    from scipy.stats import gamma
    from subject_info import info
    import numpy as np
    import sys
    from scipy.misc import imsave
    from scipy.linalg import toeplitz
    from scipy.interpolate import interp1d
    import os
    import re

    #defaults taken from SPM's spm_hrf.m
    repitition_time = 2.0
    hrf_duration = 32.0
    time_bins_per_scan = 16.0
    dt = repitition_time/time_bins_per_scan

    
    hrf_domain = np.arange(0,hrf_duration,dt)
    hrf = gamma(6.).pdf(hrf_domain) - gamma(16.).pdf(hrf_domain)/6
    normalized_hrf = hrf/np.sum(hrf)

    WM_runs = info[subject_id][1][0]
    num_runs = len(WM_runs)

    pattern = '('
    for run in range(num_runs):
        pattern += WM_runs[run]+'|'
    pattern = pattern.strip('|')
    pattern += ')'
    timeseries_files = [file for file in timeseries_files if re.search(pattern, file)]
    motion_param_files = [file for file in motion_param_files if re.search(pattern, file)]
    art_outlier_files = [file for file in art_outlier_files if re.search(pattern, file)]
    #onsets_files = [file for file in onsets_files if re.search(pattern,file)]
    timeseries_files.sort()
    onsets_files.sort()
    motion_param_files.sort()
    art_outlier_files.sort()



    if not num_runs == len(timeseries_files):
        sys.exit('%s: number of timeseries files does not match number of runs!'%(subject_id))
    if not num_runs == len(onsets_files):
        sys.exit('%s: number of onsets files does not match number of runs!'%(subject_id))
    if not num_runs == len(motion_param_files):
        sys.exit('%s: number of motion parameter files does not match number of runs!'%(subject_id))
    if not num_runs == len(art_outlier_files):
        sys.exit('%s: number of outlier files does not match number of runs!'%(subject_id))

    
    tmp_design = []
    names = []
    tmp_signal = []
    for run in range(num_runs):        

        timeseries = np.genfromtxt(timeseries_files[run])
        
        tmp_signal.append(np.array([]))
        interp_timeseries = interp1d(np.arange(0,len(timeseries)*repitition_time,2),timeseries,kind='linear',axis=0,bounds_error=False,fill_value=0)
        tmp_signal[run] = interp_timeseries(np.arange(0,len(timeseries)*repitition_time,dt))        

        try:
            outlier_timepoints = np.genfromtxt(art_outlier_files[run])
        except IOError:
            print "no timepoints found!"
            outlier_timepoints = np.array([])

        #one thing to check will be whether the motion parameter files are absolute in the sense that they are rigid body transforms from volume 1,
        #or whether they are differentiated w/ respect to the previous volume
        motion_parameters = np.genfromtxt(motion_param_files[run])
        mat = loadmat(onsets_files[run])
        
        nam = mat['names'][0]
        ons = mat['onsets'][0]
        dur = mat['durations'][0]
        pmods = mat['pmod'][0]


        tmp_design.append(np.array([]))
        run_names = ['intercept_run%d'%(run)]
        
        tmp_design[run] = np.ones(len(timeseries)*time_bins_per_scan)

        for col_idx,event_name in enumerate(nam):
            tmp_col = None
            tmp_nam = None
            tmp_param_col = None

            if event_name[0] == modeled_event_name:

                tmp_col = np.zeros(len(timeseries)*time_bins_per_scan)

                for ons_idx,onset in enumerate(ons[col_idx][0]):
                    tmp_col[int(np.round((onset-2*repitition_time)/dt))] = 1

                    
                num_shifts = len(np.arange(-2*repitition_time,10*repitition_time,dt))
                tmp_nam = []
                for fir_idx,time in enumerate(np.arange(-2*repitition_time,10*repitition_time,dt)):
                    tmp_nam.append(event_name[0]+'_fir_%d@%.2fs'%(fir_idx,time))
                tmp_design[run] = np.column_stack([tmp_design[run],toeplitz(tmp_col,np.zeros(num_shifts))])

                for col_name in tmp_nam:
                    run_names.append(col_name+'_run%d'%(run))
                

            else:            
                tmp_col = np.zeros(len(timeseries)*time_bins_per_scan)
                tmp_nam = [event_name[0]]

                if pmods[col_idx].any():
                    tmp_param_col = np.zeros(len(timeseries)*time_bins_per_scan)
                    #using the name "load" is specific to our Working Memory paradigm. In the more general case, one should use a name that makes sense for whatever is being parametrically
                    #modulated. 
                    tmp_nam.append(event_name[0]+'xload^1')
                
                for ons_idx,onset in enumerate(ons[col_idx][0]):
                    tmp_col[int(np.round(onset/dt)):int(np.round(onset/dt+dur[col_idx][0][ons_idx]/dt))] = 1
                    if pmods[col_idx].any():
                        tmp_param_col[int(np.round(onset/dt)):int(np.round(onset/dt+dur[col_idx][0][ons_idx]/dt))] = pmods[col_idx][0][ons_idx]

                for col_name in tmp_nam:
                    run_names.append(col_name+'_run%d'%(run))
                #using "Same" may introduce boundary effects for convolution. Perhaps we should guarantee the convolution is only defined for complete overlaps and fill in boundaries with zeros?
                tmp_design[run] = np.column_stack([tmp_design[run],np.convolve(tmp_col,normalized_hrf,'same')])
                if not tmp_param_col is None:
                    #Not sure whether parametric modulators should be concolved. Probably not! They should modulate the associated block of bold uniformly across the duraiton of the HRF, no?
                    tmp_design[run] = np.column_stack([tmp_design[run],np.convolve(tmp_param_col,normalized_hrf,'same')])

        #done with "modeled" events; move on to movement, outliers, and lagrange
        for realign_param in ["x_trans","y_trans","z_trans","pitch","roll","yaw"]:
            run_names.append('realign_%s_run%d'%(realign_param,run))
        interp_motion_func = interp1d(np.arange(0,len(timeseries)*repitition_time,2),motion_parameters,kind='linear',axis=0,bounds_error=False,fill_value=0)
        tmp_design[run] = np.column_stack([tmp_design[run],interp_motion_func(np.arange(0,len(timeseries)*repitition_time,dt))])

        if outlier_timepoints.size>1:
            for outlier_idx, outlier_event in enumerate(outlier_timepoints):
                run_names.append('art_outlier%d_run%d'%(outlier_idx,run))
                tmp_col = np.zeros(len(timeseries)*time_bins_per_scan)
                tmp_col[int(np.round(outlier_event/dt)):int(np.round((outlier_event+1)/dt))-1] = 1
                tmp_design[run] = np.column_stack([tmp_design[run],tmp_col])
        elif outlier_timepoints.size == 1:
            run_names.append('art_outlier%d_run%d'%(1,run))
            tmp_col = np.zeros(len(timeseries)*time_bins_per_scan)
            tmp_col[int(np.round(outlier_timepoints/dt)):int(np.round((outlier_timepoints+1)/dt))] = 1
            tmp_design[run] = np.column_stack([tmp_design[run],tmp_col])

        for exp in range(5)[1:]:
            coeffs = np.zeros(exp).tolist()
            coeffs.append(1)
            poly = np.polynomial.legendre.Legendre(coeffs)
            tmp_col = poly.linspace(len(timeseries)*time_bins_per_scan)[1]
            tmp_design[run] = np.column_stack([tmp_design[run],tmp_col])
            run_names.append('legendre^%d_run%d'%(exp,run))
        
        for name in run_names:
            names.append(name)

    rows = 0
    columns = 0
    timeseries_length = 0
    for run_index,run_design in enumerate(tmp_design):
        print "design for run %d has (%d,%d) shape"%(run_index,run_design.shape[0],run_design.shape[1])
        rows+=run_design.shape[0]
        columns+=run_design.shape[1]
        timeseries_length += len(tmp_signal[run_index])

    curr_row = 0
    curr_col = 0
    print rows
    print columns
    design_matrix = np.zeros((rows,columns))
    signal_timeseries = np.zeros(timeseries_length)
    print design_matrix.shape
    for run_index,run_design in enumerate(tmp_design):
        print "made it to design for run %d"%(run_index)
        design_matrix[curr_row:curr_row+run_design.shape[0],curr_col:curr_col+run_design.shape[1]] = run_design
        signal_timeseries[curr_row:curr_row+len(tmp_signal[run_index])] = tmp_signal[run_index]
        curr_row += run_design.shape[0]
        curr_col += run_design.shape[1]

    design_path = '%s/design.txt'%(os.path.abspath(os.path.curdir))
    np.savetxt(design_path,design_matrix)

    names_path = '%s/names.txt'%(os.path.abspath(os.path.curdir))
    name_file = open(names_path,'w')
    for name in names:
        print >> name_file,name
    name_file.close()

    betas = np.dot(np.linalg.inv(np.dot(np.transpose(design_matrix),design_matrix)),np.dot(np.transpose(design_matrix),signal_timeseries))
    betas_path = '%s/betas.txt'%(os.path.abspath(os.path.curdir))
    np.savetxt(betas_path,betas)

    hrf_matrix = []
    for run, _garbage in enumerate(tmp_signal):
        
        hrf_indexes = [idx for idx,name in enumerate(names) if re.search('.*fir.*run%d.*'%(run),name)]
        hrf_matrix.append(betas[hrf_indexes])

    hrf_matrix = np.squeeze(np.array(hrf_matrix))
    hrf_matrix_path = '%s/hrf_by_run.txt'%(os.path.abspath(os.path.curdir))
    np.savetxt(hrf_matrix_path,hrf_matrix)

    hrf = np.average(hrf_matrix,axis = 0)
    interp_hrf = interp1d(np.arange(-2*repitition_time,10*repitition_time,dt),hrf,kind='linear',axis=0,bounds_error=False,fill_value=hrf[-1])
    downsampled_hrf = interp_hrf(np.arange(-2*repitition_time,10*repitition_time,1))
    #U,E,V = np.linalg.svd(hrf); consider using SVD/PCA to get HRF?
    fir_hrf_path = '%s/fir_avg_hrf.txt'%(os.path.abspath(os.path.curdir))
    np.savetxt(fir_hrf_path,downsampled_hrf)

    return fir_hrf_path, hrf_matrix_path, design_path, names_path, betas_path
    
    
        
        

func_aseg = pe.Workflow(name='func_aseg')
func_aseg.base_dir = '/mindhive/scratch/jsalva/mcnab/mcnab_func_aseg/'


infosource = pe.Node(interface = util.IdentityInterface(fields=['subject_id']), name="infosource")
infosource.iterables = ('subject_id', subject_list)


datasource = pe.Node(interface=nio.DataGrabber(
        infields = ['subject_id'],
        outfields = ['meanfunc','aparc_aseg','brainmask','coreg_matrix']),
    name = 'datasource')
datasource.inputs.base_directory = '/mindhive/gablab/GATES/'
datasource.inputs.template = '*'
#MAKE ABSOLUTELY SURE THESE COME OUT SORTED.
datasource.inputs.field_template = dict(
    meanfunc = 'Analysis/MCNAB/l1_mcnab/%s/subj_anat/rWM*.nii',
    aparc_aseg = 'data/%s/mri/aparc+aseg.mgz',
    brainmask = 'data/%s/mri/brainmask.mgz',
    coreg_matrix = 'Analysis/MCNAB/l1_mcnab/%s/surfreg/*.mat'
)
datasource.inputs.template_args = dict(
    meanfunc = [['subject_id']],
    aparc_aseg = [['subject_id']],  
    brainmask = [['subject_id']],
    coreg_matrix = [['subject_id']]
)


func_aseg.connect( [(infosource, datasource, [( 'subject_id' , 'subject_id' )] )] )


aparcaseg_2nii = pe.Node(fs.preprocess.MRIConvert(),
    name = 'aparcaseg_2nii')
aparcaseg_2nii.inputs.out_type = 'nii'


func_aseg.connect( [(datasource, aparcaseg_2nii, [( 'aparc_aseg' , 'in_file' )] )] )


mgz_2nii = pe.Node(fs.preprocess.MRIConvert(),
    name='mgz_2nii')
mgz_2nii.inputs.out_type = 'nii'


func_aseg.connect( [(datasource, mgz_2nii, [( 'brainmask' , 'in_file' )] )] )


fsl_reg_2_itk = pe.Node(util.Function(
        input_names = ['unwarped_brain', 'mean_func', 'out_fsl_file'],
        output_names = ['fsl2antsAffine'],function=convert_affine),
    name = 'fsl_reg_2_itk')


func_aseg.connect( [(mgz_2nii, fsl_reg_2_itk, [( 'out_file' , 'unwarped_brain' )] )] )
func_aseg.connect( [(datasource, fsl_reg_2_itk, [( 'meanfunc' , 'mean_func' )] )] )
func_aseg.connect( [(datasource, fsl_reg_2_itk, [( 'coreg_matrix' , 'out_fsl_file' )] )] )


native2func = pe.Node(ants.WarpImageMultiTransform(),
    name='native2func')
native2func.inputs.invert_affine = [1]
native2func.inputs.use_nearest = True
native2func.inputs.reslice_by_header = True


func_aseg.connect( [(aparcaseg_2nii, native2func, [( 'out_file' , 'moving_image' )] )] )
func_aseg.connect( [(datasource, native2func, [( 'meanfunc' , 'reference_image' )] )] )
func_aseg.connect( [(fsl_reg_2_itk, native2func, [( 'fsl2antsAffine' , 'transformation_series' )] )] )


preproc_sink = pe.Node(nio.DataSink(), name = 'preproc_sink')
preproc_sink.inputs.base_directory = '/mindhive/gablab/GATES/Analysis/WM_func_aparc_aseg/'


func_aseg.connect( [(native2func, preproc_sink, [( 'output_image' , 'func_seg_file' )] )] )

func_aseg.run(plugin = 'MultiProc', plugin_args = {'n_procs' : 30})


segmentation_regions = [11,12,13,50,51,52]
lut = dict()
lut[str(11)] = 'lh_caudate'
lut[str(12)] = 'lh_putamen'
lut[str(13)] = 'lh_pallidum'
lut[str(50)] = 'rh_caudate'
lut[str(51)] = 'rh_putamen'
lut[str(52)] = 'rh_pallidum'


for seg_val in segmentation_regions:

    fir_workflow = None
    fir_infosource = None
    FIR_datagrabber = None
    extract_roi_from_seg_value = None
    mean_ts_at_peak = None
    fir_model = None
    fir_datasink = None


    fir_workflow = pe.Workflow(name='wf_%s'%(lut[str(seg_val)]))
    fir_workflow.base_dir = '/mindhive/scratch/jsalva/mcnab/fir/'


    fir_infosource = pe.Node(interface = util.IdentityInterface(fields=['subject_id']), name="fir_infosource")
    fir_infosource.iterables = ('subject_id', subject_list)

    
    FIR_datagrabber = pe.Node(nio.DataGrabber(
            infields = ['subject_id'],
            outfields = ['realigned_funcs','func_aparc_aseg','filtering_contrast','motion_params','art_outlier_files','onsets']),
        name = 'fir_datagrabber')
    FIR_datagrabber.inputs.template = '*'
    FIR_datagrabber.inputs.base_directory = '/mindhive/gablab/GATES/'
    FIR_datagrabber.inputs.field_template = dict(
        realigned_funcs = 'data/%s/niftis/WM*.nii',
        func_aparc_aseg = 'Analysis/WM_func_aparc_aseg/func_seg_file/_subject_id_%s/*.nii',
        filtering_contrast = 'Analysis/MCNAB/l1_mcnab/%s/subj_contrasts/con_0001.img',
        motion_params = 'Analysis/MCNAB/l1_mcnab/%s/qc_realign/*.txt',
        art_outlier_files = 'Analysis/MCNAB/l1_mcnab/%s/qc_art/art*.txt',
        onsets = 'data/onsets/MCNAB/%s_WMrun*_onsets_mcnab.mat'
    )
    FIR_datagrabber.inputs.template_args = dict(
        realigned_funcs = [['subject_id']],
        func_aparc_aseg = [['subject_id']],
        filtering_contrast = [['subject_id']],
        motion_params = [['subject_id']],
        art_outlier_files = [['subject_id']],
        onsets = [['subject_id']]
    )


    fir_workflow.connect( [(fir_infosource, FIR_datagrabber, [( 'subject_id' , 'subject_id' )] )] )


    extract_roi_from_seg_value = pe.Node(util.Function(
            input_names = ['func_aparc_aseg','segmentation_value'],
            output_names = ['mask'], function = extract_segmentation_roi),
        name = 'extract_roi_from_seg_value')
    extract_roi_from_seg_value.inputs.segmentation_value = seg_val


    fir_workflow.connect( [(FIR_datagrabber, extract_roi_from_seg_value, [( 'func_aparc_aseg' , 'func_aparc_aseg' )] )] )


    find_peak_within_anat_mask = pe.Node(util.Function(
            input_names = ['anat_mask','func_activation'],
            output_names = ['peak_mask'], function = gen_peak_mask),
        name = 'find_peak_within_anat_mask')


    fir_workflow.connect( [(FIR_datagrabber, find_peak_within_anat_mask, [( 'filtering_contrast' , 'func_activation' )] )] )
    fir_workflow.connect( [(extract_roi_from_seg_value, find_peak_within_anat_mask, [( 'mask' , 'anat_mask' )] )] )


    mean_ts_at_peak = pe.MapNode(util.Function(
            input_names = ['func','mask'],
            output_names = ['mean_signal'], function = get_timeseries_from_mask),
        name = 'mean_ts_at_peak',
        iterfield = ['func'])


    fir_workflow.connect( [(find_peak_within_anat_mask, mean_ts_at_peak, [( 'peak_mask' , 'mask' )] )] )
    fir_workflow.connect( [(FIR_datagrabber, mean_ts_at_peak, [( 'realigned_funcs' , 'func' )] )] )


    fir_model = pe.MapNode(util.Function(
            input_names = ['subject_id','modeled_event_name','timeseries_files', 'onsets_files', 'motion_param_files', 'art_outlier_files'],
            output_names = ['HRF','HRF_matrix','design_matrix','col_names','betas'],
            function = model_fir_from_onsets),
        name = 'fir_model',
        iterfield = ['modeled_event_name'])
    fir_model.inputs.modeled_event_name = ['inst_dist','inst_nodist']


    fir_workflow.connect( [( mean_ts_at_peak, fir_model, [( 'mean_signal' , 'timeseries_files' )] )] )
    fir_workflow.connect( [( fir_infosource, fir_model, [( 'subject_id' , 'subject_id' )] )] )
    fir_workflow.connect( [( FIR_datagrabber, fir_model, [( 'onsets' , 'onsets_files' )] )] )
    fir_workflow.connect( [( FIR_datagrabber, fir_model, [( 'motion_params' , 'motion_param_files' )] )] )
    fir_workflow.connect( [( FIR_datagrabber, fir_model, [( 'art_outlier_files' , 'art_outlier_files' )] )] ) 

    fir_datasink = pe.Node(nio.DataSink(), name="fir_datasink")
    fir_datasink.inputs.base_directory = '/mindhive/gablab/GATES/Analysis/MCNAB_FIR/%s'%(lut[str(seg_val)])
    fir_datasink.inputs.substitutions = [('_subject_id_%s/'%(subject_id),'%s/'%(subject_id)) for subject_id in subject_list]

    fir_workflow.connect( [( fir_model, fir_datasink, [( 'HRF' , 'fir_model.@hrf' )] )] )
    fir_workflow.connect( [( fir_model, fir_datasink, [( 'HRF_matrix' , 'fir_model.@hrf_mat' )] )] )
    fir_workflow.connect( [( fir_model, fir_datasink, [( 'col_names' , 'fir_model.@names' )] )] )
    fir_workflow.connect( [( fir_model, fir_datasink, [( 'betas' , 'fir_model.@betas' )] )] )


    fir_workflow.run(plugin = 'MultiProc', plugin_args = {'n_procs' : 30})







