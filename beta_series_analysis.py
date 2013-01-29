"""
BOO YA
"""

import os                                    # system functions
import nipype.algorithms.modelgen as model   # model generation
import nipype.algorithms.rapidart as ra      # artifact detection
import nipype.interfaces.freesurfer as fs    # freesurfer
import nipype.interfaces.fsl as fsl          # fsl
import nipype.interfaces.io as nio           # i/o routines
import nipype.interfaces.matlab as mlab      # how to run matlab
import nipype.interfaces.spm as spm          # spm
import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine
import nipype.interfaces.ants as ants
from nipype.utils.filemanip import loadflat     # some useful stuff for debugging
import scipy.io as sio
import numpy as np
from nipype.interfaces.base import Bunch
from copy import deepcopy
import sys
from nipype import config


qsubargs = '-V -q max30 -e /dev/null -o /dev/null'

config.set('execution','keep_inputs','true')
config.set('execution','remove_unnecessary_outputs','false')

#Initialize paradigm-specific paths
study_dir = os.path.abspath('/mindhive/gablab/GATES/Analysis/beta_series_analysis/')
data_dir = os.path.abspath('/mindhive/gablab/GATES/data/')
preproc_tmpl = 'Analysis/beta_series_preprocessing/%s/preproc/output/fwhm_5.0/*_r[0-9][0-9]_bandpassed.nii.gz'
aparc_tmpl = 'data/%s/mri/aparc+aseg.mgz'
onsets_dir = os.path.join(data_dir, 'onsets/WM/')
subjects_dir = os.path.abspath('/mindhive/xnat/surfaces/GATES/')
os.environ["SUBJECTS_DIR"] = subjects_dir
#Set freesurfer default subj directory
fs.FSCommand.set_default_subjects_dir(subjects_dir)


# Set the way matlab should be called
mlab.MatlabCommand.set_default_matlab_cmd("/software/matlab_versions/2010b/bin/matlab -nodesktop -nosplash")
# If SPM is not in your MATLAB path you should add it here
mlab.MatlabCommand.set_default_paths('/software/spm8')
# Set up how FSL should write nifti files:5       Firefox_wallpaper.png  pypeline
fsl.FSLCommand.set_default_output_type('NIFTI')


#subject_list = ['401','402','403','404']#
subject_list=['300','301','302','303','304','305','306','307','308','309','310','311','312','313','314','315']

# Infosource Node (get subject specific run information)
infosource = pe.Node(interface=util.IdentityInterface(fields=['subject_id']), name="infosource")
infosource.iterables = ('subject_id', subject_list)



# Datasource node
datasource = pe.Node(interface=nio.DataGrabber(infields = ['subject_id'],outfields = ['func','brainmask','meanfunc','reg_file']), name = 'datasource')
datasource.inputs.base_directory = '/mindhive/gablab/GATES'
datasource.inputs.template = '*'
datasource.inputs.field_template = dict(func=preproc_tmpl, brainmask='data/%s/mri/brainmask.mgz',meanfunc='Analysis/beta_series_preprocessing/%s/preproc/mean/*', reg_file='Analysis/beta_series_preprocessing/%s/preproc/bbreg/*_register.dat')
datasource.inputs.template_args = dict(func=[['subject_id']], brainmask = [['subject_id']], meanfunc = [['subject_id']],reg_file=[['subject_id']])



def timeseries_design(subject_id,whatParadigm,onsets_dir):
    import scipy.signal
    import scipy.special as sp
    import numpy as np
    import math
    from nipype.interfaces.base import Bunch
    from copy import deepcopy
    from scipy.io.matlab import loadmat
    import glob
    import os
    #from Facematch import onsets_dir
    print "Entered timeseries_design once with arguments SUBID = "+subject_id+", paradigm = "+whatParadigm+", and onsets dir = "+onsets_dir+"."
    output = []
    regressor_names = None
    regressors = None
    onsets_temp = os.path.join(onsets_dir, subject_id+'*onsets.mat')
    onsets_files = sorted(glob.glob(onsets_temp))
    testmat = loadmat(onsets_files[0], struct_as_record=False)
    testnames = testmat['names'][0]
    names_count_vec = np.zeros(len(testnames))

    for r in range(len(onsets_files)):
        mat = loadmat(onsets_files[r], struct_as_record=False)
        ons = mat['onsets'][0]
        nam = mat['names'][0]
        dur = mat['durations'][0]

        names = []
        durations = []
        run_onsets = []
        for condition in range(len(nam)):
            for onset in range(len(ons[condition][0])): 
                names_count_vec[condition] += 1          
                names.append(str(nam[condition][0])+'_%d'%(names_count_vec[condition]))


                run_onsets.append([ons[condition][0][onset]])
                durations.append(dur[condition][0])
  

        print run_onsets
        print names
        print durations
        output.insert(r,
            Bunch(conditions=deepcopy(names),
                onsets=deepcopy(run_onsets),
                durations=deepcopy(durations),
                amplitudes=None,
                tmod=None,
                pmod=None,
                regressor_names=None,
                regressors=regressors)) #here is where we can do linear, quad, etc detrending
        
    return output



def sub_names(subject_id,whatParadigm,onsets_dir):
    import scipy.signal
    import scipy.special as sp
    import numpy as np
    import math
    from nipype.interfaces.base import Bunch
    from copy import deepcopy
    from scipy.io.matlab import loadmat
    import glob
    import os
    #from Facematch import onsets_dir
    print "Entered sub_names once with arguments SUBID = "+subject_id+", paradigm = "+whatParadigm+", and onsets dir = "+onsets_dir+"."

    onsets_temp = os.path.join(onsets_dir, subject_id+'*onsets.mat')
    onsets_files = sorted(glob.glob(onsets_temp))

    subs = []


    testmat = loadmat(onsets_files[0], struct_as_record=False)
    testnames = testmat['names'][0]
    names_count_vec = np.zeros(len(testnames))

    for r in range(len(onsets_files)):
        mat = loadmat(onsets_files[r], struct_as_record=False)
        ons = mat['onsets'][0]
        nam = mat['names'][0]
        dur = mat['durations'][0]
        names = []
        durations = []
        run_onsets = []
        for condition in range(len(nam)):
            for onset in range(len(ons[condition][0])):             
                names_count_vec[condition] += 1     
                names.append(str(nam[condition][0])+'_%d'%(names_count_vec[condition]))

                run_onsets.append([ons[condition][0][onset]])
                durations.append(dur[condition][0])
                subs.append(('_estimate_model%d/pe%d.nii'%(r,condition*len(ons[condition][0])+onset+1),str(nam[condition][0])+'_%04d.nii'%(names_count_vec[condition])))


    return subs


func="""
def gunzip(compressed_files):
    import gzip
    import string
    import os
    decompressed_files=[]
    print compressed_files
    for i in range(len(compressed_files)):
        print str(i)
        r_file = gzip.GzipFile(compressed_files[i],'r')
        decompressed_files.append(string.rstrip(compressed_files[i],'.gz'))
        w_file = open(decompressed_files[i],'w')
        w_file.write(r_file.read())
        w_file.close()
        r_file.close()
    print decompressed_files
    return decompressed_files
"""
gunzip = pe.Node(util.Function(input_names=['compressed_files'],output_names=['decompressed_files']),name='gunzip')
gunzip.inputs.function_str = func


# Model Specification (NiPype) Node
modelspec = pe.Node(interface=model.SpecifyModel(), name="modelspec")
modelspec.inputs.input_units = 'secs'
modelspec.inputs.time_repetition = 2.0
#128 in original script
modelspec.inputs.high_pass_filter_cutoff = 0.0001 #128#because of linear / quad regressors - otherwise ~160



level1design = pe.Node(interface=fsl.Level1Design(), 
                       name="level1design")
level1design.inputs.interscan_interval = 2.0
level1design.inputs.bases = {'dgamma':{'derivs':False}}
level1design.inputs.model_serial_correlations = False

modelgen = pe.MapNode(interface=fsl.FEATModel(), 
                      name='generate_model',
                      iterfield = ['fsf_file', 
                                   'ev_files'])
    
modelestimate = pe.MapNode(interface=fsl.FILMGLS(smooth_autocorr=True,
                                                 mask_size=5),
                           name='estimate_model',
                           iterfield = ['design_file','in_file'])




# Datasink node for saving output of the pipeline
betasink = pe.Node(interface=nio.DataSink(), name="betasink")
betasink.inputs.base_directory = os.path.join(study_dir, 'betasink')


pickfirst = lambda x: x[0]
beta_series_pipeline = pe.Workflow(name="beta_series_pipeline")
beta_series_pipeline.base_dir = os.path.abspath('/mindhive/scratch/jsalva/GATES/working/beta_series/')
beta_series_pipeline.connect([
    (infosource,datasource,[('subject_id','subject_id')]),
    (infosource,betasink,[(('subject_id', sub_names, 'WM',onsets_dir),'substitutions')]),

    (datasource,modelspec,[('func','functional_runs')]),
    (infosource,modelspec, [(('subject_id',timeseries_design,'WM',onsets_dir),'subject_info')]),
    (modelspec,level1design,[('session_info','session_info')]),
    (level1design,modelgen,[('fsf_files','fsf_file'),
                            ('ev_files','ev_files')]),
    (modelgen, modelestimate,[('design_file','design_file')]),
    (datasource,modelestimate,[('func','in_file')]),
    (modelestimate, betasink,[('param_estimates','betas.@param_estimates'),
                                  ('sigmasquareds','sigmasquared.@sigmasquareds'),
                                  #('corrections','corrections.@corrections'),
                                  ('dof_file','dof.@dof_file')])
    ])


#beta_series_pipeline.write_graph(graph2use="exec");



beta_series_pipeline.run()#plugin='PBS',plugin_args={'qsub_args':qsubargs})


#####DUMMY L2######


def getsubs(subject):
    subs = []
    for i in range(10):
        for j in range(10):
            subs.append(('conn/_subject_id_%s/_calc_fischer_z_img%d%d'%(subject,i,j),'%s'%(subject)))
    print subs
    #exit()
    return subs




def corr_image(resting_image, aparc_aseg_file, seed_region):
    import numpy as np
    import nibabel as nb
    import matplotlib.pyplot as plt
    from surfer import Brain, Surface
    import os
    import string
    import math
    LUT=dict()
    for line in open("/software/Freesurfer/current/FreeSurferColorLUT.txt"):
        if not "#" in line:
            if len(line.split()) == 6:
                LUT[int(line.split()[0])] = line.split()[1]

    aparc_aseg = nb.load(aparc_aseg_file)
    indeces=np.where(aparc_aseg.get_data()==seed_region)
    img = nb.load(resting_image)
    data=img.get_data()
    header=aparc_aseg.get_header()
    affine=aparc_aseg.get_affine()
    corr_data = np.zeros(aparc_aseg.shape[0:3])
    sum_signal = np.zeros(data.shape[3])
    for idx in range(len(indeces[0])):
        sum_signal = sum_signal + data[indeces[0][idx]][indeces[1][idx]][indeces[2][idx]]

    seed_signal = sum_signal/(len(indeces[0]))

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            for k in range(data.shape[2]):
                val = np.arctanh(np.corrcoef(seed_signal,data[i][j][k])[0,1])*np.sqrt(len(data[i][j][k])-3)
                corr_data[i][j][k] = val if not math.isnan(val) else 0
    filename = os.path.split(resting_image)[1]
    out_img = nb.Nifti1Image(corr_data,affine,header)
    corr_image = os.path.abspath('%s_z_corr_'%LUT[seed_region]+filename)
    out_img.to_filename(corr_image) 
    return corr_image

inputnode = pe.Node(util.IdentityInterface(fields=['subject_id']),
                        name='inputnode')
inputnode.iterables = ('subject_id', subject_list)

betagrabber = pe.Node(interface=nio.DataGrabber(infields = ['subject_id'], outfields = ['inst_DRY','inst_4RY', 'stim_4RY','stim_DRY','stim_2R']), name = 'betagrabber')
betagrabber.inputs.sort_filelist = True
betagrabber.inputs.base_directory = os.path.join(study_dir, 'betasink')
betagrabber.inputs.template = '*'
betagrabber.inputs.field_template = dict(inst_DRY='betas/_subject_id_%s/inst_DRY*.nii', inst_4RY='betas/_subject_id_%s/inst_4RY*.nii', stim_2R='betas/_subject_id_%s/stim_2R*.nii', stim_4RY='betas/_subject_id_%s/stim_4RY*.nii',stim_DRY='betas/_subject_id_%s/stim_DRY*.nii')
betagrabber.inputs.template_args = dict(inst_DRY=[['subject_id']], inst_4RY = [['subject_id']],stim_2R=[['subject_id']],stim_4RY=[['subject_id']],stim_DRY=[['subject_id']])

static_grabber = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],outfields=['aparc_aseg','fsl_mat','meanfunc']),name='static_grabber')
static_grabber.inputs.base_directory = '/mindhive/gablab/GATES'
static_grabber.inputs.template='*'
static_grabber.inputs.field_template = dict(aparc_aseg='data/%s/mri/aparc+aseg.mgz',fsl_mat='Analysis/WM/l1output/%s/norm_anat/fsl2antsAffine.txt',meanfunc = 'Analysis/beta_series_preprocessing/%s/preproc/mean/*mean*')
static_grabber.inputs.template_args = dict(aparc_aseg = [['subject_id']], fsl_mat=[['subject_id']],meanfunc=[['subject_id']])


#converts freesurfer segmentation into nii
aparcaseg_2nii = pe.Node(
    fs.preprocess.MRIConvert(),
    name='aparcaseg_2nii')
aparcaseg_2nii.inputs.out_type = 'nii'

merge1 = pe.Node(fsl.Merge(),name='merge1')
merge1.inputs.dimension = 't'
merge1.inputs.merged_file = 'bseries_inst_DRY.nii'

merge2 = pe.Node(fsl.Merge(),name='merge2')
merge2.inputs.dimension = 't'
merge2.inputs.merged_file = 'bseries_inst_4RY.nii'

merge3 = pe.Node(fsl.Merge(),name='merge3')
merge3.inputs.dimension = 't'
merge3.inputs.merged_file = 'bseries_stim_2R.nii'

merge4 = pe.Node(fsl.Merge(),name='merge4')
merge4.inputs.merged_file = 'bseries_stim_DRY.nii'
merge4.inputs.dimension = 't'

merge5 = pe.Node(fsl.Merge(),name='merge5')
merge5.inputs.merged_file = 'bseries_stim_4RY.nii'
merge5.inputs.dimension = 't'

threeDreg = pe.Node(ants.WarpTimeSeriesImageMultiTransform(),name='threeDreg')
threeDreg.inputs.invert_affine = [1]
threeDreg.inputs.use_nearest = True
threeDreg.inputs.dimension = 3
tolist = lambda x: [x]



calc_fischer_z_img1 = pe.MapNode(util.Function(input_names=['resting_image','aparc_aseg_file','seed_region'],
                               output_names=["corr_image"],
                               function=corr_image),
                    name="calc_fischer_z_img1",iterfield=['seed_region'])
calc_fischer_z_img1.inputs.seed_region = [11,12,13,1003,1027,50,51,52,2003,2027]
calc_fischer_z_img2 = calc_fischer_z_img1.clone(name='calc_fischer_z_img2')
calc_fischer_z_img3 = calc_fischer_z_img1.clone(name='calc_fischer_z_img3')
calc_fischer_z_img4 = calc_fischer_z_img1.clone(name='calc_fischer_z_img4')
calc_fischer_z_img5 = calc_fischer_z_img1.clone(name='calc_fischer_z_img5')

beta_series_sinker =  pe.Node(nio.DataSink(), name='beta_series_sinker')

beta_series_sinker.inputs.base_directory = os.path.join(study_dir, 'beta_series_sinker')

beta_series_pipeline2 = pe.Workflow(name="beta_series_pipeline2")
beta_series_pipeline2.base_dir = os.path.abspath('/mindhive/scratch/jsalva/GATES/working/beta_series_4d/')
beta_series_pipeline2.connect([
    (inputnode,betagrabber,[('subject_id','subject_id')]),  
    (inputnode,static_grabber,[('subject_id','subject_id')]), 
    (inputnode,beta_series_sinker,[(('subject_id', getsubs),'substitutions')]), 
    (betagrabber,merge1,[('inst_DRY','in_files')]),
    (betagrabber,merge2,[('inst_4RY','in_files')]),
    (betagrabber,merge3,[('stim_2R','in_files')]),
    (betagrabber,merge4,[('stim_DRY','in_files')]),
    (betagrabber,merge5,[('stim_4RY','in_files')]),
    (merge1,beta_series_sinker,[('merged_file','inst_DRY.beta_series.@img')]),
    (merge2,beta_series_sinker,[('merged_file','inst_4RY.beta_series.@img')]),
    (merge3,beta_series_sinker,[('merged_file','stim_2R.beta_series.@img')]),
    (merge4,beta_series_sinker,[('merged_file','stim_DRY.beta_series.@img')]),
    (merge5,beta_series_sinker,[('merged_file','stim_4RY.beta_series.@img')]),

    (static_grabber,aparcaseg_2nii,[('aparc_aseg','in_file')]),
    (aparcaseg_2nii, threeDreg,[('out_file','moving_image')]),

    (static_grabber,threeDreg,[('fsl_mat','transformation_series')]),

    (static_grabber,threeDreg,[('meanfunc','reference_image')]),

    (threeDreg,beta_series_sinker,[('output_image','segment.@aparcaseg')]),


    (threeDreg,calc_fischer_z_img1,[('output_image','aparc_aseg_file')]),
    (threeDreg,calc_fischer_z_img2,[('output_image','aparc_aseg_file')]),
    (threeDreg,calc_fischer_z_img3,[('output_image','aparc_aseg_file')]),
    (threeDreg,calc_fischer_z_img4,[('output_image','aparc_aseg_file')]),
    (threeDreg,calc_fischer_z_img5,[('output_image','aparc_aseg_file')]),


    (merge1,calc_fischer_z_img1,[('merged_file','resting_image')]),
    (merge2,calc_fischer_z_img2,[('merged_file','resting_image')]),
    (merge3,calc_fischer_z_img3,[('merged_file','resting_image')]),
    (merge4,calc_fischer_z_img4,[('merged_file','resting_image')]),
    (merge5,calc_fischer_z_img5,[('merged_file','resting_image')]),





    (calc_fischer_z_img1,beta_series_sinker,[('corr_image','inst_DRY.conn.@img')]),
    (calc_fischer_z_img2,beta_series_sinker,[('corr_image','inst_4RY.conn.@img')]),
    (calc_fischer_z_img3,beta_series_sinker,[('corr_image','stim_2R.conn.@img')]),
    (calc_fischer_z_img4,beta_series_sinker,[('corr_image','stim_DRY.conn.@img')]),
    (calc_fischer_z_img5,beta_series_sinker,[('corr_image','stim_4RY.conn.@img')])

])

beta_series_pipeline2.run()#plugin='PBS',plugin_args={'qsub_args':qsubargs})








   
