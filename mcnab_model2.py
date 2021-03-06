"""
Facematch paradigm analysis pipeline
Written for CARD study by John Salvatore
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
from nipype.utils.filemanip import loadflat     # some useful stuff for debugging
import scipy.io as sio
import numpy as np
from nipype.interfaces.base import Bunch
from copy import deepcopy
import sys
from nipype import config
from compcor_workflow import (create_compcorr,extract_noise_components)
from subject_info import (info,subject_list)

config.set('execution','keep_inputs','true')
config.set('execution','remove_unnecessary_outputs','false')
##############################################################################
#                               ARGUMENTS
##############################################################################
#perhaps implement with argparse
print sys.argv
if not len(sys.argv) == 4:
    sys.stderr.write("The paradigm and analysis level must be provided on the command line")
    sys.exit("The paradigm and analysis level must be provided on the command line")
elif sys.argv[1] in ['WMSTAT','WM','Nback_spatial','Nback_letters','MCNAB','MCNAB2']:
    #originally 'condition'
    paradigm = sys.argv[1]
    if sys.argv[2] in ['l1','l2']:
        levelToRun = sys.argv[2]
    else:
        sys.stderr.write("You must indicate whether this is l1 or l2 on the command line")
        sys.exit("You must indicate whether this is l1 or l2 on the command line")        
    if sys.argv[3] in ['torque','PBS']:
        numCores = "PBS"
    else:
        sys.stderr.write("Unacceptable number of processor cores specified, running in series instead...")
        numCores = 5
else:
    sys.exit("Paradigm must be one of the following: 'WM','WMSTAT','Nback_letters','MCNAB','MCNAB2'")

##############################################################################
#                      Participant Info and Directories
##############################################################################
#Initialize paradigm-specific paths
study_dir = os.path.abspath('/mindhive/gablab/GATES/Analysis/%s/'%(paradigm))
data_dir = os.path.abspath('/mindhive/gablab/GATES/data/')
onsets_dir = os.path.join(data_dir, 'onsets/%s/'%(paradigm))
subjects_dir = os.path.abspath('/mindhive/xnat/surfaces/GATES/')
os.environ["SUBJECTS_DIR"] = subjects_dir
#Set freesurfer default subj directory
fs.FSCommand.set_default_subjects_dir(subjects_dir)

#Test on subject that has been unpacked and reconned for WMSTAT


#Missing onsets:'325','326','327','328'

# Set the way matlab should be called
mlab.MatlabCommand.set_default_matlab_cmd("/software/matlab_versions/2010b/bin/matlab -nodesktop -nosplash")
# If SPM is not in your MATLAB path you should add it here
mlab.MatlabCommand.set_default_paths('/software/spm8_4290')
# Set up how FSL should write nifti files:5       Firefox_wallpaper.png  pypeline
fsl.FSLCommand.set_default_output_type('NIFTI')

# Infosource Node (get subject specific run information)
infosource = pe.Node(interface=util.IdentityInterface(fields=['subject_id']), name="infosource")
infosource.iterables = ('subject_id', subject_list)

##############################################################################
#                           Preprocessing Nodes
##############################################################################
def get_first(_list):
    if isinstance(_list,list):
        return _list[0]
    else:
        return _list
#Parses subj info and returns the necessary dict for datasource.template_args
def getsubjinfo(subject_id,info,whatParadigm):
    #intially whatCond
    if whatParadigm == 'WM':
        outfield = 1
    elif whatParadigm == 'Nback_letters':
        outfield = 2
    elif whatParadigm == 'MCNAB' or whatParadigm == 'MCNAB2':
        outfield = 1
    return dict(
                func=[['subject_id',info[subject_id][outfield][0]]],
                struct=[['subject_id',info[subject_id][0][0]]],
                mag=[['subject_id',info[subject_id][3][0]]],
                phase=[['subject_id',info[subject_id][4][0]]]
                )

def weight_mean(image, art_file):
    """Calculates the weighted mean of a 4d image, where 
    
    the weight of outlier timpoints is = 0.
    
    Parameters
    ----------
    image : File to take mean
    art_file : text file specifying outlier timepoints
    
    Returns
    -------
    File : weighted mean image
    """
    import nibabel as nib
    import numpy as np
    from nipype.utils.filemanip import split_filename
    import os
    import nipype.interfaces.freesurfer as fs
    
    if not isinstance(image,list):
        image = [image]
    if not isinstance(art_file,list):
        art_file = [art_file]
    
    def try_import(fname):
        try:
            a = np.genfromtxt(fname)
            return np.atleast_1d(a).astype(int)
        except:
            return np.array([]).astype(int)
    
    mean_image_fname = os.path.abspath(split_filename(image[0])[1])
    
    total_weights = []
    meanimage = []
    for i, im in enumerate(image):
        img = nib.load(im)
        weights=np.ones(img.shape[3])
        weights[try_import(art_file[i])] = 1
        meanimage.append(img.shape[3]*np.average(img.get_data(), axis=3, weights=weights))
        total_weights.append(img.shape[3])
    mean_all = np.average(meanimage, weights=total_weights, axis=0)

    final_image = nib.Nifti1Image(mean_all, img.get_affine(), img.get_header()) 
    final_image.to_filename(mean_image_fname+'.nii.gz') 

    return mean_image_fname+'.nii.gz'


def trad_mot(subinfo,files):
    print "MADE IT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    # modified to work with only one regressor at a time...
    import numpy as np
    motion_params = []
    mot_par_names = ['Pitch (rad)','Roll (rad)','Yaw (rad)','Tx (mm)','Ty (mm)','Tz (mm)']
    if not isinstance(files,list):
        files = [files]
    if not isinstance(subinfo,list):
        subinfo = [subinfo]
    for j,i in enumerate(files):
        motion_params.append([[],[],[],[],[],[]])
        #k = map(lambda x: float(x), filter(lambda y: y!='',open(i,'r').read().replace('\n',' ').split(' ')))
        print i
        a = np.genfromtxt(i)
        for z in range(6):
            motion_params[j][z] = a[:,z].tolist()#k[z:len(k):6]
        
    for j,i in enumerate(subinfo):
        if i.regressor_names == None: i.regressor_names = []
        if i.regressors == None: i.regressors = []
        for j3, i3 in enumerate(motion_params[j]):
            i.regressor_names.append(mot_par_names[j3])
            i.regressors.append(i3)
    return subinfo

def noise_mot(subinfo,files,num_noise_components):
    noi_reg_names = map(lambda x: 'noise_comp_'+str(x+1),range(num_noise_components))
    noise_regressors = []
    if not isinstance(files,list):
        files = [files]
    for j,i in enumerate(files):
        noise_regressors.append([[]]*num_noise_components)
        k = map(lambda x: float(x), filter(lambda y: y!='',open(i,'r').read().replace('\n',' ').split(' ')))
        for z in range(num_noise_components):
            noise_regressors[j][z] = k[z:len(k):num_noise_components]
    for j,i in enumerate(subinfo):
        if i.regressor_names == None: i.regressor_names = []
        if i.regressors == None: i.regressors = []
        for j3,i3 in enumerate(noise_regressors[j]):
            i.regressor_names.append(noi_reg_names[j3])
            i.regressors.append(i3)
    return subinfo



datasource = pe.Node(interface=nio.DataGrabber(), name = 'datasource')
datasource.infields = ['subject_id']
datasource.outfields = ['func','struct','mag','phase']
datasource.inputs.base_directory = data_dir
datasource.inputs.template = '%s/niftis/%s.nii'

realign = pe.Node(interface=spm.Realign(), name="realign")
realign.inputs.register_to_mean = True

trad_motn = pe.Node(util.Function(input_names=['subinfo',
                                                   'files'],
                                      output_names=['subinfo'],
                                      function=trad_mot),
                        name='trad_motn')

ad = pe.Node(interface=ra.ArtifactDetect(), name="ad")
ad.inputs.zintensity_threshold = 2
ad.inputs.norm_threshold = 0.5
ad.inputs.use_differences = [True, False]
ad.inputs.mask_type = 'spm_global'
ad.inputs.parameter_source = 'SPM'

meanimg = pe.Node(util.Function(input_names=['image','art_file'],
                                       output_names=['mean_image'],
                                       function=weight_mean),
                                       name='meanimg')


"""
fieldmap = pe.Node(interface=fsl.utils.EPIDeWarp(), name='fieldmap_unwarp')
dewarper = pe.MapNode(interface=fsl.FUGUE(),iterfield=['in_file'],name='dewarper')


compcorr = create_compcorr()
compcorr.inputs.inputspec.num_components = 6
compcorr.inputs.inputspec.selector = [True,True]


noise_motn = pe.Node(util.Function(input_names=['subinfo',
                                                    'files',
                                                    'num_noise_components'],
                                       output_names=['subinfo'],
                                       function=noise_mot),
                         name='noise_motn')
noise_motn.inputs.num_noise_components = 6
"""
art = pe.Node(interface=ra.ArtifactDetect(), name="art")
art.inputs.use_differences = [True,False]
art.inputs.use_norm = True
art.inputs.norm_threshold = 1.0
art.inputs.zintensity_threshold = 3.0
art.inputs.mask_type = 'file'
art.inputs.parameter_source = 'SPM'

stimcor = pe.Node(interface=ra.StimulusCorrelation(), name="stimcor")
stimcor.inputs.concatenated_design = False

volsmooth = pe.Node(interface=spm.Smooth(), name="volsmooth")
volsmooth.inputs.fwhm = [6,6,6]

calcSurfReg = pe.Node(interface=fs.BBRegister(),name='calcSurfReg')
calcSurfReg.inputs.init = 'fsl'
calcSurfReg.inputs.contrast_type = 't2'
calcSurfReg.inputs.registered_file = True
calcSurfReg.inputs.out_fsl_file = True

makeImgNiiT = pe.MapNode(interface=fs.MRIConvert(),name='makeImgNiiT', iterfield=['in_file'])
makeImgNiiT.inputs.in_type = 'nifti1'
makeImgNiiT.inputs.out_type = 'niigz'


makeImgNiiCon = pe.MapNode(interface=fs.MRIConvert(),name='makeImgNiiCon', iterfield=['in_file'])
makeImgNiiCon.inputs.in_type = 'nifti1'
makeImgNiiCon.inputs.out_type = 'niigz'

calcVarMap = pe.MapNode(interface=fsl.maths.MultiImageMaths(),name='calcVarMap',iterfield=['in_file','operand_files'])
calcVarMap.inputs.op_string = '-div %s -sqr'
calcVarMap.inputs.output_type = 'NIFTI_GZ'


applySurfRegT = pe.MapNode(interface=fs.ApplyVolTransform(),name='applySurfRegT', iterfield = ['source_file'])

applySurfRegCon = pe.MapNode(interface=fs.ApplyVolTransform(),name='applySurfRegCon', iterfield = ['source_file'])

applySurfRegVar = pe.MapNode(interface=fs.ApplyVolTransform(),name='applySurfRegVar', iterfield = ['source_file'])

FreeSurferSource = pe.Node(interface=nio.FreeSurferSource(), name='FreeSurferSource')
FreeSurferSource.inputs.subjects_dir = subjects_dir

ApplyVolTransform = pe.Node(interface=fs.ApplyVolTransform(), name='applyreg')
ApplyVolTransform.inputs.inverse = True

Threshold = pe.Node(interface=fs.Binarize(),name='threshold')
Threshold.inputs.min = 1
Threshold.inputs.out_type = 'nii'

##############################################################################
#                           L1 Analysis Nodes
##############################################################################

# function to build the predicted timeseries for each subject:
#originally whatcond

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
    onsets_temp = os.path.join(onsets_dir, subject_id+'_*onsets_mcnab2.mat')
    onsets_files = sorted(glob.glob(onsets_temp))
    print onsets_files
    for r in range(len(onsets_files)):
        print "Run %d"%(r)
        mat = loadmat(onsets_files[r], struct_as_record=False)
        ons = mat['onsets'][0]
        nam = mat['names'][0]
        dur = mat['durations'][0]
        pmods = mat['pmod'][0]

        #Paradigm-specifics
        if whatParadigm == 'WMSTAT': 
            #24 types...
            names = ['inst_2r','inst_4r','inst_dry','inst_4ry','stim_2r','stim_4r','stim_dry','stim_4ry','probe_2r','probe_4r','probe_dry','probe_4ry']
            durations = []
            run_onsets = []
            for i in range(len(names)):
                print names[i]+": "
                durations.append([0])                
                run_onsets.append(ons[i][0])
        else:
            names = []
            durations = []
            run_onsets = []
            pmod = []
            for i in range(len(nam)):

                if not pmods[i].any():
                    pmod.append(None)
                else:
                    pmod.append(Bunch(name=["load"],param=[pmods[i][0].tolist()],poly=[1]))
                
                names.append(str(nam[i][0]))
                run_onsets.append(ons[i][0].tolist())
                durations.append(dur[i][0].tolist())


        output.insert(r,
            Bunch(conditions=deepcopy(names),
                onsets=deepcopy(run_onsets),
                durations=deepcopy(durations),
                amplitudes=None,
                tmod=None,
                pmod=None,
                regressor_names=regressor_names,
                regressors=regressors)) #here is where we can do linear, quad, etc detrending

    return output

# Model Specification (NiPype) Node
modelspec = pe.Node(interface=model.SpecifyModel(), name="modelspec")
modelspec.inputs.input_units = 'secs'
modelspec.inputs.time_repetition = 2.0
#128 in original script
modelspec.inputs.high_pass_filter_cutoff = 128 #np.inf #because of linear / quad regressors - otherwise ~160

# Level 1 Design (SPM) Node
level1design = pe.Node(interface=spm.Level1Design(), name= "level1design")
level1design.inputs.timing_units = 'secs' #modelspec.inputs.output_units
level1design.inputs.interscan_interval = 2.0
level1design.inputs.bases = {'hrf':{'derivs':[0,0]}}
level1design.inputs.model_serial_correlations = 'none' #'AR(1)'

# Level 1 Estimation node
level1estimate = pe.Node(interface=spm.EstimateModel(), name="level1estimate")
level1estimate.inputs.estimation_method = {'Classical' : 1}

# Constrast Estimation node

contrastestimate = pe.Node(interface = spm.EstimateContrast(), name="contrastestimate")
if paradigm == 'WMSTAT':

#INST = instruction and delay period before stim
    #INST: no contrasts
    c1 = ('inst_2r','T',['inst_2r'],[1])
    c2 = ('inst_4r','T',['inst_4r'],[1])
    c3 = ('inst_dry','T',['inst_dry'],[1])
    c4 = ('inst_4ry','T',['inst_4ry'],[1])
#preparing for stuff
    #INST: all; 
    c5 = ('inst: all','T',['inst_dry','inst_2r','inst_4r','inst_4ry'],[.25,.25,.25,.25])
#preparing for non-inhibition
    #INST: all - dist
    c6 = ('inst: all-dist','T',['inst_2r','inst_4r','inst_4ry'],[1./3,1./3,1./3])
#preparing for inhibition vs no inhibition
    #INST: distract > all; look for basal ganglia (caudate, putamen, etc)?
    c7 = ('inst: dist > nodist','T',['inst_dry', 'inst_2r', 'inst_4r', 'inst_4ry'],[1,-1./3,-1./3,-1./3])
    #INST: distract > nodist RY 
    c8 = ('inst: dry > 4ry','T',['inst_dry','inst_4ry'],[1,-1])
#preparing w/ LOAD; THIS IS NOT YET KNOWN TO PARTICIPANT, so EXPECT NO DIFFERENCE!!!
    #INST: Load (no distraction)
    c9 = ('inst: 4r > 2r','T',['inst_4r', 'inst_2r'],[1,-1])
    #INST: Load (4 things)
    c10 = ('inst: 4things','T',['inst_4r','inst_4ry'],[.5,.5])
    #INST: Load (2 things)
    c11 = ('inst: 2things','T',['inst_2r','inst_dry'],[.5,.5])
    #INST: Load (4 vs 2 things)
    c12 = ('inst: 4things > 2things','T',['inst_4r','inst_4ry','inst_2r'],[.5,.5,-1])
    #INST: Load (is distract the same as load of 2?)
    c13 = ('inst: dry > 2r','T',['inst_dry','inst_2r'],[1,-1])
#Holding shit in memory w/ COLOR @ const load; ALSO NOT KNOWN SO EXPECT NO DIFFERENCE!!!
    #STIM: Red and yellow > Red
    c14 = ('inst: 4ry > 4r','T',['inst_4ry','inst_4r'],[1,-1])


#STIM = stimulus presentation and delay before response
    #STIM: no contrasts
    c15 = ('stim_2r','T',['stim_2r'],[1])
    c16 = ('stim_4r','T',['stim_4r'],[1])
    c17 = ('stim_dry','T',['stim_dry'],[1])
    c18 = ('stim_4ry','T',['stim_4ry'],[1])
#Encoding/holding shit in memory
    #STIM: all; hope to see FRONTAL: ask amy about particular region! I forgot the name ??
    c19 = ('stim: all','T',['stim_dry','stim_2r','stim_4r','stim_4ry'],[.25,.25,.25,.25])
    #STIM: all - dist
    c20 = ('stim: all-dist','T',['stim_2r','stim_4r','stim_4ry'],[1./3,1./3,1./3])
#Holding shit in memory with inhibition vs without inhibition
    #STIM: distract > all; 
    c21 = ('stim: dist > nodist','T',['stim_dry', 'stim_2r', 'stim_4r', 'stim_4ry'],[1,-1./3,-1./3,-1./3])
    #STIM: distract > nodist RY; 
    c22 = ('stim: dry > 4ry','T',['stim_dry','stim_4ry'],[1,-1])
#Holding shit in memory w/ LOAD
    #STIM: Load (no distraction)
    c23 = ('stim: 4r > 2r','T',['stim_4r', 'stim_2r'],[1,-1])
    #STIM: Load (4 things)
    c24 = ('stim: 4things','T',['stim_4r','stim_4ry'],[.5,.5])
    #STIM: Load (2 things)
    c25 = ('stim: 2things','T',['stim_2r','stim_dry'],[.5,.5])
    #STIM: Load (4 vs 2 things)
    c26 = ('stim: 4things > 2things','T',['stim_4r','stim_4ry','stim_2r'],[.5,.5,-1])
    #STIM: Load (is distract the same as load of 2?)
    c27 = ('stim: dry > 2r','T',['stim_dry','stim_2r'],[1,-1])
#Holding shit in memory w/ COLOR @ const load
    #STIM: Red and yellow > Red
    c28 = ('stim: 4ry > 4r','T',['stim_4ry','stim_4r'],[1,-1])
    

#PROBE:
    #PROBE: no contrasts
    c29 = ('probe_2r','T',['probe_2r'],[1])
    c30 = ('probe_4r','T',['probe_4r'],[1])
    c31 = ('probe_dry','T',['probe_dry'],[1])
    c32 = ('probe_4ry','T',['probe_4ry'],[1])
#Responding to shit
    #PROBE: all; 
    c33 = ('probe: all','T',['probe_dry','probe_2r','probe_4r','probe_4ry'],[.25,.25,.25,.25])
    #PROBE: all - dist
    c34 = ('probe: all-dist','T',['probe_2r','probe_4r','probe_4ry'],[1./3,1./3,1./3])
#responding with inhibition vs responding without inhibition
    #PROBE: distract > all; 
    c35 = ('probe: dist > nodist','T',['probe_dry', 'probe_2r', 'probe_4r', 'probe_4ry'],[1,-1./3,-1./3,-1./3])
    #PROBE: distract > nodist RY; 
    c36 = ('probe: dry > 4ry','T',['probe_dry','probe_4ry'],[1,-1])
#responding w/ LOAD
    #PROBE: Load (no distraction)
    c37 = ('probe: 4r > 2r','T',['probe_4r', 'probe_2r'],[1,-1])
    #PROBE: Load (4 things)
    c38 = ('probe: 4things','T',['probe_4r','probe_4ry'],[.5,.5])
    #PROBE: Load (2 things)
    c39 = ('probe: 2things','T',['probe_2r','probe_dry'],[.5,.5])
    #PROBE: Load (4 > 2 things)
    c40 = ('probe: 4things > 2things','T',['probe_4r','probe_4ry','probe_2r'],[.5,.5,-1])    
    #PROBE: Load (is distract the same as load of 2?)
    c41 = ('probe: dry > 2r','T',['probe_dry','probe_2r'],[1,-1])
#responding w/ COLOR @ const load
    #PROBE: Red and yellow > Red
    c42 = ('probe: 4ry > 4r','T',['probe_4ry','probe_4r'],[1,-1])
    
    contrastestimate.inputs.contrasts = [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,c35,c36,c37,c38,c39,c40,c41,c42]
elif paradigm =='WM':
    #distract vs no distract
    c1 = ('inst_filt','T',['inst_DRY','inst_4RY','inst_2R'],[1,-.5,-.5])
    c2 = ('inst_weakfilt','T',['inst_DRY','inst_4RY'],[.5,-.5])
    
    c3 = ('stim_filt_stimcntrl','T',['stim_DRY','stim_4RY'],[.5,-.5])
    #load
    c4 = ('stim_load','T',['stim_4RY','stim_2R'],[.5,-.5])
    #dist/load?
    c5 = ('stim_filt_efficacy','T',['stim_DRY','stim_2R'],[.5,-.5])

    #just 4RY
    c6 = ('stim4RY','T',['stim_4RY'],[1])
    #just 2r
    c7 = ('stim2R','T',['stim_2R'],[1])
    #4ry + 2R
    c8 = ('stim4RY2R','T',['stim_4RY','stim_2R'],[.5,.5])

    contrastestimate.inputs.contrasts = [c1,c2,c3,c4,c5,c6,c7,c8]

elif paradigm =='Nback_spatial' or paradigm == 'Nback_letters':
    c1 = ('3_gt_1','T',['3BACK','1BACK'],[.5,-.5])
    c2 = ('3_gt_2','T',['3BACK','2BACK'],[.5,-.5])    
    c3 = ('3_gt_0','T',['3BACK','0BACK'],[.5,-.5])

    c4 = ('2_gt_1','T',['2BACK','1BACK'],[.5,-.5])    
    c5 = ('2_gt_0','T',['2BACK','0BACK'],[.5,-.5])

    c6 = ('1_gt_0','T',['1BACK','0BACK'],[.5,-.5])    
    c7 = ('0back','T',['0BACK'],[1])    
    c8 = ('1back','T',['1BACK'],[1])    
    c9 = ('2back','T',['2BACK'],[1])    
    c10 = ('3back','T',['3BACK'],[1])

    contrastestimate.inputs.contrasts = [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10]

elif paradigm =='MCNAB':
    #inst 2R 4RY DRY
    #stim 2R 4RY DRY
    #param_stim 2R 4RY DRY
    #delay 2R 4RY DRY
    #param_delay 2R 4RY DRY
    #probe 2R 4RY DRY

    c1 = ('inst_filt','T',['inst_dist','inst_nodist'],[1,-1])

    contrastestimate.inputs.contrasts = [c1]
elif paradigm =='MCNAB2':
    #inst 2R 4RY DRY
    #stim 2R 4RY DRY
    #param_stim 2R 4RY DRY
    #delay 2R 4RY DRY
    #param_delay 2R 4RY DRY
    #probe 2R 4RY DRY

    c1 = ('inst_filt','T',['inst_dist','inst_nodist'],[1,-1])
    c2 = ('load','T',['delay_4RY','delay_2R'],[1,-1])
    c3 = ('leaky-load','T',['delay_DRY','delay_2R'],[1,-1])
    contrastestimate.inputs.contrasts = [c1,c2,c3]

# Merge node to allow registration of both con and spmT images
mergenode = pe.Node(interface=util.Merge(2, axis='hstack'),name='merge')

#normalization workflow
from ants_norm_base import get_full_norm_workflow
normflow = get_full_norm_workflow('ANTS_norm')
normflow.inputs.inputspec.template_file = '/software/atlas/mmrr-21/MMRR21_template_to_MNI152.nii.gz'

# Datasink node for saving output of the pipeline
datasink = pe.Node(interface=nio.DataSink(), name="datasink")
datasink.inputs.base_directory = os.path.join(study_dir, 'l1_mcnab2')

#returns a large list of string tuples corresponding to these strings related to each contrast index
#to be used for substitutions node of dataSink
def getsubs(subject_id,contrast_list):
    subs = [('_subject_id_%s/'%subject_id,'')]
    for i in range(len(contrast_list),0,-1):
        subs.append(('_applySurfRegCon%d/'%(i-1),''))
        subs.append(('_applySurfRegT%d/'%(i-1),''))
        subs.append(('_applySurfRegVar%d/'%(i-1),''))
        subs.append(('_warp_images%d/'%(i-1),''))
        subs.append(('con_%04d_out_maths_warped'%(i),'var_%04d_out_warped'%(i)))
        subs.append(('spmT_%04d_out_warped'%(i),'spmT_%d_%s_out_warped'%(i,contrast_list[i-1][0])))
    return subs

# Declare and connect up the pipeline for preproc + level1 analysis
#Roughly understood by JSALVA

l1pipeline = pe.Workflow(name="l1_mcnab")
#l1pipeline.base_dir = os.path.abspath('/mindhive/gablab/lex/Domain-AAA/Analysis/pipeline/results/working/%s/'%(condition))
#JSALVA: change base_dir to one relevant to CARD
l1pipeline.base_dir = os.path.abspath('/mindhive/gablab/GATES/tmp_working/mcnab/%s/'%(paradigm))
l1pipeline.connect([
                    (infosource,datasource, [('subject_id','subject_id'),
                                                (('subject_id',getsubjinfo,info,paradigm),'template_args')]),
                    (datasource,realign,[('func','in_files')]),
                    (infosource,calcSurfReg,[('subject_id','subject_id')]),                                        
                    
                    (infosource,FreeSurferSource,[('subject_id','subject_id')]),
                    (FreeSurferSource,normflow,[(('aparc_aseg',get_first),'inputspec.segmentation')]),
                    (FreeSurferSource,normflow,[('orig','inputspec.brain')]),    
                    (FreeSurferSource,ApplyVolTransform,[('brainmask','target_file')]),

                    (calcSurfReg,ApplyVolTransform,[('out_reg_file','reg_file')]),
                    (ApplyVolTransform,Threshold,[('transformed_file','in_file')]),

                    (realign,art,[('realignment_parameters','realignment_parameters'),
                                    ('realigned_files','realigned_files')]),
                    (realign,ad,[('realignment_parameters','realignment_parameters'),
                                    ('realigned_files','realigned_files')]),
                    (ad,meanimg,[('outlier_files','art_file')]),
                    (realign,meanimg,[('realigned_files','image')]),

                    (meanimg,calcSurfReg,[('mean_image','source_file')]),
                    (meanimg,ApplyVolTransform,[('mean_image','source_file')]),

                    (realign,volsmooth,[('realigned_files','in_files')]),
                    (realign,trad_motn,[('realignment_parameters','files')]),
                    (Threshold,art, [('binary_file','mask_file')]),
                    (infosource,trad_motn,[(('subject_id',timeseries_design,paradigm,onsets_dir),'subinfo')]),
                    (trad_motn,modelspec,[('subinfo','subject_info')]),
                    (volsmooth,modelspec,[('smoothed_files','functional_runs')]),
                    (art,modelspec,[('outlier_files','outlier_files')]),
                    (modelspec,level1design,[('session_info','session_info')]),
                    (Threshold,level1design,[('binary_file','mask_image')]),
                    (level1design,level1estimate,[('spm_mat_file','spm_mat_file')]),
                    (level1estimate,contrastestimate,[('spm_mat_file','spm_mat_file'),
                                                        ('beta_images','beta_images'),
                                                        ('residual_image','residual_image')]),
                    (art,stimcor,[('intensity_files','intensity_values')]),
                    (realign,stimcor,[('realignment_parameters','realignment_parameters')]),
                    (level1design,stimcor,[('spm_mat_file','spm_mat_file')]),
                    (contrastestimate,makeImgNiiCon,[('con_images','in_file')]),
                    (makeImgNiiCon,normflow,[('out_file','inputspec.input_image')]),
                    (makeImgNiiCon,applySurfRegCon,[('out_file','source_file')]),
                    (calcSurfReg,applySurfRegCon,[('out_reg_file','reg_file'),
                                                ('registered_file','target_file')]),
                    (contrastestimate,makeImgNiiT,[('spmT_images','in_file')]),
                    (makeImgNiiT,applySurfRegT,[('out_file','source_file')]),
                    (calcSurfReg,applySurfRegT,[('out_reg_file','reg_file'),
                                                ('registered_file','target_file')]),
                    (calcSurfReg,normflow,[('out_fsl_file','inputspec.out_fsl_file')]),
                    (meanimg,normflow,[('mean_image','inputspec.mean_func')]),
                    (makeImgNiiCon,calcVarMap,[('out_file','in_file')]),
                    (makeImgNiiT,calcVarMap,[('out_file','operand_files')]),
                    (calcVarMap,applySurfRegVar,[('out_file','source_file')]),
                    (calcSurfReg,applySurfRegVar,[('out_reg_file','reg_file'),
                                                ('registered_file','target_file')]),
                    (normflow,datasink,[('outputspec.warped_image','norm_func.@warped_image')]),
                    (normflow,datasink,[('outputspec.warp_field','norm_anat.@warp_field')]),
                    (normflow,datasink,[('outputspec.affine_transformation','norm_anat.@affine')]),
                    (normflow,datasink,[('outputspec.inverse_warp','norm_anat.@inverse_warp')]),
                    (normflow,datasink,[('outputspec.unwarped_brain','norm_anat.@unwarped_brain')]),
                    (normflow,datasink,[('outputspec.warped_brain','norm_anat.@warped_brain')]),
                    (normflow,datasink,[('normalize_post_struct.fsl_reg_2_itk.fsl2antsAffine','norm_anat.@fslaffine')]),
                    (infosource,datasink,[('subject_id','container'),
                                            #JSALVA: substitutions takes list of 2-tuples that replaces default
                                            #path names; can be used to eliminate entire subdirectories created by default
                                            #by including '' as the second tuple
                                            (('subject_id',getsubs,contrastestimate.inputs.contrasts),'substitutions')]),
                    (FreeSurferSource,datasink,[('brain','subj_anat.@brain')]),
                    (meanimg,datasink,[('mean_image','subj_anat.@mean')]),
                    (realign,datasink,[('realignment_parameters','qc_realign.@realign_params')]),
                    (contrastestimate,datasink,[('con_images','subj_contrasts.@con'),
                                                ('spmT_images','subj_contrasts.@tmap')]),
                    (calcVarMap,datasink,[('out_file','subj_contrasts.@var')]),
                    (calcSurfReg,datasink,[('out_reg_file','surfreg'),
                                            ('out_fsl_file','surfreg.@fslmat'),
                                            ('registered_file','subj_anat.@reg_mean')]),
                    (applySurfRegCon,datasink,[('transformed_file','reg_cons.@con')]),
                    (applySurfRegT,datasink,[('transformed_file','reg_cons.@T')]),
                    (applySurfRegVar,datasink,[('transformed_file','reg_cons.@var')]),
                    (level1estimate,datasink,[('beta_images','model.@beta'),
                                                ('spm_mat_file','model.@spm'),
                                                ('mask_image','model.@mask'),
                                                ('residual_image','model.@res'),
                                                ('RPVimage','model.@rpv')]),
                    (art,datasink,[('outlier_files','qc_art.@outliers'),
                                    ('plot_files','qc_art.@motionplots'),
                                    ('statistic_files','qc_art.@statfiles'),
                                    ]),
                    (stimcor,datasink,[('stimcorr_files','qc_stimcor')]),
                    ])

##############################################################################
#                           L2 Analysis Nodes
##############################################################################

if paradigm == 'WM':
    n_contrasts = 8
elif paradigm == 'Nback_letters':
    n_contrasts = 6
elif paradigm == 'MCNAB':
    n_contrasts = 1
elif paradigm == 'MCNAB2':
    n_contrasts = 3

def selectsubs(file_list,subject_list):
    returnlist = [] 
    for _file in file_list:
        for subj in subject_list:
            if subj in _file:
                returnlist.append(_file)
    return returnlist

contrast_list = range(1,n_contrasts+1)
func="""
def gunzip(compressed_file):
    import gzip
    import string
    import os
    r_file = gzip.GzipFile(compressed_file,'r')
    decompressed_file = string.rstrip(compressed_file,'.gz')
    w_file = open(decompressed_file,'w')
    w_file.write(r_file.read())
    w_file.close()
    r_file.close()
    os.unlink(compressed_file)
    return decompressed_file
"""
# Infosource Node (get subject specific run information)
l2infosource = pe.Node(interface=util.IdentityInterface(fields=['contrasts']), name="l2infosource")
l2infosource.iterables = ('contrasts', contrast_list)

l2datasource = pe.Node(nio.DataGrabber(infields=['contrast']), name='l2datasource')
l2datasource.inputs.base_directory = '/mindhive/gablab/GATES/Analysis/%s/'%(paradigm)
l2datasource.inputs.template = 'l1_mcnab2/*/norm_func/con_%04d*'

gunzip = pe.MapNode(util.Function(input_names=['compressed_file'],output_names=['decompressed_file']),name='gunzip', iterfield=['compressed_file'])
gunzip.inputs.function_str = func

# setup a 1-sample t-test node
onesamplettestdes = pe.Node(interface=spm.OneSampleTTestDesign(), name="onesampttestdes")
onesamplettestdes.inputs.explicit_mask_file = os.path.abspath('/software/spm8/apriori/brainmask.nii')

l2estimate = pe.Node(interface=spm.EstimateModel(), name="onesampleestimate")
l2estimate.inputs.estimation_method = {'Classical' : 1}

l2conestimate = pe.Node(interface = spm.EstimateContrast(), name="l2conestimate")
# 'mean' is the name SPM gives to this value in the design matrix
# format is as follows:  [name of contrast, type of test, name given in the SPM design matrix, weight]
con1 = ('Group(pos)','T', ['mean'],[1])
con2 = ('Group(neg)','T', ['mean'],[-1])

l2conestimate.inputs.contrasts = [con1,con2]
l2conestimate.inputs.group_contrast = True

l2datasink = pe.Node(interface=nio.DataSink(), name="l2_mcnab_model2_adults")
l2datasink.inputs.base_directory = os.path.join(study_dir, 'l2_mcnab_model2_adults/')

l2pipeline = pe.Workflow(name='l2_mcnab_storage_adults')
l2pipeline.base_dir = os.path.abspath('/mindhive/gablab/GATES/tmp_working/mcnabl2/%s/'%(paradigm))
l2pipeline.connect([(l2infosource,l2datasource,[('contrasts','contrast')]),
                    (l2datasource,onesamplettestdes,[(('outfiles',selectsubs, subject_list),'in_files')]),#gunzip,[('outfiles','compressed_file')]),
#                    (gunzip,
                    (onesamplettestdes,l2estimate,[('spm_mat_file','spm_mat_file')]),
                    (l2estimate,l2conestimate,[('spm_mat_file','spm_mat_file'),
                                                ('beta_images','beta_images'),
                                                ('residual_image','residual_image')]),
                    (l2conestimate,l2datasink,[('con_images','group_contrasts.@con'),
                                                ('spmT_images','group_contrasts.@tmap')]),

                ])

##############################################################################
#                           Run The Damn Show
##############################################################################
def save_design_image():
    import os
    import nipype.interfaces.matlab as matlab
    import glob
    disp_model = matlab.MatlabCommand()
    script='addpath \'/mindhive/gablab/users/jsalva/scripts/GATES/\';'
    model_dirs = glob.glob(os.path.join(study_dir,'l1output_fldmap','*','model'))
    for model_dir in model_dirs:
        script=script+'designmatrix(\'%s\');'%(model_dir)
    disp_model.inputs.script=script
    disp_model.run()
    return




"""
Make the magic happen!
"""
qsubargs = '-V -q max30 -l host=ba5'

if levelToRun in ['l1','L1']:
    if numCores == "PBS":
        l1pipeline.run(plugin='PBS',plugin_args={'qsub_args':qsubargs})
#        l1pipeline.write_graph();
        
    elif numCores <= 1:
        l1pipeline.run()        
#        l1pipeline.write_graph();

    else:
        l1pipeline.run(plugin='MultiProc',plugin_args={'n_procs':numCores})
#        l1pipeline.write_graph();
#    save_design_image()

elif levelToRun in ['l2','L2']:
    if numCores == "PBS":
        l2pipeline.run(plugin='PBS',plugin_args={'qsub_args':qsubargs})
    elif numCores <= 1:
        l2pipeline.run()
    else:
        l2pipeline.run(plugin='MultiProc',plugin_args={'n_procs':numCores})

