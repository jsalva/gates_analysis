import nipype.interfaces.ants as ants
import nipype.pipeline.engine as pe
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as util 
import numpy as np
import os
import nipype.interfaces.io as nio   

getfirst = lambda x: x[0]

def create_inverse_warpflow(name='inverse_warpflow'):
    import nipype.interfaces.ants as ants
    import nipype.pipeline.engine as pe
    import nipype.interfaces.freesurfer as fs
    import nipype.interfaces.utility as util 

    getfirst = lambda x: x[0]
    #inverse warp workflow to take normalized files to native space
    inverse_warpflow = pe.Workflow(name=name)

    inputspec = pe.Node(util.IdentityInterface(
        fields=['images_to_warp','reference_anatomical','inverse_warp','affine']),
        name='inputspec')

    #convert the SPM img files to niftis
    SPMimg2nii = pe.MapNode(fs.preprocess.MRIConvert(),
            name='SPMimg2nii',iterfield=['in_file'])
    SPMimg2nii.inputs.out_type = 'nii'
    inverse_warpflow.connect([(inputspec,SPMimg2nii,
        [('images_to_warp','in_file')])])

    #convert Freesurfer mgz files to niftis
    FSmgz2nii = pe.MapNode(fs.preprocess.MRIConvert(),
            name='FSmgz2nii',iterfield=['in_file'])
    FSmgz2nii.inputs.out_type = 'nii'
    
    inverse_warpflow.connect([(inputspec,FSmgz2nii,
        [('reference_anatomical','in_file')])])

    #put the transforms into a list in the proper order to go from normalized to native space
    collect_transforms = pe.Node(util.Merge(2),
        name='collect_transforms')
    inverse_warpflow.connect([(inputspec,collect_transforms,
        [('affine','in1'),
        ('inverse_warp','in2')])])

    #collect the moving, reference, and transformation series to warp via ANTS WIMT
    norm2native = pe.MapNode(ants.WarpImageMultiTransform(),
        name='norm2native',iterfield=['moving_image'])
    norm2native.inputs.invert_affine = [1]
    norm2native.inputs.use_nearest = True
    norm2native.inputs.reslice_by_header = True
    inverse_warpflow.connect([(SPMimg2nii,norm2native,
        [('out_file','moving_image')])])
    inverse_warpflow.connect([(FSmgz2nii,norm2native,
        [(('out_file',getfirst),'reference_image')])])
    inverse_warpflow.connect([(collect_transforms,norm2native,
        [('out','transformation_series')])])

    outputspec = pe.Node(util.IdentityInterface(
        fields=['warped_images']),
        name='outputspec')
    inverse_warpflow.connect([(norm2native,outputspec,
        [('output_image','warped_images')])])
    return inverse_warpflow

def create_segstats_workflow(name='segstats_workflow',regions = None):
    import nipype.pipeline.engine as pe
    import nipype.interfaces.freesurfer as fs
    import nipype.interfaces.utility as util 
    
    segstats_workflow = pe.Workflow(name=name)
    
    inputspec = pe.Node(util.IdentityInterface(
        fields=['segmentation_file','images_to_extract_stats_from','subjects_dir','subject']),
        name='inputspec')
    
    segstats = pe.MapNode(fs.SegStats(),name="segstats", iterfield=['in_file'])
    segstats.inputs.color_table_file = '/software/Freesurfer/5.1.0/FreeSurferColorLUT.txt'

    if not regions:
        print "no regions specified; doing all..."
    else:
        segstats.inputs.segment_id = regions

    segstats_workflow.connect([(inputspec,segstats,
        [('segmentation_file','segmentation_file'),
        ('images_to_extract_stats_from','in_file'),
        ('subjects_dir','subjects_dir')])])

    outputspec = pe.Node(util.IdentityInterface(
        fields=['summary_files']),
        name='outputspec')
    segstats_workflow.connect([(segstats,outputspec,
        [('summary_file','summary_files')])])

    return segstats_workflow


#project back into native from normalized space

#re-label ROI mask by identity determined from majority rule in segmentation file
def label_cluster_by_identity(binary_mask,segmentation_file):
    import os
    import nibabel as nb
    import numpy as np
    from scipy import stats
    mask = nb.load(binary_mask)
    mask_data = mask.get_data()
    header = mask.get_header()
    affine = mask.get_affine()

    
    labeled_brain = nb.load(segmentation_file)
    labeled_data = labeled_brain.get_data()

    indeces = np.where((mask_data == 1)&(labeled_data!=0))


    majority_label = int(stats.mode(labeled_data[indeces],axis=None)[0][0])
    path = '%s/%s_%s'%(os.path.split(binary_mask)[0],str(majority_label),os.path.split(binary_mask)[1])
    
    new_mask_data = np.zeros(mask.shape[0:3])
    new_mask_data[indeces] = majority_label
    new_roi_mask = nb.Nifti1Image(new_mask_data,affine,header)
    new_roi_mask.to_filename(path) 

    return path


#combine all stats into 1 table >> append things to the bottom of the table rather than overwrite it; make a different table for each measure
def asegstats2table(stats_list,measure,subject_list):
    import os
    import numpy as np
    sorted_list = np.array(sorted(map(os.path.abspath,stats_list)))
    firstrow = ['subject']
    secondrow = ['segmentation']

    cmd = "asegstats2table --all-segs -m %s -t allsubs.%s.table --transpose "%(measure,measure)

    for subject in subject_list:
        tmp_list = sorted_list[np.array(map(lambda x: x.find(subject),sorted_list))!=-1]
        
        for statsfile in tmp_list:
            firstrow.append(subject)
            secondrow.append(os.path.split(statsfile)[-1])
            cmd = cmd+' -i '+str(statsfile)

    print cmd

    os.system(cmd)
   
    readfilehandle = open(os.path.abspath('allsubs.%s.table'%(measure)),'r')
    writefilehandle = open(os.path.abspath('allsubs_formatted.%s.table'%(measure)),'w')
    content_storage = []
    for num, line in enumerate(readfilehandle):
        content_storage.insert(num,line.split())
    readfilehandle.close()
    content_storage.insert(0,tuple(secondrow))
    content_storage.insert(0,tuple(firstrow))
    print content_storage


    content_transpose = zip(*content_storage)

    formatstr = '%s\t'
    for counter in range(len(content_transpose[0])-1):
        formatstr=formatstr+'%s\t'

    writefilehandle.write('\n'.join(formatstr%line for line in content_transpose))
    writefilehandle.close()
    return os.path.abspath('allsubs_formatted.%s.table'%(measure))


#####################################################
#Project specific shit
#####################################################
def getsubs(subject_id):
    subs = [('_subject_id_%s/'%subject_id,'')]
    subs.append(('/_segstats',''))
    subs.append(('/summary',''))
    return subs
    

#where your data live, yo
subjects_dir= '/mindhive/gablab/GATES/data'
#brains you want to extract stats from
coregistered_functional_template = 'Analysis/WM/l1sink_*/%s/reg_cons/con*.nii*'
#Freesurfer aparc+aseg
segmentation_file_template = 'data/%s/mri/aparc+aseg.mgz'
#location of ROIs
roi_file_template = 'Analysis/ROI/Norm_ROIs_101112/*.img'
#ANTS stuff
reference_template = 'Analysis/WM/l1sink_*/%s/norm_anat/orig_out_masked*'
inverse_warp_template = 'Analysis/WM/l1sink_*/%s/norm_anat/ants_InverseWarp.nii*'
affine_template = 'Analysis/WM/l1sink_*/%s/norm_anat/ants_Affine.txt'

datasink_dir = '/mindhive/gablab/GATES/Analysis/roi_analysis_poster'.
from subject_info import *
"""
subject_list=['300','301','302','303','304','305','306','307','308','309','310','311','312',
'313','314','315','316','317','318','319','320','321','322','323','325','326','327',
'328','329','330','332','333','334','335','336','337','338','339','340','341','342',
'343','344','345','346','347','348','349','350','351','352','353','354','355','356',
'357','358','359','360','361','362','363','364','365','366','367','368','369','370',
'401','402','403','404','405','406','407','408','409','410','411','412',
'413','414',
'500','501','502','503','504','505','b500_2','b501_2','b502_2','b503_2']"""
regions = [11,12,13,17,1003,1027,1028,1029,50,51,52,53,2003,2027,2028,2029,4,43,31,63]



ROI = pe.Workflow(name='ROI')
ROI.base_dir = '/mindhive/scratch/jsalva/ROI/poster/'

infosource = pe.Node(util.IdentityInterface(fields=['subject_id']),
                         name='infosource')
infosource.iterables = ('subject_id', subject_list)

datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                         outfields=['brain','seg_file','func_files','inverse_warp_file', 'affine_file']),
                         name = 'datasource')
datasource.inputs.base_directory = '/mindhive/gablab/GATES'
datasource.inputs.template ='*'
datasource.inputs.field_template = dict(brain=reference_template, 
    seg_file=segmentation_file_template, 
    func_files =coregistered_functional_template, 
    inverse_warp_file=inverse_warp_template, 
    affine_file=affine_template)
datasource.inputs.template_args = dict(brain=[['subject_id']], 
    seg_file=[['subject_id']],
    func_files=[['subject_id']],
    inverse_warp_file=[['subject_id']],
    affine_file=[['subject_id']])
ROI.connect([(infosource,datasource,
    [('subject_id','subject_id')])])
"""
ROIsource = pe.Node(interface=nio.DataGrabber(outfields=['func_roi_files']),
                         name = 'ROIsource')
ROIsource.inputs.base_directory = '/mindhive/gablab/GATES'
ROIsource.inputs.template = roi_file_template


unwarp_rois = create_inverse_warpflow('unwarp_rois')
ROI.connect([(ROIsource,unwarp_rois,
    [('func_roi_files','inputspec.images_to_warp')])])
ROI.connect([(datasource,unwarp_rois,
    [('brain','inputspec.reference_anatomical'),
    ('inverse_warp_file','inputspec.inverse_warp'),
    ('affine_file','inputspec.affine')])])

label_rois = pe.MapNode(util.Function(
        input_names=['binary_mask','segmentation_file'],
        output_names=['labeled_rois'],
        function=label_cluster_by_identity),
        name='label_rois', iterfield=['binary_mask'])
ROI.connect([(unwarp_rois,label_rois,
    [('outputspec.warped_images','binary_mask')])])
ROI.connect([(datasource,label_rois,
    [('seg_file','segmentation_file')])])

functional_segstats = create_segstats_workflow('fucntional_segstats')
f_segstats_inputspec = functional_segstats.get_node('inputspec')
f_segstats_inputspec.inputs.subjects_dir = subjects_dir
ROI.connect([(datasource,f_segstats_inputspec,
    [('func_files','images_to_extract_stats_from')])])
ROI.connect([(infosource,f_segstats_inputspec,
    [('subject_id','subject')])])
"""

structural_segstats = create_segstats_workflow('structural_segstats',regions)
s_segstats_inputspec = structural_segstats.get_node('inputspec')
s_segstats_inputspec.inputs.subjects_dir = subjects_dir
ROI.connect([(datasource,s_segstats_inputspec,
    [('seg_file','segmentation_file'),
    ('func_files','images_to_extract_stats_from')])])
ROI.connect([(infosource,s_segstats_inputspec,
    [('subject_id','subject')])])

"""
collect_stat_files = pe.MapNode(util.IdentityInterface(
    fields=['structural_summary_files','functional_summary_files','subject_id']),   
    name='collect_stat_files', iterfield=['structural_summary_files','functional_summary_files','subject_id'])
ROI.connect([(structural_segstats,collect_stat_files,
        [('outputspec.summary_files','structural_summary_files')])])
ROI.connect([(functional_segstats,collect_stat_files,
        [('outputspec.summary_files','functional_summary_files')])])
ROI.connect([(infosource,collect_stat_files,
        [('subject_id','subject_id')])])
"""

datasink = pe.Node(interface=nio.DataSink(),name='datasink')
datasink.inputs.base_directory = datasink_dir
ROI.connect([(infosource,datasink,[('subject_id','container'),
    (('subject_id',getsubs),'substitutions')])])
ROI.connect([(structural_segstats,datasink,#collect_stat_files,datasink,
    [#('functional_summary_files','functional'),
    #('structural_summary_files','structural')])])
    ('outputspec.summary_files','structural')])])

ROI.run(plugin='MultiProc',plugin_args={'n_procs':32})



STATS = pe.Workflow(name='STATS')
STATS.base_dir = '/mindhive/scratch/jsalva/ROI/poster/STATS'

statsgrabber = pe.Node(interface=nio.DataGrabber(outfields=['stats']),
                         name = 'statsgrabber')
statsgrabber.inputs.base_directory = datasink_dir
statsgrabber.inputs.template ='*/structural*.stats'


combine_stat_files = pe.MapNode(util.Function(
   input_names=['stats_list','measure','subject_list'],
   output_names=['combined_summary_file'],
   function=asegstats2table),
   name='combine_stat_files', iterfield=['measure'])
combine_stat_files.inputs.measure = ['mean','std','volume']
combine_stat_files.inputs.subject_list = subject_list
STATS.connect([(statsgrabber,combine_stat_files,
    [('stats','stats_list')])])

subs = lambda x: [('_combine_stat_files./','')]
statsink = pe.Node(interface=nio.DataSink(),name='statsink')
statsink.inputs.base_directory = datasink_dir
STATS.connect([(statsgrabber,statsink,
    [(('stats',subs),'regexp_substitutions')])])
STATS.connect([(combine_stat_files,statsink,
    [('combined_summary_file','summary_file')])])


STATS.run(plugin='MultiProc',plugin_args={'n_procs':32})






