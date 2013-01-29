from bips.workflows.base import get_config

uuid = '7757e3168af611e1b9d5001e4fb1404c' 

c = get_config(uuid)

"""surf_dir : Freesurfer subjects directory"""

c.surf_dir = '/mindhive/xnat/surfaces/GATES'

"""compcor_select : The first value in the list corresponds to applying                                        t-compcor, and the second value to a-compcor. Note:                                        both can be true"""

c.compcor_select = (True, True)

"""crash_dir : Location to store crash files"""

c.crash_dir = '/mindhive/scratch/jsalva/preproc_betaseries_bips'

"""fwhm : Full width at half max. The data will be smoothed at all values                              specified in this list."""

c.fwhm = [0.0, 5.0]

"""sink_dir : Location where the BIP will store the results"""

c.sink_dir = '/mindhive/gablab/GATES/Analysis/preproc_betaseries_bips'

"""reg_params : None"""

c.reg_params = (True, True, True, True, True)

"""use_fieldmap : True to include fieldmap distortion correction. Note: field_dir                                      must be specified"""

c.use_fieldmap = False

"""z_thresh : z thresh for art"""

c.z_thresh = 3.0

"""magnitude_template : None"""

c.magnitude_template = '%s/magnitude.nii.gz'

"""TR : TR of functional"""

c.TR = 2.0

"""field_dir : Base directory of field-map data (Should be subject-independent)                                                  Set this value to None if you don't want fieldmap distortion correction"""

c.field_dir = ''

"""working_dir : Location of the Nipype working directory"""

c.working_dir = '/mindhive/scratch/jsalva/preproc_betaseries_bips'

"""subjects : Subject id's. Note: These MUST match the subject id's in the                                 Freesurfer directory. For simplicity, the subject id's should                                 also match with the location of individual functional files."""

c.subjects = ['300', '301', '302', '303', '304', '305', '306', '307', '308', '309', '310', '311', '312', '313', '314', '315']

"""SliceOrder : None"""

c.SliceOrder = ''

"""hpcutoff : highpass cutoff"""

c.hpcutoff = 128.0

"""plugin_args : Plugin arguments."""

c.plugin_args = {'qsub_args': '-q many'}

"""Interleaved : True for Interleaved"""

c.Interleaved = False

"""test_mode : Affects whether where and if the workflow keeps its                             intermediary files. True to keep intermediary files. """

c.test_mode = False

"""num_noise_components : number of principle components of the noise to use"""

c.num_noise_components = 6

"""phase_template : None"""

c.phase_template = '%s/phase.nii.gz'

"""TE_diff : difference in B0 field map TEs"""

c.TE_diff = 0.0

"""do_slicetiming : Perform slice timing correction"""

c.do_slicetiming = False

"""run_using_plugin : True to run pipeline with plugin, False to run serially"""

c.run_using_plugin = True

"""echospacing : EPI echo spacing"""

c.echospacing = 0.0

"""json_sink : Location to store json_files"""

c.json_sink = '/mindhive/gablab/GATES/Analysis/preproc_betaseries_bips'

"""norm_thresh : norm thresh for art"""

c.norm_thresh = 1.0

"""plugin : plugin to use, if run_using_plugin=True"""

c.plugin = 'PBS'

"""use_advanced_options : None"""

c.use_advanced_options = False

"""lowpass_freq : None"""

c.lowpass_freq = 1.0

"""base_dir : Base directory of data. (Should be subject-independent)"""

c.base_dir = '/mindhive/gablab/GATES/data'

"""sigma : 2D spatial gaussing smoothing stdev (default = 2mm)"""

c.sigma = 2

"""highpass_freq : None"""

c.highpass_freq = 0.0079

"""func_template : None"""

c.func_template = '%s/niftis/WM1.nii'

"""advanced_script : None"""

c.advanced_script = """  """ 

