#!/usr/bin/env python

"""
ATTENTION: This script assumes that when running it, you are in the directory where all subject directories
are located or should be created (when '-f' or '--fetch' is specified). Also if you've already fetched the data from 
sigma, make sure that the dicoms are in SUBJECT_ID/dicoms! :)
"""

from optparse import OptionParser
import os

def main(options, args):
    # check if in the right directory
    workingdir = checkdir(os.getcwd())
    
    ### SPECIFY SCANS: you can do this right here in the script or later when asked
    #########
    # Userinput required: define scans, nr of volumes and name under which they should be saved 
    # when more than one is found, 2.nii etc will automatically be appended.
    # assumes that scans starting with 'T1_MPRAGE' are the anatomical scans and
    # will be saved as .nii and .mgz
    # --> make sure this is correct (for your own sake...) 
 #   ['T1_MPRAGE_2530_256x176_4e_p2_1mm_iso',1,'mprag'],
   
    scans = [ 

    ['T1_MPRAGE_2530_220x176_4e_p2_1mm_iso',1,'mprag'],
    ['T1_MPRAGE_1mm_iso',1,'mprag_wrongprotocol'],
    ['T1_MPRAGE_2530_256x176_4e_p2_1mm_iso',1,'mprag_wrongprotocol'],
    ['ge_functional_WM_imp',267,'Nback'],
    ['ge_functional_WM_imp_2meas',276,'Nback_276tr'],
    ['ge_functional_WM_imp',199,'Nback_199tr'],
    ['ep2d_func_pace_NBACK',267,'Nback'],
    ['ep2d_func_pace_Vigilance',190,'vigilance'],
    ['ep2d_func_pace_HRF',75,'hrf'],
    ['ge_functional_WM_imp',203,'WM'],
    ['ge_functional_WM_imp',204,'WM'],
    ['ep2d_func_pace_WM_filt',204,'WM'],
   
   
    ['field_mapping_WM',1,'fieldmap_func'],
    ['field_mapping_2x2x2',1,'fieldmap_resting'],
 

    
    ['ep2d_Resting_2x2x2_32Ch',62,'Resting'],
    ['DIFFUSION_HighRes',70,'DTI']]
    

    #########
    scans_correct = False
    while (scans == [] or not scans_correct):
        print "you have currently specified the following scans for conversion:"
        print scans
        print "if this is NOT CORRECT, you can either terminate the script right now (<ctrl-c>) and ";
        print "specify the scans in the script itself and run it again OR continue and do it interactively"
        #correct = raw_input("Are the scans correct? (yes/no): ")
        correct = 'yes'
        if correct == 'yes' or correct == 'y':
            scans_correct = True
        else:
            scans = specify_scans()
    
    print "\nOK, WE'RE ALL SET!!\n"
    
    ### go through all subjects given as arguments when the script was called              
    for subjectid in args:
        print "\n PROCESSING SUBJECT: %s" %subjectid
        
        subject_dir = os.path.abspath(subjectid)
        dicom_dir = os.path.abspath(os.path.join(subjectid,'dicoms'))
        ### FETCH DATA FROM SIGMA
        if options.sigma_fetch:
            if not os.path.isdir(dicom_dir):
                os.makedirs(dicom_dir)
            # Call external python script to copy over dicom files from sigma
            fetch_cmd = "fetch_dicoms -s %s -l -q -d %s"%(subjectid, dicom_dir)
            print "[==        ]    ... fetching dicoms from sigma ..."
            #print fetch_cmd # use this for debugging instead of running the actuall command
            os.system(fetch_cmd)
        
        ### CONVERT DICOMS TO NIFTIES
        ### dicoms must be in SUBJECT_ID/dicoms!
        if not os.path.isdir(dicom_dir):
            print "ERROR: Skipping subject %s because the subject's dicom directory is not found!" %subjectid
            print "-->    make sure you're currently in the right directory (has to contain the subject directories)"
            print "-->    make sure the directory containing the subject's dicoms is labeled 'dicoms'"
            import sys
            sys.exit()
            continue
            
        ### 1. create scanonly file calling external command
        scanonlyfile = '%s/run.info' %(subject_dir)
        if not os.path.isfile(scanonlyfile):
            info_cmd = 'unpacksdcmdir -src %s -targ %s -scanonly %s' %(dicom_dir, subject_dir, scanonlyfile)
            print "[====      ]     ... creating %s in subject dir ..." %scanonlyfile
            #print info_cmd # use this for debugging instead of running the actuall command
            os.system(info_cmd)
        
        ### 2. from the file, extract further details for every scan
        scanonly_output = open(os.path.abspath(os.path.join(subject_dir,scanonlyfile)), 'r')
        # dictionary: dicoms_info['dicom_name'] = ('scan_name',TRs)
        # --> extracts this from each line of the scanoutput file 
        dicoms_info = dict()
        for line in scanonly_output.readlines():
            scan_details = line.split()
            dicoms_info[scan_details[-1]] = (scan_details[1],int(scan_details[-2]))
            
        ### 3. process all specified scans
        print "[=====     ]     ... processing specified scans ..."
        for scan in scans:
            dicoms = []
            # for each dicom look if it match the specified scans (i.e. name & # of TRs), only take nomoco scans
            dicom_list = dicoms_info.keys()
            sorted_dicoms = dict()
            for dicom in dicom_list:
                print dicom
                sorted_dicoms[int(dicom.split('-')[1])] = dicom
            for val,dicom in sorted(sorted_dicoms.iteritems()):
                info = dicoms_info[dicom]
                if info == (scan[0],scan[1]) and not is_moco(os.path.join(dicom_dir,dicom)):
                    dicoms.append(dicom)
            print "[======    ] The following nomoco dicoms match scan %s and will be converted now:" %scan[0]
            print dicoms
            # create niftis and mri directory
            niftis_dir = os.path.abspath(os.path.join(subjectid,'niftis'))
            mri_dir = os.path.abspath(os.path.join(subjectid,'mri'))
            fieldmap_dir = os.path.abspath(os.path.join(subjectid,'fieldmaps'))
            if not os.path.isdir(niftis_dir):
                os.makedirs(niftis_dir)
            if not os.path.isdir(mri_dir):
                os.makedirs(mri_dir)
            if not os.path.isdir(fieldmap_dir):
                os.makedirs(fieldmap_dir)
                
            ## convert each dicom
            for i, dicom in enumerate(dicoms, start=1):
                input_dicom = os.path.abspath(os.path.join(dicom_dir, dicom))
                # .nii files get saved as specified in scans, however there gets always a number added
                output_nifti = os.path.abspath(os.path.join(niftis_dir,scan[2]+str(i)+'.nii'))
                # call external command (freesurfer) for conversion
                convert_cmd = 'mri_convert %s %s' %(input_dicom, output_nifti)
                #print convert_cmd # use this for debugging instead of running the actuall command
                os.system(convert_cmd)
                print "--> created %s" %output_nifti
                
                # if T1_MPRAGE scan: also do anatomical conversion
                # --> usually this is the only one (only needed for recon-all), if applicable, also for flash scans
                if scan[0].startswith("T1_MPRAGE"):
                    output_mgz = os.path.join(mri_dir, '001.mgz') 
                    mriconvert_cmd = 'mri_convert %s %s' %(input_dicom, output_mgz)
                    #print mriconvert_cmd # use this for debugging instead of running the actuall command
                    os.system(mriconvert_cmd)
                    print "--> created %s" %output_mgz

                if scan[0].startswith("field_mapping"):
                    if i == 1:
                        fieldmap_nifti = os.path.abspath(os.path.join(fieldmap_dir,scan[2]+'_mag.nii'))
                    else:
                        fieldmap_nifti = os.path.abspath(os.path.join(fieldmap_dir,scan[2]+'_phase.nii'))
                    # call external command (freesurfer) for conversion
                    fieldmap_convert_cmd = 'mri_convert %s %s' %(input_dicom, fieldmap_nifti)
                    #print mriconvert_cmd # use this for debugging instead of running the actuall command
                    os.system(fieldmap_convert_cmd)
                    
        print "[==========]  Congraz, all specified scans for the subject were successfully converted!"
    
    print "\n\n DONE WITH ALL GIVEN SUBJECTS. Have a great day!" 
        
            
def checkdir(directory):
    """ make sure you are in the right directory where either all subject folders
    are located or will be created when --fetch is specified 
    returns the verified directory"""
    while True:
        print "Just double checking: Is this the directory where all subject directories are located or should be created? :"
        print directory
        print "(Even if you're just running this with only 1 subject, make sure you are not IN the subject's directory)"
        #check = raw_input("(yes/no): ")
        check = 'no'
        if check == 'yes' or check == 'y':
            return directory
        else:
            worked = False
            while not worked:
                #changedir = raw_input("please enter absolut or relative path to the right directory:\n")
                changedir = '/gablab/p/GATES/data/'
                try:
                    os.chdir(changedir)
                    directory = os.getcwd()
                    return directory
                    worked = True
                except Exception, e:
                    print "this does not seem to be a valid directory...try again =)"
                    #print e # for debugging

def specify_scans():
    """ Interactivly enter the scans that should be processed for the subjects """
    scans = []
    done = False
    i = 0
    while not done:
        scans.append([raw_input('Please enter scan name: ')])
        scans[i].append(raw_input('Please enter expected number of TRs: '))
        scans[i].append(raw_input('Please enter name under which the scan(s) should be saved (without adding a run number or .nii): '))             
        print "you've currently specified %i scans" %len(scans)
        more = raw_input('do you want to add another scan? (yes/no): ')
        if not (more == 'yes' or more == 'y'):
            done = True
        i += 1
    return scans

def is_moco(dcmfile):
    """Determine if a run has on-line motion correction. Thanks Todd :)"""
    import sys
    import subprocess
    cmd = ['mri_probedicom', '--i', dcmfile, '--t', '8', '103e']
    proc  = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    stdout = proc.communicate()[0]
    return stdout.strip().startswith('MoCoSeries')

if __name__ == "__main__":
    # parse arguments
    parser = OptionParser(usage="usage: %prog [OPTIONS] SUBJECT_ID1 [SUBJECT_ID2 ...]")
    parser.add_option("-f", "--fetch", action="store_true", dest="sigma_fetch", default=False,
                      help="fetch dicoms from sigma")

    (options, args) = parser.parse_args()
    # must be called with at least 1 subject id, if so, call main()
    if len(args) < 1:
        parser.print_help()
    else:
        main(options, args)
