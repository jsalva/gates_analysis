import numpy as np
import os
from glob import glob
from matplotlib.pyplot import *
import sys
#313, 317 outlier
kids = ["300", "301", "302", "303", "304", "305", "306", "307", "308", "309", "310", "311", "312",  "314", "315", "316","317" "318", "319", "320", "321", "322", "323", "325", "326", "327", "328", "329", "330", "332", "333", "334", "335", "336", "337", "338", "339", "340", "341", "342", "343", "344", "345", "346", "347", "348", "349", "350", "351", "352", "353", "354", "355", "356", "357", "358", "359", "360", "361", "362", "363", "364", "365", "366", "367", "368", "369", "370", "371", "372", "373", "374", "375", "376", "377", "378", "379", "380", "381", "382", "383", "384", "385", "386", "387", "388", "389", "390", "391", "392", "393"]
adults = ["400", "401", "402", "403", "404", "405", "406", "407", "408", "409", "410", "411", "412", "413", "414", "415", "416"]

base_dir = '/gablab/p/GATES/Analysis/MCNAB_FIR'
task_periods = {'inst_dist':'_fir_model0/','inst_nodist':'_fir_model1/'}
for roi in glob(os.path.join(base_dir,'*')):
    roi_name = os.path.split(roi)[-1]
    print roi_name
    for contrast in sorted(task_periods.keys()):
        print contrast
        contrast_dir = task_periods[contrast]
        directory = os.path.join(os.path.join(roi,'fir_model/*'),contrast_dir)
        beta_files = glob(os.path.join(directory,'fir_avg_hrf.txt'))
        hrf = None
        kid_beta_files = []
        adult_beta_files = []
        for kid in kids:
            kid_files = [beta_file for beta_file in beta_files if beta_file.find(kid) != -1]
            kid_beta_files.extend(kid_files)
        for adult in adults:
            adult_files = [beta_file for beta_file in beta_files if beta_file.find(adult) != -1]
            adult_beta_files.extend(adult_files)

        #print kid HRF
        hrf = None
        for hrf_file in kid_beta_files:
            curve = np.genfromtxt(hrf_file)
            if max(curve) >10:
                print hrf_file
                print max(curve)

            #norm_curve = curve/np.sum(np.abs(curve))
            plot(np.arange(-1,10,1),curve,marker='o',linestyle='None', color = 'b', alpha=0.1)
            try:
                hrf = np.column_stack([hrf,curve])
            except ValueError:                
                hrf = curve
                print 'initialized hrf'
        group_hrf = np.median(hrf,axis=1)
        #norm_group_hrf = group_hrf/np.sum(np.abs(group_hrf))
        plot(np.arange(-1,10,1),group_hrf, linestyle = '-', color='b', label='kids',linewidth=2)

        #print adult HRF
        hrf = None
        for hrf_file in adult_beta_files:
            curve = np.genfromtxt(hrf_file)
            if max(curve) >10:
                print hrf_file
                print max(curve)
            #norm_curve = curve/np.sum(np.abs(curve))
            plot(np.arange(-1,10,1),curve, marker='o', linestyle='None', color = 'r', alpha=0.1)
            try:
                hrf = np.column_stack([hrf,curve])
            except ValueError:                
                hrf = curve
                print 'initialized hrf'
        group_hrf = np.median(hrf,axis=1)
        #norm_group_hrf = group_hrf/np.sum(np.abs(group_hrf))
        plot(np.arange(-1,10,1),group_hrf, linestyle = '-', color='r', label='adults',linewidth=2)


        xlabel("TR")
        ylabel("B Val")
        
        title("%s %s"%(roi_name,contrast))
        legend()
        show()
