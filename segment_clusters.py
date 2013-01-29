import os
import shutil
import time
from nipype.utils.logger import logging, logger, fmlogger, iflogger
import numpy as np
from scipy.stats.stats import pearsonr, spearmanr
from scipy.stats import wilcoxon
import sklearn as sk
from sklearn.linear_model.base import BaseEstimator, RegressorMixin
import sklearn.metrics as skm
import sklearn.cross_validation as cv
import matplotlib
import nipype.pipeline.engine as pe
import sys

print sys.argv[1]
print sys.argv[2]
print sys.argv[3]
print sys.argv[4]
print len(sys.argv)
if not len(sys.argv) == 6:
    sys.stderr.write("How to call function: python report_clusters.py [BASE_DIR] [P_THRESHOLD] [Df] [MIN_CLUSTER_SIZE] [SEGMENTATION_FILE]")
    sys.exit("; Needs these five arguments")
#originally 'condition'
base_dir = sys.argv[1]
pthresh = float(sys.argv[2])
df = int(sys.argv[3])
min_cluster_size = int(sys.argv[4])
seg_file = sys.argv[5]

def get_coords(img, affine):
    coords = []
    labels = np.setdiff1d(np.unique(img.ravel()), [0])
    cs = []
    for label in labels:
        cs.append(np.sum(img==label))
    for label in labels[np.argsort(cs)[::-1]]:
        coords.append(np.dot(affine,
                             np.hstack((np.mean(np.asarray(np.nonzero(img==label)),
                                                axis = 1),
                                        1)))[:3].tolist())
    return coords

"""
def show_slices(img, coords=None, threshold=0.1, cmap=None, prefix=None,
                show_colorbar=None, formatter='%.2f'):

    from nipy.labs import viz
    from nibabel import load
    import pylab

    if cmap is None:
        cmap = pylab.cm.hot

    data, aff = img.get_data(), img.get_affine()
    anatimg = load('/mindhive/gablab/GATES/data/WMSTAT/201/niftis/mprag2.nii')
    anatdata, anataff = anatimg.get_data(), anatimg.get_affine()
    anatdata = anatdata.astype(np.float)
    anatdata[anatdata<10.] = np.nan
    outfile = 'cluster.svg'
    if prefix:
        outfile = '_'.join((prefix, outfile))

    outfile = os.path.join('figures', outfile)
    if coords is None:
        osl = viz.plot_map(np.asarray(data), aff, threshold=threshold,
                           cmap=cmap, black_bg=False)
        osl.frame_axes.figure.savefig(outfile, transparent=True)
    else:
        for idx,coord in enumerate(coords):
            outfile = 'cluster%02d' % idx
            if prefix:
                outfile = '_'.join((prefix, outfile))
            outfile = os.path.join('figures', outfile)
            osl = viz.plot_map(np.asarray(data), aff, anat=anatdata, anat_affine=anataff,
                               threshold=threshold, cmap=cmap,
                               black_bg=False, cut_coords=coord)
            if show_colorbar:
                cb = pylab.colorbar(pylab.gca().get_images()[1], cax=pylab.axes([0.4, 0.075, 0.2, 0.025]),
                         orientation='horizontal', format=formatter)
                cb.set_ticks([cb._values.min(), cb._values.max()])
                pylab.show()
            osl.frame_axes.figure.savefig(outfile+'.svg', bbox_inches='tight', transparent=True)
            osl.frame_axes.figure.savefig(outfile+'.png', dpi=600, bbox_inches='tight', transparent=True)

# <codecell>

def plot_regression_line(x,y, xlim, color='r'):
    model=sk.linear_model.LinearRegression().fit(x[:,None],y)
    xplot = np.arange(xlim[0], xlim[1])[:,None]
    plot(xplot, model.predict(xplot), color=color)
"""
# <codecell>

import os
from scipy.ndimage import label
import scipy.stats as ss
def get_labels(data, image, min_extent,seg_vol,lut):
    labels, nlabels = label(data)
    peaks = []
    N = []
    rois = []
    nums = []
    for idx in range(1, nlabels+1):
        if sum(sum(sum(labels==idx)))<min_extent:
            labels[labels==idx] = 0
        else:
            peaks.append(max(image[labels==idx]))
            N.append(len(image[labels==idx]))
            cluster_seg_vals = seg_vol[labels==idx]
            roilist = []
            numlist = []
            for region in np.unique(cluster_seg_vals):
                numlist.append(len(cluster_seg_vals[cluster_seg_vals==region]))
                roilist.append(lut[region])
            rois.append(roilist)
            nums.append(numlist)
    return labels, nlabels, peaks, N, rois, nums

def parse_lut(lut_file='/software/Freesurfer/5.1.0/FreeSurferColorLUT.txt'):
    import re
    a=re.compile('\\S*[^ ]')
    b = open(lut_file)
    lut=dict()
    try:
        while(1):
            c=a.findall(b.next())
            if len(c)>2 and c[0][0]!= '#':
                lut[int(c[0])]=c[1]
    except StopIteration:
        print 'done'
    return lut



# <codecell>
import glob
from nibabel import load
file_list = glob.glob(os.path.join(base_dir,'spmT*.img')) 
lut = parse_lut()
#print lut.keys()
segfile = load(seg_file)
seg_vol=segfile.get_data()
for filename in file_list:
    print 'Getting coords for '+filename
    img=load(filename)
    data = img.get_data()
    max_labels, max_nlabels, max_peaks, max_N ,max_rois,max_nums= get_labels(data>ss.t.ppf(1-pthresh,df), data, min_cluster_size,seg_vol,lut)
    data[max_labels==0] = 0
    #cmeans = get_clustermeans(X, labels, nlabels)
    max_coords = get_coords(max_labels, img.get_affine())
    print 'Maxima:'
    string = ''
    for i in range(len(max_coords)):
        string+=str(np.array(max_coords[i]).round(2))+'\t'
        string+=str(max_peaks[i])+'\t'
        string+=str(max_N[i])+'\t'
        print string
        string ='\t\t\t'
        for j in range(len(max_rois[i])):
            string+=str(max_rois[i][j])+'\t'
            string+=str(max_nums[i][j])+'\t'
            print string
            string = '\t\t\t'
        string = ''
#    show_slices(img, max_coords, threshold=0.5, prefix='uncorrected', show_colorbar=True)
#    print labels


    seg_vol = segfile.get_data()
    img=load(filename)
    inv_data = -1*img.get_data()
    min_labels, min_nlabels, min_peaks, min_N, min_rois,min_nums = get_labels(inv_data>ss.t.ppf(1-pthresh,df), inv_data, min_cluster_size,seg_vol,lut)
    inv_data[min_labels==0] = 0

    min_coords = get_coords(min_labels, img.get_affine())

    print 'Minima:'
    string = ''
    for i in range(len(min_coords)):
        string+=str(np.array(min_coords[i]).round(2))+'\t'
        string+=str(min_peaks[i])+'\t'
        string+=str(min_N[i])+'\t'
        print string
        string = '\t\t\t'
        for j in range(len(min_rois[i])):
            string+=str(min_rois[i][j])+'\t'
            string+=str(min_nums[i][j])+'\t'
            print string
            string = '\t\t\t'
        string = ''
