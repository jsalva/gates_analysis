def corr_image(resting_image, aparc_aseg_file,fwhm, seed_region):
    """This function makes correlation image on brain surface"""
    import numpy as np
    import nibabel as nb
    import matplotlib.pyplot as plt
    from surfer import Brain, Surface
    import os
    import string
    aparc_aseg = nb.load(aparc_aseg_file)
    img = nb.load(resting_image)
    corrmat = np.corrcoef(np.squeeze(img.get_data()))
    corrmat[np.isnan(corrmat)] = 0
    corrmat_npz = os.path.abspath('corrmat.npz')
    np.savez(corrmat_npz,corrmat=corrmat)

#    br = Brain('fsaverage5', 'lh', 'smoothwm')

    #br.add_overlay(corrmat[0,:], min=0.2, name=0, visible=True)
    #values = nb.freesurfer.read_annot('/software/Freesurfer/5.1.0/subjects/fsaverage5/label/lh.aparc.annot')
#    values = open('/software/Freesurfer/current/FreeSurferColorLUT.txt').read()
#    values = string.split(values,'\n')
#    values = filter(None,map(string.strip,values))


    #br.add_overlay(np.mean(corrmat[values[0]==5,:], axis=0), min=0.8, name='mean', visible=True)

    aparc_aseg_data = np.squeeze(aparc_aseg.get_data())
#    data = img.get_data()

    data = np.squeeze(img.get_data())
    

    
    seed_signal = np.mean(data[aparc_aseg_data==seed_region], axis=0)
    seed = np.corrcoef(seed_signal, data)

    plt.hist(seed[0,1:], 128)
    plt.savefig(os.path.abspath("histogram_%d.png"%seed_region))
    plt.close()

        #corr_image = os.path.abspath("corr_image%s.png"%fwhm)
        #br.save_montage(corr_image)
        #ims = br.save_imageset(prefix=os.path.abspath('fwhm_%s'%str(fwhm)),views=['medial','lateral','caudal','rostral','dorsal','ventral'])
        #br.close()
        #print ims
        #precuneus[np.isnan(precuneus)] = 0
        #plt.hist(precuneus[0,1:])

    roitable = [['Region','Mean Correlation']]
    for i, roi in enumerate(np.unique(aparc_aseg_data)):
        roitable.append([roi,np.mean(seed[aparc_aseg_data==seed_region])])

        #images = [corr_image]+ims+[os.path.abspath("histogram.png"), roitable]
    roitable=[roitable]
    histogram = os.path.abspath("histogram_%d.png"%seed_region)

    return corr_image, ims, roitable, histogram, corrmat_npz
