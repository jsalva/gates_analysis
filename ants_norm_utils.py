# Utility Functions ---------------------------------------------------------


def convert_affine(unwarped_brain, mean_func, out_fsl_file):
    """Converts fsl-style Affine registration into ANTS compatible itk format

    Parameters
    ----------
    unwarped_brain : structural reference image
    mean_func : image that was coregistered
    out_fsl_file : fsl-style coregistration matrix

    Returns
    -------
    file : returns the filename corresponding to the converted registration
    """
    import os
    cmd = "c3d_affine_tool -ref %s -src %s %s -fsl2ras \
-oitk fsl2antsAffine.txt" % (unwarped_brain, mean_func, out_fsl_file)
    os.system(cmd)
    return os.path.abspath('fsl2antsAffine.txt')


def get_image_dimensions(images):
    """Return dimensions of list of images

    Parameters
    ----------
    images : list of filenames

    Returns
    -------
    list : returns dimensions of input image list
    """
    import nibabel as nb

    if isinstance(images, list):
        dims = []
        for image in images:
            dims.append(len(nb.load(image).get_shape()))
    else:
        dims = len(nb.load(images).get_shape())
    return dims
