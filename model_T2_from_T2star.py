from os import path, scandir
from os.path import join
import nibabel as nib
import argparse
import numpy as np
from scipy.optimize import curve_fit

# custom imports
from utils import is_image_file, find_replace_R2_T2_inconsistencies

'''
Find and resolve inconsistencies between T2 fitting and MPM data.
'''

# parse in and out folders through the command line
p = argparse.ArgumentParser(description='Script for finding and resolving inconsistencies between T2 fitting and MPM data.')
p.add_argument('-i', '--input_dir', type=str, default="./processed_data/travelling_head", help='input directory where all the data is stored in a BIDS file structure.') 
args = p.parse_args()

# get files from the input directory
input_dir=args.input_dir

# get sites
sites = [f.name for f in scandir(input_dir)]
# get subject and sessions for each site
for site in sites:
    print('  working on site:     ', site)
    # get subjects
    subjects = ["sub-phy001"] #[f.name for f in scandir(join(input_dir, site))]
    
    for subject in subjects:
        print('    analyzing subject: ', subject)
        # get in sessions
        sessions = ["ses-001"] #[f.name for f in scandir(join(input_dir, site, subject))]
        
        for session in sessions:
            print('        session:       ', session)
            input_path = join(input_dir, site, subject, session)
            # get files
            files = [f.path for f in scandir(input_path)]
                        
            # init all the images for error tracking
            T1_image     = None
            T2star_image = None
            T2_image     = None
            T2_affine    = None
            
            # get all the processed data
            for file in files:
                # if the file is a nifti image
                if is_image_file(file):
                    # load MPMs (T1+T2*) + T2 
                    if file.endswith('T1.nii'):
                        T1_image  = nib.load(file).get_fdata()
                    elif file.endswith('T2s_OLS.nii'):
                        T2star_image  = nib.load(file).get_fdata()
                    elif file.endswith('T2_map.nii.gz'):
                        T2_image  = nib.load(file).get_fdata()
                        T2_affine = nib.load(file).affine
                    elif file.endswith('R2s_OLS.nii'):
                        R2star_image  = nib.load(file).get_fdata()
                    elif file.endswith('R2_map.nii.gz'):
                        R2_image  = nib.load(file).get_fdata()
                        R2_affine = nib.load(file).affine
                    elif file.endswith('brain_segmentation.nii'):
                        brain_seg = nib.load(file).get_fdata()
            
            # check whether the images exist and have the same size
            assert type(T1_image)     != type(None), 'No T1 image found in {} !!'.format(input_path)
            assert type(T2star_image) != type(None), 'No T2* image found in {} !!'.format(input_path)
            assert type(T2_image)     != type(None), 'No T2 image found in {} !!'.format(input_path)
            assert type(R2star_image) != type(None), 'No R2* image found in {} !!'.format(input_path)
            assert type(R2_image)     != type(None), 'No R2 image found in {} !!'.format(input_path)
            assert T1_image.shape == T2star_image.shape == T2_image.shape == R2star_image.shape == R2_image.shape, 'Images have different sizes!!'

            # find and replace inconsistent datapoints for T2 and R2
            T2_image, R2_image = find_replace_R2_T2_inconsistencies(T2_image,T2star_image,R2_image,R2star_image,brain_seg,0.02,4)

            """
            # save difference map for analysis
            save_vol  = nib.Nifti1Image(diff_T2_T2star, T2_affine) 
            save_name = join(input_path, "{}_{}_difference_T2_T2*.nii".format(subject,session))
            nib.save(save_vol, save_name)
            """

            """
            # save image shape and flatten T2 and T2* images
            image_shape = T2_image.shape
            
            T2_image     = np.ndarray.flatten(T2_image)
            T2star_image = np.ndarray.flatten(T2star_image)
            brain_mask   = np.ndarray.flatten(brain_mask)
            
            

             
            # use T2* for linear model to get T2 estimate
            M = np.vstack([brain_mask,T2star_image]).T  # get matrix with mask and T2*
            
            # call function 
            a, b = np.linalg.lstsq(M, T2_image)[0]      # perform optimization
            T2_image_synth = a + b*T2star_image         # get estimated image
            
            # use T2* for linear model to get T2 estimate
            def linear_model(x,a,b):
                return a+b*x 
            model_coeff, _ = curve_fit(linear_model,T2_image,T2star_image)
            
            # apply model to T2* image and revert back to the corrrect image size
            T2_image_synth = linear_model(T2star_image, *model_coeff)

            # resize back to original image size
            brain_mask      = np.resize(brain_mask, image_shape)
            T2_image        = np.resize(T2_image, image_shape)
            T2_image_synth  = np.resize(T2_image_synth, image_shape)

            # mask out brain from T2*_adjusted 
            T2_image_synth = T2_image_synth*brain_mask

            # replace parts of the T2 image where data inconsistencies are
            T2_image[diff_T2_T2star==1] = T2_image_synth[diff_T2_T2star==1]
            """

            # save difference map for analysis
            save_vol  = nib.Nifti1Image(T2_image, T2_affine) 
            save_name = join(input_path, "{}_{}_T2_synth.nii".format(subject,session))
            nib.save(save_vol, save_name)
            save_vol  = nib.Nifti1Image(R2_image, R2_affine) 
            save_name = join(input_path, "{}_{}_R2_synth.nii".format(subject,session))
            nib.save(save_vol, save_name)
