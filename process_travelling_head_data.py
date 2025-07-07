from os import path, scandir, makedirs
from os.path import join, isdir, basename
import nibabel as nib
import argparse
from numpy import inf
import numpy as np
import time
import SimpleITK

# custom imports
from utils import is_image_file, process_and_convert_R_to_T_travelling_head, clip_outliers, find_replace_R2_T2_inconsistencies

'''
Script for pre-processing the MPM and T2 data from the travelling head study
'''

# parse in and out folders through the command line
p = argparse.ArgumentParser(description='Script for pre-processing the MPM and T2 data from the travelling head study.')
p.add_argument('-i', '--input_dir', type=str, required=True, help='input directory where all the data is stored in a BIDS file structure.') 
p.add_argument('-o', '--output_dir', type=str, default="./processed_data/travelling_head", help='output directory where the BIDS file structure will be established.') 
args = p.parse_args()

# get files from the input directory and recreate BIDS structure in the output directory
input_dir= args.input_dir #"/projects/crunchie/Jan/Daten/DataTravellingHeadStudy/MPMs"
output_dir=args.output_dir

# read and recreating the BIDS structure of the travelling head data, then process the found nifti images
#print('Started processing at ', time.strftime("%H:%M:%S", time.localtime()))  

# get sites
sites = [f.name for f in scandir(input_dir)]
# get subject and sessions for each site
for site in sites:
    print('  working on site:     ', site)
    site_path = join(output_dir, site)
    # check whether directory already exists else create new sites in the output directory
    if not isdir(site_path):
        makedirs(site_path)
    # and get subjects
    subjects = [f.name for f in scandir(join(input_dir, site))]
    
    for subject in subjects:
        print('    analyzing subject: ', subject)
        subject_path = join(output_dir, site, subject)
        # check whether directory already exists else create folders for subjects
        if not isdir(subject_path):
            makedirs(subject_path)
        # and get in sessions
        sessions = [f.name for f in scandir(join(input_dir, site, subject))]
        
        for session in sessions:
            print('        session:       ', session)
            session_path = join(output_dir, site, subject, session)
            # check whether directory already exists else create folders for subjects
            if not isdir(session_path):
                makedirs(session_path)
            # and get in sessions
            input_path = join(input_dir, site, subject, session, "Results")
            files = [f.path for f in scandir(input_path)]
            
            # init all the images and the brain mask for error tracking
            PD_image      = None
            PD_affine     = None
            PD_path       = None
            R1_image      = None
            R1_affine     = None
            R1_path       = None
            R2star_image  = None
            R2star_affine = None
            R2star_path   = None
            T2_image      = None
            T2_affine     = None
            T2_path       = None
            brain_mask    = None
            brain_seg     = None
            
            for file in files:
                # if the file is a nifti image
                if is_image_file(file):
                    # check whether it is brainmask, PD, R1, R2* or T2 and load accordingly
                    if file.endswith('PD.nii'):
                        PD_image  = nib.load(file).get_fdata()
                        PD_affine = nib.load(file).affine
                        PD_path   = file
                    elif file.endswith('R1.nii'):
                        R1_image  = nib.load(file).get_fdata()    
                        # correct for ms in Bonn-data
                        if site == "Bonn_Skyra_3T_LowRes":
                            R1_image = R1_image*1000
                        R1_affine = nib.load(file).affine
                        R1_path   = file
                    elif file.endswith('R2s_OLS.nii'):
                        R2star_image  = nib.load(file).get_fdata()
                        R2star_affine = nib.load(file).affine
                        R2star_path   = file
                    elif file.endswith('brainmask.nii'):
                        brain_mask = nib.load(file).get_fdata()
                    elif file.endswith('brain_segmentation.nii'):
                        brain_seg  = nib.load(file).get_fdata()                        
                    # overwrite generated T2 and R2 maps with SimpleITK as they have a problem when trying to open them with nibabel
                    elif file.endswith('T2_map.nii.gz'):
                        tmp = SimpleITK.ReadImage(file)
                        SimpleITK.WriteImage(tmp, file)
                        T2_image  = nib.load(file).get_fdata()
                        T2_affine = nib.load(file).affine
                        T2_path   = file
                    elif file.endswith('R2_map.nii.gz'):
                        tmp = SimpleITK.ReadImage(file)
                        SimpleITK.WriteImage(tmp, file)
                        R2_image  = nib.load(file).get_fdata()
                        R2_affine = nib.load(file).affine
                        R2_path   = file
            
            # check whether images and brain mask exist
            if type(PD_image) == type(None):
                print('Warning: No PD image found in ',input_path,' !!') 
            if type(R1_image) == type(None):
                print('Warning: No R1 image found in ',input_path,' !!')
            if type(R2star_image) == type(None):
                print('Warning: No R2* image found in ',input_path,' !!') 
            if type(T2_image) == type(None):
                print('Warning: No T2 image found in ',input_path,' !!')  
            if type(R2_image) == type(None):
                print('Warning: No R2 image found in ',input_path,' !!')  
            if type(brain_seg) == type(None):
                print('Warning: No brain segmentation found in ',input_path,' !!')  
            assert type(brain_mask) != type(None), 'No brain mask found in {} !!'.format(input_path)

            # check if mask is not zero for all voxels 
            if np.sum(np.sum(brain_mask))==0:
                print('Warning: Empty brain mask in ',input_path,' --> Check your data and masks again!!')
            
            # process PD
            if type(PD_image) != type(None):
                # set negative, inf and NaN values to 0
                PD_image[PD_image<0] = 0
                PD_image[PD_image==inf] = 0
                PD_image[np.isnan(PD_image)] = 0
                
                # extract brain from PD image and clip outliers
                PD_image[brain_mask==0] = 0
                PD_image = clip_outliers(vol=PD_image,threshold=0.999)

                # save the PD map
                save_vol = nib.Nifti1Image(PD_image, PD_affine) 
                save_name = join(output_dir, site, subject, session, basename(PD_path))
                nib.save(save_vol, save_name)

            # process R1 image and create T1 just in case
            if type(R1_image) != type(None):
                # get processed R1 image and converted T1 image
                R1_image, T1_image = process_and_convert_R_to_T_travelling_head(R_image=R1_image, brain_mask=brain_mask, cutoff=0.075)

                # clip outliers (empirical values)
                R1_image = clip_outliers(vol=R1_image,threshold=0.9975)
                T1_image = clip_outliers(vol=T1_image,threshold=0.99)
                
                # save the R1 map
                save_vol  = nib.Nifti1Image(R1_image, R1_affine) 
                save_name = join(output_dir, site, subject, session, basename(R1_path))
                nib.save(save_vol, save_name)
                
                # save the T1 map (just use the affine from R1...)
                save_vol  = nib.Nifti1Image(T1_image, R1_affine) 
                save_name = join(output_dir, site, subject, session, basename(R1_path).replace("R1.nii", "T1.nii"))
                nib.save(save_vol, save_name)
                
            # process R2* image and create T2* just in case
            if type(R2star_image) != type(None):
                # get processed R2* image and converted T2* image
                R2star_image, T2star_image = process_and_convert_R_to_T_travelling_head(R_image=R2star_image, brain_mask=brain_mask, cutoff=0.025)
                
                # clip outliers 
                R2star_image = clip_outliers(vol=R2star_image,threshold=0.985)
                T2star_image = clip_outliers(vol=T2star_image,threshold=0.98)

                # save the R2* map
                save_vol  = nib.Nifti1Image(R2star_image, R2star_affine) 
                save_name = join(output_dir, site, subject, session, basename(R2star_path))
                nib.save(save_vol, save_name)
                
                # save the T2* map (just use the affine from R2*...)
                save_vol  = nib.Nifti1Image(T2star_image, R2star_affine) 
                save_name = join(output_dir, site, subject, session, basename(R2star_path).replace("R2s_OLS.nii", "T2s_OLS.nii"))
                nib.save(save_vol, save_name)
            
            # process and save T2/R2 image (both are generated using the same script so we assume both exist)
            if type(R2star_image) != type(None) and type(T2_image) != type(None) and type(brain_seg) != type(None):
                # clean images, but no masking required as they are already skull-stripped
                T2_image[T2_image<0] = 0
                T2_image[T2_image==inf] = 0
                T2_image[np.isnan(T2_image)] = 0
                R2_image[R2_image<0] = 0
                R2_image[R2_image==inf] = 0
                R2_image[np.isnan(R2_image)] = 0

                # clip outliers
                T2_image = clip_outliers(vol=T2_image,threshold=0.99)
                R2_image = clip_outliers(vol=R2_image,threshold=0.99)

                # find and replace inconsistent datapoints for T2 and R2
                T2_image, R2_image = find_replace_R2_T2_inconsistencies(T2_image,T2star_image,R2_image,R2star_image,brain_seg,0.02,4)

                # save the T2 map
                save_vol  = nib.Nifti1Image(T2_image, T2_affine) 
                save_name = join(output_dir, site, subject, session, "{}_{}_".format(subject,session)+basename(T2_path))
                nib.save(save_vol, save_name)

                # save the R2 map 
                save_vol  = nib.Nifti1Image(R2_image, R2_affine) 
                save_name = join(output_dir, site, subject, session, "{}_{}_".format(subject,session)+basename(R2_path))
                nib.save(save_vol, save_name)

#print('Finished processing at ', time.strftime("%H:%M:%S", time.localtime()))  