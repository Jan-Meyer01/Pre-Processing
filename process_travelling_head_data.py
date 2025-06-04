from os import path, scandir, makedirs
from os.path import join, isdir, basename
import nibabel as nib
import argparse
from numpy import inf
import numpy as np
import time

# custom imports
from utils import is_image_file, process_and_convert_R_to_T_travelling_head, process_and_convert_T_to_R_travelling_head

'''
Script for pre-processing the MPM and T2 data from the travelling head study
'''

# parse in and out folders through the command line
p = argparse.ArgumentParser(description='Script for pre-processing the MPM and T2 data from the travelling head study.')
p.add_argument('-i', '--input_dir', type=str, default="/projects/crunchie/Jan/Daten/DataTravellingHeadStudy/MPMs", help='input directory where all the data is stored in a BIDS file structure.') # required=True,
args = p.parse_args()

# get files from the input directory and define output directory to recreate BIDS structure under ./processed_data/travelling_head
input_dir=args.input_dir
output_dir=join(".","processed_data","travelling_head")

# read and recreating the BIDS structure of the travelling head data, then process the found nifti images
print('Started processing at ', time.strftime("%H:%M:%S", time.localtime()))  

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
                    elif file.endswith('T2.nii'):
                        T2_image  = nib.load(file).get_fdata()
                        T2_affine = nib.load(file).affine
                        T2_path   = file
                    elif file.endswith('brainmask.nii'):
                        brain_mask = nib.load(file).get_fdata()
            
            # check whether images and brain mask exist
            if type(PD_image) == type(None):
                print('Warning: No PD image found in ',input_path) 
            if type(R1_image) == type(None):
                print('Warning: No R1 image found in ',input_path)
            if type(R2star_image) == type(None):
                print('Warning: No R2* image found in ',input_path) 
            if type(T2_image) == type(None):
                print('Warning: No T2 image found in ',input_path) 
            assert type(brain_mask) != type(None); 'Error: No brain mask found in {}'.format(input_path)
            
            # process PD
            if type(PD_image) != type(None):
                # set negative, inf and NaN values to 0
                PD_image[PD_image<0] = 0
                PD_image[PD_image==inf] = 0
                PD_image[np.isnan(PD_image)] = 0
                # extract brain from PD image and save it
                PD_image[brain_mask==0] = 0
                save_vol = nib.Nifti1Image(PD_image, PD_affine) 
                save_name = join(output_dir, site, subject, session, basename(PD_path))
                nib.save(save_vol, save_name)

            # process R1 image and create T1 just in case
            if type(R1_image) != type(None):
                # get processed R1 image and converted T1 image
                R1_image, T1_image = process_and_convert_R_to_T_travelling_head(R_image=R1_image, brain_mask=brain_mask, cutoff=0.075)

                # TODO: remove outliers from images
                
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
                
                # TODO: remove outliers from images

                # save the R2* map
                save_vol  = nib.Nifti1Image(R2star_image, R2star_affine) 
                save_name = join(output_dir, site, subject, session, basename(R2star_path))
                nib.save(save_vol, save_name)
                
                # save the T2* map (just use the affine from R2*...)
                save_vol  = nib.Nifti1Image(T2star_image, R2star_affine) 
                save_name = join(output_dir, site, subject, session, basename(R2star_path).replace("R2s_OLS.nii", "T2s_OLS.nii"))
                nib.save(save_vol, save_name)
            
            # process T2 image and create R2 just in case
            if type(T2_image) != type(None):
                # get processed T2 image and converted R2 image
                T2_image, R2_image = process_and_convert_T_to_R_travelling_head(T_image=T2_image, brain_mask=brain_mask)
            
                # TODO: remove outliers from images
                
                # save the T2 map
                save_vol  = nib.Nifti1Image(T2_image, T2_affine) 
                save_name = join(output_dir, site, subject, session, basename(T2_path))
                nib.save(save_vol, save_name)
                
                # save the R2 map (just use the affine from T2...)
                save_vol  = nib.Nifti1Image(R2_image, T2_affine) 
                save_name = join(output_dir, site, subject, session, basename(T2_path).replace("T2.nii", "R2.nii"))
                nib.save(save_vol, save_name)

print('Finished processing at ', time.strftime("%H:%M:%S", time.localtime()))  