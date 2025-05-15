from os import path, mkdir
from os.path import join, isdir, basename
import nibabel as nib
import argparse
from numpy import inf
import numpy as np

# custom imports
from utils import get_image_paths, data_cleaning, remove_outliers

'''
Script for processing R maps. It combines the scripts `data_cleaner`, `apply_mask` and `outlier_remover` to reduce the number of commands that need to be executed. 
Note that both R1 and R2star data will be processed at the same time.
'''

# parse in and out folders through the command line
p = argparse.ArgumentParser(description='Process R maps.')
p.add_argument('-r', '--root', required=True, help='root folder')
p.add_argument('-f', '--folder', required=True, help='folder')
p.add_argument('-m', '--mask', required=True, help='name of the brain mask for skull-stripping (without the subject number)')
p.add_argument('-o', '--output', required=True, help='target folder for the processed T maps')
p.add_argument('-s', '--subject', required=True, nargs='*', help='subject number(s)')
args = p.parse_args()

# define target folder where the processed data will be stored as well as folder for met data
target_folder = args.output

# make sure that directory exists
if not isdir(target_folder):
    mkdir(target_folder)

for subject in args.subject:
    # get folders and paths
    input_folder = args.root+subject+args.folder
    
    # get all nifti files from the input folder
    nifti_files = get_image_paths(input_folder, ['R1','R2s']) 

    # get mask from given path
    mask_path = args.mask+subject+'.nii'
    mask = nib.load(mask_path).get_fdata()

    # read in nifti files and process R maps
    for nifti_file_path in nifti_files:
        # get image data for the volume and affine transformation for saving later
        nifti_file = nib.load(nifti_file_path)
        volume = nifti_file.get_fdata()

        # extract brain using the mask
        volume[mask == 0] = 0

        # remove NaN, zeros, etc. from volume
        volume = data_cleaning(volume)

        # remove outliers from the volume using different thresholds for R1 and R2star
        if nifti_file_path.find('R1') != -1:
            volume = remove_outliers(volume,threshold=0.999)
        elif nifti_file_path.find('R2s') != -1:
            volume = remove_outliers(volume,threshold=0.995)
        
        # turn the volume into a Nifti image using the affine transform and save the data
        save_vol = nib.Nifti1Image(volume, nifti_file.affine)
        print('Saving ',basename(nifti_file_path))
        nib.save(save_vol, join(target_folder, basename(nifti_file_path)))

