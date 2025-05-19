from os import path, mkdir
from os.path import join, isdir, basename
import nibabel as nib
import argparse
from numpy import inf
import numpy as np

# custom imports
from utils import get_image_paths, data_cleaning, clip_outliers

'''
Script for generating T maps from R maps. It combines the scripts `R_to_T_converter`, `data_cleaner`, `apply_mask` and `outlier_remover` to reduce the number of commands that need to be executed. 
Note that both T1 and T2star data will be processed at the same time.
'''

# parse in and out folders through the command line
p = argparse.ArgumentParser(description='Convert R maps to T maps.')
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

# init arrays for T maps
T1_maps = []
T2star_maps = []

for subject in args.subject:
    # get folders and paths
    input_folder = args.root+subject+args.folder
    
    # get all nifti files from the input folder
    nifti_files = get_image_paths(input_folder, ['R1','R2s']) 

    # get mask from given path
    mask_path = args.mask+subject+'.nii'
    mask = nib.load(mask_path).get_fdata()

    # read in nifti files and convert R to T maps
    for nifti_file_path in nifti_files:
        # get image data for the volume and affine transformation for saving later
        nifti_file = nib.load(nifti_file_path)
        volume = nifti_file.get_fdata()

        # clean the volume
        volume = data_cleaning(vol=volume,replacement_value=1)

        # thresholding for feasable values (we know which T values to expect in human tissue)
        cutoff = 0.1
        volume[volume<cutoff] = np.nan
        
        # invert the map
        volume = 1/volume

        # extract brain using the mask
        volume[mask == 0] = 0

        # revert NaN values back with limited range
        volume[np.isnan(volume)] = 1/cutoff
        
        # scale because values in ms are needed for JEMRIS simulation
        volume = volume * 1000

        # store volume data with file name and affine transform
        if nifti_file_path.find('R1') != -1:
            T1_maps.append([volume,nifti_file_path,nifti_file.affine])
        elif nifti_file_path.find('R2s') != -1:
            T2star_maps.append([volume,nifti_file_path,nifti_file.affine])    

assert len(T1_maps) == len(T2star_maps)

# clean data, remove outliers and save T maps
for num in range(len(T1_maps)):
    # for T1
    save_vol = nib.Nifti1Image(clip_outliers(T1_maps[num][0],threshold=0.99), T1_maps[num][2])  
    save_name = basename(T1_maps[num][1]).replace("R1", "T1") 
    print('Saving processed {} file as {}'.format(basename(T1_maps[num][1]),save_name))
    nib.save(save_vol, join(target_folder, save_name))

    # and for T2 star
    save_vol = nib.Nifti1Image(clip_outliers(T2star_maps[num][0],threshold=0.99), T2star_maps[num][2]) 
    save_name = basename(T2star_maps[num][1]).replace("R2s", "T2s") 
    print('Saving processed {} file as {}'.format(basename(T2star_maps[num][1]),save_name))
    nib.save(save_vol, join(target_folder, save_name))
