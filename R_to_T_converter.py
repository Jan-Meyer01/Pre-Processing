from os import path, mkdir
from os.path import join, isdir, basename
import nibabel as nib
import argparse
from numpy import inf
import numpy as np

# custom imports
from utils import get_image_paths

# parse in and out folders through the command line
p = argparse.ArgumentParser(description='Convert R maps to T maps.')
p.add_argument('-i', '--input', required=True, help='input nifti image e.g. R1 map')
p.add_argument('-o', '--output', required=True, help='target folder for e.g. T1 map')
args = p.parse_args()

# define target folder where the processed data will be stored
target_folder = args.output

# make sure that directory exists
if not isdir(target_folder):
    mkdir(target_folder)

# get all nifti files from the input folder
nifti_files = get_image_paths(args.input, ['R1','R2s']) 

# read in nifti files and process them
for nifti_file_path in nifti_files:
    #if nifti_file_path.find('cleaned_stripped_outlierRemoved') != -1:
    # get image data for the volume and affine transform
    nifti_file = nib.load(nifti_file_path)
    volume = nifti_file.get_fdata()
    affine = nifti_file.affine

    # thresholding for feasable values (we know which T values to expect in human tissue)
    volume[volume<0.1] = np.nan
    
    # invert the map
    volume = 1/volume
    
    # scale because values in ms are needed for JEMRIS simulation
    volume = volume * 1000
    
    # save the processed volume to the target folder
    volume_save = nib.Nifti1Image(volume, affine)
    
    # set save name
    if nifti_file_path.find('R1') != -1:
        save_name = basename(nifti_file_path).replace("R1", "T1")   
    elif nifti_file_path.find('R2s') != -1:
        save_name = basename(nifti_file_path).replace("R2s", "T2s")   
    print('Saving {} as {}'.format(basename(nifti_file_path),save_name))
    nib.save(volume_save, join(target_folder, save_name))
