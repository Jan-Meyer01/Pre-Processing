from os import path, mkdir
from os.path import join, isdir, basename
import nibabel as nib
import argparse

# custom imports
from utils import get_image_paths, data_cleaning

# parse in and out folders through the command line
p = argparse.ArgumentParser(description='data cleaning for 3D nifti images.')
p.add_argument('-r', '--root', required=True, help='root folder')
p.add_argument('-f', '--folder', required=True, help='folder')
p.add_argument('-o', '--output', required=True, help='output folder to save processed nifti images')
p.add_argument('-c', '--contrast', required=True, help='choose the contrast for which all images in the folder will be processed (only one at a time)')
p.add_argument('-s', '--subject', required=True, nargs='*', help='subject number(s)')
args = p.parse_args()

# define target folder where the processed data will be stored
target_folder = args.output

# make sure that directory exists
if not isdir(target_folder):
    mkdir(target_folder)

for subject in args.subject:
    # get folders and paths
    input_folder = args.root+subject+args.folder
    
    # get all nifti files from the input folder
    nifti_files = get_image_paths(input_folder, [args.contrast])
    
    # read in nifti files and process them
    for nifti_file_path in nifti_files:
        # get image data for the volume and affine transform
        nifti_file = nib.load(nifti_file_path)
        volume = nifti_file.get_fdata()
        affine = nifti_file.affine

        # clean data
        volume = data_cleaning(volume)

        # save the processed volume to the target folder
        volume_save = nib.Nifti1Image(volume, affine)

        # set name for saving depending on the pre-processing
        save_name = basename(nifti_file_path)#.replace(".nii", "_cleaned.nii")   
        print('Saving {} as {}'.format(basename(nifti_file_path),save_name))
        nib.save(volume_save, join(target_folder, save_name))
