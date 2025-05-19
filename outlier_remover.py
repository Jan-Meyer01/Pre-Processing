from os import path, mkdir
from os.path import join, isdir, basename
import nibabel as nib
import argparse

# custom imports
from utils import get_image_paths, clip_outliers

# parse in and out folders through the command line
p = argparse.ArgumentParser(description='Pre-processing pipeline (data cleaning, skull-stripping) for 3D nifti images.')
p.add_argument('-i', '--input', required=True, help='input folder with nifti images for outlier removal')
p.add_argument('-o', '--output', required=True, help='output folder to save processed nifti images to')
p.add_argument('-c', '--contrast', required=True, help='choose the contrast for which all images in the folder will be processed (only one at a time)')
p.add_argument('-t', '--threshold', default=0.9975, type=float, help='threshold in interval [0,1] for removing outliers (if selected in mode)')
args = p.parse_args()

# define target folder where the processed data will be stored
target_folder = args.output

# make sure that directory exists
if not isdir(target_folder):
    mkdir(target_folder)

# get all nifti files from the input folder
nifti_files = get_image_paths(args.input, [args.contrast]) 

# read in nifti files and process them
for nifti_file_path in nifti_files:
    if nifti_file_path.find('cleaned') != -1:    
        # get image data for the volume and affine transform
        nifti_file = nib.load(nifti_file_path)
        volume = nifti_file.get_fdata()
        affine = nifti_file.affine

        # remove outliers
        assert args.threshold < 1 and args.threshold > 0
        volume = clip_outliers(volume, args.threshold)

        # save the processed volume to the target folder
        volume_save = nib.Nifti1Image(volume, affine)
        
        # set save name
        save_name = basename(nifti_file_path).replace(".nii", "_outlierRemoved.nii")   
        print('Saving {} as {}'.format(basename(nifti_file_path),save_name))
        nib.save(volume_save, join(target_folder, save_name))
