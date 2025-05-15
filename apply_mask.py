from os import path, mkdir
from os.path import join, isdir, basename
import nibabel as nib
import argparse

# custom imports
from utils import get_image_paths

# parse in and out folders through the command line
p = argparse.ArgumentParser(description='Pre-processing pipeline (data cleaning, skull-stripping) for 3D nifti images.')
p.add_argument('-i', '--input', required=True, help='input folder with nifti images to process')
p.add_argument('-o', '--output', required=True, help='output folder to save processed nifti images')
p.add_argument('-m', '--mask', required=True, help='brain mask for skull-stripping')
p.add_argument('-c', '--contrast', required=True, nargs='*', help='contrast(s)')
p.add_argument('-s', '--subject', required=True, nargs='*', help='subject number(s)')
args = p.parse_args()

# define target folder where the processed data will be stored
if type(args.output) == type(None):
    # if no output folder is given define target folder = input folder
    target_folder = args.input
else:    
    # else use given output folder as target folder
    target_folder = args.output

# make sure that directory exists
if not isdir(target_folder):
    mkdir(target_folder)

for subject in args.subject:
    # get all nifti files from the input folder
    nifti_files = get_image_paths(args.input, args.contrast) 

    # get mask from given path
    mask_path = args.mask+subject+'.nii'
    mask = nib.load(mask_path).get_fdata()

    # read in nifti files and extract the brain
    for nifti_file_path in nifti_files:
        # get image data for the volume and affine transformation for saving later
        nifti_file = nib.load(nifti_file_path)
        volume = nifti_file.get_fdata()

        # extract brain using the mask
        volume[mask == 0] = 0
        
        # turn the volume into a Nifti image using the affine transform and save the data
        save_vol = nib.Nifti1Image(volume, nifti_file.affine)
        print('Saving ',basename(nifti_file_path))
        nib.save(save_vol, join(target_folder, basename(nifti_file_path)))

