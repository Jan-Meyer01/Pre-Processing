from os import path, mkdir
from os.path import join, isdir, basename
import nibabel as nib
import argparse

# custom imports
from utils import get_image_paths, data_cleaning, remove_outliers, skull_stripping, process_volume_data

# parse in and out folders through the command line
p = argparse.ArgumentParser(description='Pre-processing pipeline (data cleaning, skull-stripping) for 3D nifti images.')
p.add_argument('-i', '--input', default=join('/projects','crunchie','Jan','Daten','DataSuperresolution','DWI'), help='input folder with nifti images to process')
p.add_argument('-o', '--output', default=join('.','processed_data'), help='output folder to save processed nifti images')
p.add_argument('-t', '--threshold', default=0.9975, type=float, help='threshold in interval [0,1] for removing outliers (if selected in mode)')
p.add_argument('-m', '--mode', default=0, help='mode of preprocessing: everything (0), just cleaning (1), just outlier-removal (2), just skull-stripping (3), cleaning + skull-stripping (4)')
args = p.parse_args()

# define target folder where the processed data will be stored
target_folder = args.output

# make sure that directory exists
if not isdir(target_folder):
    mkdir(target_folder)

# get all nifti files from the input folder
nifti_files = get_image_paths(args.input, ['PD']) #,'R1','R2s'

# select pre-processing steps
mode = int(args.mode)
if mode == 0:
    functions = [data_cleaning,remove_outliers,skull_stripping]
elif mode == 1:
    functions = [data_cleaning]
elif mode == 2:
    functions = [remove_outliers]
elif mode == 3:
    functions = [skull_stripping]
elif mode == 4:
    functions = [data_cleaning,skull_stripping]
else:
    raise ValueError('No mode, aka pre-processing steps, chosen!!')   

# read in nifti files and process them
for nifti_file_path in nifti_files:
    # get image data for the volume and affine transform
    nifti_file = nib.load(nifti_file_path)
    volume = nifti_file.get_fdata()
    affine = nifti_file.affine

    # apply processing
    for function in functions:
        if function == remove_outliers:
            assert args.threshold < 1 and args.threshold > 0
            volume = function(volume, args.threshold)
        else:
            volume = function(volume)
    #vol_processed = process_volume_data(volume, mode=args.mode, debug=False)

    # save the processed volume to the target folder
    volume_save = nib.Nifti1Image(volume, affine)

    # set name for saving depending on the pre-processing
    if mode == 1:
        save_name = basename(nifti_file_path).replace(".nii", "_cleaned.nii")
    elif mode == 2:
        save_name = basename(nifti_file_path).replace(".nii", "_outlierRemoved.nii")
    elif mode == 3:
        save_name = basename(nifti_file_path).replace(".nii", "_stripped.nii")
    else:
        save_name = basename(nifti_file_path).replace(".nii", "_processed.nii")    
    print('Saving {} as {}'.format(basename(nifti_file_path),save_name))
    nib.save(volume_save, join(target_folder, save_name))
