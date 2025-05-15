from os import path, mkdir
from os.path import join, basename, isdir
import nibabel as nib

# custom imports
from utils import process_volume_data, test_processing_pipelines

# define nifti file for test
nifti_file_path = '/projects/crunchie/Jan/Daten/DataSuperresolution/DWI Clean/smi_DTI_FA.nii.gz'

# define target folder where the processed data will be stored
target_folder = join('.','processed_data')

# make sure that directory exists
if not isdir(target_folder):
    mkdir(target_folder)

# get image data for the volume and affine transform
nifti_file = nib.load(nifti_file_path)
volume = nifti_file.get_fdata()
affine = nifti_file.affine

# print min and max value for debugging and checking the values later
print('Min unprocessed: ',volume.min()) 
print('Max unprocessed: ',volume.max()) 

# apply processing in debugging mode
vol_clean, vol_skull_stripped, vol_no_outliers = process_volume_data(volume, debug=True)

# save cleaned volume to the target folder
vol_clean_save = nib.Nifti1Image(vol_clean, affine)
nib.save(vol_clean_save, join(target_folder, basename(nifti_file_path).replace(".nii", f"_clean.nii")))

# print min and max value for debugging and checking the values later
print('Min cleaned: ',volume.min()) 
print('Max cleaned: ',volume.max()) 

# save skull-stripped volume to the target folder
vol_skull_stripped_save = nib.Nifti1Image(vol_skull_stripped, affine)
nib.save(vol_skull_stripped_save, join(target_folder, basename(nifti_file_path).replace(".nii", f"_skull_stripped.nii")))

# print min and max value for debugging and checking the values later
print('Min skull-stripped: ',volume.min()) 
print('Max skull-stripped: ',volume.max()) 

# save the volume with not outliers to the target folder
vol_no_outliers_save = nib.Nifti1Image(vol_no_outliers, affine)
nib.save(vol_no_outliers_save, join(target_folder, basename(nifti_file_path).replace(".nii", f"_no_outliers.nii")))

# print min and max value for debugging and checking the values later
print('Min no-outliers: ',volume.min()) 
print('Max no-outliers: ',volume.max()) 