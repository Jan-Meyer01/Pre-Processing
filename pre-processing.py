from os import path, mkdir, system
from os.path import join, isdir, basename
import nibabel as nib
import argparse
from numpy import inf
import numpy as np
from scipy.io import savemat
from subprocess import call

# custom imports
from utils import get_image_paths, data_cleaning, clip_outliers, run_freesurfer

'''
Script for the complete pre-processing pipeline
'''

# parse in and out folders through the command line
p = argparse.ArgumentParser(description='Preprocess PD and R maps.')
p.add_argument('-r', '--root', required=True, help='root folder')
p.add_argument('-f', '--folder', required=True, help='folder')
p.add_argument('-oF', '--output_files', required=True, help='target folder for the processed maps')
p.add_argument('-oM', '--output_mat', required=True, help='target folder for the matlab files')
p.add_argument('-s', '--subject', required=True, nargs='*', help='subject number(s)')
args = p.parse_args()

# define target folder where the processed data will be stored as well as folder for mat data
target_folder_files = args.output_files
target_folder_mat = args.output_mat

# make sure that directories exist
if not isdir(target_folder_files):
    mkdir(target_folder_files)
if not isdir(target_folder_mat):
    mkdir(target_folder_mat)

# init arrays for PD, T1 and T2star maps (contain [volume, file_name, affine]) as well as resolution for mat files
PD_maps = []
T1_maps = []
T2star_maps = []
RES = []

# first clean the PD images
for subject in args.subject:
    # get folders and paths
    input_folder = args.root+subject+args.folder

    # get all PD nifti files from the input folder
    nifti_files = get_image_paths(input_folder, ['PD']) 

    # read in nifti files and convert R to T maps
    for nifti_file_path in nifti_files:
        # get image data for the volume and affine transformation for saving later
        nifti_file = nib.load(nifti_file_path)
        volume = nifti_file.get_fdata()
        PD_maps.append([data_cleaning(volume),nifti_file_path,nifti_file.affine])
        RES.append(int(nifti_file.header.get_zooms()[0]))

        # save PD file to use for synthstrip
        save_vol = nib.Nifti1Image(volume, nifti_file.affine) 
        save_name = join(target_folder_files, basename(nifti_file_path))
        nib.save(save_vol, save_name)

    # now call freesurfers synthstrip model to extract the brain mask and save it to the target folder
    mask_name = join(target_folder_files, 'brain_mask_sub-tle'+subject+'.nii')
    command = 'mri_synthstrip -i {} -o {} -m {}'.format(save_name,save_name,mask_name)
    call_output = run_freesurfer(command)

# get T maps
#for subject in args.subject:
    # get folders and paths
    #input_folder = args.root+subject+args.folder
    
    # get all R1 and R2s nifti files from the input folder
    nifti_files = get_image_paths(input_folder, ['R1','R2s']) 

    # get mask from given path
    #mask_name = args.mask+subject+'.nii'
    mask = nib.load(mask_name).get_fdata()

    # read in nifti files and convert R to T maps
    for nifti_file_path in nifti_files:
        # get image data for the volume and affine transformation for saving later
        nifti_file = nib.load(nifti_file_path)
        volume = nifti_file.get_fdata()

        # thresholding for feasable values (we know which T values to expect in human tissue)
        volume[volume<0.1] = np.nan
        
        # invert the map
        volume = 1/volume

        # extract brain using the mask
        volume[mask == 0] = 0
        
        # scale because values in ms are needed for JEMRIS simulation
        volume = volume * 1000

        # store volume data with file name and affine transform
        if nifti_file_path.find('R1') != -1:
            T1_maps.append([volume,nifti_file_path,nifti_file.affine])
        elif nifti_file_path.find('R2s') != -1:
            T2star_maps.append([volume,nifti_file_path,nifti_file.affine])    

# ensure that you have all the maps 
assert len(PD_maps) == len(T1_maps) == len(T2star_maps)

# set other parameters for mat files
db = np.zeros_like(T1_maps[0][0])               # zeros for db
offset = [0,0,0]                                # and offset

# clean data, remove outliers and save T maps as well as mat files
for num in range(len(T1_maps)):
    # for T1
    T1_maps[num][0] = clip_outliers(data_cleaning(T1_maps[num][0]),threshold=0.999)
    save_vol = nib.Nifti1Image(T1_maps[num][0], T1_maps[num][2]) 
    save_name = basename(T1_maps[num][1]).replace("R1", "T1") 
    print('Saving processed {} file as {}'.format(basename(T1_maps[num][1]),save_name))
    nib.save(save_vol, join(target_folder_files, save_name))

    # and for T2 star
    T2star_maps[num][0] = clip_outliers(data_cleaning(T2star_maps[num][0]),threshold=0.99)
    save_vol = nib.Nifti1Image(T2star_maps[num][0], T2star_maps[num][2])
    save_name = basename(T2star_maps[num][1]).replace("R2s", "T2s") 
    print('Saving processed {} file as {}'.format(basename(T2star_maps[num][1]),save_name))
    nib.save(save_vol, join(target_folder_files, save_name))

    # save all the MPMs in a single .mat file
    T2 = np.zeros_like(T1_maps[num][0])
    save_name = 'sub-tle{}_MPM.mat'.format(args.subject[num])
    print('saved values as {}'.format(save_name))
    savemat(join(target_folder_mat, save_name), {'M0':PD_maps[num][0],'T1':T1_maps[num][0],'T2S':T2star_maps[num][0],'T2':T2,'DB':db,'OFFSET':offset,'RES':RES}) 