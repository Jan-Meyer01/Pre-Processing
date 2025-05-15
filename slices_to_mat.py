from os import path, makedirs
from os.path import join, isdir, basename
import nibabel as nib
import argparse
from scipy.io import savemat
import numpy as np
from numpy import inf
import sys

# custom imports
from utils import get_image_paths

"""
Script for turning MPMs into .mat files as phantoms for the JEMRIS simulation toolbox.
Make sure that the input folder only contains PD, R1, R2 and R2star images or volumes with the same size.
"""

# parse in and out folders through the command line
p = argparse.ArgumentParser(description='Convert nifti volumes to .mat dict with slices.')
p.add_argument('-i', '--input', required=True, help='input folder with nifti images to process')
p.add_argument('-o', '--output', required=True, help='output folder to save processed nifti images')
p.add_argument('-s', '--subject', required=True, nargs='*', help='subject number(s)')
args = p.parse_args()

# define folder where the pre-processed input data is stored
data_folder = args.input

# define target folder where the .mat slices will be stored
target_folder = args.output

# make sure that directory exists
if not isdir(target_folder):
    makedirs(target_folder)

# get all nifti files from the folder
nifti_files = get_image_paths(data_folder, ['PD','T1','T2s'])

# for every subject
for subject in args.subject:
    # init values for error handling
    M0 = None
    T1 = None
    T2 = None
    T2star = None
    
    # read in nifti files and process them
    for nifti_file_path in nifti_files:
        if basename(nifti_file_path).find('sub-tle{}'.format(subject)) != -1:
            # get image data for the volume and affine transform
            nifti_file = nib.load(nifti_file_path)
            volume = nifti_file.get_fdata()

            if basename(nifti_file_path).find('PD') != -1: #_cleaned_stripped
                M0 = volume
            elif basename(nifti_file_path).find('T1') != -1: #_cleaned_stripped_outlierRemoved
                T1 = volume
            elif basename(nifti_file_path).find('T2s') != -1: #_OLS_cleaned_stripped_outlierRemoved
                T2star = volume  
                T2 = volume 
            #elif basename(nifti_file_path).find('T2') != -1: #_cleaned_stripped_outlierRemoved
            #    T2 = volume 

    # set other parameters
    db = np.zeros_like(M0)                  # zeros for db
    offset = [0,0,0]                        # and offset

    RESOS = nifti_file.header.get_zooms()   # resolution in all direction 
    RES = int(RESOS[0])                     # as well as resolution

    if RES == 0:
        RES = 1
    
    if type(M0) == type(None):
        M0 = np.zeros_like(T1)
    elif type(T1) == type(None):
        T1 = np.zeros_like(M0)
    elif type(T2) == type(None):
        T2 = np.zeros_like(T1)
    elif type(T2star) == type(None):
        T2star = np.zeros_like(T1)

    assert M0.shape == T1.shape == T2star.shape == T2.shape

    for slice_num in range(M0.shape[2]):
        # save all the MPMs in a single .mat file
        save_name = 'sub-tle{0}_MPM_Slice{1:03d}.mat'.format(subject,slice_num+1)
        sys.stdout.write("\rsub-tle{} - Slice{}".format(subject,slice_num+1))
        sys.stdout.flush()
        savemat(join(target_folder, save_name), {'M0':M0[:,:,slice_num],'T1':T1[:,:,slice_num],'T2s':T2star[:,:,slice_num],'T2':T2[:,:,slice_num],'DB':db[:,:,slice_num],'OFFSET':offset,'RES':RES,'RESOS':RESOS,'Shape':[M0.shape[0],M0.shape[1]]})
print(" ")       