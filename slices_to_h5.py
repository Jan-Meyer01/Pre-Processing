from os import path, makedirs
from os.path import join, isdir, basename
import nibabel as nib
import argparse
from scipy.io import savemat
import numpy as np
from numpy import inf
import sys
import h5py

# custom imports
from utils import get_image_paths

"""
Script for turning MPMs into .h5 files as phantoms for the JEMRIS simulation toolbox.
Make sure that the input folder only contains PD, R1, R2 and R2star images or volumes with the same size.
"""

# parse in and out folders through the command line
p = argparse.ArgumentParser(description='Convert nifti volumes to .h5 dict with slices for JEMRIS.')
p.add_argument('-i', '--input', required=True, help='input folder with nifti images to process')
p.add_argument('-o', '--output', required=True, help='output folder to save processed nifti images')
p.add_argument('-R', '--Rmaps', required=True, type=bool, help='include R maps in the phantom (True) or only T maps (False)')
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
nifti_files = get_image_paths(data_folder, ['PD','T1','T2s','R1','R2s'])

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

            if basename(nifti_file_path).find('PD') != -1: 
                M0 = volume
            elif basename(nifti_file_path).find('T1') != -1: 
                T1 = volume
            elif basename(nifti_file_path).find('T2s') != -1: 
                T2star = volume  
                T2 = volume+10      # add 10ms to the T2 sample
            elif basename(nifti_file_path).find('T2') != -1: 
                T2 = volume 
            elif basename(nifti_file_path).find('R1') != -1: 
                R1 = volume
            elif basename(nifti_file_path).find('R2s') != -1: 
                R2star = volume  
                R2 = volume+0.1     # add 1/10ms to R2 sample
            elif basename(nifti_file_path).find('R2') != -1: 
                R2 = volume 

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
    elif type(R1) == type(None):
        R1 = np.zeros_like(M0)
    elif type(R2) == type(None):
        R2 = np.zeros_like(R1)
    elif type(R2star) == type(None):
        R2star = np.zeros_like(R1)

    assert M0.shape == T1.shape == T2star.shape == T2.shape

    for slice_num in range(M0.shape[2]):
        # save all the MPMs in a single .h5 file
        save_name = 'sub-tle{0}_MPM_Slice{1:03d}.h5'.format(subject,slice_num+1)
        sys.stdout.write("\rsub-tle{} - Slice{}".format(subject,slice_num+1))
        sys.stdout.flush()

        # format data for the phantom
        data = np.zeros((M0.shape[1],M0.shape[0],5))            # first init data
        data[:,:,0] = M0[:,:,slice_num].transpose()             # add M0 (PD)  
        data[:,:,4] = db[:,:,slice_num].transpose()             # add DB (chemical shift --> set to zero)
        
        # if specified add R maps else use T maps
        if args.Rmaps:
            data[:,:,1] = R1[:,:,slice_num].transpose()         # add R1
            data[:,:,2] = R2[:,:,slice_num].transpose()         # add R2
            data[:,:,3] = R2star[:,:,slice_num].transpose()     # add R2star    
        else:
            data[:,:,1] = T1[:,:,slice_num].transpose()         # add T1
            data[:,:,2] = T2[:,:,slice_num].transpose()         # add T2
            data[:,:,3] = T2star[:,:,slice_num].transpose()     # add T2star            
        
        # create groups and datasets structure for JEMRIS phantom
        f = h5py.File(join(target_folder, save_name), "w")
        f.require_group(name='sample')
        f.create_dataset(name='sample/data', data=data)
        f.create_dataset(name='sample/offset', data=offset)
        f.create_dataset(name='sample/resolution', data=RESOS)
        f.close()
print(" ")       