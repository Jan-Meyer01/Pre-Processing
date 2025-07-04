from os import path, scandir, makedirs
from os.path import join, isdir
import nibabel as nib
import argparse
from numpy import inf
import numpy as np
import time
from scipy.io import savemat

# custom imports
from utils import is_image_file, save_HDF5

'''
Script for turning the MPM and T2 data from the travelling head study into phantoms for JEMRIS.
You should either choose Tmaps+mat for use in the GUI or Rmaps+h5 for a script.
'''

# parse in and out folders through the command line
p = argparse.ArgumentParser(description='Script for turning the MPM and T2 data from the travelling head study into phantoms for JEMRIS.')
p.add_argument('-i', '--input_dir', type=str, default="./processed_data/travelling_head", help='input directory where all the data is stored in a BIDS file structure.') 
p.add_argument('-f', '--format', type=str, default="mat", help='Choose the files format in which the phantoms will be saved: mat, h5.')
p.add_argument('-s', '--slice', action='store_true', default=False, help='Use this argument when you want to get slices instead of volumes.') 
p.add_argument('-R', '--Rmaps', action='store_true', default=False, help='Use this argument when you want to use R maps instead of T maps.')  
args = p.parse_args()

# get files from the input directory
input_dir=args.input_dir

if args.Rmaps:
    print('Using R Maps...')
else:   
    print('Using T Maps...')

# read and recreating the BIDS structure of the travelling head data, then process the found nifti images
#print('Started processing at ', time.strftime("%H:%M:%S", time.localtime()))  

# get sites
sites = [f.name for f in scandir(input_dir)]
# get subject and sessions for each site
for site in sites:
    print('  working on site:     ', site)
    # get subjects
    subjects = [f.name for f in scandir(join(input_dir, site))]
    
    for subject in subjects:
        print('    analyzing subject: ', subject)
        # get in sessions
        sessions = [f.name for f in scandir(join(input_dir, site, subject))]
        
        for session in sessions:
            print('        session:       ', session)
            input_path = join(input_dir, site, subject, session)
            # get files
            files = [f.path for f in scandir(input_path)]
            
            # define target folder for phantoms to go to
            target_folder = join(input_path,'phantoms_jemris')
            # if the folder does not yet exist create it
            if not isdir(target_folder):
                makedirs(target_folder)
            
            # init all the images for error tracking
            PD_image      = None
            R1_image      = None
            T1_image      = None
            R2star_image  = None
            T2star_image  = None
            T2_image      = None
            R2_image      = None
            
            # get all the processed data
            for file in files:
                # if the file is a nifti image
                if is_image_file(file):
                    # load MPMs + T2 and correct for ms as required by JEMRIS
                    if file.endswith('PD.nii'):
                        PD_image  = nib.load(file).get_fdata()
                        PD_file = file
                        # rescale PD image to intensity range [0,1]
                        PD_image = PD_image/np.max(np.max(PD_image))
                    elif file.endswith('R1.nii'):
                        R1_image  = nib.load(file).get_fdata()/1000
                    elif file.endswith('T1.nii'):
                        T1_image  = nib.load(file).get_fdata()*1000
                    elif file.endswith('R2s_OLS.nii'):
                        R2star_image  = nib.load(file).get_fdata()/1000
                    elif file.endswith('T2s_OLS.nii'):
                        T2star_image  = nib.load(file).get_fdata()*1000
                    elif file.endswith('T2_map.nii.gz'):
                        T2_image  = nib.load(file).get_fdata()*1000
                    elif file.endswith('R2_map.nii.gz'):
                        R2_image  = nib.load(file).get_fdata()/1000
            
            # check whether the images exist and have the same size
            assert type(PD_image) != type(None), 'No PD image found in {} !!'.format(input_path)
            if args.Rmaps:
                print('')
                assert type(R1_image)     != type(None), 'No R1 image found in {} !!'.format(input_path)
                assert type(R2star_image) != type(None), 'No R2* image found in {} !!'.format(input_path)
                assert type(R2_image)     != type(None), 'No R2 image found in {} !!'.format(input_path)
                assert PD_image.shape == R1_image.shape == R2star_image.shape == R2_image.shape, 'Images have different sizes!!'
            else:    
                assert type(T1_image)     != type(None), 'No T1 image found in {} !!'.format(input_path)
                assert type(T2star_image) != type(None), 'No T2* image found in {} !!'.format(input_path)
                assert type(T2_image)     != type(None), 'No T2 image found in {} !!'.format(input_path)
                assert PD_image.shape == T1_image.shape == T2star_image.shape == T2_image.shape, 'Images have different sizes!!'

            # create phantom out of slices
            db = np.zeros_like(PD_image)            # zeros for db
            offset = [0,0,0]                        # and offset

            RESOS = nib.load(PD_file).header.get_zooms()      # resolution in all direction 
            RES = int(RESOS[0])                     # as well as resolution

            if RES == 0:
                RES = 1
            
            if args.format == 'mat':
                # if slices are selected
                if args.slice:
                    for slice_num in range(PD_image.shape[1]):
                        # save all slices .mat files
                        save_name = '{0}_{1}_phantom_jemris_slice{2:03d}.mat'.format(subject,session,slice_num+1)
                        # if specified add R maps else use T maps
                        if args.Rmaps:
                            data = {'M0':np.flip(PD_image[:,slice_num,:].transpose(), axis=1),'R1':np.flip(R1_image[:,slice_num,:].transpose(), axis=1),'R2s':np.flip(R2star_image[:,slice_num,:].transpose(), axis=1),'R2':np.flip(R2_image[:,slice_num,:].transpose(), axis=1),'DB':db[:,slice_num,:].transpose(),'OFFSET':offset,'RES':RES,'Shape':[PD_image.shape[0],PD_image.shape[2]]}
                        else:
                            data = {'M0':np.flip(PD_image[:,slice_num,:].transpose(), axis=1),'T1':np.flip(T1_image[:,slice_num,:].transpose(), axis=1),'T2s':np.flip(T2star_image[:,slice_num,:].transpose(), axis=1),'T2':np.flip(T2_image[:,slice_num,:].transpose(), axis=1),'DB':db[:,slice_num,:].transpose(),'OFFSET':offset,'RES':RES,'Shape':[PD_image.shape[0],PD_image.shape[2]]}
                        # save .mat files to target folder
                        savemat(join(target_folder, save_name), data)
                # if there is nothing specified
                else:
                    # save the volumes to a single .mat file
                    save_name = '{0}_{1}_phantom_jemris_vol.mat'.format(subject,session)
                    # if specified add R maps else use T maps
                    if args.Rmaps:
                        data = {'M0':PD_image,'R1':R1_image,'R2s':R2star_image,'R2':R2_image,'DB':db,'OFFSET':offset,'RES':RES,'Shape':[PD_image.shape[0],PD_image.shape[1],PD_image.shape[2]]}
                    else:
                        data = {'M0':PD_image,'T1':T1_image,'T2s':T2star_image,'T2':T2_image,'DB':db,'OFFSET':offset,'RES':RES,'Shape':[PD_image.shape[0],PD_image.shape[1],PD_image.shape[2]]}
                    # save .mat file to target folder
                    savemat(join(target_folder, save_name), data)
            elif args.format == 'h5':
                # if slices are selected
                if args.slice:
                    for slice_num in range(PD_image.shape[1]):
                        # save all slices as .h5 files
                        save_name = '{0}_{1}_phantom_jemris_slice{2:03d}.h5'.format(subject,session,slice_num+1)
                        # format data for the phantom
                        data = np.zeros((PD_image.shape[2],PD_image.shape[0],5))            # first init data
                        data[:,:,0] = np.flip(PD_image[:,slice_num,:].transpose(), axis=1)  # add M0 (PD)  
                        data[:,:,4] = db[:,slice_num,:].transpose()                         # add DB (chemical shift --> set to zero)                    
                        # if specified add R maps else use T maps
                        if args.Rmaps:
                            data[:,:,1] = np.flip(R1_image[:,slice_num,:].transpose(), axis=1)         # add R1
                            data[:,:,2] = np.flip(R2_image[:,slice_num,:].transpose(), axis=1)         # add R2
                            data[:,:,3] = np.flip(R2star_image[:,slice_num,:].transpose(), axis=1)     # add R2*  
                        else:
                            data[:,:,1] = np.flip(T1_image[:,slice_num,:].transpose(), axis=1)         # add T1
                            data[:,:,2] = np.flip(T2_image[:,slice_num,:].transpose(), axis=1)         # add T2
                            data[:,:,3] = np.flip(T2star_image[:,slice_num,:].transpose(), axis=1)     # add T2* 
                        # create groups and datasets structure for JEMRIS phantom
                        save_HDF5(save_dir=join(target_folder, save_name), data=data, offset=offset, resolution=RESOS)      
                # if there is nothing specified
                else:
                    # save the volumes to a single .h5 file
                    save_name = '{0}_{1}_phantom_jemris_vol.h5'.format(subject,session)
                    # format data for the phantom
                    data = np.zeros((PD_image.shape[2],PD_image.shape[1],PD_image.shape[0],5))            # first init data
                    data[:,:,:,0] = PD_image.transpose()             # add M0 (PD)  
                    data[:,:,:,4] = db.transpose()             # add DB (chemical shift --> set to zero)                    
                    # if specified add R maps else use T maps
                    if args.Rmaps:
                        data[:,:,:,1] = R1_image.transpose()         # add R1
                        data[:,:,:,2] = R2_image.transpose()         # add R2
                        data[:,:,:,3] = R2star_image.transpose()     # add R2*    
                    else:
                        data[:,:,:,1] = T1_image.transpose()         # add T1
                        data[:,:,:,2] = T2_image.transpose()         # add T2
                        data[:,:,:,3] = T2star_image.transpose()     # add T2*  
                
                    # create groups and datasets structure for JEMRIS phantom
                    save_HDF5(save_dir=join(target_folder, save_name), data=data, offset=offset, resolution=RESOS)
        
#print('Finished processing at ', time.strftime("%H:%M:%S", time.localtime()))  