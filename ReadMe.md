# Pre-processing for nifti 3D images
TODO: All pre-processing steps can be done with a single command. This includes data cleaning, converting the R maps to T maps, skull-stripping and outlier removal. In order to use the processed images as phantoms for [JEMRIS](https://www.jemris.org/) the individual volumes need to be `.mat` files. Thus, after executing the command, you will have the PD, T1, T2star files in the specified folder and the `.mat` files of all the data for a single subject in potentially another. Please adjust the file paths to your system.

'''
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/pre-processing.py -r /projects/crunchie/Jan/Daten/DataSubjects/sub-tle -f /ses-preop/anat/Results -oF /home/janmeyer/Pre-Processing/processed_data/TLE/MPM  -oM /home/janmeyer/Pre-Processing/volumes_mat/TLE/MPM -s 001 002 005
'''

If you want to convert the phantoms from `.mat` to `.phantom` files to use in [Koma](https://juliahealth.org/KomaMRI.jl/stable/) you can do so with the following scripts (note that you need to have [Julia](https://julialang.org/downloads/) installed and change the first path to your installation). This can be either done for a whole volume (might take multiple GB of GPU to run later in the simulations) or just some slices.

'''
/home/janmeyer/.julia/juliaup/julia-1.11.4+0.x64.linux.gnu/bin/julia volume_mat_to_phantom.jl -i /home/janmeyer/Pre-Processing/volumes_mat/TLE/MPM/sub-tle -o /home/janmeyer/KomaSimulations/Phantoms/Volumes/sub-tle -s 001,002,005
/home/janmeyer/.julia/juliaup/julia-1.11.4+0.x64.linux.gnu/bin/julia slices_mat_to_phantom.jl -i /home/janmeyer/Pre-Processing/slices_mat/TLE/MPM/sub-tle -o /home/janmeyer/KomaSimulations/Phantoms/Slices/sub-tle -s 001,002,005
'''

## Example workflow for the TLE data
To get a better understanding of the pipeline, we can examine its main four steps (each can be executed as an independent script). 
First, we clean up all the PD image (see different subject numbers). The first command is for the TLE data and the second for the travelling head study data.

'''
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/data_cleaner.py -r /projects/crunchie/Jan/Daten/DataSubjects/sub-tle -f /ses-preop/anat/Results -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -c PD -s 001 002 005
'''

Then we extract the brains of the PDs (as these usually have the best quality) without any further processing using freesurfers synthstrip model (execute in the folder you want to save to and don't forget to save the brain masks!).

'''
export FREESURFER_HOME=~/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh

mri_synthstrip -i /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/sub-tle001_ses-preop_acq-mpmFlash1_rec-loraks_echo-1_flip-6_mt-off_part-mag_MPM_RFSC_PD.nii -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/sub-tle001_ses-preop_acq-mpmFlash1_rec-loraks_echo-1_flip-6_mt-off_part-mag_MPM_RFSC_PD_strip.nii -m /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/brain_mask_sub-tle001.nii
mri_synthstrip -i /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/sub-tle002_ses-preop_acq-mpmFlash1_rec-loraks_echo-1_flip-6_mt-off_part-mag_MPM_RFSC_PD.nii -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/sub-tle002_ses-preop_acq-mpmFlash1_rec-loraks_echo-1_flip-6_mt-off_part-mag_MPM_RFSC_PD_strip.nii -m /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/brain_mask_sub-tle002.nii
mri_synthstrip -i /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/sub-tle005_ses-preop_acq-flash_run-1_echo-1_flip-6_mt-off_part-mag_MPM_RFSC_PD.nii -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/sub-tle005_ses-preop_acq-flash_run-1_echo-1_flip-6_mt-off_part-mag_MPM_RFSC_PD_strip.nii -m /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/brain_mask_sub-tle005.nii
'''

Now delete the old PD files in the folder and remove the _strip part of the file names. Sadly you cannot overwrite the file with these commands for some reason.
Next, we convert the R maps to the T maps and process these (the script may need a couple seconds).

'''
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/T-map_preparation.py -r /projects/crunchie/Jan/Daten/DataSubjects/sub-tle -f /ses-preop/anat/Results  -m /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/brain_mask_sub-tle  -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -s 001 002 005
'''
 
In case you are interested in the individual processing step that this script executes follow the explanation in the next section. If you also want to process the R maps without converting them to T maps you can use the following command.

'''
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/R-map_processing.py -r /projects/crunchie/Jan/Daten/DataSubjects/sub-tle -f /ses-preop/anat/Results  -m /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/brain_mask_sub-tle  -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -s 001 002 005
'''

Finally, if we want to use the data for a JEMRIS phantom, we need to save all of the data to a `.mat` or `.h5` file. Note that both slices and whole volumes can be generated (the latter might be too large for JEMRIS). These files can also be used for the Koma tool, however a seperate (julia) script needs to be used to convert into their own `.phantom` format.
 
'''
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/volume_to_mat.py -i /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -o /home/janmeyer/Pre-Processing/volumes_mat/TLE/MPM/Mat -s 001 002 005
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/slices_to_mat.py -i /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -o /home/janmeyer/Pre-Processing/slices/TLE/MPM/Mat -s 001 002 005

/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/slices_to_h5.py -i /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -o /home/janmeyer/Pre-Processing/slices/TLE/MPM/h5 -R True -s 001 002 005

/home/janmeyer/.julia/juliaup/julia-1.11.4+0.x64.linux.gnu/bin/julia volume_mat_to_phantom.jl -i /home/janmeyer/Pre-Processing/volumes_mat/TLE/MPM/sub-tle -o /home/janmeyer/KomaSimulations/Phantoms/Volumes/sub-tle -s 001,002,005
/home/janmeyer/.julia/juliaup/julia-1.11.4+0.x64.linux.gnu/bin/julia slices_mat_to_phantom.jl -i /home/janmeyer/Pre-Processing/slices_mat/TLE/MPM -o /home/janmeyer/KomaSimulations/Phantoms/Slices/sub-tle -s 001,002,005
'''

## Example workflow for the travelling head study data
We start again with data cleaning of the PD maps.

'''
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/data_cleaner.py -r /projects/crunchie/Jan/Daten/DataTravellingHeadStudy/sub-phy -f /ses-001/anat -o /home/janmeyer/Pre-Processing/processed_data/TravellingHead/MPM -c PD -s 001 002 003 004
'''

We then need to skull-strip one PD volume per subject and afterwards rename these volumes and apply the mask to the other PD maps.

'''
export FREESURFER_HOME=~/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh

mri_synthstrip -i /home/janmeyer/Pre-Processing/processed_data/TravellingHead/MPM/sub-phy001_ses-001_acq-dznbneep3dUP_run-1_echo-1_part-phase_PDw.nii.gz -o /home/janmeyer/Pre-Processing/processed_data/TravellingHead/MPM/sub-phy001_ses-001_acq-dznbneep3dUP_run-1_echo-1_part-phase_PDw_strip.nii.gz -m /home/janmeyer/Pre-Processing/processed_data/TravellingHead/MPM/brain_mask_sub-phy001.nii

/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/apply_mask.py -i /home/janmeyer/Pre-Processing/processed_data/TravellingHead/MPM -m /home/janmeyer/Pre-Processing/processed_data/TravellingHead/MPM/brain_mask_sub-phy -c PD -s 001 002 005
'''


## Individual steps for T map preparation
In case you don't want to do every processing step the processing of `T-map_preparation.py` can be broken down using multiple different scripts. Follow the guide and enter the commands instead of the third step in the example workflow (Note that this will generate a lot of different files with all the different processing steps).
First, as we need T weighted instead of R weighted images for the JEMRIS phantom, we threshold and invert the R maps to get our T maps.

'''
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/R_to_T_converter.py -i /projects/crunchie/Jan/Daten/DataSubjects/sub-tle001/ses-preop/anat/Results -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM 
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/R_to_T_converter.py -i /projects/crunchie/Jan/Daten/DataSubjects/sub-tle002/ses-preop/anat/Results -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM 
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/R_to_T_converter.py -i /projects/crunchie/Jan/Daten/DataSubjects/sub-tle005/ses-preop/anat/Results -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM 
'''

We than clean the T maps, same as with the PD maps before.

'''
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/data_cleaner.py -i /home/janmeyer/Pre-Processing/processed_data/MPM -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -c T1
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/data_cleaner.py -i /home/janmeyer/Pre-Processing/processed_data/MPM -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -c T2s
'''

Now we can use the brain masks of each patient to skull-strip the T maps (Note that cleaning and skull-stripping can be swapped, but both need to be done before removing the outliers to ensure good results).

'''
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/apply_mask.py -i /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -m /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/brain_mask_sub-tle -c T1 T2s -s 001 002 005 
'''

Next, to improve the signal quality and image contrast, we remove some of the outliers from T1 and T2star (just T2s in the file name). The chosen contrast can be set with the `-c`argument (only one at a time). Note that R2star generally has more outlier and we thus choose a slightly lower threshold, however these values might vary so you might need to tweak it for the data to look good.

'''
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/outlier_remover.py -i /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -c T1 -t 0.999
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/outlier_remover.py -i /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -c T2s -t 0.99
'''

