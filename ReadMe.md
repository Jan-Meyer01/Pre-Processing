# Processing of the Travelling Heads Dataset and Conversion into Phantoms for the JEMRIS and KomaMRI Simulation Tools

## Introduction
The Travelling Heads Study collected data at multiple sites with the goal of comparing and validating the use of quantitative MR imaging methods. In this repository, we use the MPM and T2 data obtained in this study in order to get realistic phantoms to be used for simulations of different imaging sequences.

## Data Structure
In the study, data from four volunteers was collected.
Brain masks and segmentation also need to exist for the main pipeline.

## Main Workflow
There are two important scripts that can be executed via the command line to get the processed data and convert it into phantoms for either [JEMRIS](https://www.jemris.org/) or [KomaMRI](https://juliahealth.org/KomaMRI.jl/stable/).

First, we process the data which will create the correct folder structure in this repo:

```
python process_travelling_head_data.py -i <path_to_data>
```

There are a couple of options, which you'll need, or might want to, specify:
- `input_dir`: The input directory where your MPM and T2 data is stored (this needs to be specified!!)
- `output_dir`: The output directory where the processed data will be stord. By default this will be done in this repository and we advise against changing it as the other scripts expect this directory and additional input would be required there to use a different directory.

Once you have your processed data in your desired folder, you can start to generate the phantoms. There are two scripts for this, one for JEMRIS and one for KomaMRI. Lets start with JEMRIS:

```
python phantoms_jemris_travelling_head.py
```

This script has a lot of options, which are quickly explained here:
- `input_dir`: The input directory where your processed data from the previous step is stored (needs to be the same as the `output_dir` before!!). By default, the phantoms will automatically be stored in a subfolder called `phantoms_jemris` for each subject in this repo or your specified path.
- `format`: The file format in which the phantoms will be stored. It can either be `.mat` or `.h5`. 
- `Rmaps`: Directly tied to the format is the option of using either T or R maps. If nothing is specified, T maps will be used, however putting `--Rmaps` or `-R` will generate R maps and this option should only be used with the `.h5` file format to be called with a script as the JEMRIS GUI will complain otherwise.
- `slice`: If nothing is specified the whole volume will be used as a phantoms, but using `--slice` or `-s` slices along the z-axis can be generated in case the volumes are to computationally draining.

The script for KomaMRI is very similar, but is written in [Julia](https://julialang.org), same as the simulatinos tool. The command does not look quite as nice, but it can also be run in the terminal, you just need to change the path to Julia:

```
/home/janmeyer/.julia/juliaup/julia-1.11.4+0.x64.linux.gnu/bin/julia phantoms_koma_travelling_head.jl  
```

The options are similar to before:
- `i`: The input directory where your processed data from the previous step is stored (needs to be the same as the `output_dir` before!!). By default, the phantoms will automatically be stored in a subfolder called `phantoms_koma` for each subject in this repo or your specified path.
- `s`: If nothing is specified the whole volume will be used as a phantoms, but using the `-s` flag slices along the z-axis can be generated in case the volumes are to computationally draining.

Since there are no other file types for KomaMRI, the other option are not needed. If you followed the commands above using the default values the resulting data structure should look something like this for the first subject (others are the same):


```
├── processed_data
│   └── travelling_head
│       └── site (e.g. Bonn_Skyra_3T_LowRes)
│           ├── sub-phy001
│               ├── ses-001
│                   ├──phantoms_jemris
│                       └── sub-phy001_ses-001_phantom_jemris_vol.mat
│                   ├──phantoms_koma
│                       └── sub-phy001_ses-001_phantom_koma_vol.phantom
│                   ├── R2_map.nii.gz
│                   ├── sub-phy001_ses-001_acq-dznebnep3dPDw_run-1_echo-1_flip-4_mt-off_part-mag_MPM_PD.nii
│                   ├── sub-phy001_ses-001_acq-dznebnep3dPDw_run-1_echo-1_flip-4_mt-off_part-mag_MPM_R1.nii
│                   ├── sub-phy001_ses-001_acq-dznebnep3dPDw_run-1_echo-1_flip-4_mt-off_part-mag_MPM_R2_OLS.nii
│                   ├── sub-phy001_ses-001_acq-dznebnep3dPDw_run-1_echo-1_flip-4_mt-off_part-mag_MPM_T1.nii
│                   ├── sub-phy001_ses-001_acq-dznebnep3dPDw_run-1_echo-1_flip-4_mt-off_part-mag_MPM_T2_OLS.nii
│                   └── T2_map.nii.gz
```

While the naming for the phantoms is somewhat redundant with the folder names, it seemed usefull for a quicker overview as the file names are not that long. The structure for the other subjects is the same and only subject 1 has two sessions due to the rescan.









# Old Documentation (Remove unnecessary parts later)


## Pre-processing for nifti 3D images
This repository contains code for different processing steps of the TLE and Travelling Head datasets. Please ensure that your data is in the [NIfTI](https://nifti.nimh.nih.gov/) data format as all scripts require this. In case you need to convert your data from DICOM first use e.g. [Bidscoin](https://bidscoin.readthedocs.io/en/latest/) for a [Bids](https://bids-specification.readthedocs.io/en/stable/) conform conversion.

The folder structure used by the command/scripts in the sections below is as follows.

```
├── processed_data
│   └── TLE
│       └── MPM
│           ├── brain_mask_sub-tle001.nii
│           ├── sub-tle001_ses-preop_acq-mpmFlash1_rec-loraks_echo-1_flip-6_mt-off_part-mag_MPM_RFSC_PD.nii
│           ├── sub-tle001_ses-preop_acq-mpmFlash1_rec-loraks_echo-1_flip-6_mt-off_part-mag_MPM_RFSC_R1.nii
│           ├── sub-tle001_ses-preop_acq-mpmFlash1_rec-loraks_echo-1_flip-6_mt-off_part-mag_MPM_RFSC_R2s_OLS.nii
│           ├── sub-tle001_ses-preop_acq-mpmFlash1_rec-loraks_echo-1_flip-6_mt-off_part-mag_MPM_RFSC_T1.nii
│           └── sub-tle001_ses-preop_acq-mpmFlash1_rec-loraks_echo-1_flip-6_mt-off_part-mag_MPM_RFSC_T2s_OLS.nii
│   └── travelling_head
│       └── site
│           ├── sub-phy001
│               ├── ses-001
│                   ├── sub-phy001_ses-001_acq-dznebnep3dPDw_run-1_echo-1_flip-4_mt-off_part-mag_MPM_PD.nii
│                   ├── sub-phy001_ses-001_acq-dznebnep3dPDw_run-1_echo-1_flip-4_mt-off_part-mag_MPM_R1.nii
│                   ├── sub-phy001_ses-001_acq-dznebnep3dPDw_run-1_echo-1_flip-4_mt-off_part-mag_MPM_R2_OLS.nii
│                   ├── sub-phy001_ses-001_acq-dznebnep3dPDw_run-1_echo-1_flip-4_mt-off_part-mag_MPM_T1.nii
│                   ├── sub-phy001_ses-001_acq-dznebnep3dPDw_run-1_echo-1_flip-4_mt-off_part-mag_MPM_T2_OLS.nii
│                   └── sub-phy001_ses-001_acq-dznebnep3dPDw_run-1_echo-1_flip-4_mt-off_part-mag_MPM_T2.nii
│
├── slices
│   └── TLE
│       └── MPM
│           └── h5
│               ├── sub-tle001_MPM_Slice001.h5
│               ├── sub-tle001_MPM_Slice002.h5
│               ├── sub-tle001_MPM_Slice003.h5
│               ├── ...
│               └── sub-tle001_MPM_Slice176.h5
│           └── Mat
│               ├── sub-tle001_MPM_Slice001.mat
│               ├── sub-tle001_MPM_Slice002.mat
│               ├── sub-tle001_MPM_Slice003.mat
│               ├── ...
│               └── sub-tle001_MPM_Slice176.mat
│   └── TravellingHead
│       └── MPM
│           └── h5
│               ├── sub-phy001_MPM_Slice001.h5
│               ├── sub-phy001_MPM_Slice002.h5
│               ├── sub-phy001_MPM_Slice003.h5
│               ├── ...
│               └── sub-phy001_MPM_Slice176.h5
│           └── Mat
│               ├── sub-phy001_MPM_Slice001.mat
│               ├── sub-phy001_MPM_Slice002.mat
│               ├── sub-phy001_MPM_Slice003.mat
│               ├── ...
│               └── sub-phy001_MPM_Slice176.mat
```

Note that only example files for one subject are given and that besides MPMs you might want to process DWI or other data for a dataset, which is why the highest level is the dataset. Otherwise you might also want to add a `volumes` folder if you do not want to use slices as phantoms. You can choose your own structure easialy by changing the paths. There are no paths hidden in the scripts and everything about the folder structure should be configurable by the terminal commands. If you do not plan to create phantoms for simulations you could also skip the `processed_data` folder and move the dataset specific folder one level up. On an additional note, the folders here are not provided in the repositiory as it would make it very unreadable, however all these will be generated automatically if you run the scripts below. If you want this folder structure, which we would encurage for the code to work properly, you only need to change the external path to your unprocessed data. This path should look something like `/home/user/data/sub-tle00x/ses-preop/anat/Results` for the TLE or `/home/user/data/sub-phy00x/ses-001/anat/` for the travelling head data where you need to replace `/home/user/data/` with your systems path to the dataset containing the subjects. Thus, you just need to replace `/projects/crunchie/Jan/Daten/DataSubjects/sub-tle` and `/ses-preop/anat/Results` below and once the data is stored in the `processed_data` folder you can use the other commands without changing anything - don't worry, the original unprocessed data won't be changed just the data under `processed_data`!

# One liner
TODO: All pre-processing steps can be done with a single command. This includes data cleaning, converting the R maps to T maps, skull-stripping and outlier removal. In order to use the processed images as phantoms for [JEMRIS](https://www.jemris.org/) the individual volumes need to be `.mat` files. Thus, after executing the command, you will have the PD, T1, T2star files in the specified folder and the `.mat` files of all the data for a single subject in potentially another. Please adjust the file paths to your system.

```
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/pre-processing.py -r /projects/crunchie/Jan/Daten/DataSubjects/sub-tle -f /ses-preop/anat/Results -oF /home/janmeyer/Pre-Processing/processed_data/TLE/MPM  -oM /home/janmeyer/Pre-Processing/volumes/TLE/MPM/Mat -s 001 002 005
```

If you want to convert the phantoms from `.mat` to `.phantom` files to use in [Koma](https://juliahealth.org/KomaMRI.jl/stable/) you can do so with the following scripts (note that you need to have [Julia](https://julialang.org/downloads/) installed and change the first path to your installation). This can be either done for a whole volume (might take multiple GB of GPU to run later in the simulations) or just some slices.

```
/home/janmeyer/.julia/juliaup/julia-1.11.4+0.x64.linux.gnu/bin/julia volume_mat_to_phantom.jl -i /home/janmeyer/Pre-Processing/volumes/TLE/MPM/Mat/sub-tle -o /home/janmeyer/KomaSimulations/Phantoms/Volumes/sub-tle -s 001,002,005
/home/janmeyer/.julia/juliaup/julia-1.11.4+0.x64.linux.gnu/bin/julia slices_mat_to_phantom.jl -i /home/janmeyer/Pre-Processing/slices_mat/TLE/MPM/sub-tle -o /home/janmeyer/KomaSimulations/Phantoms/Slices/sub-tle -s 001,002,005
```

## Example workflow for the TLE data
To get a better understanding of the pipeline, we can examine its main four steps (each can be executed as an independent script). 
First, we clean up all the PD image (see different subject numbers). The first command is for the TLE data and the second for the travelling head study data.

```
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/data_cleaner.py -r /projects/crunchie/Jan/Daten/DataSubjects/sub-tle -f /ses-preop/anat/Results -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -c PD -s 001 002 005
```

Then we extract the brains of the PDs (as these usually have the best quality) without any further processing using freesurfers synthstrip model (execute in the folder you want to save to and don't forget to save the brain masks!).

```
export FREESURFER_HOME=~/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh

mri_synthstrip -i /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/sub-tle001_ses-preop_acq-mpmFlash1_rec-loraks_echo-1_flip-6_mt-off_part-mag_MPM_RFSC_PD.nii -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/sub-tle001_ses-preop_acq-mpmFlash1_rec-loraks_echo-1_flip-6_mt-off_part-mag_MPM_RFSC_PD_strip.nii -m /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/brain_mask_sub-tle001.nii
mri_synthstrip -i /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/sub-tle002_ses-preop_acq-mpmFlash1_rec-loraks_echo-1_flip-6_mt-off_part-mag_MPM_RFSC_PD.nii -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/sub-tle002_ses-preop_acq-mpmFlash1_rec-loraks_echo-1_flip-6_mt-off_part-mag_MPM_RFSC_PD_strip.nii -m /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/brain_mask_sub-tle002.nii
mri_synthstrip -i /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/sub-tle005_ses-preop_acq-flash_run-1_echo-1_flip-6_mt-off_part-mag_MPM_RFSC_PD.nii -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/sub-tle005_ses-preop_acq-flash_run-1_echo-1_flip-6_mt-off_part-mag_MPM_RFSC_PD_strip.nii -m /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/brain_mask_sub-tle005.nii
```

Now delete the old PD files in the folder and remove the _strip part of the file names. Sadly you cannot overwrite the file directrly with these commands for some reason.
Next, we convert the R maps to the T maps and process these (the script may need a couple seconds).

```
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/T-map_preparation.py -r /projects/crunchie/Jan/Daten/DataSubjects/sub-tle -f /ses-preop/anat/Results  -m /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/brain_mask_sub-tle  -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -s 001 002 005
```
 
In case you are interested in the individual processing step that this script executes follow the explanation in the next section. 
If you also want to process the R maps without converting them to T maps you can use the following command.

```
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/R-map_processing.py -r /projects/crunchie/Jan/Daten/DataSubjects/sub-tle -f /ses-preop/anat/Results  -m /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/brain_mask_sub-tle  -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -s 001 002 005
```

Finally, if we want to use the data for a JEMRIS phantom, we need to save all of the data to a `.mat` or `.h5` file. Note that both slices and whole volumes can be generated (the latter might be too large for JEMRIS). These files can also be used for the Koma tool, however a seperate (julia) script needs to be used to convert into their own `.phantom` format.
 
```
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/volume_to_mat.py -i /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -o /home/janmeyer/Pre-Processing/volumes/TLE/MPM/Mat -s 001 002 005
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/slices_to_mat.py -i /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -o /home/janmeyer/Pre-Processing/slices/TLE/MPM/Mat -s 001 002 005

/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/slices_to_h5.py -i /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -o /home/janmeyer/Pre-Processing/slices/TLE/MPM/h5 -R True -s 001 002 005

/home/janmeyer/.julia/juliaup/julia-1.11.4+0.x64.linux.gnu/bin/julia volume_mat_to_phantom.jl -i /home/janmeyer/Pre-Processing/volumes/TLE/MPM/Mat/sub-tle -o /home/janmeyer/KomaSimulations/Phantoms/Volumes/sub-tle -s 001,002,005
/home/janmeyer/.julia/juliaup/julia-1.11.4+0.x64.linux.gnu/bin/julia slices_mat_to_phantom.jl -i /home/janmeyer/Pre-Processing/slices/TLE/MPM/Mat -o /home/janmeyer/KomaSimulations/Phantoms/Slices/sub-tle -s 001,002,005
```

## Example workflow for the travelling head study data
For the travelling head study there is a joint script for the whole processing chain up to the phantom generation. You only need to input the path to your quantitative maps and the script will automatically generate the correct folder structure for you. Note that you need to be in the main folder for the following command to work properly.

```
python process_travelling_head_data.py -i /projects/crunchie/Jan/Daten/DataTravellingHeadStudy/MPMs
```


## Details on individual steps for T map preparation
In case you don't want to do every processing step the processing of `T-map_preparation.py` can be broken down using multiple different scripts. Follow the guide and enter the commands instead of the third step in the example workflow (Note that this will generate a lot of different files with all the different processing steps).
First, as we need T weighted instead of R weighted images for the JEMRIS phantom, we threshold and invert the R maps to get our T maps.

```
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/R_to_T_converter.py -i /projects/crunchie/Jan/Daten/DataSubjects/sub-tle001/ses-preop/anat/Results -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM 
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/R_to_T_converter.py -i /projects/crunchie/Jan/Daten/DataSubjects/sub-tle002/ses-preop/anat/Results -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM 
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/R_to_T_converter.py -i /projects/crunchie/Jan/Daten/DataSubjects/sub-tle005/ses-preop/anat/Results -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM 
```

We than clean the T maps, same as with the PD maps before.

```
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/data_cleaner.py -i /home/janmeyer/Pre-Processing/processed_data/MPM -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -c T1
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/data_cleaner.py -i /home/janmeyer/Pre-Processing/processed_data/MPM -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -c T2s
```

Now we can use the brain masks of each patient to skull-strip the T maps (Note that cleaning and skull-stripping can be swapped, but both need to be done before removing the outliers to ensure good results).

```
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/apply_mask.py -i /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -m /home/janmeyer/Pre-Processing/processed_data/TLE/MPM/brain_mask_sub-tle -c T1 T2s -s 001 002 005 
```

Next, to improve the signal quality and image contrast, we remove some of the outliers from T1 and T2star (just T2s in the file name). The chosen contrast can be set with the `-c`argument (only one at a time). Note that R2star generally has more outlier and we thus choose a slightly lower threshold, however these values might vary so you might need to tweak it for the data to look good.

```
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/outlier_remover.py -i /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -c T1 -t 0.999
/home/janmeyer/miniconda3/envs/PyTorch/bin/python /home/janmeyer/Pre-Processing/outlier_remover.py -i /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -o /home/janmeyer/Pre-Processing/processed_data/TLE/MPM -c T2s -t 0.99
```

