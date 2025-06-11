function MPM_3DEPI_B1_correction(analysisParameters, site, subjects, mainPath, p_out, naming_b1, naming_b0)

%% additional stuff for orientation
% 1 )Für dein Verständnis: Das 3D EPI bezieht sich hier NICHT auf das MPM Protokol (zur Erinnerung
%wir vergleichen FLASH mit 3D EPI) sondern auf die Acquirierung der B1 map

% 2)That is the errormessage that appears if you directly use the BIDS
warning('For UCL & Hamburg the BIDS converted data (b1 and b0 maps) are not enough, therefore the NIFTI converted data from SPM are used')

%converted data:

%------------ B1 MAP CALCULATION (i3D_EPI) 18-Feb-2025 16:39:51 ------------

% SE/STE EPI protocol selected ...
%
% WARNING: No B1mapNominalFAValues found in the extended header
%
% WARNING: using defaults value for nominal SE/STE flip angle values
% (115 110 105 100 95 90 85 80 75 70 65 ) instead of metadata
% 18-Feb-2025 16:39:51 - Failed  'Create hMRI maps'
% Error using hmri_create_b1map>get_b1map_params (line 807)
% Number of B1 mapping image pairs (121) does not match the number of nominal flip angles (11)!
% In file "/home/luethi/Powerslave/ISMRM2025/hMRI-toolbox-master/hmri_create_b1map.m" (???), function "get_b1map_params" at line 807.
% In file "/home/luethi/Powerslave/ISMRM2025/hMRI-toolbox-master/hmri_create_b1map.m" (???), function "hmri_create_b1map" at line 27.
% In file "/home/luethi/Powerslave/ISMRM2025/hMRI-toolbox-master/hmri_run_create.m" (???), function "hmri_create_local" at line 122.
% In file "/home/luethi/Powerslave/ISMRM2025/hMRI-toolbox-master/hmri_run_create.m" (???), function "hmri_run_create" at line 21.
%
% The following modules did not run:
% Failed: Create hMRI maps

%% Here the official code starts
fprintf('\n');
disp(['--- Now processing: ', site])

data_path = fullfile(mainPath, site);
BIDS = bids.layout(data_path);

for inx_sub = 1:length(subjects) %%Loop over the subjects
    ses = {};
    ses_label = {};

    subj = subjects{inx_sub};    %Define subject
    sessions = bids.query(BIDS, 'sessions','sub', subj); %ses reflects if its scan-rescan for phy001

    if contains(site, 'Hamburg') %Should be 1 for UCL and Leipzig, but ses-002 for Hamburg
        ses_label = 'ses-002';
        ses = sessions{2};
    else
        ses_label = 'ses-001';
        ses = sessions{1};
    end

    %Define how many runs will be done:
    if strcmp(subj, 'sub-phy001')
        nrun= 2; %For scan, rescan
    else
        nrun = 1; %Only one scan exists
    end

    if analysisParameters.useDenoising
        jobfile = {fullfile(analysisParameters.codeDir, '/JobFiles/hMRImaps_with_3DEPI_B1_and_denoising_job.m')};
        jobs = repmat(jobfile, 1, nrun);
        inputs = cell(7, nrun);
    else
        jobfile = {fullfile(analysisParameters.codeDir,'/JobFiles/hMRImaps_with_3DEPI_B1_job.m')};
        jobs = repmat(jobfile, 1, nrun);
        inputs = cell(6, nrun);
    end




    for indx_run = 1:nrun
        % Create Output directory
        outputDir = fullfile(p_out , 'MPMs', site,  subj, ses_label, sprintf('run-%03d', indx_run));

        if ~exist(outputDir, 'dir')
            mkdir(outputDir);
        else
            disp(['folder already exists ', outputDir]);
        end

        if analysisParameters.useDenoising

            outputDirDenoising = [outputDir filesep 'denoisedMPMs'];
            if ~exist(outputDirDenoising, 'dir')
                mkdir(outputDirDenoising);
            else
                disp(['folder already exists ', outputDirDenoising]);
            end

        end

        % B1 map
        NiftiDir = fullfile(p_out , 'MPMs', site,  subj, ses_label, sprintf('run-%03d', indx_run));
        %Wenn man jetzt richtig krass ist, dass könnte man die outputDir von der "DICOM_Import" function global
        %definieren und hier einfach verwenden. Leider kann ich das noch nicht :/

        b1_images=dir(fullfile(NiftiDir, naming_b1, '*.nii'));
        b1_filePaths={};
        b0_images=dir(fullfile(NiftiDir, naming_b0, '*.nii'));
        b0_filePaths={};


        %Iterate over the data
        for i=1:length(b1_images)
            if ~b1_images(i).isdir
                b1_fullFilePath=fullfile(b1_images(i).folder, b1_images(i).name);
                b1_filePaths{i,1} = b1_fullFilePath;
            end
        end

        %Iterate over the data
        for i=1:length(b0_images)
            if ~b0_images(i).isdir
                b0_fullFilePath=fullfile(b0_images(i).folder, b0_images(i).name);
                b0_filePaths{i,1} = b0_fullFilePath; %Add all DICOM images in the same array
            end
        end

        if contains(site, 'UCL')
            if strcmp(subj, 'sub-phy001') %Because for phy001 the runs are included in the BIDS Structure
                PDw_maps_path=bids.query(BIDS, 'data', 'ses', ses, 'sub', subj, 'modality', 'anat', 'extension', '.nii', 'part', 'mag', 'mt', 'off', 'flip', '6', 'run', string(indx_run), 'acq', 'pdwmfc3dflashv3kR464ch', 'rec', 'ND');
                T1w_maps_path=bids.query(BIDS, 'data', 'ses', ses, 'sub', subj, 'modality', 'anat', 'extension', '.nii', 'part', 'mag','mt', 'off', 'flip', '21', 'run', string(indx_run), 'acq', 't1wmfc3dflashv3kR464ch', 'rec', 'ND');
                MTw_maps_path=bids.query(BIDS, 'data', 'ses', ses, 'sub', subj, 'modality', 'anat', 'extension', '.nii', 'part', 'mag', 'mt', 'off', 'flip', '6', 'run', string(indx_run), 'acq', 'mtwmfc3dflashv3kR464ch', 'rec', 'ND');%Achtung, checken warum BIDS hier scheiße baut
            else %Because for phy002-004 the runs are NOT included in the BIDS Structure
                PDw_maps_path=bids.query(BIDS, 'data', 'ses', ses, 'sub', subj, 'modality', 'anat', 'extension', '.nii', 'part', 'mag', 'mt', 'off', 'flip', '6', 'acq', 'pdwmfc3dflashv3kR464ch', 'rec', 'ND');
                T1w_maps_path=bids.query(BIDS, 'data', 'ses', ses, 'sub', subj, 'modality', 'anat', 'extension', '.nii', 'part', 'mag', 'mt', 'off', 'flip', '21', 'acq', 't1wmfc3dflashv3kR464ch', 'rec', 'ND');
                MTw_maps_path=bids.query(BIDS, 'data', 'ses', ses, 'sub', subj, 'modality', 'anat', 'extension', '.nii', 'part', 'mag', 'mt', 'off', 'flip', '6', 'acq', 'mtwmfc3dflashv3kR464ch', 'rec', 'ND'); %Achtung, checken warum BIDS hier scheiße baut
            end

        elseif contains(site, 'Leipzig')
            if strcmp(subj, 'sub-phy001') %Because for phy001 the runs are included in the BIDS Structure
                PDw_maps_path=bids.query(BIDS, 'data', 'ses', ses, 'sub', subj, 'modality', 'anat', 'extension', '.nii', 'part', 'mag', 'mt', 'off', 'flip', '6', 'run', string(indx_run));
                T1w_maps_path=bids.query(BIDS, 'data', 'ses', ses, 'sub', subj, 'modality', 'anat', 'extension', '.nii', 'mt', 'off', 'flip', '21', 'run', string(indx_run));
                MTw_maps_path=bids.query(BIDS, 'data', 'ses', ses, 'sub', subj, 'modality', 'anat', 'extension', '.nii', 'part', 'mag', 'mt', 'on', 'flip', '6', 'run', string(indx_run));
            else %Because for phy002-004 the runs are NOT included in the BIDS Structure
                PDw_maps_path=bids.query(BIDS, 'data', 'ses', ses, 'sub', subj, 'modality', 'anat', 'extension', '.nii', 'part', 'mag', 'mt', 'off', 'flip', '6');
                T1w_maps_path=bids.query(BIDS, 'data', 'ses', ses, 'sub', subj, 'modality', 'anat', 'extension', '.nii', 'mt', 'off', 'flip', '21');
                MTw_maps_path=bids.query(BIDS, 'data', 'ses', ses, 'sub', subj, 'modality', 'anat', 'extension', '.nii', 'part', 'mag', 'mt', 'on', 'flip', '6');
            end
        elseif contains(site, 'Hamburg')
            if strcmp(subj, 'sub-phy001') %Because for phy001 the runs are included in the BIDS Structure
                PDw_maps_path=bids.query(BIDS, 'data', 'ses', ses, 'sub', subj, 'modality', 'anat', 'extension', '.nii', 'part', 'mag', 'mt', 'off', 'flip', '6', 'run', string(indx_run));
                T1w_maps_path=bids.query(BIDS, 'data', 'ses', ses, 'sub', subj, 'modality', 'anat', 'extension', '.nii', 'part', 'mag','mt', 'off', 'flip', '21', 'run', string(indx_run));
                MTw_maps_path=bids.query(BIDS, 'data', 'ses', ses, 'sub', subj, 'modality', 'anat', 'extension', '.nii', 'part', 'mag', 'mt', 'on', 'flip', '6', 'run', string(indx_run));%Achtung, checken warum BIDS hier scheiße baut
            else %Because for phy002-004 the runs are NOT included in the BIDS Structure
                PDw_maps_path=bids.query(BIDS, 'data', 'ses', ses, 'sub', subj, 'modality', 'anat', 'extension', '.nii', 'part', 'mag', 'mt', 'off', 'flip', '6', 'rec', 'ND');
                T1w_maps_path=bids.query(BIDS, 'data', 'ses', ses, 'sub', subj, 'modality', 'anat', 'extension', '.nii', 'part', 'mag', 'mt', 'off', 'flip', '21',  'rec', 'ND');
                MTw_maps_path=bids.query(BIDS, 'data', 'ses', ses, 'sub', subj, 'modality', 'anat', 'extension', '.nii', 'part', 'mag', 'mt', 'on', 'flip', '6', 'rec', 'ND'); %Achtung, checken warum BIDS hier scheiße baut
            end

        end

        for i=1:length(PDw_maps_path)
            PDw_maps_path{i}=[PDw_maps_path{i}, ',1'];
            T1w_maps_path{i}=[T1w_maps_path{i}, ',1'];
        end
        for i=1:length(MTw_maps_path)
            MTw_maps_path{i}=[MTw_maps_path{i}, ',1'];
        end


        if analysisParameters.useDenoising

            inputs{1, indx_run} = {outputDirDenoising}; % Denoising: Output directory - cfg_files
            if contains(site, 'UCL') %There are too many echoes for PD and T1 and they definitely have to be sorted in ascending order. That is why the sorting here
                %Clearly 10 is alphabetically before 1 :)
                inputs{2, indx_run} = PDw_maps_path([2:end, 1]); % Create hMRI maps: PD images - cfg_files
                inputs{3, indx_run} = T1w_maps_path([2:end, 1]); % Create hMRI maps: T1 images - cfg_files
            else
                inputs{2, indx_run} = PDw_maps_path; % Create hMRI maps: PD images - cfg_files
                inputs{3, indx_run} = T1w_maps_path; % Create hMRI maps: T1 images - cfg_files
            end


            inputs{4, indx_run} = MTw_maps_path; % Denoising: Magnitude input - cfg_files
            inputs{5, indx_run} = {outputDir}; % Create hMRI maps: Output directory - cfg_files
            inputs{6, indx_run} = b1_filePaths; % Create hMRI maps: B1 input - cfg_files
            inputs{7, indx_run} = b0_filePaths; % Create hMRI maps: B0 input - cfg_files


        else


            inputs{1, indx_run} = {outputDir}; % Create hMRI maps: Output directory - cfg_files
            inputs{2, indx_run} = b1_filePaths'; % Create hMRI maps: B1 input - cfg_files
            inputs{3, indx_run} = b0_filePaths'; % Create hMRI maps: B0 input - cfg_files
            inputs{4, indx_run} = MTw_maps_path; % Create hMRI maps: MT images - cfg_files

            if contains(site, 'UCL') %There are too many echoes for PD and T1 and they definitely have to be sorted in ascending order. That is why the sorting here
                %Clearly 10 is alphabetically before 1 :)
                inputs{5, indx_run} = PDw_maps_path([2:end, 1]); % Create hMRI maps: PD images - cfg_files
                inputs{6, indx_run} = T1w_maps_path([2:end, 1]); % Create hMRI maps: T1 images - cfg_files
            else
                inputs{5, indx_run} = PDw_maps_path; % Create hMRI maps: PD images - cfg_files
                inputs{6, indx_run} = T1w_maps_path; % Create hMRI maps: T1 images - cfg_files
            end



        end



    end
    spm('defaults', 'FMRI');
    spm_jobman('run', jobs, inputs{:});
end
