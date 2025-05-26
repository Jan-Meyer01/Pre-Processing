%Dicom Import for the b1 and b0 maps is currently needed for: The data from Leipzig, UCL and the
% Flash data from Hamburg (ses001)

function DICOM_Import(analysisParameters, site, subjects, mainPath, p_out, naming_b1, naming_b0, suffix)
 fprintf('\n');
    disp(['--- DICOMS from ', site, ' will be converted'])
    
    data_path = fullfile(mainPath, site);
    BIDS = bids.layout(data_path);

    for indx_sub = 1: length(subjects) %%Loop over the subjects
    
    %Define subject and session naming
        subj = subjects{indx_sub};
    
        sessions = bids.query(BIDS, 'sessions','sub', subj); %ses reflects if its scan-rescan for phy001
        
        if contains(site, 'Hamburg') %Should be 1 for UCL and Leipzig, but ses-002 for Hamburg
            ses_label = 'ses-002';
        else
            ses_label = 'ses-001';
        end
        
        %Define how many runs will be done:
        if strcmp(subj, 'sub-phy001')
            nrun= 2; %For scan, rescan
        else 
            nrun = 1; %Only one scan exists
        end
 
    jobfile = {fullfile(analysisParameters.codeDir, 'JobFiles/DICOM_Import_spm_job.m')};
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(2, nrun);

    for indx_run = 1:nrun
    % Create Output directory
    warning('The new NIFITs will be safed in the same folder as the MPMs, if you want to change this, modify it here')
        outputDir = fullfile(p_out , 'MPMs', site,  subj, ses_label, sprintf('run-%03d', indx_run));
    
        if ~exist(outputDir, 'dir')
            mkdir(outputDir);
        else
            disp(['folder already exists ', outputDir]);
        end

     %To get the B1 and B0 DICOM data
       warning('Add the path your NIFTI data are saved. Maybe needs to be modified. Currently I am expecting you to work on powerslave and the NIFTIS still being stored in the Messungen folder.')
        trunc_path = truncateToCrunchie(data_path); %Currently I am expecting you to work on powerslave and the NIFTIS are stored in the Messungen folder
        
         if strcmp(subj, 'sub-phy001') %Define the paths, because for phy001 there is an additional folde 'run-001' and for phy002-004 not
             Path_files=fullfile(filesep, trunc_path, 'Messungen', 'Travel_Head_Study', site, subj, ses_label,  sprintf('run-%03d', indx_run)); %Hinteren nehmen, da umpositioniert
        else 
            Path_files=fullfile(filesep, trunc_path, 'Messungen', 'Travel_Head_Study', site, subj, ses_label); %Hinteren nehmen, da umpositioniert
        end

        if contains(site, 'Hamburg') && strcmp(subj, 'sub-phy003') %For this specific dataset the b1 and b0 maps were acquired 2 tímes. Only the second version should be taken 
            b1_path=dir(fullfile(Path_files, naming_b1));
            b1_images=dir(fullfile(b1_path(2).folder, b1_path(2).name, suffix)); %Always take the second dataset
            b0_path=dir(fullfile(Path_files, naming_b0)); %Here are 4 folder-structures (We need 3 and 4)
            b0_selected_3=dir(fullfile(b0_path(3).folder, b0_path(3).name, suffix));
            b0_selected_4=dir(fullfile(b0_path(4).folder, b0_path(4).name, suffix));
        else 
            b1_images=dir(fullfile(Path_files, naming_b1, suffix));
            b0_images=dir(fullfile(Path_files, naming_b0, suffix));  
        end


        DICOM_filePaths={};

        %Iterate over the data
        for i=1:length(b1_images)
            if ~b1_images(i).isdir
                b1_fullFilePath=fullfile(b1_images(i).folder, b1_images(i).name);
                DICOM_filePaths{end+1} = b1_fullFilePath;
            end
        end
       
           %Iterate over the data
         if contains(site, 'Hamburg') && strcmp(subj, 'sub-phy003') %For this dataset its again special
            for i=1:length(b0_selected_3)
                if ~b0_selected_3(i).isdir
                    b0_fullFilePath=fullfile(b0_selected_3(i).folder, b0_selected_3(i).name);
                    DICOM_filePaths{end+1} = b0_fullFilePath; %Add all DICOM images in the same array
                end
            end 
            for i=1:length(b0_selected_4)
                if ~b0_selected_4(i).isdir
                    b0_fullFilePath=fullfile(b0_selected_4(i).folder, b0_selected_4(i).name);
                    DICOM_filePaths{end+1} = b0_fullFilePath; %Add all DICOM images in the same array
                end
           end   
         else
            for i=1:length(b0_images)
                if ~b0_images(i).isdir
                    b0_fullFilePath=fullfile(b0_images(i).folder, b0_images(i).name);
                    DICOM_filePaths{end+1} = b0_fullFilePath; %Add all DICOM images in the same array
                end
            end
         end


        inputs{1, indx_run} = DICOM_filePaths'; % DICOM Import: DICOM files - cfg_files
        inputs{2, indx_run} = {outputDir}; % DICOM Import: Output directory - cfg_files

    end
        spm('defaults', 'FMRI');
        spm_jobman('run', jobs, inputs{:});

        clear inputs
    end
  
