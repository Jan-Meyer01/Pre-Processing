% This is the main script for generating the MPM from the travelling head
% study data

% It contains the following steps:
%  1. Define the path to the data
%  2. Loop over all the sites and subjects
%  3. See which site and use the corresponding batch and write down the
%     problems I had with each site

% Skript to estimate the MPM datasets with the hMRI toolbox
analysisParameters = struct;
[analysisParameters.codeDir, ~, ~] = fileparts(mfilename('fullpath'));
analysisParameters.useDenoising = false; %true; % denoising causes a strange error to occur

%%  0. Setup


% change the following five paths to your setup!!
% add path to SPM installation (for linux https://en.wikibooks.org/wiki/SPM/Installation_on_64bit_Linux)
addpath('/home/janmeyer/spm12-r7771/');
% add path to hMRI installation (https://github.com/hMRI-group/hMRI-toolbox/wiki/GetStarted)
addpath('/home/janmeyer/hMRI-toolbox-0.6.1/');
% add path to BIDS installation (clone https://github.com/bids-standard/bids-matlab)
addpath(genpath('/home/janmeyer/opt/bids-matlab/')); 
% select main directory where the raw data is stored
mainPath = '/projects/crunchie/Jan/Daten/DataTravellingHeadStudy/Raw/';
% select output directory for storing the MPMs (automatically adds subfolder MPMs)
p_out = '/projects/crunchie/Jan/Daten/DataTravellingHeadStudy';


% add path to helper functions
addpath(genpath([analysisParameters.codeDir filesep 'functions']));

%mainPath = uigetdir('', 'Select the folder in which all datasets of the sites are stored');
if isempty(mainPath)
    error('No folder selected. Script will be aborted.');
end
fprintf('\n')
disp(['Selected data path: ', mainPath]);
fprintf('\n')
disp('It is assumed that this path contains all data for which MPM maps are to be estimated. The data should be stored in the order 1. site, 2. subject')

% Find the subjects and sites located in the main directory
siteList = dir(mainPath);           % List all subdirectories (Sites)
siteList = siteList([siteList.isdir] & ~startsWith({siteList.name}, '.'));  
dataStruct = struct();              % Structure to store Sites and Subjects

% Find all the subjects
for i = 1:length(siteList)          % Iterate through all Sites
    siteName = siteList(i).name;
    sitePath = fullfile(mainPath, siteName);
    
    % Find subject folders containing the subject-ID
    subjects = findSubjectFolders(sitePath, 'sub-phy'); 
    dataStruct.(siteName) = subjects; % Store in structure
end

% Display detected Sites and Subjects
fprintf('\n')
disp('Detected Sites and Subjects:');
disp(dataStruct);
fprintf('\n')

%p_out = uigetdir('', 'Select Output Directory');
if isempty(p_out)
    disp('No folder was selected, therefore the files will be safed in following structure: /mainPath/MPMs/site/subject')
    p_out = mainPath;
else
    disp(['Selected Output Directory: ', p_out]);
    disp(['The files will be safed in following structure: /' p_out '/MPMs/site/subject'])
end
warning('so far the Outputstructure for the MPM datasets is not ideal, e.g. for phy001:(SPM saves the first run as Results and the second run as Run_02). Maybe you have to change that for later')

%% 1. Estimate the MPM maps 
siteNames = fieldnames(dataStruct);

% Loop over the sites:
for indx_site = 1:length(siteNames) %[5,2]
    % get site name and subjects that were scanned there
    site = siteNames{indx_site};
    subjects = dataStruct.(site);
    
    % !!Noch etwas fishy geschrieben, da die Worte der Site im Namen
    % enthalten sein müssen
    
    if contains(site, 'Kings')
       MPM_preprocessed_B1(analysisParameters, site, subjects, mainPath, p_out, 80) %Flip angle was 80°
    elseif contains(site, 'Bonn')
       MPM_preprocessed_B1(analysisParameters, site, subjects, mainPath, p_out, 60)
    elseif contains(site, 'UCL')
        % Because of a problem with missing header information the b0 and b1 maps have to be converted with spm (DICOM to NIFTI)
        %DICOM_Import(analysisParameters, site, subjects, mainPath, p_out, '0*mfc_seste_b1map_v1h', '0*gre_field_mapping*', '*.ima')
        % Now the official MPM estimation can be performed:
        MPM_3DEPI_B1_correction(analysisParameters, site, subjects, mainPath, p_out, '*mfc_seste_b1map_v1h*', '*gre_field_mapping*')
    elseif contains(site, 'Leipzig')
        % Because of a problem with missing header information the b0 and b1 maps have to be converted with spm (DICOM to NIFTI)
        %DICOM_Import(analysisParameters, site, subjects, mainPath, p_out, '*al_B1mapping_v2d*', '*gre_field_mapping*', '*.dcm')
        % Now the official MPM estimation can be performed:   
        MPM_3DEPI_B1_correction(analysisParameters, site, subjects, mainPath, p_out, '*al_B1mapping_v2d*', '*gre_field_mapping*')
    elseif contains(site, 'Hamburg')
        warning('In my case both protocols are saved in the same folder: ses-001 = 3D EPI and ses-002 = FLASH, maybe you store them differently, then you have to adapt the code')
        % ses-001: 3D EPI - process the same as Bonn
        % ses-002: FLASH  - process like UCL/Leipzig

        % Way to go for the 3D EPI data (ses-001):
        MPM_preprocessed_B1(analysisParameters, site, subjects, mainPath, p_out, 0)

        % Way to go for the FLASH data (ses-002):
        % Because of a problem with missing header information the b0 and b1 maps have to be converted with spm (DICOM to NIFTI)
        %DICOM_Import(analysisParameters, site, subjects, mainPath, p_out, '0*mfc_seste_b1map_v1e*', '0*gre_field_map*', '*MR*')
        % Now the official MPM estimation can be performed:   
        MPM_3DEPI_B1_correction(analysisParameters, site, subjects, mainPath, p_out, '*mfc_seste_b1map_v1e*', '*gre_field_map*')
    else
       fprintf('\n')
       disp(['For site ', site, ' exists no hMRI Batch']);
    end
end

fprintf('\n')
disp('Everything worked out smoothly. Only one last warning:');
fprintf('\n')
warning('The PD maps are often different by a factor of 1000 (e.g. 3D EPI & FLASH from Hamburg. That is because of the different units from TR in the headers(s or ms). Decide which unit you want and multiply the factor accordinlgy')
