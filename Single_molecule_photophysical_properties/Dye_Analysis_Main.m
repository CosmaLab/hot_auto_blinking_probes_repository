%% Dye Analysis Main Script
% This script coordinates the analysis of multiple dye localization datasets
% It calls modular functions to perform localization, blinking analysis
%
% Required User Configuration:
% 1. Camera Parameters:
%    - params.camera.gain      : Camera gain (electrons/ADU)
%    - params.camera.QE        : Quantum efficiency
%    - params.camera.exposure_time : Exposure time (seconds)
%    - params.camera.frame_interval: Time between frames (seconds)
%    - params.EM_wave_FP       : Emission wavelength (nm)
%
% 2. Analysis Parameters:
%    - params.threshold        : Spot detection threshold
%    - params.mask_sz         : Mask size for localization
%    - params.threshold_distance: Min distance between spots
%    - params.halfsquare      : ROI size for blinking analysis
%    - params.signalChange    : Intensity change threshold
%
% 3. Path Configuration:
%    - params.totalfn_path    : Path to utility functions
%    - params.baseDir         : Base directory for data
close all
clear 
%% Check for Bio-Formats (needed for ND2 files)
if ~exist('bfGetReader', 'file')
    % Check if Bio-Formats is in the MATLAB path
    bfPath = fullfile(fileparts(mfilename('fullpath')), 'bfmatlab');
     
    % If Bio-Formats folder exists, add it to path
    if exist(bfPath, 'dir')
        addpath(bfPath);
        fprintf('Added Bio-Formats to MATLAB path: %s\n', bfPath);
    else
        warning(['Bio-Formats package not found. ND2 files will not be supported. ', ...
                'Download from https://www.openmicroscopy.org/bio-formats/downloads/ ', ...
                'and place in a folder named "bfmatlab" in the current directory.']);
    end
end

%% Global Parameters
% Analysis parameters
params = struct();
% parameters for localization
params.threshold = 5;           % Spot detection threshold
params.mask_sz = 3;            % Mask size for localization, iQ_pkfnd_V6 create a square matrix of ones with size sz×sz to find peak
params.dil_sz = 3;              % Dilation size
% params.bkg_percent = 10;        % Percentage of end frames for background
params.threshold_distance = 2.5; % Minimum distance to consider two spots as separate
params.frames_to_skip = 0;    % Number of initial frames to skip
% parameters for blinking analysis
params.halfsquare = 1;          % Half size of ROI around each point (pixels)
params.signalChange = 6;       % Threshold for intensity change detection 5-7

% Camera parameters of 95b sCMOS
% params.camera.gain = 0.89;      % electrons/ADU
% params.camera.QE = 0.92;        % Quantum efficiency
% Camera parameters of OCRA-Flash4.0 v3
params.camera.gain = 0.46;      % electrons/ADU
params.camera.QE = 0.75;        % Quantum efficiency

%Camera settings of this experiment
params.camera.exposure_time = 0.02;  % Exposure time in seconds
params.camera.frame_interval = 0.02; % Time between frames in seconds
params.EM_wave_FP = 660;

% Pixel size information
params.IM_info = [110, 110];    % Pixel size in nm [x, y]

% Paths
params.totalfn_path = 'D:\CityU\3.self-blinking dye\code\Dye_code-main_0520\Dye_code-main_0520\totalfn';
params.baseDir = 'D:\CityU\3.self-blinking dye\20251123\JF646-Hoechst\mv';  % Add baseDir to params

% Generate timestamp for this analysis run
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
params.timestamp = timestamp;
% Create a unique run ID for this analysis session
params.runID = ['Run_', timestamp];

% Load utility functions
utils = DyeAnalysisUtils();

%% Initial User Input - Select Analysis Mode and Datasets
%% Select analysis mode
analysisMode = questdlg('Select analysis mode:', 'Analysis Mode', ...
    'Localization Only', 'Bleach Curve Analysis Only', 'Localization Only');

if isempty(analysisMode)
    error('No analysis mode selected. Analysis canceled.');
end

%% Get datasets to process - Different options based on analysis mode
datasetType = '';
if strcmp(analysisMode, 'Localization Only')
    % For Localization Only, only allow single dataset
    datasetType = 'Single dataset';
    fprintf('Localization mode selected - processing single dataset only\n');
else
    % For Bleach Curve Analysis, allow multiple datasets
    datasetType = questdlg('Select dataset type:', 'Dataset Type', ...
        'Search in folder', 'Single dataset', 'Search in folder');
    
    if isempty(datasetType)
        error('No dataset type selected. Analysis canceled.');
    end
end

switch datasetType
    case 'Search in folder'
        % Select parent folder to search
        searchDir = uigetdir(params.baseDir, 'Select folder to search for datasets');
        if searchDir == 0
            error('No folder selected. Analysis canceled.');
        end
        
        % Find all datasets recursively
        foundDatasets = utils.findDatasets(searchDir);
        
        if isempty(foundDatasets)
            error('No datasets found in the selected directory and its subdirectories.');
        end
        
        % Create display names for selection
        displayNames = cell(length(foundDatasets), 1);
        for i = 1:length(foundDatasets)
            [~, name, ext] = fileparts(foundDatasets{i});
            if isdir(foundDatasets{i})
                displayNames{i} = sprintf('TIFF Sequence: %s', name);
            else
                displayNames{i} = sprintf('ND2 File: %s%s', name, ext);
            end
        end
        
        % Let user select which datasets to process
        [selection, ok] = listdlg('ListString', displayNames, ...
            'SelectionMode', 'multiple', ...
            'Name', 'Select Datasets', ...
            'PromptString', 'Select datasets to analyze:', ...
            'ListSize', [400 300]);
        
        if ~ok || isempty(selection)
            error('No datasets selected. Analysis canceled.');
        end
        
        % Store selected datasets
        datasetFolders = foundDatasets(selection);
        params.baseDir = searchDir;
        
    case 'Single dataset'
        % Original single dataset selection code
        datasetType = questdlg('Select dataset type:', 'Dataset Type', ...
            'Folder with TIFF sequence', 'ND2 file', 'ND2 file');
        
        switch datasetType
            case 'Folder with TIFF sequence'
                folderPath = uigetdir(params.baseDir, 'Select folder containing TIFF sequence');
                if folderPath == 0
                    error('No folder selected. Analysis canceled.');
                end
                datasetFolders = {folderPath};
                params.baseDir = fileparts(folderPath);
                
            case 'ND2 file'
                [fileName, filePath] = uigetfile({'*.nd2', 'Nikon ND2 Files (*.nd2)'}, ...
                    'Select ND2 file', params.baseDir);
                if fileName == 0
                    error('No file selected. Analysis canceled.');
                end
                nd2FullPath = fullfile(filePath, fileName);
                datasetFolders = {nd2FullPath};
                params.baseDir = filePath;
                
            otherwise
                error('Invalid dataset type selection. Analysis canceled.');
        end
end

% Define output base directory
% for multiple datasets will be the user selected folder
% for single dataset will be the parent folder of the dataset
outputBaseDir = fullfile(params.baseDir, 'Analysis_Results');
if ~exist(outputBaseDir, 'dir')
    mkdir(outputBaseDir);
    fprintf('Created output base directory: %s\n', outputBaseDir);
end

%% Process each dataset with SIMPLIFIED DATA FLOW
results = struct();
allResults = struct();

% Add totalfn path to MATLAB path
addpath(params.totalfn_path);

for i = 1:length(datasetFolders)
    % Current dataset path - points to the raw data location
    currentDataPath = datasetFolders{i};
    
    % Initialize results for this dataset
    currentResults = struct();
    
    % Determine output directory <<name>> based on input type
    [currentDataParentPath, currentDataName, ext] = fileparts(currentDataPath);
    % for .nd2, outputFolderName prefix will be currentDataName
    if ~isempty(ext) && strcmpi(ext, '.nd2')
        % For ND2 files, use the filename without extension
        outputFolderName = currentDataName;
    elseif contains(lower(currentDataName), 'sliced')
        % For sliced folders, use the parent folder name
        [~, parentFolder] = fileparts(currentDataParentPath);
        outputFolderName = parentFolder;
    else
        % For regular folders, use the folder name
        outputFolderName = currentDataName;
    end
    % Clean up folder name to remove invalid characters for filesystem
    outputFolderName = regexprep(outputFolderName, '[\\/:*?"<>|]', '_');
    % Use cleaned folder name for display and organization
    datasetName = outputFolderName;

    % Create analysis tag based on selected mode
    analysisTag = '';
    switch analysisMode
        case 'Localization Only'
            analysisTag = sprintf('_LOC_MSK%d_THR%d', params.mask_sz, params.threshold);
        case 'Bleach Curve Analysis Only'
            analysisTag = '_BCA';
    end

    analysisTag = [analysisTag, '_HFSQR', num2str(params.halfsquare),'_THR',num2str(params.signalChange)];

    % Define parent directory Analysis_Results path (for loading)
    parentAnalysisDir = fullfile(currentDataParentPath, 'Analysis_Results');
    
    % Define base directory Analysis_Results path (for saving in multiple datasets mode)
    baseAnalysisDir = fullfile(params.baseDir, 'Analysis_Results');
    
    % Determine if we're in multiple datasets mode
    isMultipleDatasets = strcmp(datasetType, 'Search in folder');

    % Set output directory based on dataset mode
    if isMultipleDatasets
        % For multiple datasets, save results to base directory
        outputDir = fullfile(baseAnalysisDir, [outputFolderName, analysisTag, '_', params.runID]);
    else
        % For single dataset, save results to parent directory
        outputDir = fullfile(parentAnalysisDir, [outputFolderName, analysisTag, '_', params.runID]);
    end
    
    % Store these values in params for use in helper functions
    params.outputBaseDir = baseAnalysisDir; %base directory Analysis_Results path (for saving in multiple datasets mode)
    params.parentAnalysisDir = parentAnalysisDir;% Define parent directory Analysis_Results path (for loading)
    params.outputFolderName = outputFolderName;
    params.currentDataPath = currentDataPath;
    params.currentDataParentPath = currentDataParentPath;
    params.isMultipleDatasets = isMultipleDatasets;
    % Store directory structure in params for consistent access
    params.outputDir = outputDir; % with tag and timestamp
    params.localizationDir = fullfile(outputDir, '2_Localization');
    params.blinkingDir = fullfile(outputDir, '3_BlinkingAnalysis');
    params.visualizationDir = fullfile(outputDir, '4_Visualizations');
    % params.statsDir = fullfile(params.blinkingDir, 'Statistics');
    % Define common data directory in parent folder (for loading)
    params.parentCommonDataDir = fullfile(parentAnalysisDir, [outputFolderName, '_CommonData']);
    params.parentCommonImageDataDir = fullfile(params.parentCommonDataDir, 'ImageData');

    % Create main output directory
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
        fprintf('Created output directory: %s\n', outputDir);
    end

    fprintf('Processing dataset %d/%d: %s\n', i, length(datasetFolders), datasetName);
    
    % Create progress log file
    progressLogFile = fullfile(outputDir, 'analysis_progress.log');
    params.progressLogFile = progressLogFile;
    
    fid = fopen(progressLogFile, 'w');
    fprintf(fid, 'Analysis started at: %s\n', datestr(now));
    fprintf(fid, 'Dataset: %s\n', datasetName);
    fprintf(fid, 'Analysis Mode: %s\n', analysisMode);
    fprintf(fid, 'Dataset Mode: %s\n', datasetType);
    fprintf(fid, 'Dataset Path: %s\n', currentDataPath);
    fprintf(fid, 'Parent Path: %s\n', currentDataParentPath);
    
    % Log all parameters in a compact format
    fprintf(fid, '=== Analysis Parameters ===\n');
    fprintf(fid, 'Localization: thr=%d, mask=%d, dil=%d, dist=%.1f, skip=%d\n', ...
        params.threshold, params.mask_sz, params.dil_sz, params.threshold_distance, ...
        params.frames_to_skip);
    fprintf(fid, 'Blinking: ROI=%d, sigChange=%d\n', ...
        params.halfsquare, params.signalChange);
    fprintf(fid, 'Camera: gain=%.2f, QE=%.2f, exp=%.3fs, int=%.3fs, λ=%dnm\n', ...
        params.camera.gain, params.camera.QE, params.camera.exposure_time, ...
        params.camera.frame_interval, params.EM_wave_FP);
    fprintf(fid, 'Pixel: x=%d nm, y=%d nm\n', params.IM_info(1), params.IM_info(2));
    fprintf(fid, 'Run ID: %s\n', params.runID);
    fprintf(fid, '======================\n\n');
    fclose(fid);

    %% SIMPLIFIED DATA FLOW - Stage 1: Image Data
    % Only load image data if needed for the selected analysis mode
    imageData = [];
    metadata = [];
    
    if strcmp(analysisMode, 'Localization Only') || ...
        (strcmp(analysisMode, 'Bleach Curve Analysis Only'))
         
         % Set flag to automatically use latest image data
         params.autoUseLatest = true;
         
         % Try to load existing image data
         utils.updateProgress(progressLogFile, 'Image Data', 'Searching for image data...');
        % For Bleach Curve Analysis, we need to load from parent directory
        utils.updateProgress(progressLogFile, 'Image Data', sprintf('Looking in parent directory: %s', params.parentCommonImageDataDir));
        % Try to load image data from parent directory's CommonData folder
        if exist(params.parentCommonImageDataDir, 'dir')
            % Find the latest image data file
            imageDataFiles = dir(fullfile(params.parentCommonImageDataDir, 'image_data*.mat'));
            
            if ~isempty(imageDataFiles)
                % Sort by date to get the latest
                [~, idx] = sort([imageDataFiles.datenum], 'descend');
                latestImageFile = fullfile(imageDataFiles(idx(1)).folder, imageDataFiles(idx(1)).name);
                
                try
                    utils.updateProgress(progressLogFile, 'Image Data', sprintf('Loading from %s', latestImageFile));
                    loadedData = load(latestImageFile);
                    if isfield(loadedData, 'imageData')
                        imageData = loadedData.imageData;
                        if isfield(loadedData, 'metadata')
                            metadata = loadedData.metadata;
                        end
                        utils.updateProgress(progressLogFile, 'Image Data', sprintf('Successfully loaded image data with dimensions: %dx%dx%d', ...
                            size(imageData, 1), size(imageData, 2), size(imageData, 3)));
                    end
                catch ME
                    utils.updateProgress(progressLogFile, 'Warning', sprintf('Error loading image data: %s', ME.message));
                end
            else
                utils.updateProgress(progressLogFile, 'Warning', 'No image data files found in parent CommonData directory');
            end
        else
            utils.updateProgress(progressLogFile, 'Warning', sprintf('Parent CommonData directory not found: %s', params.parentCommonImageDataDir));
        end
        
        % Only process new image data for Localization Only mode, not for Bleach Curve Analysis
        if isempty(imageData) && strcmp(analysisMode, 'Localization Only')
            % If not loaded, process new image data
            utils.updateProgress(progressLogFile, 'Image Data', 'Processing new image data...');
            [imageData, metadata] = utils.manageImageData(params, progressLogFile);
        elseif isempty(imageData) && strcmp(analysisMode, 'Bleach Curve Analysis Only')
            utils.updateProgress(progressLogFile, 'Error', 'Required image data not found for Bleach Curve Analysis. Run Localization Only first.');
            fprintf('Error: No image data found for dataset: %s. Skipping analysis.\n', datasetName);
            continue;
        else
            utils.updateProgress(progressLogFile, 'Image Data', 'Using existing image data');
        end
    end
    
    %% SIMPLIFIED DATA FLOW - Stage 2: Localization
    locResults = [];
    
    if strcmp(analysisMode, 'Localization Only') || strcmp(analysisMode, 'Bleach Curve Analysis Only')
        % For Localization Only, always perform new localization
        % For Bleach Curve Analysis, try to load existing localization results


        if strcmp(analysisMode, 'Localization Only')
            % In Localization Only mode, always perform new localization
            utils.updateProgress(progressLogFile, 'Localization', 'Performing new localization...');
            
            % Create localization directory if it doesn't exist
            if ~exist(params.localizationDir, 'dir')
                mkdir(params.localizationDir);
                utils.updateProgress(progressLogFile, 'Directory', 'Created localization directory');
            end
            
            % Perform localization
            [locResults, localizationFrames] = utils.handleLocalization(imageData, params, outputDir, progressLogFile);
            
            currentResults.localization = locResults;
            
            % Save and continue to next dataset
            utils.updateProgress(progressLogFile, 'Analysis', 'Localization analysis completed');
            utils.saveAnalysisResults(currentResults, outputDir, params, allResults, i, progressLogFile);
            fprintf('Completed localization analysis for dataset: %s\n', datasetName);
            continue;
        else
            % For Bleach Curve Analysis, we need to load from parent directory
            utils.updateProgress(progressLogFile, 'Localization', 'Searching for localization results...');
            
            % Search for localization results in parent directory
            locResultsFound = false;
            
            % Look for folders with LOC tag in parent Analysis_Results
            locFolders = dir(fullfile(parentAnalysisDir, [outputFolderName, '_LOC_*']));

            if ~isempty(locFolders)
                % Sort by date to get the latest
                [~, idx] = sort([locFolders.datenum], 'descend');
                
                % Try each folder until we find valid results
                for folderIdx = 1:length(idx)
                    latestLocFolder = fullfile(locFolders(idx(folderIdx)).folder, locFolders(idx(folderIdx)).name);
                    locResultsFile = fullfile(latestLocFolder, '2_Localization', 'localization_results.mat');
                    
                    if exist(locResultsFile, 'file')
                        try
                            utils.updateProgress(progressLogFile, 'Localization', sprintf('Loading from %s', locResultsFile));
                            loadedData = load(locResultsFile);
                            if isfield(loadedData, 'locResults')
                                locResults = loadedData.locResults;
                                locResultsFound = true;
                                utils.updateProgress(progressLogFile, 'Localization', 'Successfully loaded localization results');
                                break;
                            end
                        catch ME
                            utils.updateProgress(progressLogFile, 'Warning', sprintf('Error loading localization results: %s', ME.message));
                        end
                    end
                end
            else
                utils.updateProgress(progressLogFile, 'Warning', 'No localization folders found in parent Analysis_Results directory');
            end
            
            if ~locResultsFound
                utils.updateProgress(progressLogFile, 'Error', 'No localization data found. Cannot perform blinking analysis.');
                fprintf('Error: No localization data found for dataset: %s. Skipping analysis.\n', datasetName);
                continue;
            else
                utils.updateProgress(progressLogFile, 'Localization', 'Using existing localization results');
            end
        end
    end

    %% SIMPLIFIED DATA FLOW - Stage 3: Blinking Analysis
    if strcmp(analysisMode, 'Bleach Curve Analysis Only')
        % Create blinking directory if it doesn't exist
        if ~exist(params.blinkingDir, 'dir')
            mkdir(params.blinkingDir);
            utils.updateProgress(progressLogFile, 'Directory', 'Created blinking analysis directory');
        end
        
        % Create statistics directory if it doesn't exist
        params.statsDir = fullfile(params.blinkingDir, 'Statistics');
        if ~exist(params.statsDir, 'dir')
            mkdir(params.statsDir);
            utils.updateProgress(progressLogFile, 'Directory', 'Created statistics directory');
        end
    
        % Perform new blinking analysis
        utils.updateProgress(progressLogFile, 'Blinking Analysis', 'Performing new blinking analysis...');
        
        % Perform blinking analysis with selected method
        [blinkResults, ~] = utils.handleBlinkingAnalysis(imageData, locResults, params, outputDir, progressLogFile);
    
        currentResults.blinking = blinkResults;
    
        % Visualize results
        if isfield(blinkResults, 'validatedPoints') && ~isempty(blinkResults.validatedPoints)
            utils.updateProgress(progressLogFile, 'Visualization', 'Generating visualizations...');
            
            % Create visualization directory if it doesn't exist
            if ~exist(params.visualizationDir, 'dir')
                mkdir(params.visualizationDir);
                utils.updateProgress(progressLogFile, 'Directory', 'Created visualization directory');
            end
            
            % Determine results file path
            blinkResultsFile = fullfile(params.blinkingDir, 'blinking_results.mat');
            
            % Check if the file exists before trying to visualize
            if exist(blinkResultsFile, 'file')
                try
                    % Call visualization function with error handling
                    utils.handleErrorWithFallback(progressLogFile, @() visualizeResults(blinkResultsFile), ...
                        'Visualization', 'Error in visualization', 'Will attempt basic visualization');
                    
                    utils.updateProgress(progressLogFile, 'Visualization', 'Visualizations complete');
                catch ME
                    utils.updateProgress(progressLogFile, 'Warning', sprintf('Visualization error: %s', ME.message));
                end
            else
                utils.updateProgress(progressLogFile, 'Warning', sprintf('Results file not found: %s', blinkResultsFile));
            end
        else
            utils.updateProgress(progressLogFile, 'Warning', 'No points were validated for blinking analysis');
        end
    end
    
    %% Save final results for this dataset
    utils.saveAnalysisResults(currentResults, outputDir, params, allResults, i, progressLogFile);

    % Store results in the main results structure
    results.(genvarname(['dataset_' num2str(i)])) = currentResults;

    fprintf('Completed analysis of dataset: %s\n', datasetName);
end

fprintf('Analysis complete for all datasets.\n');

% Remove totalfn path from MATLAB path
rmpath(params.totalfn_path);
