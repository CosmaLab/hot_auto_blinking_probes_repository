function [FinalImage, TotalFrame, metadata] = iQ_imgload3frames_file_input(input_data, frames_to_skip)
%% Ma, 2025-03-19
% Modify the code to image sequences from various sources including:
% - Directory paths containing TIFF files
% - Direct image arrays
% - ND2 files (if Bio-Formats is available)
% - File lists (struct arrays from dir function)
%
% Inputs:
%   input_data - Can be:
%                1. Path to directory containing TIFF files
%                2. Path to ND2 file
%                3. Direct image array [height x width x frames]
%                4. File list (struct array from dir function)
%   frames_to_skip - (Optional) Number of initial frames to skip
%
% Outputs:
%   FinalImage - 3D image stack [height x width x frames]
%   TotalFrame - Number of frames in the stack
%   metadata - Structure containing metadata about the loaded images

% Initialize metadata
metadata = struct();
metadata.source = '';
metadata.sourceType = '';
metadata.dimensions = [];
metadata.originalFrames = 0;
metadata.loadedFrames = 0;
metadata.skippedFrames = 0;
metadata.loadTime = 0;

% Set default for frames_to_skip
if nargin < 2 || isempty(frames_to_skip)
    frames_to_skip = 0;
end

% Record start time
tic;

% Check input type
if isnumeric(input_data)
    % Direct image array input
    FinalImage = input_data;
    TotalFrame = size(FinalImage, 3);
    
    % Update metadata
    metadata.sourceType = 'ImageArray';
    metadata.dimensions = size(FinalImage);
    metadata.originalFrames = TotalFrame;
    metadata.loadedFrames = TotalFrame;
    metadata.loadTime = toc;
    return;
end

% Handle ND2 files
if ischar(input_data) || isstring(input_data)
    [~, ~, ext] = fileparts(input_data);
    if strcmpi(ext, '.nd2')
        % Check if Bio-Formats is available
        if ~exist('bfGetReader', 'file')
            error('Bio-Formats package is required to read ND2 files but was not found in the MATLAB path.');
        end
        
        try
            % Load ND2 file using Bio-Formats
            fprintf('Loading ND2 file: %s\n', input_data);
            reader = bfGetReader(input_data);
            omeMeta = reader.getMetadataStore();
            
            % Get dimensions
            width = omeMeta.getPixelsSizeX(0).getValue();
            height = omeMeta.getPixelsSizeY(0).getValue();
            totalFrames = omeMeta.getPixelsSizeT(0).getValue();
            
            % Apply frames_to_skip
            startFrame = frames_to_skip + 1;
            framesToLoad = totalFrames - frames_to_skip;
            
            if framesToLoad <= 0
                error('No frames to load after skipping %d frames from a total of %d frames.', frames_to_skip, totalFrames);
            end
            
            % Initialize output array
            FinalImage = zeros(height, width, framesToLoad, 'uint16');
            
            % Load frames
            fprintf('Loading %d frames (skipping first %d)...\n', framesToLoad, frames_to_skip);
            for i = 1:framesToLoad
                frameIdx = startFrame + i - 1;
                reader.setSeries(0);
                iPlane = reader.getIndex(0, 0, frameIdx - 1) + 1;
                FinalImage(:,:,i) = bfGetPlane(reader, iPlane);
            end
            
            % Close reader
            reader.close();
            
            % Set output variables
            TotalFrame = framesToLoad;
            
            % Update metadata
            metadata.sourceType = 'ND2File';
            metadata.source = input_data;
            metadata.dimensions = size(FinalImage);
            metadata.originalFrames = totalFrames;
            metadata.loadedFrames = framesToLoad;
            metadata.skippedFrames = frames_to_skip;
            metadata.loadTime = toc;
            
            fprintf('Successfully loaded %d frames from ND2 file.\n', TotalFrame);
            return;
        catch ME
            error('Error loading ND2 file: %s', ME.message);
        end
    end
end

% Handle directory or file list input for TIFF files
if isstruct(input_data)
    fileList = input_data;
    if ~isempty(fileList)
        [dir_path, ~, ~] = fileparts(fullfile(fileList(1).folder, fileList(1).name));
        metadata.source = dir_path;
    else
        error('Empty file list provided');
    end
else
    if iscell(input_data)
        input_data = input_data{1};
    end
    dir_path = input_data;
    metadata.source = dir_path;
    fileList = dir(fullfile(dir_path, '*.tif'));
    
    if isempty(fileList)
        % Try alternative extensions
        fileList = dir(fullfile(dir_path, '*.tiff'));
        if isempty(fileList)
            error('No TIFF files found in directory: %s', dir_path);
        end
    end
end

% Apply frames_to_skip
totalAvailableFrames = length(fileList);
if frames_to_skip >= totalAvailableFrames
    error('Cannot skip %d frames when only %d frames are available.', frames_to_skip, totalAvailableFrames);
end

% Adjust file list to skip frames
fileList = fileList(frames_to_skip+1:end);
TotalFrame = length(fileList);

fprintf('Loading %d TIFF frames from %s (skipping first %d)...\n', TotalFrame, metadata.source, frames_to_skip);

% Load first image to get dimensions
firstImage = imread(fullfile(fileList(1).folder, fileList(1).name));
[nImage, mImage, ~] = size(firstImage);
FinalImage = zeros(nImage, mImage, TotalFrame, 'uint16');

% Load all images
warning('off', 'all');
parfor j = 1:TotalFrame
    curr_frm = imread(fullfile(fileList(j).folder, fileList(j).name));
    if size(curr_frm, 3) > 1
        FinalImage(:,:,j) = curr_frm(:,:,1);
    else
        FinalImage(:,:,j) = curr_frm;
    end
end
warning('on', 'all');

% Update metadata
metadata.sourceType = 'TIFFSequence';
metadata.dimensions = size(FinalImage);
metadata.originalFrames = totalAvailableFrames;
metadata.loadedFrames = TotalFrame;
metadata.skippedFrames = frames_to_skip;
metadata.loadTime = toc;

fprintf('Successfully loaded %d frames from TIFF sequence.\n', TotalFrame);
end