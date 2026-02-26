function Tiff2Mat(FolderDir,SaveDir)
    % Get list of TIFF files
    % The fullfile() function is used to construct the full path 
    % for each TIFF file before calling imread(). 
    fileInfo = dir(fullfile(FolderDir, '*.tif')); 
    if isempty(fileInfo)
        error('No TIFF files found in directory');
    end

    % Get number of frames
    N_Frames = length(fileInfo);
    % fileInfo = dir(FolderDir); 
    % TifName = {fileInfo(~[fileInfo.isdir]).name}';
    % [N_Frames,~] = size(TifName);    

    % Get image dimensions from the first frame
    firstFramePath = fullfile(FolderDir, fileInfo(1).name);
    firstFrame = imread(firstFramePath);
    [H, W] = size(firstFrame);

    % Initialize output array
    mov = zeros(H, W, N_Frames, 'uint16');

    % Temp_Dir = append(FolderDir,'/',TifName(1));
    % Temp_Frame = imread(Temp_Dir);
    % [H,W] = size(Temp_Frame);
    % mov = zeros(H,W,N_Frames); mov = uint16(mov);
    
    % Load all frames
    parfor i = 1:N_Frames
        currentFramePath = fullfile(FolderDir, fileInfo(i).name);  % Ensure full path
        mov(:,:,i) = imread(currentFramePath);  % Read the current frame
    end
    % parfor i = 1 : N_Frames
    %     fprintf('\r%d', i);
    %     Temp_Dir = append(FolderDir,'/',TifName(i));
    %     Temp_Frame = imread(Temp_Dir);
    %     mov(:,:,i) = Temp_Frame;
    % end
    

% Removed the check for the metadata file 'metadata.txt' for the following reasons:
% 1. **Unnecessary Complexity**: The presence of a metadata file may not be relevant to the processing of TIFF files, 
%    and checking for it adds unnecessary complexity to the code.
% 2. **Potential Errors**: If the variable TifName is not properly defined or populated before this check, 
%    it could lead to runtime errors.
% 3. **Clarity**: The logic for handling the metadata file may not be clear to all users, especially if they are 
%    unfamiliar with its purpose. Removing this check simplifies the code and makes it easier to understand.
% 4. **Path Handling**: The original method of constructing the file path using 'append' could lead to issues 
%    with cross-platform compatibility. Using 'fullfile' is a more robust approach.
%
% Therefore, the check for 'metadata.txt' has been removed to streamline the code and focus on processing 
% only the relevant TIFF files.
    % if exist(append(FolderDir,"/metadata.txt"),"file")
    %     RemoveSth = 'metadata.txt';
    %     TifName = string(TifName);
    %     TifName = setdiff (TifName,RemoveSth);
    %     N_Frames = N_Frames - 1;
    % end
    

    fprintf('\n');
    disp("Saving .mat file, please wait.");
    save(SaveDir,"mov","-v7.3","-nocompression");
    disp("MV saved as .mat file");
    clear all;
end