% DyeAnalysisUtils.m
% Utility functions for Dye Analysis Main script

function utils = DyeAnalysisUtils()
    % Return a structure with function handles to all utility functions
    utils.updateProgress = @updateProgress;
    utils.findResultFile = @findResultFile; % New unified file search function
    utils.loadDataStage = @loadDataStage; % New unified data loading function
    utils.selectFrameRange = @selectFrameRange;
    utils.handleLocalization = @handleLocalization;
    utils.handleBleachCurveAnalysis = @handleBleachCurveAnalysis;
    utils.handleBlinkingAnalysis = @handleBlinkingAnalysis;  % Add new function handle
    utils.manageImageData = @manageImageData;
    utils.importBrightFieldROI = @importBrightFieldROI;
    utils.findDatasets = @findDatasets;
    utils.saveAnalysisResults = @saveAnalysisResults;
    utils.handleErrorWithFallback = @handleErrorWithFallback;  % New standardized error handling
end

% Progress update function
function updateProgress(logFile, step, details)
    fid = fopen(logFile, 'a');
    fprintf(fid, '[%s] %s: %s\n', char(datetime('now', 'Format', 'HH:mm:ss')), step, details);
    fclose(fid);
    fprintf('[%s] %s: %s\n', char(datetime('now', 'Format', 'HH:mm:ss')), step, details);
end

% Function to handle frame range selection
function [startFrame, endFrame] = selectFrameRange(imageData, defaultStart, defaultEnd, prompt, title)
    frameSelectionChoice = questdlg('How would you like to select frames?', ...
        'Frame Selection', 'All frames', 'Specific range', 'All frames');
    
    startFrame = defaultStart;
    endFrame = defaultEnd;
    
    if strcmp(frameSelectionChoice, 'Specific range')
        prompt = {sprintf('Start frame (min: 1):'), ...
                  sprintf('End frame (max: %d):', size(imageData, 3))};
        dims = [1 50];
        defaultInput = {num2str(defaultStart), num2str(defaultEnd)};
        answer = inputdlg(prompt, title, dims, defaultInput);
        
        if ~isempty(answer)
            userStartFrame = str2double(answer{1});
            if ~isnan(userStartFrame) && userStartFrame >= 1 && userStartFrame <= size(imageData, 3)
                startFrame = userStartFrame;
            end
            
            userEndFrame = str2double(answer{2});
            if ~isnan(userEndFrame) && userEndFrame >= startFrame && userEndFrame <= size(imageData, 3)
                endFrame = userEndFrame;
            end
        end
    end
end

function [imageData, metadata] = manageImageData(params, progressLogFile)
    % Function to load or process image data
    % Updated to work with the new directory structure
    
    imageData = [];
    metadata = struct();
    
    % Ensure frames_to_skip is set
    if ~isfield(params, 'frames_to_skip')
        params.frames_to_skip = 0;
        updateProgress(progressLogFile, 'Parameter', 'Setting default frames_to_skip=0');
    end
    
    % Load image data from source
    updateProgress(progressLogFile, 'Data Loading', 'Loading image sequence from source...');
    [imageData, metadata] = iQ_imgload3frames_file_input(params.currentDataPath, params.frames_to_skip);
    
    if isempty(imageData) || numel(imageData) == 0
        updateProgress(progressLogFile, 'Error', 'Failed to load image data from source');
        return;
    end
    
    % Determine where to save the data - in parent directory's CommonData folder
    parentCommonDataDir = fullfile(params.currentDataParentPath, 'Analysis_Results', [params.outputFolderName, '_CommonData']);
    parentImageDataDir = fullfile(parentCommonDataDir, 'ImageData');
    parentImageDataFile = fullfile(parentImageDataDir, ['image_data_', params.timestamp, '.mat']);
    
    % Create directories if they don't exist
    if ~exist(parentCommonDataDir, 'dir')
        mkdir(parentCommonDataDir);
        updateProgress(progressLogFile, 'Directory', sprintf('Created parent common data directory: %s', parentCommonDataDir));
    end
    
    if ~exist(parentImageDataDir, 'dir')
        mkdir(parentImageDataDir);
        updateProgress(progressLogFile, 'Directory', sprintf('Created parent image data directory: %s', parentImageDataDir));
    end
    
    % Save to parent directory
    try
        updateProgress(progressLogFile, 'Data Saving', sprintf('Saving image data to: %s', parentImageDataFile));
        save(parentImageDataFile, 'imageData', 'metadata', '-v7.3');
        
        % Also save a standard version for backward compatibility
        standardImageDataFile = fullfile(parentImageDataDir, 'image_data.mat');
        save(standardImageDataFile, 'imageData', 'metadata', '-v7.3');
        
        updateProgress(progressLogFile, 'Data Saving', 'Saved image data to parent directory for future use');
    catch ME
        updateProgress(progressLogFile, 'Warning', sprintf('Could not save to parent directory: %s', ME.message));
    end
    
    updateProgress(progressLogFile, 'Data Loading', sprintf('Successfully loaded image data with dimensions: %dx%dx%d', ...
        size(imageData, 1), size(imageData, 2), size(imageData, 3)));
end

% Function to handle localization process
function [locResults, localizationFrames] = handleLocalization(imageData, params, outputDir, progressLogFile)
    
    % Create localization directory only if needed
    if ~exist(params.localizationDir, 'dir')
        mkdir(params.localizationDir);
        updateProgress(progressLogFile, 'Directory', 'Created localization directory');
    end
    
    % Set output directory in params
    params.outputDir = outputDir;
    
    % Check for existing MRB in common data folder
    commonDataDir = fullfile(params.outputBaseDir, [params.outputFolderName, '_CommonData']);
    commonMRBDir = fullfile(commonDataDir, 'ROI_Data');
    commonMRBFile = fullfile(commonMRBDir, 'localization_ROIs.mat');
    commonROIMaskFile = fullfile(commonMRBDir, 'roi_mask.mat');

    MRB = {};
    mrbLoaded = false;
    
    % SIMPLIFIED ROI HANDLING LOGIC
    % 1. First check if ROI mask already exists in params
    if isfield(params, 'roiMask') && ~isempty(params.roiMask)
        updateProgress(progressLogFile, 'ROI Selection', 'Using ROI mask already in params');
        mrbLoaded = true;
        
        % If we have the mask but not MRB, try to load MRB
        if (~isfield(params, 'MRB') || isempty(params.MRB)) && exist(commonMRBFile, 'file')
            try
                loadedData = load(commonMRBFile);
                if isfield(loadedData, 'MRB') && ~isempty(loadedData.MRB)
                    params.MRB = loadedData.MRB;
                    updateProgress(progressLogFile, 'ROI Selection', 'Loaded MRB to complement existing ROI mask');
                end
            catch
                % Continue with just the ROI mask
            end
        end
    % 2. If no ROI mask in params, try to load from common data folder
    elseif exist(commonROIMaskFile, 'file')
        try
            loadedData = load(commonROIMaskFile);
            if isfield(loadedData, 'roiMask') && ~isempty(loadedData.roiMask)
                params.roiMask = loadedData.roiMask;
                updateProgress(progressLogFile, 'ROI Selection', 'Loaded ROI mask from common data folder');
                mrbLoaded = true;
                
                % Also try to load MRB if available
                if exist(commonMRBFile, 'file')
                    try
                        mrbData = load(commonMRBFile);
                        if isfield(mrbData, 'MRB') && ~isempty(mrbData.MRB)
                            params.MRB = mrbData.MRB;
                            updateProgress(progressLogFile, 'ROI Selection', 'Loaded MRB from common data folder');
                        end
                    catch
                        % Continue with just the ROI mask
                    end
                end
            end
        catch
            % Failed to load ROI mask, continue to next option
        end
    % 3. If no ROI mask, try to load MRB and create mask
    elseif exist(commonMRBFile, 'file')
        % Ask user if they want to use existing MRB or create a new one
        choice = questdlg('ROI data found. Would you like to use existing ROIs or draw new ones?', ...
            'ROI Selection', 'Use existing ROIs', 'Draw new ROIs', 'Use existing ROIs');
        
        if strcmp(choice, 'Use existing ROIs')
            try
                loadedData = load(commonMRBFile);
                if isfield(loadedData, 'MRB') && ~isempty(loadedData.MRB)
                    MRB = loadedData.MRB;
                    params.MRB = MRB;
                    mrbLoaded = true;
                    
                    % Create ROI mask from MRB
                    params.roiMask = createRoiMaskFromMRB(MRB, size(imageData), params.dil_sz);
                    
                    % Save ROI mask for future use
                    if ~exist(commonMRBDir, 'dir')
                        mkdir(commonMRBDir);
                    end
                    roiMask = params.roiMask;
                    save(commonROIMaskFile, 'roiMask', '-v7.3');
                    
                    updateProgress(progressLogFile, 'ROI Selection', ...
                        sprintf('Created and saved ROI mask from existing MRB with %d ROIs', length(MRB)));
                else
                    updateProgress(progressLogFile, 'Warning', 'MRB file exists but contains no valid ROI data.');
                    mrbLoaded = false;
                end
            catch ME
                updateProgress(progressLogFile, 'Warning', ...
                    sprintf('Error loading MRB from common data: %s', ME.message));
                mrbLoaded = false;
            end
        else
            updateProgress(progressLogFile, 'ROI Selection', 'User chose to draw new ROIs.');
            mrbLoaded = false;
        end
    end
    
    % 4. If still no ROI data, prompt user to create new ones
    if ~mrbLoaded
        % Ask if user wants to define ROIs using a bright field image
        useROI = questdlg('Do you want to define regions of interest (ROIs) for localization?', ...
            'ROI Selection', 'Yes, import bright field image', 'No, analyze entire image', 'No, analyze entire image');
        
        if strcmp(useROI, 'Yes, import bright field image')
            updateProgress(progressLogFile, 'ROI Selection', 'Importing bright field image for ROI selection...');
            [MRB, brightFieldImage, allMasks] = importBrightFieldROI(params);
            
            if ~isempty(MRB)
                % Store MRB in params
                params.MRB = MRB;
                
                % Create ROI mask from MRB
                params.roiMask = allMasks;
                
                % Save both MRB and ROI mask
                if ~exist(commonMRBDir, 'dir')
                    mkdir(commonMRBDir);
                end
                
                % Save MRB
                save(commonMRBFile, 'brightFieldImage', 'MRB', '-v7.3');
                
                % Save ROI mask separately
                roiMask = params.roiMask;
                save(commonROIMaskFile, 'roiMask', '-v7.3');
                
                updateProgress(progressLogFile, 'ROI Selection', ...
                    sprintf('Using %d user-defined ROIs for localization.', length(MRB)));
            else
                updateProgress(progressLogFile, 'ROI Selection', 'No ROIs defined. Will analyze entire image.');
                params.MRB = [];
                
                % Create full-image mask
                params.roiMask = true(size(imageData, 1), size(imageData, 2));
            end
        else
            % User chose to analyze entire image
            params.MRB = [];
            params.roiMask = true(size(imageData, 1), size(imageData, 2));
            updateProgress(progressLogFile, 'ROI Selection', 'User chose to analyze entire image.');
        end
    end
    
    % Select frame range for localization
    [startFrame, endFrame] = selectFrameRange(imageData, ...
        params.frames_to_skip + 1, size(imageData, 3), ...
        'Frame Range Selection', 'Localization Frame Range');
    
    localizationFrames = startFrame:endFrame;
    updateProgress(progressLogFile, 'Frame Selection', ...
        sprintf('Using frame range %d to %d for localization', startFrame, endFrame));
    
    % Perform localization
    updateProgress(progressLogFile, 'Step 2', 'Performing localization...');
    locResults = performLocalization(imageData, params, localizationFrames);
    
    % Save results
    locResultsFile = fullfile(outputDir, '2_Localization', 'localization_results.mat');
    save(locResultsFile, 'locResults', 'localizationFrames', '-v7.3');
    updateProgress(progressLogFile, 'Step 2', 'Saved localization results');
    
    % Log spot count
    if isfield(locResults, 'combinedPoints')
        spotCount = size(locResults.combinedPoints, 1);
    elseif isfield(locResults, 'allFilteredPoints')
        spotCount = size(locResults.allFilteredPoints, 1);
    else
        spotCount = 0;
    end
    updateProgress(progressLogFile, 'Step 2', sprintf('Localized %d spots', spotCount));
end

% centralize handle for MRB:
    % always saving both MRB and roiMask together in the same file.
    % When loading, checking if roiMask is available and creating it only if needed
% Helper function to create ROI mask from MRB
function roiMask = createRoiMaskFromMRB(MRB, imageSize, dil_sz)
    % Initialize mask with zeros (no ROIs)
    roiMask = false(imageSize(1:2));
    
    % Add each ROI to the mask
    for r = 1:size(MRB, 1)
        if ~isempty(MRB{r,1})
            % Create individual ROI mask
            gridd = false(size(roiMask));
            gridd(sub2ind(size(gridd), MRB{r,1}(:,2), MRB{r,1}(:,1))) = 1;
            
            % Fill holes and dilate
            grid1 = imfill(gridd, 'holes');
            se1 = strel('disk', dil_sz);
            grid2 = imdilate(grid1, se1);
            
            % Add to combined mask
            roiMask = roiMask | grid2;
        end
    end
    
    % Print coverage information
    roiCoverage = sum(roiMask(:)) / numel(roiMask) * 100;
    fprintf('Created ROI mask covering %.2f%% of the image area\n', roiCoverage);
end

% Helper function to save ROI data
function saveROIData(MRB, brightFieldImage, roiMask, outputDir, commonMRBDir, params, progressLogFile)
    % Save to localization directory
    locDir = fullfile(outputDir, '2_Localization');
    if ~exist(locDir, 'dir')
        mkdir(locDir);
    end
    
    % Save ROI data with timestamp
    mrbFile = fullfile(locDir, ['MRB_', params.timestamp, '.mat']);
    standardMrbFile = fullfile(locDir, 'MRB.mat');
    
    % Save to localization directory
    try
        save(mrbFile, 'brightFieldImage', 'roiMask', 'MRB', '-v7.3');
        save(standardMrbFile, 'brightFieldImage', 'roiMask', 'MRB', '-v7.3');
        updateProgress(progressLogFile, 'ROI Selection', 'Saved ROIs to localization directory');
    catch ME
        updateProgress(progressLogFile, 'Warning', ...
            sprintf('Failed to save ROIs to localization directory: %s', ME.message));
    end
    
    % Save to common data folder
    if ~exist(commonMRBDir, 'dir')
        mkdir(commonMRBDir);
        updateProgress(progressLogFile, 'Directory', 'Created common ROI data directory');
    end
    
    commonMRBFile = fullfile(commonMRBDir, 'localization_ROIs.mat');
    try
        save(commonMRBFile, 'brightFieldImage', 'roiMask', 'MRB', '-v7.3');
        updateProgress(progressLogFile, 'ROI Selection', 'Saved ROIs to common data folder for future analyses');
    catch ME
        updateProgress(progressLogFile, 'Warning', ...
            sprintf('Failed to save ROIs to common data folder: %s', ME.message));
    end
    
    % Verify the saved file
    if exist(mrbFile, 'file')
        try
            testLoad = load(mrbFile, 'MRB');
            if isfield(testLoad, 'MRB') && ~isempty(testLoad.MRB)
                updateProgress(progressLogFile, 'ROI Selection', ...
                    sprintf('Successfully saved and verified %d user-defined ROIs', length(MRB)));
            else
                updateProgress(progressLogFile, 'Warning', ...
                    'MRB file was saved but verification failed. ROIs may not be available for future runs.');
            end
        catch
            updateProgress(progressLogFile, 'Warning', ...
                'Failed to verify saved MRB file. ROIs may not be available for future runs.');
        end
    else
        updateProgress(progressLogFile, 'Warning', ...
            'Failed to save MRB file. ROIs may not be available for future runs.');
    end
end

% Function to save final results
function saveAnalysisResults(results, outputDir, params, allResults, i, progressLogFile)
    updateProgress(progressLogFile, 'Final Step', 'Saving final results...');
    
    % Save timestamped results
    finalResultsFile = fullfile(outputDir, ['complete_results_', params.timestamp, '.mat']);
    save(finalResultsFile, 'results', '-v7.3');
    
    % Store in all results
    allResults.(genvarname(['dataset_' num2str(i)])) = results;
    
    % Extract and save statistics
    saveStatistics(results, outputDir);

    updateProgress(progressLogFile, 'Complete', 'Analysis successfully completed');
end

function saveStatistics(results, outputDir)
    % Extract key statistics for CSV export
    stats = struct();
    
    % Extract statistics from results structure
    if isfield(results, 'localization')
        locData = results.localization;
        if isfield(locData, 'allFilteredPoints')
            stats.total_filtered_points = size(locData.allFilteredPoints, 1);
        end
        if isfield(locData, 'combinedPoints')
            stats.combined_points = size(locData.combinedPoints, 1);
        end
    end
    
    % Blinking stats
    if isfield(results, 'blinking')
        blinkData = results.blinking;
        if isfield(blinkData, 'numValidatedPoints')
            stats.validated_points = blinkData.numValidatedPoints;
        end
        if isfield(blinkData, 'blinking_on_duration_ms') && ~isempty(blinkData.blinking_on_duration_ms)
            stats.mean_on_time_ms = mean(blinkData.blinking_on_duration_ms);
            stats.median_on_time_ms = median(blinkData.blinking_on_duration_ms);
        end
        if isfield(blinkData, 'blinking_off_duration_ms') && ~isempty(blinkData.blinking_off_duration_ms)
            stats.mean_off_time_ms = mean(blinkData.blinking_off_duration_ms);
            stats.median_off_time_ms = median(blinkData.blinking_off_duration_ms);
        end
    end
    
    % Duty cycle stats
    if isfield(results, 'dutyCycle')
        dutyData = results.dutyCycle;
        if isfield(dutyData, 'duty_cycle') && ~isempty(dutyData.duty_cycle)
            stats.mean_duty_cycle = mean(dutyData.duty_cycle);
            stats.median_duty_cycle = median(dutyData.duty_cycle);
        end
    end
    
    % Convert to table and save as CSV
    if ~isempty(fieldnames(stats))
        statsTable = struct2table(stats, 'AsArray', true);
        statsFile = fullfile(outputDir, 'analysis_statistics.csv');
        writetable(statsTable, statsFile);
    end
end

function [MRB, brightFieldImage, allMasks] = importBrightFieldROI(params)
    % Import a bright field image and allow user to draw ROIs for localization
    % Returns:
    %   MRB - cell array of ROI coordinates compatible with localization functions
    %   brightFieldImage - the loaded bright field image
    %   allMasks - combined binary mask of all ROIs
    
    % Ask user to select a bright field image
    [fileName, filePath] = uigetfile({'*.tif;*.tiff', 'TIFF Images (*.tif, *.tiff)'}, ...
        'Select Bright Field Image', params.baseDir);
    
    if fileName == 0
        % User canceled, return empty MRB
        MRB = {};
        return;
    end
    
    % Load the bright field image
    brightFieldImage = imread(fullfile(filePath, fileName));
    
    % Display the image for ROI selection
    hFig = figure('Name', 'Draw ROIs for Localization', 'Position', [100, 100, 800, 600]);
    hAx = axes('Parent', hFig);
    imshow(brightFieldImage, [], 'Parent', hAx);
    title('Draw ROIs for localization. Double-click to finish each ROI.');
    
    % Initialize ROI storage
    roiCount = 0;
    MRB = cell(1,1);
    allMasks = false(size(brightFieldImage));
    roiHandles = struct('line', {}, 'label', {}, 'mask', {}, 'active', {});
    
    % Flag to control when to close the figure
    isFinished = false;
    
    % Add buttons for control
    uicontrol('Style', 'pushbutton', 'String', 'Add ROI', ...
        'Position', [10 10 100 30], ...
        'Callback', @addROI);
    uicontrol('Style', 'pushbutton', 'String', 'Finish', ...
        'Position', [120 10 100 30], ...
        'Callback', @finishROIs);
    uicontrol('Style', 'pushbutton', 'String', 'Clear All', ...
        'Position', [230 10 100 30], ...
        'Callback', @clearAllROIs);
    
    % Callback for adding new ROI
    function addROI(~,~)
        try
            h = imfreehand(hAx);
            if isempty(h)
                return;
            end
            
            % Create a position callback that will be triggered when the ROI is modified
            addNewPositionCallback(h, @(p) updatePosition(h));
            
            % Create a custom double-click listener
            set(hFig, 'WindowButtonDownFcn', @(src,evt)checkDoubleClick(h));
            
            % Wait without using the problematic wait() function
            while isvalid(h) && ~isempty(h)
                drawnow limitrate;
                pause(0.05);
            end
            
        catch ME
            fprintf('Error in ROI creation: %s\n', ME.message);
            if exist('h', 'var') && isvalid(h)
                delete(h);
            end
        end
    end
    
    % Helper function to handle double-click events
    function checkDoubleClick(h)
        persistent lastClick
        persistent clickCount
        
        if isempty(lastClick)
            lastClick = tic;
            clickCount = 1;
        else
            timeDiff = toc(lastClick);
            if timeDiff < 0.3  % 300ms threshold for double-click
                clickCount = clickCount + 1;
                if clickCount == 2
                    % Double-click detected, process the ROI
                    processROI(h);
                    % Reset click tracking
                    lastClick = [];
                    clickCount = 0;
                    % Reset the window button function
                    set(hFig, 'WindowButtonDownFcn', []);
                end
            else
                % Too much time has passed, reset counter
                clickCount = 1;
            end
            lastClick = tic;
        end
    end
    
    % Helper function to process ROI after completion
    function processROI(h)
        try
            if ~isvalid(h)
                return;
            end
            
            % Get the position
            pos = getPosition(h);
            
            % Check if ROI is valid (at least 3 points)
            if size(pos, 1) < 3
                delete(h);
                return;
            end
            
            % Create mask using poly2mask for robustness
            mask = poly2mask(pos(:,1), pos(:,2), size(brightFieldImage,1), size(brightFieldImage,2));
            [y, x] = find(mask);
            
            % Check if mask is valid
            if isempty(x) || isempty(y)
                delete(h);
                return;
            end
            
            roiCount = roiCount + 1;
            
            % Store coordinates in MRB
            MRB{roiCount,1} = [x, y];
            allMasks = allMasks | mask;
            
            % Draw boundary
            hold(hAx, 'on');
            roiLine = plot(hAx, pos(:,1), pos(:,2), 'g', 'LineWidth', 2);
            
            % Add label
            centroid = mean(pos);
            roiLabel = text(hAx, centroid(1), centroid(2), sprintf('ROI %d', roiCount), ...
                'Color', 'green', 'FontWeight', 'bold', 'FontSize', 12);
            
            % Create context menu for deletion
            cmenu = uicontextmenu(hFig);
            uimenu(cmenu, 'Label', 'Delete ROI', 'Callback', {@deleteROI, roiCount});
            roiLine.UIContextMenu = cmenu;
            
            % Store handles for later management
            roiHandles(roiCount).line = roiLine;
            roiHandles(roiCount).label = roiLabel;
            roiHandles(roiCount).mask = mask;
            roiHandles(roiCount).active = true;
            
            % Clean up
            delete(h);
            
        catch ME
            fprintf('Error processing ROI: %s\n', ME.message);
            if exist('h', 'var') && isvalid(h)
                delete(h);
            end
        end
    end
    
    % Helper function to update ROI position
    function updatePosition(h)
        % This function will be called when the ROI is being modified
        if isvalid(h)
            % You can add additional handling here if needed
            drawnow limitrate;
        end
    end
    
    % Callback for deleting ROI
    function deleteROI(~, ~, roiIdx)
        try
            if roiIdx <= length(roiHandles) && isfield(roiHandles, 'active') && roiHandles(roiIdx).active
                % Delete graphics objects if they're valid
                if isfield(roiHandles(roiIdx), 'line') && isvalid(roiHandles(roiIdx).line)
                    delete(roiHandles(roiIdx).line);
                end
                
                if isfield(roiHandles(roiIdx), 'label') && isvalid(roiHandles(roiIdx).label)
                    delete(roiHandles(roiIdx).label);
                end
                
                % Mark as inactive
                roiHandles(roiIdx).active = false;
                
                % Update combined mask
                allMasks = false(size(brightFieldImage));
                for i = 1:length(roiHandles)
                    if isfield(roiHandles(i), 'active') && roiHandles(i).active
                        allMasks = allMasks | roiHandles(i).mask;
                    end
                end
                
                % Update MRB cell array
                MRB{roiIdx} = [];
            end
        catch ME
            fprintf('Error deleting ROI: %s\n', ME.message);
        end
    end
    
    % Callback for clearing all ROIs
    function clearAllROIs(~,~)
        try
            % Delete all ROI graphics
            for i = 1:length(roiHandles)
                if isfield(roiHandles(i), 'active') && roiHandles(i).active
                    if isfield(roiHandles(i), 'line') && isvalid(roiHandles(i).line)
                        delete(roiHandles(i).line);
                    end
                    if isfield(roiHandles(i), 'label') && isvalid(roiHandles(i).label)
                        delete(roiHandles(i).label);
                    end
                    roiHandles(i).active = false;
                end
            end
            
            % Reset variables
            roiCount = 0;
            MRB = cell(1,1);
            allMasks = false(size(brightFieldImage));
            roiHandles = struct('line', {}, 'label', {}, 'mask', {}, 'active', {});
            
            % Refresh display
            cla(hAx);
            imshow(brightFieldImage, [], 'Parent', hAx);
            title('Draw ROIs for localization. Double-click to finish each ROI.');
        catch ME
            fprintf('Error clearing ROIs: %s\n', ME.message);
        end
    end
    
    % Callback for finishing ROI drawing
    function finishROIs(~,~)
        try
            % Count active ROIs
            activeCount = 0;
            for i = 1:length(roiHandles)
                if isfield(roiHandles(i), 'active') && roiHandles(i).active
                    activeCount = activeCount + 1;
                end
            end
            
            % Create new cell array for valid ROIs
            cleanMRB = cell(activeCount, 1);
            validIdx = 1;
            
            % Copy only active ROIs
            for i = 1:length(roiHandles)
                if isfield(roiHandles(i), 'active') && roiHandles(i).active && ~isempty(MRB{i})
                    cleanMRB{validIdx} = MRB{i};
                    validIdx = validIdx + 1;
                end
            end
            
            MRB = cleanMRB;
            
            % Close figure
            delete(hFig);
        catch ME
            fprintf('Error finishing ROIs: %s\n', ME.message);
            % Ensure figure is closed even if there's an error
            if isvalid(hFig)
                delete(hFig);
            end
        end
    end
    
    % Set figure close request function to use the Finish button
    set(hFig, 'CloseRequestFcn', @finishROIs);
    
    % Wait for user to finish by clicking the Finish button
    waitfor(hFig);
    
    % If no ROIs were drawn or all were deleted, create default ROI
    if isempty(MRB) || all(cellfun(@isempty, MRB))
        [height, width] = size(brightFieldImage);
        [X, Y] = meshgrid(1:width, 1:height);
        MRB = {[X(:), Y(:)]};
        allMasks = true(size(brightFieldImage)); % Full image mask
        fprintf('No valid ROIs remain. Created default ROI covering the entire image.\n');
    else
        fprintf('Created %d ROIs for localization.\n', length(MRB));
    end
end

% Function to handle blinking analysis with different methods
function [blinkResults, performNewAnalysis] = handleBlinkingAnalysis(imageData, locResults, params, outputDir, progressLogFile)
    % Create blinking directory if it doesn't exist
    if ~exist(params.blinkingDir, 'dir')
        mkdir(params.blinkingDir);
        updateProgress(progressLogFile, 'Directory', 'Created blinking analysis directory');
    end
    
    % Set frame range for analysis
    if ~isempty(imageData) && numel(imageData) > 0
        % Automatically select all frames for analysis if we have image data
        startFrame = 1;
        endFrame = size(imageData, 3);
        params.analysisStartFrame = startFrame;
        params.analysisEndFrame = endFrame;
        params.analysisFrameRange = startFrame:endFrame;
    else
        params.analysisStartFrame = 1;
        if isfield(params, 'totalFrames') && ~isempty(params.totalFrames)
            params.analysisEndFrame = params.totalFrames;
        end
        params.analysisFrameRange = params.analysisStartFrame:params.analysisEndFrame;
    end
    
    updateProgress(progressLogFile, 'Frame Selection', ...
        sprintf('Using frame range %d to %d for bleach curve analysis', ...
        params.analysisStartFrame, params.analysisEndFrame));
    
    % Try to load existing blinking results
    [blinkResults, loadedFromFile] = loadDataStage(params.outputBaseDir, ...
        params.outputFolderName, '3_BlinkingAnalysis', 'blinking_results', params);
    performNewAnalysis = ~loadedFromFile;
    
    if performNewAnalysis
        updateProgress(progressLogFile, 'Blinking Analysis', 'Performing new analysis...');
        
        % Ensure we have localization data
        if isempty(locResults) || ~isfield(locResults, 'combinedPoints') || isempty(locResults.combinedPoints)
            % Try to load localization data
            [locResults, locLoaded] = loadDataStage(params.outputBaseDir, ...
                params.outputFolderName, '2_Localization', 'localization_results', params);
            
            if ~locLoaded || isempty(locResults) || ~isfield(locResults, 'combinedPoints') || isempty(locResults.combinedPoints)
                error('No localization data found. Cannot perform blinking analysis.');
            end
        end
        
        if isempty(imageData) || numel(imageData) == 0
            % Try to load image data
            [imageData, imgLoaded] = loadDataStage(params.outputBaseDir, ...
                params.outputFolderName, 'ImageData', 'image_data', params);
            
            if ~imgLoaded || isempty(imageData) || numel(imageData) == 0
                error('No image data found. Cannot perform image-based blinking analysis.');
            end
        end

        % Choose analysis method based on user selection
        try

            blinkResults = analyzeBlinkingFromImage(locResults, params, imageData);
            % Save results with consistent naming
            blinkResultsFile = fullfile(outputDir, '3_BlinkingAnalysis', ...
                ['blinking_results_', params.timestamp, '.mat']);
            save(blinkResultsFile, 'blinkResults', '-v7.3');
            
            % Also save a standard version for compatibility
            standardResultsFile = fullfile(outputDir, '3_BlinkingAnalysis', 'blinking_results.mat');
            save(standardResultsFile, 'blinkResults', '-v7.3');
            
            updateProgress(progressLogFile, 'Blinking Analysis', 'Saved blinking results');
        catch ME
            handleErrorWithFallback(progressLogFile, ME, 'Error in blinking analysis', ...
                'Analysis failed, returning empty results');
            blinkResults = struct();
        end      
    end
end

function datasetPaths = findDatasets(baseDir)
    % Initialize cell array for dataset paths
    datasetPaths = {};
    
    % Get all items in the directory
    items = dir(baseDir);
    
    % Process each item
    for i = 1:length(items)
        if strcmp(items(i).name, '.') || strcmp(items(i).name, '..')
            continue;
        end
        
        fullPath = fullfile(items(i).folder, items(i).name);
        
        if items(i).isdir
            % Check if this directory contains analysis results
            if contains(items(i).name, '_Run_') || ...
               contains(items(i).name, '_CommonData')
                continue; % Skip analysis results directories
            end
            
            % Recursively search subdirectories
            subPaths = findDatasets(fullPath);
            if ~isempty(subPaths)
                if ~iscell(subPaths)
                    subPaths = {subPaths};
                end
                datasetPaths = [datasetPaths(:); subPaths(:)];
            end
        else
            % Check for ND2 files
            [~, ~, ext] = fileparts(items(i).name);
            if strcmpi(ext, '.nd2')
                datasetPaths = [datasetPaths(:); {fullPath}];
            end
        end
    end
    
    % Look for TIFF sequences in current directory
    tiffFiles = dir(fullfile(baseDir, '*.tif'));
    tiffFiles = [tiffFiles; dir(fullfile(baseDir, '*.tiff'))];
    
    if ~isempty(tiffFiles)
        % If there are TIFF files, add the directory as a dataset
        datasetPaths = [datasetPaths(:); {baseDir}];
    end
    
    % Ensure output is a column cell array
    datasetPaths = datasetPaths(:);
end

% Function to find the most recent result file
function resultFile = findResultFile(outputBaseDir, datasetName, stagePath, dataName, params)
    resultFile = '';
    
    % Check if we should automatically use the latest file
    autoUseLatest = isfield(params, 'autoUseLatest') && params.autoUseLatest;
    
    % Define possible locations based on dataset structure
    possibleLocations = {};
    
    % 1. First check in the dataset-specific results directory
    if isfield(params, 'outputDir') && ~isempty(params.outputDir)
        possibleLocations{end+1} = fullfile(params.outputDir, stagePath);
    end
    
    % 2. Check in the common data directory for this dataset
    possibleLocations{end+1} = fullfile(outputBaseDir, [datasetName, '_CommonData'], stagePath);
    
    % 3. For image data, also check directly in the ImageData folder
    if strcmp(dataName, 'image_data') || strcmp(stagePath, 'ImageData')
        possibleLocations{end+1} = fullfile(outputBaseDir, [datasetName, '_CommonData'], 'ImageData');
    end
    
    % 4. For localization data, check in various analysis result folders
    if strcmp(dataName, 'localization_results') || strcmp(stagePath, '2_Localization')
        % Find all analysis result folders for this dataset
        analysisPattern = fullfile(outputBaseDir, [datasetName, '_*_Run_*']);
        analysisDirs = dir(analysisPattern);
        for i = 1:length(analysisDirs)
            if analysisDirs(i).isdir
                possibleLocations{end+1} = fullfile(analysisDirs(i).folder, analysisDirs(i).name, '2_Localization');
            end
        end
    end
    
    % Search all possible locations
    for i = 1:length(possibleLocations)
        location = possibleLocations{i};
        
        % Look for files matching the pattern
        filePattern = fullfile(location, [dataName, '*.mat']);
        files = dir(filePattern);
        
        % Also try without the dataName prefix (for standard files)
        if isempty(files) && ~strcmp(dataName, 'image_data')
            filePattern = fullfile(location, '*.mat');
            files = dir(filePattern);
        end
        
        if ~isempty(files)
            if autoUseLatest
                % Sort by date and use the most recent
                [~, idx] = sort([files.datenum], 'descend');
                resultFile = fullfile(files(idx(1)).folder, files(idx(1)).name);
                return; % Found a file, return immediately
            else
                % Let user select from available files
                fileNames = {files.name};
                [selection, ok] = listdlg('ListString', fileNames, ...
                    'SelectionMode', 'single', ...
                    'Name', 'Select Result File', ...
                    'PromptString', ['Select ', dataName, ' file:'], ...
                    'ListSize', [300 200]);
                
                if ok && ~isempty(selection)
                    resultFile = fullfile(files(selection).folder, files(selection).name);
                    return; % Found a file, return immediately
                end
            end
        end
    end
    
    % If we get here, no file was found in any location
    resultFile = '';
end

% Unified data loading function that handles all cases
% Unified data loading function that handles all cases
function [data, loadedFromFile] = loadDataStage(outputBaseDir, datasetName, stagePath, dataName, params)
    % Initialize
    data = [];
    loadedFromFile = false;
    
    % Log the search
    if isfield(params, 'progressLogFile')
        updateProgress(params.progressLogFile, 'Data Search', ...
            sprintf('Searching for %s in %s for dataset %s', dataName, stagePath, datasetName));
    end
    
    % Define possible locations based on dataset structure
    possibleLocations = {};
    
    % Get the dataset parent folder path (for ND2 files or TIFF sequences)
    datasetParentPath = '';
    if isfield(params, 'datasetParentPath') && ~isempty(params.datasetParentPath)
        datasetParentPath = params.datasetParentPath;
    elseif isfield(params, 'currentDataPath') && ~isempty(params.currentDataPath)
        datasetParentPath = fileparts(params.currentDataPath);
    end
    
    % Extract SR number from dataset name for matching
    srNumber = '';
    if isfield(params, 'rawDatasetName') && ~isempty(params.rawDatasetName)
        srMatch = regexp(params.rawDatasetName, 'SR(\d+)', 'tokens', 'once');
        if ~isempty(srMatch)
            srNumber = srMatch{1};
        end
    elseif ~isempty(datasetName)
        srMatch = regexp(datasetName, 'SR(\d+)', 'tokens', 'once');
        if ~isempty(srMatch)
            srNumber = srMatch{1};
        end
    end

    % Check if we're looking for image data, localization data, or MRB data
    isImageDataSearch = strcmp(stagePath, 'ImageData') || strcmp(dataName, 'image_data');
    isLocalizationSearch = strcmp(stagePath, '2_Localization') || strcmp(dataName, 'localization_results');
    isMRBSearch = strcmp(dataName, 'MRB');
    
    % Check for search tags
    hasImageDataTag = isfield(params, 'imageDataSearchTag') && strcmp(params.imageDataSearchTag, 'CommonData');
    hasLocalizationTag = isfield(params, 'localizationSearchTag') && strcmp(params.localizationSearchTag, 'LOC');
    
    % for single dataset
    % 1. First check in the dataset-specific results directory
    if isfield(params, 'outputDir') && ~isempty(params.outputDir)
        if (isImageDataSearch && hasImageDataTag && contains(params.outputDir, 'CommonData')) || ...
           (isLocalizationSearch && hasLocalizationTag && contains(params.outputDir, 'LOC')) || ...
           isMRBSearch
            possibleLocations{end+1} = fullfile(params.outputDir, stagePath);
        end
    end
    
    % for multiple dataset
    % 2. Check in the common data directory for this dataset (only for image data)
    if isImageDataSearch && hasImageDataTag
        possibleLocations{end+1} = fullfile(outputBaseDir, [datasetName, '_CommonData'], stagePath);
    end
    
    % 3. For localization data or MRB data, search in all LOC folders in the Analysis_Results directory
    if (isLocalizationSearch && hasLocalizationTag) || isMRBSearch
        % Get the raw dataset name without extension for better matching
        rawDatasetName = '';
        if isfield(params, 'rawDatasetName') && ~isempty(params.rawDatasetName)
            rawDatasetName = params.rawDatasetName;
        end
        
        % Search in parent directory's Analysis_Results folder
        if ~isempty(datasetParentPath)
            parentAnalysisDir = fullfile(datasetParentPath, 'Analysis_Results');
            
            % If parent Analysis_Results exists, search for LOC folders
            if exist(parentAnalysisDir, 'dir')
                % Look for folders with LOC tag
                locFolders = dir(fullfile(parentAnalysisDir, [datasetName, '_LOC_*']));
                
                if ~isempty(locFolders)
                    % Sort by date to get the latest
                    [~, idx] = sort([locFolders.datenum], 'descend');
                    
                    % Add each LOC folder to possible locations
                    for i = 1:length(idx)
                        locFolder = fullfile(locFolders(idx(i)).folder, locFolders(idx(i)).name);
                        possibleLocations{end+1} = fullfile(locFolder, stagePath);
                    end
                end
            end
        end
    end
    
    % Search all possible locations
    for i = 1:length(possibleLocations)
        location = possibleLocations{i};
        
        % Skip if location doesn't exist
        if ~exist(location, 'dir')
            continue;
        end
        
        % Look for files matching the pattern
        if isMRBSearch
            % For MRB, look for both timestamped and standard files
            filePattern = fullfile(location, 'MRB*.mat');
        else
            filePattern = fullfile(location, [dataName, '*.mat']);
        end
        
        files = dir(filePattern);
        
        % Also try without the dataName prefix (for standard files)
        if isempty(files) && ~isImageDataSearch
            filePattern = fullfile(location, '*.mat');
            files = dir(filePattern);
        end
        
        if ~isempty(files)
            % Check if we should automatically use the latest file
            autoUseLatest = isfield(params, 'autoUseLatest') && params.autoUseLatest;
            
            if autoUseLatest
                % Sort by date and use the most recent
                [~, idx] = sort([files.datenum], 'descend');
                resultFile = fullfile(files(idx(1)).folder, files(idx(1)).name);
                
                try
                    loadedData = load(resultFile);
                    
                    % Handle different data types
                    if isMRBSearch
                        if isfield(loadedData, 'MRB')
                            data = loadedData.MRB;
                            loadedFromFile = true;
                            
                            if isfield(params, 'progressLogFile')
                                updateProgress(params.progressLogFile, 'Data Loading', ...
                                    sprintf('Successfully loaded MRB from %s', resultFile));
                            end
                            return;
                        end
                    elseif isfield(loadedData, dataName)
                        data = loadedData.(dataName);
                        loadedFromFile = true;
                        
                        if isfield(params, 'progressLogFile')
                            updateProgress(params.progressLogFile, 'Data Loading', ...
                                sprintf('Successfully loaded %s from %s', dataName, resultFile));
                        end
                        return;
                    elseif isImageDataSearch && isfield(loadedData, 'imageData')
                        data = loadedData.imageData;
                        loadedFromFile = true;
                        
                        if isfield(params, 'progressLogFile')
                            updateProgress(params.progressLogFile, 'Data Loading', ...
                                sprintf('Successfully loaded image data from %s', resultFile));
                        end
                        return;
                    end
                catch ME
                    if isfield(params, 'progressLogFile')
                        updateProgress(params.progressLogFile, 'Warning', ...
                            sprintf('Error loading from %s: %s', resultFile, ME.message));
                    end
                end
            else
                % Let user select from available files
                fileNames = {files.name};
                [selection, ok] = listdlg('ListString', fileNames, ...
                    'SelectionMode', 'single', ...
                    'Name', 'Select Result File', ...
                    'PromptString', ['Select ', dataName, ' file:'], ...
                    'ListSize', [300 200]);
                
                if ok && ~isempty(selection)
                    resultFile = fullfile(files(selection).folder, files(selection).name);
                    
                    try
                        loadedData = load(resultFile);
                        
                        % Handle different data types
                        if isMRBSearch
                            if isfield(loadedData, 'MRB')
                                data = loadedData.MRB;
                                loadedFromFile = true;
                                
                                if isfield(params, 'progressLogFile')
                                    updateProgress(params.progressLogFile, 'Data Loading', ...
                                        sprintf('Successfully loaded MRB from %s', resultFile));
                                end
                                return;
                            end
                        elseif isfield(loadedData, dataName)
                            data = loadedData.(dataName);
                            loadedFromFile = true;
                            
                            if isfield(params, 'progressLogFile')
                                updateProgress(params.progressLogFile, 'Data Loading', ...
                                    sprintf('Successfully loaded %s from %s', dataName, resultFile));
                            end
                            return;
                        elseif isImageDataSearch && isfield(loadedData, 'imageData')
                            data = loadedData.imageData;
                            loadedFromFile = true;
                            
                            if isfield(params, 'progressLogFile')
                                updateProgress(params.progressLogFile, 'Data Loading', ...
                                    sprintf('Successfully loaded image data from %s', resultFile));
                            end
                            return;
                        end
                    catch ME
                        if isfield(params, 'progressLogFile')
                            updateProgress(params.progressLogFile, 'Warning', ...
                                sprintf('Error loading from %s: %s', resultFile, ME.message));
                        end
                    end
                end
            end
        end
    end
    
    % If we get here, no file was found or loaded in any location
    if isfield(params, 'progressLogFile')
        updateProgress(params.progressLogFile, 'Data Loading', ...
            sprintf('No %s found in any location', dataName));
    end
    
    data = [];
    loadedFromFile = false;
end

function handleErrorWithFallback(progressLogFile, funcHandle, step, errorMessage, fallbackMessage)
    try
        % Execute the function
        funcHandle();
    catch ME
        % Log the error
        updateProgress(progressLogFile, step, [errorMessage ': ' ME.message]);
        
        % Display error in console
        fprintf('ERROR: %s: %s\n', errorMessage, ME.message);
        
        % Log fallback action if provided
        if nargin > 4 && ~isempty(fallbackMessage)
            updateProgress(progressLogFile, step, fallbackMessage);
            fprintf('RECOVERY: %s\n', fallbackMessage);
        end
        
        % Re-throw the error if needed
        % rethrow(ME);
    end
end
