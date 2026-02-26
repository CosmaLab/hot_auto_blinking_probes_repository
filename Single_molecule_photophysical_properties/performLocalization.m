function results = performLocalization(imageData, params, frameRange)
    % Check if frameRange is provided and valid
    if nargin < 3 || isempty(frameRange)
        frameRange = 1:size(imageData, 3);
    end
    
    % Ensure frameRange is within valid bounds
    frameRange = frameRange(frameRange >= 1 & frameRange <= size(imageData, 3));
    
    % Check if we have any valid frames to process
    if isempty(frameRange)
        error('No valid frames to process for localization');
    end
    
    % Extract parameters
    threshold = params.threshold;
    mask_sz = params.mask_sz;
    IM_info = params.IM_info;
    CAM_set = [params.camera.gain, params.camera.QE];  % Create camera params array
    EM_wave_FP = params.EM_wave_FP;
    dil_sz = params.dil_sz;
    threshold_distance = params.threshold_distance;

    % Handle optional frameRange parameter
    if nargin < 3 || isempty(frameRange)
        % Default: use all available frames
        frameRange = 1:size(imageData, 3);
        fprintf('No frame range specified, using all %d frames\n', length(frameRange));
    else
        % Ensure frameRange is valid - but don't modify user's selection based on frames_to_skip
        % Just check if it's within bounds of the image data
        if min(frameRange) < 1 || max(frameRange) > size(imageData, 3)
            warning('Frame range contains invalid frames. Adjusting to valid range.');
            frameRange = frameRange(frameRange >= 1 & frameRange <= size(imageData, 3));
            if isempty(frameRange)
                error('No valid frames in the specified range');
            end
        end
        fprintf('Using user-specified frame range: %d to %d (%d frames total)\n', ...
            min(frameRange), max(frameRange), length(frameRange));
    end
    
    % Create a unique filename based on the frame range
    frameRangeStr = sprintf('%d-%d', min(frameRange), max(frameRange));
    matFilePath = fullfile(params.outputDir, sprintf('selected_frames_range_%s.mat', frameRangeStr));
    
    % Create directory if it doesn't exist
    if ~exist(fileparts(matFilePath), 'dir')
        mkdir(fileparts(matFilePath));
    end
    
    % Save the selected frames and parameters with metadata
    selectedFrames = imageData(:,:,frameRange);
    frameMetadata = struct();
    frameMetadata.frameRange = frameRange;
    frameMetadata.totalFrames = size(imageData, 3);
    frameMetadata.analysisDate = datestr(now);
    frameMetadata.frameRangeDescription = frameRangeStr;
    
    save(matFilePath, 'selectedFrames', 'frameRange', 'params', 'frameMetadata', '-v7.3');
    fprintf('Saved selected frames (range %s) to: %s\n', frameRangeStr, matFilePath);
    
    % Initialize results structure
    results = struct();
    
    % Get number of frames to analyze - use all frames in the range
    numFramesToAnalyze = length(frameRange);
    fprintf('Analyzing %d frames in range %s\n', numFramesToAnalyze, frameRangeStr);
    
    % Create subset of frames for analysis
    imageSubset = imageData(:,:,frameRange);
    
    % Determine if we should show visualizations based on frame count
    showfit = 1;  % Default to showing visualizations
    % Get MRB from params if available
    MRB = {};
    if isfield(params, 'MRB')
        MRB = params.MRB;
    end
    
    % Ensure we have a ROI mask - but don't recreate if it already exists
    if ~isfield(params, 'roiMask') || isempty(params.roiMask)
        fprintf('No ROI mask found in params. ');
        
        if ~isempty(MRB)
            fprintf('Creating ROI mask from MRB data...\n');
            % Create ROI mask from MRB
            roiMask = false(size(imageData, 1), size(imageData, 2));
            
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
            
            % Add to params for use in localization
            params.roiMask = roiMask;
            
            % Visualize the ROI mask for debugging
            fig = figure('Name', 'ROI Mask Visualization', 'Position', [100, 100, 800, 600]);
            imagesc(roiMask);
            colormap('gray');
            title('ROI Mask');
            saveas(fig, fullfile(params.outputDir, 'roi_mask.png'));
            close(fig);
        else
            fprintf('No MRB data available. Using entire image as ROI.\n');
            params.roiMask = true(size(imageData, 1), size(imageData, 2));
        end
    else
        fprintf('Using ROI mask already provided in params (covers %.2f%% of image)...\n', ...
            sum(params.roiMask(:))/numel(params.roiMask)*100);
    end

    % Call localization with camera parameters as array and explicitly pass numFramesToAnalyze
    RCA = [0,0,0,0,0,0];
    try
        % Pass the imageSubset directly (already contains only the selected frames)
        [~, ~, ~, XYZ, PC] = iQ_posfit3frames_forDyeLocalization_file_input(imageSubset, threshold, mask_sz, IM_info, CAM_set, EM_wave_FP, MRB, dil_sz, RCA, [], showfit,params);
        %[XYZ, ~] = iQ_localization_irregular_ROI(imageSubset, threshold, mask_sz, IM_info, CAM_set, EM_wave_FP, MRB, dil_sz, params)
        fprintf('Localization completed successfully for frames %s.\n', frameRangeStr);
    catch ME
        fprintf('Error during localization: %s\n', ME.message);
        % Create a more detailed error report
        errorReport = getReport(ME, 'extended');
        errorLogFile = fullfile(params.outputDir, 'localization_error.log');
        fid = fopen(errorLogFile, 'w');
        fprintf(fid, 'Localization Error Report\n');
        fprintf(fid, 'Date: %s\n', datestr(now));
        fprintf(fid, 'Frame Range: %s\n', frameRangeStr);
        fprintf(fid, 'Error Message: %s\n\n', ME.message);
        fprintf(fid, 'Detailed Error Report:\n%s\n', errorReport);
        fclose(fid);
        fprintf('Detailed error report saved to: %s\n', errorLogFile);

        % Return empty results with error information
        results.error = ME.message;
        results.errorDetails = errorReport;
        results.frameRange = frameRange;
        results.frameRangeStr = frameRangeStr;
        return;
    end
    
    % Store results
    results.XYZ = XYZ;
    results.PC = PC;

    % % Combine all ROI points into one array
    % allROIPoints = [];
    % for i = 1:length(PC)
    %     if ~isempty(PC{i})
    %         allROIPoints = [allROIPoints; PC{i}];
    %     end
    % end

    % % Use XYZ directly instead of recombining PC
    % % XYZ already contains all points within ROIs
    allROIPoints = XYZ;

    % Print initial points information
    fprintf('\nInitial Points Statistics:\n');
    fprintf('Total points detected across all ROIs: %d\n', size(allROIPoints, 1));
    
    % % Print initial points information
    % fprintf('\nInitial Points Statistics:\n');
    % fprintf('Total points detected: %d\n', size(XYZ, 1));

    % Filter localization results
    FIL = [ 1,    1,      0,      0,      0,      0,        0,      0;...
            80,   80,     0,      0,      0,      10000,   0,      600;...
            350,  350,    0.3,    60,     60,     1000000, 90000,  10000];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TOP of the routine
    FIL_sgmx=FIL(1,1);      FIL_sgmy=FIL(1,2);    FIL_ellip=FIL(1,3);
    FIL_errx=FIL(1,4);      FIL_erry=FIL(1,5);    FIL_int=FIL(1,6);
    FIL_sgmxy=FIL(1,7);     FIL_amp=FIL(1,8);
    %% filter the data
    % flt_new=XYZ;
    flt_new = allROIPoints;  % Use XYZ directly

    % sigmax filter
    if FIL_sgmx == 1; MIN_sgmx=FIL(2,1); MAX_sgmx=FIL(3,1);
        sgmx_raw=flt_new;
        sgmxf= sgmx_raw(:,8)>MIN_sgmx & sgmx_raw(:,8)<MAX_sgmx;
        sgmx_new=sgmx_raw(sgmxf,:);
        flt_new=sgmx_new;
    end
    % sigmay filter
    if FIL_sgmy == 1; MIN_sgmy=FIL(2,2); MAX_sgmy=FIL(3,2);
        sgmy_raw=flt_new;
        sgmyf= sgmy_raw(:,9)>MIN_sgmy & sgmy_raw(:,9)<MAX_sgmy;
        sgmy_new=sgmy_raw(sgmyf,:);
        flt_new=sgmy_new;
    end
    % Elli filter
    if FIL_ellip == 1; MIN_ellip=FIL(2,3); MAX_ellip=FIL(3,3);
        ell_raw=flt_new;
        ell_f=ell_raw(:,10)>MIN_ellip &ell_raw(:,10)<MAX_ellip ;
        ell_new=ell_raw(ell_f,:);
        flt_new=ell_new;
    end
    % Errx filter
    if FIL_errx == 1; MIN_errx=FIL(2,4); MAX_errx=FIL(3,4);
        errx_raw=flt_new;
        errxf= errx_raw(:,17)>MIN_errx & errx_raw(:,17)<MAX_errx;
        errx_new=errx_raw(errxf,:);
        flt_new=errx_new;
    end
    % Erry filter
    if FIL_erry == 1; MIN_erry=FIL(2,5); MAX_erry=FIL(3,5);
        erry_raw=flt_new;
        erryf= erry_raw(:,18)>MIN_erry & erry_raw(:,18)<MAX_erry;
        erry_new=erry_raw(erryf,:);
        flt_new=erry_new;
    end
    % Int filter
    if FIL_int == 1; MIN_int=FIL(2,6); MAX_int=FIL(3,6);
        int_raw=flt_new;
        int_f=int_raw(:,15)>MIN_int &int_raw(:,15)<MAX_int;
        int_new=int_raw(int_f,:);
        flt_new=int_new;
    end
    
    % sgmxy filter
    if FIL_sgmxy == 1; MIN_sgmxy=FIL(2,7); MAX_sgmxy=FIL(3,7);
        sgmxy_raw=flt_new(:,8).*flt_new(:,9);
        sgmxy_f=sgmxy_raw >MIN_sgmxy & sgmxy_raw < MAX_sgmxy;
        sgmxy_new=flt_new(sgmxy_f,:);
        flt_new=sgmxy_new;
    end
    
    % amp filter
    if FIL_amp == 1; MIN_amp=FIL(2,8); MAX_amp=FIL(3,8);
        amp_raw=flt_new(:,13);
        amp_f=amp_raw >MIN_amp & amp_raw < MAX_amp;
        amp_new=flt_new(amp_f,:);
        flt_new=amp_new;
    end

    PC_new_flt=flt_new;
    
    % Extract points from each frame - initialize with first frame
    frameIndices = unique(PC_new_flt(:, 12));
    
    % Check if we have any points
    if isempty(frameIndices)
        fprintf('Warning: No points detected in any frame after filtering.\n');
        combinedPoints = [];
    else
        % Start with points from the first frame
        firstFrameIdx = frameIndices(1);
        combinedPoints = PC_new_flt(PC_new_flt(:, 12) == firstFrameIdx, :);
        
        % Process all remaining frames
        for frameIdx = 2:length(frameIndices)
            currentFramePoints = PC_new_flt(PC_new_flt(:, 12) == frameIndices(frameIdx), :);
            
            % Add unique points from this frame
            for j = 1:size(currentFramePoints, 1)
                distances = sqrt(sum((combinedPoints(:,1:2) - currentFramePoints(j, 1:2)).^2, 2)); 
                if isempty(find(distances < threshold_distance, 1))
                    combinedPoints = [combinedPoints; currentFramePoints(j, :)]; 
                end
            end
        end
    end
    
    fprintf('\nFinal Statistics:\n');
    fprintf('Total filtered points: %d\n', size(PC_new_flt, 1));
    fprintf('Unique points after combining frames: %d\n', size(combinedPoints, 1));
    fprintf('Points removed by distance threshold: %d\n', ...
        size(PC_new_flt, 1) - size(combinedPoints, 1));
    
    % Store results with more descriptive variable names
    % results.allLocalizedPoints = XYZ;           % Raw localization results
    results.allLocalizedPoints = allROIPoints;    % Raw localization results from all ROIs
    results.allFilteredPoints = PC_new_flt;     % Points after filtering
    results.combinedPoints = combinedPoints;    % Unique points after combining frames
    
    results.numFramesAnalyzed = numFramesToAnalyze;
    results.frameRange = frameRange;  % Store the frame range used
    results.frameRangeStr = frameRangeStr;
    
    % Create a directory for sample frame visualizations
    sampleFramesDir = fullfile(params.outputDir, sprintf('sample_frames_%s', params.runID));
    if ~exist(sampleFramesDir, 'dir')
        mkdir(sampleFramesDir);
        fprintf('Created directory for sample frames: %s\n', sampleFramesDir);
    end
    
%% Visualize and save sample frames with localization results
    numSampleFrames = min(5, numFramesToAnalyze); % Limit to 5 sample frames
    sampleFrameIndices = round(linspace(1, numFramesToAnalyze, numSampleFrames));
    
    for i = 1:length(sampleFrameIndices)
        frameIdx = sampleFrameIndices(i);
        actualFrameIdx = frameRange(frameIdx);
        
        % Create figure for this frame
        h = figure('Name', sprintf('Frame %d Localization', actualFrameIdx), 'Position', [100, 100, 800, 600], 'Visible', 'off');
        
        % Display the image
        imagesc(imageData(:,:,actualFrameIdx));
        colormap('gray');
        axis image;
        hold on;
        
        % Get points detected in this frame (adjust frame index for PC_new_flt)
        framePoints = PC_new_flt(PC_new_flt(:, 12) == frameIdx, :);
        % Plot detected points
        if ~isempty(framePoints)            
            % Draw boxes around detected points to show mask size
            half_box = floor(mask_sz/2);
            for j = 1:size(framePoints, 1)
                x = framePoints(j, 1);
                y = framePoints(j, 2);
                rectangle('Position', [x-half_box, y-half_box, mask_sz, mask_sz], ...
                          'EdgeColor', 'g', 'LineWidth', 0.1, 'LineStyle', ':');
            end
        end
        % Add title and labels
        title(sprintf('Frame %d - Detected Points: %d', actualFrameIdx, size(framePoints, 1)));
        xlabel('X (pixels)');
        ylabel('Y (pixels)');
        
        % Add colorbar
        colorbar;
        
        % Save the figure as PNG
        framePngFile = fullfile(sampleFramesDir, sprintf('frame_%d_localization.png', actualFrameIdx));
        saveas(h, framePngFile);
        
        % Also save as fig for later editing
        frameFigFile = fullfile(sampleFramesDir, sprintf('frame_%d_localization.fig', actualFrameIdx));
        saveas(h, frameFigFile);
        
        % Close the figure
        close(h);
        
        fprintf('Saved visualization for frame %d\n', actualFrameIdx);
    end
    
    % Create a summary visualization of all detected points
    h = figure('Name', 'All Detected Points', 'Position', [100, 100, 800, 800], 'Visible', 'off');
    
    % Display the maximum intensity projection if available
    if size(imageData, 3) > 1
        maxIntensityProj = max(imageData(:,:,frameRange), [], 3);
        imagesc(maxIntensityProj);
    else
        imagesc(imageData(:,:,frameRange(1)));
    end
    colormap('gray');
    axis image;
    hold on;
    
    % Overlay ROI boundaries if available
    if isfield(params, 'MRB') && ~isempty(params.MRB)
        for r = 1:length(params.MRB)  % Add the missing loop here
            % Get ROI coordinates
            x = params.MRB{r}(:,1);
            y = params.MRB{r}(:,2);
            
            % Create the same mask that's used in the localization function
            % This exactly matches the processing in iQ_posfit3frames_forDyeLocalization_file_input.m
            grid = false(size(imageData, 1), size(imageData, 2));
            grid(sub2ind(size(grid), y, x)) = 1;
            
            % Fill holes and dilate - exactly as done in the localization function
            grid1 = imfill(grid, 'holes');
            se1 = strel('disk', dil_sz);  % Use the same dilation size parameter
            grid2 = imdilate(grid1, se1);
            
            % Get the boundary of the processed mask
            B = bwboundaries(grid2);
            
            % Plot each boundary segment
            for k = 1:length(B)
                boundary = B{k};
                plot(boundary(:,2), boundary(:,1), 'g-', 'LineWidth', 2);
            end
            
            % Add ROI label at centroid of the processed mask
            [row, col] = find(grid2);
            if ~isempty(row) && ~isempty(col)
                centroid = [mean(col), mean(row)];
                text(centroid(1), centroid(2), sprintf('ROI %d', r), ...
                    'Color', 'green', 'FontWeight', 'bold', 'FontSize', 12);
            end
        end
    end

    % Draw boxes for all filtered and unique points to show mask size
    half_box = floor(mask_sz/2);
    numBoxesToShow = size(combinedPoints, 1);  % Show all points
    if numBoxesToShow > 0
        boxIndices = 1:numBoxesToShow;  % Use all indices

        for j = 1:length(boxIndices)
            idx = boxIndices(j);
            x = combinedPoints(idx, 1);
            y = combinedPoints(idx, 2);
            
            % Draw box with size based on mask_sz
            rectangle('Position', [x-half_box, y-half_box, mask_sz, mask_sz], ...
                      'EdgeColor', 'g', 'LineWidth', 0.1, 'LineStyle', ':');
        end
    end
    
    % Add title and labels
    title(sprintf('All Detected Points (%d) - Mask Size: %d px', size(combinedPoints, 1), mask_sz));
    
    % Add colorbar
    colorbar;
    
    % Save the figure as PNG
    summaryPngFile = fullfile(sampleFramesDir, 'all_detected_points.png');
    saveas(h, summaryPngFile);
    
    % Also save as fig for later editing (with a different variable name)
    summaryFigFile = fullfile(sampleFramesDir, 'all_detected_points.fig');
    saveas(h, summaryFigFile);
    
    % Close the figure
    close(h);
    
    fprintf('Saved summary visualization of all detected points\n');

    % Get output directory from params or use current directory as fallback
    if isfield(params, 'outputDir')
        outputDir = params.outputDir;
    else
        outputDir = pwd;
        warning('No output directory specified in params. Using current directory: %s', outputDir);
    end
    
    % Return results without saving - let the caller handle saving
    fprintf('Localization completed. Found %d unique points.\n', size(results.combinedPoints, 1));
end

