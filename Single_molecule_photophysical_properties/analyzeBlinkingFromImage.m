function results = analyzeBlinkingFromImage(locResults, params, imageData)
    % ANALYZEBLINKINGFROMIMAGE Analyze blinking behavior using raw image data
    %   RESULTS = ANALYZEBLINKINGFROMIMAGE(LOCRESULTS, PARAMS, IMAGEDATA) analyzes 
    %   single-molecule blinking behavior from raw image data and generates
    %   comprehensive statistics and visualizations.
    
    % Initialize results structure
    results = struct();
    results.analysisTimestamp = datestr(now, 'yyyymmdd_HHMMSS');
    results.analysisParams = params;
    
    % Check inputs
    if ~isfield(locResults, 'combinedPoints')
        error('Localization results must contain combinedPoints');
    end
    
    % Get parameters
    halfsquare = params.halfsquare;
    signalChange = params.signalChange;
    points = locResults.combinedPoints;
    
    % Get total number of frames
    totalFrames = size(imageData, 3);
    
    % Initialize arrays
    numMolecules = size(points, 1);
    integrated_roi_intensity = zeros(numMolecules, totalFrames);
    mean_bkg_int = zeros(numMolecules, 1);
    std_dev_bkg_int = zeros(numMolecules, 1);
    
    % Process each molecule
    fprintf('Extracting intensity traces for %d molecules...\n', numMolecules);
    
   % Store the full point data for reference
    results.fullPointData = points;
    
    % Pre-allocate ROI coordinates for all molecules
    x_coords = round(points(:, 1));
    y_coords = round(points(:, 2));

    for n = 1:numMolecules
        % Pre-allocate for this molecule
        molecule_intensity = zeros(1, totalFrames);
        
        % Extract intensity trace
        for j = 1:totalFrames
            roi_x_min = max(1, x_coords(n) - halfsquare);
            roi_x_max = min(size(imageData, 2), x_coords(n) + halfsquare);
            roi_y_min = max(1, y_coords(n) - halfsquare);
            roi_y_max = min(size(imageData, 1), y_coords(n) + halfsquare);
            
            box_region = imageData(roi_y_min:roi_y_max, roi_x_min:roi_x_max, j);
            molecule_intensity(j) = sum(box_region(:));
        end
        % Store intensity trace
        integrated_roi_intensity(n,:) = molecule_intensity;
        
        % Calculate background statistics
        % take the average of the lowest 95% of values as estimate of background
        intensityTrace = integrated_roi_intensity(n, :);
        sortedIntensity = sort(intensityTrace);
        numBgFrames = max(3, round(0.95 * totalFrames));
        bgValues = sortedIntensity(1:numBgFrames);
        
        mean_bkg_int(n) = mean(bgValues);
        std_dev_bkg_int(n) = std(bgValues);
        
        if mod(n, 100) == 0
            fprintf('Processed %d/%d molecules\n', n, numMolecules);
        end
    end
    
    %% Visualize sampled molecules' intensity traces
    % Create visualization directory
    if isfield(params, 'visualizationDir') && ~exist(params.visualizationDir, 'dir')
        mkdir(params.visualizationDir);
    elseif ~isfield(params, 'visualizationDir') && isfield(params, 'outputDir')
        params.visualizationDir = fullfile(params.outputDir, '4_Visualizations');
        if ~exist(params.visualizationDir, 'dir')
            mkdir(params.visualizationDir);
        end
    end
    
    % Create blinking directory if it doesn't exist
    if isfield(params, 'blinkingDir') && ~exist(params.blinkingDir, 'dir')
        mkdir(params.blinkingDir);
    elseif ~isfield(params, 'blinkingDir') && isfield(params, 'outputDir')
        params.blinkingDir = fullfile(params.outputDir, '3_BlinkingAnalysis');
        if ~exist(params.blinkingDir, 'dir')
            mkdir(params.blinkingDir);
        end
    end
    
    % Create statistics directory
    statsDir = fullfile(params.outputDir, '5_Statistics');
    if ~exist(statsDir, 'dir')
        mkdir(statsDir);
    end

    numSamples = min(50, numMolecules); % Show up to 50 molecules
    moleculesPerPlot = 12; % Number of molecules to show in each figure
    numPlots = ceil(numSamples/moleculesPerPlot);
    
    % Layout for each plot
    numRows = 3;
    numCols = 4;
    
    for plotIdx = 1:numPlots
        % Calculate molecule range for this plot
        startMol = (plotIdx-1)*moleculesPerPlot + 1;
        endMol = min(plotIdx*moleculesPerPlot, numSamples);
        
        % Create figure with proper size
        figure('Name', sprintf('Molecule Intensity Traces (Set %d/%d)', plotIdx, numPlots), ...
               'Position', [50 50 1200 800]);
        
        for i = startMol:endMol
            subplot(numRows, numCols, i-startMol+1);
            
            % Plot intensity trace
            plot(1:totalFrames, integrated_roi_intensity(i,:), 'b-', 'LineWidth', 1);
            hold on;
            
            % Plot mean background level and threshold
            plot([1 totalFrames], [mean_bkg_int(i) mean_bkg_int(i)], 'r--', 'LineWidth', 1);
            threshold = mean_bkg_int(i) + signalChange * std_dev_bkg_int(i);
            plot([1 totalFrames], [threshold threshold], 'g--', 'LineWidth', 1);
            
            % Highlight ON frames
            is_on = integrated_roi_intensity(i,:) > threshold;
            on_frames = find(is_on);
            if ~isempty(on_frames)
                % Create patches for ON regions
                on_regions = bwconncomp(is_on);
                for j = 1:on_regions.NumObjects
                    frames = on_regions.PixelIdxList{j};
                    y_min = min(integrated_roi_intensity(i,frames));
                    y_max = max(integrated_roi_intensity(i,frames));
                    patch([frames(1) frames(end) frames(end) frames(1)], ...
                          [y_min y_min y_max y_max], ...
                          [0.8 0.8 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                end
            end

            title(sprintf('Molecule %d', i), 'FontSize', 10);
            ylabel('Intensity', 'FontSize', 10);
            xlabel('Frame', 'FontSize', 10);
            legend('Intensity', 'Background', 'Threshold', 'Location', 'best', 'FontSize', 8);
            set(gca, 'FontSize', 9);
            grid on;
        end
        
        % Adjust subplot spacing
        sgtitle(sprintf('Intensity Traces (Molecules %d-%d)', startMol, endMol), 'FontSize', 12);
        set(gcf, 'Units', 'normalized');
        set(gcf, 'Position', [0.1 0.1 0.8 0.8]);
        
        % Save the figure
        if isfield(params, 'outputDir')
            visDir = fullfile(params.outputDir, '4_Visualizations');
            if ~exist(visDir, 'dir')
                mkdir(visDir);
            end
            saveas(gcf, fullfile(visDir, sprintf('intensity_traces_set%d.png', plotIdx)));
            saveas(gcf, fullfile(visDir, sprintf('intensity_traces_set%d.fig', plotIdx)));
        end
        close(gcf);
    end
    
    %% Convert intensities to photons and calculate additional metrics
    photon_traces = zeros(size(integrated_roi_intensity));
    photons_per_detection = zeros(numMolecules, 1);
    photons_per_cycle = zeros(numMolecules, 1);

    for n = 1:numMolecules
        % Convert to photons
        bg_subtracted = integrated_roi_intensity(n, :) - mean_bkg_int(n);
        bg_subtracted(bg_subtracted < 0) = 0;
        photon_traces(n,:) = bg_subtracted * params.camera.gain / params.camera.QE;
        
        % Calculate photons per detection (average of ON frames)
        onFrames = find(integrated_roi_intensity(n,:) > (mean_bkg_int(n) + signalChange * std_dev_bkg_int(n)));
        if ~isempty(onFrames)
            photons_per_detection(n) = mean(photon_traces(n, onFrames));
        end
    end

    %% Store results in compatible format with analyzeBlinkingFromFit
    results.validatedPoints = points(:,1:2);
    results.numValidatedPoints = numMolecules;
    results.validatedIntensityTraces = cell(numMolecules, 1);
    results.validatedPhotonTraces = cell(numMolecules, 1);
    results.validatedOnFrames = cell(numMolecules, 1);

    % Calculate blinking statistics using the same approach as analyzeBlinkingFromFit
    [blinking_counts, on_durations, off_durations, on_frames] = calculateBlinkingEvents(integrated_roi_intensity, mean_bkg_int, std_dev_bkg_int, params);
    
    % Store traces in cell arrays to match analyzeBlinkingFromFit format
    for n = 1:numMolecules
        % Store as row vectors to match analyzeBlinkingFromFit format
        results.validatedIntensityTraces{n} = integrated_roi_intensity(n,:);  % Store as row vector
        results.validatedPhotonTraces{n} = photon_traces(n,:);                % Store as row vector
        results.validatedOnFrames{n} = on_frames{n};                          % Store as column vector
    end
    
    % Store blinking statistics
    results.blinking_counts = blinking_counts;
    results.blinking_on_durations = on_durations;
    results.blinking_off_durations = off_durations;
    
    % Calculate photons per cycle
    for n = 1:numMolecules
        if blinking_counts(n) > 0
            total_photons = sum(photon_traces(n, on_frames{n}));
            photons_per_cycle(n) = total_photons / blinking_counts(n);
        end
    end

    % Store additional required fields
    results.photons_per_detection = photons_per_detection;
    results.photons_per_cycle = photons_per_cycle;
    results.duty_cycle = calculateDutyCycle(on_frames, totalFrames);
    results.active_duty_cycle = calculateActiveDutyCycle(on_durations, off_durations);
    % Calculate windowed duty cycles with 100-second windows
    [results.windowed_duty_cycles, results.time_windows] = calculateWindowedDutyCycle(on_frames, totalFrames, params.camera.frame_interval);
    % Calculate window statistics
    results.window_stats = calculateWindowStats(results.windowed_duty_cycles, results.time_windows);
    % Visualize windowed duty cycles
    visualizeSampleWindowedDutyCycle(results.windowed_duty_cycles, results.time_windows, params);
    
    if isfield(params.camera, 'frame_interval')
        frame_interval_ms = params.camera.frame_interval * 1000;
        % Convert durations to milliseconds as numeric arrays
        results.blinking_on_duration_ms = convertToMs(on_durations, frame_interval_ms);
        results.blinking_off_duration_ms = convertToMs(off_durations, frame_interval_ms);
    end
    
    % Add additional fields to match analyzeBlinkingFromFit
    if ~isfield(results, 'analysisMethod')
        results.analysisMethod = 'image_based';
    end
    
    %% Export statistics to CSV files
    % Create analysis tag for file naming
    analysisTag = sprintf('_HFSQR%d_THR%d', params.halfsquare, params.signalChange);
    
    % Export molecule-level statistics
    molecule_stats = table();
    molecule_stats.MoleculeID = (1:numMolecules)';
    molecule_stats.X = points(:,1);
    molecule_stats.Y = points(:,2);
    molecule_stats.BlinkingEvents = blinking_counts;
    molecule_stats.DutyCycle = results.duty_cycle;
    molecule_stats.ActiveDutyCycle = results.active_duty_cycle;
    molecule_stats.PhotonsPerDetection = results.photons_per_detection;
    molecule_stats.PhotonsPerCycle = results.photons_per_cycle;

    % Add mean on/off times for each molecule
    molecule_stats.MeanOnTime_ms = nan(numMolecules, 1);
    molecule_stats.MeanOffTime_ms = nan(numMolecules, 1);

    % Calculate mean on/off times for each molecule
    for n = 1:numMolecules
        if ~isempty(results.blinking_on_duration_ms)
            molecule_indices = find(results.blinking_on_duration_ms > 0);
            if ~isempty(molecule_indices)
                molecule_stats.MeanOnTime_ms(n) = mean(results.blinking_on_duration_ms(molecule_indices));
            end
        end
        
        if ~isempty(results.blinking_off_duration_ms)
            molecule_indices = find(results.blinking_off_duration_ms > 0);
            if ~isempty(molecule_indices)
                molecule_stats.MeanOffTime_ms(n) = mean(results.blinking_off_duration_ms(molecule_indices));
            end
        end
    end
    
    % Add windowed duty cycle metrics to the statistics table
    if isfield(results, 'windowed_duty_cycles') && isfield(results, 'time_windows')
        for w = 1:length(results.time_windows)
            colName = sprintf('DutyCycle_t%d', w);  % Simplified column naming
            molecule_stats.(colName) = results.windowed_duty_cycles(:, w);
        end
    end

    for i = 1:numMolecules
        if ~isempty(on_durations{i})
            molecule_stats.MeanOnTime_frames(i) = mean(on_durations{i});
        end
        if ~isempty(off_durations{i})
            molecule_stats.MeanOffTime_frames(i) = mean(off_durations{i});
        end
        if ~isempty(on_frames{i})  % Using on_frames from calculateBlinkingEvents
            molecule_stats.DutyCycle(i) = length(on_frames{i}) / totalFrames;
        end
    end

    % Convert to milliseconds if frame interval is available
    if isfield(params.camera, 'frame_interval')
        frame_interval_ms = params.camera.frame_interval * 1000;
        molecule_stats.MeanOnTime_ms = molecule_stats.MeanOnTime_frames * frame_interval_ms;
        molecule_stats.MeanOffTime_ms = molecule_stats.MeanOffTime_frames * frame_interval_ms;
    end
    
    % Save molecule statistics with analysis tag
    writetable(molecule_stats, fullfile(statsDir, ['molecule_statistics', analysisTag, '.csv']));
    
    % Create summary statistics table
    summary_stats = table();
    
    % Variables to summarize
    variables = {'BlinkingEvents', 'DutyCycle', 'ActiveDutyCycle', 'PhotonsPerDetection', 'PhotonsPerCycle'};
    
    % Add windowed duty cycle variables
    if isfield(results, 'windowed_duty_cycles')
        for w = 1:size(results.windowed_duty_cycles, 2)
            variables{end+1} = sprintf('DutyCycle_t%d', w);
        end
    end
    
    if isfield(params.camera, 'frame_interval')
        variables = [variables, {'MeanOnTime_ms', 'MeanOffTime_ms'}];
    end
    
    % Calculate summary statistics for each variable
    for i = 1:length(variables)
        var_name = variables{i};
        var_data = molecule_stats.(var_name);
        
        % Only filter out NaNs, keep zeros if they are meaningful
        valid_data = var_data(~isnan(var_data));
        
        if ~isempty(valid_data)
            row = {var_name, mean(valid_data), std(valid_data), min(valid_data), ...
                   max(valid_data), median(valid_data), length(valid_data)};
        else
            row = {var_name, NaN, NaN, NaN, NaN, NaN, 0};
        end
        
        % Add to summary table
        if i == 1
            summary_stats = cell2table(row, 'VariableNames', ...
                {'Variable', 'Mean', 'Std', 'Min', 'Max', 'Median', 'Count'});
        else
            summary_stats = [summary_stats; cell2table(row, 'VariableNames', ...
                {'Variable', 'Mean', 'Std', 'Min', 'Max', 'Median', 'Count'})];
        end
    end

    % Add windowed statistics to the summary table
    for w = 1:length(results.time_windows)
        time_point = results.time_windows(w);
        window_stats_row = {sprintf('Window_%ds', time_point), ...
                           results.window_stats.mean_duty_cycle(w), ...
                           results.window_stats.std_duty_cycle(w), ...
                           min(results.windowed_duty_cycles(:,w)), ...
                           max(results.windowed_duty_cycles(:,w)), ...
                           median(results.windowed_duty_cycles(:,w)), ...
                           results.window_stats.num_molecules(w)};
        
        summary_stats = [summary_stats; cell2table(window_stats_row, 'VariableNames', ...
            {'Variable', 'Mean', 'Std', 'Min', 'Max', 'Median', 'Count'})];
    end
    
    % Save summary statistics with analysis tag
    writetable(summary_stats, fullfile(statsDir, ['molecule_statistics_summary', analysisTag, '.csv']));
    
    % Also save a version without the tag for compatibility with visualization scripts
    writetable(summary_stats, fullfile(statsDir, 'molecule_statistics_summary.csv'));
    
    % Save the complete results structure
    save(fullfile(params.blinkingDir, ['blinking_results', analysisTag, '.mat']), 'results');
    
    % Also save a version without the tag for compatibility with visualization scripts
    save(fullfile(params.blinkingDir, 'blinking_results.mat'), 'results');

    fprintf('Completed image-based blinking analysis for %d molecules\n', numMolecules);
end

function [blinking_counts, on_durations, off_durations, on_frames] = calculateBlinkingEvents(intensity_traces, bg_mean, bg_std, params)
    numMolecules = size(intensity_traces, 1);
    blinking_counts = zeros(numMolecules, 1);
    on_durations = cell(numMolecules, 1);
    off_durations = cell(numMolecules, 1);
    on_frames = cell(numMolecules, 1);
    
    % Gap-closing threshold (in frames) - match the value used in analyzeBlinkingFromFit
    gapClosingThreshold = 1;  % Same as in analyzeBlinkingFromFit

    for i = 1:numMolecules
        trace = intensity_traces(i,:); 
        threshold = bg_mean(i) + params.signalChange * bg_std(i);
        is_on = trace > threshold;
        
        % Find frames where molecule is ON
        on_frame_indices = find(is_on);
        on_frames{i} = on_frame_indices(:);  % Force column vector
        
        if ~isempty(on_frame_indices)
            % Find gaps between ON frames
            gaps = diff(on_frame_indices);
            
            % Find gaps that are larger than the threshold
            gapIndices = find(gaps > gapClosingThreshold);
            
            if isempty(gapIndices)
                % Single event (all gaps are small enough to be closed)
                event_starts = 1;
                event_ends = length(on_frame_indices);
                % Ensure column vectors (like in analyzeBlinkingFromFit.m)
                event_starts = event_starts(:);
                event_ends = event_ends(:);
            else
               % Multiple events - ensure column vectors
               event_starts = vertcat(1, gapIndices(:) + 1);
               event_ends = vertcat(gapIndices(:), length(on_frame_indices));
            end
            
            % Count blinking events
            blinking_counts(i) = length(event_starts);
            
            % Calculate ON durations (in frames)
            on_durations{i} = zeros(length(event_starts), 1);
            for j = 1:length(event_starts)
                % Calculate duration including the gaps that were closed
                start_idx = event_starts(j);
                end_idx = event_ends(j);
                on_durations{i}(j) = on_frame_indices(end_idx) - on_frame_indices(start_idx) + 1;
            end
            
            % Calculate OFF durations (in frames)
            if length(event_starts) > 1
                off_durations{i} = zeros(length(event_starts)-1, 1);
                for j = 1:length(event_starts)-1
                    off_durations{i}(j) = on_frame_indices(event_starts(j+1)) - on_frame_indices(event_ends(j)) - 1;
                end
            else
                off_durations{i} = [];
            end
        else
            blinking_counts(i) = 0;
            on_durations{i} = [];
            off_durations{i} = [];
        end
    end
end

% Convert durations to milliseconds
function ms_durations = convertToMs(durations, frame_interval_ms)
    % Initialize output as a numeric array instead of cell array
    if iscell(durations)
        % Count total number of durations
        total_durations = 0;
        for i = 1:length(durations)
            if ~isempty(durations{i})
                total_durations = total_durations + length(durations{i});
            end
        end
        
        % Create a single numeric array
        ms_durations = zeros(total_durations, 1);
        current_idx = 1;
        
        % Convert and concatenate all durations
        for i = 1:length(durations)
            if ~isempty(durations{i})
                num_durations = length(durations{i});
                ms_durations(current_idx:current_idx+num_durations-1) = durations{i} * frame_interval_ms;
                current_idx = current_idx + num_durations;
            end
        end
    else
        % If input is already numeric, just multiply by frame interval
        ms_durations = durations * frame_interval_ms;
    end
end

function duty_cycle = calculateDutyCycle(onFrames, totalFrames)
    numMolecules = length(onFrames);
    duty_cycle = zeros(numMolecules, 1);
    
    for i = 1:numMolecules
        duty_cycle(i) = length(onFrames{i}) / totalFrames;
    end
end

% Add new function to calculate active duty cycle (ton/(ton+toff))
function active_duty_cycle = calculateActiveDutyCycle(on_durations, off_durations)
    numMolecules = length(on_durations);
    active_duty_cycle = zeros(numMolecules, 1);
    
    for i = 1:numMolecules
        total_on_time = sum(on_durations{i});
        total_off_time = sum(off_durations{i});
        
        if (total_on_time + total_off_time) > 0
            active_duty_cycle(i) = total_on_time / (total_on_time + total_off_time);
        end
    end
end

% Add new function to calculate time-windowed duty cycles
function [windowed_duty_cycles, time_windows] = calculateWindowedDutyCycle(on_frames, totalFrames, frame_interval)
    numMolecules = length(on_frames);
    
    % Fixed window size of 100 seconds
    window_size_seconds = 100;
    frames_per_window = round(window_size_seconds / frame_interval);
    
    % Calculate number of complete windows
    num_windows = floor(totalFrames * frame_interval / window_size_seconds);
    if num_windows < 1
        num_windows = 1;
    end
    
    % Initialize outputs
    windowed_duty_cycles = zeros(numMolecules, num_windows);
    time_windows = (1:num_windows) * window_size_seconds;
    
    % Calculate duty cycles for each window
    for i = 1:numMolecules
        if ~isempty(on_frames{i})
            for w = 1:num_windows
                start_frame = 1 + (w-1) * frames_per_window;
                end_frame = min(w * frames_per_window, totalFrames);
                window_frames = start_frame:end_frame;
                
                on_in_window = ismember(on_frames{i}, window_frames);
                windowed_duty_cycles(i, w) = sum(on_in_window) / length(window_frames);
            end
        end
    end
end

function window_stats = calculateWindowStats(windowed_duty_cycles, time_windows)
    num_windows = length(time_windows);
    window_stats = struct();
    
    % Initialize arrays for statistics
    window_stats.time_points = time_windows;
    window_stats.mean_duty_cycle = zeros(num_windows, 1);
    window_stats.std_duty_cycle = zeros(num_windows, 1);
    window_stats.num_molecules = zeros(num_windows, 1);
    
    % Calculate statistics for each window
    for w = 1:num_windows
        window_data = windowed_duty_cycles(:, w);
        valid_data = window_data(~isnan(window_data));
        
        window_stats.mean_duty_cycle(w) = mean(valid_data, 'omitnan');
        window_stats.std_duty_cycle(w) = std(valid_data, 'omitnan');
        window_stats.num_molecules(w) = sum(~isnan(window_data));
    end
end

function visualizeSampleWindowedDutyCycle(windowed_duty_cycles, time_windows, params)
    figure('Name', 'Windowed Duty Cycles', 'Position', [100, 100, 1000, 800]);
    
    % Create subplot for individual traces
    subplot(2,1,1);
    numSamples = min(10, size(windowed_duty_cycles, 1));
    
    % Plot individual molecules
    for i = 1:numSamples
        plot(time_windows, windowed_duty_cycles(i,:), 'o-', 'LineWidth', 1, ...
            'DisplayName', sprintf('Molecule %d', i));
        hold on;
    end
    
    % Calculate and plot statistics
    mean_dc = mean(windowed_duty_cycles, 1, 'omitnan');
    std_dc = std(windowed_duty_cycles, 0, 1, 'omitnan');
    
    % Plot mean with error bars
    subplot(2,1,2);
    errorbar(time_windows, mean_dc, std_dc, 'k-o', 'LineWidth', 2, ...
        'MarkerSize', 8, 'DisplayName', 'Mean ± Std');
    
    % Format both subplots
    for sp = 1:2
        subplot(2,1,sp);
        xlabel('Time (seconds)', 'FontSize', 14);
        ylabel('Duty Cycle', 'FontSize', 14);
        grid on;
        if sp == 1
            title('Individual Molecule Traces', 'FontSize', 16);
            legend('Location', 'best');
        else
            title('Population Statistics', 'FontSize', 16);
        end
    end
    
    % Save figures
    if isfield(params, 'visualizationDir')
        saveas(gcf, fullfile(params.visualizationDir, 'windowed_duty_cycles.png'));
        saveas(gcf, fullfile(params.visualizationDir, 'windowed_duty_cycles.fig'));
    end
    close(gcf);
end
