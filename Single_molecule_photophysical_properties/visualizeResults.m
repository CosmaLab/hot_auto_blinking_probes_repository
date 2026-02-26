function visualizeResults(resultsFilePath)
% VISUALIZERESULTS Generate visualizations for single-molecule blinking analysis
%   VISUALIZERESULTS(resultsFilePath) loads single-molecule blinking analysis 
%   results and generates comprehensive visualizations including intensity 
%   traces and statistical distributions.

% Load the results file
if ~exist(resultsFilePath, 'file')
    error('Results file not found: %s', resultsFilePath);
end

fprintf('Loading blinking results from: %s\n', resultsFilePath);
loadedData = load(resultsFilePath);

% Check if the loaded data contains blinkResults or results
if isfield(loadedData, 'blinkResults')
    results = loadedData.blinkResults;
elseif isfield(loadedData, 'results')
    results = loadedData.results;
else
    error('No valid results found in the file. Expected fields: blinkResults or results');
end

% Get the correct visualization directory
[resultsDir, ~, ~] = fileparts(resultsFilePath);
[analysisDir, ~, ~] = fileparts(resultsDir);
visualizationDir = fullfile(analysisDir, '4_Visualizations');

% Create directory if it doesn't exist
if ~exist(visualizationDir, 'dir')
    mkdir(visualizationDir);
end

% Check if we have validated points
if ~isfield(results, 'validatedPoints') || isempty(results.validatedPoints)
    warning('No validated points found in results');
    return;
end

% Create a figure for intensity and photon traces
try
    createTracePlot(results, visualizationDir);
catch ME
    warning(ME.identifier, '%s', ME.message);
end

% Create basic blinking statistics plots
try
    createBlinkingStatsPlot(results, visualizationDir);
catch ME
    warning(ME.identifier, '%s', ME.message);
end

try
    createOnOffVisualization(results, visualizationDir);
catch ME
    warning(ME.identifier, '%s', ME.message);
end

end

function createTracePlot(results, visualizationDir)
    % Create intensity and photon traces plot
    h = figure('Name', 'Molecule Traces', 'Position', [100, 100, 1200, 800], 'Visible', 'on');

    % Determine subplot layout based on number of points
    numPoints = results.numValidatedPoints;

    % Sample points if there are too many
    maxPointsToShow = 8; % Maximum number of points to display (reduced to accommodate both traces)
    if numPoints > maxPointsToShow
        fprintf('Large dataset detected (%d points). Sampling %d points for visualization.\n', numPoints, maxPointsToShow);
        sampledIndices = round(linspace(1, numPoints, maxPointsToShow));
    else
        sampledIndices = 1:numPoints;
    end

    numSampledPoints = length(sampledIndices);
    
    % Check if we have photon traces
    hasPhotonTraces = isfield(results, 'validatedPhotonTraces') && ...
                     ~isempty(results.validatedPhotonTraces);
    
    % Create a figure for all sampled points
    for idx = 1:numSampledPoints
        i = sampledIndices(idx);
        
        % Intensity trace subplot
        subplot(numSampledPoints, 2, 2*idx-1);
        
        % Get intensity trace
        if isfield(results, 'validatedIntensityTraces')
            intensityTrace = results.validatedIntensityTraces{i};
        else
            warning('No intensity traces found');
            break;
        end
        
        % Plot intensity trace with better formatting
        plot(intensityTrace, 'b-', 'LineWidth', 1);
        title(sprintf('Molecule %d - Intensity', i), 'FontSize', 10);
        xlabel('Frame', 'FontSize', 9);
        ylabel('Intensity (counts)', 'FontSize', 9);
        
        % Highlight on frames if available
        if isfield(results, 'validatedOnFrames') && length(results.validatedOnFrames) >= i
            onFrames = results.validatedOnFrames{i};
            if ~isempty(onFrames)
                hold on;
                plot(onFrames, intensityTrace(onFrames), 'ro', 'MarkerSize', 2);
            end
        end
        
        % Improve plot appearance
        grid on;
        box on;
        set(gca, 'FontSize', 8);
        
        % Set reasonable y-axis limits
        if any(intensityTrace > 0)
            ylim([0, max(intensityTrace) * 1.2]);
        end
        
        % Photon trace subplot
        subplot(numSampledPoints, 2, 2*idx);
        
        % Get photon trace
        if hasPhotonTraces && length(results.validatedPhotonTraces) >= i
            photonTrace = results.validatedPhotonTraces{i};
            
            % Plot photon trace
            plot(photonTrace, 'g-', 'LineWidth', 1);
            title(sprintf('Molecule %d - Photons', i), 'FontSize', 10);
            xlabel('Frame', 'FontSize', 9);
            ylabel('Photons', 'FontSize', 9);
            
            % Highlight on frames if available
            if isfield(results, 'validatedOnFrames') && length(results.validatedOnFrames) >= i
                onFrames = results.validatedOnFrames{i};
                if ~isempty(onFrames)
                    hold on;
                    plot(onFrames, photonTrace(onFrames), 'ro', 'MarkerSize', 2);
                end
            end
            
            % Improve plot appearance
            grid on;
            box on;
            set(gca, 'FontSize', 8);
            
            % Set reasonable y-axis limits
            if any(photonTrace > 0)
                ylim([0, max(photonTrace) * 1.2]);
            end
        else
            % If no photon traces, create a blank subplot with a message
            text(0.5, 0.5, 'No photon trace data available', ...
                'HorizontalAlignment', 'center', 'FontSize', 9);
            axis off;
        end
    end

    % Adjust layout
    sgtitle('Single-Molecule Intensity and Photon Traces', 'FontSize', 14);
    set(gcf, 'Color', 'w');

    % Save figure and close
    saveas(h, fullfile(visualizationDir, 'molecule_traces.fig'));
    saveas(h, fullfile(visualizationDir, 'molecule_traces.png'));
    fprintf('Saved molecule traces plot to %s\n', visualizationDir);
    
    % Close the figure
    if ishandle(h)
        close(h);
    end
end

function createBlinkingStatsPlot(results, visualizationDir)
    % Create individual plots for each statistic
    createBlinkingEventsPlot(results, visualizationDir);
    createOnTimesPlot(results, visualizationDir);
    createOffTimesPlot(results, visualizationDir);
    createDutyCyclePlot(results, visualizationDir);
    createPhotonStatsPlot(results, visualizationDir);
end

function createBlinkingEventsPlot(results, visualizationDir)
    if ~isfield(results, 'blinking_counts')
        warning('No blinking counts data available');
        return;
    end
    
    h = figure('Name', 'Blinking Events', 'Position', [100, 100, 600, 500], 'Visible', 'on');
    data = results.blinking_counts;
    histObj = histogram(data, 20, 'Normalization', 'probability');
    counts = histObj.Values;
    edges = histObj.BinEdges;
    title('Blinking Events Distribution');
    xlabel('Number of Events');
    ylabel('Probability');
    grid on;
    
    % Add statistics
    hold on;
    meanVal = mean(data);
    stdVal = std(data);
    xline(meanVal, 'r-', sprintf('Mean: %.1f', meanVal), 'LineWidth', 1.5);
    xline(meanVal + stdVal, 'g--', sprintf('Std: %.1f', stdVal));
    hold off;
    
    saveFigure(h, visualizationDir, 'blinking_events');
end

function createOnTimesPlot(results, visualizationDir)
    if ~isfield(results, 'blinking_on_duration_ms')
        warning('No ON duration data available');
        return;
    end
    
    h = figure('Name', 'ON Times', 'Position', [100, 100, 600, 500], 'Visible', 'on');
    data = results.blinking_on_duration_ms/1000; % Convert to seconds
    histObj = histogram(data, 20, 'Normalization', 'probability');
    counts = histObj.Values;
    edges = histObj.BinEdges;
    title('ON Times Distribution');
    xlabel('Time (s)');
    ylabel('Probability');
    grid on;
    
    % Fit exponential decay if possible
    fitExponentialDecay(data, counts, edges);
    
    saveFigure(h, visualizationDir, 'on_times');
end

function createOffTimesPlot(results, visualizationDir)
    if ~isfield(results, 'blinking_off_duration_ms')
        warning('No OFF duration data available');
        return;
    end
    
    h = figure('Name', 'OFF Times', 'Position', [100, 100, 600, 500], 'Visible', 'on');
    data = results.blinking_off_duration_ms/1000; % Convert to seconds
    histObj = histogram(data, 20, 'Normalization', 'probability');
    counts = histObj.Values;
    edges = histObj.BinEdges;
    title('OFF Times Distribution');
    xlabel('Time (s)');
    ylabel('Probability');
    grid on;
    
    % Fit exponential decay if possible
    fitExponentialDecay(data, counts, edges);
    
    saveFigure(h, visualizationDir, 'off_times');
end

function createDutyCyclePlot(results, visualizationDir)
    if ~isfield(results, 'duty_cycle')
        warning('No duty cycle data available');
        return;
    end
    
    h = figure('Name', 'Duty Cycle', 'Position', [100, 100, 600, 500], 'Visible', 'on');
    data = results.duty_cycle;
    histogram(data, 20, 'Normalization', 'probability');
    title('Duty Cycle Distribution');
    xlabel('Duty Cycle');
    ylabel('Probability');
    grid on;
    
    % Add statistics with proper formatting for small numbers
    hold on;
    meanVal = mean(data);
    stdVal = std(data);
    
    % Format the numbers based on their magnitude
    if meanVal < 0.01
        meanStr = sprintf('%.4f', meanVal);
    elseif meanVal < 0.1
        meanStr = sprintf('%.3f', meanVal);
    else
        meanStr = sprintf('%.2f', meanVal);
    end
    
    if stdVal < 0.01
        stdStr = sprintf('%.4f', stdVal);
    elseif stdVal < 0.1
        stdStr = sprintf('%.3f', stdVal);
    else
        stdStr = sprintf('%.2f', stdVal);
    end
    
    text(0.7, 0.8, sprintf('Mean = %s\nStd = %s', meanStr, stdStr), ...
        'Units', 'normalized', 'FontSize', 9);
    hold off;
    
    saveFigure(h, visualizationDir, 'duty_cycle');
end

function createPhotonStatsPlot(results, visualizationDir)
    h = figure('Name', 'Photon Statistics', 'Position', [100, 100, 1200, 500], 'Visible', 'on');
    
    % Plot 1: Photons per detection
    subplot(1,2,1);
    if isfield(results, 'photons_per_detection')
        data = results.photons_per_detection;
        histObj = histogram(data, 20, 'Normalization', 'probability');
        counts = histObj.Values;
        edges = histObj.BinEdges;
        title('Photon Counts per Detection');
        xlabel('Photons per ON Frame');
        ylabel('Probability');
        grid on;
        fitExponentialDecay(data, counts, edges);
    end
    
    % Plot 2: Photons per cycle
    subplot(1,2,2);
    if isfield(results, 'photons_per_cycle')
        data = results.photons_per_cycle(results.photons_per_cycle > 0);
        histObj = histogram(data, 20, 'Normalization', 'probability');
        counts = histObj.Values;
        edges = histObj.BinEdges;
        title('Photons per Switching Cycle');
        xlabel('Photons per Cycle');
        ylabel('Probability');
        grid on;
        fitExponentialDecay(data, counts, edges);
    end
    
    sgtitle('Photon Statistics', 'FontSize', 14);
    saveFigure(h, visualizationDir, 'photon_statistics');
end

function fitExponentialDecay(data, counts, edges)
    centers = (edges(1:end-1) + edges(2:end))/2;
    validIdx = counts > 0;
    if sum(validIdx) > 2
        try
            % Use fitdist for exponential distribution like in Dye_bleachcurve.m
            pd = fitdist(data(:), 'Exponential');
            x = linspace(0, max(data), 100);
            y = exppdf(x, pd.mu);
            
            % Scale the PDF to match histogram height
            scaleFactor = max(counts) / max(y);
            y = y * scaleFactor;
            
            hold on;
            plot(x, y, 'r-', 'LineWidth', 2, 'DisplayName', 'Fit (Exp)');
            
            % Display k = 1/mu as in Dye_bleachcurve.m
            k = 1/pd.mu;
            meanVal = mean(data);
            stdVal = std(data);
            text(0.7*max(centers), 0.8*max(counts), ...
                sprintf('k = %.3f\nMean = %.2f\nStd = %.2f', ...
                k, meanVal, stdVal), 'FontSize', 9);
            hold off;
            
        catch
            % If fitting fails, just show basic statistics
            hold on;
            meanVal = mean(data);
            stdVal = std(data);
            text(0.7, 0.8, sprintf('Mean = %.2f\nStd = %.2f', meanVal, stdVal), ...
                'Units', 'normalized', 'FontSize', 9);
            hold off;
        end
    end
end

function saveFigure(h, visualizationDir, name)
    set(h, 'Color', 'w');
    saveas(h, fullfile(visualizationDir, [name '.fig']));
    saveas(h, fullfile(visualizationDir, [name '.png']));
    if ishandle(h)
        close(h);
    end
end

function createOnOffVisualization(results, visualizationDir)
    % Create on/off time visualization for sampled points
    h = figure('Name', 'On/Off Time Visualization', 'Position', [100, 100, 1200, 800], 'Visible', 'on');
    
    % Determine number of points to sample
    numPoints = results.numValidatedPoints;
    numSampledPoints = min(5, numPoints); % Sample up to 5 points
    
    if numSampledPoints == 0
        warning('No points available for on/off visualization');
        return;
    end
    
    % Sample points evenly across the dataset
    sampledIndices = round(linspace(1, numPoints, numSampledPoints));
    
    % Check if we have the required fields
    if ~isfield(results, 'validatedIntensityTraces') || ~isfield(results, 'validatedOnFrames')
        warning('Missing required fields for on/off visualization');
        return;
    end
    
    % Create a 1x2 subplot for each sampled point
    for i = 1:numSampledPoints
        idx = sampledIndices(i);
        
        % Get intensity trace and on frames
        intensityTrace = results.validatedIntensityTraces{idx};
        onFrames = results.validatedOnFrames{idx};
        
        % Create binary on/off trace
        binaryTrace = zeros(size(intensityTrace));
        binaryTrace(onFrames) = 1;
        
        % Calculate duty cycle
        dutyCycle = sum(binaryTrace) / length(binaryTrace);
        
        % Use pre-calculated durations and throw error if not available
        if isfield(results, 'blinking_on_durations') && length(results.blinking_on_durations) >= idx
            onDurations = results.blinking_on_durations{idx};
        else
            error('Pre-calculated ON durations not found in results. Field blinking_on_durations is required.');
        end
        
        % Use pre-calculated off durations and throw error if not available
        if isfield(results, 'blinking_off_durations') && length(results.blinking_off_durations) >= idx
            offDurations = results.blinking_off_durations{idx};
        else
            error('Pre-calculated OFF durations not found in results. Field blinking_off_durations is required.');
        end

        % Ensure onDurations and offDurations are column vectors
        onDurations = onDurations(:);
        offDurations = offDurations(:);
        
        % Create subplot for intensity trace
        subplot(numSampledPoints, 2, 2*i-1);
        plot(intensityTrace, 'b-', 'LineWidth', 1);
        hold on;
        plot(onFrames, intensityTrace(onFrames), 'ro', 'MarkerSize', 2);
        hold off;
        title(sprintf('Molecule %d - Intensity', idx), 'FontSize', 10);
        xlabel('Frame', 'FontSize', 9);
        ylabel('Intensity (counts)', 'FontSize', 9);
        grid on;
        
        % Create subplot for binary on/off trace
        subplot(numSampledPoints, 2, 2*i);
        stairs(1:length(binaryTrace), binaryTrace, 'LineWidth', 1.5);
        ylim([-0.1, 1.1]);
        title(sprintf('Molecule %d - On/Off States', idx), 'FontSize', 10);
        xlabel('Frame', 'FontSize', 9);
        ylabel('State (1=On, 0=Off)', 'FontSize', 9);
        grid on;

        % Add statistics
        meanOnDuration = mean(onDurations);
        stdOnDuration = std(onDurations);
        meanOffDuration = mean(offDurations);
        stdOffDuration = std(offDurations);
        
        statsText = sprintf('Duty Cycle: %.3f\nOn Time: %.1f ± %.1f frames\nOff Time: %.1f ± %.1f frames', ...
            dutyCycle, meanOnDuration, stdOnDuration, meanOffDuration, stdOffDuration);
        
        text(0.7, 0.8, statsText, 'Units', 'normalized', 'FontSize', 8, ...
            'BackgroundColor', [1 1 1 0.7]);
    end
    
    % Adjust layout and save
    sgtitle('Single-Molecule On/Off Time Analysis', 'FontSize', 14);
    set(gcf, 'Color', 'w');
    saveas(h, fullfile(visualizationDir, 'on_off_visualization.fig'));
    saveas(h, fullfile(visualizationDir, 'on_off_visualization.png'));
    
    if ishandle(h)
        close(h);
    end
end