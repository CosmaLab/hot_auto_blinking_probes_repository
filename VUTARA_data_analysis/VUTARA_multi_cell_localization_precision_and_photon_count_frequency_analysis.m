% VUTARA_multi_cell_localization_precision_and_photon_count_frequency_analysis.m
% 
% DESCRIPTION:
%   This script processes multiple Vutara CSV files (one per cell) to generate 
%   frequency distribution histograms of photon counts and localization precision 
%   (X, Y, Z) across all cells. It calculates mean and standard error of the mean 
%   (SEM) for each bin across cells, and produces publication-quality figures 
%   with customizable parameters. Output includes PDF, TIFF, FIG files, and 
%   summary CSV tables.
%
% Clear workspace and close all figures
clear all; close all; clc;

% ============================================
% USER-CUSTOMIZABLE PARAMETERS SECTION
% ============================================

% Define file path
filePath = 'E:\Project\8 Code\Vutara\Vutara 3D\CF3';

% 1. CUSTOMIZE BIN WIDTHS FOR EACH PARAMETER
photonBinWidth = 100;         % Bin width for photon counts (photons)
precisionBinWidthX = 4;       % Bin width for X precision (nm)
precisionBinWidthY = 4;       % Bin width for Y precision (nm)
precisionBinWidthZ = 8;       % Bin width for Z precision (nm)

% 2. CUSTOMIZE ERROR BAR HORIZONTAL LINE LENGTH
% This controls how long the horizontal lines on error bars are
% Value is a fraction of the bin width (e.g., 0.2 = 20% of bin width)
errorbarLineScale = 0.2;      % Default: 0.2 (20% of bin width)

% 3. CUSTOMIZE BAR WIDTH (applies to all bar plots)
barWidth = 0.7;               % Width of histogram bars (0 to 1)

% 4. CUSTOMIZE COLORS
photonColor = [0.2 0.6 0.8];  % Color for photon count histogram
xColor = [0.2 0.6 0.8];       % Color for X precision histogram
yColor = [0.2 0.6 0.8];       % Color for Y precision histogram
zColor = [0.2 0.6 0.8];       % Color for Z precision histogram

% 5. CUSTOMIZE PLOT RANGES (optional - auto-calculated if not specified)
% Set to [] for automatic calculation based on percentiles
% Example: photonRange = [0, 2000]; % Force range from 0 to 2000 photons
photonRange = [];             % Auto-calculate photon range
precisionRangeX = [];         % Auto-calculate X precision range
precisionRangeY = [];         % Auto-calculate Y precision range
precisionRangeZ = [];         % Auto-calculate Z precision range

% 6. CUSTOMIZE MAXIMUM DISPLAY RANGE FOR PRECISION PLOTS (for visualization only)
maxDisplayPrecisionX = 100;   % Maximum X value to display (nm)
maxDisplayPrecisionY = 100;   % Maximum Y value to display (nm)
maxDisplayPrecisionZ = 200;   % Maximum Z value to display (nm)

% ============================================
% END OF USER-CUSTOMIZABLE PARAMETERS
% ============================================

% Display parameter settings
fprintf('=== CUSTOMIZABLE PARAMETERS ===\n');
fprintf('Photon bin width: %d photons\n', photonBinWidth);
fprintf('X precision bin width: %d nm\n', precisionBinWidthX);
fprintf('Y precision bin width: %d nm\n', precisionBinWidthY);
fprintf('Z precision bin width: %d nm\n', precisionBinWidthZ);
fprintf('Error bar line scale: %.2f (%.0f%% of bin width)\n', errorbarLineScale, errorbarLineScale*100);
fprintf('Bar width: %.2f\n', barWidth);
fprintf('===============================\n\n');

% List all CSV files
csvFiles = dir(fullfile(filePath, '*.csv'));

% Check if CSV files exist
if isempty(csvFiles)
    error('No CSV files found in the specified path!');
end

fprintf('Found %d CSV files\n', length(csvFiles));

% Initialize arrays to store frequency distributions
cellNames = {};
cellFrequencyPhotons = [];  % Will store frequency counts for each cell
cellFrequencyX = [];        % Frequency counts for X precision
cellFrequencyY = [];        % Frequency counts for Y precision
cellFrequencyZ = [];        % Frequency counts for Z precision

% Define common bin edges for all cells (to make distributions comparable)
fprintf('\nDefining common bin edges for all cells...\n');

% First pass: find global ranges and collect all data for overall median
allPhotons = [];
allXPrecision = [];
allYPrecision = [];
allZPrecision = [];

% Process each CSV file to find global ranges
for i = 1:length(csvFiles)
    fileName = csvFiles(i).name;
    fullPath = fullfile(filePath, fileName);
    
    % Skip summary files using strfind for MATLAB 2016a compatibility
    skipFile = false;
    if ~isempty(strfind(fileName, 'Summary_Statistics')) || ...
       ~isempty(strfind(fileName, 'Summary_')) || ...
       ~isempty(strfind(fileName, 'SummaryPhoton')) || ...
       ~isempty(strfind(fileName, 'SummaryPrecision')) || ...
       ~isempty(strfind(fileName, 'SummaryOverview')) || ...
       ~isempty(strfind(fileName, 'Combined_')) || ...
       ~isempty(strfind(fileName, 'CellLevel_')) || ...
       ~isempty(strfind(fileName, 'Frequency_'))
        fprintf('Skipping summary file: %s\n', fileName);
        skipFile = true;
    end
    
    if skipFile
        continue;
    end
    
    try
        % Read CSV data
        data = readtable(fullPath);
        
        % Check if required columns exist
        requiredColumns = [14, 29, 30, 31];
        if size(data, 2) < max(requiredColumns)
            continue;
        end
        
        % Extract required columns
        photons = data{:, 14};
        xPrecision = data{:, 29};
        yPrecision = data{:, 30};
        zPrecision = data{:, 31};
        
        % Remove NaN values
        validIdx = ~isnan(photons) & ~isnan(xPrecision) & ...
                   ~isnan(yPrecision) & ~isnan(zPrecision);
        photons = photons(validIdx);
        xPrecision = xPrecision(validIdx);
        yPrecision = yPrecision(validIdx);
        zPrecision = zPrecision(validIdx);
        
        % Collect for global range calculation
        allPhotons = [allPhotons; photons];
        allXPrecision = [allXPrecision; xPrecision];
        allYPrecision = [allYPrecision; yPrecision];
        allZPrecision = [allZPrecision; zPrecision];
        
    catch ME
        continue;
    end
end

% Check if we have data
if isempty(allPhotons)
    error('No valid data was processed!');
end

% Calculate overall medians for combined data
overallMedianPhotons = median(allPhotons);
overallMedianX = median(allXPrecision);
overallMedianY = median(allYPrecision);
overallMedianZ = median(allZPrecision);

fprintf('Overall medians from combined data:\n');
fprintf('  Photons: %.1f\n', overallMedianPhotons);
fprintf('  X Precision: %.2f nm\n', overallMedianX);
fprintf('  Y Precision: %.2f nm\n', overallMedianY);
fprintf('  Z Precision: %.2f nm\n', overallMedianZ);

% Define common bin edges for all parameters using CUSTOM bin widths
% For photons
if isempty(photonRange)
    photonMin = floor(prctile(allPhotons, 1)/photonBinWidth)*photonBinWidth;  % 1st percentile
    photonMax = ceil(prctile(allPhotons, 99)/photonBinWidth)*photonBinWidth;  % 99th percentile
else
    photonMin = photonRange(1);
    photonMax = photonRange(2);
end

photonEdges = photonMin:photonBinWidth:photonMax;
photonCenters = (photonEdges(1:end-1) + photonEdges(2:end))/2;
numPhotonBins = length(photonCenters);

% For X precision
if isempty(precisionRangeX)
    precisionMinX = 0;  % Typically precision starts at 0
    precisionMaxX = ceil(prctile(allXPrecision, 99)/precisionBinWidthX)*precisionBinWidthX;
else
    precisionMinX = precisionRangeX(1);
    precisionMaxX = precisionRangeX(2);
end

precisionEdgesX = precisionMinX:precisionBinWidthX:precisionMaxX;
precisionCentersX = (precisionEdgesX(1:end-1) + precisionEdgesX(2:end))/2;
numPrecisionBinsX = length(precisionCentersX);

% For Y precision
if isempty(precisionRangeY)
    precisionMinY = 0;
    precisionMaxY = ceil(prctile(allYPrecision, 99)/precisionBinWidthY)*precisionBinWidthY;
else
    precisionMinY = precisionRangeY(1);
    precisionMaxY = precisionRangeY(2);
end

precisionEdgesY = precisionMinY:precisionBinWidthY:precisionMaxY;
precisionCentersY = (precisionEdgesY(1:end-1) + precisionEdgesY(2:end))/2;
numPrecisionBinsY = length(precisionCentersY);

% For Z precision
if isempty(precisionRangeZ)
    precisionMinZ = 0;
    precisionMaxZ = ceil(prctile(allZPrecision, 99)/precisionBinWidthZ)*precisionBinWidthZ;
else
    precisionMinZ = precisionRangeZ(1);
    precisionMaxZ = precisionRangeZ(2);
end

precisionEdgesZ = precisionMinZ:precisionBinWidthZ:precisionMaxZ;
precisionCentersZ = (precisionEdgesZ(1:end-1) + precisionEdgesZ(2:end))/2;
numPrecisionBinsZ = length(precisionCentersZ);

% Find maximum precision bins for display
maxPrecisionBins = max([numPrecisionBinsX, numPrecisionBinsY, numPrecisionBinsZ]);

fprintf('\nDefined bin edges:\n');
fprintf('  Photons: %d bins (%.0f-%.0f, width=%d)\n', numPhotonBins, photonMin, photonMax, photonBinWidth);
fprintf('  X Precision: %d bins (%.1f-%.1f nm, width=%d nm)\n', numPrecisionBinsX, precisionMinX, precisionMaxX, precisionBinWidthX);
fprintf('  Y Precision: %d bins (%.1f-%.1f nm, width=%d nm)\n', numPrecisionBinsY, precisionMinY, precisionMaxY, precisionBinWidthY);
fprintf('  Z Precision: %d bins (%.1f-%.1f nm, width=%d nm)\n', numPrecisionBinsZ, precisionMinZ, precisionMaxZ, precisionBinWidthZ);

% Second pass: calculate frequency distributions for each cell
fprintf('\nCalculating frequency distributions for each cell...\n');

validCellCount = 0;
for i = 1:length(csvFiles)
    fileName = csvFiles(i).name;
    fullPath = fullfile(filePath, fileName);
    
    % Skip summary files using strfind for MATLAB 2016a compatibility
    skipFile = false;
    if ~isempty(strfind(fileName, 'Summary_Statistics')) || ...
       ~isempty(strfind(fileName, 'Summary_')) || ...
       ~isempty(strfind(fileName, 'SummaryPhoton')) || ...
       ~isempty(strfind(fileName, 'SummaryPrecision')) || ...
       ~isempty(strfind(fileName, 'SummaryOverview')) || ...
       ~isempty(strfind(fileName, 'Combined_')) || ...
       ~isempty(strfind(fileName, 'CellLevel_')) || ...
       ~isempty(strfind(fileName, 'Frequency_'))
        skipFile = true;
    end
    
    if skipFile
        continue;
    end
    
    fprintf('Processing file: %s (%d/%d)\n', fileName, i, length(csvFiles));
    
    try
        % Read CSV data
        data = readtable(fullPath);
        
        % Check if required columns exist
        requiredColumns = [14, 29, 30, 31];
        if size(data, 2) < max(requiredColumns)
            fprintf('  Warning: File %s has insufficient columns, skipping\n', fileName);
            continue;
        end
        
        % Extract required columns
        photons = data{:, 14};
        xPrecision = data{:, 29};
        yPrecision = data{:, 30};
        zPrecision = data{:, 31};
        
        % Remove NaN values
        validIdx = ~isnan(photons) & ~isnan(xPrecision) & ...
                   ~isnan(yPrecision) & ~isnan(zPrecision);
        photons = photons(validIdx);
        xPrecision = xPrecision(validIdx);
        yPrecision = yPrecision(validIdx);
        zPrecision = zPrecision(validIdx);
        
        % Calculate frequency distributions for this cell using CUSTOM bin edges
        cellFreqPhotons = histcounts(photons, photonEdges);
        cellFreqX = histcounts(xPrecision, precisionEdgesX);
        cellFreqY = histcounts(yPrecision, precisionEdgesY);
        cellFreqZ = histcounts(zPrecision, precisionEdgesZ);
        
        % Store results
        validCellCount = validCellCount + 1;
        cellNames{validCellCount} = fileName;
        
        if validCellCount == 1
            % Initialize arrays
            cellFrequencyPhotons = cellFreqPhotons;
            cellFrequencyX = cellFreqX;
            cellFrequencyY = cellFreqY;
            cellFrequencyZ = cellFreqZ;
        else
            % Add as new rows
            cellFrequencyPhotons = [cellFrequencyPhotons; cellFreqPhotons];
            cellFrequencyX = [cellFrequencyX; cellFreqX];
            cellFrequencyY = [cellFrequencyY; cellFreqY];
            cellFrequencyZ = [cellFrequencyZ; cellFreqZ];
        end
        
        fprintf('  Cell %d: %d localizations\n', validCellCount, length(photons));
        
    catch ME
        fprintf('  Error processing file %s: %s\n', fileName, ME.message);
    end
end

numCells = validCellCount;
fprintf('\n=== FREQUENCY DISTRIBUTION SUMMARY ===\n');
fprintf('Total cells processed: %d\n', numCells);

% Calculate statistics for each bin across cells
fprintf('\nCalculating statistics for each bin across cells...\n');

% For each bin, calculate mean and SEM across cells
meanFreqPhotons = mean(cellFrequencyPhotons, 1);
semFreqPhotons = std(cellFrequencyPhotons, 0, 1) / sqrt(numCells);

meanFreqX = mean(cellFrequencyX, 1);
semFreqX = std(cellFrequencyX, 0, 1) / sqrt(numCells);

meanFreqY = mean(cellFrequencyY, 1);
semFreqY = std(cellFrequencyY, 0, 1) / sqrt(numCells);

meanFreqZ = mean(cellFrequencyZ, 1);
semFreqZ = std(cellFrequencyZ, 0, 1) / sqrt(numCells);

% ============================================
% 1. PHOTON COUNT FREQUENCY DISTRIBUTION
% ============================================
fprintf('\nCreating photon count frequency distribution plot...\n');

fig1 = figure('Visible', 'off');
set(gcf, 'Position', [100, 100, 1400, 600]);

% Main plot: mean frequency with SEM error bars
subplot(1, 2, 1);
hold on;

% Plot mean frequencies as bars with CUSTOM bar width
bar(photonCenters, meanFreqPhotons, barWidth, 'FaceColor', photonColor, 'EdgeColor', 'none');

% Add SEM error bars with CUSTOM horizontal line length
for i = 1:numPhotonBins
    if meanFreqPhotons(i) > 0
        xPos = photonCenters(i);
        yPos = meanFreqPhotons(i);
        semVal = semFreqPhotons(i);
        
        % Draw error bar manually
        line([xPos, xPos], [yPos - semVal, yPos + semVal], 'Color', 'k', 'LineWidth', 1);
        % Use CUSTOM error bar line scale
        line([xPos - photonBinWidth*errorbarLineScale, xPos + photonBinWidth*errorbarLineScale], ...
             [yPos - semVal, yPos - semVal], 'Color', 'k', 'LineWidth', 1);
        line([xPos - photonBinWidth*errorbarLineScale, xPos + photonBinWidth*errorbarLineScale], ...
             [yPos + semVal, yPos + semVal], 'Color', 'k', 'LineWidth', 1);
    end
end

% Add red vertical line for overall median from combined data
ylimits = ylim;
line([overallMedianPhotons, overallMedianPhotons], [0, ylimits(2)], ...
    'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
text(overallMedianPhotons + photonBinWidth*0.05, ylimits(2)*0.9, ...
    sprintf('Median: %.0f', overallMedianPhotons), ...
    'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');

xlabel('Photons');
ylabel('Mean Frequency � SEM (across cells)');
title(sprintf('Photon Count Frequency Distribution (n=%d cells)', numCells));
grid on;

% Add overall statistics
totalMean = mean(sum(cellFrequencyPhotons, 2));
totalSEM = std(sum(cellFrequencyPhotons, 2)) / sqrt(numCells);
text(0.05, 0.95, sprintf('Total per cell: %.0f � %.0f', totalMean, totalSEM), ...
    'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', ...
    'BackgroundColor', 'white', 'EdgeColor', 'black');

% Add bin width annotation
text(0.05, 0.85, sprintf('Bin width: %d photons', photonBinWidth), ...
    'Units', 'normalized', 'FontSize', 9, 'BackgroundColor', 'white', 'EdgeColor', 'black');

% Subplot 2: Individual cell distributions
subplot(1, 2, 2);
hold on;

% Plot each cell's distribution as thin lines
for i = 1:min(numCells, 10)  % Plot first 10 cells for clarity
    plot(photonCenters, cellFrequencyPhotons(i, :), ...
        'Color', [photonColor, 0.3], 'LineWidth', 0.5);
end

% Plot mean distribution as thick line
plot(photonCenters, meanFreqPhotons, 'r-', 'LineWidth', 2);

xlabel('Photons');
ylabel('Frequency');
title('Individual Cell Distributions (first 10 cells)');
legend({'Individual Cells', 'Mean across cells'}, 'Location', 'best');
grid on;

% Add overall title
annotation(fig1, 'textbox', [0.3, 0.93, 0.4, 0.05], 'String', ...
    'Photon Count Frequency Distribution Analysis', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none');

% Save figure
outputBase1 = fullfile(filePath, 'Frequency_Photon_Distribution');
print([outputBase1 '.pdf'], '-dpdf', '-bestfit');
saveas(gcf, [outputBase1 '.fig']);
print([outputBase1 '.tiff'], '-dtiff', '-r300');
close(fig1);
fprintf('Photon frequency distribution saved as: %s\n', outputBase1);

% ============================================
% 2. LOCALIZATION PRECISION FREQUENCY DISTRIBUTIONS
% ============================================
fprintf('\nCreating localization precision frequency distribution plots...\n');

fig2 = figure('Visible', 'off');
set(gcf, 'Position', [100, 100, 1400, 800]);

% X precision frequency distribution
subplot(2, 2, 1);
hold on;
% Use CUSTOM bar width and color
bar(precisionCentersX, meanFreqX, barWidth, 'FaceColor', xColor, 'EdgeColor', 'none');

% Add SEM error bars for X precision with CUSTOM horizontal line length
numBinsToShowX = min(numPrecisionBinsX, ceil(maxDisplayPrecisionX/precisionBinWidthX));
for i = 1:numBinsToShowX
    if meanFreqX(i) > 0
        xPos = precisionCentersX(i);
        yPos = meanFreqX(i);
        semVal = semFreqX(i);
        
        % Draw error bar
        line([xPos, xPos], [yPos - semVal, yPos + semVal], 'Color', 'k', 'LineWidth', 1);
        % Use CUSTOM error bar line scale
        line([xPos - precisionBinWidthX*errorbarLineScale, xPos + precisionBinWidthX*errorbarLineScale], ...
             [yPos - semVal, yPos - semVal], 'Color', 'k', 'LineWidth', 1);
        line([xPos - precisionBinWidthX*errorbarLineScale, xPos + precisionBinWidthX*errorbarLineScale], ...
             [yPos + semVal, yPos + semVal], 'Color', 'k', 'LineWidth', 1);
    end
end

% Add red vertical line for overall median from combined data
ylimits = ylim;
line([overallMedianX, overallMedianX], [0, ylimits(2)], ...
    'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
text(overallMedianX + precisionBinWidthX*0.5, ylimits(2)*0.9, ...
    sprintf('Median: %.1f nm', overallMedianX), ...
    'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');

xlabel('X Precision (nm)');
ylabel('Mean Frequency � SEM');
title(sprintf('X Precision Distribution (n=%d cells)', numCells));
xlim([0 maxDisplayPrecisionX]);  % Use custom display range
grid on;

% Add bin width annotation
text(0.05, 0.95, sprintf('Bin width: %d nm', precisionBinWidthX), ...
    'Units', 'normalized', 'FontSize', 9, 'BackgroundColor', 'white', 'EdgeColor', 'black');

% Y precision frequency distribution
subplot(2, 2, 2);
hold on;
% Use CUSTOM bar width and color
bar(precisionCentersY, meanFreqY, barWidth, 'FaceColor', yColor, 'EdgeColor', 'none');

% Add SEM error bars for Y precision with CUSTOM horizontal line length
numBinsToShowY = min(numPrecisionBinsY, ceil(maxDisplayPrecisionY/precisionBinWidthY));
for i = 1:numBinsToShowY
    if meanFreqY(i) > 0
        xPos = precisionCentersY(i);
        yPos = meanFreqY(i);
        semVal = semFreqY(i);
        
        % Draw error bar
        line([xPos, xPos], [yPos - semVal, yPos + semVal], 'Color', 'k', 'LineWidth', 1);
        % Use CUSTOM error bar line scale
        line([xPos - precisionBinWidthY*errorbarLineScale, xPos + precisionBinWidthY*errorbarLineScale], ...
             [yPos - semVal, yPos - semVal], 'Color', 'k', 'LineWidth', 1);
        line([xPos - precisionBinWidthY*errorbarLineScale, xPos + precisionBinWidthY*errorbarLineScale], ...
             [yPos + semVal, yPos + semVal], 'Color', 'k', 'LineWidth', 1);
    end
end

% Add red vertical line for overall median from combined data
ylimits = ylim;
line([overallMedianY, overallMedianY], [0, ylimits(2)], ...
    'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
text(overallMedianY + precisionBinWidthY*0.5, ylimits(2)*0.9, ...
    sprintf('Median: %.1f nm', overallMedianY), ...
    'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');

xlabel('Y Precision (nm)');
ylabel('Mean Frequency � SEM');
title(sprintf('Y Precision Distribution (n=%d cells)', numCells));
xlim([0 maxDisplayPrecisionY]);  % Use custom display range
grid on;

% Add bin width annotation
text(0.05, 0.95, sprintf('Bin width: %d nm', precisionBinWidthY), ...
    'Units', 'normalized', 'FontSize', 9, 'BackgroundColor', 'white', 'EdgeColor', 'black');

% Z precision frequency distribution
subplot(2, 2, 3);
hold on;
% Use CUSTOM bar width and color
bar(precisionCentersZ, meanFreqZ, barWidth, 'FaceColor', zColor, 'EdgeColor', 'none');

% Add SEM error bars for Z precision with CUSTOM horizontal line length
numBinsToShowZ = min(numPrecisionBinsZ, ceil(maxDisplayPrecisionZ/precisionBinWidthZ));
for i = 1:numBinsToShowZ
    if meanFreqZ(i) > 0
        xPos = precisionCentersZ(i);
        yPos = meanFreqZ(i);
        semVal = semFreqZ(i);
        
        % Draw error bar
        line([xPos, xPos], [yPos - semVal, yPos + semVal], 'Color', 'k', 'LineWidth', 1);
        % Use CUSTOM error bar line scale
        line([xPos - precisionBinWidthZ*errorbarLineScale, xPos + precisionBinWidthZ*errorbarLineScale], ...
             [yPos - semVal, yPos - semVal], 'Color', 'k', 'LineWidth', 1);
        line([xPos - precisionBinWidthZ*errorbarLineScale, xPos + precisionBinWidthZ*errorbarLineScale], ...
             [yPos + semVal, yPos + semVal], 'Color', 'k', 'LineWidth', 1);
    end
end

% Add red vertical line for overall median from combined data
ylimits = ylim;
line([overallMedianZ, overallMedianZ], [0, ylimits(2)], ...
    'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
text(overallMedianZ + precisionBinWidthZ*0.5, ylimits(2)*0.9, ...
    sprintf('Median: %.1f nm', overallMedianZ), ...
    'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');

xlabel('Z Precision (nm)');
ylabel('Mean Frequency � SEM');
title(sprintf('Z Precision Distribution (n=%d cells)', numCells));
xlim([0 maxDisplayPrecisionZ]);  % Use custom display range
grid on;

% Add bin width annotation
text(0.05, 0.95, sprintf('Bin width: %d nm', precisionBinWidthZ), ...
    'Units', 'normalized', 'FontSize', 9, 'BackgroundColor', 'white', 'EdgeColor', 'black');

% Simplified comparison of peak frequencies
subplot(2, 2, 4);
% Find peak positions for each precision type
[~, peakX] = max(meanFreqX);
[~, peakY] = max(meanFreqY);
[~, peakZ] = max(meanFreqZ);

peakPositions = [precisionCentersX(peakX), precisionCentersY(peakY), precisionCentersZ(peakZ)];
peakFrequencies = [meanFreqX(peakX), meanFreqY(peakY), meanFreqZ(peakZ)];
peakSEMs = [semFreqX(peakX), semFreqY(peakY), semFreqZ(peakZ)];

% Create bars with error bars - simplified for MATLAB 2016a
barPositions = 1:3;
hold on;

% Plot each bar separately with CUSTOM width
h1 = bar(barPositions(1), peakFrequencies(1), barWidth, 'FaceColor', xColor);
h2 = bar(barPositions(2), peakFrequencies(2), barWidth, 'FaceColor', yColor);
h3 = bar(barPositions(3), peakFrequencies(3), barWidth, 'FaceColor', zColor);

% Add error bars manually with CUSTOM horizontal line length
for i = 1:3
    xPos = barPositions(i);
    yPos = peakFrequencies(i);
    semVal = peakSEMs(i);
    
    line([xPos, xPos], [yPos - semVal, yPos + semVal], ...
        'Color', 'k', 'LineWidth', 2);
    % Use CUSTOM error bar line scale (relative to bar width)
    line([xPos - barWidth*errorbarLineScale, xPos + barWidth*errorbarLineScale], ...
         [yPos - semVal, yPos - semVal], 'Color', 'k', 'LineWidth', 2);
    line([xPos - barWidth*errorbarLineScale, xPos + barWidth*errorbarLineScale], ...
         [yPos + semVal, yPos + semVal], 'Color', 'k', 'LineWidth', 2);
end

set(gca, 'XTick', 1:3, 'XTickLabel', {'X', 'Y', 'Z'});
xlabel('Precision Type');
ylabel('Peak Mean Frequency � SEM');
title('Peak Frequency Comparison');
grid on;

% Add legend
legend([h1, h2, h3], {'X Peak', 'Y Peak', 'Z Peak'}, 'Location', 'best');

% Add overall title
annotation(fig2, 'textbox', [0.3, 0.93, 0.4, 0.05], 'String', ...
    'Localization Precision Frequency Distribution Analysis', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', 14, 'FontWeight', 'bold', 'EdgeColor', 'none');

% Save figure
outputBase2 = fullfile(filePath, 'Frequency_Precision_Distribution');
print([outputBase2 '.pdf'], '-dpdf', '-bestfit');
saveas(gcf, [outputBase2 '.fig']);
print([outputBase2 '.tiff'], '-dtiff', '-r300');
close(fig2);
fprintf('Precision frequency distribution saved as: %s\n', outputBase2);

% ============================================
% 3. CREATE SUMMARY TABLES
% ============================================
fprintf('\nCreating frequency distribution summary tables...\n');

% Create table for photon frequency statistics (Mean � SEM)
photonStats = cell(numPhotonBins+1, 5);
photonStats{1, 1} = 'Bin_Center';
photonStats{1, 2} = 'Bin_Start';
photonStats{1, 3} = 'Bin_End';
photonStats{1, 4} = 'Mean_Frequency';
photonStats{1, 5} = 'SEM_Frequency';

for i = 1:numPhotonBins
    photonStats{i+1, 1} = photonCenters(i);
    photonStats{i+1, 2} = photonEdges(i);
    photonStats{i+1, 3} = photonEdges(i+1);
    photonStats{i+1, 4} = meanFreqPhotons(i);
    photonStats{i+1, 5} = semFreqPhotons(i);
end

% Create table for precision frequency statistics (Mean � SEM)
% We'll create separate tables for X, Y, Z or combine them

% X precision statistics
xPrecisionStats = cell(numPrecisionBinsX+1, 5);
xPrecisionStats{1, 1} = 'Bin_Center_nm';
xPrecisionStats{1, 2} = 'Bin_Start_nm';
xPrecisionStats{1, 3} = 'Bin_End_nm';
xPrecisionStats{1, 4} = 'Mean_Freq_X';
xPrecisionStats{1, 5} = 'SEM_Freq_X';

for i = 1:numPrecisionBinsX
    xPrecisionStats{i+1, 1} = precisionCentersX(i);
    xPrecisionStats{i+1, 2} = precisionEdgesX(i);
    xPrecisionStats{i+1, 3} = precisionEdgesX(i+1);
    xPrecisionStats{i+1, 4} = meanFreqX(i);
    xPrecisionStats{i+1, 5} = semFreqX(i);
end

% Y precision statistics
yPrecisionStats = cell(numPrecisionBinsY+1, 5);
yPrecisionStats{1, 1} = 'Bin_Center_nm';
yPrecisionStats{1, 2} = 'Bin_Start_nm';
yPrecisionStats{1, 3} = 'Bin_End_nm';
yPrecisionStats{1, 4} = 'Mean_Freq_Y';
yPrecisionStats{1, 5} = 'SEM_Freq_Y';

for i = 1:numPrecisionBinsY
    yPrecisionStats{i+1, 1} = precisionCentersY(i);
    yPrecisionStats{i+1, 2} = precisionEdgesY(i);
    yPrecisionStats{i+1, 3} = precisionEdgesY(i+1);
    yPrecisionStats{i+1, 4} = meanFreqY(i);
    yPrecisionStats{i+1, 5} = semFreqY(i);
end

% Z precision statistics
zPrecisionStats = cell(numPrecisionBinsZ+1, 5);
zPrecisionStats{1, 1} = 'Bin_Center_nm';
zPrecisionStats{1, 2} = 'Bin_Start_nm';
zPrecisionStats{1, 3} = 'Bin_End_nm';
zPrecisionStats{1, 4} = 'Mean_Freq_Z';
zPrecisionStats{1, 5} = 'SEM_Freq_Z';

for i = 1:numPrecisionBinsZ
    zPrecisionStats{i+1, 1} = precisionCentersZ(i);
    zPrecisionStats{i+1, 2} = precisionEdgesZ(i);
    zPrecisionStats{i+1, 3} = precisionEdgesZ(i+1);
    zPrecisionStats{i+1, 4} = meanFreqZ(i);
    zPrecisionStats{i+1, 5} = semFreqZ(i);
end

% Save tables
photonStatsPath = fullfile(filePath, 'Frequency_Photon_Statistics.csv');
xPrecisionStatsPath = fullfile(filePath, 'Frequency_X_Precision_Statistics.csv');
yPrecisionStatsPath = fullfile(filePath, 'Frequency_Y_Precision_Statistics.csv');
zPrecisionStatsPath = fullfile(filePath, 'Frequency_Z_Precision_Statistics.csv');

% Write photon statistics
fid = fopen(photonStatsPath, 'w');
fprintf(fid, '%s,%s,%s,%s,%s\n', photonStats{1, :});
for i = 2:size(photonStats, 1)
    fprintf(fid, '%.0f,%.0f,%.0f,%.2f,%.2f\n', photonStats{i, :});
end
fclose(fid);

% Write X precision statistics
fid = fopen(xPrecisionStatsPath, 'w');
fprintf(fid, '%s,%s,%s,%s,%s\n', xPrecisionStats{1, :});
for i = 2:size(xPrecisionStats, 1)
    fprintf(fid, '%.1f,%.1f,%.1f,%.2f,%.2f\n', xPrecisionStats{i, :});
end
fclose(fid);

% Write Y precision statistics
fid = fopen(yPrecisionStatsPath, 'w');
fprintf(fid, '%s,%s,%s,%s,%s\n', yPrecisionStats{1, :});
for i = 2:size(yPrecisionStats, 1)
    fprintf(fid, '%.1f,%.1f,%.1f,%.2f,%.2f\n', yPrecisionStats{i, :});
end
fclose(fid);

% Write Z precision statistics
fid = fopen(zPrecisionStatsPath, 'w');
fprintf(fid, '%s,%s,%s,%s,%s\n', zPrecisionStats{1, :});
for i = 2:size(zPrecisionStats, 1)
    fprintf(fid, '%.1f,%.1f,%.1f,%.2f,%.2f\n', zPrecisionStats{i, :});
end
fclose(fid);

fprintf('Photon frequency statistics saved to: %s\n', photonStatsPath);
fprintf('X precision frequency statistics saved to: %s\n', xPrecisionStatsPath);
fprintf('Y precision frequency statistics saved to: %s\n', yPrecisionStatsPath);
fprintf('Z precision frequency statistics saved to: %s\n', zPrecisionStatsPath);

% Display final summary
fprintf('\n=== FINAL FREQUENCY DISTRIBUTION ANALYSIS ===\n');
fprintf('Number of cells analyzed: %d\n', numCells);
fprintf('Photon bins: %d (range: %.0f-%.0f, width: %d photons)\n', numPhotonBins, photonMin, photonMax, photonBinWidth);
fprintf('X precision bins: %d (range: %.1f-%.1f nm, width: %d nm)\n', numPrecisionBinsX, precisionMinX, precisionMaxX, precisionBinWidthX);
fprintf('Y precision bins: %d (range: %.1f-%.1f nm, width: %d nm)\n', numPrecisionBinsY, precisionMinY, precisionMaxY, precisionBinWidthY);
fprintf('Z precision bins: %d (range: %.1f-%.1f nm, width: %d nm)\n', numPrecisionBinsZ, precisionMinZ, precisionMaxZ, precisionBinWidthZ);

% Find and report peak frequencies (mean)
[peakFreqPhotons, peakIdxPhotons] = max(meanFreqPhotons);
peakPhotonBin = photonCenters(peakIdxPhotons);

fprintf('\nPeak photon frequency: %.0f � %.0f at %.0f photons\n', ...
    peakFreqPhotons, semFreqPhotons(peakIdxPhotons), peakPhotonBin);
fprintf('Overall median (combined data): %.1f photons\n', overallMedianPhotons);
fprintf('Overall medians for precision (combined data):\n');
fprintf('  X: %.2f nm, Y: %.2f nm, Z: %.2f nm\n', overallMedianX, overallMedianY, overallMedianZ);

% Display user parameters summary
fprintf('\n=== USER PARAMETERS USED ===\n');
fprintf('Error bar line scale: %.2f (horizontal lines = %.0f%% of bin width)\n', errorbarLineScale, errorbarLineScale*100);
fprintf('Bar width in plots: %.2f\n', barWidth);
fprintf('Display ranges: X=0-%d nm, Y=0-%d nm, Z=0-%d nm\n', maxDisplayPrecisionX, maxDisplayPrecisionY, maxDisplayPrecisionZ);

disp('All frequency distribution analyses completed successfully!');