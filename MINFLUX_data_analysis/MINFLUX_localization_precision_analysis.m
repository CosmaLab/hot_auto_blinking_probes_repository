% Script Description:
% This script is compatible with MATLAB 2018b.
% This script processes all .mat files in the specified folder by:
% 1. Applying detailed filtering conditions to MINFLUX data for each cell
% 2. Calculating localization precision for each cell separately
% 3. Combining data from all cells and calculating mean ± SEM
% 4. Generating histograms with median precision indicators
% 5. Saving results to 'localization_precisions.csv'
% 6. Saving figures in PDF, TIFF, and FIG formats

% Initialize variables
resultsTable = table([], [], [], [],...
    'VariableNames', {'file_name', 'loc_precision_X', 'loc_precision_Y', 'loc_precision_Z'});

% Configure paths
dataFolder = 'C:\Users\18317\Desktop\code 20260225\Euchromatin_5cells';
matFiles = dir(fullfile(dataFolder, '*.mat'));

% Define colors
faceColor = [0.2, 0.6, 0.8];    
edgeColor = 'k';       

% Initialize cell arrays to store precision data from all cells
all_X_precisions = [];
all_Y_precisions = [];
all_Z_precisions = [];

% Cell arrays to store individual cell data for statistics
cell_X_data = {};
cell_Y_data = {};
cell_Z_data = {};

% Main processing loop - First pass: collect data from all cells
for k = 1:length(matFiles)
    fileName = matFiles(k).name;
    fprintf('Processing file: %s\n', fileName);
    data = load(fullfile(dataFolder, fileName));
    
    %% Extract and Process MINFLUX Localizations - Apply Filtering Conditions
    fprintf('Extracting and filtering MINFLUX localizations...\n');
    
    % Step 1: Get initial valid localizations (data.vld == 1)
    initial_valid = data.vld == 1;
    initial_count = sum(initial_valid);
    
    % Extract all localization data initially
    all_locations = squeeze(data.itr.loc(:, end, :)); % Last iteration positions
    all_traceIDs = data.tid;
    
    % Apply initial valid filter
    locations_step1 = all_locations(initial_valid, :);
    traceIDs_step1 = all_traceIDs(initial_valid);
    
    % Step 2: Filter by cfr (itr.cfr last column <= 0.8)
    if isfield(data.itr, 'cfr')
        cfr_data = data.itr.cfr(:, end);
        cfr_valid = cfr_data <= 0.8;
        
        all_indices = 1:length(data.vld);
        valid_indices = all_indices(initial_valid);
        cfr_valid_for_step1 = cfr_valid(valid_indices);
        
        locations_step2 = locations_step1(cfr_valid_for_step1, :);
        traceIDs_step2 = traceIDs_step1(cfr_valid_for_step1);
    else
        locations_step2 = locations_step1;
        traceIDs_step2 = traceIDs_step1;
    end
    
    % Step 3: Filter by efo (itr.efo last column <= 70000)
    if isfield(data.itr, 'efo')
        efo_data = data.itr.efo(:, end);
        efo_valid = efo_data <= 70000;
        
        all_indices = 1:length(data.vld);
        initial_valid_indices = all_indices(initial_valid);
        
        if exist('cfr_valid_for_step1', 'var')
            step2_indices = initial_valid_indices(cfr_valid_for_step1);
        else
            step2_indices = initial_valid_indices;
        end
        
        efo_valid_for_current = efo_valid(step2_indices);
        locations_step3 = locations_step2(efo_valid_for_current, :);
        traceIDs_step3 = traceIDs_step2(efo_valid_for_current);
    else
        locations_step3 = locations_step2;
        traceIDs_step3 = traceIDs_step2;
    end
    
    % Step 4: Filter traces with at least 2 points
    unique_traceIDs_step3 = unique(traceIDs_step3);
    trace_counts = zeros(length(unique_traceIDs_step3), 1);
    
    for i = 1:length(unique_traceIDs_step3)
        trace_counts(i) = sum(traceIDs_step3 == unique_traceIDs_step3(i));
    end
    
    valid_traces = unique_traceIDs_step3(trace_counts >= 2);
    trace_filter = ismember(traceIDs_step3, valid_traces);
    
    % Apply trace filter
    locations_final = locations_step3(trace_filter, :);
    traceIDs_final = traceIDs_step3(trace_filter);
    
    % Calculate statistics for this cell
    [~, ~, newGroupIdx] = unique(traceIDs_final);
    
    % Calculate standard deviation for each trace
    traceSTD = splitapply(@(x) std(x,0,1), locations_final, newGroupIdx);
    
    % Store individual cell data (convert to nm)
    cell_X_data{end+1} = traceSTD(:,1) * 1e9;
    cell_Y_data{end+1} = traceSTD(:,2) * 1e9;
    if size(traceSTD, 2) > 2
        cell_Z_data{end+1} = traceSTD(:,3) * 1e9;
    end
    
    % Compute median precision (convert to nm)
    precisionX = median(traceSTD(:,1)) * 1e9;
    precisionY = median(traceSTD(:,2)) * 1e9;
    if size(traceSTD, 2) > 2
        precisionZ = median(traceSTD(:,3)) * 1e9;
    else
        precisionZ = 0;
    end
    
    % Update results table
    newRow = table(string(fileName), precisionX, precisionY, precisionZ,...
        'VariableNames', resultsTable.Properties.VariableNames);
    resultsTable = [resultsTable; newRow];
    
    % Store all precisions for combined analysis
    all_X_precisions = [all_X_precisions; traceSTD(:,1) * 1e9];
    all_Y_precisions = [all_Y_precisions; traceSTD(:,2) * 1e9];
    if size(traceSTD, 2) > 2
        all_Z_precisions = [all_Z_precisions; traceSTD(:,3) * 1e9];
    end
    
    % Generate individual cell visualization (optional)
    % Comment out if you don't need individual cell plots
    generateIndividualCellPlot(fileName, traceSTD, precisionX, precisionY, precisionZ, dataFolder);
end

% Save results table
outputPath = fullfile(dataFolder, 'localization_precisions.csv');
writetable(resultsTable, outputPath);

% Now create combined histograms with mean ± SEM
% Check if we have 3D data
has3D = ~isempty(all_Z_precisions) && any(all_Z_precisions ~= 0);

% Create combined figure
if has3D
    fig = figure('Name', 'Combined Localization Precision (All Cells)', 'NumberTitle', 'off',...
        'Color', 'w', 'Position', [100, 100, 1500, 500],...
        'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
    
    % X axis
    subplot(1,3,1);
    createCombinedHistogram(cell_X_data, 'X', all_X_precisions, faceColor);
    ax1 = gca;
    
    % Y axis
    subplot(1,3,2);
    createCombinedHistogram(cell_Y_data, 'Y', all_Y_precisions, faceColor);
    ax2 = gca;
    
    % Z axis
    subplot(1,3,3);
    createCombinedHistogram(cell_Z_data, 'Z', all_Z_precisions, faceColor);
    ax3 = gca;
else
    fig = figure('Name', 'Combined Localization Precision (All Cells)', 'NumberTitle', 'off',...
        'Color', 'w', 'Position', [100, 100, 1200, 500],...
        'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
    
    % X axis
    subplot(1,2,1);
    createCombinedHistogram(cell_X_data, 'X', all_X_precisions, faceColor);
    ax1 = gca;
    
    % Y axis
    subplot(1,2,2);
    createCombinedHistogram(cell_Y_data, 'Y', all_Y_precisions, faceColor);
    ax2 = gca;
end

% Save combined figure
savePath = fullfile(dataFolder, 'Combined_Localization_Precision');

% Save as FIG
saveas(fig, [savePath '.fig']);

% Save as PDF (vector format)
print(fig, [savePath '.pdf'], '-dpdf', '-painters', '-r300')

% Save as TIFF (300 dpi)
print(fig, [savePath '.tiff'], '-dtiff', '-r300')

% Save each subplot as individual PDF
if has3D
    saveIndividualSubplotPDF(ax1, dataFolder, 'Combined_X_Axis_Precision', 'X');
    saveIndividualSubplotPDF(ax2, dataFolder, 'Combined_Y_Axis_Precision', 'Y');
    saveIndividualSubplotPDF(ax3, dataFolder, 'Combined_Z_Axis_Precision', 'Z');
else
    saveIndividualSubplotPDF(ax1, dataFolder, 'Combined_X_Axis_Precision', 'X');
    saveIndividualSubplotPDF(ax2, dataFolder, 'Combined_Y_Axis_Precision', 'Y');
end

% Also save individual cell data for further analysis
save(fullfile(dataFolder, 'cell_precision_data.mat'),...
    'cell_X_data', 'cell_Y_data', 'cell_Z_data',...
    'all_X_precisions', 'all_Y_precisions', 'all_Z_precisions');

fprintf('Processing complete. Results saved to: %s\n', dataFolder);

%% ========== Helper Functions ==========

function createCombinedHistogram(cell_data, axisLabel, all_data, faceColor)
    % Create histogram of localization precision (std) for multiple cells
    % Y-axis shows mean frequency count per bin ± SEM across cells
    
    % Determine bin edges based on all data
    all_combined = [];
    for i = 1:length(cell_data)
        all_combined = [all_combined; cell_data{i}];
    end
    
    % Set bin parameters
    binWidth = 0.3;  % in nm
    maxValue = min(20, prctile(all_combined, 95)); % Cap at 95th percentile or 20nm
    binEdges = 0:binWidth:maxValue;
    binCenters = binEdges(1:end-1) + binWidth/2;
    
    % Initialize matrix to store histograms from all cells
    nCells = length(cell_data);
    nBins = length(binCenters);
    cell_histograms = zeros(nCells, nBins);
    
    % Calculate histogram for each cell (frequency counts, not probability density)
    for i = 1:nCells
        if ~isempty(cell_data{i})
            counts = histcounts(cell_data{i}, binEdges, 'Normalization', 'count');
            cell_histograms(i, :) = counts;
        end
    end
    
    % Calculate mean and SEM for each bin
    mean_freq = mean(cell_histograms, 1);
    sem_freq = std(cell_histograms, 0, 1) / sqrt(nCells);
    
    % Create bar plot with error bars
    hBar = bar(binCenters, mean_freq, 'FaceColor', faceColor, 'EdgeColor', 'k',...
        'FaceAlpha', 0.7, 'BarWidth', 0.8);
    hold on;
    
    % Add error bars (SEM)
    errorbar(binCenters, mean_freq, sem_freq, sem_freq,...
        'LineStyle', 'none', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
    
    % Calculate and plot median of all combined data
    if ~isempty(all_data)
        median_all = median(all_data);
        
        % Plot red dashed line for median
        xline(median_all, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 2,...
            'DisplayName', sprintf('Median: %.2f nm', median_all));
        
        % Add text label for median value
        y_lim = ylim;
        text_y = y_lim(2) * 0.95;  % Place at 95% of the y-axis height
        text(median_all, text_y, sprintf('Median: %.2f nm', median_all),...
            'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'top',...
            'FontSize', 10,...
            'Color', 'r',...
            'FontWeight', 'bold',...
            'BackgroundColor', [1 1 1 0.7]);  % Semi-transparent white background
    end
    
    % Set plot properties
    title(sprintf('%s Axis Precision (All Cells, n=%d)', axisLabel, nCells), 'FontSize', 12);
    xlabel('Standard Deviation (nm)', 'FontSize', 11);
    ylabel('Frequency', 'FontSize', 11);   % Y-axis label changed to Frequency
    xlim([0, maxValue]);
    ylim([0, max(mean_freq + sem_freq) * 1.2]);
    grid off;
    set(gca, 'FontSize', 10, 'LineWidth', 1, 'Box', 'on');
    
    % Add legend
    if ~isempty(all_data)
        legend('Mean frequency', 'SEM', 'Median (all data)',...
            'Location', 'northeast');
    else
        legend('Mean frequency', 'SEM', 'Location', 'northeast');
    end
    
    hold off;
end

% Helper function for individual cell plots (optional)
function generateIndividualCellPlot(fileName, traceSTD, precisionX, precisionY, precisionZ, dataFolder)
    binSize = 0.3;  % in nm
    faceColor = [0.2, 0.6, 0.8];    
    edgeColor = 'k';       
    
    fig = figure('Name', fileName, 'NumberTitle', 'off',...
        'Color', 'w', 'Visible', 'off',...
        'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
    
    if precisionZ == 0
        % 2D data visualization
        subplot(1,2,1);
        createHistogram(traceSTD(:,1)*1e9, binSize, 'X', precisionX, faceColor, edgeColor);
        ax1 = gca;
        
        subplot(1,2,2);
        createHistogram(traceSTD(:,2)*1e9, binSize, 'Y', precisionY, faceColor, edgeColor);
        ax2 = gca;
    else
        % 3D data visualization
        subplot(1,3,1);
        createHistogram(traceSTD(:,1)*1e9, binSize, 'X', precisionX, faceColor, edgeColor);
        ax1 = gca;
        
        subplot(1,3,2);
        createHistogram(traceSTD(:,2)*1e9, binSize, 'Y', precisionY, faceColor, edgeColor);
        ax2 = gca;
        
        subplot(1,3,3);
        createHistogram(traceSTD(:,3)*1e9, binSize, 'Z', precisionZ, faceColor, edgeColor);
        ax3 = gca;
    end
    
    % Save individual figures
    [~, baseName] = fileparts(fileName);
    savePath = fullfile(dataFolder, 'Individual_Cells', baseName);
    
    % Create directory if it doesn't exist
    if ~exist(fullfile(dataFolder, 'Individual_Cells'), 'dir')
        mkdir(fullfile(dataFolder, 'Individual_Cells'));
    end
    
    % Save as FIG
    saveas(fig, [savePath '.fig']);
    
    % Save as PDF (vector format)
    print(fig, [savePath '.pdf'], '-dpdf', '-painters', '-r300')
    
    % Save as TIFF (300 dpi)
    print(fig, [savePath '.tiff'], '-dtiff', '-r300')
    
    % Save each subplot as individual PDF
    if precisionZ == 0
        saveIndividualSubplotPDF(ax1, fullfile(dataFolder, 'Individual_Cells'), [baseName '_X'], 'X');
        saveIndividualSubplotPDF(ax2, fullfile(dataFolder, 'Individual_Cells'), [baseName '_Y'], 'Y');
    else
        saveIndividualSubplotPDF(ax1, fullfile(dataFolder, 'Individual_Cells'), [baseName '_X'], 'X');
        saveIndividualSubplotPDF(ax2, fullfile(dataFolder, 'Individual_Cells'), [baseName '_Y'], 'Y');
        saveIndividualSubplotPDF(ax3, fullfile(dataFolder, 'Individual_Cells'), [baseName '_Z'], 'Z');
    end
    
    close(fig);
end

% Helper function for histogram creation (for individual cells)
function createHistogram(data, binWidth, axisLabel, medianValue, faceColor, edgeColor)
    histogram(data,...
        'BinWidth', binWidth,...
        'FaceColor', faceColor,...
        'EdgeColor', edgeColor,...
        'LineWidth', 0.5);
    
    hold on;
    xline(medianValue,...
        'LineStyle', '--',...
        'Color', 'r',...
        'LineWidth', 2);
    
    title(sprintf('%s Axis Precision: %.2f nm', axisLabel, medianValue));
    xlabel('Standard Deviation (nm)');
    ylabel('Frequency');   % Already 'Frequency' for individual histograms
    xlim([0, max(20, prctile(data, 95))]);
    grid off;
    set(gca, 'FontSize', 10, 'LineWidth', 1);
end

% New function to save individual subplot as PDF
function saveIndividualSubplotPDF(ax, folder, filename, axisLabel)
    % Create a new figure with the same properties as the subplot
    fig_single = figure('Visible', 'off', 'Color', 'w',...
        'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
    
    % Copy the axes content to the new figure
    new_ax = copyobj(ax, fig_single);
    
    % Adjust the axes position to fill the figure
    set(new_ax, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.75, 0.75]);
    
    % Set the figure size to match the aspect ratio of the plot
    set(fig_single, 'Units', 'inches', 'Position', [1, 1, 6, 5]);
    
    % Save the single plot as PDF
    savePath = fullfile(folder, filename);
    print(fig_single, [savePath '.pdf'], '-dpdf', '-painters', '-r300');
    
    % Also save as PNG for quick viewing
    print(fig_single, [savePath '.png'], '-dpng', '-r300');
    
    % Close the temporary figure
    close(fig_single);
    
    fprintf('Saved individual %s axis plot: %s.pdf\n', axisLabel, filename);
end