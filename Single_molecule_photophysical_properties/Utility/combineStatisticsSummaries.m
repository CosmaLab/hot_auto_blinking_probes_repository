function combineStatisticsSummaries2()
    % The script will:
    %   1. load from "molecule_statistics_summary.csv":
    %      - Variable: Name of the metric
    %        * MeanOnTime_ms: Average ON-state duration in milliseconds
    %        * MeanOffTime_ms: Average OFF-state duration in milliseconds
    %        * DutyCycle: Fraction of time molecule is in ON state
    %        * ActiveDutyCycle: Duty cycle during active periods
    %        * BlinkingEvents: Number of ON-OFF transitions
    %        * PhotonsPerCycle: Average photon count per ON-OFF cycle
    %        * PhotonsPerDetection: Average photon count per detection event
    %        * DutyCycle_t1, DutyCycle_t2, etc.: Time-windowed duty cycles
    %        * DC_100s_t100: 100-second window duty cycle at t=100s
    %      - Mean: Average value across all molecules
    %      - Std: Standard deviation across molecules
    %      - Min: Minimum value
    %      - Max: Maximum value
    %      - Median: Median value
    %      - Count: Number of molecules used in calculation
    %   2. Combine the molecule_statistics_summary files into one consolidated file
    %   3. Save the combined data as both .mat and .csv files
    
    %% load statistics
        % Ask user to select the base directory
        baseDir = uigetdir('', 'Select the base directory containing analysis results');
        if baseDir == 0
            disp('Operation cancelled by user.');
            return;
        end
        
        % Find all analysis run directories, including in subdirectories
        allDirs = [];
        
        % Recursively search for Run directories in current and subdirectories
        function dirs = findRunDirs(startPath)
            % Find only BCA analysis directories
            dirs = dir(fullfile(startPath, '*_BCA_*_Run_*'));
            dirs = dirs([dirs.isdir]);
            
            % Search in subdirectories
            subDirs = dir(startPath);
            subDirs = subDirs([subDirs.isdir]);
            
            % Remove '.' and '..' directories
            subDirs = subDirs(~ismember({subDirs.name}, {'.', '..'}));
            
            % Recursively search each subdirectory
            for i = 1:length(subDirs)
                if ~contains(subDirs(i).name, '_Run_')  % Skip if it's already a Run directory
                    newDirs = findRunDirs(fullfile(subDirs(i).folder, subDirs(i).name));
                    dirs = [dirs; newDirs];
                end
            end
        end
    
        % Find all Run directories
        allDirs = findRunDirs(baseDir);
        
        % Convert to cell array of names with full paths and extract timestamps
        runDirs = cell(length(allDirs), 1);
        runDates = NaT(length(allDirs), 1);  % Initialize datetime array
        
        for i = 1:length(allDirs)
            runDirs{i} = fullfile(allDirs(i).folder, allDirs(i).name);
            % Extract timestamp from directory name
            [~, dirName] = fileparts(runDirs{i});
            timestamp = regexp(dirName, '_Run_(\d{8}_\d{6})', 'tokens');
            if ~isempty(timestamp)
                runDates(i) = datetime(timestamp{1}{1}, 'InputFormat', 'yyyyMMdd_HHmmss');
            end
        end
        
        % Sort directories by timestamp
        [runDates, sortIdx] = sort(runDates);
        runDirs = runDirs(sortIdx);
        
        if isempty(runDirs)
            disp('No BCA analysis run directories found.');
            return;
        end
    
        % Create display names with timestamps (using sorted indices)
        displayNames = cell(length(runDirs), 1);
        for i = 1:length(runDirs)
            % Format the timestamp for display
            runDateStr = datestr(runDates(i), 'yyyy-mm-dd HH:MM:SS');
            
            % Extract analysis type and method
            [~, dirName] = fileparts(runDirs{i});
            
            % Determine analysis method
            if contains(dirName, '_MSK')
                analysisMethod = 'fit-based';
                maskParams = regexp(dirName, '_MSK(\d+)_THR(\d+)', 'tokens');
                if ~isempty(maskParams)
                    analysisMethod = sprintf('%s (MSK=%s,THR=%s)', ...
                        analysisMethod, maskParams{1}{1}, maskParams{1}{2});
                end
            elseif contains(dirName, '_HFSQR')
                analysisMethod = 'image-based';
                imgParams = regexp(dirName, '_HFSQR(\d+)_THR(\d+)', 'tokens');
                if ~isempty(imgParams)
                    analysisMethod = sprintf('%s (HFSQR=%s,THR=%s)', ...
                        analysisMethod, imgParams{1}{1}, imgParams{1}{2});
                end
            else
                analysisMethod = 'unknown';
            end
            
            % Include timestamp and analysis info in the display name
            displayNames{i} = sprintf('[%s] %s (%s-%s)', ...
                runDateStr, runDirs{i}, analysisMethod);
        end
        
        % Display chronological summary before selection
        fprintf('\nFound %d analysis runs in chronological order:\n', length(runDirs));
        for i = 1:length(runDirs)
            fprintf('%d. %s\n', i, displayNames{i});
        end
        fprintf('\n');
        
        % Let user select which runs to include
        [selectedIndices, ok] = listdlg('ListString', displayNames, ...
            'SelectionMode', 'multiple', ...
            'Name', 'Select Analysis Runs', ...
            'PromptString', sprintf('Select runs to combine (%d runs found):', length(runDirs)), ...
            'ListSize', [800 400]);  % Made wider to accommodate timestamps
        
        if ~ok || isempty(selectedIndices)
            disp('No runs selected. Operation cancelled.');
            return;
        end
        
        %% Initialize data structure for combined results
        combinedData = struct();
        combinedData.runInfo = cell(length(selectedIndices), 1);
        combinedData.statistics = cell(length(selectedIndices), 1);
        combinedTable = table();
        
        %% Process each selected run
        for i = 1:length(selectedIndices)
            runIdx = selectedIndices(i);
            runDir = runDirs{runIdx};  % Already contains full path
            statsDir = fullfile(runDir, '5_Statistics');
            
            % Check if statistics directory exists
            if ~exist(statsDir, 'dir')
                fprintf('Warning: Statistics directory not found for %s. Skipping.\n', runDirs{runIdx});
                continue;
            end
            
            %% Look for molecule_statistics_summary file
            summaryFile = fullfile(statsDir, 'molecule_statistics_summary.csv');
            
            if exist(summaryFile, 'file')
                % Read the CSV file
                statsTable = readtable(summaryFile);
                
                % Find all duty cycle columns (including windowed duty cycles)
                dcColumns = contains(statsTable.Variable, 'DutyCycle');
                dcWindowColumns = {};
                
                % Find windowed duty cycle columns with the new format (DutyCycle_t1, DutyCycle_t2, etc.)
                for j = 1:height(statsTable)
                    varName = statsTable.Variable{j};
                    if startsWith(varName, 'DutyCycle_t') && ~strcmp(varName, 'DutyCycle')
                        dcWindowColumns{end+1} = varName;
                    elseif startsWith(varName, 'DC_') 
                        % For backward compatibility with older formats
                        dcWindowColumns{end+1} = varName;
                    end
                end
                
                fprintf('Found %d regular duty cycle metrics and %d windowed duty cycle metrics in %s\n', ...
                    sum(dcColumns), length(dcWindowColumns), runDirs{runIdx});
                
                % Store run information
                runInfo = struct();
                runInfo.name = runDirs{runIdx};
                runInfo.displayName = displayNames{runIdx};
                
                % Extract timestamp
                timestamp = regexp(runDirs{runIdx}, '_Run_(\d{8}_\d{6})', 'tokens');
                if ~isempty(timestamp)
                    runInfo.timestamp = timestamp{1}{1};
                    runInfo.datetime = datetime(timestamp{1}{1}, 'InputFormat', 'yyyyMMdd_HHmmss');
                else
                    runInfo.timestamp = '';
                    runInfo.datetime = datetime('now');
                end
                
                % Extract analysis method (fit-based or image-based)
                [~, dirName] = fileparts(runDirs{runIdx});
                if contains(dirName, '_MSK')
                    runInfo.analysisMethod = 'fit_based';
                    maskParams = regexp(dirName, '_MSK(\d+)_THR(\d+)', 'tokens');
                    if ~isempty(maskParams)
                        runInfo.maskSize = str2double(maskParams{1}{1});
                        runInfo.threshold = str2double(maskParams{1}{2});
                    end
                elseif contains(dirName, '_HFSQR')
                    runInfo.analysisMethod = 'image_based';
                    imgParams = regexp(dirName, '_HFSQR(\d+)_THR(\d+)', 'tokens');
                    if ~isempty(imgParams)
                        runInfo.halfSquare = str2double(imgParams{1}{1});
                        runInfo.threshold = str2double(imgParams{1}{2});
                    end
                else
                    runInfo.analysisMethod = 'unknown';
                end
                
                % Add run information to the table
                statsTable.RunName = repmat({runInfo.name}, height(statsTable), 1);
                statsTable.Timestamp = repmat({runInfo.timestamp}, height(statsTable), 1);
                statsTable.AnalysisMethod = repmat({runInfo.analysisMethod}, height(statsTable), 1);
                
                % Add to combined table
                if isempty(combinedTable)
                    combinedTable = statsTable;
                else
                    % Ensure all columns exist in both tables
                    existingCols = combinedTable.Properties.VariableNames;
                    newCols = statsTable.Properties.VariableNames;
                    
                    % Add missing columns to combinedTable
                    for j = 1:length(newCols)
                        if ~ismember(newCols{j}, existingCols)
                            combinedTable.(newCols{j}) = NaN(height(combinedTable), 1);
                            if iscell(statsTable.(newCols{j}))
                                combinedTable.(newCols{j}) = cell(height(combinedTable), 1);
                            end
                        end
                    end
                    
                    % Add missing columns to statsTable
                    for j = 1:length(existingCols)
                        if ~ismember(existingCols{j}, newCols)
                            statsTable.(existingCols{j}) = NaN(height(statsTable), 1);
                            if iscell(combinedTable.(existingCols{j}))
                                statsTable.(existingCols{j}) = cell(height(statsTable), 1);
                            end
                        end
                    end
                    
                    % Now combine the tables
                    combinedTable = [combinedTable; statsTable];
                end
                
                fprintf('Added statistics from %s\n', runDirs{runIdx});
            else
                fprintf('Warning: No statistics summary file found for %s\n', runDirs{runIdx});
                continue;
            end
        end
        
        %% STEP 5: VALIDATE COMBINED DATA
        % Check if we have any data after processing all runs
        if isempty(combinedTable)
            disp('No valid statistics data found in the selected runs.');
            return;
        end
        
        %% STEP 6: PREPARE OUTPUT DIRECTORY
        % Extract analysis parameters from the first selected run for folder naming
        firstRunIdx = selectedIndices(1);
        [~, firstDirName] = fileparts(runDirs{firstRunIdx});
        folderAnalysisTag = '';
        
        if contains(firstDirName, '_MSK')
            maskParams = regexp(firstDirName, '_MSK(\d+)_THR(\d+)', 'tokens');
            if ~isempty(maskParams)
                folderAnalysisTag = sprintf('_MSK%s_THR%s', maskParams{1}{1}, maskParams{1}{2});
            end
        elseif contains(firstDirName, '_HFSQR')
            imgParams = regexp(firstDirName, '_HFSQR(\d+)_THR(\d+)', 'tokens');
            if ~isempty(imgParams)
                folderAnalysisTag = sprintf('_HFSQR%s_THR%s', imgParams{1}{1}, imgParams{1}{2});
            end
        end
        
        % Create output directory with analysis tag
        outputDir = fullfile(baseDir, sprintf('Combined_Statistics%s', folderAnalysisTag));
        if ~exist(outputDir, 'dir')
            mkdir(outputDir);
            fprintf('Created output directory: %s\n', outputDir);
        end
    
        %% STEP 7: EXTRACT SAMPLE METADATA
        % Extract dye type and sample indices
        dyePatterns = {'2CF3', 'CF3', 'SO2', 'C4','Sample 3', 'Sample 4', 'Sample 5', 'Sample 6','Sample 7'};  % Order matters: check 2CF3 before CF3
        sampleIndices = cell(height(combinedTable), 1);  % Change to cell array to handle both numeric and string identifiers
        dyeTypes = {};
        
        for i = 1:height(combinedTable)
            runName = combinedTable.RunName{i};
            
            % Extract SR number
            srMatch = regexp(runName, 'SR(\d+)', 'tokens');
            if ~isempty(srMatch)
                sampleIndices{i} = str2double(srMatch{1}{1});
                fprintf('Found SR number %s in run: %s\n', srMatch{1}{1}, runName);
            else
                % If no SR number, try to extract sample identifier from run name
                % Extract the part before _BCA from paths like:
                % "...\Analysis_Results\Sample 3_561_20mW_Exp_20ms_2_BCA_HFSQR2_THR10_Run_20250519_151714"
                [~, fileName] = fileparts(runName);
                bcaMatch = regexp(fileName, '(.+?)_BCA_', 'tokens');
                
                if ~isempty(bcaMatch)
                    % Use the sample identifier as a string
                    sampleId = bcaMatch{1}{1};
                    fprintf('No SR number found, using sample identifier: %s\n', sampleId);
                    sampleIndices{i} = sampleId;
                else
                    % Try alternative pattern for older naming conventions
                    sampleIdMatch = regexp(runName, '([^_\\\/]+)_Run_\d{8}_\d{6}', 'tokens');
                    if ~isempty(sampleIdMatch)
                        % Use the sample identifier as a string
                        sampleId = sampleIdMatch{1}{1};
                        fprintf('No SR number found, using sample identifier: %s\n', sampleId);
                        sampleIndices{i} = sampleId;
                    else
                        fprintf('Warning: No sample identifier found in run name: %s\n', runName);
                        sampleIndices{i} = sprintf('Unknown_%d', i);
                    end
                end
            end
            
            % Extract dye type from folder structure
            folderParts = strsplit(runName, filesep);
            dyeFound = false;
            
            % Check each folder level for dye pattern
            for j = 1:length(dyePatterns)
                pattern = dyePatterns{j};
                folderMatch = find(contains(folderParts, pattern), 1);
                
                if ~isempty(folderMatch)
                    dyeTypes{end+1} = pattern;
                    dyeFound = true;
                    fprintf('Found dye %s in folder structure: %s\n', pattern, runName);
                    break;
                end
            end
            
            if ~dyeFound
                fprintf('Warning: No dye type found in folder structure for run: %s\n', runName);
                dyeTypes{end+1} = 'Unknown';
            end
        end
        
        %% STEP 8: PREPARE FOR STATISTICS CALCULATION
        % Get all unique samples across all dye types
        allUniqueSamples = unique(sampleIndices);
        
        % Initialize base variable names with only found samples
        baseVarNames = {'Metric', 'Dye', 'All_Samples'};
        for j = 1:length(allUniqueSamples)
            if isnumeric(allUniqueSamples{j})
                baseVarNames{end+1} = sprintf('SR%d', allUniqueSamples{j});
            else
                % For string-based sample identifiers, create a valid variable name
                sampleId = regexprep(allUniqueSamples{j}, '[^a-zA-Z0-9]', '_');
                baseVarNames{end+1} = sprintf('Sample_%s', sampleId);
            end
        end
        
        %% STEP 9: CALCULATE STATISTICS FOR EACH DYE TYPE
        % Get unique dye types from the combined data
        uniqueDyeTypes = unique(dyeTypes);
    
        % Create statistics summary for each dye type
        allDyesTable = [];
        
        
    
            %% STEP 12: CALCULATE STATISTICS FOR EACH METRIC
        % Initialize summary table for all metrics and dyes
        summaryTable = table();
        
        % Create a simplified summary table for final output
        simplifiedSummary = table();
        
        for dyeIdx = 1:length(uniqueDyeTypes)
            dyeType = uniqueDyeTypes{dyeIdx};
            fprintf('Processing statistics for dye type: %s\n', dyeType);
            
            % Get data for this dye type
            dyeMask = strcmp(dyeTypes, dyeType);
            dyeData = combinedTable(dyeMask, :);
            
            % Get unique metrics for this dye
            uniqueMetrics = unique(dyeData.Variable);
            
            % Initialize summary table for this dye
            summaryTable = table();
            summaryTable.Metric = uniqueMetrics;
            summaryTable.Dye = repmat({dyeType}, length(uniqueMetrics), 1);
            summaryTable.Movie_Mean = nan(length(uniqueMetrics), 1);
            summaryTable.Movie_Std = nan(length(uniqueMetrics), 1);
            summaryTable.Movie_SEM = nan(length(uniqueMetrics), 1);
            summaryTable.Num_Movies = zeros(length(uniqueMetrics), 1);
            
            % Process each metric
            for metricIdx = 1:length(uniqueMetrics)
                metric = uniqueMetrics{metricIdx};
                
                % Get data for this metric
                metricMask = strcmp(dyeData.Variable, metric);
                metricData = dyeData(metricMask, :);
                
                % Get unique movies (run names) for this metric
                runNames = metricData.RunName;
                uniqueMovies = unique(runNames);
                numMovies = length(uniqueMovies);
                
                % Calculate mean value for each movie
                movieMeans = zeros(numMovies, 1);
                for movieIdx = 1:numMovies
                    movieMask = strcmp(runNames, uniqueMovies{movieIdx});
                    movieValues = metricData.Mean(movieMask);
                    movieMeans(movieIdx) = mean(movieValues, 'omitnan');
                end
                
                % Calculate statistics across movies
                summaryTable.Movie_Mean(metricIdx) = mean(movieMeans, 'omitnan');
                summaryTable.Movie_Std(metricIdx) = std(movieMeans, 'omitnan');
                summaryTable.Movie_SEM(metricIdx) = summaryTable.Movie_Std(metricIdx) / sqrt(numMovies);
                summaryTable.Num_Movies(metricIdx) = numMovies;
                
                % Get unique movies (run names) for this metric
                runNames = metricData.RunName;
                uniqueMovies = unique(runNames);
                numMovies = length(uniqueMovies);
                
                % Calculate mean value for each movie
                movieMeans = zeros(numMovies, 1);
                for movieIdx = 1:numMovies
                    movieMask = strcmp(runNames, uniqueMovies{movieIdx});
                    movieValues = metricData.Mean(movieMask);
                    movieMeans(movieIdx) = mean(movieValues, 'omitnan');
                end
                
                % Calculate statistics across movies
                summaryTable.Movie_Mean(metricIdx) = mean(movieMeans, 'omitnan');
                summaryTable.Movie_Std(metricIdx) = std(movieMeans, 'omitnan');
                summaryTable.Movie_SEM(metricIdx) = summaryTable.Movie_Std(metricIdx) / sqrt(numMovies);
                summaryTable.Num_Movies(metricIdx) = numMovies;
                
                % Extract SR numbers or sample identifiers for each movie
            movieSamples = cell(numMovies, 1);
            for movieIdx = 1:numMovies
                movieName = uniqueMovies{movieIdx};
                srMatch = regexp(movieName, 'SR(\d+)', 'tokens');
                
                if ~isempty(srMatch)
                    movieSamples{movieIdx} = str2double(srMatch{1}{1});
                else
                    % Try to extract sample identifier from run name
                    [~, fileName] = fileparts(movieName);
                    bcaMatch = regexp(fileName, '(.+?)_BCA_', 'tokens');
                    
                    if ~isempty(bcaMatch)
                        movieSamples{movieIdx} = bcaMatch{1}{1};
                    else
                        sampleIdMatch = regexp(movieName, '([^_\\\/]+)_Run_\d{8}_\d{6}', 'tokens');
                        if ~isempty(sampleIdMatch)
                            movieSamples{movieIdx} = sampleIdMatch{1}{1};
                        else
                            movieSamples{movieIdx} = sprintf('Unknown_%d', movieIdx);
                        end
                    end
                end
            end
            
            % Add sample info to summary table
            summaryTable.Sample_IDs{metricIdx} = unique(movieSamples);
            
            % Create a list of sample identifiers for display
            uniqueSamples = unique(movieSamples);
            sampleListStr = '';
            
            for sampleIdx = 1:length(uniqueSamples)
                if sampleIdx == 1
                    if isnumeric(uniqueSamples{sampleIdx})
                        sampleListStr = sprintf('SR%d', uniqueSamples{sampleIdx});
                    else
                        sampleListStr = uniqueSamples{sampleIdx};
                    end
                else
                    if isnumeric(uniqueSamples{sampleIdx})
                        sampleListStr = sprintf('%s, SR%d', sampleListStr, uniqueSamples{sampleIdx});
                    else
                        sampleListStr = sprintf('%s, %s', sampleListStr, uniqueSamples{sampleIdx});
                    end
                end
            end
            
            summaryTable.SR_List{metricIdx} = sampleListStr;
        end
            

            % Create formatted display table
            displayTable = table();
            displayTable.Metric = summaryTable.Metric;
            displayTable.Dye = summaryTable.Dye;
            displayTable.SR_List = summaryTable.SR_List;  % Add the SR list instead of range
            displayTable.Movie_Statistics_Mean_SEM = cell(height(summaryTable), 1);
            displayTable.Movie_Statistics_Mean_STD = cell(height(summaryTable), 1);

            % Format statistics for display
            for j = 1:height(summaryTable)
                metric = summaryTable.Metric{j};
                
                % Format based on metric type
                if contains(metric, 'Time_ms')
                    % Convert time from ms to s for display
                    mean_value = summaryTable.Movie_Mean(j) / 1000;
                    sem_value = summaryTable.Movie_SEM(j) / 1000;
                    std_value = summaryTable.Movie_Std(j) / 1000;
                    unit = 's';
                elseif contains(metric, 'Photons')
                    % Convert photon counts to thousands
                    mean_value = summaryTable.Movie_Mean(j) / 1000;
                    sem_value = summaryTable.Movie_SEM(j) / 1000;
                    std_value = summaryTable.Movie_Std(j) / 1000;
                    unit = '×10³';
                elseif contains(metric, 'DutyCycle') || startsWith(metric, 'DC_')
                    % For duty cycle metrics, use more decimal places
                    mean_value = summaryTable.Movie_Mean(j) * 100; % Convert to percentage
                    sem_value = summaryTable.Movie_SEM(j) * 100;
                    std_value = summaryTable.Movie_Std(j) * 100;
                    unit = '%';
                    displayTable.Movie_Statistics_Mean_SEM{j} = sprintf('%.2f ± %.2f %s (n=%d movies)', ...
                        mean_value, sem_value, unit, summaryTable.Num_Movies(j));
                    displayTable.Movie_Statistics_Mean_STD{j} = sprintf('%.2f ± %.2f %s (n=%d movies)', ...
                        mean_value, std_value, unit, summaryTable.Num_Movies(j));
                    continue;
                else
                    mean_value = summaryTable.Movie_Mean(j);
                    sem_value = summaryTable.Movie_SEM(j);
                    std_value = summaryTable.Movie_Std(j);
                    unit = '';
                end
                
                displayTable.Movie_Statistics_Mean_SEM{j} = sprintf('%.3f ± %.3f %s (n=%d movies)', ...
                    mean_value, sem_value, unit, summaryTable.Num_Movies(j));
                displayTable.Movie_Statistics_Mean_STD{j} = sprintf('%.3f ± %.3f %s (n=%d movies)', ...
                    mean_value, std_value, unit, summaryTable.Num_Movies(j));
            end
            
            % Add key metrics to simplified summary
            newRow = height(simplifiedSummary) + 1;
            
            % Initialize all fields for the new row to avoid warnings
            if newRow == 1
                % First row - create the structure with all fields
                simplifiedSummary = table();
                % Define all columns that will be used
                simplifiedSummary.Dye = cell(1,1);
                simplifiedSummary.SR_List = cell(1,1);
                simplifiedSummary.MeanOnTime_s_SEM = cell(1,1);
                simplifiedSummary.MeanOnTime_s_STD = cell(1,1);
                simplifiedSummary.MeanOffTime_s_SEM = cell(1,1);
                simplifiedSummary.MeanOffTime_s_STD = cell(1,1);
                simplifiedSummary.DutyCycle_pct_SEM = cell(1,1);
                simplifiedSummary.DutyCycle_pct_STD = cell(1,1);
                simplifiedSummary.BlinkingEvents_SEM = cell(1,1);
                simplifiedSummary.BlinkingEvents_STD = cell(1,1);
                simplifiedSummary.PhotonsPerCycle_SEM = cell(1,1);
                simplifiedSummary.PhotonsPerCycle_STD = cell(1,1);
                simplifiedSummary.PhotonsPerDetection_SEM = cell(1,1);
                simplifiedSummary.PhotonsPerDetection_STD = cell(1,1);
                
                % Initialize the first row with empty values
                simplifiedSummary.Dye{1} = '';
                simplifiedSummary.SR_List{1} = '';
                simplifiedSummary.MeanOnTime_s_SEM{1} = '';
                simplifiedSummary.MeanOnTime_s_STD{1} = '';
                simplifiedSummary.MeanOffTime_s_SEM{1} = '';
                simplifiedSummary.MeanOffTime_s_STD{1} = '';
                simplifiedSummary.DutyCycle_pct_SEM{1} = '';
                simplifiedSummary.DutyCycle_pct_STD{1} = '';
                simplifiedSummary.BlinkingEvents_SEM{1} = '';
                simplifiedSummary.BlinkingEvents_STD{1} = '';
                simplifiedSummary.PhotonsPerCycle_SEM{1} = '';
                simplifiedSummary.PhotonsPerCycle_STD{1} = '';
                simplifiedSummary.PhotonsPerDetection_SEM{1} = '';
                simplifiedSummary.PhotonsPerDetection_STD{1} = '';
            else
                % For subsequent rows, use a more robust approach
                % Create a temporary row with the same structure as the table
                tempRow = simplifiedSummary(1,:);
                
                % Clear all values in the temporary row
                for colIdx = 1:width(simplifiedSummary)
                    tempRow{1, colIdx} = {''};
                end
                
                % Add the empty row to the table
                simplifiedSummary = [simplifiedSummary; tempRow];
            end
            
            % Now assign the actual values
            simplifiedSummary.Dye{newRow} = dyeType;
            simplifiedSummary.SR_List{newRow} = summaryTable.SR_List{1}; % Add SR list to simplified summary
            
            % Find and add key metrics in the specific format
            metricMap = containers.Map();
            for j = 1:height(summaryTable)
                metricMap(summaryTable.Metric{j}) = j;
            end
            
            % Mean On Time (s) - Fix to properly display in seconds
        if isKey(metricMap, 'MeanOnTime_ms')
            idx = metricMap('MeanOnTime_ms');
            % Convert from ms to s and format with proper units
            simplifiedSummary.MeanOnTime_s_SEM{newRow} = sprintf('%.3f ± %.3f', ...
                summaryTable.Movie_Mean(idx)/1000, summaryTable.Movie_SEM(idx)/1000);
            simplifiedSummary.MeanOnTime_s_STD{newRow} = sprintf('%.3f ± %.3f', ...
                summaryTable.Movie_Mean(idx)/1000, summaryTable.Movie_Std(idx)/1000);
        else
            simplifiedSummary.MeanOnTime_s_SEM{newRow} = 'N/A';
            simplifiedSummary.MeanOnTime_s_STD{newRow} = 'N/A';
        end

        % Mean Off Time (s) - Add this metric to the summary
        if isKey(metricMap, 'MeanOffTime_ms')
            idx = metricMap('MeanOffTime_ms');
            % Convert from ms to s and format with proper units
            simplifiedSummary.MeanOffTime_s_SEM{newRow} = sprintf('%.3f ± %.3f', ...
                summaryTable.Movie_Mean(idx)/1000, summaryTable.Movie_SEM(idx)/1000);
            simplifiedSummary.MeanOffTime_s_STD{newRow} = sprintf('%.3f ± %.3f', ...
                summaryTable.Movie_Mean(idx)/1000, summaryTable.Movie_Std(idx)/1000);
        else
            simplifiedSummary.MeanOffTime_s_SEM{newRow} = 'N/A';
            simplifiedSummary.MeanOffTime_s_STD{newRow} = 'N/A';
        end
            
            % Duty Cycle (%)
            if isKey(metricMap, 'DutyCycle')
                idx = metricMap('DutyCycle');
                simplifiedSummary.DutyCycle_pct_SEM{newRow} = sprintf('%.2f ± %.2f', ...
                    summaryTable.Movie_Mean(idx)*100, summaryTable.Movie_SEM(idx)*100);
                simplifiedSummary.DutyCycle_pct_STD{newRow} = sprintf('%.2f ± %.2f', ...
                    summaryTable.Movie_Mean(idx)*100, summaryTable.Movie_Std(idx)*100);
            else
                simplifiedSummary.DutyCycle_pct_SEM{newRow} = 'N/A';
                simplifiedSummary.DutyCycle_pct_STD{newRow} = 'N/A';
            end
            
            % Blinking Events
            if isKey(metricMap, 'BlinkingEvents')
                idx = metricMap('BlinkingEvents');
                simplifiedSummary.BlinkingEvents_SEM{newRow} = sprintf('%.2f ± %.2f', ...
                    summaryTable.Movie_Mean(idx), summaryTable.Movie_SEM(idx));
                simplifiedSummary.BlinkingEvents_STD{newRow} = sprintf('%.2f ± %.2f', ...
                    summaryTable.Movie_Mean(idx), summaryTable.Movie_Std(idx));
            else
                simplifiedSummary.BlinkingEvents_SEM{newRow} = 'N/A';
                simplifiedSummary.BlinkingEvents_STD{newRow} = 'N/A';
            end
            
            % Photons Per Cycle (×10³)
            if isKey(metricMap, 'PhotonsPerCycle')
                idx = metricMap('PhotonsPerCycle');
                simplifiedSummary.PhotonsPerCycle_SEM{newRow} = sprintf('%.2f ± %.2f', ...
                    summaryTable.Movie_Mean(idx)/1000, summaryTable.Movie_SEM(idx)/1000);
                simplifiedSummary.PhotonsPerCycle_STD{newRow} = sprintf('%.2f ± %.2f', ...
                    summaryTable.Movie_Mean(idx)/1000, summaryTable.Movie_Std(idx)/1000);
            else
                simplifiedSummary.PhotonsPerCycle_SEM{newRow} = 'N/A';
                simplifiedSummary.PhotonsPerCycle_STD{newRow} = 'N/A';
            end
            
            % Photons Per Detection (×10³)
            if isKey(metricMap, 'PhotonsPerDetection')
                idx = metricMap('PhotonsPerDetection');
                simplifiedSummary.PhotonsPerDetection_SEM{newRow} = sprintf('%.2f ± %.2f', ...
                    summaryTable.Movie_Mean(idx)/1000, summaryTable.Movie_SEM(idx)/1000);
                simplifiedSummary.PhotonsPerDetection_STD{newRow} = sprintf('%.2f ± %.2f', ...
                    summaryTable.Movie_Mean(idx)/1000, summaryTable.Movie_Std(idx)/1000);
            else
                simplifiedSummary.PhotonsPerDetection_SEM{newRow} = 'N/A';
                simplifiedSummary.PhotonsPerDetection_STD{newRow} = 'N/A';
            end
            
            % Display the formatted table
            disp(displayTable);
            
            % Get current timestamp for file naming
            timestamp = datestr(now, 'yyyymmdd_HHMMSS');
            
            % Extract analysis tag for file naming
            analysisTag = '';
            if ~isempty(selectedIndices)
                firstRunIdx = selectedIndices(1);
                [~, firstDirName] = fileparts(runDirs{firstRunIdx});
                
                if contains(firstDirName, '_MSK')
                    maskParams = regexp(firstDirName, '_MSK(\d+)_THR(\d+)', 'tokens');
                    if ~isempty(maskParams)
                        analysisTag = sprintf('_MSK%s_THR%s', maskParams{1}{1}, maskParams{1}{2});
                    end
                elseif contains(firstDirName, '_HFSQR')
                    imgParams = regexp(firstDirName, '_HFSQR(\d+)_THR(\d+)', 'tokens');
                    if ~isempty(imgParams)
                        analysisTag = sprintf('_HFSQR%s_THR%s', imgParams{1}{1}, imgParams{1}{2});
                    end
                end
            end
            
            % Save detailed table with separate mean/std columns
            writetable(summaryTable, fullfile(outputDir, ...
                sprintf('%s%s_combined_statistics_detailed_%s.csv', dyeType, analysisTag, timestamp)));
            
            % Save formatted table with mean ± std in each cell
            writetable(displayTable, fullfile(outputDir, ...
                sprintf('%s%s_combined_statistics_summary_%s.csv', dyeType, analysisTag, timestamp)));
            
            % Add to the all dyes table
            if isempty(allDyesTable)
                allDyesTable = displayTable;
            else
                allDyesTable = [allDyesTable; displayTable];
            end
        end
        
        % Create a formatted simplified summary table for display
        if height(simplifiedSummary) > 0
            % Create a formatted table with proper headers
            formattedTable = table();
            formattedTable.Dye = simplifiedSummary.Dye;
            formattedTable.SR_List = simplifiedSummary.SR_List;
            
            % Use the new field names with _SEM suffix
            formattedTable.MeanOnTime_s = simplifiedSummary.MeanOnTime_s_SEM;
            formattedTable.MeanOffTime_s = simplifiedSummary.MeanOffTime_s_SEM;  % Add Mean Off Time
            formattedTable.DutyCycle_pct = simplifiedSummary.DutyCycle_pct_SEM;
            formattedTable.BlinkingEvents = simplifiedSummary.BlinkingEvents_SEM;
            formattedTable.PhotonsPerCycle = simplifiedSummary.PhotonsPerCycle_SEM;
            formattedTable.PhotonsPerDetection = simplifiedSummary.PhotonsPerDetection_SEM;
            
            % Also create a version with STD
            formattedTableSTD = table();
            formattedTableSTD.Dye = simplifiedSummary.Dye;
            formattedTableSTD.SR_List = simplifiedSummary.SR_List;
            formattedTableSTD.MeanOnTime_s = simplifiedSummary.MeanOnTime_s_STD;
            formattedTableSTD.MeanOffTime_s = simplifiedSummary.MeanOffTime_s_STD;  % Add Mean Off Time
            formattedTableSTD.DutyCycle_pct = simplifiedSummary.DutyCycle_pct_STD;
            formattedTableSTD.BlinkingEvents = simplifiedSummary.BlinkingEvents_STD;
            formattedTableSTD.PhotonsPerCycle = simplifiedSummary.PhotonsPerCycle_STD;
            formattedTableSTD.PhotonsPerDetection = simplifiedSummary.PhotonsPerDetection_STD;
            
            % Rename columns for display
            formattedTable.Properties.VariableNames = {'Dye', 'SR Numbers', 'Mean On Time (s)', ...
                'Mean Off Time (s)', 'Duty Cycle (%)', 'Blinking Events', ...
                'Photons/Cycle (×10³)', 'Photons/Detection (×10³)'};
            
            formattedTableSTD.Properties.VariableNames = {'Dye', 'SR Numbers', 'Mean On Time (s)', ...
                'Mean Off Time (s)', 'Duty Cycle (%)', 'Blinking Events', ...
                'Photons/Cycle (×10³)', 'Photons/Detection (×10³)'};
            
            % Display the simplified table
            disp('Summary of Key Metrics (Mean ± SEM):');
            disp(formattedTable);
            
            disp('Summary of Key Metrics (Mean ± STD):');
            disp(formattedTableSTD);
            
            % Save the simplified summaries
            writetable(formattedTable, fullfile(outputDir, ...
                sprintf('All_Dyes%s_key_metrics_summary_SEM_%s.csv', analysisTag, timestamp)));
            
            writetable(formattedTableSTD, fullfile(outputDir, ...
                sprintf('All_Dyes%s_key_metrics_summary_STD_%s.csv', analysisTag, timestamp)));
                
            fprintf('\nSaved key metrics summaries at %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        end
        
        % Save the combined table with all dye types
        if height(allDyesTable) > 0
            writetable(allDyesTable, fullfile(outputDir, ...
                sprintf('All_Dyes%s_combined_statistics_%s.csv', analysisTag, timestamp)));
            fprintf('\nSaved combined statistics for all dye types at %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        end
        
        fprintf('\nCombined statistics saved to: %s\n', outputDir);
    end
