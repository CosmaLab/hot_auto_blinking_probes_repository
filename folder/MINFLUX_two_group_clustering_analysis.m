%% MINFLUX Density-Based K-means Clustering Analysis Script 
% Multi-Cell Version - Processes two groups: Euchromatin and Heterochromatin
% Each group contains multiple cell datasets
% Parameters are customizable for optimized analysis of different datasets
% Clustering based on original points 
% Cluster area calculation - Convex Hull for clusters with >=3 points
% Two-group comparison with detailed individual cell analysis and combined analysis results

clear, clc, close all

%% 0. CUSTOM PARAMETER CONFIGURATION FOR BOTH GROUPS
% Users can customize all clustering parameters here
% Modify the values directly

% ========== CUSTOM CLUSTERING PARAMETERS ==========
density_pixel_size = 4;           % Density map pixel size (nm) - controls density map resolution
gaussian_sigma = 1;                % Gaussian smoothing coefficient (pixels) - controls smoothing level
density_threshold_factor = 0.10;   % Density threshold factor - controls peak detection sensitivity
min_peak_distance = 30;            % Minimum peak distance (nm) - controls minimum spacing between clusters
% ==================================================

% ===== VISUALIZATION: XY projection =====
xy_point_size = 30;      % trace point size
xy_center_size = 30;    % cluster center 'x' size
xy_center_linewidth = 1.5;
% =====================================================

% ========== TWO-GROUP CONFIGURATION ==========
% Specify folders for Euchromatin and Heterochromatin groups
% Each folder should contain .mat files for each cell
euchromatin_folder = 'C:\Users\18317\Desktop\code 20260225\Euchromatin_5cells';
heterochromatin_folder = 'C:\Users\18317\Desktop\code 20260225\Heterochromatin_3cells';

% Group names for display
group_names = {'Euchromatin', 'Heterochromatin'};
group_folders = {euchromatin_folder, heterochromatin_folder};
nGroups = length(group_names);

% Store custom parameters
CUSTOM_PARAMS.density_pixel_size = density_pixel_size;
CUSTOM_PARAMS.gaussian_sigma = gaussian_sigma;
CUSTOM_PARAMS.density_threshold_factor = density_threshold_factor;
CUSTOM_PARAMS.min_peak_distance = min_peak_distance;
CUSTOM_PARAMS.xy_point_size = xy_point_size;
CUSTOM_PARAMS.xy_center_size = xy_center_size;
CUSTOM_PARAMS.xy_center_linewidth = xy_center_linewidth;

% Display parameters
fprintf('==========================================\n');
fprintf('TWO-GROUP ANALYSIS PARAMETERS:\n');
fprintf('==========================================\n');
fprintf('  Density pixel size: %.1f nm\n', density_pixel_size);
fprintf('  Gaussian sigma: %.1f pixels\n', gaussian_sigma);
fprintf('  Density threshold factor: %.2f\n', density_threshold_factor);
fprintf('  Minimum peak distance: %.1f nm\n', min_peak_distance);
fprintf('  Euchromatin group folder: %s\n', euchromatin_folder);
fprintf('  Heterochromatin group folder: %s\n', heterochromatin_folder);
fprintf('==========================================\n\n');

%% Initialize storage for all groups
all_groups_results = cell(nGroups, 1);
group_summary = struct();
group_cell_stats = cell(nGroups, 1); % Store statistics for each cell in each group

% Process each group
for group_idx = 1:nGroups
    current_group = group_names{group_idx};
    current_folder = group_folders{group_idx};
    
    fprintf('\n==========================================\n');
    fprintf('PROCESSING GROUP %d/%d: %s\n', group_idx, nGroups, current_group);
    fprintf('Folder: %s\n', current_folder);
    fprintf('==========================================\n\n');
    
    %% 1. Find all .mat files in the group folder
    fprintf('Searching for .mat files in folder: %s\n', current_folder);
    
    % Find all .mat files
    matFiles = dir(fullfile(current_folder, '*.mat'));
    
    if isempty(matFiles)
        error('No .mat files found in the specified folder: %s', current_folder);
    end
    
    % Get cell names and file names
    cell_names = cell(length(matFiles), 1);
    file_names = cell(length(matFiles), 1);
    
    for i = 1:length(matFiles)
        [~, name, ~] = fileparts(matFiles(i).name);
        cell_names{i} = name;
        file_names{i} = matFiles(i).name;
    end
    
    nCells = length(cell_names);
    fprintf('Found %d .mat files in the folder:\n', nCells);
    for i = 1:nCells
        fprintf('  %d. %s\n', i, cell_names{i});
    end
    fprintf('\n');
    
    %% Initialize storage for current group
    group_results = [];
    group_all_areas = [];
    group_all_nnd = [];
    group_all_traceIDs = [];
    group_all_densities = [];
    
    % Initialize storage for cell-level statistics
    cell_stats = struct();
    cell_stats.cell_names = cell_names;
    cell_stats.areas_data = cell(nCells, 1); % Store area data for each cell
    cell_stats.nnd_data = cell(nCells, 1);   % Store NND data for each cell
    cell_stats.traceIDs_data = cell(nCells, 1); % Store trace IDs per cluster for each cell
    cell_stats.densities_data = cell(nCells, 1); % Store density data for each cell
    
    % Loop through each cell in the group
    for cell_idx = 1:nCells
        current_cell = cell_names{cell_idx};
        current_file = file_names{cell_idx};
        
        fprintf('\n==========================================\n');
        fprintf('PROCESSING CELL %d/%d: %s\n', cell_idx, nCells, current_cell);
        fprintf('File: %s\n', current_file);
        fprintf('==========================================\n\n');
        
        %% 2. Configuration and Parameters for current cell
        % Create export folder for this cell
        exportFolder = fullfile(current_folder, ['export_' current_cell]);
        if ~exist(exportFolder, 'dir')
            mkdir(exportFolder);
        end
        
        % Load MINFLUX .mat file
        filePath = fullfile(current_folder, current_file);
        fprintf('Loading MINFLUX data: %s\n', current_file);
        
        % Load data
        data = load(filePath);
        
        %% 3. Extract and Process MINFLUX Localizations - Apply Filtering Conditions
        fprintf('Extracting and filtering MINFLUX localizations...\n');
        
        % Step 1: Get initial valid localizations (data.vld == 1)
        initial_valid = data.vld == 1;
        initial_count = sum(initial_valid);
        fprintf('Step 1 - Valid localizations (vld==1): %d\n', initial_count);
        
        % Extract all localization data initially
        all_locations = squeeze(data.itr.loc(:, end, :)); % Last iteration positions
        all_traceIDs = data.tid;
        
        % Apply initial valid filter
        locations_step1 = all_locations(initial_valid, :);
        traceIDs_step1 = all_traceIDs(initial_valid);
        
        % Step 2: Filter by cfr (itr.cfr last column <= 0.8)
        fprintf('Step 2 - Applying cfr filter (<= 0.8)...\n');
        if isfield(data.itr, 'cfr')
            % Check dimensions
            cfr_data = data.itr.cfr(:, end);
            cfr_valid = cfr_data <= 0.8;
            
            % Apply cfr filter to the already valid indices
            all_indices = 1:length(data.vld);
            valid_indices = all_indices(initial_valid); % Indices that passed step 1
            
            % Now filter these by cfr condition
            cfr_valid_for_step1 = cfr_valid(valid_indices);
            
            % Apply filter
            locations_step2 = locations_step1(cfr_valid_for_step1, :);
            traceIDs_step2 = traceIDs_step1(cfr_valid_for_step1);
            
            step2_count = sum(cfr_valid_for_step1);
            fprintf('  After cfr filter: %d localizations (removed %d)\n', ...
                step2_count, initial_count - step2_count);
        else
            locations_step2 = locations_step1;
            traceIDs_step2 = traceIDs_step1;
            step2_count = initial_count;
            fprintf('  cfr field not found, skipping this filter\n');
        end
        
        % Step 3: Filter by efo (itr.efo last column <= 70000)
        fprintf('Step 3 - Applying efo filter (<= 70000)...\n');
        if isfield(data.itr, 'efo')
            % Check dimensions
            efo_data = data.itr.efo(:, end);
            efo_valid = efo_data <= 70000;
            
            % We need to track indices through all filtering steps
            all_indices = 1:length(data.vld);
            initial_valid_indices = all_indices(initial_valid);
            
            if exist('cfr_valid_for_step1', 'var')
                % Get indices that passed both step1 and step2
                step2_indices = initial_valid_indices(cfr_valid_for_step1);
            else
                step2_indices = initial_valid_indices;
            end
            
            % Apply efo filter to current indices
            efo_valid_for_current = efo_valid(step2_indices);
            
            % Apply filter
            locations_step3 = locations_step2(efo_valid_for_current, :);
            traceIDs_step3 = traceIDs_step2(efo_valid_for_current);
            
            step3_count = sum(efo_valid_for_current);
            fprintf('  After efo filter: %d localizations (removed %d)\n', ...
                step3_count, step2_count - step3_count);
        else
            locations_step3 = locations_step2;
            traceIDs_step3 = traceIDs_step2;
            step3_count = step2_count;
            fprintf('  efo field not found, skipping this filter\n');
        end
        
        % Step 4: Filter traces with at least 2 points
        fprintf('Step 4 - Filtering traces with at least 2 points...\n');
        
        % Get unique trace IDs and count points in each trace
        unique_traceIDs_step3 = unique(traceIDs_step3);
        trace_counts = zeros(length(unique_traceIDs_step3), 1);
        
        for i = 1:length(unique_traceIDs_step3)
            trace_counts(i) = sum(traceIDs_step3 == unique_traceIDs_step3(i));
        end
        
        % Identify traces with at least 2 points
        valid_traces = unique_traceIDs_step3(trace_counts >= 2);
        trace_filter = ismember(traceIDs_step3, valid_traces);
        
        % Apply trace filter
        locations_final = locations_step3(trace_filter, :);
        traceIDs_final = traceIDs_step3(trace_filter);
        
        final_count = sum(trace_filter);
        trace_count = length(valid_traces);
        fprintf('  After trace filter: %d traces with >=2 points, %d localizations\n', ...
            trace_count, final_count);
        fprintf('  Removed %d traces with <2 points\n', sum(trace_counts < 2));
        
        % Summary of filtering
        fprintf('\n==========================================\n');
        fprintf('FILTERING SUMMARY for %s:\n', current_cell);
        fprintf('==========================================\n');
        fprintf('Initial valid localizations (vld==1): %d\n', initial_count);
        fprintf('After cfr filter: %d\n', step2_count);
        fprintf('After efo filter: %d\n', step3_count);
        fprintf('After trace filter: %d traces, %d localizations\n', trace_count, final_count);
        fprintf('==========================================\n\n');
        
        %% 4. Extract All Localization Points 
        fprintf('Extracting all localization points for clustering...\n');
        
        % Convert to nanometers
        locations_final_nm = locations_final * 1e9;
        
        % Use X and Y coordinates for 2D clustering
        XY_all = locations_final_nm(:, 1:2);
        
        % Keep traceIDs for each point
        traceIDs_all = traceIDs_final;
        
        fprintf('Total localizations for clustering: %d\n', size(XY_all, 1));
        fprintf('Total unique trace IDs: %d\n', length(unique(traceIDs_all)));
        
        %% 5. Calculate Trace Centroids for Density Map
        fprintf('Calculating trace centroids for density map...\n');
        
        % Get unique trace IDs
        uniqueTraceIDs = unique(traceIDs_all);
        nTraces = length(uniqueTraceIDs);
        
        % Initialize arrays for trace centroids
        trace_centroids = zeros(nTraces, 2);
        
        % Calculate centroid (mean) for each trace
        for i = 1:nTraces
            traceID = uniqueTraceIDs(i);
            
            % Get all points of this trace
            trace_indices = find(traceIDs_all == traceID);
            trace_points = XY_all(trace_indices, :);
            
            % Calculate centroid (mean) of trace points
            trace_centroids(i, :) = mean(trace_points, 1);
        end
        
        fprintf('Calculated centroids for %d traces\n', nTraces);
        
        %% 6. Save All Localization Points Data to CSV File
        fprintf('Saving all localization points to CSV...\n');
        
        % Create localization points data table
        locPointsData = table(traceIDs_all(:), ...
                               locations_final_nm(:,1), ...
                               locations_final_nm(:,2), ...
                               locations_final_nm(:,3), ...
                               ones(size(locations_final_nm,1),1)); % Each point is a localization
        
        locPointsData.Properties.VariableNames = {'TraceID', 'X_nm', 'Y_nm', 'Z_nm', 'LocalizationIndex'};
        csvFile = fullfile(exportFolder, [current_cell '_AllLocalizations.csv']);
        writetable(locPointsData, csvFile);
        fprintf('All localization points saved to: %s\n', csvFile);
        
        %% 7. Density-Peak Based K-means Clustering on All Points
        fprintf('==========================================\n');
        fprintf('Performing Density-Peak Based K-means Clustering on ALL POINTS\n');
        fprintf('==========================================\n');
        
        X = XY_all; % ALL points for clustering
        nPoints = size(X, 1);
        fprintf('Total points to cluster: %d\n', nPoints);
        
        fprintf('Custom Parameters:\n');
        fprintf('  Density pixel size: %.1f nm\n', density_pixel_size);
        fprintf('  Gaussian sigma for smoothing: %.1f pixels\n', gaussian_sigma);
        fprintf('  Density threshold factor: %.2f\n', density_threshold_factor);
        fprintf('  Minimum peak distance: %.1f nm\n', min_peak_distance);
        
        % Step 1: Create 2D density histogram based on TRACE CENTROIDS
        fprintf('1. Creating 2D density histogram based on TRACE CENTROIDS...\n');
        x_range = range(trace_centroids(:,1));
        y_range = range(trace_centroids(:,2));
        
        % Calculate histogram edges based on trace centroids
        x_edges = floor(min(trace_centroids(:,1))/density_pixel_size)*density_pixel_size : ...
                  density_pixel_size : ...
                  ceil(max(trace_centroids(:,1))/density_pixel_size)*density_pixel_size;
        y_edges = floor(min(trace_centroids(:,2))/density_pixel_size)*density_pixel_size : ...
                  density_pixel_size : ...
                  ceil(max(trace_centroids(:,2))/density_pixel_size)*density_pixel_size;
        
        % Create histogram using trace centroids (one count per trace)
        density_hist = histcounts2(trace_centroids(:,1), trace_centroids(:,2), x_edges, y_edges);
        density_hist = density_hist'; % Transpose to match coordinate system
        
        fprintf('  Density histogram created: %d x %d pixels\n', size(density_hist, 2), size(density_hist, 1));
        fprintf('  Maximum trace density: %d traces per pixel\n', max(density_hist(:)));
        fprintf('  Pixel area: %.1f nm²\n', density_pixel_size^2);
        
        % Convert to density in traces/nm²
        pixel_area_nm2 = density_pixel_size^2;
        density_hist_traces_nm2 = density_hist / pixel_area_nm2;
        
        % Step 2: Apply Gaussian smoothing to density map to reduce noise
        fprintf('2. Applying Gaussian smoothing...\n');
        if gaussian_sigma > 0
            % Create Gaussian filter
            filter_size = ceil(3*gaussian_sigma)*2 + 1;
            gaussian_filter = fspecial('gaussian', filter_size, gaussian_sigma);
            density_smoothed = imfilter(density_hist_traces_nm2, gaussian_filter, 'replicate');
        else
            density_smoothed = density_hist_traces_nm2;
        end
        
        % Step 3: Find local density peaks
        fprintf('3. Finding local density peaks...\n');
        
        % Method 1: Use imregionalmax to find local maxima
        BW = imregionalmax(density_smoothed);
        
        % Get peak positions (pixel coordinates)
        [peak_y, peak_x] = find(BW);
        n_peaks_raw = length(peak_x);
        fprintf('  Found %d raw local maxima\n', n_peaks_raw);
        
        % Convert pixel coordinates to actual coordinates (nm)
        peak_coords_raw = zeros(n_peaks_raw, 2);
        peak_values_raw = zeros(n_peaks_raw, 1);
        
        for i = 1:n_peaks_raw
            peak_coords_raw(i, 1) = x_edges(peak_x(i)) + density_pixel_size/2;
            peak_coords_raw(i, 2) = y_edges(peak_y(i)) + density_pixel_size/2;
            peak_values_raw(i) = density_smoothed(peak_y(i), peak_x(i));
        end
        
        % Step 4: Filter peaks (based on density threshold and minimum distance)
        fprintf('4. Filtering peaks...\n');
        
        % Set density threshold
        max_density = max(peak_values_raw);
        density_threshold = density_threshold_factor * max_density;
        
        % Preliminary screening based on density threshold
        valid_density = peak_values_raw >= density_threshold;
        peak_coords_filtered = peak_coords_raw(valid_density, :);
        peak_values_filtered = peak_values_raw(valid_density);
        
        fprintf('  After density threshold (%.4f traces/nm²): %d peaks\n', ...
            density_threshold, sum(valid_density));
        
        % Further screening based on minimum distance (ensure peaks are not too close)
        if min_peak_distance > 0 && size(peak_coords_filtered, 1) > 1
            % Calculate distance matrix between peaks
            dist_matrix = pdist2(peak_coords_filtered, peak_coords_filtered);
            
            % Identify peaks to keep (based on density value and distance)
            n_filtered = size(peak_coords_filtered, 1);
            keep_peak = true(n_filtered, 1);
            
            for i = 1:n_filtered
                if keep_peak(i)
                    % Find other peaks near current peak
                    too_close = dist_matrix(i, :) < min_peak_distance;
                    too_close(i) = false; % Exclude self
                    
                    % Among peaks that are too close, keep only the highest density one
                    close_indices = find(too_close);
                    for j = 1:length(close_indices)
                        idx = close_indices(j);
                        if keep_peak(idx) && peak_values_filtered(idx) < peak_values_filtered(i)
                            keep_peak(idx) = false;
                        elseif keep_peak(idx) && peak_values_filtered(idx) >= peak_values_filtered(i)
                            keep_peak(i) = false;
                            break;
                        end
                    end
                end
            end
            
            peak_coords_final = peak_coords_filtered(keep_peak, :);
            peak_values_final = peak_values_filtered(keep_peak);
        else
            peak_coords_final = peak_coords_filtered;
            peak_values_final = peak_values_filtered;
        end
        
        n_peaks_final = size(peak_coords_final, 1);
        fprintf('  After distance filtering: %d final peaks\n', n_peaks_final);
        
        % Step 5: Perform K-means clustering using density peaks as initial centers
        fprintf('5. Performing K-means with %d initial centers...\n', n_peaks_final);
        
        if n_peaks_final < 2
            error('Not enough density peaks found. Please adjust parameters (especially density_threshold_factor).');
        end
        
        % K-means parameter settings
        kmeans_options = statset('MaxIter', 1000, 'Display', 'iter');
        
        % Run K-means using density peaks as initial centers
        try
            [labels, centroids, sumd, distances] = kmeans(X, n_peaks_final, ...
                'Start', peak_coords_final, ...
                'Replicates', 1, ...
                'Options', kmeans_options, ...
                'Distance', 'sqeuclidean');
            
            fprintf('K-means converged successfully.\n');
            
        catch ME
            warning('K-means failed with initial centers. Trying random initialization...');
            % If using density peaks fails, try random initialization
            [labels, centroids] = kmeans(X, n_peaks_final, ...
                'Replicates', 3, ...
                'Options', kmeans_options);
        end
        
        % Step 6: Reassign points of the same trace ID to the same cluster
        fprintf('6. Reassigning points of the same trace ID to the same cluster...\n');
        
        % Get unique trace IDs
        uniqueTraceIDs = unique(traceIDs_all);
        nTraces = length(uniqueTraceIDs);
        
        % Initialize reassigned labels
        reassigned_labels = zeros(size(labels));
        
        % For each trace ID
        for i = 1:nTraces
            traceID = uniqueTraceIDs(i);
            
            % Find all points of this trace
            trace_indices = find(traceIDs_all == traceID);
            
            if length(trace_indices) > 0
                % Get current cluster assignments for this trace
                trace_labels = labels(trace_indices);
                
                % Find the most common cluster (mode)
                if length(trace_labels) == 1
                    % Single point, keep its cluster
                    reassigned_labels(trace_indices) = trace_labels;
                else
                    % Multiple points, find the most common cluster
                    [unique_labels, ~, ic] = unique(trace_labels);
                    counts = accumarray(ic, 1);
                    [max_count, max_idx] = max(counts);
                    
                    % Check if there is a tie for most common
                    ties = find(counts == max_count);
                    if length(ties) > 1
                        % Tie: assign to the cluster with the nearest centroid to trace center
                        trace_points = X(trace_indices, :);
                        trace_center = mean(trace_points, 1);
                        
                        % Calculate distances to each tied cluster centroid
                        tie_distances = zeros(length(ties), 1);
                        for t = 1:length(ties)
                            tie_distances(t) = norm(trace_center - centroids(unique_labels(ties(t)), :));
                        end
                        
                        % Choose the closest cluster
                        [~, closest_idx] = min(tie_distances);
                        chosen_cluster = unique_labels(ties(closest_idx));
                    else
                        % Clear winner
                        chosen_cluster = unique_labels(max_idx);
                    end
                    
                    % Assign all points of this trace to the chosen cluster
                    reassigned_labels(trace_indices) = chosen_cluster;
                end
            end
        end
        
        % Update labels with reassigned labels
        labels = reassigned_labels;
        
        % Statistical analysis of clustering results
        uniqueClusters = unique(labels);
        nClusters = length(uniqueClusters);
        
        % Check for empty clusters
        cluster_sizes = zeros(nClusters, 1);
        for k = 1:nClusters
            cluster_sizes(k) = sum(labels == k);
        end
        empty_clusters = sum(cluster_sizes == 0);
        
        if empty_clusters > 0
            fprintf('Warning: %d clusters are empty\n', empty_clusters);
            % Renumber non-empty clusters
            non_empty = cluster_sizes > 0;
            new_labels = zeros(size(labels));
            new_cluster_id = 0;
            
            for k = 1:nClusters
                if cluster_sizes(k) > 0
                    new_cluster_id = new_cluster_id + 1;
                    new_labels(labels == k) = new_cluster_id;
                end
            end
            
            labels = new_labels;
            centroids = centroids(non_empty, :);
            nClusters = new_cluster_id;
            fprintf('Renumbered clusters: now have %d non-empty clusters\n', nClusters);
        end
        
        fprintf('==========================================\n');
        fprintf('Density-Peak K-means Clustering Complete for %s:\n', current_cell);
        fprintf('  Found %d density peaks\n', n_peaks_final);
        fprintf('  Result: %d clusters from %d localization points\n', nClusters, nPoints);
        fprintf('  After reassignment: %d unique trace IDs distributed in clusters\n', nTraces);
        fprintf('==========================================\n');
        
        % Save density map and peak information for visualization
        density_data.x_edges = x_edges;
        density_data.y_edges = y_edges;
        density_data.density_map = density_smoothed; % This is now in traces/nm²
        density_data.density_map_raw = density_hist_traces_nm2; % Raw density map
        density_data.peak_coords = peak_coords_final;
        density_data.peak_values = peak_values_final;
        density_data.trace_centroids = trace_centroids; % Save trace centroids
        
        %% 8. Calculate Cluster Statistics (Based on Original Points, Density Based on Trace IDs)
        if nClusters > 0
            % Extract cluster information
            clusterStats = struct();
            clusterStats.centroids = [];
            clusterStats.nTraces = []; % Number of unique trace IDs in each cluster
            clusterStats.nPoints = []; % Number of localization points in each cluster
            clusterStats.sigmaX = [];
            clusterStats.sigmaY = [];
            clusterStats.areas = [];
            clusterStats.densities = [];
            clusterStats.area_method = cell(nClusters, 1); % Store area calculation method
            clusterStats.labels = labels; % Store cluster labels
            
            % Calculate statistics for each cluster
            for k = 1:nClusters
                % Get all points in this cluster
                clusterPoints = X(labels == k, :);
                clusterTraceIDs = traceIDs_all(labels == k);
                
                nPointsInCluster = size(clusterPoints, 1);
                nTracesInCluster = length(unique(clusterTraceIDs));
                
                % Calculate centroid (mean of all points)
                centroid = mean(clusterPoints, 1);
                
                % Calculate standard deviation (as a measure of cluster spread)
                if nPointsInCluster > 1
                    sigmaX = std(clusterPoints(:, 1));
                    sigmaY = std(clusterPoints(:, 2));
                else
                    sigmaX = 5; % Default value for single point (nm)
                    sigmaY = 5; % Default value for single point (nm)
                end
                
                % Calculate cluster area based on number of points
                if nPointsInCluster >= 3
                    % For clusters with 3 or more points: use Convex Hull area
                    try
                        % Calculate convex hull
                        [K, area] = convhull(clusterPoints(:,1), clusterPoints(:,2));
                        area_method = 'Convex Hull';
                        
                        % Optional: Display convex hull for debugging
                        if k <= 5 && nPointsInCluster < 100  % Only show first 5 small clusters
                            fprintf('  Cluster %d: %d points, Convex Hull area = %.1f nm²\n', ...
                                k, nPointsInCluster, area);
                        end
                    catch ME
                        % If convex hull fails (e.g., collinear points), use elliptical area
                        area = pi * sigmaX * sigmaY;
                        area_method = 'Elliptical (Convex Hull failed)';
                        fprintf('  Warning: Convex Hull failed for cluster %d (%d points). Using elliptical area.\n', ...
                            k, nPointsInCluster);
                    end
                else
                    % For 1 or 2 points: use elliptical area (original method)
                    area = pi * sigmaX * sigmaY;
                    if nPointsInCluster == 1
                        area_method = 'Elliptical (single point)';
                    else
                        area_method = 'Elliptical (2 points)';
                    end
                end
                
                % Calculate density: nTraces / area (traces per nm^2)
                if area > 0
                    density = nTracesInCluster / area;
                else
                    density = NaN; % Handle zero area case
                end
                
                % Store results
                clusterStats.centroids = [clusterStats.centroids; centroid];
                clusterStats.nTraces = [clusterStats.nTraces; nTracesInCluster];
                clusterStats.nPoints = [clusterStats.nPoints; nPointsInCluster];
                clusterStats.sigmaX = [clusterStats.sigmaX; sigmaX];
                clusterStats.sigmaY = [clusterStats.sigmaY; sigmaY];
                clusterStats.areas = [clusterStats.areas; area];
                clusterStats.densities = [clusterStats.densities; density];
                clusterStats.area_method{k} = area_method;
            end
            
            % Calculate Nearest Neighbor Distance (NND) for clusters
            if nClusters > 1
                nndValues = zeros(nClusters, 1);
                
                for i = 1:nClusters
                    % Calculate distances to all other cluster centers
                    distances = zeros(nClusters, 1);
                    for j = 1:nClusters
                        if i == j
                            distances(j) = Inf; % Exclude self
                        else
                            % Calculate Euclidean distance
                            distances(j) = sqrt(sum((clusterStats.centroids(j, :) - clusterStats.centroids(i, :)).^2));
                        end
                    end
                    
                    % Find minimum distance
                    nndValues(i) = min(distances);
                end
                
                clusterStats.nnd = nndValues;
            else
                clusterStats.nnd = [];
            end
        else
            clusterStats = [];
            fprintf('No clusters found!\n');
            % Continue to next cell if no clusters found
            continue;
        end
        
        %% 9. Calculate Overall Cluster Statistics for current cell
        % Overall statistics
        medianClusterArea = median(clusterStats.areas);
        meanClusterArea = mean(clusterStats.areas);
        stdClusterArea = std(clusterStats.areas);
        
        % Density statistics
        validDensities = clusterStats.densities(~isnan(clusterStats.densities));
        if ~isempty(validDensities)
            medianDensity = median(validDensities);
            meanDensity = mean(validDensities);
            stdDensity = std(validDensities);
        else
            medianDensity = NaN;
            meanDensity = NaN;
            stdDensity = NaN;
        end
        
        if nClusters > 1
            medianNND = median(clusterStats.nnd);
            meanNND = mean(clusterStats.nnd);
            stdNND = std(clusterStats.nnd);
        else
            medianNND = NaN;
            meanNND = NaN;
            stdNND = NaN;
        end
        
        % Count clusters by area calculation method
        convex_hull_clusters = sum(strcmp(clusterStats.area_method, 'Convex Hull'));
        
        % Use strfind for compatibility with older MATLAB versions
        elliptical_clusters = 0;
        convex_hull_failed = 0;
        elliptical_1_2_points = 0;
        
        for i = 1:length(clusterStats.area_method)
            if strcmp(clusterStats.area_method{i}, 'Elliptical (single point)') || ...
               strcmp(clusterStats.area_method{i}, 'Elliptical (2 points)')
                elliptical_clusters = elliptical_clusters + 1;
                elliptical_1_2_points = elliptical_1_2_points + 1;
            elseif strcmp(clusterStats.area_method{i}, 'Elliptical (Convex Hull failed)')
                elliptical_clusters = elliptical_clusters + 1;
                convex_hull_failed = convex_hull_failed + 1;
            end
        end
        
        %% Store results for this cell
        cell_result = struct();
        cell_result.cell_name = current_cell;
        cell_result.file_name = current_file;
        cell_result.clusterStats = clusterStats;
        cell_result.density_data = density_data;
        cell_result.XY_all = XY_all;
        cell_result.labels = labels;
        cell_result.traceIDs_all = traceIDs_all;
        cell_result.centroids = centroids;
        cell_result.nClusters = nClusters;
        cell_result.nPoints = nPoints;
        cell_result.nTraces = nTraces;
        cell_result.medianClusterArea = medianClusterArea;
        cell_result.meanClusterArea = meanClusterArea;
        cell_result.stdClusterArea = stdClusterArea;
        cell_result.medianDensity = medianDensity;
        cell_result.meanDensity = meanDensity;
        cell_result.stdDensity = stdDensity;
        cell_result.medianNND = medianNND;
        cell_result.meanNND = meanNND;
        cell_result.stdNND = stdNND;
        cell_result.convex_hull_clusters = convex_hull_clusters;
        cell_result.elliptical_clusters = elliptical_clusters;
        cell_result.exportFolder = exportFolder;
        
        % Store cell-level statistics
        cell_stats.areas_data{cell_idx} = clusterStats.areas;
        if ~isempty(clusterStats.nnd)
            cell_stats.nnd_data{cell_idx} = clusterStats.nnd;
        else
            cell_stats.nnd_data{cell_idx} = [];
        end
        cell_stats.traceIDs_data{cell_idx} = clusterStats.nTraces;
        valid_densities_cell = clusterStats.densities(~isnan(clusterStats.densities));
        cell_stats.densities_data{cell_idx} = valid_densities_cell;
        
        % Add to group results
        % FIX: Check if group_results is empty, initialize properly
        if isempty(group_results)
            group_results = cell_result;
        else
            group_results(end+1) = cell_result;
        end
        
        % Collect data for group-level statistics
        group_all_areas = [group_all_areas; clusterStats.areas];
        if ~isempty(clusterStats.nnd)
            group_all_nnd = [group_all_nnd; clusterStats.nnd];
        end
        group_all_traceIDs = [group_all_traceIDs; clusterStats.nTraces];
        valid_densities = clusterStats.densities(~isnan(clusterStats.densities));
        group_all_densities = [group_all_densities; valid_densities];
        
        %% 10. Save All Clustering Results to Files for current cell
        % Create detailed clustering results table
        detailedClusterResults = table((1:nClusters)', ...
                                       clusterStats.centroids(:,1), ...
                                       clusterStats.centroids(:,2), ...
                                       clusterStats.nTraces, ...
                                       clusterStats.nPoints, ...
                                       clusterStats.sigmaX, ...
                                       clusterStats.sigmaY, ...
                                       clusterStats.areas, ...
                                       clusterStats.densities);
        
        % Add area method as a cell array to the table
        area_method_cell = clusterStats.area_method;
        detailedClusterResults.AreaMethod = area_method_cell;
        
        if nClusters > 1 && ~isempty(clusterStats.nnd)
            detailedClusterResults.NND_nm = clusterStats.nnd;
        else
            detailedClusterResults.NND_nm = NaN(nClusters, 1);
        end
        
        detailedClusterResults.Properties.VariableNames = {'ClusterID', 'CentroidX_nm', 'CentroidY_nm', ...
            'nTraceIDs', 'nPoints', 'SigmaX_nm', 'SigmaY_nm', 'Area_nm2', 'Density_nm2', 'AreaMethod', 'NND_nm'};
        
        % Save to Excel
        excelFile = fullfile(exportFolder, [current_cell '_DensityPeak_Kmeans_Results.xlsx']);
        writetable(detailedClusterResults, excelFile);
        fprintf('\nDetailed cluster results saved to: %s\n', excelFile);
        
        % Save all localization points data with cluster assignments
        locPointsDataWithClusters = table(traceIDs_all(:), ...
                                           locations_final_nm(:,1), ...
                                           locations_final_nm(:,2), ...
                                           locations_final_nm(:,3), ...
                                           ones(size(locations_final_nm,1),1), ... % LocalizationIndex
                                           labels(:));
        
        locPointsDataWithClusters.Properties.VariableNames = {'TraceID', 'X_nm', 'Y_nm', 'Z_nm', ...
            'LocalizationIndex', 'ClusterID'};
        csvFileWithClusters = fullfile(exportFolder, [current_cell '_AllLocalizations_DensityKmeans.csv']);
        writetable(locPointsDataWithClusters, csvFileWithClusters);
        fprintf('All localization points with cluster assignments saved to: %s\n', csvFileWithClusters);
        

        
        %% 11. Create Visualization Charts for current cell 
        fprintf('Creating enhanced visualization figures for %s...\n', current_cell);
        
        fig1 = figure('Position', [50, 50, 1600, 900], 'Color', 'w', ...
            'Name', sprintf('Clustering Results - %s - ALL POINTS', current_cell));
        
        % Plot parameters
        barFaceColor = [0.2, 0.6, 0.8];
        barEdgeColor = 'k';
        medianLineColor = 'r';
        medianLineStyle = '--';
        numBins = 30;
        
        % --- Subplot 1: XY Projection colored by cluster (ALL POINTS) ---
        subplot(2, 3, 1);
        if nClusters > 0
            hsv_colors = jet(nClusters);
            random_idx = randperm(nClusters);
            colors = hsv_colors(random_idx, :);
        
            hold on;
        
            for k = 1:nClusters
                clusterPoints = XY_all(labels == k, :);
                scatter(clusterPoints(:,1), clusterPoints(:,2), ...
                    xy_point_size, colors(k,:), ...
                    'filled', 'MarkerFaceAlpha', 0.7);
            end
        
            scatter(clusterStats.centroids(:,1), clusterStats.centroids(:,2), ...
                xy_center_size, 'k', 'x', ...
                'LineWidth', xy_center_linewidth);
        
            xlabel('X (nm)');
            ylabel('Y (nm)');
            title(sprintf('XY Projection - %d Clusters\nALL %d Points', nClusters, nPoints));
            axis equal tight;
            set(gca, 'YDir', 'reverse');
            grid off;
            hold off;
        end
        
        % --- Subplot 2: Cluster trace ID distribution ---
        subplot(2, 3, 2);
        if nClusters > 0
            h = histogram(clusterStats.nTraces, numBins, ...
                'FaceColor', barFaceColor, 'EdgeColor', barEdgeColor);
            
            maxCount = max(h.Values);
            if maxCount > 0
                ylim([0, maxCount * 1.1]);
            end
            
            xlabel('Trace IDs per Cluster');
            ylabel('Frequency');
            title(sprintf('Cluster Trace ID Distribution\nMedian = %.1f trace IDs', median(clusterStats.nTraces)));
            grid off;
            hold on;
            yLimits = ylim;
            plot([median(clusterStats.nTraces), median(clusterStats.nTraces)], yLimits, ...
                medianLineColor, 'LineStyle', medianLineStyle, 'LineWidth', 2);
            hold off;
        end
        
        % --- Subplot 3: Cluster area histogram ---
        subplot(2, 3, 3);
        if nClusters > 0 && ~all(isnan(clusterStats.areas))
            validAreas = clusterStats.areas(~isnan(clusterStats.areas));
            if ~isempty(validAreas)
                h = histogram(validAreas, numBins, 'FaceColor', barFaceColor, 'EdgeColor', barEdgeColor);
                
                maxCount = max(h.Values);
                if maxCount > 0
                    ylim([0, maxCount * 1.1]);
                end
                
                xlabel('Cluster Area (nm^2)');
                ylabel('Frequency');
                title(sprintf('Cluster Areas\nMedian = %.1f nm^2', medianClusterArea));
                grid off;
                hold on;
                if ~isnan(medianClusterArea)
                    yLimits = ylim;
                    plot([medianClusterArea, medianClusterArea], yLimits, ...
                        medianLineColor, 'LineStyle', medianLineStyle, 'LineWidth', 2);
                end
                hold off;
            end
        end
        
        % --- Subplot 4: NND histogram ---
        subplot(2, 3, 4);
        if nClusters > 1 && ~isempty(clusterStats.nnd) && ~all(isnan(clusterStats.nnd))
            validNND = clusterStats.nnd(~isnan(clusterStats.nnd));
            if ~isempty(validNND)
                h = histogram(validNND, numBins, 'FaceColor', barFaceColor, 'EdgeColor', barEdgeColor);
                
                maxCount = max(h.Values);
                if maxCount > 0
                    ylim([0, maxCount * 1.1]);
                end
                
                xlabel('Nearest Neighbor Distance (nm)');
                ylabel('Frequency');
                title(sprintf('Nearest Neighbor Distances\nMedian = %.1f nm', medianNND));
                grid off;
                hold on;
                if ~isnan(medianNND)
                    yLimits = ylim;
                    plot([medianNND, medianNND], yLimits, ...
                        medianLineColor, 'LineStyle', medianLineStyle, 'LineWidth', 2);
                end
                hold off;
            end
        end
        
        % --- Subplot 5: Trace size distribution (points per trace) ---
        subplot(2, 3, 5);
        % Calculate points per trace
        uniqueTraceIDs = unique(traceIDs_all);
        points_per_trace = zeros(length(uniqueTraceIDs), 1);
        for i = 1:length(uniqueTraceIDs)
            points_per_trace(i) = sum(traceIDs_all == uniqueTraceIDs(i));
        end
        
        maxTraceSizeToShow = min(100, max(points_per_trace));
        traceSizesFiltered = points_per_trace(points_per_trace <= maxTraceSizeToShow);
        
        if ~isempty(traceSizesFiltered)
            h = histogram(traceSizesFiltered, 20, ...
                'BinLimits', [0, maxTraceSizeToShow], ...
                'FaceColor', barFaceColor, 'EdgeColor', barEdgeColor);
            
            maxCount = max(h.Values);
            if maxCount > 0
                ylim([0, maxCount * 1.1]);
            end
            
            xlabel('Points per Trace');
            ylabel('Frequency');
            title(sprintf('Trace Size Distribution\nMedian = %.1f points/trace', median(points_per_trace)));
            grid off;
            hold on;
            yLimits = ylim;
            plot([median(points_per_trace), median(points_per_trace)], yLimits, ...
                medianLineColor, 'LineStyle', medianLineStyle, 'LineWidth', 2);
            hold off;
        else
            text(0.5, 0.5, 'No trace size data', ...
                'HorizontalAlignment', 'center', 'Units', 'normalized');
            xlabel('Trace Size');
            ylabel('Frequency');
            title('Trace Size Distribution');
        end
        
        % Save Figure 1 in multiple formats
        fig1File = fullfile(exportFolder, [current_cell '_ClusteringResults_ALLPOINTS.fig']);
        saveas(fig1, fig1File);
        png1File = fullfile(exportFolder, [current_cell '_ClusteringResults_ALLPOINTS.png']);
        saveas(fig1, png1File, 'png');
        fprintf('Main figure saved as: %s\n', png1File);
        
        %% 12. Export PNG subplots to PDF with identical appearance for current cell
        fprintf('\nCreating PDF subplots for %s...\n', current_cell);
        fprintf('All subplots will match PNG exactly, including color, size, and style\n');
        
        % Reference the main figure
        main_fig = fig1;
        
        % Find all axes objects in the main figure
        subplot_axes = findobj(main_fig, 'Type', 'axes');
        n_subplots = length(subplot_axes);
        
        % Titles for PDF subplots (5 subplots)
        subplot_titles = {
            'XY_Projection_by_Cluster_ALL_POINTS';
            'Cluster_Trace_ID_Distribution';
            'Cluster_Areas_Histogram';
            'NND_Histogram';
            'Trace_Size_Distribution'
        };
        
        % Ensure we have titles for all subplots
        if length(subplot_titles) < n_subplots
            for i = length(subplot_titles)+1:n_subplots
                subplot_titles{i} = sprintf('Subplot_%d', i);
            end
        end
        
        % Loop through each subplot and create a PDF
        for i = 1:n_subplots
            % Create a new invisible figure for each subplot
            fig_sub = figure('Position', [100, 100, 800, 600], 'Color', 'w', 'Visible', 'off');
            
            % Copy the axes from the main figure
            ax_copy = copyobj(subplot_axes(i), fig_sub);
            set(ax_copy, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.75, 0.75]);
                
            % --- Histogram subplots (2-5) ---
            if i >= 2 && i <= 5
                % Set histogram bars face and edge color
                hist_objs = findobj(ax_copy, 'Type', 'histogram');
                for hh = 1:length(hist_objs)
                    set(hist_objs(hh), 'FaceColor', barFaceColor, 'EdgeColor', barEdgeColor);
                end
                
                % Set median lines style
                line_objs = findobj(ax_copy, 'Type', 'line');
                for ll = 1:length(line_objs)
                    if strcmp(get(line_objs(ll), 'Color'), medianLineColor) && ...
                       strcmp(get(line_objs(ll), 'LineStyle'), medianLineStyle)
                        set(line_objs(ll), 'LineWidth', 2);
                    end
                end
            end
            
            % --- Set title safely ---
            title_text = get(get(ax_copy, 'Title'), 'String');
            if isempty(title_text)
                title_text = subplot_titles{i};
            elseif iscell(title_text)
                title_text = title_text{1};
            end
            safe_title = regexprep(title_text, '[\\/:*?"<>|]', '_');
            safe_title = regexprep(safe_title, '\s+', '_');
            safe_title = regexprep(safe_title, '\n', '_');
            
            % --- Save PDF ---
            subplot_pdf_file = fullfile(exportFolder, sprintf('%s_Subplot%d_%s.pdf', current_cell, i, safe_title));
            set(fig_sub, 'PaperPositionMode', 'auto');
            print(fig_sub, subplot_pdf_file, '-dpdf', '-bestfit');
            
            fprintf('  Created PDF subplot %d: %s\n', i, safe_title);
            
            % Close temporary figure
            close(fig_sub);
        end
        
        fprintf('\nPDF subplots created successfully for %s!\n', current_cell);
        
        % Close the main figure for this cell
        close(fig1);
        
        %% 13. Create Summary Report with Parameter Documentation for current cell
        summaryFile = fullfile(exportFolder, [current_cell '_DensityPeak_Analysis_Summary_ALLPOINTS.txt']);
        fid = fopen(summaryFile, 'w');
        fprintf(fid, '==========================================\n');
        fprintf(fid, 'DENSITY-PEAK K-MEANS CLUSTERING ANALYSIS SUMMARY\n');
        fprintf(fid, 'CUSTOM PARAMETERS VERSION - ALL POINTS CLUSTERING\n');
        fprintf(fid, '==========================================\n\n');
        fprintf(fid, 'Cell: %s\n', current_cell);
        fprintf(fid, 'Data file: %s\n', current_file);
        fprintf(fid, 'Analysis date: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        fprintf(fid, 'Clustering method: Density-Peak Based K-means on ALL POINTS\n');
        fprintf(fid, 'Parameter mode: Custom parameters (direct modification in script)\n\n');
        
        fprintf(fid, 'CUSTOM PARAMETERS USED:\n');
        fprintf(fid, '  Density pixel size: %.1f nm\n', density_pixel_size);
        fprintf(fid, '  Gaussian sigma for smoothing: %.1f pixels\n', gaussian_sigma);
        fprintf(fid, '  Density threshold factor: %.2f\n', density_threshold_factor);
        fprintf(fid, '  Minimum peak distance: %.1f nm\n', min_peak_distance);
       
        
        fprintf(fid, 'FILTERING RESULTS:\n');
        fprintf(fid, '  Step 1 - Valid localizations (vld==1): %d\n', initial_count);
        fprintf(fid, '  Step 2 - After cfr filter: %d\n', step2_count);
        fprintf(fid, '  Step 3 - After efo filter: %d\n', step3_count);
        fprintf(fid, '  Step 4 - Final: %d traces, %d localizations\n\n', trace_count, final_count);
        
        fprintf(fid, 'CLUSTERING RESULTS:\n');
        fprintf(fid, '  Total localizations clustered: %d\n', nPoints);
        fprintf(fid, '  Total unique trace IDs: %d\n', nTraces);
        fprintf(fid, '  Initial density peaks found: %d\n', n_peaks_final);
        fprintf(fid, '  Final clusters identified: %d\n', nClusters);
        fprintf(fid, '  Average trace IDs per cluster: %.1f\n', mean(clusterStats.nTraces));
        fprintf(fid, '  Average points per cluster: %.1f\n\n', mean(clusterStats.nPoints));
        
        fprintf(fid, 'CLUSTER AREA STATISTICS (nm^2):\n');
        fprintf(fid, '  Median: %.2f\n', medianClusterArea);
        fprintf(fid, '  Mean: %.2f ± %.2f\n', meanClusterArea, stdClusterArea);
        fprintf(fid, '  Range: %.2f to %.2f\n', min(clusterStats.areas), max(clusterStats.areas));
        fprintf(fid, '  Convex Hull clusters only - Median: %.2f, Mean: %.2f\n\n', ...
            median(clusterStats.areas(strcmp(clusterStats.area_method, 'Convex Hull'))), ...
            mean(clusterStats.areas(strcmp(clusterStats.area_method, 'Convex Hull'))));
        
        fprintf(fid, 'CLUSTER DENSITY STATISTICS (trace IDs/nm^2):\n');
        if ~isnan(medianDensity)
            fprintf(fid, '  Median: %.4f\n', medianDensity);
            fprintf(fid, '  Mean: %.4f ± %.4f\n', meanDensity, stdDensity);
        else
            fprintf(fid, '  Density calculation not available\n');
        end
        fprintf(fid, '\n');
        
        fprintf(fid, 'NEAREST NEIGHBOR DISTANCE STATISTICS (nm):\n');
        if nClusters > 1
            fprintf(fid, '  Median: %.2f\n', medianNND);
            fprintf(fid, '  Mean: %.2f ± %.2f\n', meanNND, stdNND);
            fprintf(fid, '  Range: %.2f to %.2f\n\n', min(clusterStats.nnd), max(clusterStats.nnd));
        else
            fprintf(fid, '  Only one cluster found, NND not applicable\n\n');
        end
        
        fprintf(fid, 'OUTPUT FILES:\n');
        fprintf(fid, '  1. All localization points: %s\n', csvFile);
        fprintf(fid, '  2. All localization points with cluster assignments: %s\n', csvFileWithClusters);
        fprintf(fid, '  3. Detailed clustering results (including area method): %s\n', excelFile);
        fprintf(fid, '  5. Main visualization figure (PNG): %s\n', png1File);
        fprintf(fid, '  6. PDF subplots (identical to PNG subplots): %s\n', exportFolder);
        fprintf(fid, '  7. This summary report: %s\n\n', summaryFile);
        
        fclose(fid);
        fprintf('Summary report saved to: %s\n', summaryFile);
        
        fprintf('\n==========================================\n');
        fprintf('COMPLETED PROCESSING FOR CELL %d/%d: %s\n', cell_idx, nCells, current_cell);
        fprintf('==========================================\n\n');
    end
    
    %% Store group results and cell statistics
    all_groups_results{group_idx} = group_results;
    group_cell_stats{group_idx} = cell_stats;
    
    % Store group summary statistics
    group_summary(group_idx).group_name = current_group;
    group_summary(group_idx).nCells = nCells;
    group_summary(group_idx).total_clusters = sum([group_results.nClusters]);
    group_summary(group_idx).total_points = sum([group_results.nPoints]);
    group_summary(group_idx).total_traces = sum([group_results.nTraces]);
    group_summary(group_idx).median_area = median(group_all_areas);
    group_summary(group_idx).mean_area = mean(group_all_areas);
    group_summary(group_idx).std_area = std(group_all_areas);
    group_summary(group_idx).median_nnd = median(group_all_nnd);
    group_summary(group_idx).mean_nnd = mean(group_all_nnd);
    group_summary(group_idx).std_nnd = std(group_all_nnd);
    group_summary(group_idx).median_traceIDs = median(group_all_traceIDs);
    group_summary(group_idx).mean_traceIDs = mean(group_all_traceIDs);
    group_summary(group_idx).std_traceIDs = std(group_all_traceIDs);
    group_summary(group_idx).median_density = median(group_all_densities);
    group_summary(group_idx).mean_density = mean(group_all_densities);
    group_summary(group_idx).std_density = std(group_all_densities);
    group_summary(group_idx).all_areas = group_all_areas;
    group_summary(group_idx).all_nnd = group_all_nnd;
    group_summary(group_idx).all_traceIDs = group_all_traceIDs;
    group_summary(group_idx).all_densities = group_all_densities;
    
    fprintf('\n==========================================\n');
    fprintf('COMPLETED PROCESSING FOR GROUP %d/%d: %s\n', group_idx, nGroups, current_group);
    fprintf('==========================================\n\n');
end

%% 14. Create Group-Level Analysis: Individual Cell Statistics and Combined Histograms 
fprintf('\n==========================================\n');
fprintf('CREATING GROUP-LEVEL ANALYSIS\n');
fprintf('==========================================\n');

% Create a folder for group-level analysis in the specified path
base_path = 'C:\Users\18317\Desktop\code 20260225';
groupAnalysisFolder = fullfile(base_path, 'export_GroupLevel_Analysis');
if ~exist(groupAnalysisFolder, 'dir')
    mkdir(groupAnalysisFolder);
end

% Define bin configurations for different metrics
bin_configs = struct();
bin_configs.areas = struct('n_bins', 20);
bin_configs.nnd = struct('n_bins', 20);
bin_configs.traceIDs = struct('n_bins', 20);
% bin_configs.densities removed (no longer plotted)

% Initialize combined data storage
combined_data = struct();
combined_data.areas = cell(nGroups, 1);
combined_data.nnd = cell(nGroups, 1);
combined_data.traceIDs = cell(nGroups, 1);
combined_data.densities = cell(nGroups, 1); % still stored for summary, but not plotted

%% Process each group separately for individual group analysis
for group_idx = 1:nGroups
    current_group = group_names{group_idx};
    cell_stats = group_cell_stats{group_idx};
    nCells = length(cell_stats.cell_names);
    
    fprintf('\nCreating group-level analysis for %s (%d cells)...\n', current_group, nCells);
    
    % Create a figure for this group's analysis (3 metrics in one figure)
    fig_group = figure('Position', [50, 50, 1600, 600], 'Color', 'w', ...
        'Name', sprintf('%s Group Analysis - Individual Cell Statistics', current_group));
    
    % 1. Cluster Areas Analysis
    subplot(1, 3, 1);
    hold on;
    
    % Collect all area data for this group
    all_areas_group = [];
    for cell_idx = 1:nCells
        if ~isempty(cell_stats.areas_data{cell_idx})
            all_areas_group = [all_areas_group; cell_stats.areas_data{cell_idx}];
        end
    end
    
    if ~isempty(all_areas_group)
        % Determine bin edges
        max_area = max(all_areas_group);
        min_area = min(all_areas_group);
        bin_edges = linspace(min_area, max_area, bin_configs.areas.n_bins + 1);
        bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
        
        % Calculate histogram for each cell
        cell_hist_counts = zeros(length(bin_centers), nCells);
        for cell_idx = 1:nCells
            areas_data = cell_stats.areas_data{cell_idx};
            if ~isempty(areas_data)
                [counts, ~] = histcounts(areas_data, bin_edges);
                cell_hist_counts(:, cell_idx) = counts;
            end
        end
        
        % Calculate mean and SEM
        mean_counts = mean(cell_hist_counts, 2, 'omitnan');
        sem_counts = std(cell_hist_counts, 0, 2, 'omitnan') / sqrt(nCells);
        
        % Create bar plot with unified color
        bar_positions = 1:length(bin_centers);
        bar_width = 0.8;
        
        % Create bar plot with unified color - changed to [0.2, 0.6, 0.8]
        h_bar = bar(bar_positions, mean_counts, bar_width, ...
            'FaceColor', [0.2, 0.6, 0.8], ...  % Modified to match first script
            'EdgeColor', 'k', ...
            'LineWidth', 1.0);
        
        % Add error bars manually (without CapSize)
        for i = 1:length(bar_positions)
            x_pos = bar_positions(i);
            y_pos = mean_counts(i);
            error_length = sem_counts(i);
            
            % Draw vertical error bar
            line([x_pos, x_pos], [y_pos - error_length, y_pos + error_length], ...
                'Color', 'k', 'LineWidth', 1.5);
            
            % Draw horizontal caps
            cap_width = bar_width * 0.2;
            line([x_pos - cap_width, x_pos + cap_width], ...
                [y_pos - error_length, y_pos - error_length], ...
                'Color', 'k', 'LineWidth', 1.5);
            line([x_pos - cap_width, x_pos + cap_width], ...
                [y_pos + error_length, y_pos + error_length], ...
                'Color', 'k', 'LineWidth', 1.5);
        end
        
        % Format x-axis with integer values
        x_ticks = [];
        x_tick_labels = {};
        
        % Find appropriate tick positions based on data range
        if max_area <= 100
            tick_step = ceil(max_area/10/10)*10;
            if tick_step < 10
                tick_step = 10;
            end
            x_ticks = 0:tick_step:ceil(max_area);
        elseif max_area <= 1000
            tick_step = ceil(max_area/10/50)*50;
            if tick_step < 50
                tick_step = 50;
            end
            x_ticks = 0:tick_step:ceil(max_area);
        else
            if max_area <= 2000
                tick_step = ceil(max_area/10/100)*100;
                if tick_step < 100
                    tick_step = 100;
                end
                x_ticks = 0:tick_step:ceil(max_area);
            else
                tick_step = ceil(max_area/10/500)*500;
                x_ticks = 0:tick_step:ceil(max_area);
            end
        end
        
        % Convert x_ticks to bin positions
        x_tick_positions = zeros(1, length(x_ticks));
        x_tick_labels = cell(1, length(x_ticks));
        
        for i = 1:length(x_ticks)
            tick_value = x_ticks(i);
            [~, closest_idx] = min(abs(bin_centers - tick_value));
            x_tick_positions(i) = bar_positions(closest_idx);
            x_tick_labels{i} = sprintf('%d', tick_value);
        end
        
        set(gca, 'XTick', x_tick_positions, 'XTickLabel', x_tick_labels);
        
        % Set Y axis to show only >=0
        y_lim = ylim();
        ylim([0, y_lim(2)]);
        
        xlabel('Cluster Area (nm²)', 'FontSize', 12);
        ylabel('Frequency', 'FontSize', 12);
        title(sprintf('%s: Cluster Area Distribution (n=%d cells)', current_group, nCells), 'FontSize', 14);
        grid off;  % Turn off grid
        box on;    % Add box
        
        % Add red dashed line for median of all combined data
        if ~isempty(all_areas_group)
            combined_median = median(all_areas_group);
            [~, median_bin_idx] = min(abs(bin_centers - combined_median));
            median_x_pos = bar_positions(median_bin_idx);
            
            y_limits = ylim();
            plot([median_x_pos, median_x_pos], [0, y_limits(2)*0.95], ...
                'r--', 'LineWidth', 2);
            
            text(median_x_pos, y_limits(2)*0.97, sprintf('Median: %.1f', combined_median), ...
                'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'BackgroundColor', 'w');
        end
        
        hold off;
    else
        text(0.5, 0.5, 'No area data available', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 12);
        xlabel('Cluster Area (nm²)');
        ylabel('Frequency');
        title(sprintf('%s: Cluster Area Distribution', current_group));
        grid off;
        box on;
    end
    
    % 2. NND Analysis
    subplot(1, 3, 2);
    hold on;
    
    % Collect all NND data for this group
    all_nnd_group = [];
    for cell_idx = 1:nCells
        nnd_data = cell_stats.nnd_data{cell_idx};
        if ~isempty(nnd_data)
            all_nnd_group = [all_nnd_group; nnd_data];
        end
    end
    
    if ~isempty(all_nnd_group)
        % Determine bin edges
        max_nnd = max(all_nnd_group);
        min_nnd = min(all_nnd_group);
        bin_edges = linspace(min_nnd, max_nnd, bin_configs.nnd.n_bins + 1);
        bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
        
        % Calculate histogram for each cell
        cell_hist_counts = zeros(length(bin_centers), nCells);
        for cell_idx = 1:nCells
            nnd_data = cell_stats.nnd_data{cell_idx};
            if ~isempty(nnd_data)
                [counts, ~] = histcounts(nnd_data, bin_edges);
                cell_hist_counts(:, cell_idx) = counts;
            end
        end
        
        % Calculate mean and SEM
        mean_counts = mean(cell_hist_counts, 2, 'omitnan');
        sem_counts = std(cell_hist_counts, 0, 2, 'omitnan') / sqrt(nCells);
        
        % Create bar plot with unified color
        bar_positions = 1:length(bin_centers);
        bar_width = 0.8;
        
        h_bar = bar(bar_positions, mean_counts, bar_width, ...
            'FaceColor', [0.2, 0.6, 0.8], ...  % Modified to match first script
            'EdgeColor', 'k', ...
            'LineWidth', 1.0);
        
        % Add error bars manually
        for i = 1:length(bar_positions)
            x_pos = bar_positions(i);
            y_pos = mean_counts(i);
            error_length = sem_counts(i);
            
            line([x_pos, x_pos], [y_pos - error_length, y_pos + error_length], ...
                'Color', 'k', 'LineWidth', 1.5);
            
            cap_width = bar_width * 0.2;
            line([x_pos - cap_width, x_pos + cap_width], ...
                [y_pos - error_length, y_pos - error_length], ...
                'Color', 'k', 'LineWidth', 1.5);
            line([x_pos - cap_width, x_pos + cap_width], ...
                [y_pos + error_length, y_pos + error_length], ...
                'Color', 'k', 'LineWidth', 1.5);
        end
        
        % Format x-axis with integer values
        x_ticks = [];
        x_tick_labels = {};
        
        if max_nnd <= 100
            tick_step = ceil(max_nnd/10/10)*10;
            if tick_step < 10
                tick_step = 10;
            end
            x_ticks = 0:tick_step:ceil(max_nnd);
        elseif max_nnd <= 1000
            tick_step = ceil(max_nnd/10/50)*50;
            if tick_step < 50
                tick_step = 50;
            end
            x_ticks = 0:tick_step:ceil(max_nnd);
        else
            if max_nnd <= 2000
                tick_step = ceil(max_nnd/10/100)*100;
                if tick_step < 100
                    tick_step = 100;
                end
                x_ticks = 0:tick_step:ceil(max_nnd);
            else
                tick_step = ceil(max_nnd/10/500)*500;
                x_ticks = 0:tick_step:ceil(max_nnd);
            end
        end
        
        x_tick_positions = zeros(1, length(x_ticks));
        x_tick_labels = cell(1, length(x_ticks));
        
        for i = 1:length(x_ticks)
            tick_value = x_ticks(i);
            [~, closest_idx] = min(abs(bin_centers - tick_value));
            x_tick_positions(i) = bar_positions(closest_idx);
            x_tick_labels{i} = sprintf('%d', tick_value);
        end
        
        set(gca, 'XTick', x_tick_positions, 'XTickLabel', x_tick_labels);
        
        y_lim = ylim();
        ylim([0, y_lim(2)]);
        
        xlabel('Nearest Neighbor Distance (nm)', 'FontSize', 12);
        ylabel('Frequency', 'FontSize', 12);
        title(sprintf('%s: NND Distribution (n=%d cells)', current_group, nCells), 'FontSize', 14);
        grid off;  % Turn off grid
        box on;    % Add box
        
        if ~isempty(all_nnd_group)
            combined_median = median(all_nnd_group);
            [~, median_bin_idx] = min(abs(bin_centers - combined_median));
            median_x_pos = bar_positions(median_bin_idx);
            
            y_limits = ylim();
            plot([median_x_pos, median_x_pos], [0, y_limits(2)*0.95], ...
                'r--', 'LineWidth', 2);
            
            text(median_x_pos, y_limits(2)*0.97, sprintf('Median: %.1f', combined_median), ...
                'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'BackgroundColor', 'w');
        end
        
        hold off;
    else
        text(0.5, 0.5, 'No NND data available', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 12);
        xlabel('Nearest Neighbor Distance (nm)');
        ylabel('Frequency');
        title(sprintf('%s: NND Distribution', current_group));
        grid off;
        box on;
    end
    
    % 3. Cluster Size (Trace IDs per Cluster) Analysis
    subplot(1, 3, 3);
    hold on;
    
    % Collect all traceIDs data for this group
    all_traceIDs_group = [];
    for cell_idx = 1:nCells
        traceIDs_data = cell_stats.traceIDs_data{cell_idx};
        if ~isempty(traceIDs_data)
            all_traceIDs_group = [all_traceIDs_group; traceIDs_data];
        end
    end
    
    if ~isempty(all_traceIDs_group)
        % Determine bin edges
        max_traceIDs = max(all_traceIDs_group);
        min_traceIDs = min(all_traceIDs_group);
        bin_edges = linspace(min_traceIDs, max_traceIDs, bin_configs.traceIDs.n_bins + 1);
        bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
        
        % Calculate histogram for each cell
        cell_hist_counts = zeros(length(bin_centers), nCells);
        for cell_idx = 1:nCells
            traceIDs_data = cell_stats.traceIDs_data{cell_idx};
            if ~isempty(traceIDs_data)
                [counts, ~] = histcounts(traceIDs_data, bin_edges);
                cell_hist_counts(:, cell_idx) = counts;
            end
        end
        
        % Calculate mean and SEM
        mean_counts = mean(cell_hist_counts, 2, 'omitnan');
        sem_counts = std(cell_hist_counts, 0, 2, 'omitnan') / sqrt(nCells);
        
        % Create bar plot with unified color
        bar_positions = 1:length(bin_centers);
        bar_width = 0.8;
        
        h_bar = bar(bar_positions, mean_counts, bar_width, ...
            'FaceColor', [0.2, 0.6, 0.8], ...  % Modified to match first script
            'EdgeColor', 'k', ...
            'LineWidth', 1.0);
        
        % Add error bars manually
        for i = 1:length(bar_positions)
            x_pos = bar_positions(i);
            y_pos = mean_counts(i);
            error_length = sem_counts(i);
            
            line([x_pos, x_pos], [y_pos - error_length, y_pos + error_length], ...
                'Color', 'k', 'LineWidth', 1.5);
            
            cap_width = bar_width * 0.2;
            line([x_pos - cap_width, x_pos + cap_width], ...
                [y_pos - error_length, y_pos - error_length], ...
                'Color', 'k', 'LineWidth', 1.5);
            line([x_pos - cap_width, x_pos + cap_width], ...
                [y_pos + error_length, y_pos + error_length], ...
                'Color', 'k', 'LineWidth', 1.5);
        end
        
        % Format x-axis with integer values
        x_ticks = [];
        x_tick_labels = {};
        
        if max_traceIDs <= 20
            tick_step = ceil(max_traceIDs/10/2)*2;
            if tick_step < 2
                tick_step = 2;
            end
            x_ticks = 0:tick_step:ceil(max_traceIDs);
        elseif max_traceIDs <= 100
            tick_step = ceil(max_traceIDs/10/10)*10;
            if tick_step < 10
                tick_step = 10;
            end
            x_ticks = 0:tick_step:ceil(max_traceIDs);
        else
            tick_step = ceil(max_traceIDs/10/50)*50;
            x_ticks = 0:tick_step:ceil(max_traceIDs);
        end
        
        x_tick_positions = zeros(1, length(x_ticks));
        x_tick_labels = cell(1, length(x_ticks));
        
        for i = 1:length(x_ticks)
            tick_value = x_ticks(i);
            [~, closest_idx] = min(abs(bin_centers - tick_value));
            x_tick_positions(i) = bar_positions(closest_idx);
            x_tick_labels{i} = sprintf('%d', tick_value);
        end
        
        set(gca, 'XTick', x_tick_positions, 'XTickLabel', x_tick_labels);
        
        y_lim = ylim();
        ylim([0, y_lim(2)]);
        
        xlabel('Trace IDs per Cluster', 'FontSize', 12);
        ylabel('Frequency', 'FontSize', 12);
        title(sprintf('%s: Cluster Size Distribution (n=%d cells)', current_group, nCells), 'FontSize', 14);
        grid off;  % Turn off grid
        box on;    % Add box
        
        if ~isempty(all_traceIDs_group)
            combined_median = median(all_traceIDs_group);
            [~, median_bin_idx] = min(abs(bin_centers - combined_median));
            median_x_pos = bar_positions(median_bin_idx);
            
            y_limits = ylim();
            plot([median_x_pos, median_x_pos], [0, y_limits(2)*0.95], ...
                'r--', 'LineWidth', 2);
            
            text(median_x_pos, y_limits(2)*0.97, sprintf('Median: %.1f', combined_median), ...
                'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'BackgroundColor', 'w');
        end
        
        hold off;
    else
        text(0.5, 0.5, 'No cluster size data available', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 12);
        xlabel('Trace IDs per Cluster');
        ylabel('Frequency');
        title(sprintf('%s: Cluster Size Distribution', current_group));
        grid off;
        box on;
    end
    
    % Save the group analysis figure
    group_png = fullfile(groupAnalysisFolder, sprintf('%s_Group_Analysis.png', current_group));
    saveas(fig_group, group_png, 'png');
    group_pdf = fullfile(groupAnalysisFolder, sprintf('%s_Group_Analysis.pdf', current_group));
    saveas(fig_group, group_pdf, 'pdf');
    fprintf('Group analysis figure saved for %s: %s\n', current_group, group_png);
    
    %% ===== Create PDF subplots for each metric =====
    fprintf('Creating PDF subplots for %s...\n', current_group);
    
    % Create folder for PDF subplots
    pdfSubplotsFolder = fullfile(groupAnalysisFolder, sprintf('%s_PDF_Subplots', current_group));
    if ~exist(pdfSubplotsFolder, 'dir')
        mkdir(pdfSubplotsFolder);
    end
    
    % Get all axes from the group figure
    all_axes = findobj(fig_group, 'Type', 'axes');
    
    % Titles for PDF subplots (3 subplots)
    subplot_titles = {
        'Cluster_Area_Distribution';
        'NND_Distribution';
        'Cluster_Size_Distribution';
    };
    
    % Ensure we have titles for all subplots
    if length(subplot_titles) < length(all_axes)
        for i = length(subplot_titles)+1:length(all_axes)
            subplot_titles{i} = sprintf('Subplot_%d', i);
        end
    end
    
    % Create PDF for each subplot
    for subplot_idx = 1:length(all_axes)
        % Create new figure
        fig_single = figure('Position', [100, 100, 800, 600], ...
                           'Color', 'w', ...
                           'Visible', 'off');
        
        % Copy axes
        ax_copy = copyobj(all_axes(subplot_idx), fig_single);
        
        % Adjust position
        set(ax_copy, 'Units', 'normalized', ...
                     'Position', [0.15, 0.15, 0.75, 0.75]);
        
        % Improve appearance
        % Title
        original_title = get(get(ax_copy, 'Title'), 'String');
        if ~isempty(original_title)
            title_obj = get(ax_copy, 'Title');
            set(title_obj, 'FontSize', 14, 'FontWeight', 'bold');
        end
        
        % Axis labels
        xlabel_obj = get(ax_copy, 'XLabel');
        if ~isempty(xlabel_obj)
            set(xlabel_obj, 'FontSize', 12);
        end
        
        ylabel_obj = get(ax_copy, 'YLabel');
        if ~isempty(ylabel_obj)
            set(ylabel_obj, 'FontSize', 12);
        end
        
        % Grid (copy will retain original settings, but we can enforce)
        % No need to modify here since original already has grid off/box on
        
        % Determine title
        if isempty(original_title)
            title_text = subplot_titles{subplot_idx};
        elseif iscell(original_title)
            title_text = original_title{1};
        else
            title_text = original_title;
        end
        
        % Sanitize title for filename
        safe_title = regexprep(title_text, '[\\/:*?"<>|]', '_');
        safe_title = regexprep(safe_title, '\s+', '_');
        safe_title = regexprep(safe_title, '\n', '_');
        
        % Set figure name
        set(fig_single, 'Name', safe_title);
        
        % Save PDF
        pdf_file = fullfile(pdfSubplotsFolder, ...
                           sprintf('%s_Subplot%d_%s.pdf', current_group, subplot_idx, safe_title));
        
        set(fig_single, 'PaperPositionMode', 'auto');
        print(fig_single, pdf_file, '-dpdf', '-bestfit');
        
        % Save PNG as well
        png_file = fullfile(pdfSubplotsFolder, ...
                           sprintf('%s_Subplot%d_%s.png', current_group, subplot_idx, safe_title));
        print(fig_single, png_file, '-dpng', '-r300');
        
        fprintf('  Created PDF subplot %d: %s\n', subplot_idx, safe_title);
        
        % Close temporary figure
        close(fig_single);
    end
    
    fprintf('PDF subplots saved in: %s\n\n', pdfSubplotsFolder);
    
    % Close the group figure
    close(fig_group);
    
    % Store combined data for two-group comparison
    combined_data.areas{group_idx} = all_areas_group;
    combined_data.nnd{group_idx} = all_nnd_group;
    combined_data.traceIDs{group_idx} = all_traceIDs_group;
    combined_data.densities{group_idx} = group_all_densities; % for summary only
end
%% 15. Create CSV Table with Cell-Level Median Statistics
fprintf('\n==========================================\n');
fprintf('CREATING CSV TABLE WITH CELL-LEVEL MEDIAN STATISTICS\n');
fprintf('==========================================\n');

% Initialize cell array for storing statistics
cell_stats_data = [];

for group_idx = 1:nGroups
    current_group = group_names{group_idx};
    group_results = all_groups_results{group_idx};
    
    % Check if group_results is valid
    if isempty(group_results) || ~isstruct(group_results)
        fprintf('Warning: No valid group results for %s\n', current_group);
        continue;
    end
    
    % Number of cells in this group
    nCells = length(group_results);
    
    for cell_idx = 1:nCells
        if isempty(group_results(cell_idx))
            continue;
        end
        
        try
            cell_result = group_results(cell_idx);
            
            % Check for required fields
            if ~isfield(cell_result, 'cell_name') || ~isfield(cell_result, 'medianClusterArea') || ...
               ~isfield(cell_result, 'medianNND') || ~isfield(cell_result, 'clusterStats')
                fprintf('Warning: Missing fields in cell %d of group %s\n', cell_idx, current_group);
                continue;
            end
            
            % Extract statistics
            median_area = cell_result.medianClusterArea;
            median_nnd = cell_result.medianNND;
            
            % Median trace IDs
            if isfield(cell_result.clusterStats, 'nTraces') && ~isempty(cell_result.clusterStats.nTraces)
                median_traceIDs = median(cell_result.clusterStats.nTraces);
            else
                median_traceIDs = NaN;
            end
            
            % Handle NaN
            if isempty(median_area) || isnan(median_area)
                median_area = NaN;
            end
            
            if isempty(median_nnd) || isnan(median_nnd)
                median_nnd = NaN;
            end
            
            if isempty(median_traceIDs) || isnan(median_traceIDs)
                median_traceIDs = NaN;
            end
            
            % New row - ensure numeric type
            new_row = {current_group, ...
                      cell_result.cell_name, ...
                      double(median_area), ...
                      double(median_nnd), ...
                      double(median_traceIDs)};
            
            % Append to data
            if isempty(cell_stats_data)
                cell_stats_data = new_row;
            else
                cell_stats_data = [cell_stats_data; new_row];
            end
            
        catch ME
            fprintf('Warning: Error processing cell %d in group %s: %s\n', ...
                cell_idx, current_group, ME.message);
            continue;
        end
    end
end

% Create table and save CSV
if ~isempty(cell_stats_data)
    % Check number of columns
    if size(cell_stats_data, 2) == 5
        cell_stats_table = cell2table(cell_stats_data, ...
            'VariableNames', {'Group', 'CellName', 'MedianClusterArea_nm2', 'MedianNND_nm', 'MedianTraceIDs'});
        
        % Save CSV
        csv_filename = fullfile(groupAnalysisFolder, 'CellLevel_Median_Statistics.csv');
        writetable(cell_stats_table, csv_filename);
        fprintf('Cell-level median statistics saved to: %s\n', csv_filename);
        
        % Display table
        fprintf('\nCell-Level Median Statistics:\n');
        disp(cell_stats_table);
    else
        fprintf('Warning: Data has incorrect number of columns (%d instead of 5).\n', size(cell_stats_data, 2));
        fprintf('Data content:\n');
        disp(cell_stats_data);
    end
else
    fprintf('Warning: No cell-level statistics data available for CSV export.\n');
end

%% 16. Create Two-Group Statistical Summary
fprintf('\nCreating two-group statistical summary...\n');

summary_file = fullfile(groupAnalysisFolder, 'TwoGroup_Statistical_Summary.txt');
fid = fopen(summary_file, 'w');

fprintf(fid, '==========================================\n');
fprintf(fid, 'TWO-GROUP STATISTICAL SUMMARY\n');
fprintf(fid, 'Euchromatin vs Heterochromatin\n');
fprintf(fid, '==========================================\n\n');

fprintf(fid, 'Analysis Date: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, 'Custom Parameters Used:\n');
fprintf(fid, '  Density pixel size: %.1f nm\n', density_pixel_size);
fprintf(fid, '  Gaussian sigma for smoothing: %.1f pixels\n', gaussian_sigma);
fprintf(fid, '  Density threshold factor: %.2f\n', density_threshold_factor);
fprintf(fid, '  Minimum peak distance: %.1f nm\n\n', min_peak_distance);

% Group information
fprintf(fid, 'GROUP INFORMATION:\n');
fprintf(fid, '%-15s %-10s %-15s %-15s %-15s\n', ...
    'Group', 'Cells', 'Total Clusters', 'Total Points', 'Total Traces');
fprintf(fid, '%s\n', repmat('-', 70, 1));

for group_idx = 1:nGroups
    fprintf(fid, '%-15s %-10d %-15d %-15d %-15d\n', ...
        group_summary(group_idx).group_name, ...
        group_summary(group_idx).nCells, ...
        group_summary(group_idx).total_clusters, ...
        group_summary(group_idx).total_points, ...
        group_summary(group_idx).total_traces);
end
fprintf(fid, '\n');

% Cluster Area Statistics
fprintf(fid, 'CLUSTER AREA STATISTICS (nm²):\n');
fprintf(fid, '%-15s %-10s %-15s %-15s %-15s\n', ...
    'Group', 'Median', 'Mean', 'Std', 'Range');
fprintf(fid, '%s\n', repmat('-', 70, 1));

for group_idx = 1:nGroups
    data = combined_data.areas{group_idx};
    if ~isempty(data)
        fprintf(fid, '%-15s %-10.1f %-15.1f %-15.1f %-15s\n', ...
            group_names{group_idx}, ...
            median(data), ...
            mean(data), ...
            std(data), ...
            sprintf('[%.1f, %.1f]', min(data), max(data)));
    else
        fprintf(fid, '%-15s %-10s %-15s %-15s %-15s\n', ...
            group_names{group_idx}, 'N/A', 'N/A', 'N/A', 'N/A');
    end
end
fprintf(fid, '\n');

% NND Statistics
fprintf(fid, 'NEAREST NEIGHBOR DISTANCE STATISTICS (nm):\n');
fprintf(fid, '%-15s %-10s %-15s %-15s %-15s\n', ...
    'Group', 'Median', 'Mean', 'Std', 'Range');
fprintf(fid, '%s\n', repmat('-', 70, 1));

for group_idx = 1:nGroups
    data = combined_data.nnd{group_idx};
    if ~isempty(data)
        fprintf(fid, '%-15s %-10.1f %-15.1f %-15.1f %-15s\n', ...
            group_names{group_idx}, ...
            median(data), ...
            mean(data), ...
            std(data), ...
            sprintf('[%.1f, %.1f]', min(data), max(data)));
    else
        fprintf(fid, '%-15s %-10s %-15s %-15s %-15s\n', ...
            group_names{group_idx}, 'N/A', 'N/A', 'N/A', 'N/A');
    end
end
fprintf(fid, '\n');

% Cluster Size Statistics
fprintf(fid, 'CLUSTER SIZE STATISTICS (Trace IDs per Cluster):\n');
fprintf(fid, '%-15s %-10s %-15s %-15s %-15s\n', ...
    'Group', 'Median', 'Mean', 'Std', 'Range');
fprintf(fid, '%s\n', repmat('-', 70, 1));

for group_idx = 1:nGroups
    data = combined_data.traceIDs{group_idx};
    if ~isempty(data)
        fprintf(fid, '%-15s %-10.1f %-15.1f %-15.1f %-15s\n', ...
            group_names{group_idx}, ...
            median(data), ...
            mean(data), ...
            std(data), ...
            sprintf('[%.1f, %.1f]', min(data), max(data)));
    else
        fprintf(fid, '%-15s %-10s %-15s %-15s %-15s\n', ...
            group_names{group_idx}, 'N/A', 'N/A', 'N/A', 'N/A');
    end
end
fprintf(fid, '\n');

% Density Statistics
fprintf(fid, 'CLUSTER DENSITY STATISTICS (Trace IDs/nm²):\n');
fprintf(fid, '%-15s %-10s %-15s %-15s %-15s\n', ...
    'Group', 'Median', 'Mean', 'Std', 'Range');
fprintf(fid, '%s\n', repmat('-', 70, 1));

for group_idx = 1:nGroups
    data = combined_data.densities{group_idx};
    if ~isempty(data)
        fprintf(fid, '%-15s %-10.4f %-15.4f %-15.4f %-15s\n', ...
            group_names{group_idx}, ...
            median(data), ...
            mean(data), ...
            std(data), ...
            sprintf('[%.4f, %.4f]', min(data), max(data)));
    else
        fprintf(fid, '%-15s %-10s %-15s %-15s %-15s\n', ...
            group_names{group_idx}, 'N/A', 'N/A', 'N/A', 'N/A');
    end
end
fprintf(fid, '\n');

% Combined Statistics (Both Groups)
fprintf(fid, 'COMBINED STATISTICS (Both Groups):\n');

% Combine all areas
all_areas_combined = [];
for group_idx = 1:nGroups
    if ~isempty(combined_data.areas{group_idx})
        all_areas_combined = [all_areas_combined; combined_data.areas{group_idx}];
    end
end

% Combine all NND
all_nnd_combined = [];
for group_idx = 1:nGroups
    if ~isempty(combined_data.nnd{group_idx})
        all_nnd_combined = [all_nnd_combined; combined_data.nnd{group_idx}];
    end
end

% Combine all trace IDs
all_traceIDs_combined = [];
for group_idx = 1:nGroups
    if ~isempty(combined_data.traceIDs{group_idx})
        all_traceIDs_combined = [all_traceIDs_combined; combined_data.traceIDs{group_idx}];
    end
end

% Combine all densities
all_densities_combined = [];
for group_idx = 1:nGroups
    if ~isempty(combined_data.densities{group_idx})
        all_densities_combined = [all_densities_combined; combined_data.densities{group_idx}];
    end
end

if ~isempty(all_areas_combined)
    fprintf(fid, '  Cluster Area Combined Median: %.1f nm²\n', median(all_areas_combined));
    fprintf(fid, '  Cluster Area Combined Mean: %.1f nm²\n', mean(all_areas_combined));
end
if ~isempty(all_nnd_combined)
    fprintf(fid, '  NND Combined Median: %.1f nm\n', median(all_nnd_combined));
    fprintf(fid, '  NND Combined Mean: %.1f nm\n', mean(all_nnd_combined));
end
if ~isempty(all_traceIDs_combined)
    fprintf(fid, '  Cluster Size Combined Median: %.1f trace IDs\n', median(all_traceIDs_combined));
    fprintf(fid, '  Cluster Size Combined Mean: %.1f trace IDs\n', mean(all_traceIDs_combined));
end
if ~isempty(all_densities_combined)
    fprintf(fid, '  Cluster Density Combined Median: %.4f trace IDs/nm²\n', median(all_densities_combined));
    fprintf(fid, '  Cluster Density Combined Mean: %.4f trace IDs/nm²\n', mean(all_densities_combined));
end
fprintf(fid, '\n');

fprintf(fid, '\n');

fprintf(fid, 'OUTPUT FILES:\n');
for group_idx = 1:nGroups
    fprintf(fid, '  %s group analysis: %s_Group_Analysis.png\n', group_names{group_idx}, group_names{group_idx});
end
fprintf(fid, '  Cell-level median statistics: CellLevel_Median_Statistics.csv\n');
fprintf(fid, '  Individual PDF subplots: %s\n', groupAnalysisFolder);
fprintf(fid, '  This statistical summary: %s\n', summary_file);

fclose(fid);
fprintf('Two-group statistical summary saved to: %s\n', summary_file);

%% 17. Display Final Summary
fprintf('\n==========================================\n');
fprintf('ANALYSIS COMPLETE FOR BOTH GROUPS!\n');
fprintf('==========================================\n');

for group_idx = 1:nGroups
    fprintf('\n%s Group Summary:\n', group_names{group_idx});
    fprintf('  Cells analyzed: %d\n', group_summary(group_idx).nCells);
    fprintf('  Total clusters: %d\n', group_summary(group_idx).total_clusters);
    fprintf('  Total points: %d\n', group_summary(group_idx).total_points);
    fprintf('  Total traces: %d\n', group_summary(group_idx).total_traces);
    fprintf('  Median cluster area: %.1f nm²\n', group_summary(group_idx).median_area);
    fprintf('  Median NND: %.1f nm\n', group_summary(group_idx).median_nnd);
    fprintf('  Median trace IDs per cluster: %.1f\n', group_summary(group_idx).median_traceIDs);
    fprintf('  Median cluster density: %.4f trace IDs/nm²\n', group_summary(group_idx).median_density);
end

fprintf('\nAnalysis Results:\n');
fprintf('  Group-level analysis folder: %s\n', groupAnalysisFolder);
for group_idx = 1:nGroups
    fprintf('  %s group analysis: %s_Group_Analysis.png\n', group_names{group_idx}, group_names{group_idx});
end
fprintf('  Cell-level median statistics: CellLevel_Median_Statistics.csv\n');
fprintf('  Individual PDF subplots saved in: %s\n', groupAnalysisFolder);
fprintf('  Statistical summary: %s\n', summary_file);

fprintf('\n==========================================\n');
fprintf('TWO-GROUP ANALYSIS COMPLETED SUCCESSFULLY!\n');
fprintf('==========================================\n');