%% Ma, 2025-03-19
% Modify the code to handle both the directory input and add the discard frames option
% change: iQ_imgload3frames_file_input
% Modify the photon count method.
% checks if all frames have no detected spots and handles this case by returning empty matrices with the correct dimensions.
% checks if XYZ is empty after filtering out zero entries.
% adds a check for empty MRB and creates a default one if needed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ DFT,DFT_dr,DFT_drave,XYZ,PC,tEnd1,tEnd2] = iQ_posfit3frames_forDyeLocalization_file_input(dir_mv,threshold,mask_sz,IM_info,CAM_set,EM_wave_FP,MRB,dil_sz,refcoor,Z_cal_curve,showfit,params)
%% load the mv file
    tStart = tic;
    % Check if dir_mv is already an image array (for user-selected frames)
    if isnumeric(dir_mv) && ndims(dir_mv) == 3
        fprintf('Using pre-loaded image data with %d frames\n', size(dir_mv, 3));
        FinalImage = dir_mv;
        TotalFrame = size(FinalImage, 3);
    else
        % Load from file
        [FinalImage, TotalFrame] = iQ_imgload3frames_file_input(dir_mv);
    end
    % cd(setpath);
    tEnd1 = toc(tStart);

%% Initialize variables
    DFT = cell(TotalFrame,1);
    DFT_dr = cell(TotalFrame,1);
    DFT_drave = cell(TotalFrame,1);
    XYZ = cell(TotalFrame,1);
    BKinfo = cell(TotalFrame,1);

%% loop through all frames and do Gaussian fitting
    for i = 1:TotalFrame
        %% find peaks with background removal
        im = double(FinalImage(:,:,i));  % load curent frame
        % Check if ROI mask is available in params - no need to recreate if already exists
        if ~isfield(params, 'roiMask') || isempty(params.roiMask)
            % Only create ROI mask if MRB is provided and mask doesn't exist
            if exist('MRB', 'var') && ~isempty(MRB)
                fprintf('Creating ROI mask from MRB in frame %d...\n', i);
                % Initialize mask with zeros (no ROIs)
                roiMask = false(size(im));
                
                % Add each ROI to the mask
                for r = 1:size(MRB, 1)
                    if ~isempty(MRB{r,1})
                        % Create individual ROI mask
                        gridd = false(size(im));
                        gridd(sub2ind(size(gridd), MRB{r,1}(:,2), MRB{r,1}(:,1))) = 1;
                        
                        % Fill holes and dilate
                        grid1 = imfill(gridd, 'holes');
                        se1 = strel('disk', dil_sz);
                        grid2 = imdilate(grid1, se1);
                        
                        % Add to combined mask
                        roiMask = roiMask | grid2;
                    end
                end
                
                % Add to params for later use
                params.roiMask = roiMask;
                fprintf('ROI mask created and stored in params (covers %.2f%% of image)\n', ...
                    sum(roiMask(:))/numel(roiMask)*100);
            else
                % Create a default mask covering the entire image
                params.roiMask = true(size(im));
                fprintf('No MRB provided. Using entire image as ROI.\n');
            end
        end
        
        % Check if ROI mask is available
        if isfield(params, 'roiMask') && ~isempty(params.roiMask)
            % Apply ROI mask to restrict peak finding to ROI regions
            [im_bk, peak_fnd, ~, bkinfo] = iQ_pkfnd_V6_roi(im, threshold, mask_sz, params.roiMask);
        else
            % Use original peak finding without ROI restriction
            [im_bk,peak_fnd,~,bkinfo]=iQ_pkfnd_V6(im,threshold,mask_sz);   % find out the local maximum position in the current frame
        end
        
        BKinfo{i,1} = bkinfo;
        
        %% Visualize and save background removal and peak finding process (for first frame or every 100th frame)
        if i == 1 || mod(i, 1000) == 0
            % Create a directory for sample frame visualizations if it doesn't exist
            SampleFramesDir = fullfile(params.outputDir, '2_Localization', 'sample_frames');
            if ~exist(SampleFramesDir, 'dir')
                mkdir(SampleFramesDir);
                fprintf('Created directory for sample frames: %s\n', SampleFramesDir);
            end

            fig = figure('Name', sprintf('Frame %d - Background Removal and Peak Finding', i), ...
                        'Position', [100, 100, 1200, 500]);
            
            % Original image
            subplot(1, 3, 1);
            imagesc(im);
            colormap('gray');
            title('Original Image');
            colorbar;
            axis image;
            
            % Background-subtracted image
            subplot(1, 3, 2);
            imagesc(im_bk);
            colormap('gray');
            title({sprintf('Background-Subtracted (mask size: %d)', mask_sz), ...
                   sprintf('Mean: %.2f, Std: %.2f', bkinfo(1), bkinfo(2)), ...
                   sprintf('Threshold: %.2f', bkinfo(3))});
            colorbar;
            axis image;
            
            % Detected peaks overlaid on background-subtracted image
            subplot(1, 3, 3);
            imagesc(im_bk);
            colormap('gray');
            hold on;

            % Overlay ROI boundaries if available
            if isfield(params, 'roiMask') && ~isempty(params.roiMask)
                % Get the boundary of the ROI mask
                B = bwboundaries(params.roiMask);
                
                % Plot each boundary segment with more visible line
                for k = 1:length(B)
                    boundary = B{k};
                    plot(boundary(:,2), boundary(:,1), 'g-', 'LineWidth', 1);
                end
            end

            % Only visualize peaks that are inside the ROI
            if ~isempty(peak_fnd)
                % Calculate half box size for visualization
                half_box = floor(mask_sz/2);
                
                % Filter peaks to only show those inside ROI if ROI mask exists
                if isfield(params, 'roiMask') && ~isempty(params.roiMask)
                    % Convert peak coordinates to indices
                    validPeakIndices = [];
                    for p = 1:size(peak_fnd, 1)
                        % Check if peak is within image bounds
                        if peak_fnd(p,2) > 0 && peak_fnd(p,1) > 0 && ...
                           peak_fnd(p,2) <= size(im, 1) && peak_fnd(p,1) <= size(im, 2)
                            % Check if peak is inside ROI
                            if params.roiMask(round(peak_fnd(p,2)), round(peak_fnd(p,1)))
                                validPeakIndices = [validPeakIndices; p];
                            end
                        end
                    end
                    
                    % Only visualize peaks inside ROI
                    peak_fnd_vis = peak_fnd(validPeakIndices, :);
                    
                    % Add note about filtering
                    if size(peak_fnd_vis, 1) < size(peak_fnd, 1)
                        fprintf('Note: %d peaks were found but only %d are inside ROI and shown in visualization\n', ...
                            size(peak_fnd, 1), size(peak_fnd_vis, 1));
                    end
                else
                    % If no ROI mask, show all peaks
                    peak_fnd_vis = peak_fnd;
                end
                
                % Draw boxes around each detected peak
                for k = 1:size(peak_fnd_vis, 1)
                    x = peak_fnd_vis(k, 1);
                    y = peak_fnd_vis(k, 2);
                    
                    % Draw box with size based on mask_sz
                    rectangle('Position', [x-half_box, y-half_box, mask_sz, mask_sz], ...
                            'EdgeColor', 'r', 'LineWidth', 0.5, 'LineStyle', ':');
                end
                
                title({sprintf('%d Peaks Detected', size(peak_fnd, 1)), ...
                       sprintf('%d Peaks Inside ROI', size(peak_fnd_vis, 1))});
            else
                title('No Peaks Detected');
            end
            colorbar;
            axis image;
            zoom on;
            pan on;
            hold off;

            % Save figure as PNG with frame number in filename
            figFilename = fullfile(SampleFramesDir, sprintf('frame_%04d_visualization.png', i));
            saveas(fig, figFilename);
            
            % Optional: save as MATLAB .fig file for later editing
            figMatFilename = fullfile(SampleFramesDir, sprintf('frame_%04d_visualization.fig', i));
            saveas(fig, figMatFilename);     
            fprintf('Saved visualization for frame %d to: %s\n', i, figFilename);
        end

        %% Drift correction
        if refcoor == [0,0,0,0,0,0] % no drift correction excution
            DftSptCoor = [];
            DftSptDr = [];
            DftSptDrAve = [0 0 0 0 0];
        else % with drift correction excution
            if ~isempty(peak_fnd)
                DftSpt = cell2mat(iQ_nearestPoints(refcoor,peak_fnd));
            else
                if i == 1
                    DftSpt = round(refcoor);
                end
            end
            
            if ~isempty(Z_cal_curve) % 3D case
                DftSptCoor = iQ_Gaussian2DFit( DftSpt,im_bk,mask_sz,i,IM_info,Z_cal_curve )';
            else % 2D case
                DftSptCoor = iQ_Gaussian2DFit( DftSpt,im_bk,mask_sz,i,IM_info)';
            end
            if i == 1
                DftsptCoor1 = DftSptCoor;
            end
            
            DftSptDr = cellfun(@(x,y) x(1,1:5)-y(1,1:5),DftSptCoor,DftsptCoor1,'uni',false);
            DftSptDrAve = mean(cell2mat(DftSptDr'),1);
            % if the drift is larger than 3 pixel, probably the marker is not
            % good, the fitting give incorrect values, if this happens, the
            % drift correction will be deactived for this particular frame.
            if DftSptDrAve(1) > 3||DftSptDrAve(2) > 3
                DftSptDr = DFT_dr{i-1,1};
                DftSptDrAve = DFT_drave{i-1,1};
            end
        end
        DFT{i,1} = DftSptCoor;
        DFT_dr{i,1} = DftSptDr;
        DFT_drave{i,1} = DftSptDrAve;
        % get the spot fitting results
        % ----------------------------------------------------------------------
        if ~isempty(peak_fnd)
            % since MRG is determined by the first frame, any peaks found
            % in the following frames need to drift corrected first then be judged
            % the peak is keeped/discarded(inside/outside the mask).
            DFTxy = round(DftSptDrAve(ones(size(peak_fnd,1),1),:)); % this makes drift vector a matrix that has the same dimension as peak_fnd.
            DFTxy = DFTxy(:,1:2);
            peak_fnddc = peak_fnd-DFTxy; % drift correct the peak first
    
            % IntSpt = intersect(peak_fnddc,MRG,'rows'); % pick up only those peaks are inside the cell.
            DFTxy = DFTxy(1:size(peak_fnddc,1),:);
            peak_fnddc = peak_fnddc+DFTxy; % get back to the original peak position to start the fitting
            clear DFTxy
        else
            peak_fnddc=[];
        end
        %% Gaussian fit
        if isempty(peak_fnddc)
            XYZ{i,1}=[0 0 0 0 0 0 0 0 0 0 0 i 0 0 0 0 ];
            % Fix the escaped character issue by using disp() instead of display()
            % and properly formatting the string without escape sequences
            disp(['Frame no. ', int2str(i), '   Spots fitted:  nothing inside the cell']);
        else
            if ~isempty(Z_cal_curve);
                IntSptCoor = iQ_Gaussian2DFit( peak_fnddc,im_bk,mask_sz,i,IM_info,Z_cal_curve );   % 3D Gaussian fit all spots inside the cell
            else
                IntSptCoor = iQ_Gaussian2DFit( peak_fnddc,im_bk,mask_sz,i,IM_info);
            end
            IntSptCoor = cell2mat(IntSptCoor);
            IntSptCoor(:,1:5) = IntSptCoor(:,1:5)-DftSptDrAve(ones(size(IntSptCoor(:,1:5),1),1),:);
            XYZ{i,1} = IntSptCoor;
            % Fix the escaped character issue by using disp() instead of display()
            % and properly formatting the string without escape sequences
            disp(['Frame no. ', int2str(i), '   Spots fitted: ', int2str(size(peak_fnddc,1)), '  spots']);
        end
    end

%% Handle the case where no spots were detected in any frame
    % Check if all XYZ cells are empty or contain only placeholder values
    allEmpty = true;
    for i = 1:length(XYZ)
        if ~isempty(XYZ{i,1}) && ~(size(XYZ{i,1},1) == 1 && XYZ{i,1}(1,1) == 0)
            allEmpty = false;
            break;
        end
    end
    
    
    if allEmpty
        fprintf('Warning: No spots detected in any of the analyzed frames.\n');
        % Return empty matrices with correct dimensions
        XYZ = zeros(0, 16);  % Empty matrix with 16 columns
        DFT_dr = zeros(0, 5);
        DFT_drave = zeros(0, 5);
        PC = {zeros(0, 22)};  % Empty cell with one element containing empty matrix with 22 columns
        return;
    end
    
%% Continue with normal processing if spots were detected, update variables
    XYZ = cell2mat(XYZ);
    DFT_dr = cell2mat(vertcat(DFT_dr{:}));
    DFT_drave = cell2mat(DFT_drave);
    XYZ(XYZ(:,1)==0,:)=[];
    BKinfo = cell2mat(BKinfo);
    
%% Group the spots according cell mask and estimate num of photons, errxy
    %======================================================================
    % calculate the Errx and Erry for each spot
    %----------------------------------------------------------------------
    
    % h = 6.626e-34;    % Planck Constant (m^2 kg / s)
    % c = 299792458;    % Speed of light (m / s)
    % ev = 1.602e-19;    % Energy of 1eV (kg m^2 / s^2)
    % E_ph= h*(c/(EM_wave_FP*1e-9)); % calculate the photon energy of emission light
    % g = CAM_set(1); S = CAM_set(2); QE = CAM_set(3);
    camera_params.gain = CAM_set(1);          % electrons/ADU
    camera_params.QE = CAM_set(2);          % Quantum efficiency
    
    Pixel_x = IM_info(1); Pixel_y = IM_info(2);
    
    % cts_n = XYZ(:,15);
    % N_e = (cts_n/g)*(S/QE)*3.65;
    % N_p = N_e*ev/E_ph;
    cts_n = XYZ(:,15);
    N_p = intensity2photon(cts_n,camera_params);
    
    % cts_b = XYZ(:,16);
    % b_e = (cts_b/g)*(S/QE)*3.65;
    % b_p = b_e*ev/E_ph;
    % Get the global background statistics for each spot's frame
    bkgmean = BKinfo(XYZ(:,12), 1);  % Global background mean
    bkgstd = BKinfo(XYZ(:,12), 2);   % Global background std ← THIS IS WHAT WE NEED
    bkgth = BKinfo(XYZ(:,12), 3);    % Threshold

     % add Err_x, Err_y as placeholders first (will be overwritten)
    XYZ = [XYZ(:,1:16), ...              % Original 16 columns
                zeros(size(XYZ,1), 1), ...    % Column 17: Err_x (placeholder)
                zeros(size(XYZ,1), 1), ...    # Column 18: Err_y (placeholder)
                N_p, ...                      % Column 19: N_p (photon count)
                bkgmean, bkgstd, bkgth];      % Columns 20-22: BKinfo

    % measure bkg by local bkgstd
    cts_b = XYZ(:,21);            
    cts_b = max(cts_b, 0);
    b_p = intensity2photon(cts_b, camera_params);

    
    % Calculate Thompson localization precision using the formula
    sgm_x = XYZ(:,8);  % XYZ(:,8) gives the sigmax in nm
    Err_x = ((sgm_x.^2)./N_p+...
        Pixel_x^2./(12*N_p)+...
        (8*pi*(sgm_x.^4).*b_p.^2)./(Pixel_x*N_p).^2).^0.5;
    
    sgm_y = XYZ(:,9);  % XYZ(:,9) gives the sigmay in nm
    Err_y = ((sgm_y.^2)./N_p+...
        Pixel_y^2./(12*N_p)+...
        (8*pi*(sgm_y.^4).*b_p.^2)./(Pixel_y*N_p).^2).^0.5;
    
    % PC_new=[x(px) y(px) x(nm) y(nm) z(nm)
    %        sigma_x(px) sigma_y(px) sigma_x(nm) sigma_y(nm) ellip PSFW_H(px)
    %        frame# amp(cts) background(cts) Int(cts) std_b(cts)  Err_x(nm) Err_y(nm) Intint(cts)]
    bkgmean = BKinfo(XYZ(:,12),1); bkgstd = BKinfo(XYZ(:,12),2);
    bkgth = BKinfo(XYZ(:,12),3);
    XYZ = [XYZ(:,1:16) Err_x Err_y,N_p,bkgmean,bkgstd,bkgth];% XYZ(:,1:17)];
    
    %======================================================================
    % find out spots inside the region of interest (ROI)
    %----------------------------------------------------------------------

    % Check if XYZ is empty after filtering
    if isempty(XYZ)
        fprintf('Warning: No valid spots remain after filtering.\n');
        % Return empty matrices with correct dimensions
        PC = {zeros(0, 22)};  % Empty cell with one element containing empty matrix with 22 columns
        return;
    end
    % Initialize ROI_grid
    ROI_grid = cell(size(MRB, 1), 1);

    % find out and index the ROI and create the grid for each ROIs
    for j = 1:size(MRB,1)
        gridd = false(size(im));
        if ~isempty(MRB{j,1})
            gridd(sub2ind(size(gridd),MRB{j,1}(:,2),MRB{j,1}(:,1))) = 1;
            
            grid1 = imfill(gridd, 'holes');
            se1 = strel('disk',dil_sz);
            grid2 = imdilate(grid1,se1);
            
            [row_cc,col_cc] = find(grid2>0);
            ROI_grid{j,1} = [col_cc,row_cc];
        else
            ROI_grid{j,1} = [];
        end
    end

    % relocate spots into corresponding ROIs.
    XYZ_round = round(XYZ(:,1:2));
    PC = cell(length(ROI_grid),1);
    for m = 1:length(ROI_grid)  % loop through all cells
        if ~isempty(ROI_grid{m,1})
        spt_member = ismember(XYZ_round,ROI_grid{m,1},'rows');
        XYZ_ROI = XYZ(spt_member,:);
        PC{m,1} = XYZ_ROI;      % protein coordinates inside each cell
        else
            PC{m,1} = zeros(0, size(XYZ, 2));  % Empty matrix with same number of columns as XYZ
        end
    end
    
    
%% report elapsed time
    tEnd2 = toc(tStart);
    fprintf('Total mv loading time: %d minutes and %f seconds\n',floor(tEnd1/60),rem(tEnd1,60));
    fprintf('Total spot fitting time: %d minutes and %f seconds\n',floor(tEnd2/60),rem(tEnd2,60));
end
    
    