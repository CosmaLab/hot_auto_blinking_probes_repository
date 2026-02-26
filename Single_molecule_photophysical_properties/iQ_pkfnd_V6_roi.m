%% Ma, 2025-03-19
% iQ_pkfnd_V6_roi - Modify the code to Find local maxima in an image restricted to ROI regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [im_bk,peak,intensity,bk] = iQ_pkfnd_V6_roi(im, threshold, sz, roiMask)

    if nargin < 3, sz = 11; end
    if nargin < 2, threshold = 4; end
    if nargin < 4 || isempty(roiMask)
        % If no ROI mask provided, create a default mask covering the entire image
        roiMask = true(size(im));
        fprintf('No ROI mask provided. Using entire image.\n');
    end
    
    % Ensure ROI mask is binary and same size as image
    if ~islogical(roiMask)
        roiMask = logical(roiMask);
    end
    
    if ~isequal(size(roiMask), size(im))
        error('ROI mask dimensions must match image dimensions');
    end

    % Print ROI coverage information
    roiCoverage = sum(roiMask(:)) / numel(roiMask) * 100;
    fprintf('ROI covers %.2f%% of the image area\n', roiCoverage);

    im = double(im);   
    [nr,nc] = size(im);       % find out the size of image
    
    %% STEP 1: BACKGROUND REMOVAL(subtracting the locally estimated background, low-pass filtering)
    bk_bxr = ones(20,20)/400; %original (15,15)/225, modified:(20,20)/400
    bk_fit = imfilter(im,bk_bxr,'replicate');
    
    im_bk = im-bk_fit;% remove the background from original image
    
    %% STEP 2: LOCAL MAXIMA DETECTION
    mask = ones(sz);                                      % define the mask for one PSF spot   
    mask(ceil(sz/2),ceil(sz/2))=0;   
    % A pixel is a local maximum if its value is greater than the maximum value in its neighborhood
    bw = im_bk>imdilate(im_bk,mask);% imdilate finds max value in neighborhood defined by mask
    
    % Apply ROI mask to restrict peaks to ROI regions - CRITICAL STEP
    bw = bw & roiMask;
    
    % Count potential peaks before and after ROI masking for debugging
    totalPeaks = sum(im_bk>imdilate(im_bk,mask), 'all');
    roiPeaks = sum(bw, 'all');
    fprintf('Found %d potential peaks, %d inside ROI (%.1f%%)\n',...
            totalPeaks, roiPeaks, (roiPeaks/totalPeaks*100));
    
    %% STEP 3: EXTRACT COORDINATES OF LOCAL MAXIMA
    [r,c] = find(bw>0);  rc=[r,c];                        % get the coordinates of local maximum spots
    rc((rc(:,1)<2) | (rc(:,2)<2) | (rc(:,1)>nr-1)|(rc(:,2)>nc-1),:) = [];
    
    % Double-check that all peaks are within ROI (redundant but safe)
    if ~isempty(rc)
        validIndices = false(size(rc,1),1);
        for i = 1:size(rc,1)
            validIndices(i) = roiMask(rc(i,1), rc(i,2));
        end
        rc = rc(validIndices,:);
    end
    
    %% STEP 4: CALCULATE BACKGROUND STATISTICS
    % use bottom 95% spots to estimate the mean background value and the
    % standard deviation ( this is useful when your ROI is small, basically
    % eliminates all the bright spots and only consider the real background.
    int_im_bk = im_bk(:);                                
    vs = sort(int_im_bk,'ascend');
    n95 = round(numel(int_im_bk)*95/100)+1;
    v95 = vs(1:n95);
    
    im_bk_mean = mean(v95);                       % find the average value of background removed image
    im_bk_std = std(v95);                         % find the standard deviation of background removed image
    th = im_bk_mean+threshold*im_bk_std;                  % set up the threshold for peak determination
    bk = [im_bk_mean,im_bk_std,th];                       % generate the bk info: [im_bk_mean,im_bk_std,th]=[average of image,standard deviation of image,threshold for peak finding
    
    %% STEP 5: FILTER PEAKS BASED ON INTENSITY THRESHOLD
    if isempty(rc)
        fprintf('No local maxima found within ROI\n');
        peak = []; intensity = [];
        return;
    else
        rc_th = [];
        r1 = 1;                                               % obtain the coordinates which shows intensity larger than the threshold
        [rc_row,~] = size(rc);
    for q = 1:rc_row;
        % Extract 3x3 neighborhood around the peak
        pkint = im_bk(rc(q,1)-1:rc(q,1)+1,rc(q,2)-1:rc(q,2)+1);
        pkint = sort(pkint(:),'descend'); % Sort neighborhood values in descending order
        pkint = mean(pkint(1:3)); % Average the top 3 values (robust peak intensity)
        % Keep peak if its intensity exceeds the threshold
        if pkint > th;
            rc_th(r1,:) = rc(q,:);
            r1 = r1+1;
        end
    end
    
    %% STEP 6: HANDLE CASE WHERE NO PEAKS EXCEED THRESHOLD
    if isempty(rc_th)
        fprintf('No peaks above threshold within ROI\n');
        peak = []; intensity = [];
        return;
    else
        %% STEP 7: SUPPLEMENT WITH ADDITIONAL PEAKS FROM GUESSING ALGORITHM
        rc_th1 = rc_th;
        gues = Guessing_v2(im,size(im),7,97,7,false,false,false); 
        rc_th2 = gues(:,2:3);
        
        % Filter guessed peaks to only include those within the ROI - CRITICAL STEP
        if ~isempty(rc_th2)
            validGuessIndices = [];
            for i = 1:size(rc_th2, 1)
                % Check if the guessed peak is within the ROI
                % Note: rc_th2 is [x,y] but roiMask is indexed as [row,col]
                if rc_th2(i,2) > 0 && rc_th2(i,1) > 0 && ...
                   rc_th2(i,2) <= nr && rc_th2(i,1) <= nc && ...
                   roiMask(rc_th2(i,2), rc_th2(i,1))
                    validGuessIndices = [validGuessIndices; i];
                end
            end
            rc_th2 = rc_th2(validGuessIndices, :);
            fprintf('Guessing algorithm found %d additional peaks within ROI\n', size(rc_th2, 1));
        end
        
        rc_th = rc_th1; % Start with original peaks
    
        % Add peaks from guessing algorithm if they're not too close to existing peaks
        for ii = 1:size(rc_th2,1)
            cur = rc_th2(ii,:);
            dists = sqrt((rc_th1(:,1)-cur(1,1)).^2+(rc_th1(:,2)-cur(1,2)).^2);
            if find(dists <= 2) % close than 2 px
                % Skip this peak as it's too close to an existing one
            else
                rc_th = [rc_th;cur];
            end
        end
    
        peak(:,1)=rc_th(:,2);                               % generate the peak output: [column1,column2]= [x,y]
        peak(:,2)=rc_th(:,1);
    
    %% STEP 9: REMOVE PEAKS TOO CLOSE TO IMAGE BOUNDARIES, Keep only peaks that are at least half_sz away from any boundary
        ind_boundary= peak(:,1)>ceil(sz/2) & peak(:,1)<nc-ceil(sz/2) & peak(:,2)>ceil(sz/2) & peak(:,2)<nr-ceil(sz/2);
        peak=peak(ind_boundary,:); 
    
    %% STEP 10: FINAL ROI CHECK - Ensure all peaks are within ROI
        if ~isempty(peak)
            finalValidIndices = false(size(peak,1),1);
            for i = 1:size(peak,1)
                % Note: peak is [x,y] but roiMask is indexed as [row,col]
                if peak(i,2) > 0 && peak(i,1) > 0 && ...
                   peak(i,2) <= nr && peak(i,1) <= nc && ...
                   roiMask(peak(i,2), peak(i,1))
                    finalValidIndices(i) = true;
                end
            end
            peak = peak(finalValidIndices,:);
        end
    
    %% STEP 11: Calculate sum of intensities in sz×sz window around each peak
        if isempty(peak)
            fprintf('No valid peaks remain after all filtering steps\n');
            intensity = [];
            return;
        end
        
        ind_th=sub2ind(size(im_bk),peak(:,2),peak(:,1));    % convert new coordinates to index
        im_bk_sums=imfilter(im_bk,ones(sz));                % calculate the intergrated intensity over specified area
        im_bk_int=im_bk_sums(ind_th);
        intensity=im_bk_int;       
        
        % Print summary of results
        fprintf('Found %d peaks within ROI after all filtering steps\n', size(peak, 1));         
    end
    end
end