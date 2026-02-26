function [ DFT,DFT_dr,DFT_drave,XYZ,PC,tEnd1,tEnd2] = iQ_posfit3frames_forDyeLocalization( dir_mv,threshold,mask_sz,IM_info,CAM_set,EM_wave_FP,MRB,dil_sz,refcoor,Z_cal_curve,showfit)
%% load the mv file
tStart = tic;
[FinalImage, TotalFrame ]= iQ_imgload3frames(dir_mv);
% cd(setpath);
tEnd1 = toc(tStart);
%%%%%%%%%%
%%  loop through the mv and do 2D Gaussian fittings for spots
DFT = cell(TotalFrame,1);
DFT_dr = cell(TotalFrame,1);
DFT_drave = cell(TotalFrame,1);
XYZ = cell(TotalFrame,1);
BKinfo = cell(TotalFrame,1);

% for i = 1:TotalFrame
for i = 1:5
    im = double(FinalImage(:,:,i));  % load first frame
    [im_bk,peak_fnd,~,bkinfo]=iQ_pkfnd_V6(im,threshold,mask_sz);   % find out the local maximum position in the current frame
    BKinfo{i,1} = bkinfo;
    
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
    if isempty(peak_fnddc)
        XYZ{i,1}=[0 0 0 0 0 0 0 0 0 0 0 i 0 0 0 0 ];
        display(['Frame no. ',int2str(i),'   Sptos fitted: ',' nothing inside the cell'])
    else
        if ~isempty(Z_cal_curve);
            IntSptCoor = iQ_Gaussian2DFit( peak_fnddc,im_bk,mask_sz,i,IM_info,Z_cal_curve );   % 3D Gaussian fit all spots inside the cell
        else
            IntSptCoor = iQ_Gaussian2DFit( peak_fnddc,im_bk,mask_sz,i,IM_info);
        end
        IntSptCoor = cell2mat(IntSptCoor);
        IntSptCoor(:,1:5) = IntSptCoor(:,1:5)-DftSptDrAve(ones(size(IntSptCoor(:,1:5),1),1),:);
        XYZ{i,1} = IntSptCoor;
        display(['Frame no. ',int2str(i),'   Sptos fitted: ',int2str(size(peak_fnddc,1)), '  spots'])
    end
    %     ======================================================================
    % indicate the fitted spots with white box
    % ----------------------------------------------------------------------
    if showfit == 1
        for ii = 1:size(XYZ{i,1},1)
            curr_x = round(XYZ{i,1}(ii,1));
            curr_y = round(XYZ{i,1}(ii,2));
            if curr_x <= 4
                curr_x = 5;
            elseif curr_x > size(im,2)-4
                curr_x = size(im,2);
            end
       
            if curr_y <= 4
                curr_y = 5;
            elseif curr_y > size(im,1)-4
                curr_y = size(im,1);
            end
            im(curr_y-4,curr_x-4:curr_x+4) = max(max(im));
            im(curr_y+4,curr_x-4:curr_x+4) = max(max(im));
            im(curr_y-4:curr_y+4,curr_x-4) = max(max(im));
            im(curr_y-4:curr_y+4,curr_x+4) = max(max(im));
        end
        imshow(im,[]);
    end
end

BKinfo = cell2mat(BKinfo);
XYZ = cell2mat(XYZ);
DFT_dr = cell2mat(vertcat(DFT_dr{:}));
DFT_drave = cell2mat(DFT_drave);
XYZ(XYZ(:,1)==0,:)=[];

%% Group the spots according cell mask and estimate num of photons, errxy
%======================================================================
% calculate the Errx and Erry for each spot
%----------------------------------------------------------------------

h = 6.626e-34;    % Planck Constant (m^2 kg / s)
c = 299792458;    % Speed of light (m / s)
ev = 1.602e-19;    % Energy of 1eV (kg m^2 / s^2)
E_ph= h*(c/(EM_wave_FP*1e-9)); % calculate the photon energy of emission light
g = CAM_set(1); S = CAM_set(2); QE = CAM_set(3);

Pixel_x = IM_info(1); Pixel_y = IM_info(2);

cts_n = XYZ(:,15);
N_e = (cts_n/g)*(S/QE)*3.65;
N_p = N_e*ev/E_ph;

cts_b = XYZ(:,16);
b_e = (cts_b/g)*(S/QE)*3.65;
b_p = b_e*ev/E_ph;

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

% find out and index the ROI and create the grid for each ROIs

ROI_grid = cell(size(MRB,1),1);
for j = 1:size(MRB,1);
    grid = false(size(im));
    grid(sub2ind(size(im),MRB{j,1}(:,2),MRB{j,1}(:,1))) = 1;
    
    grid1 = imfill(grid, 'holes');
    se1 = strel('disk',dil_sz);
    grid2 = imdilate(grid1,se1);
    
    [row_cc,col_cc] = find(grid2>0);
    ROI_grid{j,1} = [col_cc,row_cc];
end
% relocate spots into corresponding ROIs.
XYZ_round = round(XYZ(:,1:2));
PC = cell(length(ROI_grid),1);
for m = 1:length(ROI_grid)  % loop through all cells
    spt_member = ismember(XYZ_round,ROI_grid{m,1},'rows');
    XYZ_ROI = XYZ(spt_member,:);
    PC{m,1} = XYZ_ROI;      % protein coordinates inside each cell
end


%% report elapsed time
tEnd2 = toc(tStart);
fprintf('Total mv loading time: %d minutes and %f seconds\n',floor(tEnd1/60),rem(tEnd1,60));
fprintf('Total spot fitting time: %d minutes and %f seconds\n',floor(tEnd2/60),rem(tEnd2,60));
end

