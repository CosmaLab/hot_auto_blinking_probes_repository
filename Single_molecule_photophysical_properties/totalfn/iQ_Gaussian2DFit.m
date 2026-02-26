 function [ GaussianFitResult ] = iQ_Gaussian2DFit( pick,im_bk,mask_sz,n,IM_info,Z_cal_curve )
%UNTITLED2 Summary of this function goes here
% Parameters:
Pixel_x = IM_info(1); 
corr_cylin = IM_info(1)/IM_info(2);
sz = (mask_sz-1)/2;
num_spt = size(pick,1);        
GaussianFitResult = cell(num_spt,1);

for j = 1:num_spt 
    
   [xi,yi] = meshgrid(pick(j,1)-sz:pick(j,1)+sz,pick(j,2)-sz:pick(j,2)+sz);
   zi = im_bk(pick(j,2)-sz:pick(j,2)+sz,pick(j,1)-sz:pick(j,1)+sz);
   
   results = iQ_autoGaussianSurfML(xi,yi,zi);
   
   % fitting results
   amp = results.a; 
   background = results.b; 
   sigma_x = results.sigmax; 
   sigma_y = results.sigmay/corr_cylin;
   mu_x = results.x0; 
   mu_y = results.y0/corr_cylin; 
   Int = sum(sum(results.G-results.b));
   residue = zi-results.G; 
   std_b = mean(std(residue(:)));
   ellip = 2*abs((sigma_x-sigma_y))/(sigma_x+sigma_y);
   PSFW_H = sigma_x-sigma_y;     % calculate the width-height vaue from PSF fitting
   
   mu_z = 0;
   
   if nargin == 6;
        x = Z_cal_curve(:,1);                                     
        y = Z_cal_curve(:,2);
        mu_z = interp1(y,x,PSFW_H); % calulate z from the width-height of psf
   end
   
   GaussianFitResult{j,1} = [...
       mu_x,     mu_y,     mu_x*Pixel_x,     mu_y*Pixel_x,     mu_z...
       sigma_x,  sigma_y,  sigma_x*Pixel_x,  sigma_y*Pixel_x,  ellip...
       PSFW_H,   n,        amp,              background,       Int,     std_b ];

end
end