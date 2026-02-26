function [FinalImage, TotalFrame] = iQ_imgload3frames(dir_mv)
%% Bing Fu, 2024-05-11
% Modified from iQ_imgload3.m, corresponding to iQ_posfit3frames.m, works
% for unstacked individual SMT frames
%% Bing Fu, Oct. 4, 2017
%% Corresponding to iQ_rodcellfnd3, use TIFFStack to load movie fast
%% This function load the mv into matlab super-fast!!!
% Input: directory path of the folder contains the movie frames


%% Start of the routine
tic
cd(dir_mv); D = dir('*.tif');
TotalFrame=length(D);

% copy the tifflib to target dir
%dir_tiflib = 'C:\Program Files\MATLAB\R2012a\imagesci\private';
%cd(dir_tiflib)
cd(dir_mv);

% load the movie
FileTif=D(1).name;
InfoImage=imfinfo(FileTif)
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
FinalImage=zeros(nImage,mImage,TotalFrame,'uint16');
warning('off','all');
parfor j = 1:TotalFrame
    cd(dir_mv);
    FileTif=D(j).name;
    InfoImage=imfinfo(FileTif);
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages = numel(InfoImage);
    curr_frm = imread([dir_mv filesep FileTif]);
    FinalImage(:,:,j) = curr_frm(:,:,:);
end

toc
end