clear all
close all

%% flat field correction with gaussian blurring
% read in 
raw_image = imread("/Users/andyzjc/Dropbox (Princeton)/Polarizatino Engineered Aberration Robust Adaptive Light Sheet Microscope/Figures/Figure 3/Fig3_3/Data/Panel A/Panel B/collagen_mip.png");
% raw_image = uint16(raw_image);
% raw_image = raw_image./max(raw_image,[],'all');

imagesc(raw_image)
axis image

%% gaussian blur
sigma = 100;
dark_field = imflatfield(raw_image,sigma);
% dark_field = dark_field./max(dark_field,[],'all');
figure; 
imagesc(uint16(dark_field))
colormap("hot")
axis image

imwrite(uint16(dark_field),'corrected.png')
% close all
