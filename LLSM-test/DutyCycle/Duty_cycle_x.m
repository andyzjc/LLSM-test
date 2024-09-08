%% test 
clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

clc
NA1 = 0.5; %0.6 for gaussian
deltaNA = 0.1;
LatticeType = 'hex';
ProfileType = 'tophat';
SWweighting = 1; %4/3 for equal OTF V2 LLS, 1/sqrt(2) for V1 LLS
Latticeweighting = [1,1,1,1,1,1]; % 1.9 for V2 LLS
center = (N+1)/2;

[LatticePupil,LatticeMask,LatticeMetaData] = GetLatticePupil(LatticeType,ProfileType, ...
                                                             NA1,deltaNA, ...
                                                             0.8,0.0,...
                                                             Latticeweighting);

[SWPupil,SWMask,SWPupilMetaData] = GetSWPairPupil(ProfileType,NA1,NA1/2,...
                                                           deltaNA,deltaNA*2,...
                                                           SWweighting);

% simulate, and find dither range for LLS, determine one period to compare
xzLatticePSF = abs( fftshift( ifft2(ifftshift(LatticePupil)) ) ).^2;
xzLLSmax = max(xzLatticePSF,[],'all')
xzLLSmean = mean(xzLatticePSF(:))

xzSW1PSF = abs( fftshift( ifft2(ifftshift(SWPupil(:,:,1))) ) ).^2;
xzSW2PSF = abs( fftshift( ifft2(ifftshift(SWPupil(:,:,2))) ) ).^2;
xzPEARLSPSF = xzSW1PSF + xzSW2PSF;
xzPEARLSmax = max(xzPEARLSPSF,[],'all')
xzPEARLSmean = mean(xzPEARLSPSF(:))

[~,LOCS] = findpeaks(xzLatticePSF(center,center:end),'MinPeakHeight',0.8 .* max(xzLatticePSF,[],'all'));
dither_range = LOCS(1,1)-1 %pixel 
area_of_interest = center: (center+dither_range);

% xzPEARLSPSF = xzPEARLSPSF/max(xzPEARLSPSF,[],'all');
% xzLatticePSF = xzLatticePSF./max(xzLatticePSF,[],'all');

clc
xLLSPSF = xzLatticePSF(center,area_of_interest);
xPEARLSPSF = xzPEARLSPSF(center,area_of_interest);

% same peak, compare mean
xPEARLSPSF_peak_scale = xzLLSmax./xzPEARLSmax .* xPEARLSPSF; % scale to same peak 
PEARLSmean_same_peak = mean(xPEARLSPSF_peak_scale);
PEARLSpeak_same_peak = max(PEARLSmean_same_peak);
LLSmean = mean(xLLSPSF);
LLSpeak = max(xLLSPSF);
same_peak_mean_ratio = LLSmean/PEARLSmean_same_peak

% same mean, compare peak 
xPEARLSmean = mean(xPEARLSPSF);
xPEALRSpeak = max(xPEARLSPSF);
% xPEALRS_ratio = xPEALRSpeak./xPEARLSmean;

xLLSmean = mean(xLLSPSF);
xLLSpeak = max(xLLSPSF);
% xLLSpeak_ratio = xLLSpeak./xLLSmean;

xPEARLSPSF_mean_scale = xLLSmean./xPEARLSmean .* xPEARLSPSF; % scale to same mean 
same_mean_LLS_to_PEARLS_peak_ratio = xLLSpeak./max(xPEARLSPSF_mean_scale)


%% Create heat map, same mean, compare peak, 

clear all
close all
clc

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

% same peak, compare mean
NA_min = 0.3;
NA_max = 0.6;
NA_min_k_pixel = round(NA_min/(deltak .* n / k_wave));  % just to make sure everything is align to pixel
NA_max_k_pixel = round(NA_max/(deltak .* n / k_wave));
NA1_pixel = NA_min_k_pixel: 1 : NA_max_k_pixel;
NA1 = NA1_pixel*deltak .* n /k_wave;
deltaNA1_pixel = 1:1:round(0.2/(deltak .* n / k_wave)); 
deltaNA1 = deltaNA1_pixel .*deltak .* n /k_wave;
center = (N+1)/2;

LatticeType = 'hex'; 
ProfileType = 'tophat';
SWweighting = 1; %4/3 for equal OTF V2 LLS, 1/sqrt(2) for V1 LLS
Latticeweighting = [1,1,1,1,1,1]; % 1.9 for V2 LLS

% ratio maps
PEALRS_peak_mean_ratio_map = zeros(length(NA1), length(deltaNA1));
LLS_peak_mean_ratio_map = PEALRS_peak_mean_ratio_map;
same_mean_LLS_to_PEARLS_peak_ratio_map = PEALRS_peak_mean_ratio_map;
% same_peak_LLS_to_PEARLS_mean_ratio_map = PEALRS_peak_mean_ratio_map;
just_overall_LLS_to_PEARLS_peak_ratio_map = PEALRS_peak_mean_ratio_map;
just_overall_LLS_to_PEARLS_mean_ratio_map = PEALRS_peak_mean_ratio_map;

PEALRS_peak_mean_ratio_array = [];
LLS_peak_mean_ratio_array = [];
same_mean_LLS_to_PEARLS_peak_ratio_array = [];
% same_peak_LLS_to_PEARLS_mean_ratio_array = [];
just_overall_LLS_to_PEARLS_peak_ratio_array = [];
just_overall_LLS_to_PEARLS_mean_ratio_array = [];

counter = 1;
for i = 1:length(NA1)
    for j = 1:length(deltaNA1)

        % determine if PEARLS can have equivalent LLS 

        %%%%%%%%%%%%%%%%%
        % deltaNA2 = sqrt(2*NA1(i)*deltaNA1(j)); NA2 = deltaNA2./2; % square 
        deltaNA2 = deltaNA1(j)*2; NA2 = NA1(i) ./ 2; % hex
        %%%%%%%%%%%%%%%%%

        NA1min = NA1(i) - deltaNA1(j)/2;
        NA1max = NA1(i) + deltaNA1(j)/2;
        NA2min = NA2 - deltaNA2/2
        NA2max = NA2 + deltaNA2/2

        if NA2max < NA1min
            NA1Array(counter,1) = NA1(i);
            deltaNA1Array(counter,1) = deltaNA1(j);

            % generate Pupil 
            [LatticePupil,LatticeMask,LatticeMetaData] = GetLatticePupil(LatticeType,ProfileType, ...
                                                 NA1(i),deltaNA1(j), ...
                                                 0.8,0.0,...
                                                 Latticeweighting);

            if isequal(LatticeType,'hex')
                [SWPupil,~,~] = GetSWPairPupil(ProfileType,NA1(i),NA2,...
                                                        deltaNA1(j),deltaNA2,...
                                                        SWweighting);
                NA2Array(counter,1) = NA2;
                deltaNA2Array(counter,1) = deltaNA2;
            else
                [SWPupil,~,~] = GetSWPairPupil(ProfileType,NA1(i),NA2,...
                                                        deltaNA1(j),deltaNA2,...
                                                        SWweighting);
                NA2Array(counter,1) = deltaNA2/2;
                deltaNA2Array(counter,1) = deltaNA2;
            end

            % simulate, and find dither range for LLS, determine one period to compare
            xzLatticePSF = abs( fftshift( ifft2(ifftshift(LatticePupil)) ) ).^2;
            xzLLSmax = max(xzLatticePSF,[],'all');
            xzLLSmean = mean(xzLatticePSF(:));
            
            xzSW1PSF = abs( fftshift( ifft2(ifftshift(SWPupil(:,:,1))) ) ).^2;
            xzSW2PSF = abs( fftshift( ifft2(ifftshift(SWPupil(:,:,2))) ) ).^2;
            xzPEARLSPSF = xzSW1PSF + xzSW2PSF;
            xzPEARLSmax = max(xzPEARLSPSF,[],'all');
            xzPEARLSmean = mean(xzPEARLSPSF(:));

            just_overall_LLS_to_PEARLS_peak_ratio_map(i,j) = xzLLSmax/xzPEARLSmax; % just curious
            just_overall_LLS_to_PEARLS_mean_ratio_map(i,j) = xzLLSmean/xzPEARLSmean; % just curious
            
            [~,LOCS] = findpeaks(xzLatticePSF(center,center:end),'MinPeakHeight',0.8 .* max(xzLatticePSF,[],'all'));
            dither_range = LOCS(1,1)-1; %pixel 
            area_of_interest = center: (center+dither_range);
            
            xLLSPSF = xzLatticePSF(center,area_of_interest);
            xPEARLSPSF = xzPEARLSPSF(center,area_of_interest);
            
            % same mean, compare peak 
            xPEARLSmean = mean(xPEARLSPSF);
            xPEALRSpeak = max(xPEARLSPSF);
            PEALRS_peak_mean_ratio_map(i,j) = xPEALRSpeak./xPEARLSmean; % this is what the task wants, should always be 1
            
            xLLSmean = mean(xLLSPSF);
            xLLSpeak = max(xLLSPSF);
            LLS_peak_mean_ratio_map(i,j) = xLLSpeak./xLLSmean; % this is what the task wants
            
            xPEARLSPSF_mean_scale = xLLSmean./xPEARLSmean .* xPEARLSPSF; % scale to same mean 
            same_mean_LLS_to_PEARLS_peak_ratio_map(i,j) = xLLSpeak./max(xPEARLSPSF_mean_scale); % this is what I understand, should be same as above

            % same peak, compare mean (just for fun) 
            % xPEARLSPSF_peak_scale = xzLLSmax./xzPEARLSmax .* xPEARLSPSF; % scale to same peak 
            % PEARLSmean_same_peak = mean(xPEARLSPSF_peak_scale);
            % PEARLSpeak_same_peak = max(PEARLSmean_same_peak);
            % same_peak_LLS_to_PEARLS_mean_ratio_map(i,j) = xLLSmean/PEARLSmean_same_peak; % I thought this should be closest to experiment, but result is anti-intuitive

            % for ploting purpose
            PEALRS_peak_mean_ratio_array(counter,1) = PEALRS_peak_mean_ratio_map(i,j);
            LLS_peak_mean_ratio_array(counter,1) = LLS_peak_mean_ratio_map(i,j);
            same_mean_LLS_to_PEARLS_peak_ratio_array(counter,1) = same_mean_LLS_to_PEARLS_peak_ratio_map(i,j);
            % same_peak_LLS_to_PEARLS_mean_ratio_array(counter,1) = xLLSmean/PEARLSmean_same_peak;
            just_overall_LLS_to_PEARLS_peak_ratio_array(counter,1) = just_overall_LLS_to_PEARLS_peak_ratio_map(i,j);
            just_overall_LLS_to_PEARLS_mean_ratio_array(counter,1) = just_overall_LLS_to_PEARLS_mean_ratio_map(i,j);

            counter = counter + 1;
        end
    end
end

% example heat map showing peak intensity ratio
figure;imagesc(NA1,deltaNA1, LLS_peak_mean_ratio_map);colorbar;axis square;xlabel("NA");ylabel("deltaNA")

% example plot as on your notes
figure;scatter(deltaNA1Array./NA1Array, LLS_peak_mean_ratio_array);axis square;xlabel("deltaNA/NA");ylabel("Ratio")