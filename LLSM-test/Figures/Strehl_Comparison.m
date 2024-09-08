%% 
clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

PSFdet = getDetectionPSF;
PSFdet = PSFdet./(max(max(max(PSFdet))));

MinRadialOrder = 2;
MaxRadialOrder = 6;
PhaseAmplitude = 6*wavelength_exc/(2*pi); 
PupilNA = 0.65;

%% Compare 2D and 1D gaussian
NA1 = 0.6; %0.6 for gaussian
SWweighting = 1; %4/3 for equal OTF V2 LLS, 1/sqrt(2) for V1 LLS
Latticeweighting = [1,1,1,1,1,1]; % 1.9 for V2 LLS
k_apertureNA = NA1 * k_wave / n;

gassianPupil2D = zeros(N,N);
gaussian_mask = zeros(N,N);
gaussian_mask= (k_apertureNA).^2 >= abs((kx_exc.^2 + kz_exc.^2));
gassianPupil2D = exp( -(kx_exc.^2 + kz_exc.^2)/ ((k_apertureNA).^2) );
gassianPupil2D = gassianPupil2D.* gaussian_mask;

gassianPupil1D = zeros(N,N);
gaussian_mask = zeros(N,N);
gaussian_mask(:,(N+1)/2) = (k_apertureNA) >= abs(kz_exc(:,1));
gassianPupil1D(:,(N+1)/2) = exp( -(kz_exc(:,1).^2)/ ((k_apertureNA).^2) );
gassianPupil1D = gassianPupil1D.* gaussian_mask;

[gaussian1D_SRatio,gaussian2D_SRatio,RadioOrderArray,AngularFrequencyArray] = ...
    Strehl(gassianPupil1D,gassianPupil2D,MinRadialOrder,MaxRadialOrder,PhaseAmplitude,PSFdet,PupilNA);
legend("1d Gaussian","2d Guassian")
%% Compare PEARLS and 0.6 NA 1D Gaussian
NA1 = 0.58; %0.6 for gaussian
deltaNA = 0.04;
SWweighting = 1/sqrt(2); %4/3 for equal OTF V2 LLS, 1/sqrt(2) for V1 LLS

SWPupil = zeros(N,N,2);
[SWPupil,~,~] = GetSWPairPupil('gaussian',NA1,NA1/2,...
deltaNA,2*deltaNA,...
SWweighting);

NA1 = 0.6; %0.6 for gaussian
SWweighting = 1; %4/3 for equal OTF V2 LLS, 1/sqrt(2) for V1 LLS
Latticeweighting = [1,1,1,1,1,1]; % 1.9 for V2 LLS
k_apertureNA = NA1 * k_wave / n;

gassianPupil1D = zeros(N,N);
gaussian_mask = zeros(N,N);
gaussian_mask(:,(N+1)/2) = (k_apertureNA) >= abs(kz_exc(:,1));
gassianPupil1D(:,(N+1)/2) = exp( -(kz_exc(:,1).^2)/ ((k_apertureNA).^2) );
gassianPupil1D = gassianPupil1D.* gaussian_mask;

[SWPupil_SRatio,gaussian1D_SRatio,RadioOrderArray,AngularFrequencyArray] = ...
    Strehl(SWPupil,gassianPupil1D,MinRadialOrder,MaxRadialOrder,PhaseAmplitude,PSFdet,PupilNA);
legend("0.6NA hex PEARLS","1d Gaussian")

%% Compare PEARLS and Bessel
NA1 = 0.58; %0.6 for gaussian
deltaNA = 0.04;
SWweighting = 4/3; %4/3 for equal OTF V2 LLS, 1/sqrt(2) for V1 LLS

SWPupil = zeros(N,N,2);
[SWPupil,~,~] = GetSWPairPupil('gaussian',NA1,NA1/2,...
deltaNA,2*deltaNA,...
SWweighting);

NA1 = 0.575; %0.6 for gaussian
deltaNA = 0.05;
SWweighting = 1; %4/3 for equal OTF V2 LLS, 1/sqrt(2) for V1 LLS
Latticeweighting = [1,1,1,1,1,1]; % 1.9 for V2 LLS
k_apertureNAmax = (NA1+ deltaNA/2) * k_wave / n;
k_apertureNAmin = (NA1-deltaNA/2) * k_wave / n;

BesselPupil = zeros(N,N);
BesselPupil = (k_apertureNAmax) > sqrt(kx_exc.^2 + kz_exc.^2) & (k_apertureNAmin) < sqrt(kx_exc.^2 + kz_exc.^2);

[SWPupil_SRatio,Bessel_SRatio,RadioOrderArray,AngularFrequencyArray] = ...
    Strehl(SWPupil,BesselPupil,MinRadialOrder,MaxRadialOrder,PhaseAmplitude,PSFdet,PupilNA);
legend("0.6NA hex PEARLS","Bessel")

%% Compare PEARLS and Double SW 
NA1 = 0.58; %0.6 for gaussian
deltaNA = 0.04;
SWweighting = 4/3; %4/3 for equal OTF V2 LLS, 1/sqrt(2) for V1 LLS

SWPupil = zeros(N,N,2);
[SWPupil,~,~] = GetSWPairPupil('gaussian',NA1,NA1/2,...
deltaNA,2*deltaNA,...
SWweighting);

DoubleSWPupil = zeros(N,N,2);
DoubleSWPupil = SWPupil(:,:,1) + SWPupil(:,:,2);

[SWPupil_SRatio,DoubleSW_SRatio,RadioOrderArray,AngularFrequencyArray] = ...
    Strehl(SWPupil,DoubleSWPupil,MinRadialOrder,MaxRadialOrder,PhaseAmplitude,PSFdet,PupilNA);
legend("0.6NA hex PEARLS","0.6NA Double SW")
