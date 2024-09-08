%% 
clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

% %% detection 
% % detection
% PSFdet = getDetectionPSF;
% PSFdet = PSFdet./(max(max(max(PSFdet))));
% 
% % xzPSFdet = PSFdet(:,:,(N+1)/2);
% % yzPSFdet = squeeze(PSFdet(:,(N+1)/2,:)); 
% % xzOTFdet = fftshift(fft2(ifftshift(xzPSFdet)));
% % zOTFdet = real(xzOTFdet(:,(N+1)/2));

%% 
NA1 = 0.58;
deltaNA = 0.04;
LatticeType = 'hex';
ProfileType = 'tophat';
PhaseAmplitude = 6*wavelength_exc/(2*pi); 
FSWeighting1 = [1 0 0 0 0 1];
FSWeighting2 = [0 1 0 0 1 0]; % 1.9 for V2 LLS
FSWeighting3 = [0 0 1 1 0 0]; 

%% FS 
FSPupil = zeros(N,N,3);
[FSPupil(:,:,1),~,~] = GetLatticePupil(LatticeType,ProfileType, ...
NA1,deltaNA, ...
0.8,0,...
FSWeighting1);

[FSPupil(:,:,2),~,~] = GetLatticePupil(LatticeType,ProfileType, ...
NA1,deltaNA, ...
0.8,0,...
FSWeighting2);

[FSPupil(:,:,3),~,~] = GetLatticePupil(LatticeType,ProfileType, ...
NA1,deltaNA, ...
0.8,0,...
FSWeighting3);

[Beam1PSF,~,~] = SimulateLattice(FSPupil(:,:,1));
[Beam2PSF,~,~] = SimulateLattice(FSPupil(:,:,2));
[Beam3PSF,~,~] = SimulateLattice(FSPupil(:,:,3));
FSPSF = Beam1PSF + Beam2PSF + Beam3PSF;
FScenter = max(FSPSF,[],'all');
FSPSF = FSPSF/max(FSPSF,[],'all');

Phase = GetSingleZmodePupil(2,-2,0.65);
ComplexPhase = exp(PhaseAmplitude .* 1i .* Phase);

AberratedFSPupil = FSPupil .* ComplexPhase;
[AberratedBeam1PSF,~,~] = SimulateLattice(AberratedFSPupil(:,:,1));
[AberratedBeam2PSF,~,~] = SimulateLattice(AberratedFSPupil(:,:,2));
[AberratedBeam3PSF,~,~] = SimulateLattice(AberratedFSPupil(:,:,3));
AberratedFSPSF = AberratedBeam1PSF + AberratedBeam2PSF + AberratedBeam3PSF;
AberratedFSPSF = AberratedFSPSF/FScenter;

% [SW_SRatio,Lattice_SRatio,RadioOrderArray,AngularFrequencyArray] = Strehl(FSPupil,SWPupil(:,:,1)+SWPupil(:,:,2),0,6,PhaseAmplitude,PSFdet,0.65);

%% plot
savingdir = [LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_' ProfileType '/'];
mkdir(savingdir) 
FSsavingdir = [savingdir 'FS_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_WR' num2str(1) '/'];

mkdir(FSsavingdir) 

grapthSW(NA1,deltaNA,LatticeType,1,AberratedFSPSF,AberratedFSPSF,FSPupil,PSFdet,FSsavingdir )



