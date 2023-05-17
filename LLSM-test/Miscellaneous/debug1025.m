%% Simulation without aberrration
clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

%% detection PSF
PSFdet = getDetectionPSF;
PSFdet = PSFdet./(max(max(max(PSFdet))));

%%
[SWPupil,SWMask,SWPupilMetaData] = GetSWPairPupil('tophat',0.58,0.29,...
                                                           0.04,0.08,...
                                                           7/10);

[LatticePupil,LatticeMask,LatticeMetaData] = GetLatticePupil('hex','tophat', ...
                                                             0.58,0.04, ...
                                                             0.8,0.0,...
                                                             1);
SWPSF = abs(fftshift(ifft2(ifftshift(SWPupil(:,:,1))))).^2 + abs(fftshift(ifft2(ifftshift(SWPupil(:,:,2))))).^2;
LatticePSF = abs(fftshift(ifft2(ifftshift(LatticePupil)))).^2;
LatticePSF = meshgrid(mean(LatticePSF,2))';
[SWint,~] = max(SWPSF,[],'all');
[Latticint,~] = max(LatticePSF,[],'all');

[ComplexPhase,Phase] = GetSingleZmodePupil(0,0,6);


AberratedSWPupil = SWPupil.*ComplexPhase;
AberratedLatticePupil= LatticePupil.*ComplexPhase;


AberratedSWPSF = abs(fftshift(ifft2(ifftshift(AberratedSWPupil(:,:,1))))).^2 + abs(fftshift(ifft2(ifftshift(AberratedSWPupil(:,:,2))))).^2;
AberratedSWPSF = AberratedSWPSF/SWint;

AberratedLatticePSF = abs(fftshift(ifft2(ifftshift(AberratedLatticePupil)))).^2;
AberratedLatticePSF = meshgrid(mean(AberratedLatticePSF,2))';
AberratedLatticePSF = AberratedLatticePSF/Latticint;

% OverallAberratedSWPSF = AberratedSWPSF.* PSFdet(:,:,(N+1)/2);
% OverallAberratedLatticePSF = AberratedLatticePSF.* PSFdet(:,:,(N+1)/2);

OverallAberratedSWPSF = AberratedSWPSF;
OverallAberratedLatticePSF = AberratedLatticePSF;

% figure
% imagesc(SWPSF)
% colorbar
% colormap(fire(256))
% 
% figure
% imagesc(LatticePSF)
% colorbar
% colormap(fire(256))
% 
% figure
% imagesc(AberratedSWPSF)
% colorbar
% colormap(fire(256))
% clim([0,1])
% 
% figure
% imagesc(AberratedLatticePSF)
% colorbar
% colormap(fire(256))
% clim([0,1])
% 
% figure
% imagesc(OverallAberratedSWPSF)
% [SWint,SWindex]= max(OverallAberratedSWPSF(:,(N+1)/2),[],'all');
% colorbar
% colormap(fire(256))
% clim([0,1])
% 
% figure
% imagesc(OverallAberratedLatticePSF)
% [Latticeint,Latticeindex] = max(OverallAberratedLatticePSF(:,(N+1)/2),[],'all');
% colorbar
% colormap(fire(256))
% clim([0,1])
% 
% figure
% imagesc(imfuse(Phase/2/pi,SWPupil(:,:,1)+SWPupil(:,:,2)))
% 
% figure
% imagesc(imfuse(Phase/2/pi,LatticePupil))

AberratedSWPSF((N+1)/2,(N+1)/2)
AberratedLatticePSF((N+1)/2,(N+1)/2)

figure
plot(Z_exc,AberratedSWPSF(:,(N+1)/2))
hold on
plot(Z_exc,AberratedLatticePSF(:,(N+1)/2))
legend("A.SW","A.LLS")
xlim([-1,1])

