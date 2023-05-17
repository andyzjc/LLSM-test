%% 
clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

%% Getting phase
PhaseAmplitude1 = 6;
PhaseAmplitude2 = 3;
[~,Phase1] = GetSingleZmodePupil(3,-1,1);
[~,Phase2] = GetSingleZmodePupil(4,0,1);
Phase = PhaseAmplitude1* Phase1+ PhaseAmplitude2*Phase2;
ComplexPhase = exp(1i .* Phase);

[SWPupil1,~,~] = GetSWPairPupil('gaussian',0.58,0,...
                                           0.04,0,...
                                           1);
AberratedSWPupil1 = SWPupil1 .* ComplexPhase;

[SWPupil2,~,~] = GetSWPairPupil('gaussian',0,0.29,...
                                           0,0.08,...
                                           1);
AberratedSWPupil2 = SWPupil2 .* ComplexPhase;

% Coherent/Incoherent Propagation of SW pairs
% [PSFCoherent,PSFIncoherent] = SimulateSWPair(SWPupil);
[~,PSFIncoherent1,SWcenter1] = SimulateSWPair(AberratedSWPupil1);
PSFIncoherent1 = PSFIncoherent1/max(PSFIncoherent1,[],'all');

[~,PSFIncoherent2,SWcenter2] = SimulateSWPair(AberratedSWPupil2);
PSFIncoherent2 = PSFIncoherent2/max(PSFIncoherent2,[],'all');

AberratedPSFIncoherent = PSFIncoherent1+PSFIncoherent2;
AberratedPSFIncoherent = AberratedPSFIncoherent/max(AberratedPSFIncoherent,[],'all');

[~,PSFIncoherent,SWcenter2] = SimulateSWPair(SWPupil1+SWPupil2);
PSFIncoherent = PSFIncoherent/max(PSFIncoherent,[],'all');

figure
imagesc(Y_exc,Z_exc,squeeze(PSFIncoherent1(:,(N+1)/2,:)))
colormap(fire(256))
axis image

figure
imagesc(Y_exc,Z_exc,squeeze(PSFIncoherent2(:,(N+1)/2,:)))
colormap(fire(256))
axis image

figure
imagesc(Y_exc,Z_exc,squeeze(AberratedPSFIncoherent(:,(N+1)/2,:)))
colormap(fire(256))
axis image

figure
imagesc(Y_exc,Z_exc,squeeze(PSFIncoherent(:,(N+1)/2,:)))
colormap(fire(256))
axis image


