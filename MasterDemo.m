clear all
close all

%% Simulation without aberrration
addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters;
CalculatePhysics;

% generate pupil 
% [SWPupil,SWMask] = GetSWPairPupil(ProfileType,NA1Ideal,NA2Ideal,deltaNA1,deltaNA2,...
%                          Weighting,WeightRatio=I(NA1)/I(NA2));
[SWPupil,SWMask] = GetSWPairPupil('Tophat',0.4,0.2,0.08,0.16,1,1);

% Coherent/Incoherent Propagation of SW pairs
% [PSFCoherent,PSFIncoherent] = SimulateSWPair(SWPupil,Plot=1/0,Save=1/0);
[SWPairPSF_Coherent,SWPairPSF_Incoherent] = SimulateSWPair(SWPupil,1,1);

% Create general Lattice
% [LatticePupil,LatticeMask] = GetSWPairPupil(LatticeType,ProfileType,NAIdeal,deltaNA,Weighting);
[LatticePupil,LatticeMask] = GetLatticePupil('Hex','Tophat',0.4,0.08,1);

% Simulate lattice
% [LatticePSF_3D] = SimulateLattice(LatticePupil,Plot=1/0,Save=1/0)
[LatticePSF] = SimulateLattice(LatticePupil,1,1);

%% Simulation with aberrration

% Get zernike mode pupil function 
% [ZmodePupil] = GetZmodePupil(ZModeOrder);
[ZmodePupil] = GetZmodePupil(1);

% Simulate aberration for SW Pair
% [PSFCoherent,PSFIncoherent] = SimulateSWPairAberrated(SWPupil,ZModePupil,Plot=1/0,Save=1/0)
[AberratedSWPSF_Coherent,AberratedSWPSF_Incoherent] = SimulateSWPairAberrated(SWPupil,ZModePupil,1,1);

% Simulate aberration for Lattice
% [AberratedPSF] = SimulateLattice(LatticePupil,ZModePupil,Plot=1/0,Save=1/0)
[AberratedLatticePSF] = SimulateLatticeAberrated(LatticePupil,ZModePupil,1,1);




