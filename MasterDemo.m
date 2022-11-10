clear all
close all

%% Simulation without aberrration
addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

%% SW Pair
% generate pupil 
% [SWPupil,SWMask] = GetSWPairPupil(ProfileType,NA1Ideal,NA2Ideal,deltaNA1,deltaNA2,...
%                          NA1Weighting,WeightRatio=I(NA1)/I(NA2));
[SWPupil,SWMask,SWPupilMetaData] = GetSWPairPupil('gaussian',0.4,0.126,0.08,0.253,1,1);

% Coherent/Incoherent Propagation of SW pairs
% [PSFCoherent,PSFIncoher
% ent] = SimulateSWPair(SWPupil);
[PSFCoherent,PSFIncoherent] = SimulateSWPair(SWPupil);

% Plot SWPairPSF
PrettyPlotSWPair(SWPupil,SWMask,SWPupilMetaData,PSFCoherent,PSFIncoherent);

%% Lattice
% Create general Lattice
% [LatticePupil,LatticeMask] = GetLatticePupil(LatticeType,ProfileType,NAIdeal,deltaNA,Weighting);
[LatticePupil,LatticeMask,LatticeMetaData] = GetLatticePupil('square','tophat',0.4,0.08,1);

% Simulate lattice
% [LatticePSF_3D] = SimulateLattice(LatticePupil)
[LatticePSF,LatticePSFDithered] = SimulateLattice(LatticePupil);

%Plot Lattice
%PrettyPlotLattice(LatticePupil,LatticeMask,LatticeMetaData,LatticePSF,LatticePSFDithered);
PrettyPlotLattice(LatticePupil,LatticeMask,LatticeMetaData,LatticePSF,LatticePSFDithered);

%% Simulation with aberrration

% Get zernike mode pupil function 
% [ZmodePupil] = GetZmodePupil(ZModeOrder);
[ZmodePupil] = GetZmodePupil(1);

% Simulate aberration for SW Pair
% [PSFCoherent,PSFIncoherent] = SimulateSW PairAberrated(SWPupil,ZModePupil,Plot=1/0,Save=1/0)
[AberratedSWPSF_Coherent,AberratedSWPSF_Incoherent] = SimulateSWPairAberrated(SWPupil,ZModePupil,1,1);

% Simulate aberration for Lattice
% [AberratedPSF] = SimulateLattice(LatticePupil,ZModePupil,Plot=1/0,Save=1/0)
[AberratedLatticePSF] = SimulateLatticeAberrated(LatticePupil,ZModePupil,1,1);




