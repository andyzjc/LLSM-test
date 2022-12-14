%% Simulation without aberrration
clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

%% Simulation of SW Pair
% generate pupil 
% [SWPupil,SWMask] = GetSWPairPupil(ProfileType,NA1Ideal,NA2Ideal,deltaNA1,deltaNA2,...
%                          NA1Weighting,WeightRatio=I(NA1)/I(NA2));
[SWPupil,SWMask,SWPupilMetaData] = GetSWPairPupil('tophat',0.4,0.2,...
                                                           0.04,0.08,...
                                                           1,1);

% Coherent/Incoherent Propagation of SW pairs
% [PSFCoherent,PSFIncoherent] = SimulateSWPair(SWPupil);
[PSFCoherent,PSFIncoherent] = SimulateSWPair(SWPupil);

% Plot SWPairPSF
PrettyPlotSWPair(SWPupil,SWMask,SWPupilMetaData,PSFCoherent,PSFIncoherent);

%% Simulation of Lattice
% Create general Lattice
% [LatticePupil,LatticeMask] = GetLatticePupil(LatticeType,ProfileType,NAIdeal,deltaNA,Weighting);
[LatticePupil,LatticeMask,LatticeMetaData] = GetLatticePupil('hex','tophat',0.4,0.04,1);

% Simulate lattice
% [LatticePSF_3D] = SimulateLattice(LatticePupil)
[LatticePSF,LatticePSFDithered] = SimulateLattice(LatticePupil);

%Plot Lattice
%PrettyPlotLattice(LatticePupil,LatticeMask,LatticeMetaData,LatticePSF,LatticePSFDithered);
PrettyPlotLattice(LatticePupil,LatticeMask,LatticeMetaData,LatticePSF,LatticePSFDithered);

%% Simulation of SW Pair with a single aberrration mode 
% Get single zernike mode function
% Phase = GetSingleZmodePupil(RadialOrder,AngularFrequency,PhaseAmplitude);
Phase = GetSingleZmodePupil(1,-1,10);

% [SWPupil,SWMask] = GetSWPairPupil(ProfileType,NA1Ideal,NA2Ideal,deltaNA1,deltaNA2,...
%                          NA1Weighting,WeightRatio=I(NA1)/I(NA2));
[SWPupil,SWMask,SWPupilMetaData] = GetSWPairPupil('tophat',0.4,0.1265,...
                                                           0.08,0.253,...
                                                           1,1); 
AberratedSWPupil = SWPupil .* Phase;

% Simulate Single mode aberrated SW Pair
% [PSFCoherent,PSFIncoherent] = SimulateSWPair(SWPupil);
[AberratedPSFCoherent,AberratedPSFIncoherent] = imulateSWPair(AberratedSWPupil);

% Plot SWPairPSF
PrettyPlotSWPair(AberratedSWPupil,SWMask,SWPupilMetaData,AberratedPSFCoherent,AberratedPSFIncoherent);

%% Simulation of Lattice with a single aberrration mode 
% Get single zernike mode function
% Phase = GetSingleZmodePupil(RadialOrder,AngularFrequency,PhaseAmplitude);
Phase = GetSingleZmodePupil(1,-1,10);

% [LatticePupil,LatticeMask] = GetLatticePupil(LatticeType,ProfileType,NAIdeal,deltaNA,Weighting);
[LatticePupil,LatticeMask,LatticeMetaData] = GetLatticePupil('hex','tophat',0.4,0.04,1);
AberratedLatticePupil = LatticePupil.*Phase;

% Simulate Single mode aberrated lattice
% [LatticePSF_3D] = SimulateLattice(LatticePupil)
[AberratedLatticePSF,AberratedLatticePSFDithered] =  SimulateLattice(AberratedLatticePupil);

%Plot Lattice
%PrettyPlotLattice(LatticePupil,LatticeMask,LatticeMetaData,LatticePSF,LatticePSFDithered);
PrettyPlotLattice(AberratedLatticePupil,LatticeMask,LatticeMetaData,AberratedLatticePSF,AberratedLatticePSFDithered);

%% Simulation with all aberration mode
% Show all zernike mode function
GetZmodePupil(6);

% Simulate aberration for SW Pair
% [AberratedPSFCoherent,AberratedPSFIncoherent] = SimulateSWPairAberration(SWPupil,MaxRadialOrder,PhaseAmplitude)
SimulateSWPairAberration(SWPupil,6,10);

% Simulate aberration for Lattice
% [PSFCoherent,PSFIncoherent] = SimulateLatticeAberration(LatticePupil,MaxRadialOrder,PhaseAmplitude)
SimulateLatticeAberration(LatticePupil,6,10);




