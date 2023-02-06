%% Simulation without aberrration
clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

%% Simulation of SW Pair
% generate pupil 
% [SWPupil,SWMask] = GetSWPairPupil(ProfileType,NA1Ideal,NA2Ideal,
%                                               deltaNA1,deltaNA2,...
%                                               WeightRatio=I(NA1)/I(NA2));
[SWPupil,SWMask,SWPupilMetaData] = GetSWPairPupil('tophat',0.3,0.15,...
                                                           0.08,0.16,...
                                                           1);

% Coherent/Incoherent Propagation of SW pairs
% [PSFCoherent,PSFIncoherent] = SimulateSWPair(SWPupil);
[PSFCoherent,PSFIncoherent,SWcenter] = SimulateSWPair(SWPupil);

% Plot SWPairPSF
PrettyPlotSWPair(SWPupil,SWMask,SWPupilMetaData,PSFCoherent,PSFIncoherent,SWcenter);

%% Simulation of Lattice
% Create general Lattice
% [LatticePupil,LatticeMask] = GetLatticePupil(LatticeType,ProfileType,
%                                                   NAIdeal,deltaNA,
%                                                   MaskNAmax,MaskNAmin,
%                                                   WeightingRatio=I(SW)/I(Lattice));
[LatticePupil,LatticeMask,LatticeMetaData] = GetLatticePupil('hex','tophat', ...
                                                             0.4,0.08, ...
                                                             0.44,0.36,...
                                                             1);

% Simulate lattice
% [LatticePSF_3D] = SimulateLattice(LatticePupil)
[LatticePSF,LatticePSFDithered,Latticecenter] = SimulateLattice(LatticePupil);

%Plot Lattice
%PrettyPlotLattice(LatticePupil,LatticeMask,LatticeMetaData,LatticePSF,LatticePSFDithered);
PrettyPlotLattice(LatticePupil,LatticeMask,LatticeMetaData,LatticePSF,LatticePSFDithered,Latticecenter);

%% Simulation of SW Pair with a single aberrration mode 
% Get single zernike mode function
% Phase = GetSingleZmodePupil(RadialOrder,AngularFrequency,PhaseAmplitude);
Phase = GetSingleZmodePupil(4,4,1);

% [SWPupil,SWMask] = GetSWPairPupil(ProfileType,NA1Ideal,NA2Ideal,
%                                               deltaNA1,deltaNA2,...
%                                               WeightRatio=I(NA1)/I(NA2));
[SWPupil,SWMask,SWPupilMetaData] = GetSWPairPupil('tophat',0.4,0.2,...
                                                           0.08,0.16,...
                                                           0.68); 
AberratedSWPupil = SWPupil .* Phase;

% Simulate Single mode aberrated SW Pair
% [PSFCoherent,PSFIncoherent] = SimulateSWPair(SWPupil);
[AberratedPSFCoherent,AberratedPSFIncoherent,AberratedSWcenter] = SimulateSWPair(AberratedSWPupil);

% Plot SWPairPSF
PrettyPlotSWPair(AberratedSWPupil,SWMask,SWPupilMetaData,AberratedPSFCoherent,AberratedPSFIncoherent,AberratedSWcenter);

%% Simulation of Lattice with a single aberrration mode 
% Get single zernike mode function
% Phase = GetSingleZmodePupil(RadialOrder,AngularFrequency,PhaseAmplitude);
Phase = GetSingleZmodePupil(4,4,1);

% [LatticePupil,LatticeMask] = GetLatticePupil(LatticeType,ProfileType,
%                                                   NAIdeal,deltaNA,
%                                                   MaskNAmax,MaskNAmin,
%                                                   WeightingRatio=I(SW)/I(Lattice));
[LatticePupil,LatticeMask,LatticeMetaData] = GetLatticePupil('hex','tophat', ...
                                                             0.4,0.08, ...
                                                             0.6,0.2,...
                                                             1);
AberratedLatticePupil = LatticePupil.*Phase;

% Simulate Single mode aberrated lattice
% [LatticePSF_3D] = SimulateLattice(LatticePupil)
[AberratedLatticePSF,AberratedLatticePSFDithered,AberratedLatticecenter] =  SimulateLattice(AberratedLatticePupil);

%Plot Lattice
%PrettyPlotLattice(LatticePupil,LatticeMask,LatticeMetaData,LatticePSF,LatticePSFDithered);
PrettyPlotLattice(AberratedLatticePupil,LatticeMask,LatticeMetaData,AberratedLatticePSF,AberratedLatticePSFDithered,AberratedLatticecenter);

%% Simulation with all aberration mode
% Show all zernike mode function
GetZmodePupil(6);

% Simulate aberration for SW Pair
% [AberratedPSFCoherent,AberratedPSFIncoherent] = SimulateSWPairAberration(SWPupil,MaxRadialOrder,PhaseAmplitude)
SimulateSWPairAberration(SWPupil,6,1);

% Simulate aberration for Lattice
% [PSFCoherent,PSFIncoherent] = SimulateLatticeAberration(LatticePupil,MaxRadialOrder,PhaseAmplitude)
SimulateLatticeAberration(LatticePupil,6,1);
       



