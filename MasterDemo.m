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

%% Simulation of SW Pair
% generate pupil 
% [SWPupil,SWMask] = GetSWPairPupil(ProfileType,NA1Ideal,NA2Ideal,
%                                               deltaNA1,deltaNA2,...
%                                               WeightRatio=I(NA1)/I(NA2));
[SWPupil,SWMask,SWPupilMetaData] = GetSWPairPupil('tophat',0.55,sqrt(2*0.04*0.55)/2,...
                                                           0.04,sqrt(2*0.04*0.55),...
                                                           4);

% Coherent/Incoherent Propagation of SW pairs
% [PSFCoherent,PSFIncoherent] = SimulateSWPair(SWPupil);
[PSFCoherent,PSFIncoherent,SWcenter] = SimulateSWPair(SWPupil);

PSFIncoherent = PSFIncoherent/max(PSFIncoherent,[],'all');
% Plot SWPairPSF
% PrettyPlotSWPair(SWPupil,SWMask,SWPupilMetaData,PSFCoherent,PSFIncoherent,SWcenter);

%% Simulation of Lattice
% Create general Lattice
% [LatticePupil,LatticeMask] = GetLatticePupil(LatticeType,ProfileType,
%                                                   NAIdeal,deltaNA,
%                                                   MaskNAmax,MaskNAmin,
%                                                   WeightingRatio=I(SW)/I(Lattice));
tic
[LatticePupil,LatticeMask,LatticeMetaData] = GetLatticePupil('hex','tophat', ...
                                                             0.435,0.07, ...
                                                             0.8,0.0,...
                                                             [0 1 0 0 1 0]);

% Simulate lattice
% [LatticePSF_3D] = SimulateLattice(LatticePupil)
[LatticePSF,LatticePSFDithered,Latticecenter] = SimulateLattice(LatticePupil);

%Plot Lattice
%PrettyPlotLattice(LatticePupil,LatticeMask,LatticeMetaData,LatticePSF,LatticePSFDithered);
PrettyPlotLattice(LatticePupil,LatticeMask,LatticeMetaData,LatticePSF,LatticePSFDithered,Latticecenter);

toc
%% Simulation of SW Pair with a single aberrration mode 
% Get single zernike mode function
% Phase = GetSingleZmodePupil(RadialOrder,AngularFrequency,PhaseAmplitude);
Phase = GetSingleZmodePupil(2,-2,1);
ComplexPhase = exp(PhaseAmplitude .* 1i .* phase);

% [SWPupil,SWMask] = GetSWPairPupil(ProfileType,NA1Ideal,NA2Ideal,
%                                               deltaNA1,deltaNA2,...
%                                               WeightRatio=I(NA1)/I(NA2));
[SWPupil,SWMask,SWPupilMetaData] = GetSWPairPupil('tophat',0.4,0.2,...
                                                           0.08,0.16,...
                                                           4/3); 
AberratedSWPupil = SWPupil .* ComplexPhase;

% Simulate Single mode aberrated SW Pair
% [PSFCoherent,PSFIncoherent] = SimulateSWPair(SWPupil);
[AberratedPSFCoherent,AberratedPSFIncoherent,AberratedSWcenter] = SimulateSWPair(AberratedSWPupil);

% Plot SWPairPSF
PrettyPlotSWPair(AberratedSWPupil,SWMask,SWPupilMetaData,AberratedPSFCoherent,AberratedPSFIncoherent,AberratedSWcenter);

%% Simulation of Lattice with a single aberrration mode 
% Get single zernike mode function
% Phase = GetSingleZmodePupil(RadialOrder,AngularFrequency,PhaseAmplitude);
Phase = GetSingleZmodePupil(2,-2,1);

% [LatticePupil,LatticeMask] = GetLatticePupil(LatticeType,ProfileType,
%                                                   NAIdeal,deltaNA,
%                                                   MaskNAmax,MaskNAmin,
%                                                   WeightingRatio=I(SW)/I(Lattice));
[LatticePupil,LatticeMask,LatticeMetaData] = GetLatticePupil('hex','tophat', ...
                                                             0.4,0.08, ...
                                                             0.6,0.2,...
                                                             [1 1 1 1 1 1]);
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
[SWCorrCoef,SWStrehlPeaks,SWStrehlCenter] ...
    = SimulateSWPairAberration(SWPupil,PSFdet,6,4*pi);

% Simulate aberration for Lattice
% [PSFCoherent,PSFIncoherent] = SimulateLatticeAberration(LatticePupil,MaxRadialOrder,PhaseAmplitude)
[LatticeCorrCoef,LatticeStrehlPeaks,LatticeStrehlCenter] ...
    = SimulateLatticeAberration(LatticePupil,PSFdet,6,4*pi);
       



