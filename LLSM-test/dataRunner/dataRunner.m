%% 
clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

%% detection 
% detection
PSFdet = getDetectionPSF;
PSFdet = PSFdet./(max(max(max(PSFdet))));

% xzPSFdet = PSFdet(:,:,(N+1)/2);
% yzPSFdet = squeeze(PSFdet(:,(N+1)/2,:)); 
% xzOTFdet = fftshift(fft2(ifftshift(xzPSFdet)));
% zOTFdet = real(xzOTFdet(:,(N+1)/2));

%% Unaberrated

NA1 = 0.58;
deltaNA = 0.04;
LatticeType = 'hex';
ProfileType = 'tophat';
SWweighting = 4/3; %4/3 for equal OTF V2 LLS, 7/10 for V1 LLS
Latticeweighting = 1.9; % 1.9 for V2 LLS
SNR = 10;
Iter = 10;
OTFthreshold = 0.001;

savingdir = [LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_' ProfileType '/'];
mkdir(savingdir) 
SWsavingdir = [savingdir 'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_WR' num2str(SWweighting) '/'];
mkdir(SWsavingdir) 
Latticesavingdir = [savingdir LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_WR' num2str(Latticeweighting) '/'];
mkdir(Latticesavingdir) 
SWDatasavingdir = [savingdir 'SWdata/'];
mkdir(SWDatasavingdir) 
LLSDatasavingdir = [savingdir 'LLSdata/'];
mkdir(LLSDatasavingdir) 

if isequal(LatticeType,'hex')
    [SWPupil,SWmask,SWPupilMetaData] = GetSWPairPupil(ProfileType,NA1,NA1/2,...
    deltaNA,2*deltaNA,...
    SWweighting);

    [LatticePupil,~,~] = GetLatticePupil(LatticeType,ProfileType, ...
    NA1,deltaNA, ...
    0.8,0,...
    Latticeweighting);
elseif isequal(LatticeType,'SW')
    [SWPupil,SWmask,SWPupilMetaData] = GetSWPairPupil(ProfileType,NA1,NA1,...
    deltaNA,deltaNA,...
    SWweighting);

    [LatticePupil,~,~] = GetLatticePupil('hex',ProfileType, ...
    NA1,deltaNA, ...
    0.8,0,...
    Latticeweighting);
elseif isequal(LatticeType,'NC')
    [SWPupil,~,~] = GetSWPairPupil(ProfileType,NA1,NA1*3/4,...
    deltaNA,deltaNA*4/3,...
    SWweighting);

    [LatticePupil,~,~] = GetLatticePupil('hex',ProfileType, ...
    NA1,deltaNA, ...
    0.8,0,...
    Latticeweighting);
elseif isequal(LatticeType,'sinc')
    [SWPupil,~,~] = GetSWPairPupil('tophat',0,0,...
    NA1*2,NA1*2,...
    SWweighting);

    LatticePupil = load('LLSM-test/Aberration/Zernike/AiryPupil513_0p6.mat');
    temp = LatticePupil.AiryPupil;
    LatticePupil = zeros(size(temp));
    LatticePupil(:,(N+1)/2) = temp(:,(N+1)/2);
elseif isequal(LatticeType,'1dgaussian')
    [SWPupil,~,~] = GetSWPairPupil('gaussian',0,0,...
    0,NA1,...
    SWweighting);

    LatticePupil = load('LLSM-test/Aberration/Zernike/1dGaussianPupil513_NA_0p26.mat');
    LatticePupil = LatticePupil.GaussianPupil;
elseif isequal(LatticeType,'2dgaussian')
    [SWPupil,~,~] = GetSWPairPupil('gaussian',0,0,...
    0,NA1,...
    SWweighting);

    LatticePupil = load('LLSM-test/Aberration/Zernike/2dGaussianPupil513_NA_0p6.mat');
    LatticePupil = LatticePupil.GaussianPupil;
elseif isequal(LatticeType,'bessel')
    [SWPupil,~,~] = GetSWPairPupil('gaussian',NA1,NA1,...
    deltaNA,deltaNA,...
    SWweighting);

    LatticePupil = load('LLSM-test/Aberration/Zernike/1dBesselPupil513_NA_0p6_0p05.mat');
    LatticePupil = LatticePupil.BesselPupil;
else % square
    [SWPupil,~,~] = GetSWPairPupil(ProfileType,NA1,sqrt(2*NA1*deltaNA)/2,...
    deltaNA,sqrt(2*NA1*deltaNA),...
    SWweighting);

    [LatticePupil,~,~] = GetLatticePupil(LatticeType,ProfileType, ...
    NA1,deltaNA, ...
    0.8,0,...
    Latticeweighting);
end
%% SW
[PSFCoherent,PSFIncoherent,SWcenter] = SimulateSWPair(SWPupil);
PSFCoherent = PSFCoherent/max(PSFCoherent,[],'all');
PSFIncoherent = PSFIncoherent/max(PSFIncoherent,[],'all');
mkdir([SWsavingdir 'Unaberrated/'])
grapthSW(NA1,deltaNA,LatticeType,SWweighting,PSFIncoherent,PSFCoherent,SWPupil,PSFdet,[SWsavingdir 'Unaberrated/'])
    yindex = 1-(squeeze(PSFIncoherent((N+1)/2,(N+1)/2,:)) <= 0.5*max(squeeze(PSFIncoherent((N+1)/2,(N+1)/2,:))));
    SWyFWHM1 = find(yindex,1,'first') ;
    SWyFWHM2 = find(yindex,1,'last');

% get OTF mask for FC3
overallOTF = fftshift(ifftn(ifftshift(PSFIncoherent.* PSFdet)));
xzOTF = abs(squeeze(overallOTF(:,(N+1)/2,:)))/max(abs(squeeze(overallOTF(:,(N+1)/2,:))),[],'all');
SWOTFmask = xzOTF((N+1)/2:end,(N+1)/2:end) >= OTFthreshold;

[SWAveragefc3,SWAveragefc3FWHM,~,~] = RunFC3(PSFIncoherent,PSFdet,SWOTFmask,Iter,SNR,SWyFWHM2);
[SWconvLines,lineSpot,~,Line_Z] = ConvRes(PSFIncoherent,PSFdet);
save([SWDatasavingdir '/PSFIncoherent.mat'],'PSFIncoherent')
save([ SWDatasavingdir '/PSFIncoherent_Center.mat'], 'SWcenter')
save([ SWDatasavingdir '/PSFIncoherent_FC3.mat'], 'SWAveragefc3')
save([ SWDatasavingdir '/PSFIncoherent_FC3FWHM.mat'], 'SWAveragefc3FWHM')
save([ SWDatasavingdir '/PSFIncoherent_LineGrating.mat'], 'SWconvLines')

%% Lattice
[LatticePSF,LatticePSFDithered,Latticecenter] = SimulateLattice(LatticePupil);
LatticePSF = LatticePSF/max(LatticePSF,[],'all');
LatticePSFDithered = LatticePSFDithered/max(LatticePSFDithered,[],'all');
mkdir([Latticesavingdir 'Unaberrated/'])
grapthLattice(NA1,deltaNA,LatticeType,Latticeweighting,LatticePSF,LatticePSFDithered,LatticePupil,PSFdet,[Latticesavingdir 'Unaberrated/'])
    yindex = 1-(squeeze(LatticePSFDithered((N+1)/2,(N+1)/2,:)) <= 0.5*max(squeeze(LatticePSFDithered((N+1)/2,(N+1)/2,:))));
    LatticeyFWHM1 = find(yindex,1,'first');
    LatticeyFWHM2 = find(yindex,1,'last');

% get OTF mask for FC3
overallOTF = fftshift(ifftn(ifftshift(LatticePSFDithered.* PSFdet)));
xzOTF = abs(squeeze(overallOTF(:,(N+1)/2,:)))/max(abs(squeeze(overallOTF(:,(N+1)/2,:))),[],'all');
LatticeOTFmask = xzOTF((N+1)/2:end,(N+1)/2:end) >= OTFthreshold;

[LatticeAveragefc3,LatticeAveragefc3FWHM,~,~] = RunFC3(LatticePSFDithered,PSFdet,LatticeOTFmask,Iter,SNR,LatticeyFWHM2);
[LatticeconvLines,~,~,~] = ConvRes(LatticePSFDithered,PSFdet);
save([ LLSDatasavingdir '/LatticePSF.mat'], 'LatticePSF')
save([ LLSDatasavingdir '/LatticePSFDithered.mat'], 'LatticePSFDithered')
save([ LLSDatasavingdir '/Latticecenter.mat'], 'Latticecenter')
save([ LLSDatasavingdir '/LatticePSFDithered_FC3.mat'], 'LatticeAveragefc3')
save([ LLSDatasavingdir '/LatticePSFDithered_FC3FWHM.mat'], 'LatticeAveragefc3FWHM')
save([ LLSDatasavingdir '/LatticePSFDithered_LineGrating.mat'], 'LatticeconvLines')

%% Aberrated
clc
[theta,r] = cart2pol(kx_exc./(0.65./n*k_wave),kz_exc./(0.65./n*k_wave));
idx = r<=1;

MinRadialOrder = 2;
MaxRadialOrder = 2;
PhaseAmplitude = 6; 

RadioOrderArray = [];
AngularFrequencyArray =[];
counter =1;
for i = MinRadialOrder:MaxRadialOrder
    AngularFrequency = -i:2:i;
    for k = 1:length(AngularFrequency)
        RadioOrderArray(1,counter) = i;
        AngularFrequencyArray(1,counter) = AngularFrequency(k);
        counter = counter +1;
    end
end

LabelArray = cell(length(RadioOrderArray),1);
for i = 1:length(RadioOrderArray)
    LabelArray{i,1} = "Z^{" + num2str(AngularFrequencyArray(i)) + "}_{" + num2str(RadioOrderArray(i)) + "}";
end

RadioOrderArray = zeros(length(RadioOrderArray),1);
AngularFrequencyArray = RadioOrderArray;

AberratedPSFCoherent = cell(length(RadioOrderArray),1);
AberratedPSFIncoherent = AberratedPSFCoherent;
AberratedSWPupil = AberratedPSFCoherent;
AberratedSWcenter = AberratedPSFCoherent;

AberratedLatticePupil = AberratedPSFCoherent;
AberratedLatticecenter = AberratedPSFCoherent;
AberratedLatticePSF = AberratedPSFCoherent;
AberratedLatticePSFDithered = AberratedPSFCoherent;
%% SW
counter = 1;
for i = MinRadialOrder:MaxRadialOrder
    AngularFrequency = -i:2:i;
    for k = 1:length(AngularFrequency)
        RadioOrderArray(counter) = i;
        AngularFrequencyArray(counter) = AngularFrequency(k);
        phase = zeros(size(kx_exc));
        phase(idx) = zernfun(i,AngularFrequency(k),r(idx),theta(idx),'norm');

        %SW 
        AberratedPupil1 = zeros(size(phase));
        AberratedPupil2 = zeros(size(phase));
        SWPupil1 = squeeze(SWPupil(:,:,1));
        SWPupil2 = squeeze(SWPupil(:,:,2));

        AberratedPupil1(idx) = SWPupil1(idx) .* exp(PhaseAmplitude.* 1i.*phase(idx));
        AberratedPupil2(idx) = SWPupil2(idx) .* exp(PhaseAmplitude.* 1i.*phase(idx));
        temp1(:,:,1) = AberratedPupil1;
        temp1(:,:,2) = AberratedPupil2;
        [temp3,temp4,temp5] = SimulateSWPair(temp1);
        AberratedSWPupil{counter,1} = temp1;
        AberratedPSFCoherent{counter,1} = temp3 /SWcenter(1,1);
        AberratedPSFIncoherent{counter,1} = temp4 /SWcenter(2,1);
        AberratedSWcenter{counter,1} = temp5;
        SWAberratedsavingdir = [SWsavingdir 'Z_' num2str(i) '_' num2str(AngularFrequency(k)) '_Amplitude_' num2str(PhaseAmplitude) '/'];
        mkdir(SWAberratedsavingdir)
        grapthSW(NA1,deltaNA,LatticeType,SWweighting,AberratedPSFIncoherent{counter,1},AberratedPSFCoherent{counter,1},AberratedSWPupil{counter,1},...
           PSFdet,[SWAberratedsavingdir 'Z_' num2str(i) '_' num2str(AngularFrequency(k)) '_']);
        counter = counter +1;
    end
end

save([SWDatasavingdir '/AberratedPSFCoherent.mat'],'AberratedPSFCoherent','-v7.3')
save([SWDatasavingdir '/AberratedPSFIncoherent.mat'],'AberratedPSFIncoherent','-v7.3')
save([SWDatasavingdir '/AberratedSWPupil.mat'],'AberratedSWPupil','-v7.3')
save([SWDatasavingdir '/AberratedSWcenter.mat'],'AberratedSWcenter','-v7.3')

%% LLS
counter = 1;
for i = MinRadialOrder:MaxRadialOrder
    AngularFrequency = -i:2:i;
    for k = 1:length(AngularFrequency)
        RadioOrderArray(counter) = i;
        AngularFrequencyArray(counter) = AngularFrequency(k);
        phase = zeros(size(kx_exc));
        phase(idx) = zernfun(i,AngularFrequency(k),r(idx),theta(idx),'norm');
        % %Lattice
        temp2 = zeros(size(phase));
        temp2(idx) = LatticePupil(idx) .* exp(PhaseAmplitude.* 1i .*phase(idx));
        [temp6,temp7,temp8] = SimulateLattice(temp2);
        AberratedLatticePupil{counter,1} = temp2;
        AberratedLatticePSF{counter,1} = temp6/Latticecenter(1,1);
        AberratedLatticePSFDithered{counter,1} = temp7/Latticecenter(2,1);
        AberratedLatticecenter{counter,1} = temp8;
        LatticeAberratedsavingdir = [Latticesavingdir 'Z_' num2str(i) '_' num2str(AngularFrequency(k)) '_Amplitude_' num2str(PhaseAmplitude) '/'];
        mkdir(LatticeAberratedsavingdir)
        grapthLattice(NA1,deltaNA,LatticeType,Latticeweighting,AberratedLatticePSF{counter,1},AberratedLatticePSFDithered{counter,1},AberratedLatticePupil{counter,1},...
            PSFdet,[LatticeAberratedsavingdir 'Z_' num2str(i) '_' num2str(AngularFrequency(k)) '_'])
        counter = counter +1;
    end
end        
save([ LLSDatasavingdir '/AberratedLatticePupil.mat'], 'AberratedLatticePupil','-v7.3')
save([ LLSDatasavingdir '/AberratedLatticePSF.mat'], 'AberratedLatticePSF','-v7.3')
save([ LLSDatasavingdir '/AberratedLatticePSFDithered.mat'], 'AberratedLatticePSFDithered','-v7.3')
save([ LLSDatasavingdir '/AberratedLatticeCenter.mat'], 'AberratedLatticecenter','-v7.3')

%% Analysis
SW_SRatio_Focal = zeros(size(AberratedPSFIncoherent));
SW_SRatio_Focal_FWHM = SW_SRatio_Focal;
SW_SRatio_3D = SW_SRatio_Focal;
SW_SRatio_yz = SW_SRatio_Focal;
SW_SRatio_yz_overall = SW_SRatio_Focal;
SW_SRatio_corrected = SW_SRatio_Focal;

AberratedSWAveragefc3 = cell(length(RadioOrderArray),1);
AberratedSWAveragefc3FWHM = AberratedSWAveragefc3;
AberratedSWconvLines = AberratedSWAveragefc3;
for i = 1:length(AberratedPSFIncoherent)
    tempPSF = AberratedPSFIncoherent{i,1};
    tempSWcenter = AberratedSWcenter{i,1};

    % Srethl ratio
    SW_SRatio_Focal(i,1) = tempPSF((N+1)/2,(N+1)/2,(N+1)/2); % focal point, excitation
    SW_SRatio_Focal_FWHM(i,1) = tempPSF((N+1)/2,(N+1)/2,SWyFWHM2); % HWHM, axis, excitation
    tempyzPSF = squeeze(tempPSF(:,(N+1)/2,:));
    SW_SRatio_yz(i,1) = max(tempyzPSF,[],'all'); % yz excitation 

    [~,col] = find(tempyzPSF==max(max(tempyzPSF)));
    SW_SRatio_yz_overall(i,1) = max( PSFdet(:,(N+1)/2,(N+1)/2) .* tempyzPSF(:,col) ); % yz overall 
    
    SW_SRatio_3D(i,1) = max(tempPSF,[],'all');

    % aberration correction (quick guess) 
    SW_SRatio_corrected(i,1) = (tempSWcenter(3,1) + tempSWcenter(4,1))/SWcenter(2,1);

    % FC3 
    [AberratedSWAveragefc3{i,1},AberratedSWAveragefc3FWHM{i,1},~,~] = RunFC3(tempPSF,PSFdet,SWOTFmask,Iter,SNR,SWyFWHM2);

    % line grating
    [AberratedSWconvLines{i,1},~,~,~] = ConvRes(tempPSF,PSFdet);
end
save([SWDatasavingdir '/SW_SRatio_Focal.mat'],'SW_SRatio_Focal')
save([SWDatasavingdir '/SW_SRatio_FWHM.mat'],'SW_SRatio_Focal_FWHM')
save([SWDatasavingdir '/SW_SRatio_3D.mat'],'SW_SRatio_3D')
save([SWDatasavingdir '/SW_SRatio_corrected.mat'],'SW_SRatio_corrected')
save([SWDatasavingdir '/SW_SRatio_yz.mat'],'SW_SRatio_yz')
save([SWDatasavingdir '/SW_SRatio_yz_overall.mat'],'SW_SRatio_yz_overall')

save([SWDatasavingdir '/AberratedSWAveragefc3.mat'],'AberratedSWAveragefc3')
save([SWDatasavingdir '/AberratedSWAveragefc3FWHM.mat'],'AberratedSWAveragefc3FWHM')
save([SWDatasavingdir '/AberratedSWconvLines.mat'],'AberratedSWconvLines')

Lattice_SRatio_Focal = zeros(size(AberratedLatticePSFDithered));
Lattice_SRatio_Focal_FWHM = Lattice_SRatio_Focal;
Lattice_SRatio_yz = Lattice_SRatio_Focal;
Lattice_SRatio_yz_overall = Lattice_SRatio_Focal;
Lattice_SRatio_3D = Lattice_SRatio_Focal;

AberratedLatticeAveragefc3 = AberratedSWAveragefc3;
AberratedLatticeAveragefc3FWHM = AberratedLatticeAveragefc3;
AberratedLatticeconvLines = AberratedLatticeAveragefc3;
for i = 1:length(AberratedLatticePSFDithered)
    tempPSF = AberratedLatticePSFDithered{i,1};

    % Srethl ratio
    Lattice_SRatio_Focal(i,1) = tempPSF((N+1)/2,(N+1)/2,(N+1)/2);
    Lattice_SRatio_Focal_FWHM(i,1) = tempPSF((N+1)/2,(N+1)/2,LatticeyFWHM2);
    Lattice_SRatio_yz(i,1) = max(squeeze(tempPSF(:,(N+1)/2,:)),[],'all'); % yz excitation  

    tempyzPSF = squeeze(tempPSF(:,(N+1)/2,:));
    [~,col] = find(tempyzPSF==max(max(tempyzPSF)));
    Lattice_SRatio_yz_overall(i,1) = max( PSFdet(:,(N+1)/2,(N+1)/2) .* tempyzPSF(:,col) ); % yz overall 

    Lattice_SRatio_3D(i,1) = max(tempPSF,[],'all');

    % FC3 
    [AberratedLatticeAveragefc3{i,1},AberratedLatticeAveragefc3FWHM{i,1},~,~] = RunFC3(tempPSF,PSFdet,LatticeOTFmask,Iter,SNR,LatticeyFWHM2);

    % line grating
    [AberratedLatticeconvLines{i,1},~,~,~] = ConvRes(tempPSF,PSFdet);
end
save([ LLSDatasavingdir '/Lattice_SRatio_Focal.mat'], 'Lattice_SRatio_Focal')
save([ LLSDatasavingdir '/Lattice_SRatio_FWHM.mat'], 'Lattice_SRatio_Focal_FWHM')
save([ LLSDatasavingdir '/Lattice_SRatio_yz.mat'], 'Lattice_SRatio_yz')
save([ LLSDatasavingdir '/Lattice_SRatio_yz_overall.mat'], 'Lattice_SRatio_yz_overall')
save([ LLSDatasavingdir '/Lattice_SRatio_3D.mat'], 'Lattice_SRatio_3D')
save([ LLSDatasavingdir '/AberratedLatticeAveragefc3.mat'], 'AberratedLatticeAveragefc3')
save([ LLSDatasavingdir '/AberratedLatticeAveragefc3FWHM.mat'], 'AberratedLatticeAveragefc3FWHM')
save([ LLSDatasavingdir '/AberratedLatticeconvLines.mat'], 'AberratedLatticeconvLines')

%% Pretty plots, Pupil & Strehl
clc

Analysis_savingdir = [savingdir 'AberrationAnalysis/'];
mkdir(Analysis_savingdir) 

Pupilsavingdir = [Analysis_savingdir 'PupilError/'];
mkdir(Pupilsavingdir)
counter = 1;
for i = MinRadialOrder:MaxRadialOrder
    AngularFrequency = -i:2:i;
    for k = 1:length(AngularFrequency)
        RadioOrderArray(1,counter) = i;
        AngularFrequencyArray(1,counter) = AngularFrequency(k);
        counter = counter +1;
        phase = zeros(size(kx_exc));
        phase(idx) = zernfun(i,AngularFrequency(k),r(idx),theta(idx),'norm');
        
        fig1 = figure;
        imagesc(KX_exc,KZ_exc,PhaseAmplitude*phase/2/pi)
        xlim([-0.5,0.5])
        ylim([-0.5,0.5])
        title("0.65NA")
        colormap(jet)
        colorbar
        print(fig1, '-dsvg', [ Pupilsavingdir  'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(k)) '_WFE_Amplitude_' num2str(PhaseAmplitude) '.SVG'],'-r300')
        print(fig1, '-dpng', [ Pupilsavingdir 'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(k)) '__WFE_Amplitude_' num2str(PhaseAmplitude) '.PNG'],'-r300')
    end
end
close all

%% 
Strehlsavingdir = [Analysis_savingdir 'StrehlRatio/'];
mkdir(Strehlsavingdir)

fig1 = figure;
    h1 = subplot(1,1,1,'Parent',fig1);
    plot(1:length(SW_SRatio_Focal),SW_SRatio_Focal,'Parent',h1,'LineStyle','-','Marker','o');
    hold on
    plot(1:length(Lattice_SRatio_Focal),Lattice_SRatio_Focal,'Parent',h1,'LineStyle','-.','Marker','o');
    lgd = legend("iSW","LLS");
    lgd.Location = 'northoutside';
    h1.XAxis.TickValues = 1:length(RadioOrderArray);
    h1.XAxis.TickLabels = LabelArray;
    h1.XAxis.FontSize = 6;
    grid on
    xlabel("Aberration Mode")
    ylabel("Strehl Ratio")
    title("Focal Plane")
    pbaspect([5 1 1])
    hold off
    print(fig1, '-dsvg', [Strehlsavingdir LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_FocalStrehl_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.SVG'],'-r300')
    print(fig1, '-dpng', [Strehlsavingdir LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_FocalStrehl_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.PNG'],'-r300')

fig2 = figure;
    h1 = subplot(1,1,1,'Parent',fig2);
    plot(1:length(SW_SRatio_Focal_FWHM),SW_SRatio_Focal_FWHM,'Parent',h1,'LineStyle','-.','Marker','o');
    hold on
    plot(1:length(Lattice_SRatio_Focal_FWHM),Lattice_SRatio_Focal_FWHM,'Parent',h1,'LineStyle','-.','Marker','o');
    % p1 = plot(1:length(SW_SRatio_Focal),SW_SRatio_Focal./Lattice_SRatio_Focal,'magenta');
    lgd = legend("iSW","LLS");
    lgd.Location = 'northoutside';
    h1.XAxis.TickValues = 1:length(RadioOrderArray);
    h1.XAxis.TickLabels = LabelArray;
    h1.XAxis.FontSize = 6;
    grid on
    xlabel("Aberration Mode")
    ylabel("Strehl Ratio")
    title("yFWHM Plane")
    pbaspect([5 1 1])
    hold off
    print(fig2, '-dsvg', [ Strehlsavingdir  LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_yFWHMStrehl_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.SVG'],'-r300')
    print(fig2, '-dpng', [ Strehlsavingdir  LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_yFWHMStrehl_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.PNG'],'-r300')

fig3 = figure;
    h1 = subplot(1,1,1,'Parent',fig3);
    plot(1:length(SW_SRatio_3D),SW_SRatio_3D,'Parent',h1,'LineStyle','-.','Marker','o');
    hold on
    plot(1:length(Lattice_SRatio_3D),Lattice_SRatio_3D,'Parent',h1,'LineStyle','-.','Marker','o');
    lgd = legend("iSW","LLS");
    lgd.Location = 'northoutside';
    h1.XAxis.TickValues = 1:length(RadioOrderArray);
    h1.XAxis.TickLabels = LabelArray;
    h1.XAxis.FontSize = 6;
    grid on
    xlabel("Aberration Mode")
    ylabel("Strehl Ratio")
    title("3D Excitation")
    pbaspect([5 1 1])
    hold off
    print(fig3, '-dsvg', [ Strehlsavingdir   LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_3DStrehl_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.SVG'],'-r300')
    print(fig3, '-dpng', [ Strehlsavingdir  LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_3DStrehl_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.PNG'],'-r300')

fig4 = figure;
    h1 = subplot(1,1,1,'Parent',fig4);
    plot(1:length(SW_SRatio_Focal),SW_SRatio_Focal,'Parent',h1,'LineStyle','-.','Marker','o');
    hold on
    plot(1:length(SW_SRatio_corrected),SW_SRatio_corrected,'Parent',h1,'LineStyle','-.','Marker','o');
    lgd = legend("Aberrated","Corrected");
    lgd.Location = 'northoutside';
    h1.XAxis.TickValues = 1:length(RadioOrderArray);
    h1.XAxis.TickLabels = LabelArray;
    h1.XAxis.FontSize = 6;
    grid on
    xlabel("Aberration Mode")
    ylabel("Strehl Ratio")
    title("Correction")
    pbaspect([5 1 1])
    hold off
    print(fig4, '-dsvg', [ Strehlsavingdir   LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_CorrectedStrehl_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.SVG'],'-r300')
    print(fig4, '-dpng', [ Strehlsavingdir  LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_CorrectedStrehl_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.PNG'],'-r300')

fig5 = figure;
    h1 = subplot(1,1,1,'Parent',fig5);
    plot(1:length(SW_SRatio_yz),SW_SRatio_yz,'Parent',h1,'LineStyle','-.','Marker','o');
    hold on
    plot(1:length(Lattice_SRatio_yz),Lattice_SRatio_yz,'Parent',h1,'LineStyle','-.','Marker','o');
    lgd = legend("iSW","LLS");
    lgd.Location = 'northoutside';
    h1.XAxis.TickValues = 1:length(RadioOrderArray);
    h1.XAxis.TickLabels = LabelArray;
    h1.XAxis.FontSize = 6;
    grid on
    xlabel("Aberration Mode")
    ylabel("Strehl Ratio")
    title("yz Excitation")
    pbaspect([5 1 1])
    hold off
    print(fig5, '-dsvg', [ Strehlsavingdir   LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_yzExcitationStrehl_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.SVG'],'-r300')
    print(fig5, '-dpng', [ Strehlsavingdir  LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_yzExcitationStrehl_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.PNG'],'-r300')

fig6 = figure;
    h1 = subplot(1,1,1,'Parent',fig6);
    plot(1:length(SW_SRatio_yz_overall),SW_SRatio_yz_overall,'Parent',h1,'LineStyle','-.','Marker','o');
    hold on
    plot(1:length(Lattice_SRatio_yz_overall),Lattice_SRatio_yz_overall,'Parent',h1,'LineStyle','-.','Marker','o');
    lgd = legend("iSW","LLS");
    lgd.Location = 'northoutside';
    h1.XAxis.TickValues = 1:length(RadioOrderArray);
    h1.XAxis.TickLabels = LabelArray;
    h1.XAxis.FontSize = 6;
    grid on
    xlabel("Aberration Mode")
    ylabel("Strehl Ratio")
    title("yz Overall")
    pbaspect([5 1 1])
    hold off
    print(fig6, '-dsvg', [ Strehlsavingdir   LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_yzOverallStrehl_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.SVG'],'-r300')
    print(fig6, '-dpng', [ Strehlsavingdir  LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_yzOverallStrehl_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.PNG'],'-r300')
close all

%% Pretty plots, FC3 
clc
SWFC3savingdir = [Analysis_savingdir 'FC3/iSW/'];
mkdir(SWFC3savingdir)
fig1 = figure;
for i = 1:length(RadioOrderArray)
    h1 = subplot(2,3,1,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),SWAveragefc3)
    title("Unaberrated")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    h1 = subplot(2,3,2,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedSWAveragefc3{i,1})
    title("Aberrated")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    Ratiomap = SWAveragefc3./AberratedSWAveragefc3{i,1};
    h1 = subplot(2,3,3,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),Ratiomap)
    title("Unaberrated/Aberrated")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,hot(256))
    colorbar
    clim([0,10])
    set(gca,'YDir','normal')
    axis image

    h1 = subplot(2,3,4,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),SWAveragefc3FWHM)
    title("yFWHM")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    h1 = subplot(2,3,5,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedSWAveragefc3FWHM{i,1})
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    Ratiomap = SWAveragefc3FWHM./AberratedSWAveragefc3FWHM{i,1};
    h1 = subplot(2,3,6,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),Ratiomap)
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,hot(256))
    colorbar
    clim([0,10])
    set(gca,'YDir','normal')
    axis image

    print(fig1, '-dsvg', [ SWFC3savingdir  'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(i)) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_SWFC3_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.SVG'],'-r300')
    print(fig1, '-dpng', [ SWFC3savingdir 'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(i)) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_SWFC3_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.PNG'],'-r300')
end
close all

LatticeFC3savingdir = [Analysis_savingdir 'FC3/LLS/'];
mkdir(LatticeFC3savingdir)
fig1 = figure;
for i = 1:length(RadioOrderArray)
    h1 = subplot(2,3,1,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),LatticeAveragefc3)
    title("Unaberrated")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    h1 = subplot(2,3,2,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedLatticeAveragefc3{i,1})
    title("Aberrated")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    Ratiomap = LatticeAveragefc3./AberratedLatticeAveragefc3{i,1};
    h1 = subplot(2,3,3,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),Ratiomap)
    title("Unaberrated/Aberrated")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,hot(256))
    colorbar
    clim([0,10])
    set(gca,'YDir','normal')
    axis image

    h1 = subplot(2,3,4,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),LatticeAveragefc3FWHM)
    title("yFWHM")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    h1 = subplot(2,3,5,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedLatticeAveragefc3FWHM{i,1})
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    Ratiomap = LatticeAveragefc3FWHM./AberratedLatticeAveragefc3FWHM{i,1};
    h1 = subplot(2,3,6,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),Ratiomap)
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,hot(256))
    colorbar
    clim([0,10])
    set(gca,'YDir','normal')
    axis image

    print(fig1, '-dsvg', [ LatticeFC3savingdir  'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(i)) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_LatticeFC3_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.SVG'],'-r300')
    print(fig1, '-dpng', [ LatticeFC3savingdir 'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(i)) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_LatticeFC3_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.PNG'],'-r300')
end
close all

RatioMapFC3savingdir = [Analysis_savingdir 'FC3/RatioMap/'];
mkdir(RatioMapFC3savingdir)
fig1 = figure;
for i = 1:length(RadioOrderArray)
    h1 = subplot(2,3,1,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedSWAveragefc3{i,1})
    title("SW")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    h1 = subplot(2,3,2,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedLatticeAveragefc3{i,1})
    title("LLS")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    Ratiomap = AberratedSWAveragefc3{i,1}./AberratedLatticeAveragefc3{i,1};
    h1 = subplot(2,3,3,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),Ratiomap)
    title("SW/LLS")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,hot(256))
    colorbar
    clim([0,10])
    set(gca,'YDir','normal')
    axis image

    h1 = subplot(2,3,4,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedSWAveragefc3{i,1})
    title("yFWHM")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    h1 = subplot(2,3,5,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedLatticeAveragefc3{i,1})
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    Ratiomap = AberratedSWAveragefc3{i,1}./AberratedLatticeAveragefc3{i,1};
    h1 = subplot(2,3,6,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),Ratiomap)
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,hot(256))
    colorbar
    clim([0,10])
    set(gca,'YDir','normal')
    axis image

    print(fig1, '-dsvg', [ RatioMapFC3savingdir  'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(i)) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_RatioMapFC3_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.SVG'],'-r300')
    print(fig1, '-dpng', [ RatioMapFC3savingdir 'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(i)) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_RatioMapFC3_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.PNG'],'-r300')
end
close all

%% line grating
clc
lineGratesavingdir = [Analysis_savingdir 'ConvRes/'];
mkdir(lineGratesavingdir)

for i = 1:length(RadioOrderArray)
    fig1 = figure;
    temp1 = AberratedSWconvLines{i,1};
    temp2 = AberratedLatticeconvLines{i,1};
    plot(Line_Z,SWconvLines((N+1)/2,:),'magenta','LineWidth',2)
    hold on
    plot(Line_Z,LatticeconvLines((N+1)/2,:),'Color','b','LineWidth',2,'LineStyle','-.')
    plot(Line_Z,temp1((N+1)/2,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
    plot(Line_Z,temp2((N+1)/2,:),'Color',[0.8500 0.3250 0.0980],'LineStyle','-.','LineWidth',2)
    xline(Line_Z(lineSpot),'k')
    xlim([0,20])
    legend("iSW","LLS","aberrated iSW","aberrated LLS")
    title("Spacing increment=220+" + num2str(30) + "nm")
    xlabel("z(um)")
    ylabel("Intensity (a.u)")
    hold off 
    grid on
    pbaspect([7 2 1])
    print(fig1, '-dsvg', [ lineGratesavingdir  'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(i)) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_LineGrating_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.SVG'],'-r300')
    print(fig1, '-dpng', [ lineGratesavingdir 'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(i)) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_LineGrating_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.PNG'],'-r300')

    fig2 = figure;
    imagesc(Line_Z,Line_Z,temp1)
    colormap(hot)
    hold on
    xline(Line_Z(lineSpot),'k','LineWidth',1)
    xlim([0,20])
    ylim([0,20])
    clim([0,100])
    colorbar
    title("Spacing increment=220+" + num2str(30) + "nm")
    xlabel("z(um)")
    ylabel("um")
    hold off
    print(fig2, '-dsvg', [ lineGratesavingdir  'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(i)) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_SWLineGratingFullImage_SWweight' num2str(SWweighting) '.SVG'],'-r300')
    print(fig2, '-dpng', [ lineGratesavingdir 'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(i)) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_SWLineGratingFullImage_SWweight' num2str(SWweighting)  '.PNG'],'-r300')
    
    fig3 = figure;
    imagesc(Line_Z,Line_Z,temp2)
    colormap(hot)
    hold on
    xline(Line_Z(lineSpot),'k','LineWidth',1)
    xlim([0,20])
    ylim([0,20])
    clim([0,100])
    colorbar
    title("Spacing increment=220+" + num2str(30) + "nm")
    xlabel("z(um)")
    ylabel("um")
    hold off
    print(fig3, '-dsvg', [ lineGratesavingdir  'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(i)) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_LatticeLineGratingFullImage_SWweight' num2str(SWweighting) '.SVG'],'-r300')
    print(fig3, '-dpng', [ lineGratesavingdir 'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(i)) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_LatticeLineGratingFullImage_SWweight' num2str(SWweighting)  '.PNG'],'-r300')
end
close all

%%
% figure
% plot(KZ_exc,log(zOTFexc2),'g','LineWidth',2,'LineStyle','-')
%     hold on
%     % axis image
% plot(KZ_exc,log(zOTFexc1),'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
%     xlabel("k_z/(4\pin/\lambda_{exc})")
%     ylabel("ln(OTF Strength)")
%     xticks([-0.5,-0.25,0,0.25,0.5]);
%     yticks([])
%     % yticks([0,0.5,1]);
%     xlim([-0.5,0.5])
%     ylim([-7,0])
%     grid on
% 
% figure
% plot(Z_exc,xzPSFexc2(:,(N+1)/2),'g','LineWidth',2,'LineStyle','-')
%     hold on
%     % axis image
% plot(Z_exc,xzPSFexc1(:,(N+1)/2),'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
%     xlabel("z/(\lambda_{exc}/n)")
%     ylabel("Intensity (a.u)")
%     xticks([-10,-5,0,5,10]);
%     yticks([0,0.25,0.5,0.75,1]);
%     xlim([-10,10])
%     ylim([0,1])
%     grid on