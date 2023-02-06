clear all
close all

getParameters; %modify image parameter here
CalculatePhysics;

%% aberration
Phase_factor = GetSingleZmodePupil(2,-2,4*pi);

%% detection 
% detection
PSFdet = getDetectionPSF;
PSFdet = PSFdet./(max(max(max(PSFdet))));

xzPSFdet = PSFdet(:,:,(N+1)/2);
xyPSFdet = squeeze(PSFdet((N+1)/2,:,:)); 
yzPSFdet = squeeze(PSFdet(:,(N+1)/2,:)); 
xzOTFdet = fftshift(fft2(ifftshift(xzPSFdet)));
zOTFdet = real(xzOTFdet(:,(N+1)/2));
xOTFdet = real(xzOTFdet((N+1)/2,:));

%% SW
[SWPupil,SWMask,SWPupilMetaData] = GetSWPairPupil('tophat',0.4,0.2,...
0.08,0.16,...
4/3);
AberratedSWPupil = SWPupil .* Phase_factor;

[PSFCoherent,PSFIncoherent,SWcenter] = SimulateSWPair(SWPupil);
% PrettyPlotSWPair(SWPupil,SWMask,SWPupilMetaData,PSFCoherent,PSFIncoherent,SWcenter);

[AberratedPSFCoherent,AberratedPSFIncoherent,AberratedSWcenter] = SimulateSWPair(AberratedSWPupil);
% normalize itself by putting AberratedSWcenter
% PrettyPlotSWPair(AberratedSWPupil,SWMask,SWPupilMetaData,AberratedPSFCoherent,AberratedPSFIncoherent,AberratedSWcenter);

%% Lattice
[LatticePupil,LatticeMask,LatticeMetaData] = GetLatticePupil('hex','tophat', ...
0.4,0.08, ...
0.44,0.36,...
2);
AberratedLatticePupil = LatticePupil .* Phase_factor;

[LatticePSF,LatticePSFDithered,Latticecenter] = SimulateLattice(LatticePupil);
% PrettyPlotLattice(LatticePupil,LatticeMask,LatticeMetaData,LatticePSF,LatticePSFDithered,Latticecenter); 

[AberratedLatticePSF,AberratedLatticePSFDithered,AberratedLatticecenter] =  SimulateLattice(AberratedLatticePupil); 
% PrettyPlotLattice(AberratedLatticePupil,LatticeMask,LatticeMetaData,AberratedLatticePSF,AberratedLatticePSFDithered,AberratedLatticecenter);

%% Excitation
% excitation
PSFIncoherent2 = PSFIncoherent/max(PSFIncoherent,[],'all');
AberratedPSFIncoherent2 = AberratedPSFIncoherent/SWcenter(2,1); % normalized to unaberrated center
% AberratedPSFIncoherent2 = AberratedPSFIncoherent/max(max(max(AberratedPSFIncoherent))); % normalized to itself

LatticePSFDithered2 = LatticePSFDithered/max(LatticePSFDithered,[],'all');
AberratedLatticePSFDithered2 = AberratedLatticePSFDithered/Latticecenter(2,1); % normalized to unaberrated center
% AberratedLatticePSFDithered2 = AberratedLatticePSFDithered/max(max(max(AberratedLatticePSFDithered))); % normalized to itself

PSFexc = AberratedLatticePSFDithered2;
xzPSFexc = PSFexc(:,:,(N+1)/2); 
xyPSFexc = squeeze(PSFexc((N+1)/2,:,:)); 
yzPSFexc = squeeze(PSFexc(:,(N+1)/2,:)); 
xzOTFexc = fftshift(fft2(ifftshift(xzPSFexc)));
zOTFexc = real(xzOTFexc(:,(N+1)/2));
xOTFexc = real(xzOTFexc((N+1)/2,:));

% overall
PSFoverall = PSFexc .* PSFdet;
PSFoverall = PSFoverall./(max(max(max(PSFoverall))));
OTFoverall = fftshift(fftn(ifftshift(PSFoverall)));
OTFoverall = OTFoverall./(max(max(max(OTFoverall))));

xzPSFoverall = PSFoverall(:,:,(N+1)/2); 
xyPSFoverall = squeeze(PSFoverall((N+1)/2,:,:)); 
yzPSFoverall = squeeze(PSFoverall(:,(N+1)/2,:)); 
xzOTFoverall = OTFoverall(:,:,(N+1)/2); 
zOTFoverall = real(xzOTFoverall(:,(N+1)/2));
xOTFoverall = real(xzOTFoverall((N+1)/2,:));

%% beads 2d
SNR = 25;
[SWBeads1,SWBeads2] = beadsSimulation(PSFIncoherent,PSFdet,SNR);
[AberratedSWBeads1,AberratedSWBeads2] = beadsSimulation(AberratedPSFIncoherent,PSFdet,SNR);

[LatticeBeads1,LatticeBeads2] = beadsSimulation(LatticePSFDithered,PSFdet,SNR);
[AberratedLatticeBeads1,AberratedLatticeBeads2] = beadsSimulation(AberratedLatticePSFDithered,PSFdet,SNR);

%% FC2
fc2SW = FC2(SWBeads1,SWBeads2);
fc2 = FC2(LatticeBeads1,LatticeBeads2);

%% beads 3d
[SWBeadsVol1,SWBeadsVol2] = beadsSimulation(AberratedPSFIncoherent);
[LatticeBeadsVol1,LatticeBeadsVol2] = beadsSimulation(AberratedLatticePSFDithered);

%% FC3
fc3SW = FC2(SWBeadsVol1,SWBeadsVol2);
fc3Lattice = FC2(LatticeBeadsVol1,LatticeBeadsVol2);

%% Iterate FC2
SNR = [5];
iteration = [100];
for j = 1:length(SNR)
    SNRloop = SNR(j)
    Averagefc2 = zeros((N+1)/2,(N+1)/2);
    for i = 1:length(iteration)
        iterationloop = iteration(i)
        for k = 1:iterationloop
            [Beads1,Beads2] = beadsSimulation(PSFexc,PSFdet,SNRloop);
            fc2 = FC2(Beads1,Beads2);
            Averagefc2 = Averagefc2 + fc2;
        end
        Averagefc2 = Averagefc2/iterationloop;
        Averagefc2 = Averagefc2/max(max(Averagefc2));

        fig1 = figure;
        fig1.Name = "SNR=" + num2str(SNRloop) + ",Iteration=" + num2str(iterationloop);
        fig1.WindowState = 'maximized';
    
        h1 = subplot(1,3,1);
        imagesc(KX_exc((N+1)/2:N),flip(KZ_exc((N+1)/2:N)),Averagefc2)
        xlabel("k_x/(4\pin/\lambda_{exc})")
        ylabel("k_z/(4\pin/\lambda_{exc})")
        colormap(hot)
        colorbar
        set(h1, 'YDir','normal')
        axis image
        caxis([0,1]);
    
        h2 = subplot(1,3,2);
        fc2contour = zeros(size(Averagefc2));
        fc2contour(Averagefc2>=1/7) = 1;
        imagesc(KX_exc((N+1)/2:N),flip(KZ_exc((N+1)/2:N)),fc2contour)
        xlabel("k_x/(4\pin/\lambda_{exc})")
        ylabel("k_z/(4\pin/\lambda_{exc})")
        title("1/7")
        set(h2, 'YDir','normal')
        axis image
        caxis([0,1]);
        colorbar
    
        subplot(1,3,3)
        plot(KX_exc((N+1)/2:N),Averagefc2(end,:))
        hold on
        grid on
        axis image
        ylim([-0.3,1])
        xlim([0 1])
        plot(flip(KZ_exc((N+1)/2:N)),Averagefc2(:,1))
        plot(KZ_exc((N+1)/2:N),zOTFoverall((N+1)/2:end));
        plot(KX_exc((N+1)/2:N),xOTFoverall((N+1)/2:end))
        yline(1/7)
        legend("KX","KZ","zOTF","xOTF",'Location','northeast')
        xlabel("k/(4\pin/\lambda_{exc})")    

%         % overall OTF
%         subplot(2,2,4)
%         imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),abs(xzOTFoverall(1:(N+1)/2,(N+1)/2:end)));
%         xlabel("k_x/(4\pin/\lambda_{exc})")
%         ylabel("k_z/(4\pin/\lambda_{exc})")
%         title("Overall OTF")
%         colorbar
%         axis image
%         caxis([0,0.1])
    end
end

    
