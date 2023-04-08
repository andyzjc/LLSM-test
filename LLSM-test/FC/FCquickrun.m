clear all
close all

getParameters; %modify image parameter here
CalculatePhysics;

%% aberration
Phase_factor = GetSingleZmodePupil(2,0,2*pi);

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
[SWPupil,SWMask,SWPupilMetaData] = GetSWPairPupil('gaussian',0,0.21,...
0,0.16,...
7/10); 
AberratedSWPupil = SWPupil .* Phase_factor;

[PSFCoherent,PSFIncoherent,SWcenter] = SimulateSWPair(SWPupil);
% PrettyPlotSWPair(SWPupil,SWMask,SWPupilMetaData,PSFCoherent,PSFIncoherent,SWcenter);

[AberratedPSFCoherent,AberratedPSFIncoherent,AberratedSWcenter] = SimulateSWPair(AberratedSWPupil);
% normalize itself by putting AberratedSWcenter
% PrettyPlotSWPair(AberratedSWPupil,SWMask,SWPupilMetaData,AberratedPSFCoherent,AberratedPSFIncoherent,AberratedSWcenter);

%% Lattice
[LatticePupil,LatticeMask,LatticeMetaData] = GetLatticePupil('square','gaussian', ...
0.4,0.08, ...
0.6,0.2,...
1);
%AberratedLatticePupil = LatticePupil .* Phase_factor;

[LatticePSF,LatticePSFDithered,Latticecenter] = SimulateLattice(LatticePupil);
% PrettyPlotLattice(LatticePupil,LatticeMask,LatticeMetaData,LatticePSF,LatticePSFDithered,Latticecenter); 

%[AberratedLatticePSF,AberratedLatticePSFDithered,AberratedLatticecenter] =  SimulateLattice(AberratedLatticePupil); 
% PrettyPlotLattice(AberratedLatticePupil,LatticeMask,LatticeMetaData,AberratedLatticePSF,AberratedLatticePSFDithered,AberratedLatticecenter);

%% Field Synthesis 
[FSPupil,FSMask,FSMetaData] = GetFSPupil('hex','tophat', ...
0.56,0.04, ...
0.42,0.38,...
2);

[FSPSF,FScenter] = SimulateFS(FSPupil);
% PrettyPlotFS(FSPupil,FSMask,FSMetaData,FSPSF,FScenter); 
FSPSF = FSPSF/max(FSPSF,[],'all');

%% Excitation
% excitation
PSFIncoherent2 = PSFIncoherent/max(PSFIncoherent,[],'all');
% AberratedPSFIncoherent2 = AberratedPSFIncoherent/SWcenter(2,1); % normalized to unaberrated center
AberratedPSFIncoherent2 = AberratedPSFIncoherent/max(max(max(AberratedPSFIncoherent))); % normalized to itself

LatticePSFDithered2 = LatticePSFDithered/max(LatticePSFDithered,[],'all');
LatticePSF2 = LatticePSF/max(LatticePSF,[],'all');
% AberratedLatticePSFDithered2 = AberratedLatticePSFDithered/Latticecenter(2,1); % normalized to unaberrated center
AberratedLatticePSFDithered2 = AberratedLatticePSFDithered/max(max(max(AberratedLatticePSFDithered))); % normalized to itself

%% Overall
PSFexc = PSFIncoherent;
PSFexc = PSFexc/max(PSFexc,[],'all');
xzPSFexc = PSFexc(:,:,(N+1)/2); 
xyPSFexc = squeeze(PSFexc((N+1)/2,:,:)); 
yzPSFexc = squeeze(PSFexc(:,(N+1)/2,:)); 
xzOTFexc = fftshift(fft2(ifftshift(xzPSFexc)));
zOTFexc = real(xzOTFexc(:,(N+1)/2)); zOTFexc = zOTFexc/max(zOTFexc);
xOTFexc = real(xzOTFexc((N+1)/2,:));

% yindex = 1-(xyPSFexc((N+1)/2,:) <= 0.5*max(xyPSFexc((N+1)/2,:)));
% yFWHM1 = find(yindex,1,'first') ;
% yFWHM2 = find(yindex,1,'last');
% 
% xindex = 1-(xyPSFexc(:,(N+1)/2) <= 0.5*max(xyPSFexc(:,(N+1)/2)));
% xFWHM1 = find(xindex,1,'first') ;
% xFWHM2 = find(xindex,1,'last');

% figure(1)
% imagesc(X_exc,Y_exc,squeeze(max(PSFexc,[],1)))
% xlabel("y/(\lambda_{exc}/n)")
% ylabel("x/(\lambda_{exc}/n)")
% colorbar
% colormap(hot)
% hold on
% rectangle('Position',[Y_exc(yFWHM1) X_exc(xFWHM1) Y_exc(yFWHM2)-Y_exc(yFWHM1) X_exc(xFWHM2)-X_exc(xFWHM1)],'EdgeColor','g','LineWidth',2)

% overall
 PSFoverall = PSFexc .* PSFdet;
% PSFoverallNOISE = PSFoverall + poissrnd(PSFoverall) .* 1/20; 
% PSFoverallNOISE = fillmissing(PSFoverallNOISE,'constant',0);
% 
% PSFdecon = LatticePSFDithered2.*PSFdet;
% PSFdecon = PSFdecon + poissrnd(PSFdecon) .* 1/20; 
% PSFdecon = fillmissing(PSFdecon,'constant',0);
% PSFoverall = deconvlucy(PSFoverallNOISE,PSFdecon,20);
% PSFoverall = PSFoverallNOISE;

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
[SWBeads1,SWBeads2] = beadsSimulation2d(PSFIncoherent,PSFdet,SNR);
[AberratedSWBeads1,AberratedSWBeads2] = beadsSimulation2d(AberratedPSFIncoherent,PSFdet,SNR);

[LatticeBeads1,LatticeBeads2] = beadsSimulation2d(LatticePSFDithered,PSFdet,SNR);
[AberratedLatticeBeads1,AberratedLatticeBeads2] = beadsSimulation2d(AberratedLatticePSFDithered,PSFdet,SNR);

%% FC2
fc2SW = FC2(SWBeads1,SWBeads2);
fc2 = FC2(LatticeBeads1,LatticeBeads2);

%% beads 3d
[SWBeadsVol1,SWBeadsVol2,~] = beadsSimulation3d(AberratedPSFIncoherent,PSFdet,SNR);
[LatticeBeadsVol1,LatticeBeadsVol2] = beadsSimulation3d(AberratedLatticePSFDithered,PSFdet,SNR);

%% FC3
fc3SW = FC3(SWBeadsVol1,SWBeadsVol2);
fc3Lattice = FC3(LatticeBeadsVol1,LatticeBeadsVol2);

%% Iterate FC2
% testOTF = real(fftshift(fft2(ifftshift(xzPSFoverall))));
% testOTF = testOTF/max(max(testOTF));

SNR = [100];
iteration = [500];
for j = 1:length(SNR)
    SNRloop = SNR(j)
    Averagefc2 = zeros((N+1)/2,(N+1)/2);
    for i = 1:length(iteration)
        iterationloop = iteration(i)
        for k = 1:iterationloop
            test1 = [];
            test2 = [];
%              test1 = xzPSFoverall + poissrnd(xzPSFoverall) * 1/SNR;
%              test2 = xzPSFoverall + poissrnd(xzPSFoverall) * 1/SNR; 
%              test1 = fillmissing(test1,'constant',0);
%              test2 = fillmissing(test2,'constant',0);
%              fc2 = FC2(test1,test2);
            [Beads1,Beads2] = beadsSimulation2d(PSFexc,PSFdet,SNRloop);
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
        plot(KZ_exc((N+1)/2:N),testOTF((N+1)/2:end,(N+1)/2))
        plot(KX_exc((N+1)/2:N),testOTF((N+1)/2,(N+1)/2:end))
        yline(1/7)
        legend("KX","KZ","correct zOTF","correct xOTF","Incorrect zOTF","Incorrect xOTF",'Location','northeast')
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

%% Iterate FC3

allPSFexc = cell(4,1);
allPSFexc{1,1} = LatticePSFDithered2;
allPSFexc{2,1} = AberratedLatticePSFDithered2;
allPSFexc{3,1} = PSFIncoherent2;
allPSFexc{4,1} = AberratedPSFIncoherent2;
FWHMIndex = [346,346,354,354];

SNR = 20;
iteration = 10;
for dd = 1:length(allPSFexc)
    PSFexc = allPSFexc{dd,1};
    for j = 1:length(SNR)
        SNRloop = SNR(j)
        Averagefc3 = zeros((N+1)/2,(N+1)/2);
        for i = 1:length(iteration)
            iterationloop = iteration(i)
            for k = 1:iterationloop
                [BeadsVol1,BeadsVol2,~] = beadsSimulation3d_focal(PSFexc,PSFdet,SNRloop);
                [FRCGraph,~,~,KR,~,~] = FC3(BeadsVol1,BeadsVol2);
                Averagefc3 = Averagefc3 + FRCGraph;
            end
            Averagefc3 = Averagefc3/iterationloop;
            Averagefc3 = Averagefc3/max(max(Averagefc3));
        end
    %         % lateral cutoff (kr)
    %         for aa = 1:size(Averagefc3,1)
    %             line = Averagefc3(aa,:);
    %             FRCcutoffLateral(aa,1) = min(KR(abs(line)<1/7));
    %         end
    %         
    %         % axial cutoff (kz)
             KZ = KZ_exc((N+1)/2:end);
    %         for bb = 1:size(Averagefc3,2)
    %                 line = Averagefc3(:,bb);
    %             if isempty(min(KZ(abs(line)<1/7)))
    %                 FRCcutoffAxial(bb,1) = nan;
    %             else
    %                 FRCcutoffAxial(bb,1) = min(KZ(abs(line)<1/7));
    %             end
    %         end
    
            % find max kz
            kr_line = Averagefc3(1,:);
            kr_Index = find((abs(kr_line)<=1/7),1);
            Rboundary = KR(kr_Index);
            if isempty(Rboundary)
                [~,kr_Index] = min(kr_line);
                Rboundary = KR(kr_Index);
            end
    
            halfIndex = round(kr_Index/2);
            kz_line = Averagefc3(:,halfIndex);
            kz_Index = find((abs(kz_line)<=1/7),1);
            Zboundary = KZ(kz_Index);
            if isempty(Zboundary)
                [~,kz_Index] = min(kr_line);
                Zboundary = KZ(kz_Index);
            end
    
            fig1 = figure;
            fig1.Name = "SNR=" + num2str(SNRloop) + ",Iteration=" + num2str(iterationloop);
            hold on
    %         fig1.WindowState = 'maximized';
        
    %         h1 = subplot(1,3,1);
            imagesc(KR,KZ_exc((N+1)/2:N),Averagefc3)
            xlabel("k_r/(4\pin/\lambda_{exc})")
            ylabel("k_z/(4\pin/\lambda_{exc})")
            title("Kr Boundary=" + num2str(Rboundary)+", Max Kz=" + num2str(Zboundary))
            colorbar
            axis image
            yline(Zboundary,'--','Color','r','LineWidth',2)
            xline(KR(halfIndex),'--','Color','r','LineWidth',2)
            xline(Rboundary,'--','Color','r','LineWidth',2)
        
    %         h2 = subplot(1,3,2); 
    %         plot(KZ,FRCcutoffLateral,'Color','r','LineWidth',2)
    %         grid on
    %         xlabel("k_z/(4\pin/\lambda_{exc})")
    %         ylabel("Cutoff k_r/(4\pin/\lambda_{exc})")
    %         yline(1/7)
    %         axis image
    %         ylim([0,1])
    %         xlim([0,1])
    %     
    %         h3 = subplot(1,3,3);
    %         plot(KR,FRCcutoffAxial,'Color','r','LineWidth',2)
    %         grid on
    %         xlabel("k_r/(4\pin/\lambda_{exc})")
    %         ylabel("Cutoff k_z/(4\pin/\lambda_{exc})")
    %         yline(1/7) 
    %         drawnow
    %         axis image
    %         ylim([0,1])
    %         xlim([0,1])
    end
end

%% Iterate FC3 FWHM
for dd = 1:length(allPSFexc)
    PSFexc = allPSFexc{dd,1};
    for j = 1:length(SNR)
        SNRloop = SNR(j)
        Averagefc3 = zeros((N+1)/2,(N+1)/2);
        for i = 1:length(iteration)
            iterationloop = iteration(i)
            for k = 1:iterationloop
                [BeadsVol1,BeadsVol2,~] = beadsSimulation3d_FWHM(PSFexc,PSFdet,SNRloop,FWHMIndex(dd));
                [FRCGraph,~,~,KR,~,~] = FC3(BeadsVol1,BeadsVol2);
                Averagefc3 = Averagefc3 + FRCGraph;
            end
            Averagefc3 = Averagefc3/iterationloop;
            Averagefc3 = Averagefc3/max(max(Averagefc3));
        end
    %         % lateral cutoff (kr)
    %         for aa = 1:size(Averagefc3,1)
    %             line = Averagefc3(aa,:);
    %             FRCcutoffLateral(aa,1) = min(KR(abs(line)<1/7));
    %         end
    %         
    %         % axial cutoff (kz)
             KZ = KZ_exc((N+1)/2:end);
    %         for bb = 1:size(Averagefc3,2)
    %                 line = Averagefc3(:,bb);
    %             if isempty(min(KZ(abs(line)<1/7)))
    %                 FRCcutoffAxial(bb,1) = nan;
    %             else
    %                 FRCcutoffAxial(bb,1) = min(KZ(abs(line)<1/7));
    %             end
    %         end
    
            % find max kz
            kr_line = Averagefc3(1,:);
            kr_Index = find((abs(kr_line)<=1/7),1);
            Rboundary = KR(kr_Index);
            if isempty(Rboundary)
                [~,kr_Index] = min(kr_line);
                Rboundary = KR(kr_Index);
            end
    
            halfIndex = round(kr_Index/2);
            kz_line = Averagefc3(:,halfIndex);
            kz_Index = find((abs(kz_line)<=1/7),1);
            Zboundary = KZ(kz_Index);
            if isempty(Zboundary)
                [~,kz_Index] = min(kr_line);
                Zboundary = KZ(kr_Index);
            end
    
            fig1 = figure;
            fig1.Name = "SNR=" + num2str(SNRloop) + ",Iteration=" + num2str(iterationloop);
            hold on
    %         fig1.WindowState = 'maximized';
        
    %         h1 = subplot(1,3,1);
            imagesc(KR,KZ_exc((N+1)/2:N),Averagefc3)
            xlabel("k_r/(4\pin/\lambda_{exc})")
            ylabel("k_z/(4\pin/\lambda_{exc})")
            title("Kr Boundary=" + num2str(Rboundary)+", Max Kz=" + num2str(Zboundary))
            colorbar
            axis image
            yline(Zboundary,'--','Color','r','LineWidth',2)
            xline(KR(halfIndex),'--','Color','r','LineWidth',2)
            xline(Rboundary,'--','Color','r','LineWidth',2)
        
    %         h2 = subplot(1,3,2); 
    %         plot(KZ,FRCcutoffLateral,'Color','r','LineWidth',2)
    %         grid on
    %         xlabel("k_z/(4\pin/\lambda_{exc})")
    %         ylabel("Cutoff k_r/(4\pin/\lambda_{exc})")
    %         yline(1/7)
    %         axis image
    %         ylim([0,1])
    %         xlim([0,1])
    %     
    %         h3 = subplot(1,3,3);
    %         plot(KR,FRCcutoffAxial,'Color','r','LineWidth',2)
    %         grid on
    %         xlabel("k_r/(4\pin/\lambda_{exc})")
    %         ylabel("Cutoff k_z/(4\pin/\lambda_{exc})")
    %         yline(1/7) 
    %         drawnow
    %         axis image
    %         ylim([0,1])
    %         xlim([0,1])
    end
end
