clear all
close all

getParameters; %modify image parameter here
CalculatePhysics;

%% Lattice
NA= [0.4,0.56];
deltaNA = [0.04,0.08];
counter = 1;
for i = 1:length(NA)
    for j = 1:length(deltaNA)
        NA(i)
        deltaNA(j)
        [LatticePupil,LatticeMask,LatticeMetaData] = GetLatticePupil('hex','tophat', ...
        NA(i),deltaNA(j), ...
        0.60,0.52,...
        1);
        [LatticePSF,~,Latticecenter] = SimulateLattice(LatticePupil);
        LatticePSF2 = LatticePSF/max(LatticePSF,[],'all');

        PSFexc = LatticePSF2;
        xyPSFexc = squeeze(max(PSFexc,[],1));
        yindex = 1-(xyPSFexc((N+1)/2,:) <= 0.5*max(xyPSFexc((N+1)/2,:)));
        yFWHM1 = find(yindex,1,'first') ;
        yFWHM2 = find(yindex,1,'last');
        
        xindex = 1-(xyPSFexc(:,(N+1)/2) <= 0.5*max(xyPSFexc(:,(N+1)/2)));
        xFWHM1 = find(xindex,1,'first') ;
        xFWHM2 = find(xindex,1,'last');
        
        line1 = flip(xyPSFexc(1:(N+1)/2,yFWHM1));
        xindex1 = 1-(line1 <= 0.5*max(xyPSFexc(:,yFWHM1)));
        line2 = xyPSFexc((N+1)/2:end,yFWHM1);
        xindex2 = 1-(line2 <= 0.5*max(xyPSFexc(:,yFWHM1)));
        xFWHM1_FWHM= (N+1)/2-find(abs(diff(xindex1)),1,'last') ;
        xFWHM2_FWHM = find(abs(diff(xindex1)),1,'last') + (N+1)/2;
        
        figure(counter)
        title("NA="+num2str(NA(i))+", deltaNA="+num2str(deltaNA(j)) )
        imagesc(X_exc,Y_exc,squeeze(max(PSFexc,[],1)))
        xlabel("y/(\lambda_{exc}/n)")
        ylabel("x/(\lambda_{exc}/n)")
        colorbar
        colormap(hot)
        hold on
        rectangle('Position',[Y_exc(yFWHM1) X_exc(xFWHM1) Y_exc(yFWHM2)-Y_exc(yFWHM1) X_exc(xFWHM2)-X_exc(xFWHM1)],'EdgeColor','g','LineWidth',2)
        rectangle('Position',[Y_exc(yFWHM1) X_exc(xFWHM1_FWHM) Y_exc(yFWHM2)-Y_exc(yFWHM1) X_exc(xFWHM2_FWHM)-X_exc(xFWHM1_FWHM)],'EdgeColor','r','LineWidth',2)
        hold off
        counter = counter + 1;
    end
end

%% SW
NA = [0.4];
deltaNA = [0.04];
counter = 1;
for i = 1:length(NA1)
    for j = 1:length(deltaNA)
        NA(i)
        deltaNA(j)
        [SWPupil,SWMask,SWPupilMetaData] = GetSWPairPupil('tophat',NA1(i),0,...
        deltaNA(j),0,...
        1);
        
        [~,PSFIncoherent,~] = SimulateSWPair(SWPupil);
        PSFIncoherent = PSFIncoherent/max(PSFIncoherent,[],'all');

        PSFexc = PSFIncoherent;
        xyPSFexc = squeeze(sum(PSFexc,1));
        yindex = 1-(xyPSFexc((N+1)/2,:) <= 0.5*max(xyPSFexc((N+1)/2,:)));
        yFWHM1 = find(yindex,1,'first') ;
        yFWHM2 = find(yindex,1,'last');
        
        xindex = 1-(xyPSFexc(:,(N+1)/2) <= 0.5*max(xyPSFexc(:,(N+1)/2)));
        xFWHM1 = find(xindex,1,'first') ;
        xFWHM2 = find(xindex,1,'last');
        
%         xindex2 = 1-(xyPSFexc(:,yFWHM1) <= 0.5*max(xyPSFexc(:,(N+1)/2)));
%         xFWHM1_FWHM= find(xindex2,1,'first') ;
%         xFWHM2_FWHM = find(xindex2,1,'last');

        line1 = flip(xyPSFexc(1:(N+1)/2,yFWHM1));
        xindex1 = 1-(line1 <= 0.5*max(xyPSFexc(:,yFWHM1)));
        line2 = xyPSFexc((N+1)/2:end,yFWHM1);
        xindex2 = 1-(line2 <= 0.5*max(xyPSFexc(:,yFWHM1)));
        xFWHM1_FWHM= (N+1)/2-find(abs(diff(xindex1)),1,'first') ;
        xFWHM2_FWHM = find(abs(diff(xindex1)),1,'first') + (N+1)/2;
        
        figure(counter)
        title("NA="+num2str(NA(i))+", deltaNA="+num2str(deltaNA(j)) )
        imagesc(X_exc,Y_exc,squeeze(sum(PSFexc,1)))
%         imagesc(X_exc,Y_exc,squeeze(PSFexc((N+1)/2,:,:)))
        xlabel("y/(\lambda_{exc}/n)")
        ylabel("x/(\lambda_{exc}/n)")
        colorbar
        colormap(hot)
        hold on
        rectangle('Position',[Y_exc(yFWHM1) X_exc(xFWHM1) Y_exc(yFWHM2)-Y_exc(yFWHM1) X_exc(xFWHM2)-X_exc(xFWHM1)],'EdgeColor','g','LineWidth',2)
        rectangle('Position',[Y_exc(yFWHM1) X_exc(xFWHM1_FWHM) Y_exc(yFWHM2)-Y_exc(yFWHM1) X_exc(xFWHM2_FWHM)-X_exc(xFWHM1_FWHM)],'EdgeColor','r','LineWidth',2)
        hold off
        counter = counter + 1;
    end
end

%% FS
NA1 = [0.4,0.56];
deltaNA = [0.04,0.08];

counter = 1;
for i = 1:length(NA)
    for j = 1:length(deltaNA)
        NA(i)
        deltaNA(j)
        [FSPupil,FSMask,FSMetaData] = GetFSPupil('hex','tophat', ...
        NA(i),deltaNA(j), ...
        0.60,0.52,...
        1);
        
        [FSPSF,~] = SimulateFS(FSPupil);
        % PrettyPlotFS(FSPupil,FSMask,FSMetaData,FSPSF,FScenter); 
        FSPSF = FSPSF/max(FSPSF,[],'all');

        PSFexc = FSPSF;
        xyPSFexc = squeeze(max(PSFexc,[],1));
        yindex = 1-(xyPSFexc((N+1)/2,:) <= 0.5*max(xyPSFexc((N+1)/2,:)));
        yFWHM1 = find(yindex,1,'first') ;
        yFWHM2 = find(yindex,1,'last');
        
        xindex = 1-(xyPSFexc(:,(N+1)/2) <= 0.5*max(xyPSFexc(:,(N+1)/2)));
        xFWHM1 = find(xindex,1,'first') ;
        xFWHM2 = find(xindex,1,'last');
        
%         xindex2 = 1-(xyPSFexc(:,yFWHM1) <= 0.5*max(xyPSFexc(:,(N+1)/2)));
%         xindex2 = 1-(xyPSFexc(:,yFWHM1) <= 0.5*xyPSFexc((N+1)/2,yFWHM1));

        line1 = flip(xyPSFexc(1:(N+1)/2,yFWHM1));
        xindex1 = 1-(line1 <= 0.5*max(xyPSFexc(:,yFWHM1)));
        line2 = xyPSFexc((N+1)/2:end,yFWHM1);
        xindex2 = 1-(line2 <= 0.5*max(xyPSFexc(:,yFWHM1)));
        xFWHM1_FWHM= (N+1)/2-find(abs(diff(xindex1)),1,'first') ;
        xFWHM2_FWHM = find(abs(diff(xindex1)),1,'first') + (N+1)/2;
        
        figure(counter)
        title("NA="+num2str(NA(i))+", deltaNA="+num2str(deltaNA(j)) )
        imagesc(X_exc,Y_exc,squeeze(max(PSFexc,[],1)))
        xlabel("y/(\lambda_{exc}/n)")
        ylabel("x/(\lambda_{exc}/n)")
        colorbar
        colormap(hot)
        hold on
        rectangle('Position',[Y_exc(yFWHM1) X_exc(xFWHM1) Y_exc(yFWHM2)-Y_exc(yFWHM1) X_exc(xFWHM2)-X_exc(xFWHM1)],'EdgeColor','g','LineWidth',2)
        rectangle('Position',[Y_exc(yFWHM1) X_exc(xFWHM1_FWHM) Y_exc(yFWHM2)-Y_exc(yFWHM1) X_exc(xFWHM2_FWHM)-X_exc(xFWHM1_FWHM)],'EdgeColor','r','LineWidth',2)
        hold off
        counter = counter + 1;
    end
end
