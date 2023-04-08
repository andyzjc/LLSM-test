function [FRCGraph,FRCUpperGraph,FRCLowerGraph,KR,FRCcutoffLateral,FRCcutoffAxial] = FC3(Vol1,Vol2)
    getParameters; %modify image parameter here
    CalculatePhysics;

    % get rid of Nan
    Vol1 = fillmissing(Vol1,'constant',0);
    Vol2 = fillmissing(Vol2,'constant',0);

    % fourier transform 
    FTVol1 = fftshift(fftn(ifftshift(Vol1)));
    FTVol2 = fftshift(fftn(ifftshift(Vol2)));    

    % loop through kz, and perform FRC
    FRCGraph = zeros((N+1)/2,(N+1)/2);
    FRCUpperGraph = FRCGraph;
    FRCLowerGraph = FRCGraph;
    for i = 0:(N+1)/2-1
        % positve kz
        kXY1upper = squeeze(FTVol1((N+1)/2+i,:,:)); 
        kXY2upper = squeeze(FTVol2((N+1)/2+i,:,:));
        [FRCValueUpper,KR,~] = FRC(kXY1upper,kXY2upper,deltax);

        % negative kz
        kXY1lower = squeeze(FTVol1((N+1)/2-i,:,:));
        kXY2lower = squeeze(FTVol2((N+1)/2-i,:,:));
        [FRCValueLower,~,~] = FRC(kXY1lower,kXY2lower,deltax);

        % average 
        % FRCUpperGraph(i+1,:) = abs(FRCValueUpper)';
        % FRCLowerGraph(i+1,:) = abs(FRCValueLower)';
        FRCGraph(i+1,:) = ((abs(FRCValueUpper)+abs(FRCValueLower))./2)'; 

    end

    KR = KR./(2*k_wave);
%     % lateral cutoff (kr)
%     for i = 1:size(FRCGraph,1)
%         line = FRCGraph(i,:);
%         FRCcutoffLateral(i,1) = min(KR(abs(line)<1/7));
%     end
% 
%     % axial cutoff (kz)
%     KZ = KZ_exc((N+1)/2:end);
%     for i = 1:size(FRCGraph,2)
%         line = FRCGraph(:,i);
%         if isempty(min(KZ(abs(line)<1/7)))
%             FRCcutoffAxial(i,1) = nan;
%         else
%             FRCcutoffAxial(i,1) = min(KZ(abs(line)<1/7));
%         end
%         
%     end
    FRCcutoffLateral = 0;
    FRCcutoffAxial = 0;
    FRCUpperGraph = 0;
    FRCLowerGraph = 0;
%         h1 = subplot(1,3,1);
%         imagesc(KR,KZ_exc((N+1)/2:N),FRCGraph)
%         xlabel("k_r/(4\pin/\lambda_{exc})")
%         ylabel("k_z/(4\pin/\lambda_{exc})")
%         colorbar
%         axis image
%     
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