function [Averagefc3,Averagefc3FWHM,KR,KZ] = RunFC3(PSFexc,PSFdet,OTFmask,Iter,SNR,yFWHM)
    getParameters; %modify image parameter here
    CalculatePhysics;

    Averagefc3 = zeros((N+1)/2,(N+1)/2);
    Averagefc3FWHM = Averagefc3;

    for i = 1:length(Iter)
        [BeadsVol1,BeadsVol2,GT] = beadsSimulation3d_plane(PSFexc,PSFdet,SNR,(N+1)/2); %focal plane 
        [FRCGraph,~,~,KR,~,~] = FC3(BeadsVol1,BeadsVol2);
        Averagefc3 = Averagefc3 + FRCGraph;
    end
    Averagefc3 = Averagefc3/Iter;
    Averagefc3 = Averagefc3/max(Averagefc3,[],'all');
    Averagefc3 = Averagefc3 .* OTFmask;

    for i = 1:length(Iter)
        [BeadsVol1,BeadsVol2,~] = beadsSimulation3d_plane(PSFexc,PSFdet,SNR,yFWHM);
        [FRCGraph,~,~,~,~,~] = FC3(BeadsVol1,BeadsVol2);
        Averagefc3FWHM = Averagefc3FWHM + FRCGraph;
    end
    Averagefc3FWHM = Averagefc3FWHM/Iter;
    Averagefc3FWHM = Averagefc3FWHM/max(Averagefc3FWHM,[],'all');
    Averagefc3FWHM = Averagefc3FWHM .* OTFmask;

    % find max kz
    KZ = KZ_exc((N+1)/2:end);
    % kr_line = Averagefc3(1,:);
    % kr_Index = find((abs(kr_line)<=1/7),1);
    % Rboundary = KR(kr_Index);
    % if isempty(Rboundary)
    %     [~,kr_Index] = min(kr_line);
    %     Rboundary = KR(kr_Index);
    % end
    % 
    % halfIndex = round(kr_Index/2);
    % kz_line = Averagefc3(:,halfIndex);
    % kz_Index = find((abs(kz_line)<=1/7),1);
    % Zboundary = KZ(kz_Index);
    % if isempty(Zboundary)
    %     [~,kz_Index] = min(kr_line);
    %     Zboundary = KZ(kz_Index);
    % end

    % fig1 = figure;
    % fig1.Name = "SNR=" + num2str(SNRloop) + ",Iteration=" + num2str(iterationloop);
    % hold on
    % imagesc(KR,KZ_exc((N+1)/2:N),Averagefc3)
    % xlabel("k_r/(4\pin/\lambda_{exc})")
    % ylabel("k_z/(4\pin/\lambda_{exc})")
    % title("Kr Boundary=" + num2str(Rboundary)+", Max Kz=" + num2str(Zboundary))
    % colorbar
    % axis image
    % yline(Zboundary,'--','Color','r','LineWidth',2)
    % xline(KR(halfIndex),'--','Color','r','LineWidth',2)
    % xline(Rboundary,'--','Color','r','LineWidth',2)
    