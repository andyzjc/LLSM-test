function grapthSW(NA1,deltaNA,LatticeType,weighting,PSFIncoherent,PSFCoherent,SWPupil,PSFdet,savingdir)
    getParameters; %modify image parameter here
    CalculatePhysics;
    firemap = fire(256);
    
    Sampleoffset = 20;
    Koffset = 0.5;

    % PSFIncoherent = PSFCoherent;

    IncoherentPSFexc = PSFIncoherent;
    % IncoherentPSFexc = IncoherentPSFexc/max(IncoherentPSFexc,[],'all');
    IncoherentxzPSFexc = IncoherentPSFexc(:,:,(N+1)/2); 
    IncoherentyzPSFexc = squeeze(IncoherentPSFexc(:,(N+1)/2,:)); 
    IncoherentxyPSFexc = squeeze(IncoherentPSFexc((N+1)/2,:,:)); 

    % IncoherentxyPSFexcDithered = zeros(size(IncoherentxyPSFexc));
    % for j = -5:5
    %     IncoherentxyPSFexcDithered = IncoherentxyPSFexcDithered + ...
    %         circshift(IncoherentxyPSFexc,j,1);
    % end
    % IncoherentxyPSFexc = IncoherentxyPSFexcDithered;
    % IncoherentxyPSFexc = IncoherentxyPSFexc/max(IncoherentxyPSFexc,[],'all');

    IncoherentzPSFexc = IncoherentyzPSFexc(:,(N+1)/2);
    IncoherentxzOTFexc = fftshift(fft2(ifftshift(IncoherentxzPSFexc))); 
    IncoherentzOTFexc = IncoherentxzOTFexc(:,(N+1)/2); 

    CoherentPSFexc = PSFCoherent;
    CoherentxzPSFexc = CoherentPSFexc(:,:,(N+1)/2); 
    CoherentyzPSFexc = squeeze(CoherentPSFexc(:,(N+1)/2,:)); 
    CoherentxyPSFexc = squeeze(CoherentPSFexc((N+1)/2,:,:)); 
    CoherentzPSFexc = CoherentyzPSFexc(:,(N+1)/2);
    CoherentxzOTFexc = fftshift(fft2(ifftshift(CoherentxzPSFexc)));
    CoherentzOTFexc = real(CoherentxzOTFexc(:,(N+1)/2));

    PSFoverall = PSFIncoherent .* PSFdet; 
    % PSFoverall = PSFoverall/max(PSFoverall,[],'all');
    OTFoverall = fftshift(fftn(ifftshift(PSFoverall)));
    xzPSFoverall = PSFoverall(:,:,(N+1)/2); 
    max(max(xzPSFoverall))
    xzOTFoverall = OTFoverall(:,:,(N+1)/2); 
    zOTFoverall = real(xzOTFoverall(:,(N+1)/2));

    IncoherentSumI = cumsum(IncoherentzPSFexc);
    IncoherentSumI = IncoherentSumI/max(IncoherentSumI,[],'all');
    [~,yzeropoint1,~,~] = intersections(1:length(IncoherentSumI),IncoherentSumI, ...
                                        [(N+1)/2,(N+1)/2],[0,1]);
    [~,yzeropoint2,~,~] = intersections(1:length(IncoherentSumI),IncoherentSumI, ...
                                        [(N+1)/2-1,(N+1)/2-1],[0,1]);
    yzeropoint = (yzeropoint1(1,1)+yzeropoint2(1,1))/2;
    absIncoherentSumI = abs(IncoherentSumI-yzeropoint);
    temp = zeros(1,length(absIncoherentSumI)+1);
    temp(1:(N+1)/2-1) = absIncoherentSumI(1:(N+1)/2-1);
    temp(1,(N+1)/2) = 0;
    temp(1,(N+1)/2+1:end) = absIncoherentSumI((N+1)/2:end);
    absIncoherentSumI = 2* temp(1:end-1);
    % find IFWHM
    [~,minindex] = min(absIncoherentSumI);
    [~,maxindex] = max(absIncoherentSumI);
    half = (absIncoherentSumI(maxindex)-absIncoherentSumI(minindex))/2;
    index1 = (absIncoherentSumI(minindex:end) <= half);
    index2 = (absIncoherentSumI(1:minindex) <= half);
    if ~isempty(index1) && ~isempty(index2)
        IncoherentIFWHM1 = Z_exc(minindex + find(index1,1,'last')-1);
        IncoherentIFWHM2 = Z_exc(find(index2,1,'first'));
        if abs(IncoherentIFWHM1) == abs(IncoherentIFWHM2)
            IncoherentIFWHM = abs(IncoherentIFWHM1) + abs(IncoherentIFWHM2);
        elseif abs(Z_exc(minindex) - IncoherentIFWHM1) > abs(Z_exc(minindex) - IncoherentIFWHM2)
            IncoherentIFWHM = abs(Z_exc(minindex) - IncoherentIFWHM1)*2;
        else
            IncoherentIFWHM = abs(Z_exc(minindex) - IncoherentIFWHM2)*2;
        end
    else
        IncoherentIFWHM = "N/A";
    end

    CoherentSumI = cumsum(CoherentzPSFexc);
    CoherentSumI = CoherentSumI/max(CoherentSumI,[],'all');
    [~,yzeropoint1,~,~] = intersections(1:length(CoherentSumI),CoherentSumI, ...
                                        [(N+1)/2,(N+1)/2],[0,1]);
    [~,yzeropoint2,~,~] = intersections(1:length(CoherentSumI),CoherentSumI, ...
                                        [(N+1)/2-1,(N+1)/2-1],[0,1]);
    yzeropoint = (yzeropoint1(1,1)+yzeropoint2(1,1))/2;
    absCoherentSumI = abs(CoherentSumI-yzeropoint);
    temp = zeros(1,length(absCoherentSumI)+1);
    temp(1:(N+1)/2-1) = absCoherentSumI(1:(N+1)/2-1);
    temp(1,(N+1)/2) = 0;
    temp(1,(N+1)/2+1:end) = absCoherentSumI((N+1)/2:end);
    absCoherentSumI = 2* temp(1:end-1);

    % fine IFWHM
    [~,minindex] = min(absCoherentSumI);
    [~,maxindex] = max(absCoherentSumI);
    half = (absCoherentSumI(maxindex)-absCoherentSumI(minindex))/2;
    index1 = (absCoherentSumI(minindex:end) <= half);
    index2 = (absCoherentSumI(1:minindex) <= half);
    if ~isempty(index1) && ~isempty(index2)
        CoherentIFWHM1 = Z_exc(minindex + find(index1,1,'last')-1);
        CoherentIFWHM2 = Z_exc(find(index2,1,'first'));
        if abs(CoherentIFWHM1) == abs(CoherentIFWHM2)
            CoherentIFWHM = abs(CoherentIFWHM1) + abs(CoherentIFWHM2);
        elseif abs(Z_exc(minindex) - CoherentIFWHM1) > abs(Z_exc(minindex) - CoherentIFWHM2)
            CoherentIFWHM = abs(Z_exc(minindex) - CoherentIFWHM1)*2;
        else
            CoherentIFWHM = abs(Z_exc(minindex) - CoherentIFWHM2)*2;
        end
    else
        CoherentIFWHM = "N/A";
    end
    
    yindex = 1-(IncoherentyzPSFexc((N+1)/2,:) <= 0.5*max(IncoherentyzPSFexc((N+1)/2,:)));
    yFWHM1 = find(yindex,1,'first') ;
    yFWHM2 = find(yindex,1,'last');
    Y_exc(yFWHM2)
    if isempty(yFWHM1)
        yFWHM1 = 1;
    end
    if isempty(yFWHM2)
        yFWHM2 = N;
    end

    [~,maxindex] = max(IncoherentyzPSFexc(:,(N+1)/2));
    index1 = (IncoherentyzPSFexc(maxindex:end,(N+1)/2) <= 0.5*max(IncoherentyzPSFexc(:,(N+1)/2)));
    index2 = (IncoherentyzPSFexc(1:maxindex,(N+1)/2) <= 0.5*max(IncoherentyzPSFexc(:,(N+1)/2)));
    if ~isempty(index1) && ~isempty(index2)
        zFWHM1 = Z_exc(maxindex + find(index1,1,'first')) ;
        zFWHM2 = Z_exc(find(index2,1,'last')-1);
        if abs(zFWHM1) == abs(zFWHM2)
            zFWHM = abs(zFWHM1) + abs(zFWHM2);
        elseif abs(Z_exc(maxindex) - zFWHM1) > abs(Z_exc(maxindex) - zFWHM2)
            zFWHM = abs(Z_exc(maxindex) - zFWHM1)*2;
        else
            zFWHM = abs(Z_exc(maxindex) - zFWHM2)*2;
        end
    else
        zFWHM = "N/A";
    end

    Pupil = SWPupil(:,:,1)+SWPupil(:,:,2);
    Pupil = Pupil/max(Pupil,[],'all');
    fig1 = figure;
    imagesc(KX_exc,KZ_exc,real(Pupil) )
    axis image
    xlabel("k_x/(4\pin/\lambda_{exc})");
    ylabel("k_z/(4\pin/\lambda_{exc})");
    xlim([-0.5,0.5])
    ylim([-0.5,0.5])
    colormap(firemap)
    colorbar
    hold on
    plot(real(Pupil(:,(N+1)/2))*Koffset,KZ_exc,'LineWidth',1,'Color','w')
    print(fig1, '-dsvg', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_PupilAmp_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig1, '-dpng', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_PupilAmp_weighting_' num2str(weighting) '.PNG'],'-r300')

    fig2 = figure;
    imagesc(KX_exc,KZ_exc,angle(Pupil) )
    axis image
    xlim([-0.5,0.5])
    ylim([-0.5,0.5])
    xlabel("k_x/(4\pin/\lambda_{exc})");
    ylabel("k_z/(4\pin/\lambda_{exc})");
    colormap(firemap)
    colorbar
    print(fig2, '-dsvg', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_PupilPhase_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig2, '-dpng', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_PupilPhase_weighting_' num2str(weighting) '.PNG'],'-r300')

    fig3 = figure;
    imagesc(X_exc,Z_exc,IncoherentxzPSFexc)
    axis image
    xlabel("x(\lambda_{exc}/n)");
    ylabel("z(\lambda_{exc}/n)");
    xlim([-20,20])
    ylim([-20,20])
    colormap(firemap)
    colorbar
    clim([0,1])
    hold on
    plot(IncoherentxzPSFexc(:,(N+1)/2)*Sampleoffset-Sampleoffset,Z_exc,'LineWidth',1,'Color','g')
    print(fig3, '-dsvg', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzPSF_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig3, '-dpng', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzPSF_weighting_' num2str(weighting) '.PNG'],'-r300')

    fig4 = figure;
    imagesc(X_exc,Z_exc,IncoherentPSFexc(:,:,yFWHM2))
    axis image
    xlabel("x(\lambda_{exc}/n)");
    ylabel("z(\lambda_{exc}/n)");
    xlim([-20,20])
    ylim([-20,20])
    colormap(firemap)
    colorbar
    clim([0,1])
    hold on
    plot(squeeze(IncoherentPSFexc(:,(N+1)/2,yFWHM2))*Sampleoffset-Sampleoffset,Z_exc,'LineWidth',1,'Color','g')
    print(fig4, '-dsvg', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzPSF_yFWHM_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig4, '-dpng', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzPSF_yFWHM_weighting_' num2str(weighting) '.PNG'],'-r300')

    fig5 = figure;
    imagesc(Y_exc,Z_exc,IncoherentyzPSFexc)
    axis image
    xlabel("y(\lambda_{exc}/n)");
    ylabel("z(\lambda_{exc}/n)");
    ylim([-40,40])
    xlim([-100,100])
    colormap(firemap)
    colorbar
    clim([0,1])
    hold on
    title("yFWHM=" + num2str(Y_exc(yFWHM2)*2) + ", zFWHM=" + num2str(zFWHM))
    xline(Y_exc(yFWHM2),'Color','magenta',"LineWidth",1,'LineStyle','--')
    xline(Y_exc(yFWHM1),'Color','magenta',"LineWidth",1,'LineStyle','--')
    plot(Y_exc,-IncoherentyzPSFexc((N+1)/2,:)*Sampleoffset+Sampleoffset,'LineWidth',1,'Color','g')
    print(fig5, '-dsvg', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_yzPSF_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig5, '-dpng', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_yzPSF_weighting_' num2str(weighting) '.PNG'],'-r300')

    fig6 = figure;
    imagesc(KX_exc,KZ_exc,real(IncoherentxzOTFexc)/max(real(IncoherentxzOTFexc),[],'all'))
    axis image
    xlabel("k_x/(4\pin/\lambda_{exc})");
    ylabel("k_z/(4\pin/\lambda_{exc})");
    xlim([-0.5,0.5])
    ylim([-0.5,0.5])
    colormap(firemap)
    colorbar
    clim([0,1])
    hold on
    plot(real(IncoherentzOTFexc)/max(real(IncoherentzOTFexc))/2,KZ_exc,'LineWidth',0.5,'Color','w')
    print(fig6, '-dsvg', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzOTFAmp_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig6, '-dpng', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzOTFAmp_weighting_' num2str(weighting) '.PNG'],'-r300')

    fig7 = figure;
    imagesc(KX_exc,KZ_exc,angle(IncoherentxzOTFexc))
    axis image
    xlabel("k_x/(4\pin/\lambda_{exc})");
    ylabel("k_z/(4\pin/\lambda_{exc})");
    xlim([-0.5,0.5])
    ylim([-0.5,0.5])
    colormap(firemap)
    colorbar
    hold on
    % plot(angle(IncoherentzOTFexc)/max(angle(IncoherentzOTFexc))/pi,KZ_exc,'LineWidth',2,'Color','green')
    print(fig7, '-dsvg', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzOTFPhase_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig7, '-dpng', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzOTFPhase_weighting_' num2str(weighting) '.PNG'],'-r300')

    fig8 = figure;
    imagesc(X_exc,Z_exc,xzPSFoverall)
    axis image
    xlabel("x(\lambda_{exc}/n)");
    ylabel("z(\lambda_{exc}/n)");
    xlim([-20,20])
    ylim([-20,20])
    colormap(firemap)
    colorbar
    clim([0,1])
    hold on
    plot(xzPSFoverall(:,(N+1)/2)*10-10,Z_exc,'LineWidth',1,'Color','g')
    print(fig8, '-dsvg', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzPSF_overall_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig8, '-dpng', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzPSF_overall_weighting_' num2str(weighting) '.PNG'],'-r300')

    fig9 = figure;
    imagesc(KX_exc,KZ_exc,real(xzOTFoverall)/max(real(xzOTFoverall),[],'all'))
    axis image
    xlabel("k_x/(4\pin/\lambda_{exc})");
    ylabel("k_z/(4\pin/\lambda_{exc})");
    xlim([-1,1])
    ylim([-1,1])
    colormap(firemap)
    colorbar
    clim([0,1])
    hold on
    plot(real(zOTFoverall)/max(real(zOTFoverall))/2,KZ_exc,'LineWidth',0.5,'Color','w')
    print(fig9, '-dsvg', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzOTF_overallAmp_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig9, '-dpng', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzOTF_overallAmp_weighting_' num2str(weighting) '.PNG'],'-r300')

    fig10 = figure;
    imagesc(KX_exc,KZ_exc,angle(xzOTFoverall))
    axis image
    xlabel("k_x/(4\pin/\lambda_{exc})");
    ylabel("k_z/(4\pin/\lambda_{exc})");
    xlim([-1,1])
    ylim([-1,1])
    colormap(firemap)
    colorbar
    hold on
    % plot(angle(IncoherentzOTFexc)/max(angle(IncoherentzOTFexc))/pi,KZ_exc,'LineWidth',2,'Color','green')
    print(fig10, '-dsvg', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzOTF_overallPhase_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig10, '-dpng', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzOTF_overallPhase_weighting_' num2str(weighting) '.PNG'],'-r300')

    fig13 = figure;
    imagesc(X_exc,Z_exc,IncoherentxyPSFexc)
    axis image
    xlabel("y(\lambda_{exc}/n)");
    ylabel("x(\lambda_{exc}/n)");
    xlim([-40,40])
    ylim([-40,40])
    colormap(firemap)
    colorbar
    clim([0,1])
    hold on
    print(fig13, '-dsvg', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xyPSF_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig13, '-dpng', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xyPSF_weighting_' num2str(weighting) '.PNG'],'-r300')

    fig14 = figure;
    plot(Z_exc,IncoherentxzPSFexc(:,(N+1)/2))
    hold on
    axis square
    plot(Z_exc,IncoherentxzPSFexc((N+1)/2,:))
    plot(Y_exc,IncoherentyzPSFexc((N+1)/2,:))
    xlabel("z(\lambda_{exc}/n)");
    ylabel("Intensity");
    xlim([-5,5])
    ylim([0,1])
    legend("z","x","y")
    grid on
    print(fig14, '-dsvg', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_lineProfile_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig14, '-dpng', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_lineProfile_weighting_' num2str(weighting) '.PNG'],'-r300')

    fig12 = figure;
    plot(Z_exc+0.125,IncoherentSumI,'Color','r','LineWidth',2)
    hold on
    grid on
    plot(Z_exc+0.125,CoherentSumI,'Color','g','LineWidth',1.5)
    plot(Z_exc,absIncoherentSumI,'LineWidth',2)
    plot(Z_exc,absCoherentSumI,'LineWidth',1.5)
    legend("Incoherent","Coherent","Incoherent","Coherent")
    axis square
    ylim([0,1])
    xlim([-40,40])
    yline(half)
    title("Incoherent IFWHM=" + num2str(IncoherentIFWHM) + ", Coherent IFWHM=" + num2str(CoherentIFWHM) )
    print(fig12, '-dsvg', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_zInt_' num2str(weighting) '.SVG'],'-r300')
    print(fig12, '-dpng', [savingdir  'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_zInt_' num2str(weighting) '.PNG'],'-r300')

    close all
    savingdir = [savingdir 'AxialProfile/'];
    mkdir(savingdir)

%     fig11 = figure;
%     counter = 1;
%     for jj = (N+1)/2-100:10:(N+1)/2+100
% 
%         PSF1 = IncoherentPSFexc(:,:,jj);
%         OTF1 = fftshift(fft2(ifftshift(PSF1))); 
% 
%         PSF2 = CoherentPSFexc(:,:,jj);
%         OTF2 = fftshift(fft2(ifftshift(PSF2))); 
% 
%         subplot(3,3,1);
%         imagesc(X_exc,Z_exc,PSF1);
%         title("Incoherent")
%         axis image;
%         xlim([-20,20])
%         ylim([-20,20])
%         xlabel("x/(\lambda_{exc}/n)")
%         ylabel("z/(\lambda_{exc}/n)")
%         colormap(firemap)
%         colorbar;
%         clim([0,1])
% 
%         subplot(3,3,2);
%         imagesc(X_exc,Z_exc,PSF2);
%         title("Coherent")
%         axis image;
%         xlim([-20,20])
%         ylim([-20,20])
%         xlabel("x/(\lambda_{exc}/n)")
%         ylabel("z/(\lambda_{exc}/n)")
%         colormap(firemap)
%         colorbar;
%         clim([0,1])
% 
%         subplot(3,3,3);
%         plot(Z_exc,PSF1(:,(N+1)/2),'Color','magenta','LineWidth',1);
%         hold on
%         plot(Z_exc,PSF2(:,(N+1)/2),'Color','g','LineWidth',1);
%         axis square
%         xlim([-10,10])
%         ylim([0,1])
%         xlabel("z/(\lambda_{exc}/n)")
%         ylabel("Intensity (a.u.)")
%         lgd = legend("Incoherent","Coherent");
%         lgd.FontSize = 3;
%         lgd.Box = 'off';
%         grid on
%         hold off
% 
%         subplot(3,3,4);
%         imagesc(KX_exc,KZ_exc,real(OTF1)/max(max(real(OTF1))));
%         axis image;
%         title("Amplitude")
%         xlabel("k_x/(4\pin/\lambda_{exc})");
%         ylabel("k_z/(4\pin/\lambda_{exc})");
%         xlim([-0.5,0.5])
%         ylim([-0.5,0.5])
%         colormap(firemap)
%         colorbar
%         clim([0,1])
% 
%         subplot(3,3,5);
%         imagesc(KX_exc,KZ_exc,real(OTF2)/max(max(real(OTF2))));
%         axis image;
%         xlabel("k_x/(4\pin/\lambda_{exc})");
%         ylabel("k_z/(4\pin/\lambda_{exc})");
%         xlim([-0.5,0.5])
%         ylim([-0.5,0.5])
%         colormap(firemap)
%         colorbar
%         clim([0,1])
% 
%         subplot(3,3,6);
%         plot(KZ_exc,real(OTF1(:,(N+1)/2))/max(real(OTF1(:,(N+1)/2))),'Color','magenta','LineWidth',1);
%         hold on
%         plot(KZ_exc,real(OTF2(:,(N+1)/2))/max(real(OTF2(:,(N+1)/2))),'Color','g','LineWidth',1);
%         axis square
%         xlabel("k_z/(4\pin/\lambda_{exc})");
%         ylabel("OTF strength")
%         xlim([-0.5,0.5])
%         ylim([-0.5,1])
%         lgd = legend("Incoherent","Coherent");
%         lgd.FontSize = 3;
%         lgd.Box = 'off';
%         grid on
%         hold off
% 
%         subplot(3,3,7);
%         imagesc(KX_exc,KZ_exc,angle(OTF1));
%         axis image;
%         title("Phase")
%         xlabel("k_x/(4\pin/\lambda_{exc})");
%         ylabel("k_z/(4\pin/\lambda_{exc})");
%         xlim([-0.5,0.5])
%         ylim([-0.5,0.5])
%         colormap(firemap)
%         colorbar
%         clim([0,1])
% 
%         subplot(3,3,8);
%         imagesc(KX_exc,KZ_exc,angle(OTF2));
%         axis image;
%         xlabel("k_x/(4\pin/\lambda_{exc})");
%         ylabel("k_z/(4\pin/\lambda_{exc})");
%         xlim([-0.5,0.5])
%         ylim([-0.5,0.5])
%         colormap(firemap)
%         colorbar
%         clim([0,1])
% 
%         subplot(3,3,9);
%         plot(KZ_exc,angle(OTF1(:,(N+1)/2)),'Color','magenta','LineWidth',1);
%         hold on
%         plot(KZ_exc,angle(OTF2(:,(N+1)/2)),'Color','g','LineWidth',1);
%         axis square
%         xlabel("k_z/(4\pin/\lambda_{exc})");
%         ylabel("OTF phase")
%         xlim([-0.5,0.5])
%         ylim([-pi,pi])
%         lgd = legend("Incoherent","Coherent");
%         lgd.FontSize = 3;
%         lgd.Box = 'off';
%         grid on
%         hold off
% 
%         print(fig11, '-dpng', [savingdir  num2str(counter) '_SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_weighting_' num2str(weighting) '_y_' num2str(Y_exc(jj)) '.PNG'],'-r300')
%         counter = counter + 1; 
%     end
% close all
end