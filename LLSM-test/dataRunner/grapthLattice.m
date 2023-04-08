function grapthLattice(NA1,deltaNA,LatticeType,weighting,LatticePSF,LatticePSFDithered,LatticePupil,PSFdet,savingdir)
    getParameters; %modify image parameter here
    CalculatePhysics;
    firemap = fire(256);

    Sampleoffset = 40;
    Koffset = 0.5;

    DitheredPSFexc = LatticePSFDithered;
    DitheredxzPSFexc = DitheredPSFexc(:,:,(N+1)/2); 
    DitheredyzPSFexc = squeeze(DitheredPSFexc(:,(N+1)/2,:)); 
    DitheredxzOTFexc = fftshift(fft2(ifftshift(DitheredxzPSFexc))); 
    DitheredzOTFexc = DitheredxzOTFexc(:,(N+1)/2); 

    PSFexc = LatticePSF;
    xzPSFexc = PSFexc(:,:,(N+1)/2); 
    yzPSFexc = squeeze(PSFexc(:,(N+1)/2,:)); 
    xzOTFexc = fftshift(fft2(ifftshift(xzPSFexc)));
    zOTFexc = real(xzOTFexc(:,(N+1)/2));

    PSFoverall = LatticePSFDithered .* PSFdet;
    OTFoverall = fftshift(fftn(ifftshift(PSFoverall)));
    xzPSFoverall = PSFoverall(:,:,(N+1)/2); 
    xzOTFoverall = OTFoverall(:,:,(N+1)/2); 
    zOTFoverall = real(xzOTFoverall(:,(N+1)/2));
    
    yindex = 1-(DitheredyzPSFexc((N+1)/2,:) <= 0.5*max(DitheredyzPSFexc((N+1)/2,:)));
    yFWHM1 = find(yindex,1,'first') ;
    yFWHM2 = find(yindex,1,'last');
    Y_exc(yFWHM2)
    if isempty(yFWHM1)
        yFWHM1 = 1;
    end
    if isempty(yFWHM2)
        yFWHM2 = N;
    end

    fig1 = figure;
    imagesc(KX_exc,KZ_exc,real(LatticePupil)./max(real(LatticePupil),[],'all') )
    axis image
    xlabel("k_x/(4\pin/\lambda_{exc})");
    ylabel("k_z/(4\pin/\lambda_{exc})");
    xlim([-0.5,0.5])
    ylim([-0.5,0.5])
    colormap(firemap)
    colorbar
    print(fig1, '-dsvg', [savingdir  LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_PupilAmp_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig1, '-dpng', [savingdir  LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_PupilAmp_weighting_' num2str(weighting) '.PNG'],'-r300')

    fig2 = figure;
    imagesc(KX_exc,KZ_exc,angle(LatticePupil) )
    axis image
    xlim([-0.5,0.5])
    ylim([-0.5,0.5])
    xlabel("k_x/(4\pin/\lambda_{exc})");
    ylabel("k_z/(4\pin/\lambda_{exc})");
    colormap(firemap)
    colorbar
    print(fig2, '-dsvg', [savingdir  LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_PupilPhase_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig2, '-dpng', [savingdir  LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_PupilPhase_weighting_' num2str(weighting) '.PNG'],'-r300')

    fig3 = figure;
    imagesc(X_exc,Z_exc,DitheredxzPSFexc)
    axis image
    xlabel("x(\lambda_{exc}/n)");
    ylabel("z(\lambda_{exc}/n)");
    xlim([-40,40])
    ylim([-40,40])
    colormap(firemap)
    colorbar
    clim([0,1])
    hold on
    plot(DitheredxzPSFexc(:,(N+1)/2)*Sampleoffset-Sampleoffset,Z_exc,'LineWidth',1,'Color','g')
    print(fig3, '-dsvg', [savingdir  LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzPSF_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig3, '-dpng', [savingdir  LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzPSF_weighting_' num2str(weighting) '.PNG'],'-r300')

    fig4 = figure;
    imagesc(X_exc,Z_exc,DitheredPSFexc(:,:,yFWHM2))
    axis image
    xlabel("x(\lambda_{exc}/n)");
    ylabel("z(\lambda_{exc}/n)");
    xlim([-40,40])
    ylim([-40,40])
    colormap(firemap)
    colorbar
    clim([0,1])
    hold on
    plot(squeeze(DitheredPSFexc(:,(N+1)/2,yFWHM2))*Sampleoffset-Sampleoffset,Z_exc,'LineWidth',1,'Color','g')
    print(fig4, '-dsvg', [savingdir   LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzPSF_yFWHM_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig4, '-dpng', [savingdir   LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzPSF_yFWHM_weighting_' num2str(weighting) '.PNG'],'-r300')

    fig5 = figure;
    imagesc(Y_exc,Z_exc,DitheredyzPSFexc)
    axis image
    xlabel("y(\lambda_{exc}/n)");
    ylabel("z(\lambda_{exc}/n)");
    ylim([-40,40])
    xlim([-60,60])
    colormap(firemap)
    colorbar
    clim([0,1])
    hold on
    title(num2str(Y_exc(yFWHM2)))
    xline(Y_exc(yFWHM2),'Color','magenta',"LineWidth",1,'LineStyle','--')
    xline(Y_exc(yFWHM1),'Color','magenta',"LineWidth",1,'LineStyle','--')
    plot(Y_exc,-DitheredyzPSFexc((N+1)/2,:)*Sampleoffset+Sampleoffset,'LineWidth',1,'Color','g')
    print(fig5, '-dsvg', [savingdir   LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_yzPSF_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig5, '-dpng', [savingdir   LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_yzPSF_weighting_' num2str(weighting) '.PNG'],'-r300')

    fig6 = figure;
    imagesc(KX_exc,KZ_exc,real(DitheredxzOTFexc)/max(real(DitheredxzOTFexc),[],'all'))
    axis image
    xlabel("k_x/(4\pin/\lambda_{exc})");
    ylabel("k_z/(4\pin/\lambda_{exc})");
    xlim([-0.5,0.5])
    ylim([-0.5,0.5])
    colormap(firemap)
    colorbar
    clim([0,1])
    hold on
    plot(real(DitheredzOTFexc)/max(real(DitheredzOTFexc))/2,KZ_exc,'LineWidth',0.5,'Color','w')
    print(fig6, '-dsvg', [savingdir   LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzOTFAmp_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig6, '-dpng', [savingdir   LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzOTFAmp_weighting_' num2str(weighting) '.PNG'],'-r300')

    fig7 = figure;
    imagesc(KX_exc,KZ_exc,angle(DitheredxzOTFexc))
    axis image
    xlabel("k_x/(4\pin/\lambda_{exc})");
    ylabel("k_z/(4\pin/\lambda_{exc})");
    xlim([-0.5,0.5])
    ylim([-0.5,0.5])
    colormap(firemap)
    colorbar
    hold on
    print(fig7, '-dsvg', [savingdir   LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzOTFPhase_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig7, '-dpng', [savingdir   LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzOTFPhase_weighting_' num2str(weighting) '.PNG'],'-r300')

    fig8 = figure;
    imagesc(X_exc,Z_exc,xzPSFoverall)
    axis image
    xlabel("x(\lambda_{exc}/n)");
    ylabel("z(\lambda_{exc}/n)");
    xlim([-10,10])
    ylim([-10,10])
    colormap(firemap)
    colorbar
    clim([0,1])
    hold on
    plot(xzPSFoverall(:,(N+1)/2)*10-10,Z_exc,'LineWidth',1,'Color','g')
    print(fig8, '-dsvg', [savingdir   LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzPSF_overall_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig8, '-dpng', [savingdir   LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzPSF_overall_weighting_' num2str(weighting) '.PNG'],'-r300')

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
    plot(real(zOTFoverall)/max(real(zOTFoverall)),KZ_exc,'LineWidth',0.5,'Color','w')
    print(fig9, '-dsvg', [savingdir   LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzOTF_overallAmp_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig9, '-dpng', [savingdir   LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzOTF_overallAmp_weighting_' num2str(weighting) '.PNG'],'-r300')

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
    print(fig10, '-dsvg', [savingdir   LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzOTF_overallPhase_weighting_' num2str(weighting) '.SVG'],'-r300')
    print(fig10, '-dpng', [savingdir   LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_xzOTF_overallPhase_weighting_' num2str(weighting) '.PNG'],'-r300')

    close all
    savingdir = [savingdir 'AxialProfile/'];
    mkdir(savingdir)

    fig11 = figure;
    counter = 1;
    for jj = (N+1)/2-100:10:(N+1)/2+100

        PSF1 = LatticePSF(:,:,jj);
        OTF1 = fftshift(fft2(ifftshift(PSF1))); 

        PSF2 = LatticePSFDithered(:,:,jj);
        OTF2 = fftshift(fft2(ifftshift(PSF2))); 
       
        subplot(3,3,1);
        imagesc(X_exc,Z_exc,PSF1);
        axis image;
        xlim([-40,40])
        ylim([-40,40])
        xlabel("x/(\lambda_{exc}/n)")
        ylabel("z/(\lambda_{exc}/n)")
        colormap(firemap)
        colorbar;
        clim([0,1])
        
        subplot(3,3,2);
        imagesc(X_exc,Z_exc,PSF2);
        axis image;
        xlim([-40,40])
        ylim([-40,40])
        xlabel("x/(\lambda_{exc}/n)")
        ylabel("z/(\lambda_{exc}/n)")
        colormap(firemap)
        colorbar;
        clim([0,1])
        
        subplot(3,3,3);
        plot(Z_exc,PSF1(:,(N+1)/2),'Color','magenta','LineWidth',1);
        hold on
        plot(Z_exc,PSF2(:,(N+1)/2),'Color','g','LineWidth',1);
        axis square
        xlim([-40,40])
        ylim([0,1])
        xlabel("z/(\lambda_{exc}/n)")
        ylabel("Intensity (a.u.)")
        lgd = legend("LLS","Dithered");
        lgd.FontSize = 3;
        lgd.Box = 'off';
        grid on
        hold off

        subplot(3,3,4);
        imagesc(KX_exc,KZ_exc,real(OTF1)/max(max(real(OTF1))));
        axis image;
        title("Amplitude")
        xlabel("k_x/(4\pin/\lambda_{exc})");
        ylabel("k_z/(4\pin/\lambda_{exc})");
        xlim([-0.5,0.5])
        ylim([-0.5,0.5])
        colormap(firemap)
        colorbar
        clim([0,1])

        subplot(3,3,5);
        imagesc(KX_exc,KZ_exc,real(OTF2)/max(max(real(OTF2))));
        axis image;
        xlabel("k_x/(4\pin/\lambda_{exc})");
        ylabel("k_z/(4\pin/\lambda_{exc})");
        xlim([-0.5,0.5])
        ylim([-0.5,0.5])
        colormap(firemap)
        colorbar
        clim([0,1])

        subplot(3,3,6);
        plot(KZ_exc,real(OTF1(:,(N+1)/2))/max(real(OTF1(:,(N+1)/2))),'Color','magenta','LineWidth',1);
        hold on
        plot(KZ_exc,real(OTF2(:,(N+1)/2))/max(real(OTF2(:,(N+1)/2))),'Color','g','LineWidth',1);
        axis square
        xlabel("k_z/(4\pin/\lambda_{exc})");
        ylabel("OTF strength")
        xlim([-0.5,0.5])
        ylim([-0.5,1])
        lgd = legend("LLS","Dithered");
        lgd.FontSize = 3;
        lgd.Box = 'off';
        grid on
        hold off

        subplot(3,3,7);
        imagesc(KX_exc,KZ_exc,angle(OTF1));
        axis image;
        title("Phase")
        xlabel("k_x/(4\pin/\lambda_{exc})");
        ylabel("k_z/(4\pin/\lambda_{exc})");
        xlim([-0.5,0.5])
        ylim([-0.5,0.5])
        colormap(firemap)
        colorbar
        clim([0,1])

        subplot(3,3,8);
        imagesc(KX_exc,KZ_exc,angle(OTF2));
        axis image;
        xlabel("k_x/(4\pin/\lambda_{exc})");
        ylabel("k_z/(4\pin/\lambda_{exc})");
        xlim([-0.5,0.5])
        ylim([-0.5,0.5])
        colormap(firemap)
        colorbar
        clim([0,1])

        subplot(3,3,9);
        plot(KZ_exc,angle(OTF1(:,(N+1)/2)),'Color','magenta','LineWidth',1);
        hold on
        plot(KZ_exc,angle(OTF2(:,(N+1)/2)),'Color','g','LineWidth',1);
        axis square
        xlabel("k_z/(4\pin/\lambda_{exc})");
        ylabel("OTF phase")
        xlim([-0.5,0.5])
        ylim([-pi,pi])
        lgd = legend("LLS","Dithered");
        lgd.FontSize = 3;
        lgd.Box = 'off';
        grid on
        hold off

        print(fig11, '-dpng', [savingdir  num2str(counter) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_weighting_' num2str(weighting) '_y_' num2str(Y_exc(jj)) '.PNG'],'-r300')
        counter = counter + 1;     
    end
    close all
end
