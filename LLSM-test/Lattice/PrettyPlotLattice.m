function PrettyPlotLattice(LatticePupil,LatticeMask,LatticeMetaData,LatticePSF,LatticePSFDithered,Latticecenter)

    getParameters; %modify image parameter here
    CalculatePhysics;
    NAmax = LatticeMetaData.NAmax;
    NAmin = LatticeMetaData.NAmin;
    MaskNAmax = LatticeMetaData.MaskNAmax;
    MaskNAmin = LatticeMetaData.MaskNAmin;
    
    % check center value and index for aberration 
    [PSFcenterInt,PSFcenter] = max(LatticePSFDithered,[],'all'); % value, index
    if isequal(PSFcenter,Latticecenter(2,2)) && isequal(PSFcenterInt,Latticecenter(2,1))
        LatticePSF = LatticePSF/max(max(max(LatticePSF)));
        LatticePSFDithered = LatticePSFDithered/max(max(max(LatticePSFDithered)));
    else
        LatticePSF = LatticePSF/Latticecenter(1,1);
        LatticePSFDithered = LatticePSFDithered/Latticecenter(2,1);
    end

    xzPSFexc = LatticePSF(:,:,(N+1)/2);
    zPSFexc = xzPSFexc(:,(N+1)/2);
    xzOTFexc = real(fftshift(fft2(ifftshift(xzPSFexc))));
    xzOTFexc = xzOTFexc/max(max(xzOTFexc));
    zOTFexc = xzOTFexc(:,(N+1)/2);

    xzPSFexc_dither = LatticePSFDithered(:,:,(N+1)/2); 
    zPSFexc_dither = xzPSFexc_dither(:,(N+1)/2);
    xzOTFexc_dither = real(fftshift(fft2(ifftshift(xzPSFexc_dither)))); 
    xzOTFexc_dither = xzOTFexc_dither./max(max(xzOTFexc_dither));
    zOTFexc_dither = xzOTFexc_dither(:,(N+1)/2);

    yzPSFexc = squeeze(LatticePSF(:,(N+1)/2,:));
    yzPSFexc_dither = squeeze(LatticePSFDithered(:,(N+1)/2,:));
    yPSFexc = yzPSFexc((N+1)/2,:);
    yPSFexc_dither = yzPSFexc_dither((N+1)/2,:);

%%
    fig1 = figure;
    fig1.Name = "Rear Pupil," + "NAmax = " + num2str(NAmax) +...
                  ", NAmin = " + num2str(NAmin); 
    fig1.WindowState = 'maximized';
    colormap(hot(256))
    
    subplot(1,2,1);
    image11 = imagesc(KX_exc,KZ_exc, real(LatticePupil) );
    title("Rear Pupil, " +...
          "N_{xz} = " + num2str(N))
    xlabel("k_x/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    axis image
    image11.Parent.XLim = [-0.5,0.5];
    image11.Parent.YLim = [-0.5,0.5];
    colorbar;    
    
    Illum_mask = imfuse(real(LatticePupil),LatticeMask,"falsecolor","ColorChannels","green-magenta");
    subplot(1,2,2)
    image15 = imagesc( KX_exc, KZ_exc,...
                  Illum_mask);
    title("Masking, " +...
          "NA_{max} = " + num2str(MaskNAmax) +...
          ", NA_{min} = " + num2str(MaskNAmin) )
    xlabel("k_x/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    axis image
    image15.Parent.XLim = [-0.5,0.5];
    image15.Parent.YLim = [-0.5,0.5];
    
    %% Figure 2 - Excitation
    fig2 = figure;
    fig2.Name = "Focal PSF/OTF";
    fig2.WindowState = 'maximized';
    colormap(hot(256)) 
    
    subplot(2,4,1)
    image21 = imagesc(X_exc, Z_exc, xzPSFexc);
    title("XZ-Excitation PSF")
    xlabel("x/(\lambda_{exc}/n)")
    ylabel("z/(\lambda_{exc}/n)")
    colorbar;
    axis image;
    image21.Parent.XLim = [-20,20];
    image21.Parent.YLim = [-20,20];
    
    subplot(2,4,2)
    image22 = plot(Z_exc,zPSFexc);
    title("Z-Excitation PSF, " + ...
        "X=0" + ",Y=0" )
    ylabel("Normalized a.u. ")
    xlabel("z/(\lambda_{exc}/n)")
    image22.LineWidth = 2;
    image22.Color = 'r';
    image22.Parent.XLim = [-10,10];
    image22.Parent.YAxis.TickValues = linspace(0,1,11);
    image22.Parent.XAxis.TickValues = linspace(-10,10,21);
    image22.Parent.YLim = [0,1];
    grid on
    axis square
    
    subplot(2,4,3)
    image18 = imagesc( KX_exc,...
                  KZ_exc,...
                 xzOTFexc) ;
    title("XZ-Excitation real(OTF), K_Y=0")
    xlabel("k_x/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colorbar;
    axis image
    image18.Parent.XLim = [-0.5,0.5];
    image18.Parent.YLim = [-0.5,0.5];
    
    subplot(2,4,4)
    image23 = plot( KZ_exc, zOTFexc);
    title("Z-Excitation real(OTF), " + "K_X=0, " + "K_Y=0")
    ylabel("Normalized a.u. ")
    xlabel("k_z/(4\pin/\lambda_{exc})")
    image23.Color = 'r';
    image23.LineWidth = 2;
    image23.Parent.XLim = [-0.3,1];
%     image23.Parent.YAxis.TickValues = linspace(0,1,11);
    image23.Parent.XAxis.TickValues = linspace(-0.5,0.5,11);
    image23.Parent.XLim = [-0.5,0.5];
    image23.Parent.YLim = [-0.3,1];
    grid on
    axis square
    
    subplot(2,4,5)
    hold on;
    image19 = imagesc(X_exc, Z_exc ,xzPSFexc_dither);
    title("Dithered XZ-Excitation PSF")
    xlabel("x/(\lambda_{exc}/n)")
    ylabel("z/(\lambda_{exc}/n)")
    colorbar;
    axis image;
    image19.Parent.XLim = [-20,20];
    image19.Parent.YLim = [-20,20];
    
    subplot(2,4,6)
    image110 = plot( Z_exc, zPSFexc_dither);
        title("Z-Excitation PSF, " + ...
        "X=0" + ",Y=0" )
    ylabel("Normalized a.u. ")
    xlabel("z/(\lambda_{exc}/n)")
    image110.LineWidth = 2;
    image110.Color = 'r';
    image110.Parent.YAxis.TickValues = linspace(0,1,11);
    image110.Parent.XAxis.TickValues = linspace(-10,10,21);
    image110.Parent.XLim = [-10,10];
    image110.Parent.YLim = [0,1];
    grid on
    axis square
    
    subplot(2,4,7)
    image111 = imagesc( KX_exc,...
                    KZ_exc,...
                    xzOTFexc_dither );
    title("Dithered XZ-Excitation real(OTF), K_Y=0")
    xlabel("k_x/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colorbar;
    axis image
    image111.Parent.XLim = [-0.5,0.5];
    image111.Parent.YLim = [-0.5,0.5];
    
    subplot(2,4,8)
    image112 = plot( KZ_exc, zOTFexc_dither);
    title("Dithered Z-Excitation real(OTF), " + "K_X=0, " + "K_Y=0")
    ylabel("Normalized a.u. ")
    xlabel("k_z/(4\pin/\lambda_{exc})")
    image112.LineWidth = 2;
    image112.Color = 'r';
%     image112.Parent.YAxis.TickValues = linspace(0,1,11);
    image112.Parent.XAxis.TickValues = linspace(-0.5,0.5,11);
    image112.Parent.XLim = [-0.5,0.5];
    image112.Parent.YLim = [-0.3,1];
    grid on
    axis square
    %% 
    fig3 = figure;
    fig3.Name = "Propagation";
    fig3.WindowState = 'maximized';
    colormap(hot(256)) 
    
    subplot(2,4,1:2);
    image24 = imagesc(Y_exc, Z_exc, yzPSFexc);
    title("YZ-Excitation-PSF, " + "X = 0" )
    axis image
    xlabel("y/(\lambda_{exc}/n)")
    ylabel("z/(\lambda_{exc}/n)")
    image24.Parent.YLim = [-40,40];
    colorbar
    
    subplot(2,4,5:6);
    image25 = imagesc(Y_exc,Z_exc, yzPSFexc_dither );
    title("Dithered YZ-Excitation-PSF, " + "X = 0" )
    xlabel("y/(\lambda_{exc}/n)")
    ylabel("z/(\lambda_{exc}/n)")
    axis image
    image25.Parent.YLim = [-40,40];
    colorbar;


    subplot(2,4,[3:4 7:8]);
    hold on
    % Calculate yFWHM
    [~,maxindex] = max(yPSFexc);
    index = 1-(yPSFexc <= 0.5*max(yPSFexc));
    if ~isempty(index)
        yFWHM1 = Y_exc(find(index,1,'first')) ;
        yFWHM2 = Y_exc(find(index,1,'last'));
        if abs(yFWHM1) == abs(yFWHM2)
            yFWHM = abs(yFWHM1) + abs(yFWHM2);
        elseif abs(Y_exc(maxindex) - yFWHM1) > abs(Y_exc(maxindex) - yFWHM2)
            yFWHM = abs(Y_exc(maxindex) - yFWHM1)*2;
        else
            yFWHM = abs(Y_exc(maxindex) - yFWHM2)*2;
        end
    else
        yFWHM = "N/A";
    end

    [~,maxindex] = max(yPSFexc_dither);
    index_dither = 1-(yPSFexc_dither <= 0.5*max(yPSFexc_dither));
    if ~isempty(index_dither)
        yFWHM1 = Y_exc(find(index_dither,1,'first')) ;
        yFWHM2 = Y_exc(find(index_dither,1,'last'));
        if abs(yFWHM1) == abs(yFWHM2)
            yFWHM_dither = abs(yFWHM1) + abs(yFWHM2);
        elseif abs(Y_exc(maxindex) - yFWHM1) > abs(Y_exc(maxindex) - yFWHM2)
            yFWHM_dither = abs(Y_exc(maxindex) - yFWHM1)*2;
        else
            yFWHM_dither = abs(Y_exc(maxindex) - yFWHM2)*2;
        end
    else
        yFWHM_dither = "N/A";
    end   
    image26_1 = plot(Y_exc, yPSFexc );
    image26_2 = plot(Y_exc, yPSFexc_dither);
    title("Y-Excitation-PSF, " + "X = 0, Z = 0, " + ...
          "yFWHM = " + num2str(yFWHM) + "  \lambda_{exc}/n, " +...
          "yFWHM_{dither} = " + num2str(yFWHM_dither) + "\lambda_{exc}/n")
    xlabel("y/(\lambda_{exc}/n)")
    ylabel("Normalized a.u. ")
    image26_1.Color = 'g';
    image26_1.LineWidth = 2;
    image26_2.Color = 'r';
    image26_2.LineWidth = 2;
    xlim([-150,150])
    ylim([0,1])
    grid on
    axis square
    legend("lattice","Dithered Lattice")