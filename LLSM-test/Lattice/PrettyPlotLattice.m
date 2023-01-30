function PrettyPlotLattice(LatticePupil,LatticeMask,LatticeMetaData,LatticePSF,LatticePSFDithered)

    getParameters; %modify image parameter here
    CalculatePhysics;
    NAmax = LatticeMetaData.NAmax;
    NAmin = LatticeMetaData.NAmin;
    MaskNAmax = LatticeMetaData.MaskNAmax;
    MaskNAmin = LatticeMetaData.MaskNAmin;

    fig1 = figure;
    fig1.Name = "Rear Pupil," + "NAmax = " + num2str(NAmax) +...
                  ", NAmin = " + num2str(NAmin); 
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
    colormap(hot(256)) 
    
    subplot(2,4,1)
    xzPSF_exc = LatticePSF(:,:,(N+1)/2);
    xzPSF_exc = xzPSF_exc/max(max(xzPSF_exc));
    image21 = imagesc(X_exc, Z_exc, xzPSF_exc);
    title("XZ-Excitation PSF")
    xlabel("x/(\lambda_{exc}/n)")
    ylabel("z/(\lambda_{exc}/n)")
    colorbar;
    axis image;
    image21.Parent.XLim = [-20,20];
    image21.Parent.YLim = [-20,20];
    
    subplot(2,4,2)
    zPSF_exc = xzPSF_exc(:,(N+1)/2);
    image22 = plot(Z_exc,zPSF_exc);
    title("Z-Excitation PSF, " + ...
        "X=0" + ",Y=0" )
    ylabel("Normalized a.u. ")
    xlabel("z/(\lambda_{exc}/n)")
    image22.LineWidth = 2;
    image22.Color = 'r';
    image22.Parent.XLim = [-10,10];
    image22.Parent.YAxis.TickValues = linspace(0,1,11);
    image22.Parent.XAxis.TickValues = linspace(-10,10,21);
    grid on
    axis square
    
    subplot(2,4,3)
    xzOTF_exc = fftshift(fft2(xzPSF_exc));
    xzOTF_exc = xzOTF_exc/max(max(xzOTF_exc));
    image18 = imagesc( KX_exc,...
                  KZ_exc,...
                 abs(xzOTF_exc)) ;
    title("XZ-Excitation OTF, K_Y=0")
    xlabel("k_x/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colorbar;
    axis image
    image18.Parent.XLim = [-0.5,0.5];
    image18.Parent.YLim = [-0.5,0.5];
    
    subplot(2,4,4)
    zOTF_exc = abs(xzOTF_exc(:,(N+1)/2));
    image23 = plot( KZ_exc, zOTF_exc);
    title("Z-Excitation-OTF, " + "K_X=0, " + "K_Y=0")
    ylabel("Normalized a.u. ")
    xlabel("k_z/(4\pin/\lambda_{exc})")
    image23.Color = 'r';
    image23.LineWidth = 2;
    image23.Parent.XLim = [-1,1];
    image23.Parent.YAxis.TickValues = linspace(0,1,11);
    image23.Parent.XAxis.TickValues = linspace(-0.5,0.5,11);
    image23.Parent.XLim = [-0.5,0.5];
    grid on
    axis square
    
    subplot(2,4,5)
    hold on;
    xzPSF_exc_dither = LatticePSFDithered(:,:,(N+1)/2); 
    xzPSF_exc_dither = xzPSF_exc_dither/max(max(xzPSF_exc_dither));
    image19 = imagesc(X_exc, Z_exc ,xzPSF_exc_dither);
    title("Dithered XZ-Excitation PSF")
    xlabel("x/(\lambda_{exc}/n)")
    ylabel("z/(\lambda_{exc}/n)")
    colorbar;
    axis image;
    image19.Parent.XLim = [-20,20];
    image19.Parent.YLim = [-20,20];
    
    subplot(2,4,6)
    zPSF_exc_dither = xzPSF_exc_dither(:,(N+1)/2);
    image110 = plot( Z_exc, zPSF_exc_dither);
        title("Z-Excitation PSF, " + ...
        "X=0" + ",Y=0" )
    ylabel("Normalized a.u. ")
    xlabel("z/(\lambda_{exc}/n)")
    image110.LineWidth = 2;
    image110.Color = 'r';
    image110.Parent.YAxis.TickValues = linspace(0,1,11);
    image110.Parent.XAxis.TickValues = linspace(-10,10,21);
    image110.Parent.XLim = [-10,10];
    grid on
    axis square
    
    subplot(2,4,7)
    xzOTF_exc_dither = fftshift(fft2(xzPSF_exc_dither));  
    image111 = imagesc( KX_exc,...
                    KZ_exc,...
                    abs(xzOTF_exc_dither)/max(max(abs(xzOTF_exc_dither))));
    title("Dithered XZ-Excitation OTF, K_Y=0")
    xlabel("k_x/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colorbar;
    axis image
    image111.Parent.XLim = [-0.5,0.5];
    image111.Parent.YLim = [-0.5,0.5];
    
    subplot(2,4,8)
    zOTF_exc_dither = abs(xzOTF_exc_dither(:,(N+1)/2));
    image112 = plot( KZ_exc, zOTF_exc_dither/max(max(zOTF_exc_dither)));
    title("Dithered Z-Excitation-OTF, " + "K_X=0, " + "K_Y=0")
    ylabel("Normalized a.u. ")
    xlabel("k_z/(4\pin/\lambda_{exc})")
    image112.LineWidth = 2;
    image112.Color = 'r';
    image112.Parent.YAxis.TickValues = linspace(0,1,11);
    image112.Parent.XAxis.TickValues = linspace(-0.5,0.5,11);
    image112.Parent.XLim = [-0.5,0.5];
    grid on
    axis square
    %% 
    fig3 = figure;
    fig3.Name = "Propagation";
    colormap(hot(256)) 
    
    subplot(2,4,1:2);
    yzPSF_exc = squeeze(LatticePSF(:,(N+1)/2,:));
    yzPSF_exc = yzPSF_exc/max(max(yzPSF_exc));
    image24 = imagesc(Y_exc, Z_exc, yzPSF_exc);
    title("YZ-Excitation-PSF, " + "X = 0" )
    axis image
    xlabel("y/(\lambda_{exc}/n)")
    ylabel("z/(\lambda_{exc}/n)")
    image24.Parent.YLim = [-40,40];
    colorbar
    
    subplot(2,4,5:6);
    yzPSF_exc_dither = squeeze(LatticePSFDithered(:,(N+1)/2,:));
    image25 = imagesc(Y_exc,Z_exc, yzPSF_exc_dither );
    title("Dithered YZ-Excitation-PSF, " + "X = 0" )
    xlabel("y/(\lambda_{exc}/n)")
    ylabel("z/(\lambda_{exc}/n)")
    axis image
    image25.Parent.YLim = [-40,40];
    
    colorbar;
    
    subplot(2,4,[3:4 7:8]);
    hold on
    yPSF_exc = squeeze(LatticePSF((N+1)/2,(N+1)/2,:));
    yPSF_exc = yPSF_exc/max(yPSF_exc);
    yPSF_exc_dither = squeeze(LatticePSFDithered((N+1)/2,(N+1)/2,:));
    yPSF_exc_dither = yPSF_exc_dither/max(yPSF_exc_dither);
    
    % Calculate yFWHM
    index = find(yPSF_exc((N+1)/2:end) < 0.5);
    index_dither = find(yPSF_exc_dither((N+1)/2:end) < 0.5);
    if ~isempty(index)
        yFWHM = 2*Y_exc((N+1)/2+index(1)-1);
    else
        yFWHM = "N/A";
    end
    if ~isempty(index_dither)
        yFWHM_dither = 2*Y_exc((N+1)/2+index_dither(1)-1);
    else
        yFWHM = "N/A";
    end   
    image26_1 = plot(Y_exc, yPSF_exc );
    image26_2 = plot(Y_exc, yPSF_exc_dither);
    title("Y-Excitation-PSF, " + "X = 0, Z = 0, " + ...
          "yFWHM = " + num2str(yFWHM) + "  \lambda_{exc}/n, " +...
          "yFWHM_{dither} = " + num2str(yFWHM_dither) + "\lambda_{exc}/n")
    xlabel("y/(\lambda_{exc}/n)")
    ylabel("Normalized a.u. ")
    image26_1.Color = 'g';
    image26_1.LineWidth = 2;
    image26_2.Color = 'r';
    image26_2.LineWidth = 2;
    grid on
    axis square