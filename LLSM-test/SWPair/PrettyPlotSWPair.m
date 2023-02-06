function PrettyPlotSWPair(SWPupil,SWMask,SWPupilMeta,PSFCoherent,PSFIncoherent,SWcenter)
    % Display Plots based on Coherent/Incoherent pair of SW propgation 
getParameters;
CalculatePhysics;

Pupil1 = squeeze(SWPupil(:,:,1));
Pupil2 = squeeze(SWPupil(:,:,2));
Pupil_sum = Pupil1 + Pupil2;
A_mask1 = squeeze(SWMask(:,:,1));
A_mask2 = squeeze(SWMask(:,:,2));
beam1NAmax = SWPupilMeta.NA1max;
beam1NAmin = SWPupilMeta.NA1min;
beam2NAmax = SWPupilMeta.NA2max;
beam2NAmin = SWPupilMeta.NA2min;
WeightingRatio = SWPupilMeta.WeightingRatio;

[PSFcenterInt,PSFcenter] = max(PSFIncoherent,[],'all'); % value, index
if isequal(PSFcenter,SWcenter(2,2)) && isequal(PSFcenterInt,SWcenter(2,1))
    PSFCoherent = PSFCoherent/max(max(max(PSFCoherent)));
    PSFIncoherent = PSFIncoherent/max(max(max(PSFIncoherent)));
else
    PSFCoherent = PSFCoherent/SWcenter(1,1);
    PSFIncoherent = PSFIncoherent/SWcenter(2,1);
end

xzPSFIncoherent = PSFIncoherent(:,:,(N+1)/2);
zPSFIncoherent = xzPSFIncoherent(:,(N+1)/2);
yzPSFIncoherent = squeeze(PSFIncoherent(:,(N+1)/2,:));
yPSFIncoherent = yzPSFIncoherent((N+1)/2,:);

xzPSFCoherent = PSFCoherent(:,:,(N+1)/2);
% zPSFCoherent = xzPSFCoherent(:,(N+1)/2);
yzPSFCoherent = squeeze(PSFCoherent(:,(N+1)/2,:));

xzOTFIncoherent = real(fftshift(fft2(ifftshift(xzPSFIncoherent))));
xzOTFIncoherent = xzOTFIncoherent./max(max(xzOTFIncoherent));
zOTFIncoherent = xzOTFIncoherent(:,(N+1)/2);
xzOTFCoherent = real(fftshift(fft2(ifftshift(xzPSFCoherent))));
xzOTFCoherent = xzOTFCoherent./max(max(xzOTFCoherent));
% zOTFCoherent = xzOTFCoherent(:,(N+1)/2);


%% Figure 1 - Rear Pupil 
    fig1 = figure;
    fig1.Name = "Rear Pupil";
    fig1.WindowState = 'maximized';
    colormap(hot(256))

    subplot(1,2,1)
    Illum_mask1 = imfuse(real(Pupil1),A_mask1,"falsecolor","ColorChannels","green-magenta");
    Illum_mask2 = imfuse(real(Pupil2),A_mask2,"falsecolor","ColorChannels","green-magenta");
    Illum_mask = Illum_mask1 + Illum_mask2;
    image15 = imagesc( KX_exc, KZ_exc,...
                      real(Illum_mask));
    title("beam1=" + num2str(beam1NAmax) + "/" + num2str(beam1NAmin)...
            +",beam2=" + num2str(beam2NAmax) + "/" + num2str(beam2NAmin) ) 
    xlabel("k_x/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    axis image
    image15.Parent.XLim = [-0.5,0.5];
    image15.Parent.YLim = [-0.5,0.5];

    subplot(1,2,2)
Pupils = imagesc(KZ_exc,KX_exc,real(Pupil_sum));
    title("Rear Pupil, I(NA1)/I(NA2)=" + num2str(WeightingRatio))
    xlabel("k_x/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colorbar
    axis image
    Pupils.Parent.XLim = [-0.5,0.5];
    Pupils.Parent.YLim = [-0.5,0.5];

%% Figure 2 - Coherent/Incoherent
    fig2 = figure;
    fig2.Name = "Incoherent/Coherent" ...
            +",beam1=" + num2str(beam1NAmax) + "/" + num2str(beam1NAmin)...
            +",beam2=" + num2str(beam2NAmax) + "/" + num2str(beam2NAmin);
    fig2.WindowState = 'maximized';
    colormap(hot(256))

    subplot(3,3,1)
    add_intensity = imagesc(X_exc,Z_exc,xzPSFIncoherent);
    title("Incoherent, Y=0")
    xlabel("x/(\lambda_{exc}/n)")
    ylabel("z /(\lambda_{exc}/n)")
    axis image
    colorbar
    add_intensity.Parent.XLim = [-20,20];
    add_intensity.Parent.YLim = [-20,20];

    subplot(3,3,2)
    add_pupil = imagesc(X_exc,Z_exc,xzPSFCoherent);
    title("Coherent, Y=0")
    xlabel("x(\lambda_{exc}/n)")
    ylabel("z(\lambda_{exc}/n)")
    axis image
    colorbar
    add_pupil.Parent.XLim = [-20,20];
    add_pupil.Parent.YLim = [-20,20];

    subplot(3,3,3)
    hold on
    add_intensity_line = plot(Z_exc,zPSFIncoherent);
        add_intensity_line.Color = 'r';
        add_intensity_line.LineWidth = 3;
%     add_pupil_line = plot(Z_exc,PSFCoherent(:,(N+1)/2,(N+1)/2));
%         add_pupil_line.Color = 'g';
%         add_pupil_line.LineWidth = 3;
    title("Focal Plane Profile,X=0,Y=0")
    xlabel("z(\lambda_{exc}/n)")
    ylabel("Normalized Intensity")
    xlim([-10,10])
    legend("Incoherent")
    grid on
    axis square

    subplot(3,3,4)
image18 = imagesc( KX_exc,...
                  KZ_exc,...
                    xzOTFIncoherent) ;
    title("Incoherent real(OTF), K_Y=0")
    xlabel("k_x/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colorbar;
    axis image
    image18.Parent.XLim = [-0.5,0.5];
    image18.Parent.YLim = [-0.5,0.5];

    subplot(3,3,5)
image19 = imagesc( KX_exc,...
                  KZ_exc,...
                    xzOTFCoherent) ;
    title("Coherent real(OTF), K_Y=0")
    xlabel("k_x/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colorbar;
    axis image
    image19.Parent.XLim = [-0.5,0.5];
    image19.Parent.YLim = [-0.5,0.5];

    subplot(3,3,6)
    hold on
    incoherent = plot( KZ_exc, zOTFIncoherent);
    incoherent.LineWidth = 3;
    incoherent.Color = 'r';
% coherent = plot( KZ_exc, coherent_OTF(:,(N+1)/2));
%     coherent.LineWidth = 2;
%     coherent.Color = 'g';
    title("Z-OTF Profile, " + "K_X=0, " + "K_Y=0")
    ylabel("Normalized a.u. ")
    xlabel("k_z/(4\pin/\lambda_{exc})")
%     incoherent.Parent.YAxis.TickValues = linspace(0,1,11);
    incoherent.Parent.XAxis.TickValues = linspace(-0.5,0.5,11);
    incoherent.Parent.XLim = [-0.5,0.5];
    incoherent.Parent.YLim = [-0.3,1];
    legend("Incoherent")
    grid on
    axis square
    hold off

    subplot(3,3,7)
    add_intensityyz = imagesc(Y_exc,Z_exc,yzPSFIncoherent);
    title("Incoherent")
    xlabel("y(\lambda_{exc}/n)")
    ylabel("z(\lambda_{exc}/n)")
    axis image
    colorbar
    add_intensityyz.Parent.YLim = [-40,40];

    subplot(3,3,8)
    add_pupilyz = imagesc(Y_exc,Z_exc,yzPSFCoherent);
    title("Coherent")
    xlabel("y(\lambda_{exc}/n)")
    ylabel("z(\lambda_{exc}/n)")
    axis image
    colorbar
    add_pupilyz.Parent.YLim = [-40,40];

    % Calculate yFWHM
    [~,maxindex] = max(yPSFIncoherent);
    index = 1-(yPSFIncoherent <= 0.5*max(yPSFIncoherent));
    if ~isempty(index)
        yFWHM1 = Y_exc(find(index,1,'first'));
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

    subplot(3,3,9)
    hold on
    add_intensity_line_yz = plot(Y_exc, yPSFIncoherent );
        add_intensity_line_yz.Color = 'r';
        add_intensity_line_yz.LineWidth = 3;
%     add_pupil_line_yz = plot(Y_exc, squeeze(PSFCoherent((N+1)/2,(N+1)/2,:)) );
%         add_pupil_line_yz.Color = 'g';
%         add_pupil_line_yz.LineWidth = 3;
    title("Z=0,X=0,Incoherent yFWHM=" + num2str(yFWHM) + "/lambda")
    xlabel("y(\lambda_{exc}/n)")
    ylabel("Normalized Intensity")
    legend("Incoherent")
    xlim([-150,150])
    ylim([0,1])
    grid on
    axis square







