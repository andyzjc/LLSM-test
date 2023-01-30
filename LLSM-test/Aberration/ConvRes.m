function Spacing = ConvRes(PSFexc,deconPSFexc)
    getParameters; %modify image parameter here
    CalculatePhysics;

    lineSpot = [0.0000,0.2202,0.4770,0.7705,1.1008,1.4677,1.8713,2.3116,...
                2.7886,3.3023,3.8526,4.4397,5.0635,5.7239,6.4211,7.1549,...
                7.9254,8.7326,9.5765,10.4571,14.0896]; %um 
%     Spacing = 0.22:0.01:0.9;
%     lineSpot = zeros(1,length(Spacing)+1);
%     for i = 1:length(Spacing)
%         lineSpot(1,i+1) = Spacing(1,i) + lineSpot(1,i);
%     end

    % define line array 
    N_line = 2049;
    deltax_line = 0.01; %um
    GTLines = zeros(N_line,N_line);
    GTLines(:,round(lineSpot/deltax_line)+1) = 1;
    Spacing = diff(round(lineSpot/deltax_line) * deltax_line * 1000); %nm
    z_exc_line = (0:N_line-1) * deltax_line;
    Z_exc_line = z_exc_line / wavelength_exc;
    
    % get overall PSF
    xzPSFOverall = getOverallPSF(PSFexc(:,:,(N+1)/2));
    rescale_factor = deltax/deltax_line;
    up_Image_size = round(size(xzPSFOverall,1) * rescale_factor);
    Image_center = (up_Image_size+1)/2;
    scaledxzPSFOverall = imresize(xzPSFOverall,[up_Image_size,up_Image_size]);
    scaledxzPSFOverall = scaledxzPSFOverall/max(max(scaledxzPSFOverall));
    scaledxzPSFOverall = scaledxzPSFOverall(Image_center-(N_line+1)/2+1:Image_center+(N_line+1)/2-1,...
                                  Image_center-(N_line+1)/2+1:Image_center+(N_line+1)/2-1);

   % get overall deconPSF
    deconxzPSFOverall = getOverallPSF(deconPSFexc(:,:,(N+1)/2));
    rescale_factor = deltax/deltax_line;
    up_Image_size = round(size(deconxzPSFOverall,1) * rescale_factor);
    Image_center = (up_Image_size+1)/2;
    scaleddeconxzPSFOverall = imresize(deconxzPSFOverall,[up_Image_size,up_Image_size]);
    scaleddeconxzPSFOverall = scaleddeconxzPSFOverall/max(max(scaleddeconxzPSFOverall));
    scaleddeconxzPSFOverall = scaleddeconxzPSFOverall(Image_center-(N_line+1)/2+1:Image_center+(N_line+1)/2-1,...
                                  Image_center-(N_line+1)/2+1:Image_center+(N_line+1)/2-1);

    % Convolution with GT lines
    convLines = conv2(GTLines,scaledxzPSFOverall','same');
    convLines = convLines/max(max(convLines));

%     % deconvlution 
    DeconvLines = deconvlucy(convLines,scaleddeconxzPSFOverall,20);
    DeconvLines = DeconvLines/max(max(DeconvLines));

    % Plotting 
    plot(z_exc_line,GTLines((N_line+1)/2,:)/max(GTLines((N_line+1)/2,:)),'k')
    hold on
    grid on
    plot(z_exc_line,convLines((N_line+1)/2,:)/max(convLines((N_line+1)/2,:)),'r')
     plot(z_exc_line,DeconvLines((N_line+1)/2,:)/max(DeconvLines((N_line+1)/2,:)),'b')
     xlim([0,20])
    legend("GT","Conv","Deconv")
    xlabel("z(um)")
    hold off

