function [convLines,lineSpot,Spacing,LineZ] = ConvRes(PSFexc,PSFdet,SNR)
    getParameters; %modify image parameter here
    CalculatePhysics;
    
    % lineSpot = [0.0000,0.2202,0.4770,0.7705,1.1008,1.4677,1.8713,2.3116,...
    %             2.7886,3.3023,3.8526,4.4397,5.0635,5.7239,6.4211,7.1549,...
    %             7.9254,8.7326,9.5765,10.4571,14.0896]; %um 
    N_line = 1025;
    deltax_line = 0.02; %um
    Spacing = 0.2:0.02:0.9; % um
    lineSpot = zeros(1,length(Spacing)+1);
    lineSpot(1,1) = 1;
    for i = 1:length(Spacing)
        lineSpot(1,i+1) = round(( Spacing(1,i))./deltax_line) + lineSpot(1,i);
    end

    % define line array 
    GTLines = zeros(N_line,N_line);
    GTLines(:,lineSpot) = 1;
    LineZ = (0:N_line-1) * deltax_line;
    LineZ = LineZ / wavelength_exc;
    
    % get overall PSF
    PSFoverall = PSFexc .* PSFdet;
    xzPSFOveralldecon = PSFoverall(:,:,(N+1)/2);
    %xzPSFOveralldecon = xzPSFOveralldecon + poissrnd(xzPSFOveralldecon) .* 1/SNR;
    
    % interpolate overall PSF
     rescale_factor = deltax/deltax_line;
    up_Image_size = round(size(xzPSFOveralldecon,1) * rescale_factor);
    Image_center = (up_Image_size+1)/2;
    scaledxzPSFOverall = imresize(xzPSFOveralldecon,[up_Image_size,up_Image_size]);
    scaledxzPSFOverall = scaledxzPSFOverall(Image_center-(N_line+1)/2+1:Image_center+(N_line+1)/2-1,...
                                  Image_center-(N_line+1)/2+1:Image_center+(N_line+1)/2-1);

   %scaledxzPSFOverall = imresize(xzPSFOveralldecon,[N_line,N_line]);
   %scaledxzPSFOverall = scaledxzPSFOverall + poissrnd(scaledxzPSFOverall) .* 1/SNR;

    % Convolution with GT lines
    convLines = conv2(GTLines,scaledxzPSFOverall','same');
    convLines = convLines + poissrnd(convLines) .* 1/SNR;
    % convLines = convLines/max(max(convLines));

    imagesc(convLines);
    % deconvlution 
    % DeconvLines = deconvlucy(convLines,xzPSFOveralldecon,deconIter);
    % DeconvLines = DeconvLines/max(max(DeconvLines));

    % Plotting 
    % plot(z_exc_line,GTLines((N+1)/2,:)/max(GTLines((N+1)/2,:)),'k')
    % hold on
    % grid on
    % plot(z_exc_line,convLines((N+1)/2,:)/max(convLines((N+1)/2,:)),'r')
    % % plot(z_exc_line,DeconvLines((N+1)/2,:)/max(DeconvLines((N+1)/2,:)),'b')
    % xlim([0,20])
    % legend("GT","Conv","deConv")
    % xlabel("z(um)")
    % hold off

