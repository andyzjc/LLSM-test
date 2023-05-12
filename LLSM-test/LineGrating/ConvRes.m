function [convLines,lineSpot,Spacing,LineZ] = ConvRes(PSFexc,PSFdet,SNR)
    getParameters; %modify image parameter here
    CalculatePhysics;
    
    % lineSpot = [0.0000,0.2202,0.4770,0.7705,1.1008,1.4677,1.8713,2.3116,...
    %             2.7886,3.3023,3.8526,4.4397,5.0635,5.7239,6.4211,7.1549,...
    %             7.9254,8.7326,9.5765,10.4571,14.0896]; %um 
    N_line = 2049;
    deltax_line = 0.01; %um
    Spacing = 0.22:0.03:0.9; % um
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
    % xzPSFOverall = xzPSFOveralldecon + poissrnd(xzPSFOveralldecon) .* 1/SNR;
    
    % interpolate overall PSF
   scaledxzPSFOverall = imresize(scaledxzPSFOverall,[N_line,N_line]);
   scaledxzPSFOverall = scaledxzPSFOverall + poissrnd(scaledxzPSFOverall) .* 1/SNR;

    % Convolution with GT lines
    convLines = conv2(GTLines,scaledxzPSFOverall','same');
    % convLines = convLines/max(max(convLines));

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

