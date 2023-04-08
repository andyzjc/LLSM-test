function [convLines,lineSpot,Spacing,LineZ] = ConvRes(PSFexc,PSFdet)
    getParameters; %modify image parameter here
    CalculatePhysics;

    % lineSpot = [0.0000,0.2202,0.4770,0.7705,1.1008,1.4677,1.8713,2.3116,...
    %             2.7886,3.3023,3.8526,4.4397,5.0635,5.7239,6.4211,7.1549,...
    %             7.9254,8.7326,9.5765,10.4571,14.0896]; %um 
    lineSpot = [1 3 6 10 15 21 28 36 45 55 66 78 91 105 120 160];

    % define line array 
    GTLines = zeros(N,N);
    GTLines(:,lineSpot) = 1;
    Spacing = diff(lineSpot * deltax * 1000); %nm
    LineZ = (0:N-1) * deltax;
    
    % get overall PSF
    PSFoverall = PSFexc .* PSFdet;
    xzPSFOveralldecon = PSFoverall(:,:,(N+1)/2);
    % xzPSFOverall = xzPSFOveralldecon + poissrnd(xzPSFOveralldecon) .* 1/SNR;

    % Convolution with GT lines
    convLines = conv2(GTLines,xzPSFOveralldecon','same');
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

