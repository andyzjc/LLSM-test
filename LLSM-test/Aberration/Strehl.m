function [SW_SRatio,Lattice_SRatio,RadioOrderArray,AngularFrequencyArray] = Strehl(SWPupil,LatticePupil,MinRadialOrder,MaxRadialOrder,PhaseAmplitude,savingdir)
    getParameters; %modify image parameter here
    CalculatePhysics;

    [theta,r] = cart2pol(kx_exc./(0.6./n*k_wave),kz_exc./(0.6./n*k_wave));
    idx = r<=1;

    SWPSF = zeros(N,N);
    for j = 1:size(SWPupil,3)
        PupilIter = SWPupil(:,:,j);
        SWPSF = SWPSF + abs(fftshift(ifft2(fftshift(PupilIter)))).^2;
    end
    [SWvalue,~] = max(SWPSF,[],'all'); % value, index
    SWPSF = SWPSF./SWvalue;


    LatticePSF = abs(fftshift(ifft2(fftshift(LatticePupil)))).^2;
    LatticePSF = meshgrid(mean(LatticePSF,2))';
    [Latticevalue,~] = max(LatticePSF,[],'all'); % value, index
    LatticePSF = LatticePSF./Latticevalue;
    
    counter = 1;
    for i = MinRadialOrder:MaxRadialOrder
        AngularFrequency = -i:2:i;
        for k = 1:length(AngularFrequency)
            RadioOrderArray(counter) = i;
            AngularFrequencyArray(counter) = AngularFrequency(k);
            phase = zeros(size(kx_exc));
            phase(idx) = zernfun(i,AngularFrequency(k),r(idx),theta(idx),'norm');
    
            AberratedSWPupil = zeros(size(phase));
            AberratedSWPSF = zeros(size(phase));
            %SW 
            for j = 1:size(SWPupil,3)
                PupilIter = SWPupil(:,:,j);
                AberratedSWPupil(idx) = PupilIter(idx) .* exp(PhaseAmplitude.* 1i.*phase(idx));
                AberratedSWPSF = AberratedSWPSF + abs(fftshift(ifft2(fftshift(AberratedSWPupil)))).^2;
            end
            AberratedSWPSF = AberratedSWPSF./SWvalue; 
            SW_SRatio(counter,1) = AberratedSWPSF((N+1)/2,(N+1)/2);

            AberratedLatticePupil = zeros(size(phase));
            AberratedLatticePSF = zeros(size(phase));
            AberratedLatticePupil(idx) = LatticePupil(idx) .* exp(PhaseAmplitude.* 1i.*phase(idx));
            AberratedLatticePSF = abs(fftshift(ifft2(fftshift(AberratedLatticePupil)))).^2;
            AberratedLatticePSF = meshgrid(mean(AberratedLatticePSF,2))';
            AberratedLatticePSF = AberratedLatticePSF/Latticevalue;
            Lattice_SRatio(counter,1) = AberratedLatticePSF((N+1)/2,(N+1)/2);
            counter = counter + 1;
        end
    end

    fig1 = figure;
    h1 = subplot(1,1,1,'Parent',fig1);
    b1 = bar(1:length(SW_SRatio),[SW_SRatio,Lattice_SRatio],0.5,'Parent',h1);
    lgd = legend(b1,"iSW","LLS");
    h1.XAxis.TickValues = 1:length(RadioOrderArray);
    LabelArray = [RadioOrderArray;AngularFrequencyArray];
    tickLabels = strtrim(sprintf('%d\\newline%d\n', LabelArray(:)));
    h1.XAxis.TickLabels = tickLabels;
    grid on
