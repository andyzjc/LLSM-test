clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

savingdir = ['NaJi_Zebrafish_motor_neuron_Aberration_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_' ProfileType '_SWWeighting_' num2str(SWweighting) '_LLSWeighting_' num2str(Latticeweighting) '_SNR_' num2str(SNR) '/'];

%% parameters 
NA1 = (0.3:0.05:0.55)';
deltaNA1 = 0.02:0.01:0.2; 
sinc_deltaNA2 = 0.1:0.01:0.4; 
weighting = 1;

NA1_diff_sum = [];
NA_Ratio = [];
for k = 1:length(NA1)
    NA1_Iter = NA1(i)

    Square_yFWHM = [];
    for i = 1:length(deltaNA1)
        Square_deltaNA2 = sqrt(2*NA1_Iter*deltaNA1(i));
        NA1min = NA1_Iter - deltaNA1(i)/2;
        NA1max = NA1_Iter + deltaNA1(i)/2;
        if Square_deltaNA2 < NA1min
            Square_NA2 = Square_deltaNA2/2;
            [SWPupil,~,~] = GetSWPairPupil('tophat',NA1_Iter,Square_NA2,...
                                                    deltaNA1(i),Square_deltaNA2,...
                                                    weighting);
            [~,PSFIncoherent,~] = SimulateSWPair(SWPupil);
            PSFIncoherent = PSFIncoherent/max(PSFIncoherent,[],'all');
            yindex = 1-(squeeze(PSFIncoherent((N+1)/2,(N+1)/2,:)) <= 0.5*max(squeeze(PSFIncoherent((N+1)/2,(N+1)/2,:))));
            SWyFWHM1 = find(yindex,1,'first') ;
            Square_yFWHM(i) = abs(Y_exc(SWyFWHM1)).*2;
        end
    end
    
    Sinc_yFWHM = [];
    for i = 1:length(sinc_deltaNA2)
        k_apertureNA = sinc_deltaNA2(i) * k_wave / n;
        sincPupil = zeros(N,N);
        sincPupil(:,(N+1)/2) = (k_apertureNA) >= abs(kz_exc(:,1));
        [~,sincPSFDithered,~] = SimulateLattice(sincPupil);
        yindex = 1-(squeeze(sincPSFDithered((N+1)/2,(N+1)/2,:)) <= 0.5*max(squeeze(sincPSFDithered((N+1)/2,(N+1)/2,:))));
        SincyFWHM1index = find(yindex,1,'first');
        Sinc_yFWHM(i) = abs(Y_exc(SincyFWHM1index)).*2;
    end
    
    beam1minyFWHM = min(Square_yFWHM);
    beam1maxyFWHM = max(Square_yFWHM);
    beam2minyFWHM = min(Sinc_yFWHM);
    beam2maxyFWHM = max(Sinc_yFWHM);
    if beam1minyFWHM >= beam2minyFWHM
        minyFWHM_equal = beam1minyFWHM;
    else
        minyFWHM_equal = beam2minyFWHM;
    end
    if beam1maxyFWHM >= beam2maxyFWHM
        maxyFWHM_equal = beam2maxyFWHM;
    else
        maxyFWHM_equal = beam1maxyFWHM;
    end
    yFWHM_equal = (minyFWHM_equal:0.1*wavelength_exc:maxyFWHM_equal)';
    
    % find intersection
    deltaNA_Square = NaN(size(yFWHM_equal));
    deltaNA_sinc = deltaNA_Square;
    for jj = 2:length(yFWHM_equal)
        [x1,~,~,~] = intersections(deltaNA1,Square_yFWHM,...
            0:0.01:0.4,yFWHM_equal(jj)*ones(1,length(0:0.01:0.4)),0);
        [x2,~,~,~] = intersections(sinc_deltaNA2,Sinc_yFWHM,...
            0:0.01:0.4,yFWHM_equal(jj)*ones(1,length(0:0.01:0.4)),0);
        deltaNA_Square(jj,1) = x1(1);
        deltaNA_sinc(jj,1) = x2(1);
    end

    % analytical curve 
    analytical_X = [];
    analytical_Y = [];
    analytical_X = deltaNA_Square ./ (2*NA1_Iter);
    analytical_Y = (NA1_Iter + (deltaNA_Square)./2 ) ./ (deltaNA_sinc);

    NA1_diff_sum = [NA1_diff_sum analytical_X];
    NA_Ratio = [NA_Ratio analytical_Y];

    fig1 = figure;
    scatter(deltaNA1,sqrt(2.*deltaNA1*NA1_Iter)) % approximation
    hold on
    A = sqrt( n.^2- (NA1_Iter - (deltaNA1/2) ).^2 );
    B = sqrt( n.^2- (NA1_Iter + (deltaNA1/2) ).^2 );
    scatter(deltaNA1, sqrt( n.^2- (n-A+B).^2 )) % strict non-diffractive
    plot(deltaNA_Square,deltaNA_sinc) % numerical simulation
    legend("Approximation","Strict","Numerical")
    xlabel("deltaNA Sqaure")
    ylabel("deltaNA sinc")
    grid on
    hold off
    print(fig1, '-dsvg', [  savingdir 'EquivalentdeltaNA_NA1_' num2str(NA1_Iter) '.SVG'],'-r300')
    print(fig1, '-dpng', [  savingdir 'EquivalentdeltaNA_NA1_' num2str(NA1_Iter) '.PNG'],'-r300') 

    fig2 = figure;
    plot(deltaNA1,Square_yFWHM)
    hold on
    plot(sinc_deltaNA2,Sinc_yFWHM)
    grid on
    legend("Square","sinc")
    print(fig2, '-dsvg', [  savingdir 'yFWHM_NA1_' num2str(NA1_Iter) '.SVG'],'-r300')
    print(fig2, '-dpng', [  savingdir 'yFWHM_NA1_' num2str(NA1_Iter) '.PNG'],'-r300') 
    close all
end

%% Analytical curve
x = 0.01:0.01:1;
y = 1/2*(sqrt(x)+sqrt(1./x));
fig3 = figure;
    plot(x,y)
    grid on
    hold on
    plot(NA1_diff_sum,NA_Ratio)
    legend("Analytical","Numerical")


