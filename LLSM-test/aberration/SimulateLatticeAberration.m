function [CorrCoef,StrehlPeaks,StrehlCenter] = SimulateLatticeAberration(LatticePupil,PSFdet,MaxRadialOrder,PhaseAmplitude)
    getParameters; %modify image parameter here
    CalculatePhysics;

    % Unaberrated
    [LatticePSF,LatticePSFDithered,Latticecenter] = SimulateLattice(LatticePupil);
    LatticePSF = LatticePSF/max(LatticePSF,[],'all');
    LatticePSFDithered = LatticePSFDithered/max(LatticePSFDithered,[],'all');

    UnaberratedxzPSF = LatticePSFDithered(:,:,(N+1)/2); 
    UnaberratedzPSF = UnaberratedxzPSF(:,(N+1)/2); 
    UnaberratedxzOTF = fftshift(fft2(ifftshift(UnaberratedxzPSF))); UnaberratedxzOTF = UnaberratedxzOTF/max(max(UnaberratedxzOTF));
    UnaberratedzOTF = UnaberratedxzOTF(:,(N+1)/2); 
    UnaberratedyzPSF = squeeze(LatticePSFDithered(:,(N+1)/2,:));
    UnaberratedyPSF = UnaberratedyzPSF((N+1)/2,:);
    
    % Overall Unaberrated 
    UnaberratedOverallPSF = LatticePSFDithered .* PSFdet;
%     UnaberratedOverallPSF = UnaberratedOverallPSF/max(UnaberratedOverallPSF,[],'all');
    UnaberratedOverallOTF = fftshift(fftn(ifftshift(UnaberratedOverallPSF)));
    UnaberratedOverallOTF = UnaberratedOverallOTF/max(UnaberratedOverallOTF,[],'all');

    UnaberratedOverallxzPSF = UnaberratedOverallPSF(:,:,(N+1)/2); 
    UnaberratedOverallzPSF = UnaberratedOverallxzPSF(:,(N+1)/2); 
    UnaberratedOverallxzOTF = UnaberratedOverallOTF(:,:,(N+1)/2); 
    UnaberratedOverallzOTF = UnaberratedOverallxzOTF(:,(N+1)/2); 

    [theta,r] = cart2pol(kx_exc./(0.6./n*k_wave),kz_exc./(0.6./n*k_wave));
    idx = r<=1;
    phase = zeros(size(kx_exc));

    fig1 = figure(1);
        fig1.Name = 'Excitation xzPSF';
        fig1.WindowState = 'maximized';
    fig2 = figure(2);
        fig2.Name = 'Excitation zPSF';
        fig2.WindowState = 'maximized';
    fig3 = figure(3);
        fig3.Name = 'Excitation xzOTF';
        fig3.WindowState = 'maximized';
    fig4 = figure(4);
        fig4.Name = 'Excitation zOTF';
        fig4.WindowState = 'maximized';
    fig5 = figure(5);
        fig5.Name = 'Excitation yzPSF';
        fig5.WindowState = 'maximized';
    fig6 = figure(6);
        fig6.Name = 'Excitation yPSF';
        fig6.WindowState = 'maximized';
    fig7 = figure(7);
        fig7.Name = 'Overall xzPSF';
        fig7.WindowState = 'maximized';
    fig8 = figure(8);
        fig8.Name = 'Overall zPSF';
        fig8.WindowState = 'maximized';
    fig9 = figure(9);
        fig9.Name = 'Overall xzOTF';
        fig9.WindowState = 'maximized';
    fig10 = figure(10);
        fig10.Name = 'Overall zOTF';
        fig10.WindowState = 'maximized';
    fig12 = figure(12);
        fig12.Name = 'Excitation xzPSF NonDither';
        fig12.WindowState = 'maximized';
    fig13 = figure(13);
        fig13.Name = 'Metric';
        fig13.WindowState = 'maximized';

    MinRadialOrder = 0;
    counter = 1;
    for i = MinRadialOrder:MaxRadialOrder
        RadialOrder = i*ones(1,i+1);
        AngularFrequency = -i:2:i;
        AngularFrequency_iteration = 1:1:length(AngularFrequency);
        for k = 1:length(AngularFrequency)
            RadioOrderArray(counter) = i;
            AngularFrequencyArray(counter) = AngularFrequency(k);

            phase(idx) = zernfun(i,AngularFrequency(k),r(idx),theta(idx),'norm');

            AberratedPupil = zeros(size(phase));
            AberratedPupil(idx) = LatticePupil(idx) .* exp(PhaseAmplitude.* 1i .*phase(idx));
            
            AberratedLatticePSF = zeros(N,N, N);
            AberratedLatticePSFDithered = zeros(N,N, N);

            % propagation
            for j = 1:length(y_exc)
                propagator_exc = exp(2*pi * 1i * ky_exc * y_exc(j));
                AberratedLatticePSF(:,:,j) = abs( fftshift( ifft2(ifftshift(AberratedPupil .* propagator_exc)) ) ).^2;
                AberratedLatticePSFDithered(:,:,j) = meshgrid(mean(AberratedLatticePSF(:,:,j),2))';
            end  

            % Normalize to unaberrated center 
            AberratedLatticePSF = AberratedLatticePSF/Latticecenter(1,1);
            AberratedLatticePSFDithered = AberratedLatticePSFDithered/Latticecenter(2,1);

            xzPSFexc = AberratedLatticePSFDithered(:,:,(N+1)/2); 
            zPSFexc = xzPSFexc(:,(N+1)/2); 
            xzOTFexc = fftshift(fft2(ifftshift(xzPSFexc))); xzOTFexc = xzOTFexc/max(max(xzOTFexc));
            zOTFexc = xzOTFexc(:,(N+1)/2); 
            yzPSFexc = squeeze(AberratedLatticePSFDithered(:,(N+1)/2,:));
            yPSFexc = yzPSFexc((N+1)/2,:); 

            OverallPSF = AberratedLatticePSFDithered .* PSFdet;
%             OverallPSF = OverallPSF/max(OverallPSF,[],'all');
            OverallOTF = fftshift(fftn(ifftshift(OverallPSF)));
            OverallOTF = OverallOTF/max(OverallOTF,[],'all');
        
            OverallxzPSF = OverallPSF(:,:,(N+1)/2); 
            OverallzPSF = OverallxzPSF(:,(N+1)/2); 
            OverallxzOTF = OverallOTF(:,:,(N+1)/2); 
            OverallzOTF = OverallxzOTF(:,(N+1)/2); 

            CoefMatrix = corrcoef(real(UnaberratedOverallzOTF),real(zOTFexc));
            CorrCoef(counter) = CoefMatrix(1,2);
            StrehlPeaks(counter) = max(AberratedLatticePSFDithered,[],'all');
            StrehlCenter(counter) = AberratedLatticePSFDithered((N+1)/2,(N+1)/2,(N+1)/2)./ LatticePSFDithered((N+1)/2,(N+1)/2,(N+1)/2);

            counter = counter + 1;

            h1 = subplot(MaxRadialOrder-MinRadialOrder+1,MaxRadialOrder+1,((MaxRadialOrder+1)*(i-MinRadialOrder))+AngularFrequency_iteration(k),'Parent',fig1);
                imagesc(h1,X_exc,Z_exc,xzPSFexc);
                h1.XAxis.Label.String = "x(\lambda_{exc}/n)";
                h1.YAxis.Label.String = "z(\lambda_{exc}/n)";
                h1.Colormap = colormap(hot(256));
                h1.XAxis.Limits = [-20,20];
                h1.YAxis.Limits = [-20,20];
                h1.Title.String = ['Z_{' num2str(RadialOrder(k)) '}^{' num2str(AngularFrequency(k)) '}'];
                h1.DataAspectRatio = [1,1,1];
                colorbar(h1)
                clim(h1,[0,1])


            h2 = subplot(MaxRadialOrder-MinRadialOrder+1,MaxRadialOrder+1,((MaxRadialOrder+1)*(i-MinRadialOrder))+AngularFrequency_iteration(k),'Parent',fig2);
                h2.NextPlot = "add";
                plot(h2,Z_exc,zPSFexc,'LineWidth',2,'Color','r')
                plot(h2,Z_exc,UnaberratedzPSF,'LineWidth',1,'Color','g')
                h2.XAxis.Label.String = "z(\lambda_{exc}/n)";
                h2.XAxis.Limits = [-20,20];
                h2.Title.String = ['Z_{' num2str(RadialOrder(k)) '}^{' num2str(AngularFrequency(k)) '}'];
                h2.XGrid = 'on';
                h2.YGrid = 'on';
                lgd = legend(h2,"Aberrated","Normal");
                lgd.FontSize = 3;
                grid on
                h2.NextPlot = "replace";

            h3 = subplot(MaxRadialOrder-MinRadialOrder+1,MaxRadialOrder+1,((MaxRadialOrder+1)*(i-MinRadialOrder))+AngularFrequency_iteration(k),'Parent',fig3);
                imagesc(h3,KX_exc,KZ_exc,real(xzOTFexc))
                h3.XAxis.Label.String = "k_x/(4\pin/\lambda_{exc})";
                h3.YAxis.Label.String = "k_z/(4\pin/\lambda_{exc})";
                h3.Colormap = colormap(hot(256));
                h3.XAxis.Limits = [-0.5,0.5];
                h3.YAxis.Limits = [-0.5,0.5];
                h3.Title.String = ['Z_{' num2str(RadialOrder(k)) '}^{' num2str(AngularFrequency(k)) '}'];
                h3.DataAspectRatio = [1,1,1];
                colorbar(h3)
                clim(h3,[0,1])

            h4 = subplot(MaxRadialOrder-MinRadialOrder+1,MaxRadialOrder+1,((MaxRadialOrder+1)*(i-MinRadialOrder))+AngularFrequency_iteration(k),'Parent',fig4);
                h4.NextPlot = "add";
                plot(h4,KZ_exc,real(zOTFexc),'LineWidth',2,'Color','r')
                plot(h4,KZ_exc,abs(zOTFexc),'LineWidth',1.5,'Color','b')
                plot(h4,KZ_exc,real(UnaberratedzOTF),'LineWidth',1,'Color','g')
                h4.XAxis.Label.String = "k_z/(4\pin/\lambda_{exc})";
                h4.XAxis.Limits = [-0.5,0.5];
                h4.YAxis.Limits = [-0.5,1];
                h4.Title.String = ['Z_{' num2str(RadialOrder(k)) '}^{' num2str(AngularFrequency(k)) '}'];
                h4.XAxis.TickValues = linspace(-0.5,0.5,11);
                h4.XGrid = 'on';
                h4.YGrid = 'on';
                lgd = legend(h4,"real(Aberrated)","abs(Aberrated)","Normal");
                lgd.FontSize = 3;
                grid on
                h4.NextPlot = "replace";

            h5 = subplot(MaxRadialOrder-MinRadialOrder+1,MaxRadialOrder+1,((MaxRadialOrder+1)*(i-MinRadialOrder))+AngularFrequency_iteration(k),'Parent',fig5);
                imagesc(h5,Y_exc,Z_exc,yzPSFexc)
                h5.XAxis.Label.String = "y(\lambda_{exc}/n)";
                h5.YAxis.Label.String = "z(\lambda_{exc}/n)";
                h5.Colormap = colormap(hot(256));
                h5.YAxis.Limits = [-40,40];
                h5.Title.String = ['Z_{' num2str(RadialOrder(k)) '}^{' num2str(AngularFrequency(k)) '}'];
                h5.DataAspectRatio = [1,1,1];
                colorbar(h5)
                clim(h5,[0,1])

            h6 = subplot(MaxRadialOrder-MinRadialOrder+1,MaxRadialOrder+1,((MaxRadialOrder+1)*(i-MinRadialOrder))+AngularFrequency_iteration(k),'Parent',fig6);
                h6.NextPlot = "add";
                plot(h6,Y_exc,yPSFexc,'LineWidth',2,'Color','r')
                plot(h6,Y_exc,UnaberratedyPSF,'LineWidth',1,'Color','g')
                h6.XAxis.Label.String = "y(\lambda_{exc}/n)";
                h6.Title.String = ['Z_{' num2str(RadialOrder(k)) '}^{' num2str(AngularFrequency(k)) '}'];
                h6.XGrid = 'on';
                h6.YGrid = 'on';
                lgd = legend(h6,"Aberrated","Normal");
                lgd.FontSize = 3;
                grid on
                h6.NextPlot = "replace";

            h7 = subplot(MaxRadialOrder-MinRadialOrder+1,MaxRadialOrder+1,((MaxRadialOrder+1)*(i-MinRadialOrder))+AngularFrequency_iteration(k),'Parent',fig7);
                imagesc(h7,X_exc,Z_exc,OverallxzPSF);
                h7.XAxis.Label.String = "x(\lambda_{exc}/n)";
                h7.YAxis.Label.String = "z(\lambda_{exc}/n)";
                h7.Colormap = colormap(hot(256));
                h7.XAxis.Limits = [-20,20];
                h7.YAxis.Limits = [-20,20];
                h7.Title.String = ['Z_{' num2str(RadialOrder(k)) '}^{' num2str(AngularFrequency(k)) '}'];
                h7.DataAspectRatio = [1,1,1];
                colorbar(h7)
                clim(h7,[0,1])

            h8 = subplot(MaxRadialOrder-MinRadialOrder+1,MaxRadialOrder+1,((MaxRadialOrder+1)*(i-MinRadialOrder))+AngularFrequency_iteration(k),'Parent',fig8);
                h8.NextPlot = "add";
                plot(h8,Z_exc,OverallzPSF,'LineWidth',2,'Color','r')
                plot(h8,Z_exc,UnaberratedOverallzPSF,'LineWidth',1,'Color','g')
                h8.XAxis.Label.String = "z(\lambda_{exc}/n)";
                h8.XAxis.Limits = [-20,20];
                h8.Title.String = ['Z_{' num2str(RadialOrder(k)) '}^{' num2str(AngularFrequency(k)) '}'];
                h8.XGrid = 'on';
                h8.YGrid = 'on';
                lgd = legend(h8,"Aberrated","Normal");
                lgd.FontSize = 3;
                grid on
                h8.NextPlot = "replace";

            h9 = subplot(MaxRadialOrder-MinRadialOrder+1,MaxRadialOrder+1,((MaxRadialOrder+1)*(i-MinRadialOrder))+AngularFrequency_iteration(k),'Parent',fig9);
                imagesc(h9,KX_exc,KZ_exc,real(OverallxzOTF))
                h9.XAxis.Label.String = "k_x/(4\pin/\lambda_{exc})";
                h9.YAxis.Label.String = "k_z/(4\pin/\lambda_{exc})";
                h9.Colormap = colormap(hot(256));
                h9.XAxis.Limits = [-0.5,0.5];
                h9.YAxis.Limits = [-0.5,0.5];
                h9.Title.String = ['Z_{' num2str(RadialOrder(k)) '}^{' num2str(AngularFrequency(k)) '}'];
                h9.DataAspectRatio = [1,1,1];
                colorbar(h9)
                clim(h9,[0,1])

             h10 = subplot(MaxRadialOrder-MinRadialOrder+1,MaxRadialOrder+1,((MaxRadialOrder+1)*(i-MinRadialOrder))+AngularFrequency_iteration(k),'Parent',fig10);
                h10.NextPlot = "add";
                plot(h10,KZ_exc,real(OverallzOTF),'LineWidth',2,'Color','r')
                plot(h10,KZ_exc,abs(OverallzOTF),'LineWidth',1.5,'Color','b')
                plot(h10,KZ_exc,real(UnaberratedOverallzOTF),'LineWidth',1,'Color','g')
                h10.XAxis.Label.String = "k_z/(4\pin/\lambda_{exc})";
                h10.XAxis.Limits = [-0.5,0.5];
                h10.YAxis.Limits = [-0.5,1];
                h10.XAxis.TickValues = linspace(-0.5,0.5,11);
                h10.Title.String = ['Z_{' num2str(RadialOrder(k)) '}^{' num2str(AngularFrequency(k)) '}'];
                h10.XGrid = 'on';
                h10.YGrid = 'on';
                lgd = legend(h10,"real(Aberrated)","abs(Aberrated)","Normal");
                lgd.FontSize = 3;
                grid on
                h10.NextPlot = "replace";           
                
            h12 = subplot(MaxRadialOrder-MinRadialOrder+1,MaxRadialOrder+1,((MaxRadialOrder+1)*(i-MinRadialOrder))+AngularFrequency_iteration(k),'Parent',fig12);
                imagesc(h12,X_exc,Z_exc,AberratedLatticePSF(:,:,(N+1)/2));
                h12.XAxis.Label.String = "x(\lambda_{exc}/n)";
                h12.YAxis.Label.String = "z(\lambda_{exc}/n)";
                h12.Colormap = colormap(hot(256));
                h12.XAxis.Limits = [-20,20];
                h12.YAxis.Limits = [-20,20];
                h12.Title.String = ['Z_{' num2str(RadialOrder(k)) '}^{' num2str(AngularFrequency(k)) '}'];
                h12.DataAspectRatio = [1,1,1];
                colorbar(h12)
                clim(h12,[0,1])
                
            for jj = (N+1)/2-100:10:(N+1)/2+100
                fig11 = figure('Name','Profile','WindowState','maximized','Visible','off');
                colormap(hot(256))

                xzPSF_exc = AberratedLatticePSF(:,:,jj); 
                xzPSF_exc_dither = AberratedLatticePSFDithered(:,:,jj); 
                zPSF_exc = xzPSF_exc(:,(N+1)/2); 
                zPSF_exc_dither = xzPSF_exc_dither(:,(N+1)/2); 
                xzOTF_exc = fftshift(fft2(ifftshift(xzPSF_exc))); xzOTF_exc = xzOTF_exc/max(max(xzOTF_exc));
                zOTF_exc = xzOTF_exc(:,(N+1)/2);

                UnaberratedxzPSF_exc = LatticePSF(:,:,jj); 
                UnaberratedxzPSF_exc_dither = LatticePSFDithered(:,:,jj); 
                UnaberratedzPSF_exc = UnaberratedxzPSF_exc(:,(N+1)/2); 
                UnaberratedzPSF_exc_dither = UnaberratedxzPSF_exc_dither(:,(N+1)/2); 
                UnaberratedxzOTF_exc = fftshift(fft2(ifftshift(UnaberratedxzPSF_exc))); UnaberratedxzOTF_exc = UnaberratedxzOTF_exc/max(max(UnaberratedxzOTF_exc));
                UnaberratedzOTF_exc = UnaberratedxzOTF_exc(:,(N+1)/2);

                h1_Unaberrated = subplot(3,3,1,'Parent',fig11);
                h1_image_unaberrated = imagesc(h1_Unaberrated,X_exc,Z_exc,UnaberratedxzPSF_exc);
                title("XZ-Excitation PSF - Y=" + num2str(Y_exc(jj),'%.2f') + "\lambda / n")
                h1_Unaberrated.Title.FontSize = 5;
                xlabel("x/(\lambda_{exc}/n)")
                ylabel("z/(\lambda_{exc}/n)")
                h1_Unaberrated.XAxis.FontSize = 5;
                h1_Unaberrated.YAxis.FontSize = 5;
                h1_Unaberrated.XAxis.FontWeight = 'bold';
                h1_Unaberrated.YAxis.FontWeight = 'bold';
                colorbar;
                clim([0,1])
                axis image;
                h1_image_unaberrated.Parent.XLim = [-40,40];
                h1_image_unaberrated.Parent.YLim = [-40,40];
                
                h1 = subplot(3,3,2,'Parent',fig11);
                h1_image = imagesc(h1,X_exc,Z_exc,xzPSF_exc);
                title("Aberrated XZ-Excitation PSF - Y=" + num2str(Y_exc(jj),'%.2f') + "\lambda / n")
                h1.Title.FontSize = 5;
                xlabel("x/(\lambda_{exc}/n)")
                ylabel("z/(\lambda_{exc}/n)")
                h1.XAxis.FontSize = 5;
                h1.YAxis.FontSize = 5;
                h1.XAxis.FontWeight = 'bold';
                h1.YAxis.FontWeight = 'bold';
                colorbar;
                clim([0,1])
                axis image;
                h1_image.Parent.XLim = [-40,40];
                h1_image.Parent.YLim = [-40,40];
        
                h2 = subplot(3,3,3,'Parent',fig11);
                h2.NextPlot = "add";
                h2_image = plot(h2,Z_exc,zPSF_exc);
                h2_image_Unaberrated = plot(h2,Z_exc,UnaberratedzPSF_exc);
                title("Z-Excitation PSF - X=0" + "\lambda / n")
                h2.Title.FontSize = 5;
                ylabel("Normalized a.u. ")
                xlabel("z/(\lambda_{exc}/n)")
                h2.XAxis.FontSize = 4;
                h2.YAxis.FontSize = 5;
                h2.XAxis.FontWeight = 'bold';
                h2.YAxis.FontWeight = 'bold';
                h2_image.LineWidth = 2;
                h2_image.Color = 'r';
                h2_image.Parent.XLim = [-20,20];
                h2_image.Parent.YAxis.TickValues = linspace(0,1,11);
                h2_image.Parent.XAxis.TickValues = linspace(-20,20,21);
                h2_image_Unaberrated.LineWidth = 1;
                h2_image_Unaberrated.Color = 'g';
                lgd = legend(h2,"Aberrated","Normal");
                lgd.FontSize = 3;
                grid on
                axis square
                h2.NextPlot = "replace";

                h1_dither_unaberrated = subplot(3,3,4,'Parent',fig11);
                h1_image_dither_unaberrated = imagesc(h1_dither_unaberrated,X_exc,Z_exc,UnaberratedxzPSF_exc_dither);
                title("Dithered XZ-Excitation PSF - Y=" + num2str(Y_exc(jj),'%.2f') + "\lambda / n")
                h1_dither_unaberrated.Title.FontSize = 5;
                xlabel("x/(\lambda_{exc}/n)")
                ylabel("z/(\lambda_{exc}/n)")
                h1_dither_unaberrated.XAxis.FontSize = 5;
                h1_dither_unaberrated.YAxis.FontSize = 5;
                h1_dither_unaberrated.XAxis.FontWeight = 'bold';
                h1_dither_unaberrated.YAxis.FontWeight = 'bold';
                colorbar;
                clim([0,1])
                axis image;
                h1_image_dither_unaberrated.Parent.XLim = [-40,40];
                h1_image_dither_unaberrated.Parent.YLim = [-40,40];

                h1_dither = subplot(3,3,5,'Parent',fig11);
                h1_image_dither = imagesc(h1_dither,X_exc,Z_exc,xzPSF_exc_dither);
                title("Dithered Aberrated XZ-Excitation PSF - Y=" + num2str(Y_exc(jj),'%.2f') + "\lambda / n")
                h1_dither.Title.FontSize = 5;
                xlabel("x/(\lambda_{exc}/n)")
                ylabel("z/(\lambda_{exc}/n)")
                h1_dither.XAxis.FontSize = 5;
                h1_dither.YAxis.FontSize = 5;
                h1_dither.XAxis.FontWeight = 'bold';
                h1_dither.YAxis.FontWeight = 'bold';
                colorbar;
                clim([0,1])
                axis image;
                h1_image_dither.Parent.XLim = [-40,40];
                h1_image_dither.Parent.YLim = [-40,40];

                h2_dither = subplot(3,3,6,'Parent',fig11);
                h2_dither.NextPlot = "add";
                h2_image_dither = plot(h2_dither,Z_exc,zPSF_exc_dither);
                h2_image_dither_unaberrated = plot(h2_dither,Z_exc,UnaberratedzPSF_exc_dither);
                title("Dithered Z-Excitation PSF - X=0" + "\lambda / n")
                h2_dither.Title.FontSize = 5;
                ylabel("Normalized a.u. ")
                xlabel("z/(\lambda_{exc}/n)")
                h2_dither.XAxis.FontSize = 4;
                h2_dither.YAxis.FontSize = 5;
                h2_dither.XAxis.FontWeight = 'bold';
                h2_dither.YAxis.FontWeight = 'bold';
                h2_image_dither.LineWidth = 2;
                h2_image_dither.Color = 'r';
                h2_image_dither.Parent.XLim = [-20,20];
                h2_image_dither.Parent.YAxis.TickValues = linspace(0,1,11);
                h2_image_dither.Parent.XAxis.TickValues = linspace(-20,20,21);
                h2_image_dither_unaberrated.LineWidth = 1;
                h2_image_dither_unaberrated.Color = 'g';
                lgd = legend(h2_dither,"Aberrated","Normal");
                lgd.FontSize = 3;
                grid on
                axis square
                h2_dither.NextPlot = "replace";

                h3_Unaberrated = subplot(3,3,7,'Parent',fig11);
                h3_image_unaberrated = imagesc(h3_Unaberrated,KX_exc,...
                          KZ_exc,...
                         real(UnaberratedxzOTF_exc)) ;
                title("real(XZ-Excitation OTF)")
                h3_Unaberrated.Title.FontSize = 5;
                xlabel("k_x/(4\pin/\lambda_{exc})")
                ylabel("k_z/(4\pin/\lambda_{exc})")
                h3_Unaberrated.XAxis.FontSize = 5;
                h3_Unaberrated.YAxis.FontSize = 5;
                h3_Unaberrated.XAxis.FontWeight = 'bold';
                h3_Unaberrated.YAxis.FontWeight = 'bold';
                colorbar;
                clim([0,1])
                axis image
                h3_image_unaberrated.Parent.XLim = [-0.5,0.5];
                h3_image_unaberrated.Parent.YLim = [-0.5,0.5];
        
                h3 = subplot(3,3,8,'Parent',fig11);
                h3_image = imagesc(h3,KX_exc,...
                          KZ_exc,...
                         abs(xzOTF_exc)) ;
                title("real(Aberrated XZ-Excitation OTF)")
                h3.Title.FontSize = 5;
                xlabel("k_x/(4\pin/\lambda_{exc})")
                ylabel("k_z/(4\pin/\lambda_{exc})")
                h3.XAxis.FontSize = 5;
                h3.YAxis.FontSize = 5;
                h3.XAxis.FontWeight = 'bold';
                h3.YAxis.FontWeight = 'bold';
                colorbar;
                clim([0,1])
                axis image
                h3_image.Parent.XLim = [-0.5,0.5];
                h3_image.Parent.YLim = [-0.5,0.5];
        
                h4 = subplot(3,3,9,'Parent',fig11);
                h4.NextPlot = "add";
                h4_image = plot(h4,KZ_exc,real(zOTF_exc));
                h4_image_abs = plot(h4,KZ_exc,abs(zOTF_exc));
                h4_image_unaberrated = plot(h4,KZ_exc,real(UnaberratedzOTF_exc));
                title("Z-Excitation OTF, K_X=0")
                h4.Title.FontSize = 5;
                xlabel("kz * \lambda / n")
                ylabel("Normalized a.u. ")
                h4.XAxis.FontSize = 5; 
                h4.YAxis.FontSize = 5; 
                h4.XAxis.FontWeight = 'bold'; 
                h4.YAxis.FontWeight = 'bold'; 
                h4_image.Color = 'r';
                h4_image.LineWidth = 2;
                h4_image.Parent.XLim = [-1,1];
                h4_image.Parent.YAxis.TickValues = linspace(-0.5,1,16);
                h4.XLim = [-0.5,0.5];
                h4.YLim = [-0.5,1];
                h4_image.Parent.XAxis.TickValues = linspace(-0.5,0.5,11);
                h4_image_unaberrated.Color = 'g';
                h4_image_unaberrated.LineWidth = 1;
                h4_image_abs.Color = 'b';
                h4_image_abs.LineWidth = 1.5;
                lgd = legend(h4,"real(Aberrated)","abs(Aberrated)","Normal");
                lgd.FontSize = 3;
                grid on
                axis square
                h4.NextPlot = "replace";
        
                exportgraphics(fig11,"Profile" + "N=" + num2str(i) + "M=" + num2str(AngularFrequency(k))...
                    + ".gif",'Append',true,'Resolution',300)

                close(fig11)
                clear fig11
            end
        end
    end

%     savefig(fig1, [pwd  '/ExcitationxzPSF.fig'])
%     savefig(fig2, [pwd  '/ExcitationzPSF.fig'])
%     savefig(fig3, [pwd  '/ExcitationxzOTF.fig'])
%     savefig(fig4, [pwd  '/ExcitationzOTF.fig'])
%     savefig(fig5, [pwd  '/ExcitationyzPSF.fig'])
%     savefig(fig6, [pwd  '/ExcitationyPSF.fig'])

    exportgraphics(fig1, [pwd  '/ExcitationxzPSF.png'],'Resolution',300)
    exportgraphics(fig2, [pwd  '/ExcitationzPSF.png'],'Resolution',300)
    exportgraphics(fig3, [pwd  '/ExcitationxzOTFreal.png'],'Resolution',300)
    exportgraphics(fig4, [pwd  '/ExcitationzOTF.png'],'Resolution',300)
    exportgraphics(fig5, [pwd  '/ExcitationyzPSF.png'],'Resolution',300)
    exportgraphics(fig6, [pwd  '/ExcitationyPSF.png'],'Resolution',300)
    exportgraphics(fig7, [pwd  '/OverallxzPSF.png'],'Resolution',300)
    exportgraphics(fig8, [pwd  '/OverallzPSF.png'],'Resolution',300)
    exportgraphics(fig9, [pwd  '/OverallxzOTFreal.png'],'Resolution',300)
    exportgraphics(fig10, [pwd  '/OverallzOTF.png'],'Resolution',300)
    exportgraphics(fig12, [pwd  '/ExcitationxzPSFNonDither.png'],'Resolution',300)
    exportgraphics(fig13, [pwd  '/Metrics.png'],'Resolution',300)

    h1 = subplot(1,1,1,'Parent',fig13);
    b1 = bar(1:length(CorrCoef),[CorrCoef;StrehlCenter;StrehlPeaks],0.5,'Parent',h1);
    lgd = legend(b1,"Corr Coef (realOTF)","Center Intensity Strehl Ratio","Peak Intensity Strehl Ratio");
    h1.XAxis.TickValues = 1:length(RadioOrderArray);
    LabelArray = [RadioOrderArray;AngularFrequencyArray];
    tickLabels = strtrim(sprintf('%d\\newline%d\n', LabelArray(:)));
    xticklabels(b1,tickLabels)
    grid on
