function SimulateLatticeAberration(LatticePupil,MaxRadialOrder,PhaseAmplitude)
    getParameters; %modify image parameter here
    CalculatePhysics;
    
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

    for i = 0:MaxRadialOrder
        RadialOrder = i*ones(1,i+1);
        AngularFrequency = -i:2:i;
        y = zernfun(RadialOrder,AngularFrequency,r(idx),theta(idx));

        AngularFrequency_iteration = 1:1:length(AngularFrequency);
        for k = 1:length(RadialOrder)
            AberratedPupil = zeros(size(phase));
            phase(idx) = y(:,k);
            AberratedPupil(idx) = LatticePupil(idx) .* exp(PhaseAmplitude.* 1i.*pi.*phase(idx));
            
            AberratedLatticePSF = zeros(N,N, N);
            AberratedLatticePSFDithered = zeros(N,N, N);
            % propagation
            for j = 1:length(y_exc)
                propagator_exc = exp(2*pi * 1i * ky_exc * y_exc(j));
                AberratedLatticePSF(:,:,j) = abs( fftshift( ifft2(AberratedPupil .* propagator_exc) ) ).^2;
                AberratedLatticePSFDithered(:,:,j) = meshgrid(mean(AberratedLatticePSF(:,:,j),2))';
            end  
            % Normalize
            AberratedLatticePSF = AberratedLatticePSF/max(max(max(AberratedLatticePSF)));
            AberratedLatticePSFDithered = AberratedLatticePSFDithered/max(max(max(AberratedLatticePSFDithered)));

            xzPSF = AberratedLatticePSF(:,:,(N+1)/2); xzPSF = xzPSF/max(max(xzPSF));
            zPSF = xzPSF(:,(N+1)/2); zPSF = zPSF/max(zPSF);
            xzOTF = abs(fftshift(fft2(xzPSF))); xzOTF = xzOTF/max(max(xzOTF));
            zOTF = xzOTF(:,(N+1)/2); zOTF = zOTF/max(zOTF);
            yzPSF = squeeze(AberratedLatticePSF(:,(N+1)/2,:)); yzPSF = yzPSF/max(max(yzPSF));
            yPSF = yzPSF((N+1)/2,:); yPSF = yPSF/max(yPSF);

            index = find(yPSF((N+1)/2:end) <= 0.5);
            if ~isempty(index)
                yFWHM = 2*Y_exc((N+1)/2+index(1)-1);
            else
                yFWHM = "N/A";
            end

            h1 = subplot(MaxRadialOrder+1,MaxRadialOrder+1,((MaxRadialOrder+1)*i)+AngularFrequency_iteration(k),'Parent',fig1);
                imagesc(h1,X_exc,Z_exc,xzPSF);
                h1.XAxis.Label.String = "x(\lambda_{exc}/n)";
                h1.YAxis.Label.String = "z(\lambda_{exc}/n)";
                h1.Colormap = colormap(hot(256));
                h1.XAxis.Limits = [-20,20];
                h1.YAxis.Limits = [-20,20];
                h1.Title.String = ['Z_{' num2str(RadialOrder(k)) '}^{' num2str(AngularFrequency(k)) '}'];
                h1.DataAspectRatio = [1,1,1];

            h2 = subplot(MaxRadialOrder+1,MaxRadialOrder+1,((MaxRadialOrder+1)*i)+AngularFrequency_iteration(k),'Parent',fig2);
                plot(h2,Z_exc,zPSF,'LineWidth',2,'Color','r')
                h2.XAxis.Label.String = "z(\lambda_{exc}/n)";
                h2.XAxis.Limits = [-20,20];
                h2.Title.String = ['Z_{' num2str(RadialOrder(k)) '}^{' num2str(AngularFrequency(k)) '}'];
                h2.XGrid = 'on';
                h2.YGrid = 'on';
                grid on

            h3 = subplot(MaxRadialOrder+1,MaxRadialOrder+1,((MaxRadialOrder+1)*i)+AngularFrequency_iteration(k),'Parent',fig3);
                imagesc(h3,KX_exc,KZ_exc,xzOTF)
                h3.XAxis.Label.String = "k_x/(4\pin/\lambda_{exc})";
                h3.YAxis.Label.String = "k_z/(4\pin/\lambda_{exc})";
                h3.Colormap = colormap(hot(256));
                h3.XAxis.Limits = [-0.5,0.5];
                h3.YAxis.Limits = [-0.5,0.5];
                h3.Title.String = ['Z_{' num2str(RadialOrder(k)) '}^{' num2str(AngularFrequency(k)) '}'];
                h3.DataAspectRatio = [1,1,1];

            h4 = subplot(MaxRadialOrder+1,MaxRadialOrder+1,((MaxRadialOrder+1)*i)+AngularFrequency_iteration(k),'Parent',fig4);
                plot(h4,KZ_exc,zOTF,'LineWidth',2,'Color','r')
                h4.XAxis.Label.String = "k_z/(4\pin/\lambda_{exc})";
                h4.XAxis.Limits = [-0.5,0.5];
                h4.Title.String = ['Z_{' num2str(RadialOrder(k)) '}^{' num2str(AngularFrequency(k)) '}'];
                grid on
                h4.XGrid = 'on';
                h4.YGrid = 'on';

            h5 = subplot(MaxRadialOrder+1,MaxRadialOrder+1,((MaxRadialOrder+1)*i)+AngularFrequency_iteration(k),'Parent',fig5);
                imagesc(h5,Y_exc((N+1)/2:end),Z_exc,yzPSF(:,(N+1)/2:end))
                h5.XAxis.Label.String = "y(\lambda_{exc}/n)";
                h5.YAxis.Label.String = "z(\lambda_{exc}/n)";
                h5.Colormap = colormap(hot(256));
                h5.YAxis.Limits = [-40,40];
                h5.Title.String = ['Z_{' num2str(RadialOrder(k)) '}^{' num2str(AngularFrequency(k)) '}'];
                h5.DataAspectRatio = [1,1,1];

            h6 = subplot(MaxRadialOrder+1,MaxRadialOrder+1,((MaxRadialOrder+1)*i)+AngularFrequency_iteration(k),'Parent',fig6);
                plot(h6,Y_exc((N+1)/2:end),yPSF((N+1)/2:end),'LineWidth',2,'Color','r')
                h6.XAxis.Label.String = "z(\lambda_{exc}/n)";
                h6.Title.String = ['Z_{' num2str(RadialOrder(k)) '}^{' num2str(AngularFrequency(k)) '}' 'yFWHM=' num2str(yFWHM) '/lambda'];
                h6.XGrid = 'on';
                h6.YGrid = 'on';
        end
    end

    savefig(fig1, [pwd  '/ExcitationxzPSF.fig'])
    savefig(fig2, [pwd  '/ExcitationzPSF.fig'])
    savefig(fig3, [pwd  '/ExcitationxzOTF.fig'])
    savefig(fig4, [pwd  '/ExcitationzOTF.fig'])
    savefig(fig5, [pwd  '/ExcitationyzPSF.fig'])
    savefig(fig6, [pwd  '/ExcitationyPSF.fig'])

    exportgraphics(fig1, [pwd  '/ExcitationxzPSF.png'],'Resolution',500)
    exportgraphics(fig2, [pwd  '/ExcitationzPSF.png'],'Resolution',500)
    exportgraphics(fig3, [pwd  '/ExcitationxzOTF.png'],'Resolution',500)
    exportgraphics(fig4, [pwd  '/ExcitationzOTF.png'],'Resolution',500)
    exportgraphics(fig5, [pwd  '/ExcitationyzPSF.png'],'Resolution',500)
    exportgraphics(fig6, [pwd  '/ExcitationyPSF.png'],'Resolution',500)
