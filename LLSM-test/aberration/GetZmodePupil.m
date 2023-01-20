function GetZmodePupil(MaxRadialOrder)
    getParameters; %modify image parameter here
    CalculatePhysics;
    
    [theta,r] = cart2pol(kx_exc./(0.6./n*k_wave),kz_exc./(0.6./n*k_wave));
    idx = r<=1;
    phase = zeros(size(kx_exc));

    fig1 = figure('Units','normalized');
    MinRadialOrder = 0;
    for i = MinRadialOrder:MaxRadialOrder
        RadialOrder = i*ones(1,i+1);
        AngularFrequency = -i:2:i;
        AngularFrequency_iteration = 1:1:length(AngularFrequency);

        for k = 1:length(AngularFrequency)
            phase(idx) = zernfun(i,AngularFrequency(k),r(idx),theta(idx),'norm');
            subplot(MaxRadialOrder-MinRadialOrder+1,MaxRadialOrder+1,((MaxRadialOrder+1)*(i-MinRadialOrder))+AngularFrequency_iteration(k),'Parent',fig1)
            pcolor(KX_exc,KZ_exc,phase), shading interp
            set(gca,'XTick',[],'YTick',[])
            axis square
            colorbar
            title(['Z_{' num2str(RadialOrder(k)) '}^{' num2str(AngularFrequency(k)) '}'])
        end
    end