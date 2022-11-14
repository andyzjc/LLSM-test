function GetZmodePupil(MaxRadialOrder)
    getParameters; %modify image parameter here
    CalculatePhysics;
    
    [theta,r] = cart2pol(kx_exc./(1.2./n*k_wave),kz_exc./(1.2./n*k_wave));
    idx = r<=1;
    z = nan(size(kx_exc));

    fig1 = figure('Units','normalized');
    for i = 0:MaxRadialOrder
        RadialOrder = i*ones(1,i+1);
        AngularFrequency = -i:2:i;
        y = zernfun(RadialOrder,AngularFrequency,r(idx),theta(idx));

        AngularFrequency_iteration = 1:1:length(AngularFrequency);
        for k = 1:length(RadialOrder)
            z(idx) = y(:,k);
            subplot(MaxRadialOrder+1,MaxRadialOrder+1,((MaxRadialOrder+1)*i)+AngularFrequency_iteration(k),'Parent',fig1)
            pcolor(KX_exc,KZ_exc,z), shading interp
            set(gca,'XTick',[],'YTick',[])
            axis square
            title(['Z_{' num2str(RadialOrder(k)) '}^{' num2str(AngularFrequency(k)) '}'])
        end
    end