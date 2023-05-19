function Phase = GetSingleZmodePupil(RadialOrder,AngularFrequency,PupilNA)
    getParameters; %modify image parameter here
    CalculatePhysics;
    
    [theta,r] = cart2pol(kx_exc./(PupilNA./n*k_wave),kz_exc./(PupilNA./n*k_wave));
    idx = r<=1;
    Zmode = zeros(size(kx_exc));
    Zmode(idx) = zernfun(RadialOrder,AngularFrequency,r(idx),theta(idx),'norm');
    Phase = (Zmode) .* 2*pi ./ wavelength_exc;

    