function Phase = GetSingleZmodePupil(RadialOrder,AngularFrequency,PhaseAmplitude)
    getParameters; %modify image parameter here
    CalculatePhysics;
    
    [theta,r] = cart2pol(kx_exc./(1.2./n*k_wave),kz_exc./(1.2./n*k_wave));
    idx = r<=1;
    Phase = zeros(size(kx_exc));
    Phase(idx) = zernfun(RadialOrder,AngularFrequency,r(idx),theta(idx));
    Phase = exp(PhaseAmplitude .* 1i .* pi .* Phase);