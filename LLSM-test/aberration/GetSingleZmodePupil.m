function Phase_factor = GetSingleZmodePupil(RadialOrder,AngularFrequency,PhaseAmplitude)
    getParameters; %modify image parameter here
    CalculatePhysics;
    
    [theta,r] = cart2pol(kx_exc./(0.6./n*k_wave),kz_exc./(0.6./n*k_wave));
    idx = r<=1;
    Phase = zeros(size(kx_exc));
    Phase(idx) = zernfun(RadialOrder,AngularFrequency,r(idx),theta(idx),'norm');
    Phase(idx) = PhaseAmplitude .* Phase(idx);
    Phase_factor = exp(1i .* Phase );