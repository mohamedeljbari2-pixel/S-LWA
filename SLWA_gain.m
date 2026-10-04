function G = SLWA_gain(phi, f, b_eff, L_eff, alpha_att)
% Physically scaled SLWA gain (Eq. 5)
    c     = 3e8;
    beta0 = 2*pi*f/c;
    arg   = 1 - (c./(2*b_eff.*f)).^2;
    beta  = beta0 .* sqrt(max(arg, 0));
    z     = (beta - 1j*alpha_att - beta0.*cos(phi)) * (L_eff/2);
    z(abs(z)<1e-12) = 1e-12;
    G     = L_eff * (sin(z)./z);          % complex pattern
end