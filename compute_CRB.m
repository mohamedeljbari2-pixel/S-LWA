function [CRB_phi, CRB_r, CRB_v] = compute_CRB(P_n, Prms, p)
    f        = Prms.f;
    phi_t    = Prms.phi_t;
    rho_t    = Prms.rho_t;
    kappa    = Prms.kappa_t;
    sigma_s2 = Prms.sigma_s2;
    alpha    = Prms.alpha_att;
    c        = 3e8;
    L_tb     = Prms.L_tb;

    b_eff = p(1);
    L_eff = p(2);

    g = zeros(length(f),1);
    for n = 1:length(f)
        Gval = SLWA_gain(phi_t, f(n), b_eff, L_eff, alpha);
        g(n) = abs(Gval)^2 * sqrt(max(P_n(n),0)) * (c/(4*pi*f(n)*rho_t^2));
    end

    divers  = 1.0 + 0.55*(std(f)/mean(f));
    SNR_eff = kappa^2 * sum(abs(g).^2) / (sigma_s2 + eps);

    CRB_phi = (1/(2*SNR_eff*L_tb + eps)) / divers * (180/pi)^2;
    CRB_r   = CRB_phi * (rho_t/12)^2;
    CRB_v   = 0.08 / (SNR_eff + eps);
end