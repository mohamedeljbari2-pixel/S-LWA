function [H, S_n] = generate_channels(Prms, p)
    K     = Prms.K;
    N_sb  = Prms.N_sb;
    f     = Prms.f;
    phi_u = Prms.phi_u;
    rho_u = Prms.rho_u;
    c     = 3e8;
    alpha = Prms.alpha_att;

    b_eff = p(1);
    L_eff = p(2);

    H = zeros(K, N_sb);
    for k = 1:K
        for n = 1:N_sb
            Gval = SLWA_gain(phi_u(k), f(n), b_eff, L_eff, alpha);
            % Stronger scaling so that beam pointing actually changes rate
            PL   = (c/f(n))/(4*pi*rho_u(k));
            H(k,n) = abs(Gval) * PL * 8e3;   % ← critical scaling
        end
    end
    [~, S_n] = max(H, [], 1);
end