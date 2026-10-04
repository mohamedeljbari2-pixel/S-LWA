function P_n = waterfilling_OFDMA(H, S_n, Pmax, sigma2)
    [~, N_sb] = size(H);
    alpha = zeros(N_sb,1);
    beta  = sigma2 * ones(N_sb,1);

    for n = 1:N_sb
        k = S_n(n);
        alpha(n) = H(k,n)^2;
        interf   = sum(H(:,n).^2) - H(k,n)^2;
        beta(n)  = beta(n) + 0.25*interf;
    end

    % Robust bisection
    nu_lo = 1e-12; nu_hi = 1e6;
    for it = 1:50
        nu = sqrt(nu_lo*nu_hi);
        P  = max(0, 1/nu - beta./(alpha+eps));
        if sum(P) > Pmax
            nu_lo = nu;
        else
            nu_hi = nu;
        end
    end
    P_n = P * (Pmax / max(sum(P), eps));
end