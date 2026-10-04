function [p_opt, P_n_opt, rate_hist, crb_hist, obj_hist] = joint_AO(Prms, lambda, Nmax)

    p         = [Prms.b_base; Prms.L_eff0];
    Pmax      = Prms.Pmax;
    rate_hist = zeros(Nmax,1);
    crb_hist  = zeros(Nmax,1);
    obj_hist  = zeros(Nmax,1);

    % Larger initial step – this is the key to visible convergence
    step_b = 0.25 * Prms.Delta_b;

    for it = 1:Nmax
        % ========== Step 1: channels + power allocation ==========
        [H, S_n] = generate_channels(Prms, p);
        P_n      = waterfilling_OFDMA(H, S_n, Pmax, Prms.sigma2);

        % Current rate
        Rsum = 0;
        for n = 1:Prms.N_sb
            k    = S_n(n);
            SINR = H(k,n)^2 * P_n(n) / (Prms.sigma2 + 0.1*sum(H(:,n).^2)*P_n(n) + eps);
            Rsum = Rsum + (Prms.B_w/Prms.N_sb)*log2(1 + SINR);
        end

        % ========== Step 2: projected gradient on b_eff ==========
        eps_b = 5e-6;
        p_try = p;
        p_try(1) = p(1) + eps_b;
        [Htry, Stry] = generate_channels(Prms, p_try);

        Rtry = 0;
        for n = 1:Prms.N_sb
            k    = Stry(n);
            SINR = Htry(k,n)^2 * P_n(n) / (Prms.sigma2 + 0.1*sum(Htry(:,n).^2)*P_n(n) + eps);
            Rtry = Rtry + (Prms.B_w/Prms.N_sb)*log2(1 + SINR);
        end

        grad_b = (Rtry - Rsum) / eps_b;

        % Adaptive & sufficiently large step
        if abs(grad_b) < 1
            step = step_b * 0.3;
        else
            step = step_b;
        end

        p(1) = p(1) + step * sign(grad_b + 1e-12);
        p(1) = max(Prms.b_base - Prms.Delta_b, ...
                   min(Prms.b_base + Prms.Delta_b, p(1)));

        % ========== Metrics ==========
        [CRB_phi, ~, ~] = compute_CRB(P_n, Prms, p);

        rate_hist(it) = Rsum / 1e9;          % Gbps
        crb_hist(it)  = CRB_phi;
        obj_hist(it)  = (1-lambda)*rate_hist(it) - lambda * min(CRB_phi, 50);

        % Mild step-size reduction after some iterations
        if it > 12
            step_b = max(step_b * 0.96, 0.04*Prms.Delta_b);
        end
    end

    p_opt   = p;
    P_n_opt = P_n;
end