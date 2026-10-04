clear; clc; close all;
addpath('functions');

%% System parameters
Prms.K        = 4;
Prms.N_sb     = 32;
Prms.B_w      = 1e9;
Prms.f        = linspace(90e9,100e9,Prms.N_sb);
Prms.phi_u    = [-60,-20,20,60]*pi/180;
Prms.rho_u    = [10 15 20 25];
Prms.phi_t    = 0;
Prms.rho_t    = 20;
Prms.kappa_t  = 1;
Prms.b_base   = 1.5e-3;
Prms.Delta_b  = 0.15e-3;
Prms.L_eff0   = 0.10;
Prms.alpha_att= 5;
Prms.sigma2   = 10^(-174/10)*1e-3*(Prms.B_w/Prms.N_sb);
Prms.sigma_s2 = 10^(-170/10)*1e-3*(Prms.B_w/Prms.N_sb);
Prms.L_tb     = 500;
Prms.Nmax     = 40;

scenarios = {
    struct('Pmax_dBm',26,'lambda',0.0,'label','P=26 dBm, \lambda=0 (comm-only)');
    struct('Pmax_dBm',23,'lambda',0.2,'label','P=23 dBm, \lambda=0.2 (ISAC)');
    struct('Pmax_dBm',20,'lambda',0.5,'label','P=20 dBm, \lambda=0.5 (sensing-heavy)');
};

rate_hist = zeros(Prms.Nmax,3);
obj_hist  = zeros(Prms.Nmax,3);
crb_hist  = zeros(Prms.Nmax,3);

for s = 1:3
    Prms.Pmax = 10^((scenarios{s}.Pmax_dBm-30)/10);
    [~,~,rate_hist(:,s),crb_hist(:,s),obj_hist(:,s)] = ...
        joint_AO(Prms, scenarios{s}.lambda, Prms.Nmax);
end

%% Plot
figure('Position',[100 100 1100 420]);
cols = [0.85 0.1 0.1; 0.1 0.55 0.1; 0.1 0.1 0.75];

subplot(1,2,1);
for s=1:3, plot(1:Prms.Nmax,rate_hist(:,s),'-','Color',cols(s,:),'LineWidth',1.9); hold on; end
grid on; xlabel('Number of iterations'); ylabel('Sum-rate (Gbps)');
legend({scenarios{1}.label,scenarios{2}.label,scenarios{3}.label},'Location','southeast');
title('(a) Sum-rate convergence');

subplot(1,2,2);
for s=1:3, plot(1:Prms.Nmax,obj_hist(:,s),'-','Color',cols(s,:),'LineWidth',1.9); hold on; end
grid on; xlabel('Number of iterations'); ylabel('Joint objective');
legend({scenarios{1}.label,scenarios{2}.label,scenarios{3}.label},'Location','southeast');
title('(b) Joint objective convergence');

sgtitle('Convergence of Joint SLWA-OFDMA-ISAC Algorithm','FontWeight','bold');