% Cellular MIMO Uplink Simulation with ASPG Algorithm
% This script evaluates multiple activity detection algorithms, including ASPG,
% in a massive MIMO uplink scenario under impulsive noise conditions.

clc; 
clear all; 
close all;
s = RandStream('mt19937ar', 'Seed', 2);
RandStream.setGlobalStream(s);

% Add paths (adjust as needed)
addpath '../CL_toolbox/SSR_algorithms/'
addpath '../algorithms'
addpath '../CL_toolbox'
addpath '../'

%% Parameters
N = 1000;
M_list = [20];
L_list = [20];
K_a_list = [20];
SNR_list = -15:4:0;      
MC_iters = 80;

% Storage for results (9 algorithms, ASPG now after CWOpt)
Pmd = zeros(length(M_list), length(L_list), length(K_a_list), length(SNR_list), 9);
Pfa = zeros(length(M_list), length(L_list), length(K_a_list), length(SNR_list), 9);
tme = zeros(length(M_list), length(L_list), length(K_a_list), length(SNR_list), 9);

%% Monte Carlo Loop
t_allstart = tic;
for snr_idx = 1:length(SNR_list)
    SNR = SNR_list(snr_idx);
    fprintf(' %%%%%%%% SNR = %.1f dB in %.2f mins %%%%%%%% \n', SNR, toc(t_allstart) / 60);
    for m_idx = 1:length(M_list)
        M = M_list(m_idx);
        for l_idx = 1:length(L_list)
            L = L_list(l_idx);
            for ka_idx = 1:length(K_a_list)
                K_a = K_a_list(ka_idx);
                % Monte Carlo iterations
                Pmd_mc = zeros(MC_iters, 9);
                Pfa_mc = zeros(MC_iters, 9);
                tme_mc = zeros(MC_iters, 9);
                % parfor mc_iter = 1:MC_iters
                for mc_iter = 1:MC_iters
                    temp_Pmd_mc = zeros(1, 9);
                    temp_Pfa_mc = zeros(1, 9);
                    temp_tme_mc = zeros(1, 9);

                    % --- Generate path loss coefficients ---
                    distance = 50 + (1000 - 50) * rand(N, 1); % Distances in [0.05, 1] km
                    beta_lsfc = 10.^((-128.1 - 37.6 .* log10(distance .* 1e-3)) / 10);
                    nBeta = beta_lsfc;

                    % --- Generate transmit signals ---
                    varpho_max = 0.2; % Max transmission power (W)
                    sup = randperm(N, K_a);
                    [~, tmp] = sort(beta_lsfc(sup), 'descend');
                    sup = sup(tmp);
                    alphaN = zeros(N, 1); alphaN(sup) = 1;
                    gamma = zeros(N, 1);
                    gamma(sup) = varpho_max * min(beta_lsfc);

                    small_scale_fading = sqrt(1/2) * (randn(K_a, M) + 1i * randn(K_a, M));
                    X = zeros(N, M);
                    X(sup, :) = diag(sqrt(gamma(sup))) * small_scale_fading;

                    % --- Generate pilot matrix ---
                    A = (1-2*binornd(1,0.5,[L, N]))/sqrt(2) + 1i*(1-2*binornd(1,0.5,[L, N]))/sqrt(2);

                    % --- Noise power ---
                    gamma_sig = min(gamma(sup));
                    sigma2_noise = gamma_sig / (10^(SNR / 10));
                    nsPower = sigma2_noise;

                    % --- Generate noise matrix ---
                    Z = generateNoise(L, M, 't', struct('nu', 2.5), nsPower);

                    % --- Received signal model ---
                    Y = A(:, sup) * X(sup, :) + Z;
                    cov_m = (1/M) * Y * Y';
                    Y_scaled = Y ./ sqrt(min(beta_lsfc));
                    nsPower_scaled = nsPower / min(beta_lsfc);
                    cov_m_scaled = cov_m / min(beta_lsfc);

                    % --- Algorithms ---
                    % 1. Noisy CAMP MMSE (Commented out for speed)
                    % tStart = tic;
                    % [temp_Pmd_mc(1), temp_Pfa_mc(1), ~] = noisyCAMPmmseforKLS(A, N, M, L, Y_scaled, X, alphaN, 150, K_a/N, nBeta, sqrt(nsPower));
                    % temp_tme_mc(1) = toc(tStart);

                    % 2. SBL Basic
                    tStart = tic;
                    [temp_Pmd_mc(2), temp_Pfa_mc(2), ~] = sblBasicFunc(Y_scaled, A, X, alphaN, nsPower_scaled, 150, 1e-4);
                    temp_tme_mc(2) = toc(tStart);

                    % 3. CLMP
                    tStart = tic;
                    sup_clmp = CLMP(A, cov_m_scaled, K_a, nsPower_scaled);
                    [temp_Pmd_mc(3), temp_Pfa_mc(3)] = computePmdPfa(sup_clmp, sup, N);
                    temp_tme_mc(3) = toc(tStart);

                    % 4. CWOpt
                    tStart = tic;
                    [sup_cwopt, ~, ~] = CWOpt(A, Y_scaled, K_a, nsPower_scaled);
                    [temp_Pmd_mc(4), temp_Pfa_mc(4)] = computePmdPfa(sup_cwopt, sup, N);
                    temp_tme_mc(4) = toc(tStart);

                    % 5. ASPG (Active Set Projected Gradient) - NEW POSITION
                    tStart = tic;
                    [sup_aspg, ~] = ASCL2(A, cov_m_scaled, K_a, M, nsPower_scaled);
                    [temp_Pmd_mc(5), temp_Pfa_mc(5)] = computePmdPfa(sup_aspg, sup, N);
                    temp_tme_mc(5) = toc(tStart);

                    % 6. SNIHT (moved down)
                    tStart = tic;
                    [~, sup_sniht] = SNIHT(Y_scaled, A, K_a);
                    [temp_Pmd_mc(6), temp_Pfa_mc(6)] = computePmdPfa(sup_sniht, sup, N);
                    temp_tme_mc(6) = toc(tStart);

                    % 7. SOMP (moved down)
                    tStart = tic;
                    [~, sup_somp] = SOMP(Y_scaled, A, K_a);
                    [temp_Pmd_mc(7), temp_Pfa_mc(7)] = computePmdPfa(sup_somp, sup, N);
                    temp_tme_mc(7) = toc(tStart);

                    % 8. RCWO (moved down)
                    tStart = tic;
                    [sup_rcwot, ~] = tCWO(A, Y_scaled, K_a, nsPower_scaled, 2.5, 50);
                    [temp_Pmd_mc(8), temp_Pfa_mc(8)] = computePmdPfa(sup_rcwot, sup, N);
                    temp_tme_mc(8) = toc(tStart);

                    % 9. RCL-MP (moved down)
                    tStart = tic;
                    [sup_rclmpt, ~] = tMP(A, Y_scaled, K_a, nsPower_scaled, 2.5);
                    [temp_Pmd_mc(9), temp_Pfa_mc(9)] = computePmdPfa(sup_rclmpt, sup, N);
                    temp_tme_mc(9) = toc(tStart);

                    % Store results
                    Pmd_mc(mc_iter, :) = temp_Pmd_mc;
                    Pfa_mc(mc_iter, :) = temp_Pfa_mc;
                    tme_mc(mc_iter, :) = temp_tme_mc;
                end
                % Average results
                Pmd(m_idx, l_idx, ka_idx, snr_idx, :) = mean(Pmd_mc);
                Pfa(m_idx, l_idx, ka_idx, snr_idx, :) = mean(Pfa_mc);
                tme(m_idx, l_idx, ka_idx, snr_idx, :) = mean(tme_mc);
            end
        end
    end
end

elapsed_time = toc(t_allstart) / 3600;
fprintf(' === All done in %.2f mins, or %.3f hours. === \n', elapsed_time*60, elapsed_time);

%% Save Results
% filename = sprintf('mat/SNRvsPMDPFA_SNR=%dto%d_K=%d_L=%d_N=%d_MC=%d_%s.mat', ...
    % SNR_list(1), SNR_list(end), K_a_list(1), L_list(1), N, MC_iters, datestr(now, 'mmdd-HHMM-SS'));
% save(filename, 'Pmd', 'Pfa', 'tme');

%% Plotting
% Define colors and markers with ASPG having similar style to CWOpt
colors = {
    [1, 0, 0],           % Red - VAMP
    [0, 0, 1],           % Blue - SBL
    [0, 0.7, 0],         % Green - CLMP
    [0, 0.7, 0.7],       % Cyan - CWOpt
    [0.2, 0.5, 0.8],     % Light Blue - ASPG (similar to CWOpt but distinguishable)
    [1, 0, 1],           % Magenta - SNIHT
    [1, 1, 0],           % Yellow - SOMP
    [0, 0, 0],           % Black - RCWO
    [0.5, 0.5, 0.5]      % Gray - RCL-MP
};

markers = {'-o', '-x', '-s', '-d', '-^', '-v', '-p', '-h', '-*'};
algorithm_names = {'VAMP', 'SBL', 'CLMP', 'CWOpt', 'ASPG', 'SNIHT', 'SOMP', 'RCWO', 'RCL-MP'};

fixed_L_idx = 1;
fixed_M_idx = 1;
fixed_Ka_idx = 1;

% Plot Pmd
figure(1); clf;
hold on;
for alg_idx = 2:9
    plot(SNR_list, squeeze(Pmd(fixed_M_idx, fixed_L_idx, fixed_Ka_idx, :, alg_idx)), ...
        markers{alg_idx}, 'Color', colors{alg_idx}, 'LineWidth', 2, ...
        'MarkerSize', 8, 'DisplayName', algorithm_names{alg_idx});
end
set(gca, 'YScale', 'log');
xlabel('SNR / dB', 'FontSize', 14);
ylabel('P_{md}', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 12);
grid on;
title(sprintf('P_{md} vs SNR, N = %d, M = %d, L = %d, K_a = %d, MC = %d', ...
    N, M_list(fixed_M_idx), L_list(fixed_L_idx), K_a_list(fixed_Ka_idx), MC_iters), ...
    'FontSize', 14);
hold off;

% Plot Runtime
figure(2); clf;
hold on;
for alg_idx = 2:9
    plot(SNR_list, squeeze(tme(fixed_M_idx, fixed_L_idx, fixed_Ka_idx, :, alg_idx)), ...
        markers{alg_idx}, 'Color', colors{alg_idx}, 'LineWidth', 2, ...
        'MarkerSize', 8, 'DisplayName', algorithm_names{alg_idx});
end
set(gca, 'YScale', 'log');
xlabel('SNR / dB', 'FontSize', 14);
ylabel('Time (s)', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 12);
grid on;
title(sprintf('Runtime vs SNR, N = %d, M = %d, L = %d, K_a = %d, MC = %d', ...
    N, M_list(fixed_M_idx), L_list(fixed_L_idx), K_a_list(fixed_Ka_idx), MC_iters), ...
    'FontSize', 14);
hold off;
