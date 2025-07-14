%% Activity Detection Performance Analysis

clc
clear all
close all

% Add paths
addpath 'CL_toolbox/SSR_algorithms/'
addpath 'CL_toolbox/'
addpath 'algorithms/'
addpath 'utils/'

%% General Settings
MC_iters = 20;   % number of MC iterations
P = 1;           % total scale (power) of X
N = 1000;        % number of MTD-s

% Fading configuration
fading.type = 'uniform';
fading.lower_limit = -15;
fading.upper_limit = 0;

pilot = 'bernoulli';

% Parameter lists
Mlist = 30:30:90;
Klist = [5, 10, 20, 40];
Llist = 30:30:90;

%% Analysis 1: Performance vs M (number of antennas)
L = 30;  % fixed number of pilots
Pmd = zeros(length(Mlist), 4, length(Klist));
Per = zeros(length(Mlist), 4, length(Klist));

fprintf('Starting analysis: Performance vs M (antennas)\n');
t_allstart = tic;

for k = 1:length(Klist)
    rng('default');
    K_a = Klist(k);
    fprintf('\n---- K=%3d ----\n', K_a);
    
    for m = 1:length(Mlist)
        M = Mlist(m);
        fprintf('- M=%3d ', M);
        
        [Pmd(m,:,k), Per(m,:,k), ~, ~] = activityDetectionPE_Robust(...
            L, N, K_a, M, P, MC_iters, 'fading', fading, 'pilot', pilot);
        
        fprintf('in %.2f mins\n', toc(t_allstart) / 60);
    end
end

elapsed_time = toc(t_allstart) / 3600;
fprintf('=== M analysis done in %.2f mins (%.3f hours) ===\n', elapsed_time*60, elapsed_time);

%% Plot results for M analysis
markers = {'-d', '-o', '-x', '->', '-^', '--v', ':d'};
algorithm_names = {'CWO', 'CL-MP', 'RCWO', 'RCL-MP'};

figure(1);
clf;
set(gcf, 'Units', 'inches', 'Position', [2, 2, 20, 4.5]);

for k = 1:length(Klist)
    K_a = Klist(k);
    subplot(1, length(Klist), k);
    hold on;
    
    for alg_idx = 1:4
        plot(Mlist, Pmd(:, alg_idx, k), markers{alg_idx}, ...
            'DisplayName', algorithm_names{alg_idx}, 'LineWidth', 3);
    end
    
    set(gca, 'YScale', 'log');
    grid on;
    xlim([min(Mlist), max(Mlist)]);
    xlabel('$M$', 'FontName', 'Times New Roman', 'FontSize', 16, 'Interpreter', 'latex');
    
    if k == 1
        ylabel('PMD', 'FontName', 'Times New Roman', 'FontWeight', 'bold', ...
            'FontSize', 16, 'Interpreter', 'latex');
        legend(algorithm_names, 'Location', 'southwest', 'FontSize', 16);
    end
    
    title(sprintf('$K$=%d', K_a), 'FontName', 'Times New Roman', ...
        'FontSize', 24, 'Interpreter', 'latex');
    
    hold off;
end

%% Analysis 2: Performance vs L (number of pilots)
M = 30;  % fixed number of antennas
Pmd = zeros(length(Llist), 4, length(Klist));
Per = zeros(length(Llist), 4, length(Klist));

fprintf('\n============ Starting L analysis ============\n');
t_allstart = tic;

for k = 1:length(Klist)
    rng('default');
    K_a = Klist(k);
    fprintf('\n---- K=%3d ----\n', K_a);
    
    for l = 1:length(Llist)
        L = Llist(l);
        fprintf('- L=%3d ', L);
        
        [Pmd(l,:,k), Per(l,:,k), ~, ~] = activityDetectionPE_Robust(...
            L, N, K_a, M, P, MC_iters, 'fading', fading, 'pilot', pilot);
        
        fprintf('in %.2f mins\n', toc(t_allstart) / 60);
    end
end

elapsed_time = toc(t_allstart) / 3600;
fprintf('=== L analysis done in %.2f mins (%.3f hours) ===\n', elapsed_time*60, elapsed_time);

%% Plot results for L analysis
figure(2);
clf;
set(gcf, 'Units', 'inches', 'Position', [2, 2, 20, 4.5]);

for k = 1:length(Klist)
    K_a = Klist(k);
    subplot(1, length(Klist), k);
    hold on;
    
    for alg_idx = 1:4
        plot(Llist, Pmd(:, alg_idx, k), markers{alg_idx}, ...
            'DisplayName', algorithm_names{alg_idx}, 'LineWidth', 3);
    end
    
    set(gca, 'YScale', 'log');
    grid on;
    xlim([min(Llist), max(Llist)]);
    xlabel('$L$', 'FontName', 'Times New Roman', 'FontSize', 16, 'Interpreter', 'latex');
    
    if k == 1
        ylabel('PMD', 'FontName', 'Times New Roman', 'FontWeight', 'bold', ...
            'FontSize', 16, 'Interpreter', 'latex');
        legend(algorithm_names, 'Location', 'southwest', 'FontSize', 16);
    end
    
    title(sprintf('$K$=%d', K_a), 'FontName', 'Times New Roman', ...
        'FontSize', 24, 'Interpreter', 'latex');
    
    hold off;
end

fprintf('\n=== All analyses completed ===\n');