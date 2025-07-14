function [Ilocs, gammaVals] = tCWO(A, Y, K, sigc2, nu, maxCycles)
% =========================================================================
% tCWO: Robust Coordinate-wise Optimization algorithm for sparse signal recovery
%
% Inputs:
%   A        - Dictionary matrix of size L x N (L: measurement vectors, N: dictionary entries)
%   Y        - Measurement matrix of size L x M (M: signal length)
%   K        - Number of non-zero sources
%   sigc2    - Noise variance
%   nu       - degree of freedom
%   maxCycles- Maximum number of iterations for coordinate-wise optimization
%
% Outputs:
%   Ilocs    - Support indices of non-zero signal powers (K-vector)
%   gammaVals- Power estimates for each signal (N-vector)
% 
% 2024.10.20 WXJ
% 2024.12.2
% =========================================================================

% Initialization
L = size(A, 1);        % Number of devices
N = size(A, 2);        % Number of dictionary entries
M = size(Y, 2);        % Number of signal length
Sigmainv = (1/sigc2) * eye(L);  % Initial inverse covariance matrix
gamma = zeros(N, 1);   % Initial power estimates

% Tolerance for convergence
deltaCWO = 5e-3;

u_func = @(t) (nu + 2*L) ./ (nu + 2*t); % Weight function u(t)

integrand = @(t) (t/2) .* (nu + 2*L) ./ (nu + t) .* chi2pdf(t, 2*L);
b = integral(integrand, 0, Inf)/L;
const = (1/(b*M));


%% Main Iteration Phase
for cycle = 1:maxCycles

    gammaNew = zeros(N, 1);  % Temporary storage for updated gamma values

    % Cyclically update each gamma_i
    for i = 1:N

        ai = A(:, i);
        ci = Sigmainv * ai;
        ciH = ci';
        ciciH = ci * ciH;
        den_aiHci = real(ciH * ai);  % ai^H * ci (scalar)

        % Update inverse covariance matrix without i-th power contribution
        Theta = Sigmainv + (gamma(i) / (1 - gamma(i) * den_aiHci)) * ciciH;


        % Precompute values used in each iteration
        bi = Theta * ai;  % L x 1 vector
        MD = real(sum(conj(Y) .* (Theta * Y)));  % y_m^H \Sigma_{\i}^-1 y_m (1 x M)
        v = abs(sum(Y .* conj(bi))).^2;  % |y_m^H b_i|^2 (1 x M)
        di = real(ai' * bi);  % Scalar value for ai^H * bi
        % Estimate the new gamma_i using the fixed-point method
        gammaNew(i) = fixed_point_iteration(MD,v,di,u_func,gamma(i),const,10,5e-3);

        % Update the covariance matrix with the new gamma_i
        delta_i = gammaNew(i) - gamma(i);
        if delta_i ~= 0
            Sigmainv = Sigmainv - (delta_i / (1 + delta_i * den_aiHci)) * ciciH;
        end
    end

    % Check convergence based on the relative error
    if norm(gammaNew - gamma, Inf) / norm(gammaNew, Inf) < deltaCWO
        break;  % Convergence achieved
    end

    gamma = gammaNew;  % Update gamma for the next cycle
end

% Select the K-largest gamma values as the final support
[~, sortedIndices] = sort(gammaNew, 'descend');
Ilocs = sortedIndices(1:K);
gammaVals = gammaNew;

end

% ======== Fixed-point iteration for gamma update ========
function gamma_t = fixed_point_iteration(MD,v,di,u_func,gamma_t, const,maxT, deltaT)

    % Fixed-point iteration
    for t = 1:maxT
        % Step 1: Compute Mahalanobis distance for m = 1,...,M
        s = MD - gamma_t / (1 + gamma_t * di) * v;  % 1 x M

        % Step 2: Compute weighted SCM with bi
        weighted_s = u_func(s);  % 1 x M
        sigmasq_est = const*sum(weighted_s .* v);  % Scalar estimation

        % Step 3: Compute new gamma
        gamma_t_new_simple = max((sigmasq_est - di) / (di^2), 0);

        % Check convergence
        if gamma_t == gamma_t_new_simple || abs(gamma_t_new_simple - gamma_t) < deltaT
            gamma_t = gamma_t_new_simple;
            break;
        end
        gamma_t = gamma_t_new_simple;
    end
end
