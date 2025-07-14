function [sup_idx, gam_vals] = ASCL2(A, Sigma_hat, K, M, nsPower, opt)
% ActiveSet: Active-set algorithm following ICASSP paper Algorithm 1
% Corrected version that properly uses omega_k and nu_k for active set selection

% ---------- Algorithm 1, Line 1: Initialize -----------------------------
[L, N] = size(A);
gamma = zeros(N, 1);  % γ^0 = 0
k = 0;                % k = 0

% Parameters initialization
if nargin < 6 || isempty(opt)
    opt = struct();
end
if ~isfield(opt, 'epsilon'), opt.epsilon = 1e-3; end

% Initialize covariance inverse
cov_inv = eye(L) / nsPower;

% Store fixed active set after warm-up
Aset_fixed = [];
grad_computed = false;

% ---------- Algorithm 1, Line 2: repeat ---------------------------------
while k < 20  % Maximum iterations
    k = k + 1;
    
    % ------ Set parameters {ω_k, ν_k, ε_k} (Paper Section 4) -----------
    omega_k = 10^(-6 - (k-1));  % Paper formula
    epsilon_k = max(10^(-(k-1)), 0.8 * 1e-3);
    
    % ------ Select active set A^k (Algorithm 1, Line 3) ----------------
    if k <= 2  % CHANGED: from 5 to 3 to 1 to 2
        % Warm-up phase: use all variables
        Aset = 1:N;
    else
        % After warm-up: use paper's selection rule
        
        % Compute gradient if needed
        if ~grad_computed || k == 3  % CHANGED: from k == 6 to k == 4 to k == 2 to k == 3
            S_inv_A = cov_inv * A;
            grad = real(sum(conj(S_inv_A) .* A, 1).' ...
                   - sum(conj(S_inv_A) .* (Sigma_hat * S_inv_A), 1).');
            grad_computed = true;
            
            % Compute nu_k based on gradient
            min_grad = min(grad);
            nu_k = min(10^(4 - (k-1)), 0.5 * abs(min_grad));
        else
            % Use stored gradient with updated nu_k
            nu_k = min(10^(4 - (k-1)), 0.5 * abs(min(grad)));
        end
        
        % Apply paper's selection rule (Equation 6)
        Aset = find(gamma > omega_k | grad < -nu_k);
        
        % Handle empty active set
        if isempty(Aset)
            % If empty, select based on gradient alone
            [~, idx] = sort(grad, 'ascend');
            Aset = idx(1:min(4*K, N));
        end
        
        % CRITICAL: Limit active set size to 3K
        if length(Aset) > 3*K
            % Score each variable in Aset
            scores = zeros(length(Aset), 1);
            for i = 1:length(Aset)
                idx = Aset(i);
                % Higher score for: large gamma, very negative gradient
                scores(i) = gamma(idx) / (max(gamma) + 1e-10) - grad(idx) / (max(abs(grad)) + 1e-10);
            end
            
            % Keep top 3K
            [~, keep_idx] = sort(scores, 'descend');
            Aset = Aset(keep_idx(1:3*K));
        end
        
        % Option: After first selection, fix the active set 
        if k == 3 && isempty(Aset_fixed)  
            Aset_fixed = Aset;
        elseif k > 3 && ~isempty(Aset_fixed) 
            Aset = Aset_fixed;
        end
    end
    
    % Debug output (optional)
    if k <= 10 && nargin > 5 && isfield(opt, 'verbose') && opt.verbose
        fprintf('Iter %d: omega_k=%.2e, nu_k=%.2e, |Aset|=%d\n', ...
            k, omega_k, exist('nu_k','var')*nu_k, length(Aset));
    end
    
    % ------ (Algorithm 1, Line 4) ----------
    if ~isempty(Aset)
        % Extract subproblem
        A_sub = A(:, Aset);
        gamma_sub_init = gamma(Aset);
        
        % Inner loop parameters
        if k <= 2  % CHANGED: from k <= 5 to k <= 3 to k <= 1 to k <= 2
            max_inner = 1;  % One round per iteration in warm-up
        else
            max_inner = 10;  % Reduced from 15 for efficiency
        end
        
        % Solver
        [gamma_sub, cov_inv_new] = solver(...
            A_sub, Sigma_hat, gamma_sub_init, cov_inv, nsPower, ...
            max_inner, epsilon_k);
        
        % Update global solution
        if k <= 2  % CHANGED: from k <= 5 to k <= 3 to k <= 1 to k <= 2
            % During warm-up: full update
            gamma = gamma_sub;
        else
            % After warm-up: update only active set
            gamma_new = zeros(N, 1);
            gamma_new(Aset) = gamma_sub;
            gamma = gamma_new;
        end
        
        % Update covariance inverse
        cov_inv = cov_inv_new;
    end
    
    % ------ Convergence check (Algorithm 1, Line 6) --------------------
    if k > 2  % CHANGED: from k > 5 to k > 3 to k > 1 to k > 2
        % Check relative change
        if exist('gamma_old', 'var')
            rel_change = norm(gamma - gamma_old) / max(norm(gamma_old), 1e-9);
            if rel_change < opt.epsilon
                break;
            end
        end
        gamma_old = gamma;
        
        % Also check if we found K-sparse solution
        if sum(gamma > max(gamma) * 1e-6) <= K && k > 10
            break;
        end
    end
end

% ---------- Algorithm 1, Line 7: Output ---------------------------------
[gam_vals, idx] = sort(gamma, 'descend');
sup_idx = idx(1:K);
gam_vals = gam_vals(1:K);

end

% ========================================================================
function [gamma_new, cov_inv_new] = solver(A_sub, Sigma_hat, ...
    gamma_init, cov_inv_init, nsPower, max_rounds, tol)

[L, n_sub] = size(A_sub);

% Initialize
if n_sub == size(cov_inv_init, 1)  % Full problem (warm-up)
    gamma = zeros(n_sub, 1);
    cov_inv = eye(L) / nsPower;
else  % Subproblem
    gamma = gamma_init;
    cov_inv = cov_inv_init;
end

% Coordinate descent
for round = 1:max_rounds
    ind_perm = randperm(n_sub);
    gamma_old = gamma;
    
    for i = ind_perm
        pilot = A_sub(:, i);
        p_vec = cov_inv * pilot;
        
        num1 = real(p_vec' * (Sigma_hat * p_vec));
        num2 = real(pilot' * p_vec);
        
        if num2 < 1e-12
            continue;
        end
        
        step = max((num1 - num2) / num2^2, -gamma(i));
        
        if abs(step) > 1e-12
            gamma(i) = gamma(i) + step;
            cov_inv = cov_inv - step * (p_vec * p_vec') / (1 + step * num2);
        end
    end
    
    % Convergence check
    if round > 1 && norm(gamma - gamma_old, 1) < tol
        break;
    end
end

gamma_new = gamma;
cov_inv_new = cov_inv;
end