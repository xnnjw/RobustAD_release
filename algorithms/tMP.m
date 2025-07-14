function [Ilocs, gamVals] = tMP(A, Y, K, sig2, nu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HubMP: CL Matching Pursuit (MP) using Huber's loss. 
%
% INPUT: 
%   A       - Dictionary of size L x N 
%   Y       - Matrix of L x M 
%   K       - number of non-zero sources
%   sig2    - Noise variance
%   v       - degrees of freedom (v > 2)
%
% OUTPUT:
%   Ilocs    - Support of non-zeros signal powers (K-vector)
%   gamVals  - Power estimates (K-vector)
%
% Esa Ollila, Aalto University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = size(A, 1);                 % Number of pilots
N = size(A, 2);                 % Number of MTDs
M = size(Y, 2);                 % Number of antennas (=measurement vecs)      
Sigmainv = (1/sig2)*eye(L);    % Initial \Sigma^(0)

Ilocs = zeros(1,K);            % A sparse active device set
gamVals = zeros(1,K);        

ufun = @(t, nu) (nu + 2*L) ./ (nu + 2*t);
rhofun = @(t, nu) (nu + 2*L) / 2 .* log(nu + 2*t);

if nargin < 5
    nu = 2.1;
end
integrand = @(t) (t/2) .* (nu + 2*L) ./ (nu + t) .* chi2pdf(t, 2*L);
b = integral(integrand, 0, Inf)/L;

const = (1/(b*M));
onetoN = 1:N;
gam = zeros(1,N);        % Initialize power estimates
I_max = 10;  % maximum number of iterations
delta_thr = 5*10^-3; % FP-iteration threshold
logdet_iSig = -L*log(sig2);

%--- Main Iteration Phase

for k = 1:K

    MD = real(sum(conj(Y).*(Sigmainv*Y))); % y_m'*Sigma^-1*y_m
    B = Sigmainv*A;    % Sigma^-1 a_n , n = 1,..,N
    D = subplus(real(sum(conj(A).*B))); % a_n'Sigma^-1 a_n
    tmp = setdiff(onetoN,Ilocs(1:k-1));
    err = zeros(1,N);
    
    f0 = const*sum(rhofun(MD,nu)) - logdet_iSig;
    
    for i = tmp

        b_i  = B(:,i);
        d_i  = D(i);
        num_i = abs(sum(Y.*conj(b_i))).^2;

        %% test if the solution is gamma(i)=0:
        gam0 = 0.001;
        c_i = 1/(1+gam0*d_i);
        eta = MD - c_i*gam0*num_i;
        f1 = const*sum(rhofun(eta,nu)) - logdet_iSig + log(1+gam0*d_i);
        if f1 > f0
            gam1 = 0;
        else
            % <-- FP iterations start
            gam0 = gam(i);
    
            for t = 1:I_max
    
                c_i = 1/(1+gam0*d_i);
                eta = MD - c_i*gam0*num_i;
                gam1 = (const*sum(ufun(eta,nu).*num_i)-d_i)/(d_i^2);
    
                if gam1==gam0 || abs(gam1 - gam0)/gam1 < delta_thr             
                    break;
                end
                gam0 = gam1;
            end
        end
             
        gam(i) = gam1;        
        % <-- FP iteratations end

        %- Compute the error 
        c_i = gam1/(1+gam1*d_i);
        eta = MD -  c_i*num_i;
        err(i) = sum(rhofun(eta,nu))*const - logdet_iSig + log(1+gam1*d_i);
           
    end

    % Find the index with minimum error 
    [~, indx] = min(err(tmp));
    ix = tmp(indx);
    gamVals(k) = gam(ix);
    Ilocs(k) = ix;

    % Update the Sigma
    if k < K 
        d = gamVals(k); 
        c = real(A(:,ix)'*B(:,ix));
        Sigmainv = Sigmainv - (d/(1+d*c))*B(:,ix)*(B(:,ix)');
    end
    logdet_iSig = -logdet_iSig + log(1+gamVals(k)*D(ix));

end
% <-- Main iteration ends
end

