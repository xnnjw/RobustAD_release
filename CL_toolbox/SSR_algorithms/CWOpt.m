function [Ilocs,gamma,j1] = CWOpt(A,Y,K,sig,max_iter)
% function [Ilocs,gamma,j1] = CWOpt(A,Y,K,sig)
%
% Sparse signal reconstruction using covariance learning and coordinatewise 
% optimization scheme described in reference [1,Algorithm 1]. Note that the
% method assumes that noise power (sig) is known and given as input to the 
% algorithm
% 
% INPUT:
%   A       - Dictionary of size N x M
%   Y       - Matrix of N x L (L is the number of measurement vectors)
%   K       - number of non-zero sources
%   sig     - noise variance (positive scalar)
%
% OUTPUT:
%   Ilocs  -  indices of the K-largest signal powers (K-vector)
%   gamma  -  M-vector containing source signal powers
%
% RERENCE: 
%   [1] Fengler et al."Non-Bayesian ac- tivity detection, large-scale 
%   fading coefficient estimation, and unsourced random access with a 
%   massive MIMO receiver," IEEE Transactions on Information Theory, 
%   vol. 67, no. 5, pp. 2925–2951, 2021
%
% (c) Esa Ollila, Aalto University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<5
    max_iter = 500; % maximum number of iterations allowed
end

%% Initialize variables
N = size(A,1);% number of rows in dictionary (= measurement dimensionality)  
M = size(A,2);% number of dictionary entries
L = size(Y,2);% number of snapshots

%-- assert: 
assert(isequal(round(K),K) && K < M && isreal(K) && K >0, ... 
    ['''K'' must be a positive integer smaller than ', num2str(M)]);
RY = (1/L)*Y*(Y'); 

%% Initial value using gamma == 0
SigmaYinv = (1/sig)*eye(N); 
gammaOld = zeros(1,M);
%%
flag = false;   % print report (false, then do not print) 
ERRTOL = 5e-4;
indices = 1:M; % We just go through the coordinates in deterministic order
% Some studies suggest that random order is preferrable for cyclic
% algorithms. This has not been tested. 
gamma = zeros(1,M);

for j1 = 1:max_iter
    
    for j2 = 1:M
      
        indx = indices(j2);
        a = A(:,indx);
        b = SigmaYinv*a;
        c = real(a'*b);
        if c < 10^-12
            c= 10^-12;
        end

        tmp = (real(b'*(RY*b)) - c)/c^2;
        d = max(tmp, -gammaOld(indx));
        gamma(indx) = gammaOld(indx) + d;
        if d~=0
            SigmaYinv = SigmaYinv - (d/(1+d*c))*b*(b');
        end
        
    end

    err1 = norm(gamma-gammaOld,Inf)/norm(gamma,Inf);
    if err1 < ERRTOL
        if flag 
            fprintf('Solution converged!\nIteration: %4d. Error: %.7f\n', j1, err1)
        end
        break; % goodbye     
    end

    gammaOld = gamma;
end

if j1 == max_iter % not convereged
   if flag
       fprintf('Solution not converged. Error: %.6f,', err1); 
   end
end

%% The choose the K-sparse support as indices of K-largest signal powers
[~,ind] = sort(gamma,'descend');
Ilocs = sort(ind(1:K)).';

end
