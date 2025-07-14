function [Ilocs,gam] = CLMP_mMTC(A,RY,K,sigc)
% function [Ilocs,gam0,sigc,gam] = CLOMP(A,Y,K)
%
% Sparse signal reconstruction (SSR) algorithm using greedy covariance 
% learning described in referOllila (2024).
%
% INPUT: 
%   A       - Dictionary of size N x M
%   RY       - Matrix of N x L (L is the number of measurement vectors)
%   K       - number of non-zero sources
%   sig    -  noise variance  (positive scalar)
%
% OUTPUT:
%   Ilocs  -  support of non-zeros signal powers (K-vector)
%   gam0   -  Kx1 vector containing non-zero source signal powers
%
% REFERENCE:
%   Esa Ollila, ""
% 
% Author: Esa Ollila, Aalto University, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize variables
N = size(A,1);% number of sensors
M = size(A,2);% number of dictionary entries
L = size(RY,2);% number of snapshots in the data covariance

%-- assert: 
assert(isequal(round(K),K) && K < M && isreal(K)    && K >0, ... 
    ['''K'' must be a positive integer smaller than ', num2str(M)]);
assert(isreal(sigc) && sigc>0,['''sig'' must be positive']);


%% Initialize
SigmaYinv = (1/sigc)*eye(N); 
SigmaY= sigc*eye(N);
Ilocs = zeros(1,K);
gam0 = zeros(K,1);

for k = 1:K
    
    %% 1. go through basis vectors
    B =  SigmaYinv*A; % \Sigma^-1 a_m , m=1,..,M
    gamma_num = subplus(real(sum(conj(B).*((RY-SigmaY)*B))));
    gamma_denum = subplus(real(sum(conj(A).*B)));
    gamma_denum(gamma_denum<=10^-12) = 10^-12; % make sure not zero
    gamma = gamma_num./gamma_denum.^2;   
    tmp  = gamma.*gamma_denum;
    fk = log(1+tmp) - tmp;
    
    %% 2. Update the index set 
    tmp = setdiff(1:M,Ilocs(1:k-1));
    [~,indx] = min(fk(tmp));
    indx = tmp(indx);
    Ilocs(k) = indx;

    %% 3. Update gamma0  
    if k==1
        gam(k) = gamma(indx)*(N/(N-1));
    else
        gam(k) = gamma(indx);
    end
    %% 4. Update the covariance matrix
    d = gam(k);
    a = A(:,indx);
    b = SigmaYinv*a;
    c = real(a'*b);
    SigmaYinv =  SigmaYinv - (d/(1+d*c))*b*(b');

    %% check ML condition
%     if check_cond
%         B =  SigmaYinv*Am; % \Sigma^-1 a_m , m=1,..,M
%         gamma_cond = real(sum(conj(B).*((RY-SigmaY)*B)));
%         %fprintf('gamma_cond = %.6f',gamma_cond);
%     end

end

%[gam2,Ilocs2] = maxk(gamma,K); 

%[~,idx] =  sort(Ilocs);
%Ilocs2 = Ilocs(idx);
% gam0 = gam0(idx).';

end

