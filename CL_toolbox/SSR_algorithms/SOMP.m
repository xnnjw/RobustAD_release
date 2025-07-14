function [X,S] = SOMP(Y,A,K)
% [X,S] = SOMP(Y,A,K)
%
% Simultaneous Orthogonal Matching Pursuit algorithm 
% 
% INPUT  
%   Y  - matrix of  L  measurement vectors 
%   A  - measurement matrix 
%   K  - the number of nonzero signals
%
% OUTPUT  
%   X  - estimated K-rowsparse signal matrix 
%   S  - estimated support set of non-zeros
% 
% Author: E. Ollila, Aalto University 
% Note: not a fast implementation of SOMP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,L]=size(Y);
X = zeros(K,L); % The values of non-zero signals 
R = Y;
S = zeros(1,K); % initial support 

for k=1:K
    E = sum(abs(A'*R),2);  % sum of absolute correlations (L1-norm of rows)
    [~,pos] = max(E);   % find index with max correlation
    S(k) = pos;         % update support 
    V = A(:,S(1:k));
    X(1:k,:) = pinv(V)*Y; % LSE
    R = Y- V*X(1:k,:);  % update residual
end


