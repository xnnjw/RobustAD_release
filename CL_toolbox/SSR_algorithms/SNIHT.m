function [Xnew, X1supp, failure, iter] = SNIHT(Y,A,k,X0supp,itermax,printitn)
% [Xnew, X1supp, failure, iter] = SNIHT(Y,A,k,X0supp,printitn)
%
% Simultaneous Normalized Iterative Hard Thresholding algorihtm proposed 
% in Blanchard et al (2014). 
% 
% INPUT  
%   Y  - matrix of measurement vectors 
%   A  - measurement matrix
%   k  - the number of nonzero coefficients 
%   printitn - print iteration number (printitn=0 (default) if not print )
%
% OUTPUT  
%   Xnew    - estimated signal K-rowsparse signal matrix 
%   X1supp  - estimated support set of non-zeros
%   failure - equal to 1 if algorithm failed
% 
% Author: E. Ollila, Aalto University 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, n] = size(A);
[~, q] = size(Y);

if nargin < 6 
    printitn=0; 
end

if nargin<5
    itermax = 500; % maximum number of iterations allowed
end
if nargin < 4 
    X0supp=[]; 
end

%%-- initial approximation is the zero matrix 
X0  = zeros(n,q); 

%%  DetectSupport
if isempty(X0supp)
    % if initial set of support is not given 
    R = A'*Y;
    [~, indx] = sort(sum(R.*conj(R),2),'descend');
    X0supp = indx(1:k);
end

ERRORTOL = 1e-8; % ERROR TOLERANCE FOR HALTING CRITERION
objold  = Inf; % a big number
failure = 0;
R = Y;

for iter = 1:itermax
    
    %% Compute the negative gradient 
    G = A'*R; 

    %%- stepsize
    mu = norm(G(X0supp,:),'fro')/norm(A(:,X0supp)*G(X0supp,:),'fro');
    mu = mu^2;

    %%-- Next proposal 
    X1 = X0 + mu*G;
    %%-- Detect Support
    [~, indx] = sort(sum(X1.*conj(X1),2),'descend');
    X1supp = indx(1:k);
    
    %%-- Threshold
    Xnew   = zeros(n,q);
    Xnew(X1supp,:)= X1(X1supp,:);
     
    R = Y - A*Xnew;
    objnew =   sum(abs(R(:).^2));  % = norm(Y-A*Xnew,'fro')^2
    
    if mod(iter,printitn) == 0
       fprintf('%3d : mu = %.2f objnew = %.2f\n',iter,mu,objnew);  
    end   
   
    if objnew > objold
       
      while objnew > objold
         
         % Reduce mu
         mu = mu/2; 
 
         % Update the proposal for b 
         X1 = X0 + mu*G;
         %%-- Detect Support
         [~, indx] = sort(sum(X1.*conj(X1),2),'descend');
         X1supp = indx(1:k);
    
         %%-- Threshold
         Xnew   = zeros(n,q);
         Xnew(X1supp,:)= X1(X1supp,:);
         
         % Update objnew
         R = Y - A*Xnew;
         objnew =   sum(abs(R(:).^2));  % = norm(Y-A*Xnew,'fro')^2
         
      end
      
      if mod(iter,printitn) == 0
          fprintf('updating mu = %.2f --> objnew = %.2f\n',mu,objnew);  
      end
    end 
      
    %%-- Stopping criteria          
    crit = norm(Xnew-X0,'fro')^2/norm(Xnew,'fro')^2;

    if mod(iter,printitn) == 0
        fprintf('%3d : mu = %.2f objnew = %.2f crit = %e\n',iter,mu,objnew,crit);  
        %    eq1 = A(:,X1supp)'*R;
        %    crit2 = norm(eq1,'fro');
        %    fprintf('%3d : mu = %.2f objnew = %.2f crit = %e crit2 = %.4f\n',iter,mu,objnew,crit,crit2);  
    end   
       
    if crit < ERRORTOL
        break;
    end 
        
    X0   = Xnew;
    X0supp = X1supp;
    objold = objnew; 
     
end

if iter == itermax
    failure=true; 
end

X1supp = sort(X1supp);
