function [AER,pfa, SSR] = sblBasicFunc(y,A,xTrue,alphaN,sig2e,itrMax,Tol)
% iterations = 100;
lambda_SBL = 1*sig2e;   % Noise variance
A = A;
y = y;
M = size(A,1);
N = size(A,2);
L = size(y,2);

x = zeros(N,L);
i = 1;
condition = 10;
alphaEst = zeros(N,1); 
gamma = .1*ones(N,1);
% SBL-EM
tic;  
while (condition>Tol&&i<itrMax)
    x_prev = x;
    Gamma = diag(gamma);
    G = diag(sqrt(gamma));
    [U,S,V] = svd(A*G,'econ');
    [d1,~] = size(S);
    if (d1 > 1)     diag_S = diag(S);
    else            diag_S = S(1);      end
    Xi = G * V * diag((diag_S./(diag_S.^2 + sig2e + 1e-16))) * U';
    mu_x = Xi * y;
    x = mu_x;

    %% 1) Traditional EM SBL
    %     Sigma_x = real( gamma - (sum(Xi.'.*(A*Gamma)))');
    %     gamma = sum(abs(mu_x).^2,2)/L + Sigma_x;

    %% 2) Fast EM SBL
    R_diag = max(0,real( (sum(Xi.'.*A)).' ));
    mu2_bar = sum(abs(mu_x).^2,2);
    gamma = max(0,sqrt( gamma.*real(mu2_bar./(L*R_diag+1e-10)) ));

    %% 3) MacKay fixed-point SBL
    %     R_diag = max(0,real( (sum(Xi.'.*A)).' ));
    %     mu2_bar = sum(abs(mu_x).^2,2);
    %     gamma = mu2_bar./(L*R_diag+1e-10);

    i = i + 1;
    
    condition = (norm(x_prev-x,'fro'))/norm(x,'fro');    
end
timeVal = toc;  
xHat = x;
%reducedNMMSEFull = (norm(xTrue-xHat,'fro'))/norm(xTrue,'fro');
%% NMSE ESTIMATE FOR SBL
[AER,SSR, pfa] = performanceMatrix(xHat,alphaN,N); 
end


function [pmd,SSR,pfa] = performanceMatrix(xHat,alphaN,N)
%% INITIALIZATION  
K = sum(alphaN); 
actSet = find(alphaN == 1);  
detectVals = vecnorm(xHat,2,2);
maxValuesPower = maxk(detectVals,K); 
[~,approxSupport] = intersect(detectVals,maxValuesPower,'stable'); 

%% EXACT RECOVERY 
if isempty(setdiff(actSet,approxSupport))
    sblSRR=1;
else
    sblSRR=0;
end
SSR = sblSRR;

%% PROBABILITY OF MISS DETECT
numberK = numel(setdiff(actSet,approxSupport)); 
pmd = numberK/K;  

%% PROBABILITY OF FALSE ALARM RATE
false_alarms = setdiff(approxSupport, actSet); % Detected active indices that are not actually active
pfa = numel(false_alarms) / (N - K); % Number of false alarms divided by the number of inactive devices

end 








