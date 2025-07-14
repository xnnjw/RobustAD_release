function [Ilocs,gam,sig,report] = CLBCD(A,Y,K,peaksparse,max_iter,print_flag)
% [Ilocs,gam,sig,report] = CLBCD(A,Y,K,peaksparse)
%
% Sparse signal reconstruction (SSR) algorithm using covariance learning 
% and Block Coordinate Descent (BDC) for updating the signal powers and the 
% noise power. For signal power, fixed point (FP) algorithm is used. See 
% reference Ollila (2024) for details. 
% 
% INPUT:
%   A          - Dictionary of size N x M
%   Y          - Matrix of N x L (L is the number of measurement vectors, 
%                   N is the dimension)
%   K          - number of non-zero sources or peaks
%   peaksparse - flag (false or true). Default = True. 
%   max_iter   - max number of iteratations. Default = 500
%   print_flag - Print some output (Default = False)
%
% OUTPUT:
%   Ilocs   -  support of non-zeros signal powers (K-vector)
%   gam     -  K x 1 vector containing non-zero source powers
%   sig     -  noise variance estimate 
%   report  -  a structure array with some stuff
%
% REFERENCE:
%   Esa Ollila, "Sparse signal recovery and source localization via 
%   covariance learning," arXiv:2401.13975 [stat.ME], Jan 2024. 
%   https://arxiv.org/abs/2401.13975
% 
% (c) Esa Ollila, Aalto University, 2024.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<6
    print_flag = false;   % print report (false, then do not print) 
end

if nargin<5
    max_iter = 500; % maximum number of iterations allowed
end

if nargin<4
    peaksparse = true;
end

%% Initialize variables
N = size(A,1);% number of rows in dictionary (= measurement dimensionality)  
M = size(A,2);% number of dictionary entries
L = size(Y,2);% number of snapshots

%-- assert: 
assert(isequal(round(K),K) && K < M && isreal(K) && K >0, ... 
    ['''K'' must be a positive integer smaller than ', num2str(M)]);
assert(islogical(print_flag) || print_flag == 0 || print_flag==1,['''print_flag'' must be ' ...
    'logical or 1 or 0']);
assert(islogical(peaksparse) || peaksparse == 0 || peaksparse==1,['''peaksparse'' must be ' ...
    'logical or 1 or 0']);

RY = (1/L)*Y*(Y'); % The sample covariance matrix (SCM)

%% Initialize
sig = real(trace(RY)/N);
SigmaYinv   = (1/sig)*eye(N); 
gammaOld = zeros(1,M);
IlocsOld = [];
%%
ERRTOL = 5e-4;

for j1 = 1:max_iter

    B =  SigmaYinv*A; % \Sigma^-1 a_m , m=1,..,M
    gamma_num = subplus(real(sum(conj(B).*(RY*B))));
    gamma_denum = subplus(real(sum(conj(A).*B)));
    gamma_denum(gamma_denum<=10^-12) = 10^-12; % make sure not zero
    gamma = (gamma_num./gamma_denum.^2) + subplus(gammaOld - 1./gamma_denum);
    
    if ~all(gamma>0) && print_flag
        fprintf('Warning: Some of the gamma-s are not positive!')
    end

    if peaksparse 
        [~, Ilocs] = find_peaks_1D(gamma,K);
    else
        [~, indx] = sort(gamma,'descend');
        Ilocs = indx(1:K);
    end

    if ~isempty(setdiff(IlocsOld,Ilocs)) || j1==1
        % We need to reupdate the noise power estimate only if support has
        % changed 
        Am = A(:,Ilocs);        
        P_N = eye(N)-Am*pinv(Am);
        sig = real(trace(P_N*RY))/(N - K); 
    end

    err1 = norm(gamma-gammaOld,Inf)/norm(gamma,Inf);
    if err1 < ERRTOL
        if print_flag 
            fprintf('Solution converged!\nIteration: %4d. Error: %.7f\n', j1, err1)
        end
        break;   
    end

    SigmaYinv   =  (sig*eye(N) + bsxfun(@times,gamma,A)*A')\eye(N); 
    gammaOld = gamma;
    IlocsOld = Ilocs;
end

if j1 == max_iter % not convereged
   if print_flag
       fprintf('Solution not converged. Error: %.6f,', err1); 
   end
end

[~,idx] =  sort(Ilocs);
Ilocs = Ilocs(idx);
report.j1 = j1;
report.gamma = gamma;
gam = gamma(Ilocs);

end


function [pks, locs] = find_peaks_1D(gamma, Nsources)
% This is the code SBLpeaks_1D with different name:
% 
% [pks, locs] = SBLpeaks_1D(gamma, Nsources)
%
% fast alternative for findpeaks in 1D case
%

% output variables
pks = zeros(Nsources,1);
locs = zeros(Nsources,1);

% zero padding on the boundary
gamma_new = zeros(length(gamma)+2,1);
gamma_new(2:end-1) = gamma;

[~, Ilocs]= sort(gamma,'descend');

% current number of peaks found
npeaks = 0;

for ii = 1:length(Ilocs)
    
    % local patch area surrounding the current array entry i.e. (r,c)
    local_patch = gamma_new(Ilocs(ii):Ilocs(ii)+2);
    
    % zero the center
    local_patch(2) = 0;
    
    if sum(sum(gamma(Ilocs(ii)) > local_patch)) == 3
        npeaks = npeaks + 1;
        
        pks(npeaks) = gamma(Ilocs(ii));
        locs(npeaks) = Ilocs(ii);
        
        % if found sufficient peaks, break
        if npeaks == Nsources
            break;
        end
    end
    
end

% if Nsources not found
if npeaks ~= Nsources
    pks(npeaks+1:Nsources) = [];
    locs(npeaks+1:Nsources) = [];
end

end