function [Pmd,Per,tme,itercwo] = activityDetectionPE(L,N,K_a,M,P,iter,varargin)
%ACTIVITYDETECTIONPE Activity Detection with the relaxed ML algorithm
%   Compute the probability of missed detection in a fading MIMO-MAC when
%   the pilots are choosen at random and the activity is detected by
%   picking the positions of the K_a largest estimated coefficients.
%
%   L:      Pilot dimension
%   N:      total users
%   K_a:    active users
%   M:      receive antennas
%   P:      input power
%   iter:   number of simulation runs
%
%   optional inputs:
%   struct fading with the following fields
%   'type' = {'no_fading'(default),'uniform','exp','pathloss','shadowing_pathloss'}
%   'lower_limit'(default = 10),'upper_limit'(default = 30): limits for the uniform case
%   'iter_ml': No. of iterations for the ML algorithm default=10
%
% Orignal version by Alexander Fengler: 
%   https://github.com/AlexFengler/mmv-activity-and-random-access
% Edited by Esa Ollila


    fading.type = 'uniform';
    fading.lower_limit = 10;
    fading.upper_limit = 30;
    pilot = 'unif_circle';
    
    iter_ml = 15;
    
    % Input parsing
    if nargin > 6
        names   = varargin(1:2:end);
        values  = varargin(2:2:end);
        for i=1:length(names)
            switch names{i}
                case 'iter_ml'
                    iter_ml = values{i};
                case 'fading'
                    fading = values{i};
                case 'pilot'
                    pilot = values{i};
            end
        end
    end
    
    % Large-scale fading
    sampleLSF = @(N) LS_fading(N,fading.type, fading.lower_limit, fading.upper_limit);
    % no fading -> all gammas = 1
    md = 0; %md2=0; %md4=0;
    md5=0; %md6=0; md7=0; %md_a =0;
    timecwo = 0; timemp = 0; itercwo=0; %timeomp=0;
    per =0; per5=0;
    onetoN = 1:N;
    for i=1:iter       
        % Create gamma
        gamma = sampleLSF(N);

         % Create random pilot matrix
         if strcmpi(pilot,'bernoulli')
            A = (1-2*binornd(1,0.5,[L, N]))/sqrt(2) + 1i*(1-2*binornd(1,0.5,[L, N]))/sqrt(2);
        else
            A = exp(2*pi*1i*rand(L,N)); % matrix with elements that are random 
            % variables on the unit circle
        end
        % Esa note: orignal code for random unit circle had scaling by 
        %  (1+1i)/sqrt(2) ( = exp(j*pi/4))
        % A(i,j)=exp(j*pi/4)*exp(j*om)=exp(j(pi/4+om)) where omega~Unif(0,2pi)
        % thus A(i,j) = exp(j*om2) where om2~Unif(pi/4,2*pi+pi/4) -> 
        % exp(j*om2) =_d exp(j*om) --> multiplier (1+1i)/sqrt(2) is not needed
        
        % Columns are normalized to L
        norms = sum(abs(A.^2));
        assert(all(abs(norms-L*ones(1,N))<1e-10));

        % Random activity
        sup = randperm(N,K_a);
        % Esa add: then organize the support from  largest to smallest gamma
        % this makes it easier to compare with outputs from CW function. 
        [gam0,tmp] = sort(gamma(sup),'descend');
        sup = sup(tmp);

        % Small-scale fading
        X           = zeros(N,M);
        X(sup,:)    = diag(sqrt(gamma(sup)))*(randn(K_a,M) + 1i*randn(K_a,M))/sqrt(2);
        X           = sqrt(P)*X;
        
        % AWGN
        Z = (randn(L,M) + 1i*randn(L,M))/sqrt(2);
        
        % MMV-MAC
        Y = A(:,sup)*X(sup,:) + Z;

        cov_m = (1/M)*Y*(Y'); 
        %% Esa: This is orignal function by A. Fengler but it is slower than
        % my function CWOpt.
%         tic;
%         [g_hat,it] = ML_coord_descent_round(cov_m, A, iter_ml, 1,[]);
%         toc;
%         [~,sup_est] = maxk(g_hat,K_a);
%         md = md + numel(setdiff(sup_est,sup))/K_a;

        %% CWOpt 
        tStart = tic;
        [sup1,gam1,j1]= CWOpt(A,cov_m,K_a,1,iter_ml); 
        tEnd = toc(tStart);
        timecwo = timecwo + tEnd;
        itercwo = itercwo + j1;
        err1=numel(setdiff(sup1,sup));
        md = md + err1/K_a;
        if err1==0
            per = per + 1;
        end

        %% MP:
        tStart = tic;
        [sup5,gam5] = CLMP_mMTC(A,cov_m,K_a,1);       
        tEnd = toc(tStart);
        timemp = timemp + tEnd;
        err5=numel(setdiff(sup5,sup));
        md5 = md5 + err5/K_a;
        if err5==0
            per5 = per5 + 1;
        end

        %%
        if mod(i,1000)==0
           % fprintf('.')
            [md,md5]/i
        end
    end
    Pmd = [md md5]/iter;
    Per = [per per5]/iter;
    tme = [timecwo timemp]/iter;
    itercwo = itercwo/iter;
end
