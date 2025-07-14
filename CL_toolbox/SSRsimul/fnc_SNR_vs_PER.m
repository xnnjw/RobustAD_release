function [meanREC,meanTP,meanTime]=fnc_SNR_vs_PER(N,M,L,K,SNR,NRSIM,GaussianS)
%
% meanREC = mean PER rates over SNR levels 
% meanTP  = mean True Positive (TP) rates over SNR levels 
% meanTime = average CPU time over SNR levels 
%
% Esa Ollila, Aalto University
% 
% Simulation set-up from reference: Esa Ollila, "Sparse signal recovery and
% source localization via covariance learning," arXiv:2401.13975 [stat.ME], 
% Jan 2024. https://arxiv.org/abs/2401.13975
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 7
    GaussianS = true;
end

nSNR = length(SNR);
sigmas = zeros(nSNR,K);

% Generate the source signal powers 2nd, 3rd, and 4th have SNR that is
% lower than 1db, 2dB and 4dB than the first source, and so on 
for ii = 1:(K/4)
    sigmas(:,(ii-1)*4+1) = 10.^(SNR.'/10);
    sigmas(:,(ii-1)*4+2) = 10.^((SNR.'-1)/10);
    sigmas(:,(ii-1)*4+3) = 10.^((SNR.'-2)/10);
    sigmas(:,(ii-1)*4+4) = 10.^((SNR.'-4)/10);
end

% NOTE: 10log10(sigmasq) = SNRdb => sigmasq = 10^SNRdb/10

timeiht = zeros(NRSIM,nSNR); 
timeomp = zeros(NRSIM,nSNR); 
timesml = zeros(NRSIM,nSNR);
timebcd = zeros(NRSIM,nSNR);
timeccd = zeros(NRSIM,nSNR);
timecwo = zeros(NRSIM,nSNR);


RECiht  = zeros(NRSIM,nSNR); 
REComp  = zeros(NRSIM,nSNR);
RECsml  = zeros(NRSIM,nSNR);
RECbcd  = zeros(NRSIM,nSNR);
RECccd  = zeros(NRSIM,nSNR);
RECcwo  = zeros(NRSIM,nSNR);

TPiht = zeros(NRSIM,nSNR); 
TPomp = zeros(NRSIM,nSNR);
TPsml = zeros(NRSIM,nSNR);
TPbcd = zeros(NRSIM,nSNR);
TPccd = zeros(NRSIM,nSNR);
TPcwo = zeros(NRSIM,nSNR);

sigmasq0 = 1; % true noise variance

%%
for  jj = 1:length(SNR)
  
    sig = sigmas(jj,:); 
    Lam = diag(sqrt(sig));    
   
    fprintf('\nStarting iterations for sig = %.2f SNR  = %.2f',sig(1), SNR(jj));    
    rng('default');
    
    for iter = 1:NRSIM
    
        %% generate sparse source S   
        if GaussianS
            S = Lam*complex(randn(K,L),randn(K,L))/sqrt(2); %  complex Gaussian data
        else
            S = Lam*exp(-1i*unifrnd(-pi, pi,[K L])); 
        end

        %% generate A, E and Y  
        A = GaussianMM(N,M);                
        E = (sqrt(sigmasq0)/sqrt(2))*complex(randn(N,L),randn(N,L)); 
        loc = sort(randsample(M,K));  % random locations of nonzeros
        Y = A(:,loc)*S + E;  
       
        %% SNIHT 
        tStart = tic;
        [~, Siht] = SNIHT(Y,A,K); 
        tEnd = toc(tStart);
        timeiht(iter,jj) = tEnd;
        TPiht(iter,jj) = K-length(setdiff(Siht,loc));
        if  isempty(setdiff(Siht,loc)) 
            RECiht(iter,jj) = 1;
        end
               
        %% SOMP 
        tStart = tic;
        [~, Somp] = SOMP(Y,A,K); 
        tEnd = toc(tStart);
        timeomp(iter,jj) = tEnd;
        TPomp(iter,jj) = K-length(setdiff(Somp,loc));
        if  isempty(setdiff(Somp,loc)) 
            REComp(iter,jj) = 1;
        end

        %% CL-OMP (proposed)
        tStart = tic;
        Ssml = CLOMP(A,Y,K);
        tEnd = toc(tStart);
        timesml(iter,jj) = tEnd;
        TPsml(iter,jj) = K-length(setdiff(Ssml.',loc));
        if  isempty(setdiff(Ssml,loc)) 
            RECsml(iter,jj) = 1;
        end      
        
        %% CL-BCD
        tStart = tic;
        [Sbcd,~,~,rep] = CLBCD(A,Y,K,false);
        tEnd = toc(tStart);
        timebcd(iter,jj) = tEnd;
        TPbcd(iter,jj)  = K-length(setdiff(Sbcd,loc));
        if  isempty(setdiff(Sbcd,loc)) 
            RECbcd(iter,jj) = 1;
        end

        %% CL-CCD (proposed)
        tStart = tic;
        [Sccd,~,~,rep2]= CLCCD(A,Y,K); 
        tEnd = toc(tStart);
        timeccd(iter,jj) = tEnd;
        TPccd(iter,jj) = K-length(setdiff(Sccd.',loc));
        if  isempty(setdiff(Sccd,loc)) 
            RECccd(iter,jj) = 1;
        end  
        %[timebcd(iter,jj),timeccd(iter,jj)]
        %[TPbcd(iter,jj),TPccd(iter,jj)]

        %% CWOpt 
        tStart = tic;
        [Scwo,~,j1]= CWOpt(A,Y,K,sigmasq0); 
        tEnd = toc(tStart);
        timecwo(iter,jj) = tEnd;
        TPcwo(iter,jj) = K-length(setdiff(Scwo,loc));
        if  isempty(setdiff(Scwo,loc)) 
             RECcwo(iter,jj) = 1;
        end

    end
    fprintf('.')

end

meanREC = [mean(RECiht);mean(REComp); mean(RECsml); mean(RECbcd);mean(RECccd);mean(RECcwo)];
meanTP = [mean(TPiht);mean(TPomp); mean(TPsml); mean(TPbcd);mean(TPccd);mean(TPcwo)]/K;
meanTime = [mean(timeiht);mean(timeomp); mean(timesml); mean(timebcd);mean(timeccd);mean(timecwo)];


