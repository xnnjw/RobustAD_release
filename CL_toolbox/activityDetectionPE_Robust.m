function [Pmd, Per, tme, itercwo] = activityDetectionPE_Robust(L, N, K_a, M, P, iter, varargin)
    %ACTIVITYDETECTIONPE Activity Detection with the relaxed ML algorithm
    % Compute the probability of missed detection in a fading MIMO-MAC when
    % the pilots are chosen at random and the activity is detected by
    % picking the positions of the K_a largest estimated coefficients.

    % Input parsing and setup
    if nargin > 6
        for i = 1:2:length(varargin)
            switch varargin{i}
                case 'iter_ml'
                    iter_ml = varargin{i+1};
                case 'fading'
                    fading = varargin{i+1};
                case 'pilot'
                    pilot = varargin{i+1};
            end
        end
    else
        fading.type = 'uniform';
        fading.lower_limit = 10;
        fading.upper_limit = 30;
        pilot = 'unif_circle';
        iter_ml = 15;
    end

    % Large-scale fading
    sampleLSF = @(N) LS_fading(N, fading.type, fading.lower_limit, fading.upper_limit);

    % Initialization of variables for parallel processing
    md = zeros(1, iter);
    md5 = zeros(1, iter);
    md6 = zeros(1, iter);
    md7 = zeros(1, iter);
    per = zeros(1, iter);
    per5 = zeros(1, iter);
    per6 = zeros(1, iter);
    per7 = zeros(1, iter);
    timecwo = zeros(1, iter);
    timemp = zeros(1, iter);
    timemprbst = zeros(1, iter);
    timercwo = zeros(1, iter);
    itercwo = zeros(1, iter);

    for i = 1:iter
        gamma = sampleLSF(N);

        % Create random pilot matrix
        if strcmpi(pilot, 'bernoulli')
            A = (1-2*binornd(1,0.5,[L, N]))/sqrt(2) + 1i*(1-2*binornd(1,0.5,[L, N]))/sqrt(2);
        else
            A = exp(2*pi*1i*rand(L,N)); % matrix with elements that are random variables on the unit circle
        end

        sup = randperm(N,K_a);
        [~,tmp] = sort(gamma(sup),'descend');
        sup = sup(tmp);

        X = zeros(N, M);
        X(sup, :) = diag(sqrt(gamma(sup))) * (randn(K_a, M) + 1i * randn(K_a, M)) / sqrt(2);
        X = sqrt(P) * X;

        % epsCon noise
        % Z = (randn(L,M) + 1i*randn(L,M))/sqrt(2);
        Z = epscont(L,M,1,0.02,0,10);
        Y = A(:, sup) * X(sup, :) + Z;

        % Covariance matrix of the received signal
        cov_m = (1/M) * Y * (Y');

        % Algorithm CWOpt
        tStart = tic;
        % [sup1, ~, j1] = CWOpt(A, cov_m, K_a, 1, 15);
        [sup1, ~, j1] = CWOpt(A, Y, K_a, 1, 50);
        tEnd = toc(tStart);
        timecwo(i) = tEnd;
        itercwo(i) = j1;
        md(i) = 1 - numel(intersect(sup1, sup))/K_a;
        per(i) = (md(i) == 0);

        % Algorithm CLMP
        tStart = tic;
        [sup5, ~] = CLMP(A, cov_m, K_a, 1);
        timemp(i) = toc(tStart);
        md5(i) = 1 - numel(intersect(sup5, sup))/K_a;
        per5(i) = (md5(i) == 0);
        
        % Algorithm RCL-MP
        tStart = tic;
        sup_rcwohub = HubMP(A, Y, K_a, 1, 0.9);
        timercwo(i) = toc(tStart);
        md7(i) = 1 - numel(intersect(sup_rcwohub, sup))/K_a;
        per7(i) = (md7(i) == 0);

        % Algorithm RCWO
        tStart = tic;
        [sup_rclmphub, ~] = HubCWO(A, Y, K_a, 1, 0.9, 20);
        timemprbst(i) = toc(tStart);
        md6(i) = 1 - numel(intersect(sup_rclmphub, sup))/K_a;
        per6(i) = (md6(i) == 0);


    end

    % Calculate probabilities and average times
    Pmd = [mean(md) mean(md5) mean(md6) mean(md7)];
    Per = [mean(per) mean(per5) mean(per6) mean(per7)];
    tme = [mean(timecwo) mean(timemp) mean(timercwo) mean(timemprbst)];
    itercwo = mean(itercwo);
end
