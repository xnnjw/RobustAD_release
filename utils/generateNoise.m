function samples = generateNoise(N, M, dist_type, params, sigma2)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generateNoise generates noise samples with a specified variance.
%
% INPUTS:
%   N         - Number of rows (users or data points)
%   M         - Number of columns (antennas or features)
%   dist_type - Type of distribution to generate:
%               'gaussian': Standard complex Gaussian noise
%               't': Student's t-distribution noise
%               'epscont': epsilon-contaminated Gaussian noise
%   params    - Structure containing distribution-specific parameters:
%               For 't': 
%                 params.nu - Degrees of freedom for t-distribution
%               For 'epscont':
%                 params.epsilon - Probability of contamination
%                 params.lambda - Outlier noise strength factor
%   sigma2    - Desired variance for the output samples
%
% OUTPUT:
%   samples   - Generated noise samples of size N x M with specified variance
% Note:
%       I want to note that, for epscont noise, the scaling is already done by
%       calculating the sigma1_sq using the equation:
%       \sigma^2 = (1 - \epsilon + \epsilon\lambda^2)\sigma_1^2
%       Thus, the samples generated are with specified sigma2.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Generate base complex normal samples
    normal_samples = (randn(N, M) + 1i * randn(N, M)) / sqrt(2);
    
    switch dist_type
        case 'gaussian'
            samples = sqrt(sigma2) * normal_samples;

        case 't'
            nu = params.nu;
            if nu <= 2
                error('nu needs to be larger than 2!');
            end
            chi2_samples = chi2rnd(nu, [N, 1]); % Chi-squared samples for each row
            scale_factors = sqrt((nu-2)./ chi2_samples); % Scale factors
            samples = normal_samples .* scale_factors; % unit variance t_v(0,1) samples           
            samples = sqrt(sigma2)*samples; % Adjust to match desired variance

        case 'epscont'
            epsilon = params.epsilon;  lambda = params.lambda; 
            % The sigma1_sq is for two Gaussian noises.
            sigma1_sq = sigma2 / (1 - epsilon + epsilon * lambda^2);

            clean_samples = normal_samples * sqrt(sigma1_sq); % Clean noise
            outlier_samples = normal_samples * sqrt(lambda^2 * sigma1_sq); % Outlier noise
            
            contamination_mask = rand(N, M) < epsilon;
            samples = clean_samples .* (1 - contamination_mask) + outlier_samples .* contamination_mask;

        otherwise
            error('Unknown distribution type');
    end
end