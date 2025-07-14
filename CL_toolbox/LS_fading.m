function [ gamma ] = LS_fading( N, type, lower_limit, upper_limit )
%LS_FADING Summary of this function goes here
%   Detailed explanation goes here
% Code by Alexander Fengler: https://github.com/AlexFengler/mmv-activity-and-random-access
    if ~exist('upper_limit','var')
        upper_limit = Inf;
    end
    if ~exist('lower_limit','var')
        lower_limit = 0;
    end
    switch type
        case 'uniform'
            lower_limit_dB = lower_limit;
            upper_limit_dB = upper_limit;
            gamma_dB = lower_limit_dB + (upper_limit_dB - lower_limit_dB)*rand(N,1);
            gamma = 10.^(gamma_dB/10);
        case 'no_fading'
            gamma = ones(N,1);
        case 'exp'
            gamma = randn(N,1).^2;
        case 'pathloss'
            alpha       = 128.1;
            beta        = 36.7;
            R           = 1; %km
            % Assume p(d) = 2d/R^2;
            d = R*sqrt(rand(N,1));               

            % path loss in DB normalized with noise density=170dBm/Hz and
            % 10Mhz BW
            z = -(alpha+beta*log10(d)) + 170 - 70+30;
            
            gamma = 10.^(z./10);

        case 'shadowing_pathloss'
            alpha       = 128.1;
            beta        = 36.7;
            R           = 1; %km
            sigma_SF2   = 8;
            % Assume p(d) = 2d/R^2;
            d = R*sqrt(rand(N,1));               
            
            %shadowing
            s = sqrt(sigma_SF2)*randn(N,1);
            
            % path loss in DB normalized with noise density=170dBm/Hz and
            % 10Mhz BW
            z = -(alpha+beta*log10(d)) - s + 170 - 70;
            
            gamma = 10.^(z./10);
    end

end

