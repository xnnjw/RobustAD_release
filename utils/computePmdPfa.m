function [pmd, pfa] = computePmdPfa(sup_test, sup, N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to compute both Probability of Missed Detection (PMD)
% and Probability of False Alarm Rate (PFA)
%
% Inputs:
%   sup_test - Indices detected as active by the algorithm
%   sup      - True indices of the actual active devices
%   N        - Total number of devices (active + inactive)
%
% Outputs:
%   pmd      - Probability of Missed Detection
%   pfa      - Probability of False Alarm Rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Find the elements in sup that are correctly detected
    common_elements = intersect(sup_test, sup);
    
    % Compute the number of missed detections
    num_missed = numel(sup) - numel(common_elements);
    % Compute the PMD
    pmd = num_missed / numel(sup);

    % Compute the number of false alarms (elements in sup_test not in sup)
    num_false_alarms = numel(sup_test) - numel(common_elements);
    % Compute the PFA
    num_inactive = N - numel(sup);
    pfa = num_false_alarms / num_inactive;

end
