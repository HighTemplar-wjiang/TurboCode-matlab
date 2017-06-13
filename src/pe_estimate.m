%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Matlab 2013b
% Author: Weiwei Jiang (wjiang@jaist.ac.jp)
% Date: 20141221
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate pe from LLRs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% app_llr_d1: a posteriori LLR of uncoded bits from decoder 1
% app_llr_d2: a posteriori LLR of uncoded bits from decoder 2
% threshold: threshold for choosing LLR pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% pe_hat: estimated Pe of two decoders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pe_hat = pe_estimate(app_llr_d1, app_llr_d2, threshold)
    
    info_len = length(app_llr_d1); % get information length
    
    % Select LLRs due to threshold
    selected_llr_index_d1 = abs(app_llr_d1) > threshold;
    selected_llr_index_d2 = abs(app_llr_d2) > threshold;
    selected_llr_index_dd = selected_llr_index_d1 .* selected_llr_index_d2;
    selected_llr_index    = ...
        selected_llr_index_d1(find(selected_llr_index_dd == 1));
    
    selected_llr_d1 = exp(app_llr_d1(selected_llr_index));
    selected_llr_d2 = exp(app_llr_d2(selected_llr_index));
    
    N = double(length(selected_llr_d1)); % get number of LLRs
    if (N == 0)
        pe_hat = -1;
    else
        pe_hat = 1.0/N * ...
            sum((selected_llr_d1 + selected_llr_d1) ./ ...
                ((1 + selected_llr_d1) .* (1 + selected_llr_d2)));
    end
    
end