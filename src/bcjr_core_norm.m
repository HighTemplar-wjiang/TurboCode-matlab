%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Matlab 2013b
% Author: Weiwei Jiang (wjiang@jaist.ac.jp)
% Date: 20141214
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BCJR recursive algorithm
% For systematic convolutoinal codes
% Adopt normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% app_llr: interleaved LLR of extrinsic LLR
% gamma_c: constant parts of gamma (Ck and Lcxy)
% transitions: 2-dimensional array to describe trellis structure
%              FromState, ToState, info, c1, ..., cn
%              Start value is always 1 instead of 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% appinfo_llr: A posterior LLR of information bits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function appinfo_llr = bcjr_core_norm(ap_llr, gamma_c, transitions)
    
    if(min(transitions(:, 1) ~= 1))
        error('Transition unmatched states must be numbered from 1');
    end
    
    info_len  = length(ap_llr); % get information length
    state_num = max(transitions(:, 1)); % get number of states
    trans_num = length(transitions(:, 1)); % get number of transitions
        
    % Allocate memories
    alpha_llr   = zeros(state_num, info_len + 1);
    beta_llr    = zeros(state_num, info_len + 1);
    
    % Initialization
    alpha_llr(1, 1)  = 1;
    beta_llr(:, end) = 1.0 / state_num;
    
    % Bound La
    ap_llr(find(ap_llr > +700)) = +700;
    ap_llr(find(ap_llr < -700)) = -700;
    
    % Get information bits in transitions
    % Expand into all time slots
    uk =         kron(transitions(:, 3), ones(1, info_len));
    app_llr_ex = kron(ap_llr,            ones(trans_num, 1));
    
    % Calculate gamma
    uk = uk * 2 - 1; % 0 -> -1, 1 -> +1
    gamma_llr = gamma_c .* exp(uk .* app_llr_ex * 0.5);
    % gamma_llr = gamma_c + uk .* app_llr_ex * 0.5;
    
    % Calculate alpha
    % Since time index and transition index are from 1
    % Time index must +1 from the equation
    for time_index = 2:info_len+1
        for transition_index = 1:trans_num
            alpha_llr(transitions(transition_index, 2), time_index) = ...
                alpha_llr(transitions(transition_index, 2), time_index) + ...
                    alpha_llr(transitions(transition_index, 1),  time_index - 1) * ...
                        gamma_llr(transition_index, time_index - 1);
        end
        normalizer_alpha = sum(alpha_llr(:, time_index));
        alpha_llr(:, time_index) = alpha_llr(:, time_index) / normalizer_alpha;
    end
    
    % Delete tail brunches
%     alpha_llr(2:end, end) = 0;
%     state_seq = [1:1:state_num]';
%     state_seq(transitions(find(transitions(:, 2) == 1), 1)) = 0;
%     alpha_llr(state_seq(find(state_seq)), end-1) = 0;
        
    % Calculate beta
    for time_rindex = 1:info_len
        time_index  = info_len + 1 - time_rindex;
        for transition_index = 1:trans_num
            beta_llr(transitions(transition_index, 1), time_index) = ...
                beta_llr(transitions(transition_index, 1), time_index) + ...
                    beta_llr(transitions(transition_index, 2), time_index + 1) * ...
                        gamma_llr(transition_index, time_index);
        end
        normalizer_beta = sum(beta_llr(:, time_index));
        beta_llr(:, time_index) = beta_llr(:, time_index) / normalizer_beta;
    end
    
    % Delete head brunches
%     beta_llr(2:end, 1) = 0;
%     state_seq = [1:1:state_num]';
%     state_seq(transitions(find(transitions(:, 1) == 1), 2)) = 0;
%     beta_llr(state_seq(find(state_seq)), 2) = 0;
        
    % Calculate LLR of information bits
    p1 = zeros(trans_num, info_len);
    p0 = zeros(trans_num, info_len);
    for time_index = 1:info_len
        for transition_index = 1:trans_num
            if transitions(transition_index, 3) == 1
                p1(transition_index, time_index) = ...
                    alpha_llr(transitions(transition_index, 1), time_index) * ...
                    gamma_llr(transition_index, time_index) * ...
                    beta_llr(transitions(transition_index, 2), time_index+1);
            else
                p0(transition_index, time_index) = ...
                    alpha_llr(transitions(transition_index, 1), time_index) * ...
                    gamma_llr(transition_index, time_index) * ...
                    beta_llr(transitions(transition_index, 2), time_index+1);
            end
        end
    end
    
    % Normalize LLR
    normalizer_p = sum(p1) + sum(p0);
    p1 = sum(p1) ./ normalizer_p;
    p0 = sum(p0) ./ normalizer_p;
    
    appinfo_llr = log(p1 ./ p0);
    
    % Bound LLR
    appinfo_llr(find(appinfo_llr > +700)) = +700.0;
    appinfo_llr(find(appinfo_llr < -700)) = -700.0;
    
end