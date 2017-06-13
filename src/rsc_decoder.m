%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Matlab 2013b
% Author: Weiwei Jiang (wjiang@jaist.ac.jp)
% Date: 20141218
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RSC encoder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% rcv_seq: Received sequence
% Lc: Channel L value
% transitions: Array of transitions of *RSC* encoders
% ap_llr: La
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% app_llr: LLR of information sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function app_llr = rsc_decoder(rcv_seq, Lc, transitions, ap_llr)

    rcv_len  = length(rcv_seq); % length of received sequence
    code_len = length(transitions(1, 4:end)); % length of a codeword
    info_len = rcv_len / code_len; % information length
    
    rcv_array = reshape(rcv_seq, code_len, info_len); % reshape received signals
    % rcv_info_seq = rcv_array(1,     :); % extract information
    % rcv_prty_seq = rcv_array(2:end, :); % extract parity check bits
    
    % Construct constant part of gamma
    gamma_c = zeros(length(transitions(:, 1)), info_len);
    xk = (transitions(:, 4:end) * 2 - 1)'; % encoder ouput, 0 -> -1, 1 -> +1
    
    for time_index = 1:info_len
        for transition_index = 1:length(transitions(:,1))
            gamma_c(transition_index, time_index) = ...
                Lc * 0.5 * sum(xk(:, transition_index) .* rcv_array(:, time_index), 1);
            % gamma_c(transition_index, time_index) = ...
            %     exp(Lc * 0.5 * sum(xk(:, transition_index) .* rcv_array(:, time_index), 1));
        end 
    end
    
    % app_llr = bcjr_core_norm(ap_llr, gamma_c, transitions);
    app_llr = bcjr_core(ap_llr, gamma_c, transitions);

end

