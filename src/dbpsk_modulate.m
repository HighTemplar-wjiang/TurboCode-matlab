%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Matlab 2013b
% Author: Weiwei Jiang (wjiang@jaist.ac.jp)
% Date: 20141223
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DBPSK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% seq: sequence to be modulated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% modulated_seq: modulated sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function modulated_seq = dbpsk_modulate(seq)

    current_phase    = seq(1) * 2 - 1; % set first phase
    modulated_seq    = zeros(1, length(seq));
    modulated_seq(1) = current_phase;    
    for bit_index = 2:length(seq)
        % Inverse phase while current bit is different from the last one
        current_phase = current_phase * (-1)^seq(bit_index);
        % Modulate bit
        modulated_seq(bit_index) = current_phase;
    end

end