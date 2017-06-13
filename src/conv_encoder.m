%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Matlab 2013b
% Author: Weiwei Jiang (wjiang@jaist.ac.jp)
% Date: 20141211
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolutional encoder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% info_seq: Information sequence
% transitions: 2-dimensional array to describe trellis structure
%              FromState, ToState, info, c1, ..., cn
%              Start value is always 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% encoded_seqence: Encoded bit sequences, one information to one column
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function encoded_seq = conv_encoder(info_seq, transitions)

    info_len = length(info_seq); % get information length
    code_len = length(transitions(1, :)) - 3; % get code length
    
    mem_state   = 1; % initialize memory state
    encoded_seq = zeros(code_len, info_len); % initialize output array
    
    % Use Trellis to encode information sequence
    for info_index = 1:info_len
        transition_index = (mem_state - 1) * 2 + 1 + info_seq(info_index);
        encoded_seq(:, info_index) = transitions(transition_index, 4:end)';
        mem_state = transitions(transition_index, 2);
    end
    
end