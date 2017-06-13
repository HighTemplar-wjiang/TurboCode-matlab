%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Matlab 2013b
% Author: Weiwei Jiang (wjiang@jaist.ac.jp)
% Date: 20141210
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transform convolution code polynomial into Trellis transitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% generator_polynomial: [[H1] ; [H2] ; ... ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% transitions: 2-dimensional array to describe trellis structure
%              FromState, ToState, input, c1, ..., cn
%              Start value is always 1 instead of 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function transitions = polynomial2trellis(generator_polynomial)
    
    code_len = length(generator_polynomial(:, 1)); % get code length
    mem_num  = floor(sqrt((max(max(generator_polynomial))))); % get number of memories
    
    % Memory allocation
    transitions   = zeros(2 ^ (mem_num + 1), 2 + code_len);
    
    % Translate generator polynomial into binary domain
    generator_bs  = fliplr(dec2bin(generator_polynomial));
    generator_bin = reshape(str2num(reshape(generator_bs, prod(size(generator_bs)), 1)), ...
                            [size(generator_polynomial) mem_num+1]);
    
    % Calculate transitions
    % Number of transitions = 2 * number of memory states
    % Memory state is numbered from 1 but 0
    for memory_state = 1:2^mem_num 
        transition_index = memory_state * 2 - 1;
        % From state
        transitions(transition_index,     1) = memory_state;
        transitions(transition_index + 1, 1) = memory_state;
        
        % Set memory state
        bit_memory = reshape(str2num(reshape(dec2bin(memory_state-1, mem_num), mem_num, 1)), ...
                             1, mem_num);
        % bit_memory = (bit_memory); % from left to right
        
        % Calculate encode bits
        for code_index = 1:code_len
            output_polynomial = reshape(generator_bin(code_index, 1, :), 1, mem_num+1);
            recsv_polynomial  = reshape(generator_bin(code_index, 2, :), 1, mem_num+1);
            
            % For input = 0
            recsv_bit  = floor(mod(sum(recsv_polynomial  .* [0 bit_memory]),         2));
            output_bit = floor(mod(sum(output_polynomial .* [recsv_bit bit_memory ]), 2));
            next_state = [recsv_bit bit_memory(1, 1:end-1)];
            
            transitions(transition_index, 2)            = bin2dec(num2str(next_state)) + 1;
            transitions(transition_index, 3)            = 0;
            transitions(transition_index, 3+code_index) = output_bit;
            
            % For input = 1
            recsv_bit  = floor(mod(sum(recsv_polynomial  .* [1 bit_memory]),         2));
            output_bit = floor(mod(sum(output_polynomial .* [recsv_bit bit_memory]), 2));
            next_state = [recsv_bit bit_memory(1, 1:end-1)];
            
            transitions(transition_index+1, 2)            = bin2dec(num2str(next_state)) + 1;
            transitions(transition_index+1, 3)            = 1;
            transitions(transition_index+1, 3+code_index) = output_bit;
        end 
    end
end