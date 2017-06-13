%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Matlab 2013b
% Author: Weiwei Jiang (wjiang@jaist.ac.jp)
% Date: 20141210
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Turbo encoder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% info_seq: Information sequence
% transitions: Array of transitions of *RSC* encoders
% interleaver: Array interleavers between encoders 
% select_matrix: Matrix for selecting parity check bits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% encoded_seq: Encoded bit sequences, one information to one column
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function encoded_seq = turbo_encoder(info_seq, transitions, interleaver, select_matrix)

    info_len    = length(info_seq); % get informatin length
    encoder_num = 2; % get number of encoders
    rsccode_len = length(transitions(1, :)) - 3; % get RSC code length
    
    if info_len ~= length(interleaver)
        error('Information and interleaver length unmatched');
    end
    
    % Encode information sequence
    code1 = conv_encoder(info_seq, transitions);
    code2 = conv_encoder(info_seq(interleaver), transitions);
    parity_seq_array = [code1(2:end, :) ; code2(2:end, :)];
    
    % Construct select sequence from select matrix
    select_seq_array = zeros((rsccode_len - 1) * encoder_num, info_len);
    for info_index = 1:info_len
        select_column = mod(info_index, length(select_matrix(1, :)));
        select_column = select_column + (select_column == 0) * length(select_matrix(1, :));
        select_seq_array(:, info_index) = select_matrix(:, select_column);
        % select_seq_array(:, info_index:encoder_num:end) = ...
        %     kron(select_matrix(:, info_index), ones(1, info_len / length(select_matrix(1, :))));
    end
    
    % Add information sequence
    select_seq_array = [ones(1, info_len) ; select_seq_array];
    code_seq_array   = [info_seq ; parity_seq_array];
    
    % Reshape to a one row sequence, remember that reshape() is low-dimension-major!
    % eg. page -> column -> row
    select_seq  = reshape(select_seq_array, 1, prod(size(select_seq_array)));
    code_seq    = reshape(code_seq_array, 1, prod(size(code_seq_array)));
    encoded_seq = code_seq(find(select_seq == 1));
    
end