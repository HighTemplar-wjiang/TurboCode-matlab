%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Matlab 2013b
% Author: Weiwei Jiang (wjiang@jaist.ac.jp)
% Date: 20141210
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Turbo encoder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% rcv_seq: Received sequence
% transitions: Array of transitions of *RSC* encoders
% interleaver: Array interleavers between encoders 
% select_matrix: Matrix for selecting parity check bits
% max_iteration: Max times of iteration
% Lc: Channel reliability measure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% info_seq_est: Estimation of information sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function info_llr = turbo_decoder(rcv_seq, transitions, interleaver, select_matrix, max_iteration, Lc)

    % warning('Turbo decoder is under debugging.');

    rcv_len     = length(rcv_seq); % get length of received sequence
    decoder_num = 2; % get number of decoders
    rsccode_len = length(transitions(1, :)) - 3; % get RSC code length
    tbcode_len  = 1 + sum(select_matrix(:, 1)); % get Turbo code length
    info_len    = rcv_len / tbcode_len; % get length of information
    
    % Construct select sequence from select matrix
    select_seq_array = zeros((rsccode_len - 1) * decoder_num, info_len);
    for info_index = 1:info_len
        select_column = mod(info_index, length(select_matrix(1, :)));
        select_column = select_column + (select_column == 0) * length(select_matrix(1, :));
        select_seq_array(:, info_index) = select_matrix(:, select_column);
%         select_seq_array(:, info_index:decoder_num:end) = ...
%             kron(select_matrix(:, info_index), ones(1, info_len / length(select_matrix(1, :))));
    end
    select_seq_array = [ones(1, info_len) ; select_seq_array]; % add information sequence
    select_seq       = reshape(select_seq_array, 1, prod(size(select_seq_array))); % transform to 1-dimension
        
    % Arrange received information for decoders
    rcv_seq_arrange = zeros(1, length(select_seq));
    rcv_seq_arrange(find(select_seq == 1)) = rcv_seq;
    rcv_seq_array = reshape(rcv_seq_arrange, size(select_seq_array));
    
    % Construct deinterleavers
    deinterleaver(interleaver) = 1:1:info_len;

    
    % Decoder properties
    Ck = 1; % Conditional LLR factor
    % fading_amp = 1; % Channel fading amplitude
    
    % Initialize decoder variables
    % Lc = 4 * fading_amp * 10 ^ (Ec_No_dB/10.0); % Channel reliability measure
    % Lc = 1.0; % test
    gamma_c_array = zeros(length(transitions(:, 1, 1)), info_len, decoder_num); % Constant part of gamma

    % Initialize gamma_c
    rcv_info_seq = rcv_seq_array(1, :);
    ys = zeros(1, info_len); % received infomation and parity check sequence
    xk = (transitions(:, 4:end) * 2 - 1)'; % encoder ouput, 0 -> -1, 1 -> +1
    ys(1, :) = rcv_info_seq;
    for decoder_index = 1:decoder_num
        ys(2, :) = rcv_seq_array(1+decoder_index:rsccode_len+decoder_index-1, :); % parity check sequence
        for time_index = 1:info_len
            for transition_index = 1:length(transitions(:,1))
                gamma_c_array(transition_index, time_index, decoder_index) = log(Ck) + ...
                    (Lc * 0.5 * sum(xk(:, transition_index) .* ys(:, time_index), 1));
                % gamma_c_array(transition_index, time_index, decoder_index) = Ck * ...
                %    exp(Lc * 0.5 * sum(xk(:, transition_index) .* ys(:, time_index), 1));
            end 
        end
        ys(1, :) = rcv_info_seq(interleaver); % interleaver for decoder 2
    end
    
    % Debug
    app_llr_array1 = zeros(max_iteration, info_len);
    app_llr_array2 = zeros(max_iteration, info_len);
    
    % Decoding process
    ap_llr2 = zeros(1, info_len); % La: L-apriori for decoder 1
    channel_info_seq = rcv_seq_array(1, :); % get received information sequence
    for iteration_index = 1:max_iteration
        
        % Decoder 1
        channel_llr1   = Lc * channel_info_seq; % channel L value
        bcjr_app_llr1  = bcjr_core(ap_llr2, gamma_c_array(:, :, 1), transitions);
        % bcjr_app_llr1  = bcjr_core_norm(ap_llr2, gamma_c_array(:, :, 1), transitions);
        extrinsic_llr1 = bcjr_app_llr1 - ap_llr2 - channel_llr1; % Calculate extrinsic LLR
        ap_llr1        = extrinsic_llr1(interleaver); % interleave for decoder 2
        % ap_llr1 = ap_llr2;
        
        % Decoder 2
        channel_llr2   = channel_llr1(interleaver);
        bcjr_app_llr2  = bcjr_core(ap_llr1, gamma_c_array(:, :, 2), transitions);
        % bcjr_app_llr2  = bcjr_core_norm(ap_llr1, gamma_c_array(:, :, 2), transitions);
        extrinsic_llr2 = bcjr_app_llr2 - ap_llr1 - channel_llr2;
        ap_llr2        = extrinsic_llr2(deinterleaver);
        
        % Debug
        app_llr_array1(iteration_index, :) = bcjr_app_llr1;
        app_llr_array2(iteration_index, :) = bcjr_app_llr2;

    end
    
    info_llr = bcjr_app_llr2(deinterleaver);
    
end