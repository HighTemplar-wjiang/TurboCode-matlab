% Doped accumulator simulation

close all;
clear all;

% Simulation parameters
info_len   = 65536; % information length
block_num  = 1000; % block numbers
ebn0_array = 0.0:0.1:5.0; % SNR array in dB
err_array  = zeros(1, length(ebn0_array)); % Bit error statistic
ber_array  = double(zeros(1, length(ebn0_array))); % BER statistic

% Encoder attributes
doping_rate = 2; % Doping rate of DACC
% transitions = [1 1 0 0 0; ...
%                1 2 1 1 0; ...
%                2 2 0 0 1; ...
%                2 1 1 1 1]; % Trellis transitions of 2/3 RSC encoder
% transitions = [1 1 0 0 0; ...
%                1 2 1 1 1; ...
%                2 1 0 0 1; ...
%                2 2 1 1 0];
% transitions = [1 1 0 0 0; 
%                1 3 1 1 1; 
%                2 1 0 1 1; 
%                2 3 1 0 0; 
%                3 2 0 1 0; 
%                3 4 1 0 1; 
%                4 2 0 0 1; 
%                4 4 1 1 0];
transitions   = polynomial2trellis([[1 1] ; [5 7]]);
select_matrix = [[1 0]; [0 1]];
interleaver   = randperm(info_len);
max_iteration = 10;

execute_time = 0; % timer
for sim_index = 1:length(ebn0_array)
    
    tic; % start timer
    
    fprintf('Simulation round %d/%d\n', sim_index, length(ebn0_array));
    
    % SNR transform
    ebn0   = ebn0_array(1, sim_index); % get SNR in dB
    ebn0_1 = 10.0 .^ (ebn0/10.0); % eb/n0 in linear scale
    ebn0_2 = ebn0_1 * 0.5; % eb/n0 * coding rate
    Lc     = 2 * ebn0_2; % channel reliability measure
    
    err_total = 0;
    
    for block_index = 1:block_num
        
        if mod(block_index, 10) == 1
            fprintf('\n\t');
        end
        
        fprintf('%3d ', block_index);
        
        % Source side
        info_seq      = randi([0 1], [1 info_len]); % generate random information
        % [encoded_seq, parity_check_seq] = dacc(info_seq, doping_rate); % encode information
        encoded_seq   = turbo_encoder(info_seq, transitions, interleaver, select_matrix);
        modulated_seq = reshape(encoded_seq, 1, info_len * 2) * 2 - 1; % 0 -> -1, 1 -> +1
        
        % AWGN channel
        n0    = 1.0 / ebn0_2;
        sigma = sqrt(n0/2); % noise var.
        noise = sigma * randn(size(modulated_seq));
        
        % Transmit
        received_seq = modulated_seq + noise; % add awgn
        % received_seq = modulated_seq; % test
        
        % Destination side
        % cv_info_seq = zeros(1, info_len);
        % rcv_code_seq = zeros(1, info_len);
        % rcv_rstr_seq = received_seq;
        % Re-arrange infomation and parity check positions
%         rcv_info_seq(    1, 1          :doping_rate:end) = ...
%             received_seq(1, 1          :doping_rate:end);
%         rcv_code_seq(    1, doping_rate:doping_rate:end) = ...
%             received_seq(1, doping_rate:doping_rate:end);
%         rcv_rstr_seq = reshape([rcv_info_seq ; rcv_code_seq], ...
%             1, info_len * 2); % re-dopped and re-constructed received sequence
        
        % log-MAP Test
        % rcv_rstr_seq = reshape(conv_encoder(info_seq, transitions), 1, info_len * 2) * 2 - 1; 
        
        info_llr = turbo_decoder(received_seq, transitions, interleaver, select_matrix, max_iteration, Lc); % BCJR decode
        info_est = int32(sign(info_llr) == 1); % hard desicion
        
        % BER statistic
        err_total = err_total + sum(info_est ~= info_seq);
        
    end
    
    err_array(1, sim_index) = err_total;
    ber_array(1, sim_index) = double(err_total) / double(info_len * block_num);
    
    fprintf('\n\tSNR:%.2fdB ERR:%d BER:%e\n', ebn0, err_total, ber_array(1, sim_index));
    
    % Time estimation
    execute_time = execute_time + toc;
    hour = floor(execute_time / 3600);
    min  = floor(mod(execute_time, 3600) / 60);
    sec  = floor(mod(execute_time, 60));
    fprintf('\tElapsed time: %d hour %d min %d sec \n', hour, min, sec);
    
    est_time = (execute_time / sim_index) * (length(ebn0_array) - sim_index);
    hour = floor(est_time / 3600);
    min  = floor(mod(est_time, 3600) / 60);
    sec  = floor(mod(est_time, 60));
    fprintf('\tEstimated time of finish: %d hour %d min %d sec \n', hour, min, sec);
    
end

% SNR-BER plot
figure;
semilogy(ebn0_array, ber_array);
