% Doped accumulator simulation

close all;
clear all;

% Simulation parameters
info_len   = 100; % information length
block_num  = 100; % block numbers
ebn0_array = 0.0:0.2:10.0; % SNR array in dB
err_array_bcjr  = zeros(1, length(ebn0_array)); % Bit error statistic of BCJR
ber_array_bcjr  = double(zeros(1, length(ebn0_array))); % BER statistic of BCJR
err_array_dd    = zeros(1, length(ebn0_array)); % Bit error statistic of differential detection
ber_array_dd    = double(zeros(1, length(ebn0_array))); % BER statistic of differential detection

% Encoder attributes
doping_rate = 2; % Doping rate of DACC
% transitions = [1 1 0 0 0; ...
%                1 2 1 1 0; ...
%                2 2 0 0 1; ...
%                2 1 1 1 1]; % Trellis transitions of 2/3 RSC encoder
transitions = polynomial2trellis([[1 3]]);
initial_phase = +1;

execute_time = 0; % timer
for sim_index = 1:length(ebn0_array)
    
    tic; % start timer
    
    fprintf('\nSimulation round %d/%d\n', sim_index, length(ebn0_array));
    
    % SNR transform
    ebn0   = ebn0_array(1, sim_index); % get SNR in dB
    ebn0_1 = 10.0 .^ (ebn0/10.0); % eb/n0 in linear scale
    ebn0_2 = ebn0_1 * 1.0; % eb/n0 * coding rate
    Lc     = 2 * ebn0_2; % channel reliability measure
    
    err_total_bcjr = 0;
    err_total_dd   = 0;
    
    fprintf('\t');
    for block_index = 1:block_num
        
        if mod(block_index, 10) == 1 && block_index ~= 1
            fprintf('\n\t');
        end
        
        fprintf('%3d ', block_index);
        
        % Source side
        info_seq      = randi([0 1], [1 info_len]); % generate random information
        % [encoded_seq, parity_check_seq] = dacc(info_seq, doping_rate); % encode information
        encoded_seq   = reshape(conv_encoder(info_seq, transitions), 1, info_len * 1);
        modulated_seq = encoded_seq * 2 - 1; % 0 -> -1, 1 -> +1
        modulated_uncoded_seq = dbpsk_modulate(info_seq);
        
        % AWGN channel
        n0    = 1.0 / ebn0_2;
        sigma = sqrt(n0/2); % noise var.
        noise = sigma * randn(size(modulated_encoded_seq));
        
        % Transmit
        received_encoded_seq = modulated_encoded_seq + noise; % add awgn
        received_uncoded_seq = modulated_uncoded_seq + noise;  
        % received_seq = modulated_seq; % test
        
        % Destination side
%         rcv_info_seq = zeros(1, info_len);
%         rcv_code_seq = zeros(1, info_len);
%         % Re-arrange infomation and parity check positions
%         rcv_info_seq(    1, 1          :doping_rate:end) = ...
%             received_seq(1, 1          :doping_rate:end);
%         rcv_code_seq(    1, doping_rate:doping_rate:end) = ...
%             received_seq(1, doping_rate:doping_rate:end);
%         rcv_rstr_seq = reshape([rcv_info_seq ; rcv_code_seq], ...
%             1, info_len * 2); % re-dopped and re-constructed received sequence
        
        % log-MAP Test
        % rcv_rstr_seq = reshape([info_seq; parity_check_seq], 1, info_len * 2) * 2 - 1 + sigma * randn([1 info_len * 2]); 
        
        info_llr = rsc_decoder(received_encoded_seq, Lc, transitions, zeros(1, info_len));
        x_dd     = differential_detection(received_uncoded_seq);
        % info_llr = logmap_dacc(rcv_rstr_seq, zeros(1, info_len), Lc); % log-MAP algorithm decoded
        info_est = int32(sign(info_llr) == 1); % hard desicion
        info_dd  = int32(sign(x_dd) == 1);
        
        % BER statistic
        err_total_bcjr = err_total_bcjr + sum(info_est ~= info_seq);
        err_total_dd   = err_total_dd   + sum(info_dd  ~= info_seq);
        
    end
    
    err_array_bcjr(1, sim_index) = err_total_bcjr;
    ber_array_bcjr(1, sim_index) = double(err_total_bcjr) / double(info_len * block_num);
    err_array_dd(1, sim_index)   = err_total_dd;
    ber_array_dd(1, sim_index)   = double(err_total_dd) / double(info_len * block_num);
    
    fprintf('\n\tSNR:%.2fdB ERR:%d BER:%e\n', ebn0, err_total_bcjr, ber_array_bcjr(1, sim_index));
    
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
semilogy(ebn0_array, [ber_array_bcjr ; ber_array_dd]);
