%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Matlab 2013b
% Author: Weiwei Jiang (wjiang@jaist.ac.jp)
% Date: 20141218
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source-Relay-Destination simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

% Simulation parameters
info_len   = 10; % information length
block_num  = 10; % block numbers
ebn0_array = 0.0:1.0:10.0; % SNR array in dB
err_array  = zeros(1, length(ebn0_array)); % Bit error statistic
ber_array  = double(zeros(1, length(ebn0_array))); % BER statistic

% Encoder attributes
doping_rate = 2; % Doping rate of DACC
transitions_dacc = polynomial2trellis([[1 1] ; [1 3]]);
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
HI_max    = 10; % max turns of horizontal iterations
VI_max    = 10; % max turns of vertical iterations
round_max = 10; % max rounds of process
pe_threshold = 1; % threshold for calculating Pe estimation
transitions_src   = polynomial2trellis([[1 1] ; [5 7]]); % encoder of source node
transitions_relay = polynomial2trellis([[1 1] ; [5 7]]); % encoder of relay node
interleaver_src   = randperm(info_len * 2); % interleaver of source node
interleaver_relay = randperm(info_len * 2); % interleaver of relay node
interleaver_sr    = randperm(info_len); % interleaver of source-relay iteration
deinterleaver_src(interleaver_src)     = 1:1:info_len * 2; % de-interleaver of source node
deinterleaver_relay(interleaver_relay) = 1:1:info_len * 2; % de-interleaver of relay node
deinterleaver_sr(interleaver_sr)       = 1:1:info_len; % de-interleaver of source-relay iteration

execute_time = 0; % timer
for sim_index = 1:length(ebn0_array)
    
    tic; % start timer
    
    fprintf('\nSimulation round %d/%d\n', sim_index, length(ebn0_array));
    
    % Source-destination SNR transform
    ebn0_sd   = ebn0_array(1, sim_index); % get SNR in dB
    ebn0_1_sd = 10.0 .^ (ebn0_sd/10.0); % eb/n0 in linear scale
    ebn0_2_sd = ebn0_1_sd * 1.0; % eb/n0 * coding rate
    Lc_sd     = 2 * ebn0_2_sd; % source-destination channel reliability measure
    
    % Source-relay SNR transform
    ebn0_sr   = ebn0_array(1, sim_index) + 0.0; % get SNR in dB
    ebn0_1_sr = 10.0 .^ (ebn0_sr/10.0); % eb/n0 in linear scale
    ebn0_2_sr = ebn0_1_sr * 1.0; % eb/n0 * coding rate
    Lc_sr     = 2 * ebn0_2_sr; % source-destination channel reliability measure
    
    % Relay-destination SNR transform
    ebn0_rd   = ebn0_array(1, sim_index) + 0.0; % get SNR in dB
    ebn0_1_rd = 10.0 .^ (ebn0_rd/10.0); % eb/n0 in linear scale
    ebn0_2_rd = ebn0_1_rd * 1.0; % eb/n0 * coding rate
    Lc_rd     = 2 * ebn0_2_rd; % source-destination channel reliability measure
    
    err_total = 0;
    
    fprintf('\t');
    for block_index = 1:block_num
        
        if mod(block_index, 10) == 1 && block_index ~= 1
            fprintf('\n\t');
        end
        
        fprintf('%3d ', block_index);
        
% Start % Source node %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        info_seq_src      = randi([0 1], [1 info_len]); % generate random information
        encoded1_seq_src  = reshape(conv_encoder(info_seq_src, transitions_src), 1, 2 * info_len);
        % encoded2_seq_src  = reshape(conv_encoder(encoded1_seq_src(interleaver_src), transitions_dacc), 1, 4 * info_len);
        dopped_seq_src    = dacc(encoded1_seq_src(interleaver_src), doping_rate);
        % dopped_seq    = dacc(encoded_seq, doping_rate);
        % dopped_seq    = reshape(conv_encoder(encoded_seq(interleaver1), dacc_transitions), 1, 4 * info_len);
        modulated_seq_src = dopped_seq_src * 2 - 1; % 0 -> -1, 1 -> +1
        % modulated_seq_src = encoded2_seq_src * 2 - 1; % 0 -> -1, 1 -> +1
        
% End % Source node %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% Start % Channel % Slot 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % AWGN channel between source and destination
        n0_sd    = 1.0 / ebn0_2_sd;
        sigma_sd = sqrt(n0_sd/2); % noise var.
        noise_sd = sigma_sd * randn(size(modulated_seq_src));
        
        % AWGN channel between source and relay
        n0_sr    = 1.0 / ebn0_2_sr;
        sigma_sr = sqrt(n0_sr/2); % noise var.
        noise_sr = sigma_sr * randn(size(modulated_seq_src));
        
        received_seq_sd = modulated_seq_src + noise_sd; % source-destination awgn channel
        received_seq_sr = modulated_seq_src + noise_sr; % source-relay awgn channel
        % received_seq_sd = modulated_seq_src; % test
        % received_seq_sr = modulated_seq_src; % test
        % Lc = 100;

% End % Channel % Slot 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start % Relay node %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        rcv_info_seq_sr = zeros(1, info_len * 2); % RSC coded sequence
        rcv_code_seq_sr = zeros(1, info_len * 2); % ACC parity check sequence
        % Re-arrange infomation and parity check positions
        rcv_info_seq_sr(    1,           1:doping_rate:end) = ...
            received_seq_sr(1,           1:doping_rate:end);
        rcv_code_seq_sr(    1, doping_rate:doping_rate:end) = ...
            received_seq_sr(1, doping_rate:doping_rate:end);
        rcv_rstr_seq_sr = reshape([rcv_info_seq_sr ; rcv_code_seq_sr], ...
            1, info_len * 4); % re-dopped and re-constructed received sequence

        % rcv_rstr_seq_sr = received_seq_sr; % None-doping ACC test
        
        % No-iteration decode
        dacc_apllr_relay = zeros(1, info_len * 2); % La for DACC decoder
        dacc_llr_relay   = rsc_decoder(rcv_rstr_seq_sr, Lc_sr, transitions_dacc, dacc_apllr_relay); % decode DACC
        dacc_exllr_relay = dacc_llr_relay - Lc_sr * rcv_rstr_seq_sr(1:2:end) - dacc_apllr_relay; % Le
        dc_apllr_relay   = dacc_exllr_relay(deinterleaver_src); % La
        dc_apllr_relay   = dc_apllr_relay(1:2:end); % La for uncoded bits
        info_llr_relay   = bcjr_core(dc_apllr_relay, zeros(length(transitions_src(:, 1)), info_len), transitions_src); % decode Ds
        info_est_relay   = (sign(info_llr_relay) == 1); % hard decision
        
        % Relay node encode
        encoded1_seq_rd  = reshape(conv_encoder(info_est_relay(interleaver_sr), transitions_relay), 1, 2 * info_len); % relay encode
        % encoded2_seq_relay  = reshape(conv_encoder(encoded1_seq_relay(interleaver_relay), transitions_dacc), 1, 4 * info_len); % ACC encode
        dacc_seq_relay   = dacc(encoded1_seq_rd, doping_rate);
        modulated_seq_rd = dacc_seq_relay * 2 - 1; % 0 -> -1, 1 -> +1
        
% End % Relay node %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start % Channel % Slot 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % AWGN channel between relay and destination
        n0_rd    = 1.0 / ebn0_2_rd;
        sigma_rd = sqrt(n0_rd/2); % noise var.
        noise_rd = sigma_rd * randn(size(modulated_seq_rd));
        
        received_seq_rd = modulated_seq_relay + noise_rd; % relay-destination awgn channel
        % received_seq_rd = modulated_seq_rd; % test
        
% End %%% Channel % Slot 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start % Destination node %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        rcv_info_seq_sd = zeros(1, info_len * 2); % RSC coded sequence
        rcv_code_seq_sd = zeros(1, info_len * 2); % ACC parity check sequence
        % Re-arrange infomation and parity check positions
        rcv_info_seq_sd(    1,           1:doping_rate:end) = ...
            received_seq_sd(1,           1:doping_rate:end);
        rcv_code_seq_sd(    1, doping_rate:doping_rate:end) = ...
            received_seq_sd(1, doping_rate:doping_rate:end);
        rcv_rstr_seq_sd = reshape([rcv_info_seq_sd ; rcv_code_seq_sd], ...
            1, info_len * 4); % re-dopped and re-constructed received sequence
        
        rcv_info_seq_rd = zeros(1, info_len * 2); % RSC coded sequence
        rcv_code_seq_rd = zeros(1, info_len * 2); % ACC parity check sequence
        % Re-arrange infomation and parity check positions
        rcv_info_seq_rd(    1,           1:doping_rate:end) = ...
            received_seq_rd(1,           1:doping_rate:end);
        rcv_code_seq_rd(    1, doping_rate:doping_rate:end) = ...
            received_seq_rd(1, doping_rate:doping_rate:end);
        rcv_rstr_seq_rd = reshape([rcv_info_seq_rd ; rcv_code_seq_rd], ...
            1, info_len * 4); % re-dopped and re-constructed received sequence
        
        % rcv_rstr_seq_sd = received_seq_sd; % None-doping ACC test
        % rcv_rstr_seq_rd = received_seq_rd; % None-doping ACC test
        
        % Process
        for round_index = 1:round_max
        
            % Destination decoder
            % Source-Destination decoding process
            dacc_apllr_sd = zeros(1, info_len * 2); % La for DACC decoder
            for HI_index = 1:HI_max
                dacc_llr_sd   = rsc_decoder(rcv_rstr_seq_sd, Lc_sd, transitions_dacc, dacc_apllr_sd); % DACC decoder
                dacc_exllr_sd = dacc_llr_sd - Lc_sd * rcv_rstr_seq_sd(1:2:end) - dacc_apllr_sd; % calculate extrinsic LLR
                % dacc_exllr_dst = dacc_llr_dst - dacc_apllr_dst; % calculate extrinsic LLR
                dc_apllr_sd   = dacc_exllr_sd(deinterleaver_src); % interleave to La of decoder
                % dc_apllr_dst   = dacc_exllr_dst; % interleave to La of decoder

                dc_apllr_sd   = dc_apllr_sd(1:2:end); % extract information LLR for decoder
                info_llr_sd   = bcjr_core(dc_apllr_sd, zeros(length(transitions_src(:, 1)), info_len), transitions_src);
                % info_llr_dst   = bcjr_core_norm(dc_apllr_dst, ones(length(transitions1(:, 1)), info_len), transitions1);
                dc_exllr_sd   = info_llr_sd - dc_apllr_sd; % extrinsic LLR
                dacc_apllr_sd = zeros(1, info_len * 2); % add zeros for DACC parity check bits
                dacc_apllr_sd(1:2:end) = dc_exllr_sd; % La for DACC decoder
                dacc_apllr_sd = dacc_apllr_sd(interleaver_src); % interleave
            end
            
            % Relay-Destination decoding process
            dacc_apllr_rd = zeros(1, info_len * 2); % La for DACC decoder
            for HI_index = 1:HI_max
                dacc_llr_rd   = rsc_decoder(rcv_rstr_seq_rd, Lc_rd, transitions_dacc, dacc_apllr_rd); % DACC decoder
                dacc_exllr_rd = dacc_llr_rd - Lc_rd * rcv_rstr_seq_rd(1:2:end) - dacc_apllr_rd; % calculate extrinsic LLR
                % dacc_exllr_dst = dacc_llr_dst - dacc_apllr_dst; % calculate extrinsic LLR
                dc_apllr_rd   = dacc_exllr_rd(deinterleaver_relay); % interleave to La of decoder
                % dc_apllr_dst   = dacc_exllr_dst; % interleave to La of decoder

                dc_apllr_rd   = dc_apllr_rd(1:2:end); % extract information LLR for decoder
                info_llr_rd   = bcjr_core(dc_apllr_rd, zeros(length(transitions_relay(:, 1)), info_len), transitions_relay);
                % info_llr_dst   = bcjr_core_norm(dc_apllr_dst, ones(length(transitions1(:, 1)), info_len), transitions1);
                dc_exllr_rd   = info_llr_rd - dc_apllr_rd; % extrinsic LLR
                dacc_apllr_rd = zeros(1, info_len * 2); % add zeros for DACC parity check bits
                dacc_apllr_rd(1:2:end) = dc_exllr_rd; % La for DACC decoder
                dacc_apllr_rd = dacc_apllr_rd(interleaver_relay); % interleave
            end
            
            % Vertial interation
            for VI_index = 1:VI_max
                pe_hat    = pe_estimate(info_llr_sd(interleaver_sr), info_llr_rd, pe_threshold); % calculate Pe_hat
                if pe_hat == -1
                    break;
                end
                dcs_apllr = fc(dc_exllr_rd(deinterleaver_sr), pe_hat); % update LLR to La of Ds
                dcr_apllr = fc(dc_exllr_sd(interleaver_sr),   pe_hat); % update LLR to La of Dr
                
                % Decode process
                info_llr_sd = bcjr_core(dcs_apllr + dc_apllr_sd, zeros(length(transitions_src(:, 1)), info_len), transitions_src);
                info_llr_rd = bcjr_core(dcr_apllr + dc_apllr_rd, zeros(length(transitions_relay(:, 1)), info_len), transitions_relay);
                
                % Calculate Le
                dc_exllr_sd = info_llr_sd - dcs_apllr - dc_apllr_sd;
                dc_exllr_rd = info_llr_rd - dcr_apllr - dc_apllr_rd;
            end
            
            assert(0 == sum(isnan(info_llr_sd)));
            
        end % end of process
        
% End % Destination node %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Hard decision
        info_est_dst = (sign(info_llr_sd) == 1);
        
        % BER statistic
        err_total = err_total + sum(info_est_dst ~= info_seq_src);
        
    end
    
    err_array(1, sim_index) = err_total;
    ber_array(1, sim_index) = double(err_total) / double(info_len * block_num);
    
    fprintf('\n\tSNR:%.2fdB ERR:%d BER:%e\n', ebn0_sd, err_total, ber_array(1, sim_index));
    
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
