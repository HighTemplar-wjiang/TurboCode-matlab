%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Matlab 2013b
% Author: Weiwei Jiang (wjiang@jaist.ac.jp)
% Date: 20141218
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MARC with beamforming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

% Simulation parameters
info_len   = 100; % information length
block_num  = 10; % block numbers
ebn0_array = 0.0:1.0:10.0; % SNR array in dB
err_array1 = zeros(1, length(ebn0_array)); % Bit error statistic, user 1
ber_array1 = double(zeros(1, length(ebn0_array))); % BER statistic, user 1
err_array2 = zeros(1, length(ebn0_array)); % Bit error statistic, user 2
ber_array2 = double(zeros(1, length(ebn0_array))); % BER statistic, user 2
frr_array1 = zeros(1, length(ebn0_array)); % Frame error statistic, user 1
fer_array1 = double(zeros(1, length(ebn0_array))); % FER statistic, user 1
frr_array2 = zeros(1, length(ebn0_array)); % Frame error statistic, user 2
fer_array2 = double(zeros(1, length(ebn0_array))); % FER statistic, user 2

% Encoder attributes
% doping_rate = 2; % Doping rate of DACC
transitions_acc   = polynomial2trellis([[1 3]]); % full accumulator
HI_max            = 2; % max turns of horizontal iterations
VI_max            = 5; % max turns of vertical iterations
round_max         = 10; % max rounds of process
pe_threshold      = 1; % threshold for calculating Pe estimation
transitions_src1  = polynomial2trellis([[1 1] ; [5 7]]); % encoder of source node 1
transitions_src2  = polynomial2trellis([[1 1] ; [5 7]]); % encoder of source node 2
transitions_relay = polynomial2trellis([[1 1] ; [5 7]]); % encoder of relay node
interleaver_src1  = randperm(info_len * 2); % interleaver of source node 1
interleaver_src2  = randperm(info_len * 2); % interleaver of source node 2
interleaver_xor   = randperm(info_len); % interleaver of relay node after xor
interleaver_relay = randperm(info_len * 2); % interleaver of relay node between encoders
interleaver_sr1   = randperm(info_len); % interleaver of source1-relay iteration
interleaver_sr2   = randperm(info_len); % interleaver of source2-relay iteration
deinterleaver_src1(interleaver_src1)   = 1:1:info_len * 2; % de-interleaver of source node 1
deinterleaver_src2(interleaver_src2)   = 1:1:info_len * 2; % de-interleaver of source node 2
interleaver_xor(interleaver_xor)       = 1:1:info_len; % de-interleaver of relay node after xor
deinterleaver_relay(interleaver_relay) = 1:1:info_len * 2; % de-interleaver of relay node between encoders
deinterleaver_sr1(interleaver_sr1)     = 1:1:info_len; % de-interleaver of source1-relay iteration
deinterleaver_sr2(interleaver_sr2)     = 1:1:info_len; % de-interleaver of source2-relay iteration

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
    
    err_total1 = 0;
    err_total2 = 0;
    frr_total1 = 0;
    frr_total2 = 0;
    
    fprintf('\t');
    for block_index = 1:block_num
        
        if mod(block_index, 10) == 1 && block_index ~= 1
            fprintf('\n\t');
        end
        
        fprintf('%3d ', block_index);
        
% Start % Source node %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Source 1
        info_seq_src1      = randi([0 1], [1 info_len]); % generate random information
        encoded1_seq_src1  = reshape(conv_encoder(info_seq_src1, transitions_src1), 1, 2 * info_len);
        encoded2_seq_src1  = reshape(conv_encoder(encoded1_seq_src1(interleaver_src1), transitions_acc), 1, 2 * info_len);
        modulated_seq_src1 = encoded2_seq_src1 * 2 - 1; % 0 -> -1, 1 -> +1
        
        % Source 2
        info_seq_src2      = randi([0 1], [1 info_len]); % generate random information
        encoded1_seq_src2  = reshape(conv_encoder(info_seq_src2, transitions_src2), 1, 2 * info_len);
        encoded2_seq_src2  = reshape(conv_encoder(encoded1_seq_src2(interleaver_src2), transitions_acc), 1, 2 * info_len);
        modulated_seq_src2 = encoded2_seq_src2 * 2 - 1; % 0 -> -1, 1 -> +1
        
% End % Source node %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% Start % Channel % Slot 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % AWGN channel between source and destination
        n0_sd    = 1.0 / ebn0_2_sd;
        sigma_sd1 = sqrt(n0_sd/2); % noise var.
        noise_sd1 = sigma_sd1 * randn(size(modulated_seq_src1));
        
        % AWGN channel between source and relay
        n0_sr    = 1.0 / ebn0_2_sr;
        sigma_sr1 = sqrt(n0_sr/2); % noise var.
        noise_sr1 = sigma_sr1 * randn(size(modulated_seq_src1));
        
        received_seq_sd1 = modulated_seq_src1 + noise_sd1; % source-destination awgn channel
        received_seq_sr1 = modulated_seq_src1 + noise_sr1; % source-relay awgn channel

% End %%% Channel % Slot 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start % Channel % Slot 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % AWGN channel between source and destination
        n0_sd    = 1.0 / ebn0_2_sd;
        sigma_sd2 = sqrt(n0_sd/2); % noise var.
        noise_sd2 = sigma_sd2 * randn(size(modulated_seq_src2));
        
        % AWGN channel between source and relay
        n0_sr    = 1.0 / ebn0_2_sr;
        sigma_sr2 = sqrt(n0_sr/2); % noise var.
        noise_sr2 = sigma_sr2 * randn(size(modulated_seq_src2));
        
        received_seq_sd2 = modulated_seq_src2 + noise_sd2; % source-destination awgn channel
        received_seq_sr2 = modulated_seq_src2 + noise_sr2; % source-relay awgn channel

% End %%% Channel % Slot 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start % Relay node %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        rcv_rstr_seq_sr1 = received_seq_sr1; % None-doping ACC
        rcv_rstr_seq_sr2 = received_seq_sr2; % None-doping ACC
        
        % No-iteration decode for user 1
        acc_apllr_relay1 = zeros(1, info_len * 2); % La for ACC decoder
        acc_llr_relay1   = rsc_decoder(rcv_rstr_seq_sr1, Lc_sr, transitions_acc, acc_apllr_relay1); % decode ACC
        acc_exllr_relay1 = acc_llr_relay1 - acc_apllr_relay1; % Le, ACC is non-systematic
        dc_apllr_relay1  = acc_exllr_relay1(deinterleaver_src1); % La
        dc_apllr_relay1  = dc_apllr_relay1(1:2:end); % La of uncoded bits
        info_llr_relay1  = bcjr_core(dc_apllr_relay1, zeros(length(transitions_src1(:, 1)), info_len), transitions_src1); % decode Ds1
        info_est_relay1  = (sign(info_llr_relay1) == 1); % hard decision
        
        % No-iteration decode for user 2
        acc_apllr_relay2 = zeros(1, info_len * 2); % La for ACC decoder
        acc_llr_relay2   = rsc_decoder(rcv_rstr_seq_sr2, Lc_sr, transitions_acc, acc_apllr_relay2); % decode ACC
        acc_exllr_relay2 = acc_llr_relay2 - acc_apllr_relay2; % Le, ACC is non-systematic
        dc_apllr_relay2  = acc_exllr_relay2(deinterleaver_src2); % La
        dc_apllr_relay2  = dc_apllr_relay2(1:2:end); % La of uncoded bits
        info_llr_relay2  = bcjr_core(dc_apllr_relay2, zeros(length(transitions_src2(:, 1)), info_len), transitions_src2); % decoder Ds2
        info_est_relay2  = (sign(info_llr_relay2) == 1); % hard decision
                
        % Relay node joint encode
        info_est_xor        = xor(info_est_relay1, info_est_relay2);
        encoded_seq_relay  = reshape(conv_encoder(info_est_xor(interleaver_sr1), transitions_relay), 1, 2 * info_len); % relay encode
        % encoded2_seq_relay  = reshape(conv_encoder(encoded1_seq_relay(interleaver_relay), transitions_acc), 1, 2 * info_len); % ACC encode
        modulated_seq_relay = encoded_seq_relay * 2 - 1; % 0 -> -1, 1 -> +1
        
% End %%% Relay node %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start % Channel % Slot 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % AWGN channel between relay and destination
        n0_rd    = 1.0 / ebn0_2_rd;
        sigma_rd = sqrt(n0_rd/2); % noise var.
        noise_rd = sigma_rd * randn(size(modulated_seq_relay));
        
        received_seq_rd = modulated_seq_relay + noise_rd; % relay-destination awgn channel
        % received_seq_rd = modulated_seq_relay; % test
        
% End %%% Channel % Slot 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start % Destination node %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Re-construct sequence
        % Full accumulator
        rcv_rstr_seq_sd1 = received_seq_sd1;
        rcv_rstr_seq_sd2 = received_seq_sd2;
        rcv_rstr_seq_rd  = received_seq_rd;
        
        % Process
        acc_apllr_sd1 = zeros(1, info_len * 2); % La for ACC decoder-user1
        acc_apllr_sd2 = zeros(1, info_len * 2); % La for ACC decoder-user2
        dacc_apllr_rd = zeros(1, info_len); % La for RSC decoder of relay-destination
        for round_index = 1:round_max
        
            % 1.Horizontal decoding
            for HI_index = 1:HI_max
                % Destination decoder
                % Source-Destination decoding process
                % User 1            
                acc_llr_sd1   = rsc_decoder(rcv_rstr_seq_sd1, Lc_sd, transitions_acc, acc_apllr_sd1); % ACC decoder
                acc_exllr_sd1 = acc_llr_sd1 - acc_apllr_sd1; % calculate extrinsic LLR
                dc_apllr_sd1  = acc_exllr_sd1(deinterleaver_src1); % interleave to La of decoder

                dc_apllr_sd1  = dc_apllr_sd1(1:2:end); % extract information LLR for RSC decoder
                info_llr_sd1  = bcjr_core(dc_apllr_sd1, zeros(length(transitions_src1(:, 1)), info_len), transitions_src1);
                % info_llr_dst  = bcjr_core_norm(dc_apllr_dst, ones(length(transitions1(:, 1)), info_len), transitions1);
                dc_exllr_sd1  = info_llr_sd1 - dc_apllr_sd1; % extrinsic LLR
                acc_apllr_sd1 = zeros(1, info_len * 2); % add zeros for ACC parity check bits
                acc_apllr_sd1(1:2:end) = dc_exllr_sd1; % La for ACC decoder
                acc_apllr_sd1 = acc_apllr_sd1(interleaver_src1); % interleave for ACC

                % User 2
                acc_llr_sd2   = rsc_decoder(rcv_rstr_seq_sd2, Lc_sd, transitions_acc, acc_apllr_sd2); % ACC decoder
                acc_exllr_sd2 = acc_llr_sd2 - acc_apllr_sd2; % Le of ACC
                dc_apllr_sd2  = acc_exllr_sd2(deinterleaver_src2); % La of RSC decoder

                dc_apllr_sd2  = dc_apllr_sd2(1:2:end); % Get information LLRs
                info_llr_sd2  = bcjr_core(dc_apllr_sd2, zeros(length(transitions_src2(:, 2)), info_len), transitions_src2);
                dc_exllr_sd2  = info_llr_sd2 - dc_apllr_sd2; % Le of RSC decoder
                acc_apllr_sd2 = zeros(1, info_len * 2); % La initialization of ACC
                acc_apllr_sd2(1:2:end) = dc_exllr_sd2; % Infomration LLR of La for ACC decoder
                acc_apllr_sd2 = acc_apllr_sd2(interleaver_src2); % interleave for ACC
            end
            
            % Error detecting (Not adopting yet)
            % info_est_dst1 = (sign(info_llr_sd1) == 1); % User 1 hard decision
            % info_est_dst2 = (sign(info_llr_sd2) == 1); % User 2 hard decision
            % err_user1 = sum(info_est_dst1 ~= info_seq_src1) ~= 0; % User 1 error flag
            % err_user2 = sum(info_est_dst2 ~= info_seq_src2) ~= 0; % User 2 error flag
            
            % Relay-Destination decoding process            
            dacc_llr_rd   = rsc_decoder(rcv_rstr_seq_rd, Lc_rd, transitions_relay, dacc_apllr_rd); % ACC decoder
            dacc_exllr_rd = dacc_llr_rd - dacc_apllr_rd; % calculate extrinsic LLR
            
            % Joint decoder process
            % TODO
            
            % Debug
            assert(0 == sum(isnan(info_llr_sd1)));
            assert(0 == sum(isnan(info_llr_sd2)));
            
        end % end of process
        
% End % Destination node %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Hard decision
        info_est_dst1 = (sign(info_llr_sd1) == 1);
        info_est_dst2 = (sign(info_llr_sd2) == 1);
        
        % BER statistic
        err_total1 = err_total1 + sum(info_est_dst1 ~= info_seq_src1);
        err_total2 = err_total2 + sum(info_est_dst2 ~= info_seq_src2);
        
        % FER statistic
        frr_total1 = frr_total1 + (sum(info_est_dst1 ~= info_seq_src1) ~= 0);
        frr_total2 = frr_total2 + (sum(info_est_dst2 ~= info_seq_src2) ~= 0);
        
    end
    
    % BER data
    err_array1(1, sim_index) = err_total1;
    ber_array1(1, sim_index) = double(err_total1) / double(info_len * block_num);
    err_array2(1, sim_index) = err_total2;
    ber_array2(1, sim_index) = double(err_total2) / double(info_len * block_num);
    
    % FER data
    frr_array1(1, sim_index) = frr_total1;
    fer_array1(1, sim_index) = double(frr_total1) / double(block_num);
    frr_array2(1, sim_index) = frr_total2;
    fer_array2(1, sim_index) = double(frr_total2) / double(block_num);
    
    fprintf('\n\tSNR1:%.2fdB ERR1:%d BER1:%e',   ebn0_sd, err_total1, ber_array1(1, sim_index));
    fprintf('\n\tSNR2:%.2fdB ERR2:%d BER2:%e\n', ebn0_sd, err_total2, ber_array2(1, sim_index));
    
    % Time estimation
    round_time   = toc;
    execute_time = execute_time + round_time;
    hour = floor(execute_time / 3600);
    min  = floor(mod(execute_time, 3600) / 60);
    sec  = floor(mod(execute_time, 60));
    fprintf('\tElapsed time: %d hour %d min %d sec \n', hour, min, sec);
    
    est_time = ((execute_time / sim_index) * 0.2 + round_time * 0.8) * (length(ebn0_array) - sim_index);
    hour = floor(est_time / 3600);
    min  = floor(mod(est_time, 3600) / 60);
    sec  = floor(mod(est_time, 60));
    fprintf('\tEstimated time of finish: %d hour %d min %d sec \n', hour, min, sec);
    
end

% SNR-BER plot
figure;
semilogy(ebn0_array, (ber_array1 + ber_array2) * 0.5);
