% transitions = [1 1 0 0 0; 
%                1 3 1 1 1; 
%                2 3 0 0 0; 
%                2 1 1 1 1;
%                3 4 0 0 1;
%                3 2 1 1 0;
%                4 2 0 0 1;
%                4 4 1 1 0];
% Variables
info_len = 9; % Frame length
max_iteration = 20; % Turbo decoder iteration upper bound
Ec_No_dB = 0.25; % Channel SNR

% Generate RSC encoders
transitions = polynomial2trellis([[1 1] ; [5 7]]);
% transitions2 = polynomial2trellis([[1 1] ; [5 7]]);
% transitions_array(:, :, 1) = transitions;
% transitions_array(:, :, 2) = transitions;

% Generate information sequence
info_seq = [rand(1, info_len) > 0.5]; % add 2 tail bits
% info_seq = [1 0 1 0 1 0 1 0 0];
% info_seq = zeros(1, info_len + 2); % All-zero sequence
% info_seq = ones(1, info_len + 2); % All-one sequence

% Generate interleaver and de-interleaver
interleaver = randperm(info_len);
% interleaver = [1 4 7 2 5 9 3 6 8];

% Generate select matrix
select_matrix = [[1 0] ; [0 1]];
% select_matrix = [[1 1] ; [1 1]];
    
% Conduct channel coding
encoded_seq = turbo_encoder(info_seq, transitions, interleaver, select_matrix);

% Modulate
modulated_seq = encoded_seq * 2 - 1; % 0 -> -1, 1 -> +1

% Mapping and transmit through channels & demapping
% received_seq = real(alamouti_2x2_stbc(Ec_No_dB, modulated_seq));
% received_seq = awgn(modulated_seq, Ec_No_dB);
received_seq = modulated_seq;
% received_seq = [2 -5 6 1 2 -1 3 -1 2 -2 -2 -2 2 1 -5 -4 -2 5 -5 -1 -6];
% received_seq = [0.3 -4.0 -1.9 -2.0 -2.4 -1.3 1.2 -1.1 0.7 -2.0 -1.0 -2.1 -0.2 -1.4 -0.3 -0.1 -1.1 0.3];
% received_seq = [0.3 0.1 -0.5 0.2 0.8 0.5 -0.5 0.3 0.1 -0.7 1.5 -0.4];

% Demodulate % decode
app_llr = turbo_decoder(received_seq, transitions, interleaver, select_matrix, max_iteration, Ec_No_dB);
info_seq_est = (sign(app_llr) == 1);
eql = (info_seq_est == info_seq);

% y = [-0.8008 1.1992 -0.8008 1.1992 1.1992 -0.8008 1.1992 1.1992 -0.8008 1.1992];
% app_llr = zeros(1, length(y)/2);
% gamma_c = zeros(length(transitions(:, 1)), length(y)/2);
% 
% lc = 5.0;
% ck = 1.0;
% 
% xk = (transitions(:, 4:5) * 2 - 1)'; % 0 -> -1, 1 -> +1
% ys = reshape(y, 2, length(y) / 2);
% 
% for time_index = 1:length(y)/2
%     for transition_index = 1:length(transitions(:,1))
%         gamma_c(transition_index, time_index) = log(ck) + ...
%             (lc * 0.5 * sum(xk(:, transition_index) .* ys(:, time_index), 1));
%     end
% end