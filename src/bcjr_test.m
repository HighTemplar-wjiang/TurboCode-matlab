info_len = 6;

transitions = [1 1 0 0 0; 
               1 3 1 1 1; 
               2 1 0 1 1; 
               2 3 1 0 0; 
               3 2 0 1 0; 
               3 4 1 0 1; 
               4 2 0 0 1; 
               4 4 1 1 0];

ap_llr = zeros(1, info_len);

% Initialize gamma_c
rcv_info_seq = [0.3 0.1 -0.5 0.2 0.8 0.5 -0.5 0.3 0.1 -0.7 1.5 -0.4];
% rcv_info_seq = [0.3 0.1 -0.5 0.2 0.8 0.5 -0.5 0.3 0.1 -0.7 0.1 -0.4];
ys = zeros(2, info_len); % received infomation and parity check sequence

xk = (transitions(:, 4:end) * 2 - 1)'; % encoder ouput, 0 -> -1, 1 -> +1
ys(1, :) = rcv_info_seq(1, 1:2:end); % information sequence
ys(2, :) = rcv_info_seq(1, 2:2:end); % parity check sequence

Ck = 1;
Lc = 5.0;
gamma_c = zeros(length(transitions(:, 1)), info_len);
for time_index = 1:info_len
    for transition_index = 1:length(transitions(:,1))
        gamma_c(transition_index, time_index) = Ck * ...
            exp(Lc * 0.5 * sum(xk(:, transition_index) .* ys(:, time_index), 1));
    end 
end

llr = bcjr_core_norm(ap_llr, gamma_c, transitions);
