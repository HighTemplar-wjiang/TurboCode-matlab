% Created on Matlab 2013b
% Author: Weiwei Jiang
% Create Date: 20140922
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit signals through a 2x2 alamouti MIMO channel

function received = alamouti_2x2_stbc(Eb_No_dB, modulated_data)

    % Lengh of data
    len = length(modulated_data);

    % Initialize random
    hStr = RandStream('mt19937ar', 'Seed', 55408); % random stream
    
    % Create Rayleigh distributed channel response matrix
    h = zeros(4, len);
    h(:, 1:2:end) = (randn(hStr, 4, len/2) / sqrt(2)) + ...
                        (1i*randn(hStr, 4, len/2) / sqrt(2));
    h(:, 2:2:end) = h(:, 1:2:end); %[h11;h21;h12;h22]
    
    % Alamouti encoding
    tx_ala = zeros(2, len);
    tx_ala(:, 1:2:end) = reshape(modulated_data, 2, len/2) / sqrt(2);
    tx_ala(:, 2:2:end) = kron(ones(1, len/2), [-1;1]) .* ...
        flipud(conj(tx_ala(:, 1:2:end)));
    
    % Transmit signals
    rx1 = awgn(sum(h(1:2:end, :).*tx_ala, 1), Eb_No_dB, 0, hStr);
    rx2 = awgn(sum(h(2:2:end, :).*tx_ala, 1), Eb_No_dB, 0, hStr);

    % Form equalization matrix
    hMod = zeros(4, len);
    hMod(:, 1:2:end) = [conj(h(1, 1:2:end)); conj(h(2,   1:2:end)); ...
                             h(3, 1:2:end);       h(end, 1:2:end)];
    hMod(:, 2:2:end) = [conj(h(3, 2:2:end)); conj(h(end, 2:2:end)); ...
                            -h(1, 2:2:end);      -h(2, 2:2:end)];
    hPower = sum(h.*conj(h), 1);
    
    % Form y matrix
    yMod = zeros(4, len);
    yMod(:, 1:2:end) = [     rx1(1:2:end);       rx2(1:2:end); ...
                        conj(rx1(2:2:end)); conj(rx2(2:2:end))];
    yMod(:, 2:2:end) = yMod(:, 1:2:end);
    
    % Equalization
    xHat = sum(hMod .* yMod, 1) ./ hPower;
    
    % Output
    received = xHat;
    
end