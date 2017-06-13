%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Matlab 2013b
% Author: Weiwei Jiang (wjiang@jaist.ac.jp)
% Date: 20150105
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CRC check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% crc_seq: sequence with CRC
% generator_polynomial: generator polynomial array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% crc_checksum: CRC checksum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function crc_checksum = crc_check(crc_seq, generator_polynomial)

    if generator_polynomial(1, end) ~= 1
        error('Generator polynomial must end with 1');
    end
    
    gp_len = length(generator_polynomial); % length of generator polynomial
    
    % Calculate remainder
    remainder = crc_seq;
    for iteration_index = 1:length(crc_seq) - gp_len + 1
        
        if remainder(1) ~= 0
            remainder = [xor(remainder(1:gp_len), generator_polynomial) remainder(gp_len+1:end)];
        end
        
        remainder = remainder(2:end);
        
    end
    
    crc_checksum = remainder;

end
