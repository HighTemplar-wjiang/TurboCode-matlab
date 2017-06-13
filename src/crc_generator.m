%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Matlab 2013b
% Author: Weiwei Jiang (wjiang@jaist.ac.jp)
% Date: 20150105
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CRC generator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% seq: sequence to be CRC-ed
% generator_ploynomial: generator polynomial array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% crc_code: CRC-ed sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function crc_code = crc_generator(seq, generator_polynomial)

    if generator_polynomial(1, end) ~= 1
        error('Generator polynomial must end with 1');
    end
    
    gp_len = length(generator_polynomial); % length of generator polynomial
    
    % Input right padded by 3 bits
    remainder = [seq zeros(1, gp_len - 1)];
    
    % Calculate remainder
    for iteration_index = 1:length(seq)
        
        if remainder(1) ~= 0
            remainder = [xor(remainder(1:gp_len), generator_polynomial) remainder(gp_len+1:end)];
        end
        
        remainder = remainder(2:end);
        
    end
    
    crc_code = [seq remainder];

end