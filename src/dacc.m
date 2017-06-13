% Tested on Matlab 2013b
% Author: Weiwei Jiang
% Accumulator encoder. doping rate is p.

function [encoded_seq, parity_seq] = dacc(uncoded_bits, p)

    conv_encoded_seq = conv_encoder(uncoded_bits, polynomial2trellis([[1 1] ; [1 3]]));
    parity_seq = conv_encoded_seq(2:2:end);
    
    encoded_seq = uncoded_bits;
    for ii=p:p:length(uncoded_bits)
        encoded_seq(ii) = parity_seq(ii); % Dopping
    end
   
end