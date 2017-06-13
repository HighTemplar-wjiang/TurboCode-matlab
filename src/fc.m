% Tested on Matlab 2013b
% Author: Xin He
% Modified by Weiwei Jiang
% LLR updating function

function LLR_new = fc(LLR, pe)

    LLR(abs(LLR) > 700) = sign(LLR(abs(LLR) > 700)) * 700;
    LLR_new = log(((1-pe).* exp(LLR)+ pe)./((1-pe)+pe.*exp(LLR)));

end