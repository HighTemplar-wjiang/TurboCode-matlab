% file name:    logmap_dacc.m
% description:  [LLR,a,b,g,Le] = logmap(r,ap,Lc) performs logmap decoding for the example (1 5/7) rsc code. 
%               the inputs are r (received sequence), ap (a priori) and Lc (channel condition).
%               the outputs include LLR (llr), a (alpha), b (beta), g (gamma) and Le (extrinsic Lllr).
%               note: it is used in 'tbcdec.m'.
% algorithm:    logmap decoding algorithm
% author:       y. jiang 
% date:         june 2010
% revision:     1.0


function [LLR,a,b,g,Le] = logmap_dacc(r,ap,Lc)

    rs = r(1:2:end-1);                                          % extract systematic bits
    r = repmat(r,4,1);

    nn = length(rs); 

    if ap == 0
        ap = zeros(1,nn);
    end

    La = ap;
    Le = ap';

    transitions = [1 1 0 0 0; ...
                   1 2 1 1 0; ...
                   2 2 0 0 1; ...
                   2 1 1 1 1];
    c = transitions(:, 4:end) * 2 - 1;       % list of codewords in a trellis (pay attention to the order)

    g = zeros(4,length(ap));                                    % g: gamma    
    k = 1;
    for i=1:2:length(r(1, :))-1                                 % compute gamma. eq.(6.30)
        g(:,k) = 0.5*Lc*(c(:,1).*r(:,i) + c(:,2).*r(:,i+1)) + 0.5*c(:,1).*Le(k);
        k = k + 1;
    end

    a = -1e10*ones(2,length(ap)+1); a(1,1) = 0;                 % a: alpha. initialization
    % a = -Inf * ones(2, length(ap)+1); a(1,1) = 0;
    for i = 2:length(a(1, :))                                   % compute alpha. eq.(6.29) and (6.36)
        a1 = g(1,i-1) + a(1,i-1); a2 = g(4,i-1) + a(2,i-1);
        a(1,i) = max(a1,a2) + log(1+exp(-abs(a1-a2)));

        a1 = g(2,i-1) + a(1,i-1); a2 = g(3,i-1) + a(2,i-1);
        a(2,i) = max(a1,a2) + log(1+exp(-abs(a1-a2)));
    end

    b = -1e10*ones(2,length(ap)+1); b(1,end) = 0;               % b: beta. initialization
    for i = length(b(1, :))-1:-1:1                              % compute beta. eq.(6.29) and (6.36)
        b1 = g(1,i) + b(1,i+1); b2 = g(2,i) + b(2,i+1);
        b(1,i) = max(b1,b2) + log(1+exp(-abs(b1-b2)));

        b1 = g(3,i) + b(2,i+1); b2 = g(4,i) + b(1,i+1);
        b(2,i) = max(b1,b2) + log(1+exp(-abs(b1-b2)));
    end

    sp = zeros(1, nn);
    sm = zeros(1, nn);
    for i=1:nn                                                  % compute LLR. eq.(6.31) and (6.36)                                                   % 
        spa = a(1,i) + g(2,i) + b(2,i+1); spb = a(2,i) + g(4,i) + b(1,i+1); 
        sp(i) = max(spa,spb) + log(1+exp(-abs(spa-spb)));

        sma = a(1,i) + g(1,i) + b(1,i+1); smb = a(2,i) + g(3,i) + b(2,i+1);
        sm(i) = max(sma,smb) + log(1+exp(-abs(sma-smb)));
    end

    LLR = sp - sm;

    Le = LLR - Lc*rs - La;                                      % comupte Le.
end