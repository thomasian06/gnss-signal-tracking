function [f, P] = singlefft(X, fs, truncate_even)
%% Function singlefft
%
% Purpose:
% To compute the single-sided FFT of X
%
% Inputs:
% X  - the input signal in the time domain
% fs - the sampling frequency of the signal
%
% Outputs:
% f  - a freqency array
% P  - the single-sided spectrum
%
% Author: Ian Thomas
% Last Modified: 11/6/2020

    L = length(X);
    if nargin > 2
        if truncate_even
            if mod(L, 2) == 1
                L = L - 1; 
                X = X(1:L);
            end
        end
    end

    Y = fft(X);
        
    PD = abs(Y/L);
    P = PD(1:floor(L/2) + 1);
    P(2:end-1) = 2*P(2:end-1);
    P = 10*log10(P);
    P(P == -Inf) = min(P(~isinf(P)));
    f = fs * (0:floor(L/2))/L;

end