function [lags, c] = cxcorr(a, b, n)
%% Circularly correlate a and b across first dimension
    c = ifft(fft(a, n, 1).*conj(fft(b, n, 1)), n, 1);
    lags = 0:n-1;
end