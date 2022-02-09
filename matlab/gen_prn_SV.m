function [c, cG1, cG2] = gen_prn_SV(sv)
%% Function gps_CA
%
% Purpose:
% To generate the course-acquisition code for a GPS signal using an LFSR
%
% Inputs:
% sv - integer from 1 - 32
%
% Output:
% c - the C/A code for the sv
%
% Author: Ian Thomas
% Last Modified: 9/4/2021

    svs = [
      2 6
      3 7
      4 8
      5 9
      1 9
      2 10
      1 8
      2 9
      3 10
      2 3
      3 4
      5 6
      6 7
      7 8
      8 9
      9 10
      1 4
      2 5
      3 6
      4 7
      5 8
      6 9
      1 3
      4 6
      5 7
      6 8
      7 9
      8 10
      1 6
      2 7
      3 8
      4 9];
    
    n = 10;
    phase_selector = svs(sv, :);
    taps1 = [3 10];
    taps2 = [2 3 6 8 9 10];
    G1 = cast(1023, 'uint16');
    G2 = cast(1023, 'uint16');
    c = zeros(1, 2^n-1);
    cG1 = zeros(1, 2^n-1);
    cG2 = zeros(1, 2^n-1);
    for i = 1:1023
        o1 = bitget(G1, n);
        cG1(i) = o1;
        o2 = bitget(sum(bitget(G2, phase_selector)), 1);
        cG2(i) = o2;
        c(i) = xor(o1, o2);
        i1 = bitget(sum(bitget(G1, taps1)), 1);
        i2 = bitget(sum(bitget(G2, taps2)), 1);
        G1 = bitset(bitshift(G1, 1), 1, i1);
        G2 = bitset(bitshift(G2, 1), 1, i2);
    end
    c = logical(c(1:i));
    cG1 = logical(cG1(1:i));
    cG2 = logical(cG2(1:i));
end