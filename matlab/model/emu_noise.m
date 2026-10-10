function [ni, nq] = emu_noise(seed, n)
%EMU_NOISE Emulator noise (spec section 1): n samples of approx. Gaussian I/Q.
%   One xorshift32 draw per sample: the current state is used, then advanced.
%   Each component is the sum of four 4-bit nibbles minus 30 (sigma ~ 9.2 LSB).
x  = uint32(seed);
ni = zeros(n, 1);
nq = zeros(n, 1);
m4 = uint32(15);
for k = 1:n
    si = 0; sq = 0;
    for b = 0:3
        si = si + double(bitand(bitshift(x, -4*b), m4));
        sq = sq + double(bitand(bitshift(x, -4*(b+4)), m4));
    end
    ni(k) = si - 30;
    nq(k) = sq - 30;
    x = bitxor(x, bitshift(x, 13));
    x = bitxor(x, bitshift(x, -17));
    x = bitxor(x, bitshift(x, 5));
end
end
