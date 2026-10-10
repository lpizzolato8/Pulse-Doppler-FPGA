function y = fx_bitrev(x, nbits)
%FX_BITREV Reverse the low nbits of each non-negative integer in x.
y = zeros(size(x));
for b = 0:nbits-1
    y = y + bitand(bitshift(x, -b), 1) * 2^(nbits-1-b);
end
end
