function y = fx_rsr(x, s)
%FX_RSR Round-half-up arithmetic right shift: floor((x + 2^(s-1)) / 2^s).
%   x must hold exact integers (|x| < 2^53).
y = floor((x + 2^(s-1)) / 2^s);
end
