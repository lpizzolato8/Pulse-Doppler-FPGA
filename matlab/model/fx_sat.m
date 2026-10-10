function y = fx_sat(x, w)
%FX_SAT Saturate to the signed w-bit range.
y = min(max(x, -2^(w-1)), 2^(w-1) - 1);
end
