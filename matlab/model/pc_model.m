function [yi, yq] = pc_model(p, L, xi, xq, taper)
%PC_MODEL Bit-accurate matched-filter pulse compression (spec section 2).
%   xi, xq : [n_fast x n_pulses] emulator output
%   taper  : 0 = matched filter, 1 = Taylor-weighted coefficients
%   yi, yq : [n_range x n_pulses], range bin r <- y[r + n_chirp - 1]
%   Every PRI is filtered on its own (the delay line is cleared at the PRI start).
%   Every partial sum is an integer < 2^53, so conv() is exact in double.

if taper
    hi = L.mft_i;  hq = L.mft_q;
else
    hi = L.mf_i;   hq = L.mf_q;
end
np = size(xi, 2);
yi = zeros(p.n_range, np);
yq = zeros(p.n_range, np);
idx = (0:p.n_range-1).' + p.n_chirp;    % 1-based index of y[r + n_chirp - 1]

for P = 1:np
    ai = conv(xi(:,P), hi) - conv(xq(:,P), hq);
    aq = conv(xi(:,P), hq) + conv(xq(:,P), hi);
    yi(:,P) = fx_sat(fx_rsr(ai(idx), p.pc_shift), p.pc_w);
    yq(:,P) = fx_sat(fx_rsr(aq(idx), p.pc_shift), p.pc_w);
end
end
