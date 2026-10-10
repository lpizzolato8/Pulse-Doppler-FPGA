function [yi, yq] = fft_r22_model(p, L, xi, xq)
%FFT_R22_MODEL Bit-accurate 64-point radix-2^2 DIF FFT (spec section 5).
%   xi, xq : [64 x M] inputs in natural order (one column per transform)
%   yi, yq : [64 x M] outputs in RTL output order: position p holds X[bitrev6(p)]

N = 64;
assert(size(xi,1) == N);
yi = xi;  yq = xq;
Lspan = N;
while Lspan >= 4
    for b = 0:Lspan:N-1
        % BF2I
        a  = b + (0:Lspan/2-1);
        bb = a + Lspan/2;
        [yi, yq] = bfly(yi, yq, a, bb);
        % BF2II (second half: multiply lower input by -j first)
        for h = 0:1
            a  = b + h*Lspan/2 + (0:Lspan/4-1);
            bb = a + Lspan/4;
            if h == 1
                t = yi(bb+1,:);
                yi(bb+1,:) = yq(bb+1,:);
                yq(bb+1,:) = -t;
            end
            [yi, yq] = bfly(yi, yq, a, bb);
        end
        % Twiddles
        if Lspan > 4
            if Lspan == 64
                tc = L.tw64_c;  ts = L.tw64_s;
            else
                tc = L.tw16_c;  ts = L.tw16_s;
            end
            for q = 0:3
                for n3 = 0:Lspan/4-1
                    e = n3 * fx_bitrev(q, 2);
                    if e == 0, continue; end
                    m = b + q*Lspan/4 + n3 + 1;
                    c = tc(e+1);  s = ts(e+1);
                    r = yi(m,:);  i = yq(m,:);
                    yi(m,:) = fx_rsr(r*c - i*s, p.tw_frac);
                    yq(m,:) = fx_rsr(r*s + i*c, p.tw_frac);
                end
            end
        end
    end
    Lspan = Lspan / 4;
end

lim = 2^(p.fft_out_w-1);
assert(all(abs(yi(:)) < lim & abs(yq(:)) < lim), 'FFT output exceeds %d bits', p.fft_out_w);
end

function [xi, xq] = bfly(xi, xq, a, b)
ai = xi(a+1,:);  aq = xq(a+1,:);
bi = xi(b+1,:);  bq = xq(b+1,:);
xi(a+1,:) = ai + bi;  xq(a+1,:) = aq + bq;
xi(b+1,:) = ai - bi;  xq(b+1,:) = aq - bq;
end
