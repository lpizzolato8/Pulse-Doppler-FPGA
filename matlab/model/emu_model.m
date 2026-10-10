function [xi, xq] = emu_model(p, L, cfg, n_pulses)
%EMU_MODEL Bit-accurate target/clutter/noise emulator (spec section 1).
%   Returns [n_fast x n_pulses] integer I and Q at the 12-bit "ADC" output,
%   starting from reset (pulse 0, xorshift seed).

nf = p.n_fast;
xi = zeros(nf, n_pulses);
xq = zeros(nf, n_pulses);

% Noise for every sample of every pulse, in time order. Generated first so
% the sequence is independent of the other settings.
[ni, nq] = emu_noise(p.noise_seed, nf * n_pulses);
ni = reshape(ni, nf, n_pulses);
nq = reshape(nq, nf, n_pulses);

ci = L.chirp_i;  cq = L.chirp_q;

for P = 0:n_pulses-1
    acc_i = zeros(nf, 1);
    acc_q = zeros(nf, 1);

    for t = 1:p.n_targets
        tg = cfg.tgt(t);
        if ~tg.en, continue; end
        phi = mod(P * tg.dphi, 2^p.phase_w);
        k   = floor(phi / 2^(p.phase_w - p.sincos_aw));
        pI  = floor(tg.amp * L.cos(k+1) / 2^p.amp_w);
        pQ  = floor(tg.amp * L.sin(k+1) / 2^p.amp_w);

        n  = tg.range + (0:p.n_chirp-1).';   % sample indices covered by the echo
        ok = n < nf;
        n  = n(ok);
        acc_i(n+1) = acc_i(n+1) + floor((ci(ok).*pI - cq(ok).*pQ) / 2^(p.coef_w-1));
        acc_q(n+1) = acc_q(n+1) + floor((ci(ok).*pQ + cq(ok).*pI) / 2^(p.coef_w-1));
    end

    if cfg.clutter_en
        acc_i = acc_i + L.clut_i;
        acc_q = acc_q + L.clut_q;
    end

    if cfg.noise_en
        acc_i = acc_i + ni(:, P+1) * 2^cfg.noise_shift;
        acc_q = acc_q + nq(:, P+1) * 2^cfg.noise_shift;
    end

    xi(:, P+1) = fx_sat(acc_i, p.adc_w);
    xq(:, P+1) = fx_sat(acc_q, p.adc_w);
end
end
