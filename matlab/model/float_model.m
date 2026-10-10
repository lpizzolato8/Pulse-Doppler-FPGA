function o = float_model(p, L, cfg, n_cpi)
%FLOAT_MODEL Ideal double-precision reference of the processing chain.
%   It uses the same scenario, the same noise sequence and the same overall
%   scaling as golden_model, but no quantization, saturation or approximations.
%   o.map is in the same units as the fixed-point map: 10*log2|X|
%   (0.602 dB per unit).

NP = n_cpi * p.n_pulse;
nf = p.n_fast;
A  = 2^(p.coef_w-1) - 1;

[ni, nq] = emu_noise(p.noise_seed, nf * NP);
noise = reshape(ni + 1j*nq, nf, NP);

x = zeros(nf, NP);
for P = 0:NP-1
    s = zeros(nf, 1);
    for t = 1:p.n_targets
        tg = cfg.tgt(t);
        if ~tg.en, continue; end
        ph = exp(1j*2*pi*P*tg.dphi/2^p.phase_w) * tg.amp/2^p.amp_w;
        n  = tg.range + (0:p.n_chirp-1).';
        ok = n < nf;
        s(n(ok)+1) = s(n(ok)+1) + A*L.chirp_f(ok)*ph;
    end
    if cfg.clutter_en, s = s + L.clut_f; end
    if cfg.noise_en,   s = s + noise(:,P+1) * 2^cfg.noise_shift; end
    x(:, P+1) = s;
end
o.emu = x;

% Pulse compression with the ideal (optionally tapered) filter, same 2^-12 scaling
h = A * conj(flipud(L.chirp_f));
if cfg.pc_taper, h = h .* flipud(L.taper_f); end
idx = (0:p.n_range-1).' + p.n_chirp;
o.pc = zeros(p.n_range, NP);
for P = 1:NP
    y = conv(x(:,P), h);
    o.pc(:,P) = y(idx) / 2^p.pc_shift;
end

if cfg.mti_en
    o.mti = o.pc - [zeros(p.n_range,1), o.pc(:,1:end-1)];
else
    o.mti = o.pc;
end

o.fft = zeros(p.n_pulse, p.n_range, n_cpi);
o.map = zeros(p.n_pulse, p.n_range, n_cpi);
for c = 1:n_cpi
    X = o.mti(:, (c-1)*p.n_pulse + (1:p.n_pulse)).';
    if cfg.win_en, X = X .* L.win_f; end
    F = fft(X);                                  % natural order
    o.fft(:,:,c) = F;
    o.map(:,:,c) = fftshift(10*log2(max(abs(F), 1)), 1);
end
end
