function o = golden_model(p, L, cfg, n_cpi)
%GOLDEN_MODEL Bit-accurate model of the full processing chain, from reset.
%   o.emu_i/q  [n_fast  x NP]          emulator output (NP = n_cpi * n_pulse)
%   o.pc_i/q   [n_range x NP]          pulse compression
%   o.mti_i/q  [n_range x NP]          MTI output (= pc when MTI is off)
%   o.win_i/q  [n_pulse x n_range x n_cpi]  windowed slow-time data (FFT input)
%   o.fft_i/q  [n_pulse x n_range x n_cpi]  FFT output, RTL output order
%   o.map      [n_pulse x n_range x n_cpi]  log magnitude, row = shifted Doppler bin
%   o.det_map, o.n_det(c)                   CFAR detections
%   o.pk_map, o.peaks{c}                    3x3 local-maximum peaks, [r d L] scan order

NP = n_cpi * p.n_pulse;
[o.emu_i, o.emu_q] = emu_model(p, L, cfg, NP);
[o.pc_i,  o.pc_q ] = pc_model(p, L, o.emu_i, o.emu_q, cfg.pc_taper);
[o.mti_i, o.mti_q] = mti_model(p, o.pc_i, o.pc_q, cfg.mti_en);

pos2d = mod(fx_bitrev((0:p.n_pulse-1).', 6) + p.n_pulse/2, p.n_pulse);  % position -> d

sz = [p.n_pulse, p.n_range, n_cpi];
[o.win_i, o.win_q, o.fft_i, o.fft_q, o.map] = deal(zeros(sz));
o.det_map = false(sz);
o.pk_map  = false(sz);
o.peaks = cell(n_cpi, 1);
o.n_det = zeros(n_cpi, 1);

for c = 1:n_cpi
    cols = (c-1)*p.n_pulse + (1:p.n_pulse);
    % Corner turn: [n_range x n_pulse] -> [n_pulse x n_range]
    xi = o.mti_i(:, cols).';
    xq = o.mti_q(:, cols).';
    if cfg.win_en
        xi = fx_rsr(xi .* L.win, p.win_w);
        xq = fx_rsr(xq .* L.win, p.win_w);
    end
    o.win_i(:,:,c) = xi;
    o.win_q(:,:,c) = xq;

    [fi, fq] = fft_r22_model(p, L, xi, xq);
    o.fft_i(:,:,c) = fi;
    o.fft_q(:,:,c) = fq;

    lg = mag_log_model(p, L, fi, fq);
    m = zeros(p.n_pulse, p.n_range);
    m(pos2d+1, :) = lg;
    o.map(:,:,c) = m;

    [o.det_map(:,:,c), o.n_det(c)] = cfar_model(p, m, cfg.cfar_off);
    [o.pk_map(:,:,c), o.peaks{c}] = peak_model(m, o.det_map(:,:,c));
end
end
