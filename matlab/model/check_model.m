function ok = check_model(make_plots)
%CHECK_MODEL Self-checks of the golden model against theory and the float model.
%   Run this before generating vectors: if the model is wrong, bit-exact RTL is
%   wrong too.
if nargin < 1, make_plots = true; end
p = radar_params();
L = radar_luts(p);
S = scenarios(p);
cfgof = @(name) S(strcmp({S.name}, name)).cfg;
ok = true;

%% 1. Fixed-point FFT vs MATLAB fft (full-scale random input)
rs = RandStream('mt19937ar', 'Seed', 1);
x  = randi(rs, [-2^15, 2^15-1], 64, 500) + 1j*randi(rs, [-2^15, 2^15-1], 64, 500);
[fi, fq] = fft_r22_model(p, L, real(x), imag(x));
br = fx_bitrev((0:63).', 6);
Xf = zeros(size(fi));
Xf(br+1, :) = fi + 1j*fq;                      % back to natural order
err = abs(Xf - fft(x));
sqnr = 20*log10(norm(fft(x), 'fro') / norm(Xf - fft(x), 'fro'));
ok = report(ok, 'FFT matches fft() (SQNR > 90 dB)', sqnr > 90, ...
            sprintf('SQNR %.1f dB, max err %.1f LSB', sqnr, max(err(:))));

%% 2. Single target lands in the expected range/Doppler cell
cfg = cfgof('single_target');
o  = golden_model(p, L, cfg, cfg.n_cpi);
of = float_model(p, L, cfg, cfg.n_cpi);
t  = cfg.tgt(1);
d_exp = mod(round(t.dphi / 2^(p.phase_w - log2(p.n_pulse))) + p.n_pulse/2, p.n_pulse);
for c = 1:cfg.n_cpi
    [~, idx] = max(reshape(o.map(:,:,c), [], 1));
    [dpk, rpk] = ind2sub([p.n_pulse, p.n_range], idx);
    ok = report(ok, sprintf('CPI %d peak at (r=%d, d=%d)', c, t.range, d_exp), ...
        rpk-1 == t.range && dpk-1 == d_exp, sprintf('found (r=%d, d=%d)', rpk-1, dpk-1));
end

% Peak level vs float model (in 0.602 dB units)
pk_fx = max(o.map(:));  pk_fl = max(of.map(:));
ok = report(ok, 'Peak level within 1 dB of float', abs(pk_fx - pk_fl)*p.db_per_lsb < 1, ...
    sprintf('fixed %.1f dB, float %.1f dB', pk_fx*p.db_per_lsb, pk_fl*p.db_per_lsb));

% Pulse-compression gain: peak / (single-sample echo amplitude)
echo_amp = (2^(p.coef_w-1)-1) * t.amp / 2^p.amp_w;
pc_peak  = max(abs(of.pc(:))) * 2^p.pc_shift;
g = 20*log10(pc_peak / echo_amp) - 20*log10(2^(p.coef_w-1)-1);
ok = report(ok, 'Compression gain = 20log10(N) = 30.1 dB (coherent)', abs(g - 20*log10(p.n_chirp)) < 0.5, ...
    sprintf('%.2f dB', g));

% Range sidelobes, matched vs Taylor-weighted filter (noise-free autocorrelation)
psl = zeros(1, 2);
for tp = 0:1
    h = conj(flipud(L.chirp_f));
    if tp, h = h .* flipud(L.taper_f); end
    y = abs(conv(L.chirp_f, h));
    [ymax, k0] = max(y);
    k = k0;
    while k < numel(y) && y(k+1) < y(k), k = k + 1; end     % first null after the peak
    psl(tp+1) = 20*log10(max(y(k:end)) / ymax);
end
snr_loss = -10*log10(sum(L.taper_f)^2 / (p.n_chirp * sum(L.taper_f.^2)));
ok = report(ok, 'Taylor taper: PSL < -24 dB, SNR loss < 1.2 dB', psl(2) < -24 && snr_loss < 1.2, ...
    sprintf('PSL %.1f -> %.1f dB, loss %.2f dB', psl(1), psl(2), snr_loss));

%% 3. Map accuracy vs float for the multi-target + noise scenario
cfg = cfgof('multi_target');
o  = golden_model(p, L, cfg, cfg.n_cpi);
of = float_model(p, L, cfg, cfg.n_cpi);
m  = of.map > median(of.map(:)) + 17;          % >= 10 dB above the noise floor
e  = (o.map(m) - of.map(m)) * p.db_per_lsb;
ok = report(ok, 'Map error vs float, >10 dB above noise (|mean|<0.5, max<1.5 dB)', ...
    abs(mean(e)) < 0.5 && max(abs(e)) < 1.5, ...
    sprintf('mean %.2f dB, std %.2f dB, max %.2f dB over %d cells', mean(e), std(e), max(abs(e)), nnz(m)));
for t = find([cfg.tgt.en])
    tg = cfg.tgt(t);
    d  = mod(round(tg.dphi / 1024) + 32, 64);
    hit = any(any(o.det_map(max(d-1,0)+1:min(d+1,63)+1, tg.range+1, 1)));
    ok = report(ok, sprintf('Target %d (r=%d) detected by CFAR', t, tg.range), hit, '');
    pl = o.peaks{1};
    dd = min(mod(pl(:,2) - d, 64), mod(d - pl(:,2), 64));
    hit = any(abs(pl(:,1) - tg.range) <= 1 & dd <= 1);
    ok = report(ok, sprintf('Target %d (r=%d) has a peak within +/-1 cell', t, tg.range), hit, '');
end

%% 4. MTI clutter suppression at zero Doppler
off = golden_model(p, L, cfgof('clutter_mti_off'), 1);
on  = golden_model(p, L, cfgof('clutter_mti_on'),  1);
z = p.n_pulse/2 + 1;
cr = (mean(off.map(z, 1:p.clutter_bins)) - mean(on.map(z, 1:p.clutter_bins))) * p.db_per_lsb;
ok = report(ok, 'MTI suppresses zero-Doppler clutter by > 20 dB', cr > 20, sprintf('%.1f dB', cr));

%% 5. CFAR false-alarm rate on noise only
cfg = cfgof('noise_only');
cfg.n_cpi = 4;
o = golden_model(p, L, cfg, cfg.n_cpi);
pfa = sum(o.n_det) / numel(o.map);
ok = report(ok, 'Noise-only Pfa < 2e-4 at 15.7 dB offset', pfa < 2e-4, sprintf('Pfa = %.2e', pfa));

%% Plots
if make_plots
    cfg = cfgof('multi_target');
    o  = golden_model(p, L, cfg, 1);
    of = float_model(p, L, cfg, 1);
    rng_m = p.r_min/1e3 + (0:p.n_range-1)*p.range_bin/1e3;
    vel   = ((0:p.n_pulse-1) - p.n_pulse/2) * p.v_res;
    figure('Name', 'Range-Doppler map: fixed vs float', 'Position', [100 100 1200 450]);
    subplot(1,2,1);
    imagesc(rng_m, vel, of.map(:,:,1)*p.db_per_lsb); axis xy; colorbar;
    title('Float reference'); xlabel('Range [km]'); ylabel('Velocity [m/s]');
    subplot(1,2,2);
    imagesc(rng_m, vel, o.map(:,:,1)*p.db_per_lsb); axis xy; colorbar; hold on;
    [dd, rr] = find(o.det_map(:,:,1));
    plot(rng_m(rr), vel(dd), 'wo', 'MarkerSize', 8);
    title('Fixed-point golden model + CFAR'); xlabel('Range [km]'); ylabel('Velocity [m/s]');
end

if ok, fprintf('\nCHECK_MODEL: ALL PASSED\n');
else,  fprintf('\nCHECK_MODEL: FAILURES\n'); end
end

function ok = report(ok, name, pass, info)
if pass, tag = 'PASS'; else, tag = 'FAIL'; end
fprintf('[%s] %-55s %s\n', tag, name, info);
ok = ok && pass;
end
