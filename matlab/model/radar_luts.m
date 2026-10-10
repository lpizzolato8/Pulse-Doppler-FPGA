function L = radar_luts(p)
%RADAR_LUTS Build every ROM used by the RTL (quantized) plus float references.
%   All definitions follow docs/fixed_point_spec.md.

A = 2^(p.coef_w-1) - 1;                      % 2047

% LFM chirp, sweeping -B/2 .. +B/2
n = (0:p.n_chirp-1).';
t = n / p.fs;
L.chirp_f = exp(1j*pi*((p.B/p.tau)*t.^2 - p.B*t));
L.chirp_i = round(A*real(L.chirp_f));
L.chirp_q = round(A*imag(L.chirp_f));

% Matched filter h[k] = conj(chirp[N-1-k])
L.mf_i =  flipud(L.chirp_i);
L.mf_q = -flipud(L.chirp_q);

% Taylor-weighted (mismatched) filter: lower range sidelobes, ~0.9 dB SNR loss
L.taper_f = taylor_win(p.n_chirp, p.taper_nbar, p.taper_sll);
h = A * conj(flipud(L.chirp_f)) .* flipud(L.taper_f);
L.mft_i = round(real(h));
L.mft_q = round(imag(h));

% Doppler phasor LUT
k = (0:2^p.sincos_aw-1).';
L.cos = round(A*cos(2*pi*k/2^p.sincos_aw));
L.sin = round(A*sin(2*pi*k/2^p.sincos_aw));

% Stationary clutter return, identical every PRI (zero Doppler)
s = RandStream('mt19937ar', 'Seed', p.clutter_seed);
rb = (0:p.clutter_bins-1).';
refl = zeros(p.n_fast, 1);
refl(rb+1) = (randn(s, numel(rb), 1) + 1j*randn(s, numel(rb), 1)) .* 10.^(-rb/80);
x = conv(refl, L.chirp_f);
x = x(1:p.n_fast);
x = x / max(abs([real(x); imag(x)])) * p.clutter_peak;
L.clut_f = x;
L.clut_i = round(real(x));
L.clut_q = round(imag(x));

% Slow-time Blackman window, unsigned Q0.16. Its -58 dB sidelobes keep a
% target's Doppler sidelobes below the noise floor inside the MTI notch,
% where Hamming's -43 dB floor would show up as false detections.
pp = (0:p.n_pulse-1).';
L.win_f = 0.42 - 0.5*cos(2*pi*pp/(p.n_pulse-1)) + 0.08*cos(4*pi*pp/(p.n_pulse-1));
L.win   = abs(round((2^p.win_w - 1) * L.win_f));   % abs: w(0) evaluates to -1e-17

% FFT twiddles W_L^e = exp(-j 2 pi e / L), Q1.14
L.tw64_c = round(2^p.tw_frac * cos(2*pi*(0:63).'/64));
L.tw64_s = round(-2^p.tw_frac * sin(2*pi*(0:63).'/64));
L.tw16_c = round(2^p.tw_frac * cos(2*pi*(0:15).'/16));
L.tw16_s = round(-2^p.tw_frac * sin(2*pi*(0:15).'/16));

% log2 mantissa correction, 10 LSB per octave
L.log_lut = round(10*log2(1 + (0:31).'/32));
end
