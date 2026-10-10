function p = radar_params()
%RADAR_PARAMS System parameters. Mirrors rtl/pkg/radar_pkg.sv -- keep in sync.

% Clocking
p.f_clk        = 100e6;
p.f_pix        = 25e6;
p.clk_per_samp = 10;
p.fs           = p.f_clk / p.clk_per_samp;   % 10 MSPS

% Frame geometry
p.n_range = 256;                             % range bins after pulse compression
p.n_chirp = 32;                              % transmit pulse length [samples]
p.n_fast  = p.n_range + p.n_chirp;           % 288 samples per PRI
p.n_pulse = 64;                              % pulses per CPI

% Emulated RF scenario
p.c         = 299792458;
p.fc        = 10e9;
p.lambda    = p.c / p.fc;
p.B         = 5e6;                           % LFM sweep bandwidth
p.tau       = p.n_chirp / p.fs;              % 3.2 us
p.pri       = p.n_fast / p.fs;               % 28.8 us
p.prf       = 1 / p.pri;                     % 34.72 kHz
p.range_bin = p.c / (2 * p.fs);              % ~15 m
p.r_min     = p.c * p.tau / 2;               % ~480 m (range of bin 0)
p.v_unamb   = p.lambda * p.prf / 4;          % +/- ~260 m/s
p.v_res     = p.lambda * p.prf / (2 * p.n_pulse);  % ~8.13 m/s

% Fixed-point formats (see docs/fixed_point_spec.md)
p.adc_w     = 12;
p.coef_w    = 12;
p.pc_shift  = 12;
p.pc_w      = 16;
p.win_w     = 16;
p.fft_out_w = 23;
p.tw_frac   = 14;
p.mag_w     = 23;
p.logmag_w  = 8;
p.db_per_lsb = 20*log10(2)/10;               % 0.602 dB

% Target emulator
p.n_targets   = 4;
p.phase_w     = 16;
p.sincos_aw   = 10;
p.amp_w       = 8;
p.noise_seed  = uint32(hex2dec('ACE12468'));
p.clutter_bins = 96;                         % extent of the clutter region
p.clutter_peak = 1500;                       % peak clutter component [ADC LSB]
p.clutter_seed = 7;

% Pulse-compression taper (Taylor, used when cfg.pc_taper = 1)
p.taper_nbar = 4;
p.taper_sll  = -35;                          % design sidelobe level [dB]

% CFAR
p.cfar_train = 16;
p.cfar_guard = 2;
p.max_dets   = 32;
end
