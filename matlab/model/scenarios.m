function S = scenarios(p)
%SCENARIOS Test scenarios used for vector generation and model checks.
%   Each one returns a cfg struct: tgt(1..4).{en,range,dphi,amp}, clutter_en,
%   noise_en, noise_shift, mti_en, win_en, cfar_off, n_cpi, pc_taper.

S = struct('name', {}, 'cfg', {});

% 1) One target over low-level noise. Checks the range/Doppler mapping and scaling.
c = base(p);
c.noise_en = 1;
c.tgt(1) = tgt(p, 100, +50, 32);
S(end+1) = struct('name', 'single_target', 'cfg', c);

% 2) Four targets over noise, including two targets 2 range bins apart
%    (resolution) and both velocity signs.
c = base(p);
c.noise_en = 1;  c.noise_shift = 2;
c.tgt(1) = tgt(p,  40, +120, 24);
c.tgt(2) = tgt(p, 150,  -80, 16);
c.tgt(3) = tgt(p, 152,  -80, 16);
c.tgt(4) = tgt(p, 230,  +16,  8);
c.pc_taper = 1;
S(end+1) = struct('name', 'multi_target', 'cfg', c);

% 3/4) A slow mover in ground clutter, MTI off and on
c = base(p);
c.noise_en = 1;  c.noise_shift = 1;  c.clutter_en = 1;
c.tgt(1) = tgt(p, 60, +30, 16);
c.pc_taper = 1;
S(end+1) = struct('name', 'clutter_mti_off', 'cfg', c);
c.mti_en = 1;
S(end+1) = struct('name', 'clutter_mti_on', 'cfg', c);

% 5) Noise only, rectangular window (exercises CFAR false alarms)
c = base(p);
c.noise_en = 1;  c.noise_shift = 3;  c.win_en = 0;
S(end+1) = struct('name', 'noise_only', 'cfg', c);

% 6) Edge cases: range bins 0 and 255, velocity aliasing, ADC saturation
%    (two full-scale targets overlapping in range)
c = base(p);
c.noise_en = 1;
c.tgt(1) = tgt(p,   0, +300, 64);   % aliases to about -220 m/s
c.tgt(2) = tgt(p, 255, -255, 64);
c.tgt(3) = tgt(p, 120,  +40, 255);
c.tgt(4) = tgt(p, 130,  -40, 255);
c.pc_taper = 1;
S(end+1) = struct('name', 'edge_cases', 'cfg', c);

% 7) Long run: 8 CPIs (covers every address rotation of the in-place corner
%    turn), clutter + MTI + Taylor taper, four targets
c = base(p);
c.noise_en = 1;  c.noise_shift = 1;  c.clutter_en = 1;
c.mti_en = 1;  c.pc_taper = 1;  c.n_cpi = 8;
c.tgt(1) = tgt(p,  30, +200, 40);
c.tgt(2) = tgt(p,  90, -150, 24);
c.tgt(3) = tgt(p, 180,  +60, 32);
c.tgt(4) = tgt(p, 250,  -20, 12);
S(end+1) = struct('name', 'long_run', 'cfg', c);
end

function c = base(p)
c.tgt = repmat(struct('en', 0, 'range', 0, 'dphi', 0, 'amp', 0), 1, p.n_targets);
c.clutter_en  = 0;
c.noise_en    = 0;
c.noise_shift = 0;
c.mti_en      = 0;
c.win_en      = 1;
c.cfar_off    = 26;      % 15.7 dB above the local mean (Pfa < 1e-4)
c.n_cpi       = 2;
c.pc_taper    = 0;       % 1 = Taylor-weighted pulse compression
end

function t = tgt(p, range_bin, v, amp)
t = struct('en', 1, 'range', range_bin, 'dphi', vel2dphi(p, v), 'amp', amp);
end
