function write_cfg(fname, p, cfg)
%WRITE_CFG Scenario configuration for the testbenches, one 16-bit field per line:
%   for t = 1..n_targets: en, range, dphi, amp
%   then clutter_en, noise_en, noise_shift, mti_en, win_en, cfar_off, n_cpi, pc_taper
v = [];
for t = 1:p.n_targets
    tg = cfg.tgt(t);
    v = [v; tg.en; tg.range; tg.dphi; tg.amp]; %#ok<AGROW>
end
v = [v; cfg.clutter_en; cfg.noise_en; cfg.noise_shift; cfg.mti_en; ...
        cfg.win_en; cfg.cfar_off; cfg.n_cpi; cfg.pc_taper];
write_hex(fname, v, 16);
end
