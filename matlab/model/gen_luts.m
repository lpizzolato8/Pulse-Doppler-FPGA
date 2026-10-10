function gen_luts()
%GEN_LUTS Write every RTL ROM initialization file (rtl/*/*.mem).
p = radar_params();
L = radar_luts(p);
root = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'rtl');
for sub = {'emu', 'dsp'}
    if ~exist(fullfile(root, sub{1}), 'dir'), mkdir(fullfile(root, sub{1})); end
end

write_hex(fullfile(root, 'emu', 'chirp_rom.mem'),   L.chirp_i, 12, L.chirp_q, 12);
write_hex(fullfile(root, 'emu', 'sincos_rom.mem'),  L.cos,     12, L.sin,     12);
write_hex(fullfile(root, 'emu', 'clutter_rom.mem'), L.clut_i,  12, L.clut_q,  12);
write_hex(fullfile(root, 'dsp', 'mf_coef.mem'),     [L.mf_i; L.mft_i], 12, [L.mf_q; L.mft_q], 12);
write_hex(fullfile(root, 'dsp', 'window_rom.mem'),  L.win,     16);
write_hex(fullfile(root, 'dsp', 'tw64_rom.mem'),    L.tw64_c,  16, L.tw64_s,  16);
write_hex(fullfile(root, 'dsp', 'tw16_rom.mem'),    L.tw16_c,  16, L.tw16_s,  16);
write_hex(fullfile(root, 'dsp', 'log_lut.mem'),     L.log_lut, 4);
fprintf('gen_luts: ROM files written under %s\n', root);
end
