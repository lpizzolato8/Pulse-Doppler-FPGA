function gen_vectors(names)
%GEN_VECTORS Run the golden model on every scenario and write testbench vectors
%   to sim/vectors/<scenario>/ (formats in docs/fixed_point_spec.md).
%   gen_vectors()                    % all scenarios
%   gen_vectors({'single_target'})   % a subset
p = radar_params();
L = radar_luts(p);
S = scenarios(p);
if nargin > 0
    S = S(ismember({S.name}, names));
end
root = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'sim', 'vectors');

for k = 1:numel(S)
    cfg = S(k).cfg;
    d = fullfile(root, S(k).name);
    if ~exist(d, 'dir'), mkdir(d); end
    tic;
    o = golden_model(p, L, cfg, cfg.n_cpi);

    write_cfg(fullfile(d, 'cfg.hex'), p, cfg);
    write_hex(fullfile(d, 'emu.hex'), o.emu_i, 12, o.emu_q, 12);
    write_hex(fullfile(d, 'pc.hex'),  o.pc_i,  16, o.pc_q,  16);
    write_hex(fullfile(d, 'mti.hex'), o.mti_i, 16, o.mti_q, 16);
    write_hex(fullfile(d, 'win.hex'), o.win_i, 16, o.win_q, 16);
    write_hex(fullfile(d, 'fft.hex'), o.fft_i, p.fft_out_w, o.fft_q, p.fft_out_w);
    % RD map word {peak, det, L8} at address d*256 + r
    write_hex(fullfile(d, 'map.hex'), double(permute(o.pk_map, [2 1 3])), 1, ...
              double(permute(o.det_map, [2 1 3])), 1, permute(o.map, [2 1 3]), 8);

    % Per CPI: n_det, n_peak, then the first max_dets peaks {r8, d6, L8}
    v = [];
    n_pk = zeros(1, cfg.n_cpi);
    for c = 1:cfg.n_cpi
        pl = o.peaks{c};
        n_pk(c) = size(pl, 1);
        pl = pl(1:min(end, p.max_dets), :);
        w  = zeros(p.max_dets, 1);
        w(1:size(pl,1)) = pl(:,1)*2^14 + pl(:,2)*2^8 + pl(:,3);
        v  = [v; min(o.n_det(c), 2^16-1); min(n_pk(c), 2^16-1); w]; %#ok<AGROW>
    end
    write_hex(fullfile(d, 'det.hex'), v, 22);

    fprintf('%-16s %5.1f s   CFAR hits per CPI: %s, peaks: %s\n', S(k).name, toc, ...
            mat2str(o.n_det.'), mat2str(n_pk));
end
end
