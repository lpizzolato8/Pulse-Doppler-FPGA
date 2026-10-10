function [det_map, n_det] = cfar_model(p, map, off)
%CFAR_MODEL Log-domain CA-CFAR along range for each Doppler row (spec section 7).
%   map     : [n_pulse x n_range] log magnitude, row = shifted Doppler bin d
%   off     : threshold offset [LSB, 0.602 dB each]
%   det_map : logical detections, same size as map
%   n_det   : total number of detections

G = p.cfar_guard;  T = p.cfar_train;  NR = p.n_range;
det_map = false(size(map));

for d = 0:p.n_pulse-1
    row = map(d+1, :);
    cs  = [0, cumsum(row)];                  % cs(k+1) = sum(row(1:k))
    wsum = @(lo, hi) cs(hi+2) - cs(lo+1);    % sum of row bins lo..hi (0-based)
    for r = 0:NR-1
        lead_ok = (r - G - T) >= 0;
        lag_ok  = (r + G + T) <= NR-1;
        if lead_ok && lag_ok
            S = wsum(r-G-T, r-G-1) + wsum(r+G+1, r+G+T);  sh = 5;
        elseif lead_ok
            S = wsum(r-G-T, r-G-1);  sh = 4;
        else
            S = wsum(r+G+1, r+G+T);  sh = 4;
        end
        det_map(d+1, r+1) = row(r+1) * 2^sh > S + off * 2^sh;
    end
end
n_det = nnz(det_map);
end
