function [pk_map, peaks] = peak_model(map, det_map)
%PEAK_MODEL 3x3 local maxima among CFAR detections (spec section 8).
%   Doppler (rows) wraps around; range (columns) is padded with 0 at the edges.
%   A detection is a peak when it is > its 4 earlier neighbours (row d-1, and
%   cell r-1) and >= its 4 later neighbours (row d+1, and cell r+1), so exactly
%   one cell of an equal pair survives.
%   pk_map : logical, same size as map
%   peaks  : [k x 3] = [r, d, L] for every peak, scan order (d outer, r inner)

Lp = [zeros(size(map,1), 1), map, zeros(size(map,1), 1)];
up = circshift(Lp,  1, 1);                % row d holds row d-1
dn = circshift(Lp, -1, 1);                % row d holds row d+1
c  = Lp(:, 2:end-1);
earlier = max(max(up(:, 1:end-2), up(:, 2:end-1)), max(up(:, 3:end), Lp(:, 1:end-2)));
later   = max(max(dn(:, 1:end-2), dn(:, 2:end-1)), max(dn(:, 3:end), Lp(:, 3:end)));
pk_map  = det_map & (c > earlier) & (c >= later);

[r1, d1] = find(pk_map.');                % column-major over [r x d]: d outer, r inner
peaks = [r1-1, d1-1, map(sub2ind(size(map), d1, r1))];
peaks = reshape(peaks, [], 3);
end
