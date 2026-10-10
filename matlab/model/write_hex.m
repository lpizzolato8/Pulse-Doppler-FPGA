function write_hex(fname, varargin)
%WRITE_HEX Write packed two's-complement words, one hex word per line.
%   write_hex(fname, x1, w1, x2, w2, ...) packs {x1[w1], x2[w2], ...} with x1
%   in the MSBs. Each x is a vector of integers (all the same length) and is
%   written in column-major order. The total width must be <= 53 bits.

fields = varargin(1:2:end);
widths = [varargin{2:2:end}];
W = sum(widths);
assert(W <= 53, 'packed width %d too large', W);

word = zeros(numel(fields{1}), 1);
for k = 1:numel(fields)
    v = fields{k}(:);
    w = widths(k);
    assert(all(v >= -2^(w-1) & v < 2^w), 'field %d out of range for %d bits', k, w);
    word = word * 2^w + mod(v, 2^w);
end

fid = fopen(fname, 'w');
assert(fid > 0, 'cannot open %s', fname);
fprintf(fid, sprintf('%%0%dX\\n', ceil(W/4)), uint64(word));
fclose(fid);
end
