function Lg = mag_log_model(p, L, xi, xq)
%MAG_LOG_MODEL Alpha-max-beta-min magnitude + 10-LSB/octave log (spec section 6).
a  = abs(xi);  b = abs(xq);
mx = max(a, b);  mn = min(a, b);
mag = mx + floor(3*mn/8);
assert(all(mag(:) < 2^p.mag_w));

Lg = zeros(size(mag));
nz = mag > 0;
[~, e2] = log2(mag(nz));          % mag = f * 2^e2, f in [0.5,1) -> exact
e = e2 - 1;                       % index of the leading one
f = floor(mag(nz) * 32 ./ 2.^e) - 32;
Lg(nz) = 10*e + L.log_lut(f+1);
end
