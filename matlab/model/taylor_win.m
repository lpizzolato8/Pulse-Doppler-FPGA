function w = taylor_win(n, nbar, sll_db)
%TAYLOR_WIN Taylor window, normalized to a peak of 1 (no toolbox needed).
%   n      : length
%   nbar   : number of nearly constant-level sidelobes
%   sll_db : design sidelobe level [dB], negative (e.g. -35)
A  = acosh(10^(-sll_db/20)) / pi;
s2 = nbar^2 / (A^2 + (nbar - 0.5)^2);
x  = ((0:n-1).' - (n-1)/2) / n;
w  = ones(n, 1);
for m = 1:nbar-1
    num = 1;
    den = 1;
    for k = 1:nbar-1
        num = num * (1 - m^2 / s2 / (A^2 + (k - 0.5)^2));
        if k ~= m
            den = den * (1 - m^2 / k^2);
        end
    end
    w = w + 2 * ((-1)^(m+1) * num / (2*den)) * cos(2*pi*m*x);
end
w = w / max(w);
end
