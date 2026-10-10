function dphi = vel2dphi(p, v)
%VEL2DPHI Radial velocity [m/s] (positive = closing) -> 16-bit Doppler phase step.
%   One Doppler bin = 2^16 / n_pulse = 1024.
fd   = 2 * v / p.lambda;
dphi = mod(round(fd * p.pri * 2^p.phase_w), 2^p.phase_w);
end
