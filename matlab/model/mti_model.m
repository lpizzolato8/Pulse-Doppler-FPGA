function [yi, yq] = mti_model(p, xi, xq, mti_en)
%MTI_MODEL Two-pulse canceller (spec section 3). Pulse -1 is taken as zero.
if ~mti_en
    yi = xi;  yq = xq;
    return;
end
pi_ = [zeros(size(xi,1),1), xi(:,1:end-1)];
pq_ = [zeros(size(xq,1),1), xq(:,1:end-1)];
yi = fx_sat(xi - pi_, p.pc_w);
yq = fx_sat(xq - pq_, p.pc_w);
end
