function vout = pwgtave(vin,ilev,dlev)
vout = zeros(size(vin,1),size(vin,2));
tt = sum(dlev(ilev)); % total atmosphere thickness
rtt = 1/tt;
for j = 1:numel(ilev)
   vout = 1/tt*dlev(ilev(j))*vin(:,:,ilev(j)) + vout;
end
