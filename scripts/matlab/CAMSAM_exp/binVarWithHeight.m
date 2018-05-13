function [pdf, cdf, bincenters] = binVarWithHeight(v,vrange,nbin,nz)
% v(1:nz,1:nsamp)
% bin the variable at each height
binedges = linspace(vrange(1),vrange(2),nbin);
for i=1:nbin-1
  bincenters(i) = mean(binedges(i:i+1));
end
for iz=1:nz
%  xx = histogram(v(iz,:),binedges); 
  pdf(:,iz) = histcounts(v(iz,:),binedges,'Normalization','probability'); 
  %npoint(1:nbin-1,iz) = xx.Values; 
  %prob(1:nbin-1,iz) = npoint(1:nbin-1,iz)/sum(npoint(1:nbin-1,iz));
  for ib = 1:nbin-1
    cdf(ib,iz) = sum(pdf(1:ib,iz));
  end
end
%prob = npoint'/sum(npoint(:));
