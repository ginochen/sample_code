function [rhomap, pmap ] = corrcoef_time_map(var1,var2,ilat,ilon,ilev)
% var2 is the variable that contains a lot of zero values, and that
% spatial index should have correlation set to zero. The reason is the
% correlation at that area must be very weak that we can neglect.
%
if ndims(var2)==4 % with vertical dimension
   rhomap = zeros(numel(ilon),numel(ilat),numel(ilev));
else % 2D spatial map
   rhomap = zeros(numel(ilon),numel(ilat),1);
end
pmap = rhomap;
for i=1:numel(ilon);
   for j=1:numel(ilat);
      for k=1:numel(ilev)
         if ndims(var2) == 3
            v1 = squeeze(var1(ilon(i),ilat(j),:));
            v2 = squeeze(var2(ilon(i),ilat(j),:));
         else 
            v1 = squeeze(var1(ilon(i),ilat(j),ilev(k),:));
            v2 = squeeze(var2(ilon(i),ilat(j),ilev(k),:));
         end
         [rhomap(i,j,k), pmap(i,j,k)] = corr(v1,v2);
      end
   end
end

