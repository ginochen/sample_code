function [rhomap, pmap ] = corrcoef_time_map_nzero_par(var1,var2,ilat,ilon,ilev)
% var2 is the variable that contains a lot of zero values, and that
% spatial index should have correlation set to zero. The reason is the
% correlation at that area must be very weak that we can neglect.
if ndims(var2)==4 % with vertical dimension
   rhomap = zeros(numel(ilon),numel(ilat),numel(ilev));
else % 2D spatial map
   rhomap = zeros(numel(ilon),numel(ilat),1);
end
pmap = rhomap;
% assign longitude indices to chunks
n=3;
idx=ceil(linspace(1,numel(ilon),n));
for ii=1:n-2; ichunk{ii}=ilon(idx(ii):idx(ii+1)-1); end; ichunk{n-1}=ilon(idx(n-1):idx(n));
% parallelize the correlation according to longitude chunks
parfor ichk=1:numel(ichunk)
   for i=1:numel(ichunk{ichk});
      for j=1:numel(ilat);
         for k=1:numel(ilev)
            if ndims(var2) == 3
               v1 = squeeze(var1(ichunk{ichk}(i),ilat(j),:));
               v2 = squeeze(var2(ichunk{ichk}(i),ilat(j),:));
            else 
               v1 = squeeze(var1(ichunk{ichk}(i),ilat(j),ilev(k),:));
               v2 = squeeze(var2(ichunk{ichk}(i),ilat(j),ilev(k),:));
            end
            if (sum(v2~=0)>10) % if there is more than 100 nonzero samples, then evaluate the corr
               inzero = find(v2~=0);
               [rhomaptmp{ichk}(i,j,k), pmaptmp{ichk}(i,j,k)] = corr(v1(inzero),v2(inzero));
            else % too many zeros basically deemed the corr to be sufficiently weak
               rhomaptmp{ichk}(i,j,:)=NaN;
               pmaptmp{ichk}(i,j,:)=NaN;
            end
         end
      end
   end
end
% glue the chunks together and output
for i=1:numel(ichunk)
   nch{i} = numel(ichunk{i});
   rhomap(1+(i-1)*nch{i}:i*nch{i},:,:) = rhomaptmp{i};
   pmap(1+(i-1)*nch{i}:i*nch{i},:,:) = pmaptmp{i};
end
