function [ix,iy,it,iz] = findnzero(var,ilon,ilat,ilev,itime)
% set ilev to 1 if level dimension doesn't exist
ii=1;
for i = 1:numel(ilon)
   for j = 1:numel(ilat)
      for k = 1:numel(ilev)
         for t = 1:numel(itime)
            if (numel(ilev)>1)
               if (var(ilon(i),ilat(j),ilev(k),itime(t))~=0)
                  ix(ii) = ilon(i);
                  iy(ii) = ilat(j);
                  iz(ii) = ilev(k);
                  it(ii) = itime(t);
                  ii=ii+1;
               end
            else
               if (var(ilon(i),ilat(j),itime(t))~=0)
                  ix(ii) = ilon(i);
                  iy(ii) = ilat(j);
                  it(ii) = itime(t);
                  ii=ii+1;
               end
            end
         end
      end
   end
end
