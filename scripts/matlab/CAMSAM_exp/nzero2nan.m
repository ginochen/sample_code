function varo = nzero2nan(var,ilon,ilat,ilev,itime)
% set ilev to 1 if level dimension doesn't exist
for i = 1:numel(ilon)
   for j = 1:numel(ilat)
      for k = 1:numel(ilev)
         for t = 1:numel(itime)
            if (numel(ilev)>1)
               if (var(ilon(i),ilat(j),ilev(k),itime(t))==0)
                  varo(ilon(i),ilat(j),ilev(k),itime(t))=NaN;
               else
                  varo(ilon(i),ilat(j),ilev(k),itime(t))=var(ilon(i),ilat(j),ilev(k),itime(t));
               end
            else
               if (var(ilon(i),ilat(j),itime(t))==0)
                  varo(ilon(i),ilat(j),itime(t))=NaN;
               else
                  varo(ilon(i),ilat(j),itime(t))=var(ilon(i),ilat(j),itime(t));
               end
            end
         end
      end
   end
end
