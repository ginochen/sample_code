function varout = nzeromean=(var,idim);
sdim=size(var);
ndim=ndims(var);
ii=1;
if (idim == 4)
   for i = 1:sdim(1)
      for j = 1:sdim(2)
         for k = 1:sdim(3)
            inzero = find(squeeze(var(i,j,k,:))~=0);
            varout(i,j,k) = mean(var(i,j,k,inzero));
         end
      end
   end
elseif (idim == 3)
   for i = 1:sdim(1)
      for j = 1:sdim(2)
         for k = 1:sdim(3)
            inzero = find(squeeze(var(i,j,k,:))~=0);
            varout(i,j,k) = mean(var(i,j,k,inzero));
         end
      end
   end
   var1 = mean(var(i,j,inzero,k))
for nd=1:ndim
   for i = 1:sdim(nd)
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
