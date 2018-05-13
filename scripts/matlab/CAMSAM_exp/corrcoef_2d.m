function [R, P, v1] = corrcoef_2d(var,ivarsall,ilat,ilon,ilev,dlev,itimes)
% 2D spatial cross-correlation coefficient matrix, with 3D variables
% pressure-weighted averaged to 2D before doing corrcoef
itt=1;
for it = itimes
   ivv = 1;
   v=zeros(numel(ilon),numel(ilat),numel(ivarsall));
   for iv = ivarsall
      if (ndims(var{iv})==4) % if iv is a 3D (spatially) variable which needed pressure-weighted vertical average
         v(:,:,ivv) = pwgtave(var{iv}(ilon,ilat,:,it),ilev,dlev);
      else
         v(:,:,ivv) = var{iv}(ilon,ilat,it); 
      end
      ivv=ivv+1;
   end
   for iv=1:numel(ivarsall)
      v1(:,iv,itt)=mv2(v(:,:,iv)); 
   end
   [R(:,:,itt), P(:,:,itt)] = corrcoef(v1(:,:,itt)); % make sure itimes is 1:700 no skipping indices
   itt=itt+1;
end
