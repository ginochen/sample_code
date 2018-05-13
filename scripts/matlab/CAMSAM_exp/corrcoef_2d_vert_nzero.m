function [rho, p] = corrcoef_2d_vert_nzero(var1,var2,ilat,ilon,itimes,iilev)
% vertical pairwise corrcoef
itt=1;
for it = itimes
   ill=1;
   for il = iilev
      if ndims(var2)==4
         var21 = mv2(var2(ilon,ilat,il,it));
      elseif ndims(var2)==3
         var21 = mv2(var2(ilon,ilat,it));
      else
         error('var2 dimension inconsistent')
      end
      inzero = find(var21~=0);
      var11 = mv2(var1(ilon,ilat,il,it));
      [rho(ill,itt),p(ill,itt)]=corr(var1(inzero),var2(inzero));
      ill=ill+1;
   end
   itt=itt+1;
end
