function [rho, p] = corrcoef_2d_vert(var1,var2,ilat,ilon,itimes,iilev)
% vertical pairwise corrcoef
itt=1;
for it = itimes
   ill=1;
   for il = iilev
      if ndims(var2)==4
         [rho(ill,itt),p(ill,itt)]=corr(mv2(var1(ilon,ilat,il,it)),mv2(var2(ilon,ilat,il,it)));
      elseif ndims(var2)==3
         [rho(ill,itt),p(ill,itt)]=corr(mv2(var1(ilon,ilat,il,it)),mv2(var2(ilon,ilat,it)));
      else
         error('var2 dimension inconsistent')
      end
      ill=ill+1;
   end
   itt=itt+1;
end
