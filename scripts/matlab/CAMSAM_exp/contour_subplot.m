function contour_z_subplot(var,iilev,irow,cmap,clim)
   ill = 1;
   icol = ceil(numel(iilev)/irow);
   for il = iilev
      subplot(irow,icol,ill);
      contourf(var(:,:,il)',100,'linestyle','none');colorbar  
      caxis(clim)
      colormap(cmap);
      ill=ill+1;
   end
