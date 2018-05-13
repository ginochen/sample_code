figure;
ilat = find(lat>=-35 & lat<=35);
ilon = 1:numel(lon);
icase = 3;
switch icase
   case 0
      title('Time-averaged P_{total}')
      vartmp = var{idv.prect}(ilon,ilat,:);
   case 1
      title('Time-averaged P_{Deep}')
      vartmp = var{idv.pzm}(ilon,ilat,:);
   case 2
      title('Time-averaged P_{Shallow}')
      vartmp = var{idv.psh}(ilon,ilat,:);
   case 3
      title('Time-averaged RH_{850-100mb}');
      vartmp = relhum_p24p30m(ilon,ilat,:);
end
lat_limits = [lat(ilat(1)) lat(ilat(end))];
lon_limits = [lon(ilon(1)) lon(ilon(end))];
m_proj('Equidistant Cylindrical','lon',lon_limits,'lat',lat_limits);
colorbar;
colormap(cmap(11));
m_coast('color',[0 0 0],'linewidth',1.5);
m_grid('box','on','fontsize',16)
hold on
m_contourf(lon(ilon),lat(ilat),mean(vartmp,3)',30,'linestyle','none'); 
caxis([0 8])
set(gca,'position',[0.05 0 0.9 1],'units','normalized');
set(gcf,'color',[1 1 1]);
set(gca,'OuterPosition',[-0.1010 -0.1350 1.1613 1.2270]);
% savefig('/projects/rsmas/kirtman/gchen/archive/matlab/figure/precip_rate/preccdzm_t-ave.fig')
