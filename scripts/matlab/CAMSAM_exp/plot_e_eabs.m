run ~/scripts/matlab/startup.m 
%load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/e_eabs_p1p30-ave.mat')
%load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/e_eabs_p14p30-ave.mat')
% load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/e_eabs_p24p30m.mat')
% load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/e_eabs_p14p24m.mat')
 load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/e_eabs_p21p30m.mat')
 load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/e_eabs_p14p21m.mat')
ie = 2
%ename={'eabs_p24p30m','eabs_p14p24m'};
ename={'eabs_p21p30m','eabs_p14p21m'};
etmp = eval(ename{ie});
figure;
colorbar;
colormap(cmap(0));
ilat = find(lat>=-35 & lat<=35);
ilon = 1:numel(lon);
lat_limits = [min(lat(ilat)) max(lat(ilat))];
lon_limits = [min(lon(ilon)) max(lon(ilon))];
m_proj('Equidistant Cylindrical','lon',lon_limits,'lat',lat_limits)
m_coast('color',[0 0 0],'linewidth',1.5);
m_grid('box','on','fontsize',10)
hold on
%m_contourf(lon(ilon),lat(ilat),mean(e_pm(ilon,ilat,:),3)',30,'linestyle','none')
% title('Q_{SPCAM} - Q_{CAM} (pressure-weighted-averaged, 14-days-averaged)')
%m_contourf(lon(ilon),lat(ilat),mean(e_p24p30m(ilon,ilat,:),3)',30,'linestyle','none')
m_contourf(lon(ilon),lat(ilat),mean(etmp(ilon,ilat,:),3)',30,'linestyle','none')
%m_contour(lon(ilon),lat(ilat),std(etmp(ilon,ilat,:),0,3)',[5,5],'r','linewidth',2)
m_contour(lon(ilon),lat(ilat),std(etmp(ilon,ilat,:),0,3)',[3,3],'r','linewidth',2)
title('|dQ|_{200-700mb}')
%title('|dQ|_{700-1000mb}')
% title('Q_{SPCAM} - Q_{CAM} (850-1000-pressure-weighted-averaged, 14-days-averaged)')
% title('|Q_{SPCAM} - Q_{CAM}|_{200-850mb} (14-days-averaged)')
%caxis([0 8])
set(gcf,'color',[1 1 1])
%savefig('/projects/rsmas/kirtman/gchen/archive/matlab/figure/e_eabs_stat/mean_std/p-ave/eabs_siggt5_p21p30m_t1t754m.fig')
savefig('/projects/rsmas/kirtman/gchen/archive/matlab/figure/e_eabs_stat/mean_std/p-ave/eabs_siggt5_p14p21m_t1t754m.fig')
% savefig('/projects/rsmas/kirtman/gchen/archive/matlab/figure/e_eabs_stat/mean/plev/e_t1t754-ave_p1p30-ave.fig')
