load var_PC1_11.mat
lat=parm.lat;
lon=parm.lon;
nlon=numel(lon);
lon_limits=[0 360];
lat_limits=[-20 20];
m_proj('Equidistant Cylindrical','lon',lon_limits,'lat',lat_limits);
h = m_coast('color',[0 0 0],'linewidth',2.5);
m_grid('box','on','fontsize',20);

plotmap(lat,lon,20,0); % plot the coast in the latitude zone
X=sparse(44,nlon);  
X(1:22,97)=1; X(22:44,81)=1; X(22,81:97)=1; 
dlon = abs(lon-180);
ilon = find(dlon==min(dlon));
X(1:44,ilon)=1;
hold on;
m_contour(lon(:),lat(75:118),X,[1 1],'color',[.5 .5 .5],'linewidth',5)
%savefig('/Users/g/archive/matlab/figure/basin_map.fig')

