function plot_mcs_pv(vo)
% US lon lat:
% lon = 235 to 300 (-125 to -60) 
% lat = 25 to 50
% plot the x-z contour of variables for a cluster 
load ~/comp
season = 'JJA';
%dlat = 'lat2050';
dlat = 'lat0055';
xpts = 50;
dpts = [num2str(xpts) 'pts'];
if strcmp(comp,'MAC')
  %cd /Users/g/archive/matlab/F_2000_SPCAM_m2005_3hrly1/atm/hist/lat2050/JJA/mcs_cluster_var/100pts
  cd(['/Users/g/archive/matlab/F_2000_SPCAM_m2005_3hrly1/atm/hist/' dlat '/' season '/mcs_cluster_var/' dpts])
  %cd /Users/g/archive/matlab/F_2000_SPCAM_m2005_3hrly1/atm/hist/lat2050/JJA/mcs_cluster_var/30pts
  load mcs_clusters.mat 
whos
  load ../../mcs_cluster_parm.mat
  load mcs_cluster_zonal_var.mat
else
  error(['plot on MAC, not on' comp])
end


%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ii=1;
for ic = 1:527; 
  for it = 1:nt4Cl(idxc(ic)); 
    iclon(it) = parm.lon(vo.mcsillt4Cl_fixlat{idxc(ic)}{it}(ceil(xpts/2),1)); 
    iclat(it) = parm.lat(vo.mcsillt4Cl_fixlat{idxc(ic)}{it}(1,2));
  end; 
  if any(iclon>=236 & iclon<=300 & iclat >=20 & iclat<=50 & island(idxc(ic))==0)
    ic_us(ii) = idxc(ic); 
    ii=ii+1;% US cluster indices
  end
  clear iclon iclat
end % all cluster central lon
disp(length(ic_us))
nt4Cl(ic_us)
for iic = 1:numel(ic_us)
ic = ic_us(iic);
uscl=1;
wscl=100;
%{
ic=527;
%}
contour(mean(cat(3,vo.thetaa{ic}{:}),3),[1:1:15],'k'); hold on; % lons over a fixed latitude
contour(mean(cat(3,vo.thetaa{ic}{:}),3),[-15:1:1],'k:'); hold on; % lons over a fixed latitude
contour(mean(cat(3,vo.pva{ic}{:}),3),[1e-6:1e-6:1e-5],'r'); hold on;
contour(mean(cat(3,vo.pva{ic}{:}),3),[-1e-5:1e-6:-1e-6],'r:'); hold on;
quiver(mean(cat(3,vo.u_stormrel{ic}{:}),3)*uscl,mean(cat(3,vo.w{ic}{:}),3)*wscl,'AutoScaleFactor',1); 
quiver(mean(cat(3,vo.v_stormrel{ic}{:}),3),mean(cat(3,vo.v_stormrel{ic}{:}),3),'AutoScaleFactor',1); hold off;
%quiver(mean(cat(3,vo.w{ic}{:}),3),mean(cat(3,vo.w{ic}{:}),3),'AutoScaleFactor',15); hold off;
pause
%%{
for it = 1:nt4Cl(ic)
  subplot(2,1,1)
  disp(['it, lat, lon = ' num2str(it) ', ' num2str(parm.lat(vo.mcsillt4Cl_fixlat{ic}{it}(1,2))) ', ' num2str(parm.lon(vo.mcsillt4Cl_fixlat{ic}{it}(25,1))'-360)])
  lons = parm.lon(vo.mcsillt4Cl_fixlat{ic}{it}(:,1))';
%%{
  contour(lons,parm.zint/1000,vo.thetaa{ic}{it},[0:1:15],'k'); hold on; % lons over a fixed latitude
  contour(lons,parm.zint/1000,vo.thetaa{ic}{it},[-15:1:0],'k:'); hold on; % lons over a fixed latitude
  contour(lons,parm.zint(1:end-1)/1000,vo.pva{ic}{it},[1e-6:1e-6:1e-5],'r'); hold on;
  contour(lons,parm.zint(1:end-1)/1000,vo.pva{ic}{it},[-1e-5:1e-6:1e-6],'r:'); hold on;
  quiver(lons,parm.zint/1000,vo.u_stormrel{ic}{it}*uscl,vo.w{ic}{it}*wscl,'AutoScaleFactor',1);
  quiver(lons,parm.zint/1000,vo.v_stormrel{ic}{it},vo.v_stormrel{ic}{it},'AutoScaleFactor',1); hold off
%  quiver(vo.w{ic}{it},vo.w{ic}{it},'AutoScaleFactor',5);
%  quiver(w{ic}{it}(:,ilons),vo.w{ic}{it}(:,ilons));
colormap(cmap(1)); 
%%}
  subplot(2,1,2)
  coast_centered(0); hold on
  scatter(parm.lon(vo.mcsillt4Cl_fixlat{ic}{it}(:,1)),parm.lat(vo.mcsillt4Cl_fixlat{ic}{it}(1,2))*ones(1,xpts),'b*'); 
  scatter(parm.lon(vo.mcsillt4Cl_fixlat{ic}{it}(ceil(xpts/2),1)),parm.lat(vo.mcsillt4Cl_fixlat{ic}{it}(1,2)),'r*'); hold off
  pause; 
end
%%}
end
