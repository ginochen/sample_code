function [lonc, latc] = plot_mcs_xy(iic,it)
% ic=20; for it=1:40; [lonc,latc] = plot_mcs_xy(ic,it); if it==1; l1=lonc([1,end]);l2=latc([1,end]);    end; xlim([l1(1)-20 l1(2)+20]); ylim([l2(1)-15 l2(2)+15]); pause; end % 11 19 20
% 20: over US, but with the weird merid convergence and also stops at 25 deg which is unreal
% 21: over China, very zonally elongated
% 22: ocean east of Taiwan, weird merid conv
% 27: Indian ocean, southward moving, convergence zone
% 28: good org
% 12: zonal shear is parallel to the propagating dir
% 9,10,11,14: zonal shear is perpendicular to the propagating dir, and lead
% 15: zonal shear is perpendicular to the propagating dir, and lag
% 8: zonal shear all over the place, but on ave it's behind the propgating dir 
% US lon lat:
% lon = 235 to 300 (-125 to -60) 
% lat = 25 to 50
% plot the x-z contour of variables for a cluster 
docorr = 0;
load ~/comp
season = 'JJA';
%dlat = 'lat2050';
dlat = 'lat2525';
xpts = 30;
dpts = [num2str(xpts) 'pts'];
dims = '2d';
casei = 'F_2000_SPCAM_m2005_3hrly2';
if strcmp(comp,'MAC')
  %cd /Users/g/archive/matlab/F_2000_SPCAM_m2005_3hrly1/atm/hist/lat2050/JJA/mcs_cluster_var/100pts
  cd(['/Users/g/archive/matlab/' casei '/atm/hist/' dlat '/' season '/mcs_cluster_var/' dpts '/' dims])
  diro = ['/Users/g/archive/matlab/' casei '/figure/' dlat '/' season '/mcsmap/w/'] % output dir for figures
  %cd /Users/g/archive/matlab/F_2000_SPCAM_m2005_3hrly1/atm/hist/lat2050/JJA/mcs_cluster_var/30pts
  load ../../../mcs_clusters.mat 
  load ../../../mcs_cluster_parm.mat
  load ../../mcs_cluster_u_storm.mat
  load ../../mcs_cluster_v_storm.mat
  load ../../mcs_cluster_u_cam.mat
  load ../../mcs_cluster_v_cam.mat
  eval(sprintf(['load mcs_cluster_xy_var_' sprintf('%02d',iic) '.mat']))
else
  error(['plot on MAC, not on' comp])
end

uscl=1;
wscl=100;
ic=idxc(iic)
disp(['cluster # = ' num2str(ic)])
npanel = 2;
nprow = 1;
yl = 1:xpts;
xl = 1:xpts;
ext = 14;
%ixc = xpts/2-5:xpts/2+5;
%iyc = xpts/2-5:xpts/2+5;
ixc = xpts/2-ext:xpts/2+ext;
iyc = xpts/2-ext:xpts/2+ext;
disp(['ntime=' num2str(Nt(iic))]);
%for it=1:1:Nt(iic)
  lonc = parm.lon(vo.ilons{ic}{it}(ixc));
  latc = parm.lat(vo.ilat{ic}{it}(iyc));
disp(['it= ' num2str(it)])
disp(['time='  t{mcsillt4Cl{ic}{it}(1,3)}]); % time=20 the wake low appears with high entropy drawn from subsidence; two kinds of subsidence, rain (low entropy) & dry air (high entropy) 
disp(['lon lat=' num2str(parm.lon(mcsilltcentroids{ic}(it,1))) ', ' num2str(parm.lat(mcsilltcentroids{ic}(it,2)))]); % time=20 the wake low appears with high entropy drawn from subsidence; two kinds of subsidence, rain (low entropy) & dry air (high entropy) 

%%%%%%%%%%%%%%%%%%%%%% CORR (spatial) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if docorr
%disp(['du vs spdt3km ' num2str(corr(vo.dU03km_cam{ic}{it}(:),reshape(vo.spdt{ic}{it}(:,:,parm.i3km),[],1),'rows','complete')) ]) % spdt at 3km is similar to precc, so corr is close
%disp(['du vs precc ' num2str(corr(reshape(vo.dU03km_cam{ic}{it}(ixc,iyc),[],1),reshape(vo.precc{ic}{it}(ixc,iyc),[],1),'rows','complete')) ]) % mcs domain corr not obvious
disp(['du vs CPI ' num2str(corr(vo.dU03km_cam{ic}{it}(:),vo.Bmax{ic}{it}(:),'rows','complete')) ])
disp(['du vs pc ' num2str(corr(vo.dU03km_cam{ic}{it}(:),vo.precc{ic}{it}(:),'rows','complete')) ]); % 0.2 to 0.5 corr on the large-scale domain (3000km * 3000km) for iic=3, but low for iic=1
disp(['CPI vs pc ' num2str(corr(vo.Bmax{ic}{it}(:),vo.precc{ic}{it+1}(:),'rows','complete'))])
disp(['pc vs CPI ' num2str(corr(vo.precc{ic}{it}(:),vo.Bmax{ic}{it+1}(:),'rows','complete'))])
%disp(['CPI vs pl ' num2str(corr(vo.Bmax{ic}{it}(:),vo.precl{ic}{it}(:),'rows','complete'))])
disp(['pc vs LI ' num2str(corr(vo.precc{ic}{it}(:),vo.LIave{ic}{it}(:),'rows','complete'))]) % no corr
end
%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
ii=1;
for ic = 1%:527; 
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
%}
%for iic = 2%:numel(ic_us)
%{
ic = ic_us(iic);
%}
%%%%%%%%%%%%%%%%%%%%%
  ip = 1;
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  title('ps')
%  [C,h]= contourf(lonc,latc,vo.ps{ic}{it}(xl,yl)'/100-mean(mean(vo.ps{ic}{it}(ixc,iyc)/100)),'linestyle','none'); hold on; %-mean(vo.ps{ic}{it}(:)));
  [C,h]= contourf(lonc,latc,vo.ps{ic}{it}(ixc,iyc)'/100-parm.p0/100,'linestyle','none'); hold on; %-mean(vo.ps{ic}{it}(:)));
  scatter(parm.lon(mcsillt4Cl{ic}{it}(:,1)),parm.lat(mcsillt4Cl{ic}{it}(:,2))); 
  coast_centered(0)
  xlim([lonc(1) lonc(end)]);
  ylim([latc(1) latc(end)]);
%  colorbar; caxis([990 1020])
  caxis([-15 15]); colorbar; colormap(cmap(24,0,30))
  h.LevelStep=0.5; % contour spaced at 5hPa
  clabel(C,h);
  axis square
  hold off
  ip = ip+1;
%{
%%%%%%%%%%%%%%%%%%%%%
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  contourf(lonc,latc,vo.dU03km_cam{ic}{it}(ixc,iyc)',150,'linestyle','none');
%load mcs_cluster_xy_var_28.mat; clear cc;c = cellfun(@mean,vo.dU03km_cam{end},'Un',0);for ii=    8.2245); cc(ii) = c{ii}(15); hold on; end; disp(max(cc))
  colorbar; caxis([0 15]); colormap(hh(ip),cmap(2)); hold on
  scatter(parm.lon(mcsillt4Cl{ic}{it}(:,1)),parm.lat(mcsillt4Cl{ic}{it}(:,2))); 
  title('0-3shear')
  coast_centered(0)
  xlim([lonc(1) lonc(end)]);
  ylim([latc(1) latc(end)]);
  axis square
  hold off
  ip=ip+1;
%%%%%%%%%%%%%%%%%%%%%
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  % current timestep
  contourf(lonc,latc,vo.Bave{ic}{it}(ixc,iyc)',300,'linestyle','none'); hold on
  colorbar; caxis([0 80]); colormap(hh(ip),cmap(2)); hold on
%  contour(lonc,latc,vo.Bmax{ic}{it}(ixc,iyc)',[300 300],'m'); hold on
%  colorbar; caxis([0 300]); colormap(hh(ip),cmap(2)); hold on
%  contour(lonc,latc,vo.dU03km_cam{ic}{it}(ixc,iyc)',[8 8],'r'); hold on
%  contour(lonc,latc,abs(diff(vo.u_cam{ic}{it}(ixc,iyc,[1 parm.i3km]),1,3)'),[3 3],'r'); hold on
  % next timestep
%  lonc1 = parm.lon(vo.ilons{ic}{it+1}(ixc));
%  latc1 = parm.lat(vo.ilat{ic}{it+1}(iyc));
%  contour(lonc1,latc1,vo.precc{ic}{it+1}(ixc,iyc)',[5 5],'b'); 
  scatter(parm.lon(mcsillt4Cl{ic}{it}(:,1)),parm.lat(mcsillt4Cl{ic}{it}(:,2))); 
  %contourf(lonc,latc,vo.LIave{ic}{it}(ixc,iyc)',150,'linestyle','none'); 
%  scatter(parm.lon(mcsillt4Cl{ic}{it+1}(:,1)),parm.lat(mcsillt4Cl{ic}{it+1}(:,2))); 
  title('CPI')
  coast_centered(0)
  xlim([lonc(1) lonc(end)]);
  ylim([latc(1) latc(end)]);
%  xlim([lonc1(xpts/2-3) lonc1(xpts/2+3)]);
%  ylim([latc1(xpts/2-3) latc1(xpts/2+3)]);
  axis square
  hold off
  ip=ip+1;
%%%%%%%%%%%%%%%%%%%%%
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  contourf(lonc,latc,vo.LIave{ic}{it}(ixc,iyc)',150,'linestyle','none'); 
  colorbar; caxis([-6 0]); colormap(hh(ip),flipud(cmap(2))); hold on
  scatter(parm.lon(mcsillt4Cl{ic}{it}(:,1)),parm.lat(mcsillt4Cl{ic}{it}(:,2))); 
  title('LI')
  coast_centered(0)
  xlim([lonc(1) lonc(end)]);
  ylim([latc(1) latc(end)]);
  axis square
  hold off
  ip=ip+1;
%%%%%%%%%%%%%%%%%
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
%  contourf(lonc,latc,vo.theta{ic}{it}(xl,yl,parm.i3km)',150,'linestyle','none'); colorbar; caxis([310 318]); colormap(hh(ip),cmap(2)); hold on
  contourf(lonc,latc,vo.thetaa{ic}{it}(ixc,iyc,parm.i3km)',150,'linestyle','none'); 
  colorbar; caxis([-3 3]); colormap(hh(ip),cmap(3)); hold on; 
  scatter(parm.lon(mcsillt4Cl{ic}{it}(:,1)),parm.lat(mcsillt4Cl{ic}{it}(:,2))); 
  title('\theta anomaly [K]')
  coast_centered(0)
  xlim([lonc(1) lonc(end)]);
  ylim([latc(1) latc(end)]);
  axis square
  hold off
  ip=ip+1;
%%%%%%%%%%%%%%%%%
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  contourf(lonc,latc,vo.pva{ic}{it}(ixc,iyc,parm.i3km)',150,'linestyle','none'); 
  colorbar; caxis([0 2e-6]); colormap(hh(ip),cmap(2)); hold on
  scatter(parm.lon(mcsillt4Cl{ic}{it}(:,1)),parm.lat(mcsillt4Cl{ic}{it}(:,2))); 
  title('potential vorticity anomaly')
  coast_centered(0)
  xlim([lonc(1) lonc(end)]);
  ylim([latc(1) latc(end)]);
  axis square
  hold off
  ip=ip+1;
%%%%%%%%%%%%%%%%%
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip); 
  contourf(lonc,latc,vo.spdt{ic}{it}(ixc,iyc,parm.i3km)',250,'linestyle','none'); colorbar; caxis([-35 35]); colormap(hh(ip),cmap(24)); hold on; 
  title('Q [K/day]')
  scatter(parm.lon(mcsillt4Cl{ic}{it}(:,1)),parm.lat(mcsillt4Cl{ic}{it}(:,2))); 
  coast_centered(0)
  xlim([lonc(1) lonc(end)]);
  ylim([latc(1) latc(end)]);
  axis square
  ip=ip+1;
  hold off;
%%%%%%%%%%%%%%%%%%%%%
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  contourf(lonc,latc,vo.precc{ic}{it+1}(ixc,iyc)',150,'linestyle','none'); 
  colorbar; caxis([0 10]); colormap(hh(ip),cmap(2)); hold on
  scatter(parm.lon(mcsillt4Cl{ic}{it}(:,1)),parm.lat(mcsillt4Cl{ic}{it}(:,2))); 
  title('convective prec')
  coast_centered(0)
  xlim([lonc(1) lonc(end)]);
  ylim([latc(1) latc(end)]);
  axis square
  hold off
  ip=ip+1;
%%%%%%%%%%%%%%%%%%%%%
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  contourf(lonc,latc,vo.precl{ic}{it}(ixc,iyc)',150,'linestyle','none'); 
  colorbar; caxis([0 4]); colormap(hh(ip),cmap(2)); hold on
  scatter(parm.lon(mcsillt4Cl{ic}{it}(:,1)),parm.lat(mcsillt4Cl{ic}{it}(:,2))); 
  title('stratiform prec')
  coast_centered(0)
  xlim([lonc(1) lonc(end)]);
  ylim([latc(1) latc(end)]);
  axis square
  hold off
  ip=ip+1;
%}



%%%%%%%%%%%%%%%%%
  hh(ip)=subplot(nprow,ceil(npanel/nprow),ip); % wake low should be at subsidence zone where high entropy air is brought in by rain induced subsidence
  du = (diff(vo.u_cam{ic}{it}(ixc,iyc,[1,parm.i3km]),1,3));
  dv = (diff(vo.v_cam{ic}{it}(ixc,iyc,[1,parm.i3km]),1,3));
  contourf(lonc,latc,du',250,'linestyle','none'); colorbar
%  contourf(lonc,latc,sqrt(vo.u_cam{ic}{it}(ixc,iyc,1)'.^2+vo.v_cam{ic}{it}(ixc,iyc,1)'.^2)); 
  colorbar; caxis([-30 30]); colormap(hh(ip),cmap(24,0,20)); hold on;
  title('shear [m/s]')
%  contourf(lonc,latc,vo.w{ic}{it}(ixc,iyc,parm.i3km)',25,'linestyle','none'); colorbar; caxis([-0.1 0.1]); colormap(hh(ip),cmap(3)); hold on; 
%  contourf(lonc,latc,vo.w{ic}{it}(ixc,iyc,1)',25,'linestyle','none'); colorbar; caxis([-0.02 0.02]); colormap(hh(ip),cmap(24,0,16)); hold on; 
  maxB = max(max(vo.Bave{ic}{it}(ixc,iyc)));
  contour(lonc,latc,vo.Bave{ic}{it}(ixc,iyc)',[10:10:80],'c','linewidth',3); 
  vo.precc{ic}{it}(isnan(vo.precc{ic}{it}(:)))=0;
  vo.precl{ic}{it}(isnan(vo.precl{ic}{it}(:)))=0;
  contour(lonc,latc,vo.precc{ic}{it}(ixc,iyc)',[5:1.5:50],'r','linewidth',3); 
  contour(lonc,latc,vo.precl{ic}{it}(ixc,iyc)',[1 1],'y','linewidth',3); 
%  contour(lonc,latc,vo.dU03km_cam{ic}{it}(ixc,iyc)',[8:15],'r'); hold on
%  contour(lonc,latc,du',[-25:2:-5],':k','linewidth',2); hold on
%  contour(lonc,latc,du',[5:2:25],'k','linewidth',2); hold on
%  contour(lonc,latc,vo.Bmax{ic}{it}(ixc,iyc)',[250:10:300],'m','linewidth',3); hold on
%  [cprec,hprec] = contour(lonc,latc,vo.prect{ic}{it}(ixc,iyc)'*3.6e6,'k','linewidth',3); hold on; 
  scatter(parm.lon(mcsillt4Cl{ic}{it}(:,1)),parm.lat(mcsillt4Cl{ic}{it}(:,2)),100,'linewidth',2); 
  coast_centered(0)
  xlim([lonc(1) lonc(end)]);
  ylim([latc(1) latc(end)]);
%  hprec.LevelStep=3;
  us = zeros(xpts,xpts);
  
  ustorm = mean(mean(mean(vo.u_cam{ic}{it}(ixc,iyc,1:parm.i6km))));
  us(xpts/2,xpts/2) = ustorm;
%  us(xpts/2,xpts/2) = vo.uclust{ic}{it};
%  us(xpts/2,xpts/2) = vo.ucl{ic}{it};
%  us(xpts/2,xpts/2) = vo.ujet{ic}{it};
%  us(xpts/2,xpts/2) = vo.u_storm{ic}{it};
%  us(xpts/2,xpts/2) = u_storm{ic}{it+1};
  vs = zeros(xpts,xpts);
  vstorm = mean(mean(mean(vo.v_cam{ic}{it}(ixc,iyc,1:parm.i6km))));
  vs(xpts/2,xpts/2) = vstorm;
%  vs(xpts/2,xpts/2) = vo.vclust{ic}{it};
%  vs(xpts/2,xpts/2) = vo.vcl{ic}{it};
%  vs(xpts/2,xpts/2) = vo.vjet{ic}{it};
%  vs(xpts/2,xpts/2) = vo.v_storm{ic}{it};
%  vs(xpts/2,xpts/2) = v_storm{ic}{it+1};
  uscl = 0.10;
%  quiver(lonc,latc,vo.u{ic}{it}(ixc,iyc,10)*uscl,vo.v{ic}{it}(ixc,iyc,10)*uscl,'k','autoscale','off'); 
%  quiver(lonc,latc,vo.u_cam{ic}{it}(ixc,iyc,1)*uscl,vo.v_cam{ic}{it}(ixc,iyc,1)*uscl,'k','autoscale','off'); 
%  quiver(lonc,latc,(vo.u_cam{ic}{it}(ixc,iyc,1)-ustorm)*uscl,(vo.v_cam{ic}{it}(ixc,iyc,1)-vstorm)*uscl,'k','autoscale','off'); 
  ncquiverref(lonc,latc,(vo.u_cam{ic}{it}(ixc,iyc,1)-ustorm)*uscl,(vo.v_cam{ic}{it}(ixc,iyc,1)-vstorm)*uscl,'m/s',1,'true',[0.3 0.3 0.3]); 
%  quiver(lonc,latc,vo.u_stormrel{ic}{it}(ixc,iyc,parm.i3km)*uscl,vo.v_stormrel{ic}{it}(ixc,iyc,parm.i3km)*uscl,'k','autoscale','off'); 
  usscl = 1;
  quiver(lonc,latc,us(ixc,iyc)*usscl,vs(ixc,iyc)*usscl,'g','autoscale','off'); 
  axis square
  hold off;
%saveas(gcf,[diro '/mcs_c' num2str(ic) '_t' sprintf('%02d',it) '.png'])
%end
%disp('done')
%pause(0.5)
%{
contour(mean(cat(3,vo.thetaa{ic}{:}),3),[1:1:15],'k'); hold on; % lons over a fixed latitude
contour(mean(cat(3,vo.thetaa{ic}{:}),3),[-15:1:1],'k:'); hold on; % lons over a fixed latitude
contour(mean(cat(3,vo.pva{ic}{:}),3),[1e-6:1e-6:1e-5],'r'); hold on;
contour(mean(cat(3,vo.pva{ic}{:}),3),[-1e-5:1e-6:-1e-6],'r:'); hold on;
quiver(mean(cat(3,vo.u_stormrel{ic}{:}),3)*uscl,mean(cat(3,vo.w{ic}{:}),3)*wscl,'AutoScaleFactor',1); 
quiver(mean(cat(3,vo.v_stormrel{ic}{:}),3),mean(cat(3,vo.v_stormrel{ic}{:}),3),'AutoScaleFactor',1); hold off;
%}
%quiver(mean(cat(3,vo.w{ic}{:}),3),mean(cat(3,vo.w{ic}{:}),3),'AutoScaleFactor',15); hold off;
%{
for it = 1:nt4Cl(ic)
  subplot(2,1,1)
%  disp(['it, lat, lon = ' num2str(it) ', ' num2str(parm.lat(vo.mcsillt4Cl_fixlat{ic}{it}(1,2))) ', ' num2str(parm.lon(vo.mcsillt4Cl_fixlat{ic}{it}(25,1))'-360)])
%  lons = parm.lon(vo.mcsillt4Cl_fixlat{ic}{it}(:,1))';
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
%  scatter(parm.lon(vo.mcsillt4Cl_fixlat{ic}{it}(:,1)),parm.lat(vo.mcsillt4Cl_fixlat{ic}{it}(1,2))*ones(1,xpts),'b*'); 
%  scatter(parm.lon(vo.mcsillt4Cl_fixlat{ic}{it}(ceil(xpts/2),1)),parm.lat(vo.mcsillt4Cl_fixlat{ic}{it}(1,2)),'r*'); hold off
end
%}
%end
