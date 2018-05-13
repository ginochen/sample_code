
% scatter plot of DTCOND, each idex of (time,lev,lat,lon) is a single
% sample
%ncfile_cam='/bkirtman4/gchen/cesm_spcam/archive/F_2000_4SPCAM2/atm/rest/DTCOND.nc';
ncfile{1} = '/bkirtman4/gchen/cesm_spcam/archive/F_2000_4SPCAM_m200501/atm/rest/F_2000_4SPCAM_m200501.cam.rh0.DTCOND.nc';
ncfile{2} = '/bkirtman4/gchen/cesm_spcam/archive/spcam_cam_rh0_m2005/atm/rest/spcam_actual_m2005_f09f09_branch.cam.rh0.DTCOND.nc';
%ncfile_spcam='/bkirtman4/gchen/cesm_spcam/archive/spcam_cam_rh0_sam1mom/DTCOND.nc';
ncfile{3} = '/bkirtman4/gchen/cesm_spcam/archive/cam_diff_m2005/m2005.diff.cam.rh0.DTCOND.nc';
%
do_scatter = 0;
do_hist = 1;
%
varname = 'DTCOND';
label       = [varname ' for cam'];
label_spcam = [varname ' for spcam'];
label_diff  = [varname ' difference between spcam and cam'];
%
spd = 86400; % seconds per day
var{1} = ncread(ncfile{1},varname)*spd;
var{2} = ncread(ncfile{2},varname)*spd;
var{3} = ncread(ncfile{3},varname)*spd;
lat = ncread(ncfile{1},'lat');
lon = ncread(ncfile{1},'lon');
lev = ncread(ncfile{1},'lev');
%
figure;
levels = [0,120,350,820,930,957,976,992];
for i=1:numel(iilev); iilev(i) = find(abs(lev-levels(i)) == min(abs(lev-levels(i)))); end
iilev(1)=0;
cutoff_lat = 35;
ilat = find(lat<=cutoff_lat & lat>=-cutoff_lat);
nbin = 100;
dtcondbin=linspace(-20,20,nbin);
for itime = 100%:700
   for j = 1:numel(levels)-1
   time = [itime];
   dtcond1 = reshape(var{1}(:,ilat,iilev(j)+1:iilev(j+1),time),1,numel(var{1}(:,ilat,iilev(j)+1:iilev(j+1),time)));
   dtcond2 = reshape(var{2}(:,ilat,iilev(j)+1:iilev(j+1),time),1,numel(var{2}(:,ilat,iilev(j)+1:iilev(j+1),time)));
%   var(:,2) = reshape(var_spcam(:,ilat,ilev,time),1,numel(var_spcam(:,ilat,ilev,time)));
%   std(var(:,1))
%   std(var(:,2))
%   mean(var(:,1))
%   mean(var(:,2))
   %
   %dscatter(dtcam(:),dtspcam(:),'PLOTTYPE',contour)
   if (do_hist)
      h(j) = figure;
      subplot(2,1,1)
      histogram(dtcond2,dtcondbin,'Normalization','probability','EdgeColor','none');
%      xlim([min(dtcond2) max(dtcond2)])
      xlim([-20 20])
      ylabel('Probability')
      title(['SPCAM P-level ' num2str(lev(iilev(j)+1)) ' to ' num2str(lev(iilev(j+1)))])
      subplot(2,1,2)
      histogram(dtcond1,dtcondbin,'Normalization','probability','EdgeColor','none');
%      xlim([min(dtcond2) max(dtcond2)])
      xlim([-20 20])
      title(['CAM P-level ' num2str(lev(iilev(j)+1)) ' to ' num2str(lev(iilev(j+1)))])
      ylabel('Probability')
      xlabel('Heating Rate (K/day)')
   end
   end
   savefig(h,['heating_rate_m2005_itime' num2str(itime)])
end
%
iloop=1
%
for itime = 100:700
   if (do_scatter)
      nbins = [100,100];
      [bin_counts bin_centers] = hist3(var,nbins);
      var_axis(:,1)=linspace(min(var(:,1)),max(var(:,1)),nbins(1));
      var_axis(:,2)=linspace(min(var(:,2)),max(var(:,2)),nbins(2));
      %plot(varcam(:),varspcam(:),'.') % scatter plot without too many
                                       % objects
      %contourf(varcam_axis,varspcam_axis,bin_counts',20,'linestyle','none');
      subplot(2,1,1)
      contourf(var_axis(:,1),var_axis(:,2),log(1+bin_counts)',20,'linestyle','none'); colorbar
      % take the log scale of bin_counts (plus one to avoid log(0)=-inf)
      % since the counts for near zero values outnumber the other values. 
      hold on
      % plot the black 1:1 line
      plot([var_axis(1,1) var_axis(end,1)],[var_axis(1,1) var_axis(end,1)],'k')
      if (iloop == 1 )
         xmin=min(var(:,1)); xmax=max(var(:,1));
         ymin=min(var(:,2)); ymax=max(var(:,2));
         cmin=min(min(log(1+bin_counts)));
         cmax=max(max(log(1+bin_counts)));
      end
      xlim([xmin xmax])
      ylim([ymin ymax])
      caxis([cmin cmax])
      xlabel(label_cam)
      ylabel(label_spcam)
      hold off
   end
   %
   subplot(2,1,2)
   contourf(var_diff(:,ilat,ilev,time)','linestyle','none')
   colorbar
   if (iloop == 1)
      cmin2=min(var_spcam(:));
      cmax2=max(var_spcam(:));
   end
   caxis([cmin2 cmax2])
   pause
   iloop=iloop+1
end
