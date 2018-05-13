
% scatter plot of DTCOND, each idex of (time,lev,lat,lon) is a single
% sample
%ncfile_cam='/bkirtman4/gchen/cesm_spcam/archive/F_2000_4SPCAM2/atm/rest/DTCOND.nc';
%ncfile_cam='/bkirtman4/gchen/cesm_spcam/archive/F_2000_4SPCAM_m200501/atm/rest/F_2000_4SPCAM_m200501.cam.rh0.DTCOND.nc';
%ncfile_spcam='/bkirtman4/gchen/cesm_spcam/archive/spcam_cam_rh0_m2005/spcam_actual_m2005_branch.cam.rh0.DTCOND.nc';
ncfile_spcamo = '/bkirtman4/gchen/cesm_spcam/archive/spcam_cam_rh0_m2005/spcam_actual_m2005_f09f09_branch.cam.rh0.0001-01-14-03600.nc'
ncfile_spcam='/bkirtman4/gchen/cesm_spcam/archive/spcam_cam_rh0_m2005/spcam_actual_m2005_f09f09_branch.cam.rh0.CLD.nc';
%ncfile_spcam='/bkirtman4/gchen/cesm_spcam/archive/spcam_cam_rh0_sam1mom/DTCOND.nc';
ncfile_diff='/bkirtman4/gchen/cesm_spcam/archive/cam_diff_m2005/m2005.diff.cam.rh0.DTCOND.nc';
%
do_scatter = 1;
do_contour = 0;
%
varname={'DTCOND','CLDLOW','CLDHGH','CLDMED'};
%
%var_cam = ncread(ncfile_cam,varname);
var_spcam = ncread(ncfile_spcam,varname{2});
var_diff = ncread(ncfile_diff,varname{1});
%label_cam = [varname ' for cam'];
label_spcam = [varname{2} ' for spcam'];
label_diff = [ varname{1} ' difference between spcam and cam'];
lat = ncread(ncfile_spcamo,'lat');
lon = ncread(ncfile_spcamo,'lon');
lev = ncread(ncfile_spcamo,'ilev');
lev_lmh = [find(abs(lev-680)==min(abs(lev-680))), find(abs(lev-440)==min(abs(lev-440)))]; % the level between low, medium, and high
%
figure;
ilev = 29;
cutoff_lat=45;
ilat = find(lat<=cutoff_lat & lat>=-cutoff_lat);
iloop=1
ii=1;
for itime = 100:700
   time = [itime];
   ind = find(var_spcam(:,ilat,time)==2); % find the cloud indices
   var1tmp=squeeze(var_spcam(:,ilat,time));
   var(:,1)=var1tmp(ind);
   jj=1;
   for ilev = 30:-1:lev_lmh(1)
      var2tmp=squeeze(var_diff(:,ilat,ilev,time));
      var3(:,jj)=var2tmp(ind);
      jj=jj+1;
   end
   hist(var3(var3(:)>=-10^-4 & var3(:)<=10^-4),100)
   xlim([-0.5*10^-3 0.5*10^-3])
%   var(:,2) = reshape(var_spcam(:,ilat,ilev,time),1,numel(var_spcam(:,ilat,ilev,time)));
%   std(var(:,1))
%   std(var(:,2))
%   mean(var(:,1))
%   mean(var(:,2))
   %
   %dscatter(dtcam(:),dtspcam(:),'PLOTTYPE',contour)
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
      xlabel(label_spcam)
      ylabel(label_diff)
      hold off
   end
   %
   if (do_contour)
      subplot(2,1,2)
      contourf(var_diff(:,ilat,ilev,time)','linestyle','none')
      colorbar
      if (iloop == 1)
         cmin2=min(var_spcam(:));
         cmax2=max(var_spcam(:));
      end
      caxis([cmin2 cmax2])
   end
   pause
   iloop=iloop+1
end
