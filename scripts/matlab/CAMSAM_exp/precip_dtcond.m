ncfile{3} = '/bkirtman4/gchen/cesm_spcam/archive/spcam_cam_rh0_m2005/atm/rest/spcam_actual_m2005_f09f09_branch.cam.rh0.0001-01-14-03600.nc'
ncfile{2} = '/bkirtman4/gchen/cesm_spcam/archive/F_2000_4SPCAM_m200501/atm/rest/F_2000_4SPCAM_m200501.cam.rh0.PRECIP.nc' % extra one time sec compared to spcam
%ncfile_spcam='/bkirtman4/gchen/cesm_spcam/archive/spcam_cam_rh0_m2005/spcam_actual_m2005_f09f09_branch.cam.rh0.CLD.nc';
%ncfile_spcam='/bkirtman4/gchen/cesm_spcam/archive/spcam_cam_rh0_sam1mom/DTCOND.nc';
ncfile{1}='/bkirtman4/gchen/cesm_spcam/archive/cam_diff_m2005/m2005.diff.cam.rh0.DTCOND.nc';
%
do_scatter = 1;
do_contour = 0;
%
varname={'DTCOND','PRECC','PRECCDZM','PRECL','PRECSH','PRECT'};
%
%var_cam = ncread(ncfile_cam,varname);
var{1} = ncread(ncfile{1},varname{1});
var{2} = ncread(ncfile{2},varname{2});
%label_cam = [varname ' for cam'];
label{1} = ['DTCOND difference between spcam and cam'];
label{2} = [varname{2} ' for cam'];
lat = ncread(ncfile{1},'lat');
lon = ncread(ncfile{1},'lon');
lev = ncread(ncfile{3},'lev');
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
%   ind = find(var{1}(:,ilat,time)r~=0); % find the cloud indices
   var1tmp=squeeze(var{1}(:,ilat,ilev,1:500));
   var2tmp=squeeze(var{2}(:,ilat,1:500));
%   vartmp(:,1)=var1tmp(ind);
   vartmp(:,1)=var1tmp(:);
   vartmp(:,2)=var2tmp(:);
%   dscatter(vartmp(:,1),vartmp(:,2),'PLOTTYPE','contour')
   if (do_scatter)
      nbins = [100,100];
      [bin_counts bin_centers] = hist3(vartmp,nbins);
      var_axis(:,1)=linspace(min(vartmp(:,1)),max(vartmp(:,1)),nbins(1));
      var_axis(:,2)=linspace(min(vartmp(:,2)),max(vartmp(:,2)),nbins(2));
      plot(vartmp(:,1),vartmp(:,2),'.') % scatter plot without too many
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
         xmin=min(vartmp(:,1)); xmax=max(vartmp(:,1));
         ymin=min(vartmp(:,2)); ymax=max(vartmp(:,2));
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
      contourf(var{2}(:,ilat,ilev,time)','linestyle','none')
      colorbar
      if (iloop == 1)
         cmin2=min(var{1}(:));
         cmax2=max(var{1}(:));
      end
      caxis([cmin2 cmax2])
   end
   pause
   iloop=iloop+1
end
