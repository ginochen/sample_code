% this script does the svd decomposition for the edtcond vertical
% profile. 
run /nethome/gchen/scripts/matlab/startup.m
ncfile{1} = '/bkirtman4/gchen/cesm_spcam/archive/cam_diff_m2005/m2005.diff.cam.rh0.DTCOND.nc';
ncfile{2} = '/bkirtman4/gchen/cesm_spcam/archive/F_2000_4SPCAM_m200501/atm/rest/F_2000_4SPCAM_m200501.cam.rh0.PRECIP.nc' % extra one time sec compared to spcam
ncfile{3} = '/bkirtman4/gchen/cesm_spcam/archive/spcam_cam_rh0_m2005/atm/rest/spcam_actual_m2005_f09f09_branch.cam.rh0.0001-01-14-03600.nc'

  
varname={'DTCOND','PRECC','PRECCDZM','PRECL','PRECSH','PRECT'};

var{1} = ncread(ncfile{1},varname{1});
var{2} = ncread(ncfile{2},varname{2});


lat = ncread(ncfile{1},'lat');
lon = ncread(ncfile{1},'lon');
lev = ncread(ncfile{3},'lev');

cutoff_lat=90;
ilat = find(lat<=cutoff_lat & lat>=-cutoff_lat);

dim = size(var{1}(:,ilat,:,:));
ii=1;
for itime = 1:100:700
   [svec{ii}, sval{ii}, pc{ii}]=svdVar(reshape(var{1}(:,ilat,:,itime),dim(1)*dim(2),dim(3))'); % decompose the vertical modes (that's the reason to transpose the matrix)
   sval_ratio{ii}=sval{ii}/sum(sval{ii})
   ii=ii+1;
end
horzcat(sval_ratio{:}) % show the sval_ratio for all itimes in one single columnwise matrix



% plotting starts here


ii=1;
for j=1:7; for i=1:12; subplot(3,4,i); plot(lev,svec{j}(:,i)); camroll(-90); ylim([-1 1]); end; figure; end % plot the eigenmodes
figure;
for j=1:7; subplot(4,2,j); contourf(reshape(pc{j}(:,1),dim(1),dim(2))');colorbar;caxis([-1.5*10^-3 1.5*10^-3]); end % plot the pc in space associated with an eigenmode
% try to project the entire time series to a given few eigenmodes, see
% if the actual DTCOND in CAM simulation could be perturbed by these
% spatial error modes
