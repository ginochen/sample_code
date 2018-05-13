function getQVTForCAPE(it)
tic
run /nethome/gchen/scripts/matlab/startup.m
%fname
%%%%%%%%%%%%%%%%%%%%%
%      Paths        %
%%%%%%%%%%%%%%%%%%%%%
% CAUTION!!! 
% 'it=1' starts at 03600 instead of 01800,
% due to *.rh0.* needed minus 1 from previous timestep.
% spcam_actual_m2005_f09f09_branch.cam.rh0.*.nc and
% F_2000_4SPCAM_m200501.cam.rh0.*.nc both start from 0001-01-14-03600 unlike 
% spcam_actual_m2005_f09f09_branch.cam.r.*.nc starts from 0001-01-14-01800 
compset='m2005'
spArchive         = ['/projects/rsmas/kirtman/gchen/cesm_spcam/archive/'];
spcamCaser.dir    = [spArchive 'spcam_actual_m2005_f09f09_branch/atm/rest/'];
%Case.dir    = '/projects/rsmas/kirtman/gchen/archive/spcam_actual_m2005_f09f09_branch/atm/rest/';
spcamdir_r        = dir([spcamCaser.dir '/spcam_actual_m2005_f09f09_branch.cam.r.*-*-*-*.nc']);
spcamCaser.name   = {spcamdir_r.name};
spcamtstep        = strrep(strrep(spcamCaser.name(it),'spcam_actual_m2005_f09f09_branch.cam.r.',''),'.nc','') % replace string with blank and obtian the timestep string
fnameo_QV             = [spcamCaser.dir '/crmave/spcam.m2005.crmave.r.' spcamtstep{1} '.QV.nc']; 
fnameo_T             = [spcamCaser.dir '/crmave/spcam.m2005.crmave.r.' spcamtstep{1} '.T.nc']; 
if ( exist(fnameo_QV,'file') | exist(fnameo_T, 'file') )
  error(['file exists for QV or T under path ' spcamCaser.dir '/crmave/'])
end
%%%%%%%%%%%%%%%%%%%%%%%%
%    Load dimensions   %
%%%%%%%%%%%%%%%%%%%%%%%%
CaseDim.name   = 'spcam_actual_m2005_f09f09_branch.cam.rh0.0001-01-14-03600'; % make sure to start at 03600 instead of 01800, and use only this file for the dimensions 
CaseDim.dir    = '/projects/rsmas/kirtman/gchen/cesm_spcam/archive/spcam_actual_m2005_f09f09_branch/atm/rest/';
lat  = loaddim('lat',CaseDim);
lon  = loaddim('lon',CaseDim);
                                % index 1 (TOA pressure) to nlev (surface pressure) 
%%%%%%%%%%%%%%%%%%%%%%%%
%     Parameters       %
%%%%%%%%%%%%%%%%%%%%%%%%
nlat = numel(lat);
nlon = numel(lon);
nx    = 32; % evenly spaced
nz    = 28; % levels of CAM from surface to nz levels (according to cam/crm_physics.F90 line 992)
tic
   qT = load_reshape(spcamCaser, it, 'CRM_QT', [1 1 1], [inf inf inf], nlon, nlat, nx, nz);
   qc = load_reshape(spcamCaser, it, 'CRM_QC', [1 1 1], [inf inf inf], nlon, nlat, nx, nz);
   if (compset == 'm2005')
      qi = 0;
   elseif (compset == 'sam1mom')
      qi = load_reshape(spcamCaser, it, 'CRM_QI', [1 1 1], [inf inf inf], nlon, nlat, nx, nz);
   end
   T  = load_reshape(spcamCaser, it, 'CRM_T',  [1 1 1], [inf inf inf], nlon, nlat, nx, nz);
   mynccreate_write(fnameo_QV,  'QV', flipdim(squeeze(mean(qT-qc-qi,3)),3), lon, lat);
   mynccreate_write(fnameo_T,   'T',  flipdim(squeeze(mean(T,3)),3), lon, lat);
toc
display(['saved under' spcamCaser.dir])


function [var]       = load_reshape(Case, it, varname,  start, count, nlon, nlat, nx, nz); 
   var  = ncread([Case.dir Case.name{it}],varname,start,count); 
   var  = reshape(var,nlon,nlat,nx,nz);

function mynccreate_write(fnameo,vname,vin,lon,lat)
   if (ndims(vin)==3)
      nccreate(fnameo, vname, 'Dimensions',{'lon',size(vin,1),'lat',size(vin,2),'lev',size(vin,3)},'Format','classic');
      ncwrite (fnameo, vname, vin);
      nccreate(fnameo, 'lon','Dimensions',{'lon',length(lon)},'Format','classic');
      ncwrite (fnameo, 'lon', lon);
      nccreate(fnameo, 'lat','Dimensions',{'lat',length(lat)},'Format','classic');
      ncwrite (fnameo, 'lat', lat);
   else
      error('input variable has > 3 dimension')
   end
